function transform_rigidBody(Session,varargin);
%function transform_rigidBody(Session,varargin);
%
% Estimates the position of the neck by locating the position relative to the head wihch
% minimizes translations of head movements.
%
% varargin:
%   display       (logical)               false
%   overwrite     (logical)               false
%   scoreFunction (struct)                struct('fun',@mad,'args',{{1}})
%   periods       (numeric Nx2)           []
%


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('display',               false,                                                 ...
                 'overwrite',             false,                                                 ...
                 'scoreFunction',         struct('fun',@mad,'args',{{1}}),                       ...
                 'periods',               []                                                     ...
);
[display,overwrite,scoreFunction,periods] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

% VALIDATE session
if ischar(Session),
% LOAD session
    Session = MTASession.validate(Session);
else
% CONVERT trial to session
    trialPeriodsByDefault = false;
    if isempty(periods), 
        periods = Session.sync.copy();         
        trialPeriodByDefault = true;
    end    
    Session = MTASession.validate([Session.name,'.',Session.maze.name,'.all']);
    if trialPeriodsByDefault,         
        periods.data = periods.data-Session.xyz.origin;
        periods.resample(Session.xyz.sampleRate);
        periods = periods.data;
    end
end

if display,hfig = figure();end

% LOAD xyz
% ADD lowpass filtered marker of the head's center of mass
xyz = Session.load('xyz');
rb = Session.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('fhcom',[0.5,1,0.5],[],ButFilter(hcom,3,[2]./(Session.xyz.sampleRate/2),'low'));
xyz.addMarker('hcom', [0.5,1,0.5],{{'head_back','head_front',[0,0,1]}},hcom);



% GENERATE orthogonal basis, origin: head's center of mass
nz = -cross(xyz(:,'head_back',:)-hcom,xyz(:,'head_left',:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3))); 
nm = nz.*20+hcom;
xyz.addMarker('htx',  [0.5,1,0.5],[],nm);

% GENERATE orthogonal basis, origin: head's center of mass
ny = cross(xyz(:,'htx',:)-hcom,xyz(:,'head_back',:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nm = ny.*20+hcom;
xyz.addMarker('hrx',  [0.5,1,0.5],[],nm);

% GENERATE orthogonal basis, origin: head's center of mass
nx = cross(xyz(:,'hrx',:)-hcom,xyz(:,'htx',:)-hcom);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx.*20+hcom;
xyz.addMarker('hbx',  [0.5,1,0.5],[],nm);

% GENERATE rotated orthogonal basis, origin: head's center of mass
xyz.addMarker('hbt',  [0.5,1,0.5],[],genRotatedMarker(xyz,'hbx',45,{'hbx','htx'}));
xyz.addMarker('hbr',  [0.5,1,0.5],[],genRotatedMarker(xyz,'hbx',45,{'hbx','hrx'}));    
xyz.addMarker('hbrt', [0.5,1,0.5],[],genRotatedMarker(xyz,'hbr',45,{'hbx','htx'}));
xyz.addMarker('hrt',  [0.5,1,0.5],[],genRotatedMarker(xyz,'hrx',45,{'hrx','htx'}));
nhm = {'hcom','hbx','hrx','htx','hbt','hbr','hbrt','hrt'};    


if display,
    clf();
    ind = xyz.get_pose_index();
    subplot2(7,7,[1:4],[1:3]);daspect([1,1,1]);hold('on');
    m = 'head_front';
    hold on,plot3(xyz(ind,m,1),xyz(ind,m,2),xyz(ind,m,3),'.b')
    hold on,scatter3(xyz(ind,m,1),xyz(ind,m,2),xyz(ind,m,3),150,'b')
    m = 'head_back';    
    hold on,plot3(xyz(ind,5,1),xyz(ind,5,2),xyz(ind,5,3),'.','Color',[0.5,0,0.5])
    hold on,scatter3(xyz(ind,5,1),xyz(ind,5,2),xyz(ind,5,3),150,[0.5,0,0.5])
    m = 'head_left';    
    hold on,plot3(xyz(ind,6,1),xyz(ind,6,2),xyz(ind,6,3),'.g')
    m = 'head_right';        
    hold on,plot3(xyz(ind,8,1),xyz(ind,8,2),xyz(ind,8,3),'.r')
    
    hold on,plot3(xyz(ind,'hbx',1),xyz(ind,'hbx',2),xyz(ind,'hbx',3),'^m')
    hold on,scatter3(xyz(ind,'hbx',1),xyz(ind,'hbx',2),xyz(ind,'hbx',3),150,'m')
    hold on,plot3(xyz(ind,'hrx',1),xyz(ind,'hrx',2),xyz(ind,'hrx',3),'^m')
    hold on,scatter3(xyz(ind,'hrx',1),xyz(ind,'hrx',2),xyz(ind,'hrx',3),150,'m')
    hold on,plot3(xyz(ind,'htx',1),xyz(ind,'htx',2),xyz(ind,'htx',3),'^m')
    hold on,scatter3(xyz(ind,'htx',1),xyz(ind,'htx',2),xyz(ind,'htx',3),150,'m')
    hold on,plot3(xyz(ind,'hcom',1),xyz(ind,'hcom',2),xyz(ind,'hcom',3),'+g')
end 



FileName_coarse = fullfile(Session.spath,[Session.filebase,'.xyz-shift_',func2str(scoreFunction.fun),'.mat']);
if ~exist(FileName_coarse,'file') || overwrite,

% SET initial search domain
    i = [-100:10:100];
    j = [-100:10:100];
    k = [-100:10:100];

% COMPUTE head speed
    vxy = xyz.vel('head_front',[1,2]);
% RESTRICT computations to periods where head speed is greater than 2cm/2        
    ind = vxy.data>2;
    
    if ~isempty(periods)
% RESTRICT computations to specified data subset
        pind = false([size(vxy,1),1]);
        for p = 1:size(periods,1),
            pind(periods(p,1):periods(p,2)) = true;
        end
        ind = ind&pind;
    else
% RESTRICT computations to good periods
        aind = Session.stc{'a'};
        if ~isempty(aind),
            %aind.resample(xyz);
            aind.cast('TimeSeries');
            aind.resample(xyz);
            aind.data(isnan(aind.data))=0;
            ind = logical(aind.data)&ind;
        end
    end


    vxyz = zeros([7,numel(i),numel(j),numel(k)]);
    txyz = xyz(:,nhm,:);
    ind = ind&nniz(txyz);
    txyz = txyz(ind,:,:);
    tnx = nx(ind,:,:);
    tny = ny(ind,:,:);   
    tnz = nz(ind,:,:);
    sampleCount = sum(ind);
    for x = 1:numel(i),tic
        disp(['x: ',num2str(x),' of ' num2str(numel(i))]);
        for y = 1:numel(j)
            for z = 1:numel(k)
                sxyz = bsxfun(@plus,tnx*i(x)+tny*j(y)+tnz*k(z),txyz);
                sxyz(:,1,:) = ButFilter(sxyz(:,1,:),3,2/(xyz.sampleRate/2),'low');
                ang = mat2cell(sq(reshape(bsxfun(@minus,sxyz(:,2:8,:),...
                                                 sxyz(:,1,:)),[],1,3)),sampleCount*7,[1,1,1]);
                [~,~,ang] = cart2sph(ang{:});
                ang = reshape(ang,[],7);
                for m = 1:7,
                    bnds = prctile(ang(:,m),[1,99]);
                    vxyz(m,x,y,z) = scoreFunction.fun(ang(bnds(1)<ang(:,m)&ang(:,m)<bnds(2),m),...
                                                      scoreFunction.args{:});
                end%for m
            end%for z
        end%for y
        
toc,end%for x

    save(FileName_coarse,'i','j','k','vxyz');
else
    load(FileName_coarse);
end




%% Figure of search
mind = [];
for m = 1:7,
    [mind(m,:),mv] = LocalMinimaN(sq(vxyz(m,:,:,:)),100,100);
end

if display,
    for m = 1:7,
        subplot2(7,7,5,m);
        imagesc(i,j,log10(sq(vxyz(m,:,:,mind(m,3))))');
        axis xy;
        title([nhm{m+1},':z: ' num2str(k(mind(m,3)))]);
        xlabel('x shift mm');
        ylabel('y shift mm');
    end
end


% SET Domain of estimate search for isosurface
ijk = cell(1,3);
[ijk{:}] = meshgrid(j,i,k); % Why did I have to flip ni and nj

% COMPUTE surfaces of thresholded region of the estimation domain
mset = [1:7];
for m = mset
    svxyz = log10(sq(vxyz(m,:,:,:)));
    iopts = [ijk,{svxyz},{prctile(svxyz(nniz(svxyz(:))),2)}];
    isos(m) = isosurface(iopts{:});
end

if display
    subplot2(7,7,[1:4],[5:7]);daspect([1,1,1]);hold('on');        
    clist = 'rbgymck';    
    hold('on');
    for m = mset
        p = patch(isos(m));
        iopts(end) = {p};
        isonormals(iopts{:});
        p.FaceColor = clist(m);
        p.EdgeColor = 'none';
        p.FaceAlpha = 0.3;
    end
    daspect([1 1 1]); 
    camlight; lighting phong
end
    
    

% FIND the center of each field
mpos = zeros([7,3]);
mind = zeros([7,3]);
for m = 1:7,
    [mind(m,:),mv] = LocalMinimaN(sq(vxyz(m,:,:,:)),100,100);
    mpos(m,:) = [j(mind(m,2)),i(mind(m,1)),k(mind(m,3))]; % Annoying flip
end

% COMPUTE the vector of the elipsoid field's long axis
varLines = zeros([7,3]);
for m = mset
    [~,S,V] = svd(bsxfun(@minus,isos(m).vertices,mpos(m,:)),0);
    varLines(m,:) = V(:,1);    
end

if display,
    quiver3(mpos(:,1),mpos(:,2),mpos(:,3),varLines(:,1),varLines(:,2),varLines(:,3),3);
end

scp = nan([7,7,3]);
tcq = nan([7,7,3]);
for p = 1:7,
    for q = p+1:7,
        u = varLines(p,:);
        v = varLines(q,:);
        w = mpos(p,:)-mpos(q,:);
        a = dot(u,u);
        b = dot(u,v);
        c = dot(v,v);
        d = dot(u,w);
        e = dot(v,w);

        sc = (b*e-c*d)/(a*c-b^2);
        tc = (a*e-b*d)/(a*c-b^2);

        scp(p,q,:) = mpos(p,:)+sc*varLines(p,:);
        tcq(p,q,:) = mpos(q,:)+tc*varLines(q,:);
    end
end

% COMPUTE cluster mean
% LOAD raw position data
% ADD estimate of neck ( new hcom )
mpoint = nanmean(nanmean(cat(3,reshape(scp,[],3),reshape(tcq,[],3)),3));
nxyz = Session.load('xyz');
nxyz.addMarker('head_neck',...
               [1,0,1],...
               {{'head_back','head_front',[1,0,1]}},...
               bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'hcom',:))...
);


% REPLACE data of each head marker with shifted position around neck
nxyz.data(:,nxyz.model.gmi('head_back'),:) =  ...
    bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'head_back',:));
nxyz.data(:,nxyz.model.gmi('head_left'),:) =  ...
    bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'head_left',:));
nxyz.data(:,nxyz.model.gmi('head_front'),:) =  ...
    bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'head_front',:));
nxyz.data(:,nxyz.model.gmi('head_right'),:) =  ...
    bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'head_right',:));

% UPDATE MTADxyz object metadata
% SAVE new MTADxyz object
nxyz.label = 'trb';
nxyz.key  = 't';
nxyz.name = 'tr_corrected_head';
nxyz.updateFilename(Session);      
nxyz.save;


if display,
    scatter3(mpoint(:,1),mpoint(:,2),mpoint(:,3),200,'filled');
    saveas(hfig,fullfile(Session.spath,'figures',[mfilename,'_diagnostic.fig']),'fig');
end


% $$$ FileName_fine = fullfile(Session.spath,[Session.filebase '.xyz-shift_fine_',func2str(scoreFunction.fun),'.mat']);
% $$$ 
% $$$ if ~exist(FileName_fine,'file')||overwrite,
% $$$ 
% $$$     nrange = [min(mind)-2;max(mind)+2]';
% $$$ 
% $$$     nrange(1,:) = clip(nrange(1,:),1,numel(i));
% $$$     nrange(2,:) = clip(nrange(2,:),1,numel(j));
% $$$     nrange(3,:) = clip(nrange(3,:),1,numel(k));
% $$$         
% $$$     ni = i(nrange(1,:));
% $$$     nj = j(nrange(2,:));
% $$$     nk = k(nrange(3,:));
% $$$ 
% $$$ % REFINE search
% $$$     ni = [ni(1):2:ni(2)];
% $$$     nj = [nj(1):2:nj(2)];
% $$$     nk = [nk(1):2:nk(2)];
% $$$ 
% $$$ % COMPUTE head speed
% $$$     vxy = xyz.vel('head_front',[1,2]);
% $$$ % RESTRICT computations to periods where head speed is greater than 2cm/2
% $$$     ind = vxy.data>2;
% $$$     
% $$$     if ~isempty(periods)
% $$$ % RESTRICT computations to specified data subset
% $$$         pind = false([size(vxy,1),1]);
% $$$         for p = 1:size(periods,1),
% $$$             pind(periods(p,1):periods(p,2)) = true;
% $$$         end
% $$$         ind = ind&pind;
% $$$     else
% $$$ % RESTRICT computations to good periods        
% $$$         aind = Session.stc{'a'};
% $$$         if ~isempty(aind),
% $$$             %aind.resample(xyz);
% $$$             aind.cast('TimeSeries');
% $$$             aind.resample(xyz);
% $$$             aind.data(isnan(aind.data))=0;        
% $$$             ind = logical(aind.data) & ind;
% $$$         end
% $$$     end    
% $$$ 
% $$$ 
% $$$     nhm = {'hcom','hbx','hrx','htx','hbt','hbr','hbrt','hrt'};    
% $$$     nvxyz = zeros([7,numel(ni),numel(nj),numel(nk)]);
% $$$     txyz = xyz(:,nhm,:);
% $$$     ind = ind&nniz(txyz);    
% $$$     txyz = txyz(ind,:,:);
% $$$     tnx = nx(ind,:,:);
% $$$     tny = ny(ind,:,:);   
% $$$     tnz = nz(ind,:,:);
% $$$     sampleCount = sum(ind);   
% $$$     for x = 1:numel(ni)
% $$$         disp(['x: ' num2str(x) ' of ' num2str(numel(ni))]),tic
% $$$         for y = 1:numel(nj)
% $$$             for z = 1:numel(nk)
% $$$                 sxyz = bsxfun(@plus,tnx*ni(x)+tny*nj(y)+tnz*nk(z),txyz);                
% $$$                 sxyz(:,1,:) = ButFilter(sxyz(:,1,:),3,2/(xyz.sampleRate/2),'low');
% $$$                 ang = mat2cell(sq(reshape(bsxfun(@minus,sxyz(:,2:8,:),...
% $$$                                                  sxyz(:,1,:)),[],1,3)),sampleCount*7,[1,1,1]);
% $$$                 [~,~,ang] = cart2sph(ang{:});
% $$$                 ang = reshape(ang,[],7);
% $$$ 
% $$$                 for m = 1:7,
% $$$                     bnds = prctile(ang(:,m),[1,99]);
% $$$                     nvxyz(m,x,y,z) = scoreFunction.fun(ang(bnds(1)<ang(:,m)&ang(:,m)<bnds(2),m),...
% $$$                                                        scoreFunction.args{:});
% $$$                 end%for m
% $$$             end%for z
% $$$         end%for y
% $$$         toc
% $$$     end%for x
% $$$     save(FileName_fine,'ni','nj','nk','nvxyz');
% $$$ else
% $$$     load(FileName_fine);
% $$$ end


% $$$ % LOCATE local minima of each variance field
% $$$ mind = [];
% $$$ for m = 1:7,
% $$$     [mind(m,:),mv] = LocalMinimaN(sq(nvxyz(m,:,:,:)),100,100);
% $$$ end
% $$$ 
% $$$ for m = 1:7,
% $$$     subplot2(7,7,6,m);
% $$$     imagesc(ni,nj,log10(sq(nvxyz(m,:,:,mind(m,3))))');
% $$$     caxis;axis xy;
% $$$     title([nhm{m+1},':z: ' num2str(nk(mind(m,3)))]);
% $$$     xlabel('x shift mm');
% $$$     ylabel('y shift mm');
% $$$ end
% $$$ 
% $$$ 
% $$$ % SET Domain of estimate search for isosurface
% $$$ ijk = cell(1,3);
% $$$ [ijk{:}] = meshgrid(nj,ni,nk); % Why did I have to flip ni and nj
% $$$ 
% $$$ % COMPUTE surfaces of thresholded region of the estimation domain
% $$$ mset = [1:7];
% $$$ for m = mset
% $$$     svxyz = log10(sq(nvxyz(m,:,:,:)));
% $$$     iopts = [ijk,{svxyz},{prctile(svxyz(nniz(svxyz(:))),3)}];
% $$$     isos(m) = isosurface(iopts{:});
% $$$ end
% $$$ 
% $$$ if display,
% $$$     subplot2(7,7,[1:4],[5:7]);daspect([1,1,1]);hold('on');    
% $$$     clist = 'rbgymck';    
% $$$     hold('on');
% $$$     for m = mset
% $$$         p = patch(isos(m));
% $$$         iopts(end) = {p};
% $$$         isonormals(iopts{:});
% $$$         p.FaceColor = clist(m);
% $$$         p.EdgeColor = 'none';
% $$$         p.FaceAlpha = 0.3;
% $$$     end
% $$$     daspect([1 1 1]); 
% $$$     camlight; lighting phong
% $$$ end
% $$$ 
% $$$ 
% $$$ % FIND the center of each field
% $$$ mpos = zeros([7,3]);
% $$$ mind = zeros([7,3]);
% $$$ for m = 1:7,
% $$$     [mind(m,:),mv] = LocalMinimaN(sq(nvxyz(m,:,:,:)),100,100);
% $$$     mpos(m,:) = [nj(mind(m,2)),ni(mind(m,1)),nk(mind(m,3))]; % Annoying flip
% $$$ end
% $$$ 
% $$$ % COMPUTE the vector of the elipsoid field's long axis
% $$$ varLines = zeros([7,3]);
% $$$ for m = mset
% $$$     [~,S,V] = svd(bsxfun(@minus,isos(m).vertices,mpos(m,:)),0);
% $$$     varLines(m,:) = V(:,1);    
% $$$ end
% $$$ 
% $$$ if display,
% $$$     quiver3(mpos(:,1),mpos(:,2),mpos(:,3),varLines(:,1),varLines(:,2),varLines(:,3),3);
% $$$ end
% $$$ 
% $$$ % COMPUTE pairwise coordinate pairs along vector orthoginal to lines
% $$$ scp = nan([7,7,3]);
% $$$ tcq = nan([7,7,3]);
% $$$ for p = 1:7,
% $$$     for q = p+1:7,
% $$$         u = varLines(p,:);
% $$$         v = varLines(q,:);
% $$$         w = mpos(p,:)-mpos(q,:);
% $$$         a = dot(u,u);
% $$$         b = dot(u,v);
% $$$         c = dot(v,v);
% $$$         d = dot(u,w);
% $$$         e = dot(v,w);
% $$$ 
% $$$         sc = (b*e-c*d)/(a*c-b^2);
% $$$         tc = (a*e-b*d)/(a*c-b^2);
% $$$ 
% $$$         scp(p,q,:) = mpos(p,:)+sc*varLines(p,:);
% $$$         tcq(p,q,:) = mpos(q,:)+tc*varLines(q,:);
% $$$     end
% $$$ end

% COMPUTE cluster mean
% LOAD raw position data
% ADD estimate of neck ( new hcom )
% $$$ mpoint = nanmean(nanmean(cat(3,reshape(scp,[],3),reshape(tcq,[],3)),3));
% $$$ nxyz = Session.load('xyz');
% $$$ nxyz.addMarker('head_neck',...
% $$$                [1,0,1],...
% $$$                {{'head_back','head_front',[1,0,1]}},...
% $$$                bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'hcom',:))...
% $$$ );
% $$$ 
% $$$ if display,
% $$$     scatter3(mpoint(:,1),mpoint(:,2),mpoint(:,3),200,'filled');
% $$$     saveas(hfig,fullfile(Session.spath,'figures',[mfilename,'_diagnostic.fig']));
% $$$     delete(hfig);
% $$$ end
% $$$ 
% $$$ 
% $$$ % REPLACE data of each head marker with shifted position around neck
% $$$ nxyz.data(:,nxyz.model.gmi('head_back'),:) =  ...
% $$$     bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'head_back',:));
% $$$ nxyz.data(:,nxyz.model.gmi('head_left'),:) =  ...
% $$$     bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'head_left',:));
% $$$ nxyz.data(:,nxyz.model.gmi('head_front'),:) =  ...
% $$$     bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'head_front',:));
% $$$ nxyz.data(:,nxyz.model.gmi('head_right'),:) =  ...
% $$$     bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'head_right',:));
% $$$ 
% $$$ % UPDATE MTADxyz object metadata
% $$$ % SAVE new MTADxyz object
% $$$ nxyz.label = 'trb';
% $$$ nxyz.key  = 't';
% $$$ nxyz.name = 'tr_corrected_head';
% $$$ nxyz.updateFilename(Session);      
% $$$ nxyz.save;
% $$$ 
% $$$ 
