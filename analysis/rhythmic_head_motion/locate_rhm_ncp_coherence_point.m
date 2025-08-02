function locate_rhm_ncp_coherence_point(Session,varargin);


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('display',false,'overwrite',false);
[display,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% DEFARGS ------------------------------------------------------------------------------------------

% $$$ if ~strcmp(Session.trialName,'all'),
% $$$     Session = MTASession.validate(Session);
% $$$ end
Session = MTASession.validate([Session.name,'.',Session.maze.name,'.all']);

ncp = fet_ncp(Session);

xyz = Session.load('xyz','trb');


rb = Session.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);

xyz.addMarker('fhcom',[128,255,128],{{'head_back','head_front',[0,0,1]}},...
               ButFilter(hcom,3,[2]./(Session.xyz.sampleRate/2),'low'));
xyz.addMarker('hcom',[128,255,128],{{'head_back','head_front',[0,0,1]}},hcom);

% GENERATE orthogonal basis, origin: head's center of mass
nz = -cross(xyz(:,'head_back',:)-hcom,xyz(:,'head_left',:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3)));
nm = nz.*20+hcom;
xyz.addMarker('htx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);
% GENERATE orthogonal basis, origin: head's center of mass
ny = cross(xyz(:,'htx',:)-hcom,xyz(:,'head_back',:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nm = ny.*20+hcom;
xyz.addMarker('hrx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);
% GENERATE orthogonal basis, origin: head's center of mass
nx = cross(xyz(:,'hrx',:)-hcom,xyz(:,'htx',:)-hcom);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx.*20+hcom;    
xyz.addMarker('hbx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);
% GENERATE rotated orthogonal basis, origin: head's center of mass
xyz.addMarker('hbt',  [0.5,1,0.5],[],genRotatedMarker(xyz,'hbx',45,{'hbx','htx'}));
xyz.addMarker('hbr',  [0.5,1,0.5],[],genRotatedMarker(xyz,'hbx',45,{'hbx','hrx'}));    
xyz.addMarker('hbrt', [0.5,1,0.5],[],genRotatedMarker(xyz,'hbr',45,{'hbx','htx'}));
xyz.addMarker('hrt',  [0.5,1,0.5],[],genRotatedMarker(xyz,'hrx',45,{'hrx','htx'}));
nhm = {'hcom','hbx','hrx','htx','hbt','hbr','hbrt','hrt'};    

fxyz = xyz.copy();
fxyz.filter('RectFilter',3,4);

%FileName_coarse = fullfile(Session.spath,[Session.filebase '.rhmXncp-shift.mat']);
FileName_coarse = fullfile(Session.spath,[Session.filebase '.rhmXncp-shift_v2.mat']);
if ~exist(FileName_coarse,'file') || overwrite,

    i = [-30:15:30];
    j = [-30:15:30];
    k = [-150:5:150];

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
    ind = MTADepoch([],[],ind,xyz.sampleRate,xyz.sync,xyz.origin,'TimeSeries');
    
    vxyz = zeros([7,numel(i),numel(j),numel(k)]);


    txyz = fxyz(:,nhm,:);
    for x = 1:numel(i),tic
        for y = 1:numel(j)
            for z = 1:numel(k)
                sxyz = bsxfun(@plus,nx*i(x)+ny*j(y)+nz*k(z),txyz);
                sxyz(nniz(sxyz),1,:) = ButFilter(sxyz(nniz(sxyz),1,:),3,2/(xyz.sampleRate/2),'low');
                ang = mat2cell(sq(reshape(bsxfun(@minus,sxyz(:,2:8,:),sxyz(:,1,:)),[],1,3)),size(sxyz,1)*7,[1,1,1]);
                [~,~,ang] = cart2sph(ang{:});
                ang = MTADfet.encapsulate(Session,...
                                          diff(reshape(ang,[],7)),...
                                          fxyz.sampleRate,'rhythmic head motion','rhm','r');
                
                ang.filter('RectFilter');
                ang.data = cat(1,zeros([1,numel(nhm)-1]),diff(ang.data),zeros([1,numel(nhm)-1]));
                ang.filter('RectFilter');
                
                for m = 1:numel(nhm)-1,
                    % Modify time stamps and spec; add padding (0's)
                    % COMPUTE cross spectra
                    [ys,fs,ts] = mtcsdglong([ang(:,m),ncp.data],2^8,fxyz.sampleRate,2^7,2^6,[],[],[],[5,14]);
                    if m ==1,
                        ts = ts+(2^6/2)/fxyz.sampleRate;
                        ssr = 1/diff(ts(1:2));
                        pad = round([ts(1),size(fxyz,1)./fxyz.sampleRate-ts(end)].*ssr)-[1,0];
                        szy = size(ys);
                        if x==1&&y==1&&z==1,
                            resample(ind,ssr);
                            if size(ys,1)>size(ind,1),
                                ind.data = cat(1,ind.data,zeros([size(ys,1)-size(ind,1),1]));
                            elseif size(ys,1)<size(ind,1),
                                ind.data = ind.data(1:size(ys,1));
                            end
                            ind = ind.data;
                        end
                    end 
                    yss = sq(mean(cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),2));
                    vxyz(m,x,y,z) = nanmean(abs(yss(ind,1,2)))./nanmean(sqrt(yss(ind,1,1).*yss(ind,2,2)));
                end%for m
            end%for z
        end%for y
        toc
    end%for x

    save(FileName_coarse,'i','j','k','vxyz');
else
    load(FileName_coarse);
end




%% Figure of search
mind = [];
for m = 1:7,
    [mind(m,:),mv] = LocalMinimaN(-sq(vxyz(m,:,:,:)),100,100);
end

% $$$ figure,
% $$$ clf
% $$$ for m = 1:7,
% $$$     subplot(2,4,m);
% $$$     imagesc(i,j,log10(sq(vxyz(m,:,:,mind(m,3))))');
% $$$     axis xy;
% $$$     title([nhm{m+1} ', with z-axis shift: ' num2str(k(mind(m,3)))]);
% $$$     xlabel('x-axis shift from origin mm');
% $$$     ylabel('y-axis shift from origin mm');
% $$$ end
% $$$ suptitle(Session.filebase)
% $$$ 


FileName_fine = fullfile(Session.spath,[Session.filebase '.xyz-shift_fine.mat']);

if ~exist(FileName_fine,'file')||overwrite,

    nrange = [min(mind)-2;max(mind)+2]';

    nrange(1,:) = clip(nrange(1,:),1,numel(i));
    nrange(2,:) = clip(nrange(2,:),1,numel(j));
    nrange(3,:) = clip(nrange(3,:),1,numel(k));
        
    ni = i(nrange(1,:));
    nj = j(nrange(2,:));
    nk = k(nrange(3,:));

% REFINE search
    ni = [ni(1):2:ni(2)];
    nj = [nj(1):2:nj(2)];
    nk = [nk(1):2:nk(2)];

    vxy = xyz.vel('head_front',[1,2]);
    ind = vxy.data>2;
    aind = Session.stc{'a'};
    if ~isempty(aind),
        %aind.resample(xyz);
        aind.cast('TimeSeries');
        aind.resample(xyz);
        aind.data(isnan(aind.data))=0;        
        ind = logical(aind.data) & ind;
    end
    
    nvxyz = zeros([7,numel(ni),numel(nj),numel(nk)]);

    nhm = {'hcom','hbx','hrx','htx','hbt','hbr','hbrt','hrt'};
    txyz = xyz.copy;
    txyz.data = xyz(:,nhm,:);
    txyz.model = xyz.model.rb(nhm);
    for x = 1:numel(ni)
        disp(['x: ' num2str(x)]),tic
        for y = 1:numel(nj)
            for z = 1:numel(nk)

                sxyz = txyz.copy;
                sxyz.data = bsxfun(@plus,nx*ni(x)+ny*nj(y)+nz*nk(z),sxyz.data);

                fhcom = zeros([sxyz.size(1),1,3]);
                fhcom(nniz(sxyz),:,:) = ButFilter(sxyz(nniz(sxyz),'hcom',:),3,2/(sxyz.sampleRate/2),'low');
                sxyz.addMarker('fhcom',[128,255,128],{{'hbx','hcom',[0,0,1]}},fhcom);

                sxyz.data = sxyz(ind,:,:);
                ang = [sxyz(:,'hbx',:)-sxyz(:,'fhcom',:);...
                       sxyz(:,'hrx',:)-sxyz(:,'fhcom',:);...
                       sxyz(:,'htx',:)-sxyz(:,'fhcom',:);...
                       ...
                       sxyz(:,'hbt',:)-sxyz(:,'fhcom',:);...
                       sxyz(:,'hbr',:)-sxyz(:,'fhcom',:);...
                       sxyz(:,'hbrt',:)-sxyz(:,'fhcom',:);...
                       sxyz(:,'hrt',:)-sxyz(:,'fhcom',:)];


                ang = mat2cell(permute(ang,[1,3,2]),size(ang,1),[1,1,1]);
                [~,~,ang] = cart2sph(ang{:});
                ang = reshape(ang,[],7);
                
                for m = 1:7,
                    bnds = prctile(ang(:,m),[.1,99.1]);
                    nvxyz(m,x,y,z) = nanvar(ang(bnds(1)<ang(:,m)&ang(:,m)<bnds(2),m));
                    
                end
            end
        end
        toc
    end
    save(FileName_fine,'ni','nj','nk','nvxyz');
else
    load(FileName_fine);
end

%save(fullfile(Session.spath,[Session.filebase '.xyz-shift_fine_a-m-s.mat']),'ni','nj','nk','nvxyz');
%load(fullfile(Session.spath,[Session.filebase '.xyz-shift_fine_a-m-s.mat']),'ni','nj','nk','nvxyz');

% LOCATE local minima of each variance field
mind = [];
for m = 1:7,
    [mind(m,:),mv] = LocalMinimaN(sq(nvxyz(m,:,:,:)),100,100);
end

%figure,
% $$$ clf
% $$$ for m = 1:7,
% $$$     subplot(2,4,m);
% $$$     imagesc(ni,nj,log10(sq(nvxyz(m,:,:,mind(m,3))))');
% $$$     caxis;axis xy;
% $$$     title([nhm{m+1} ', with z-axis shift: ' num2str(nk(mind(m,3)))]);
% $$$     xlabel('x-axis shift from origin mm');
% $$$     ylabel('y-axis shift from origin mm');
% $$$ end
% $$$ suptitle(Session.filebase)

% SET Domain of estimate search for isosurface
ijk = cell(1,3);
[ijk{:}] = meshgrid(nj,ni,nk); % Why did I have to flip ni and nj

% COMPUTE surfaces of thresholded region of the estimation domain
mset = [1:7];
for m = mset
    svxyz = log10(sq(nvxyz(m,:,:,:)));
    iopts = [ijk,{svxyz},{prctile(svxyz(nniz(svxyz(:))),2)}];
    isos(m) = isosurface(iopts{:});
end

if display,
    clist = 'rbgymck';    
    hfig = figure();
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
    [mind(m,:),mv] = LocalMinimaN(sq(nvxyz(m,:,:,:)),100,100);
    mpos(m,:) = [nj(mind(m,2)),ni(mind(m,1)),nk(mind(m,3))]; % Annoying flip
end

% COMPUTE the vector of the elipsoid field's long axis
varLines = zeros([7,3]);
for m = mset
    [~,S,V] = svd(bsxfun(@minus,isos(m).vertices,mpos(m,:)),0);
    varLines(m,:) = V(:,1);    
end

if display,
    figure(hfig);
    quiver3(mpos(:,1),mpos(:,2),mpos(:,3),varLines(:,1),varLines(:,2),varLines(:,3),3);
    saveas(hfig,fullfile(Session.spath,[mfilename,'-fine','.fig']));
    delete(hfig);
end



% COMPUTE pairwise coordinate pairs along vector orthoginal to lines
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
% $$$ plot3(scp(1),scp(2),scp(3),'*b')
% $$$ plot3(tcq(1),tcq(2),tcq(3),'*m')

% COMPUTE cluster mean
mpoint = nanmean(nanmean(cat(3,reshape(scp,[],3),reshape(tcq,[],3)),3));


% $$$ % ADD estimate of neck 
% $$$ xyz.addMarker(['head_neck'],[255,0,255],{{'head_back','head_front',[255,0,255]}},bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'hcom',:)));
% $$$ 
% $$$ % ADD estimate of head center
% $$$ i = 1;
% $$$ xyz.addMarker(['head_center'],[255,0,255],{{'head_back','head_front',[255,0,255]}},bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'hcom',:)));

% $$$ xyz.addMarker('nhb',[128,255,128],...
% $$$               {{'spine_upper','nhb',[0,1,0]},...
% $$$                {'nhb'        ,'nhr',[1,0,0]},...
% $$$                {'nhb'        ,'nhl',[0,0,1]}},...
% $$$               bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_back',:)));
% $$$ xyz.addMarker('nhl',[128,255,128],...
% $$$               {{'nhl','nhf',[0,0,1]}},...
% $$$               bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_left',:)));
% $$$ xyz.addMarker('nhf',[128,255,128],...
% $$$               {},...
% $$$               bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_front',:)));
% $$$ xyz.addMarker('nhr',[128,255,128],...
% $$$               {{'nhr','nhf',[1,0,0]}},...
% $$$               bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_right',:)));

% LOAD raw position data
nxyz = Session.load('xyz');

% $$$ nxyz.data(:,nxyz.model.gmi('head_back'),:) =  bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_back',:));
% $$$ nxyz.data(:,nxyz.model.gmi('head_left'),:) =  bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_left',:));
% $$$ nxyz.data(:,nxyz.model.gmi('head_front'),:) =  bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_front',:));
% $$$ nxyz.data(:,nxyz.model.gmi('head_right'),:) =  bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_right',:));

% ADD estimate of neck ( new hcom )
nxyz.addMarker('head_neck',...
               [255,0,255],...
               {{'head_back','head_front',[255,0,255]}},...
               bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'hcom',:))...
);

% ADD estimate of head center
% $$$ i = 1;
% $$$ nxyz.addMarker('head_center',...
% $$$                [255,0,255],...
% $$$                {{'head_back','head_front',[255,0,255]}},...
% $$$                bsxfun(@plus,nx*ni(mind(i,1))...
% $$$                       +ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'hcom',:))...
% $$$ );

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
nxyz.label = 'trb';
nxyz.key  = 't';
nxyz.name = 'tr_corrected_head';
nxyz.updateFilename(Session);      

% SAVE new MTADxyz object
nxyz.save;


