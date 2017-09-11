function transform_rigidBody(Session,varargin);
[display,overwrite] = DefaultArgs(varargin,{false,false});

if ~strcmp(Session.trialName,'all'),
    Session = MTASession.validate(Session);
end

%Session = MTASession('Ed05-20140529','all','ont');
%Session = MTASession('Ed05-20140528');
%Session = MTASession('jg05-20120317');
%Session = MTASession('jg05-20120317');
%Session = MTASession('Ed03-20140625');
%display = false;
%overwrite = false;


xyz = Session.load('xyz');
rb = Session.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);

xyz.addMarker('fhcom',[128,255,128],{{'head_back','head_front',[0,0,1]}},...
               ButFilter(hcom,3,[2]./(Session.xyz.sampleRate/2),'low'));
xyz.addMarker('hcom',[128,255,128],{{'head_back','head_front',[0,0,1]}},hcom);

nz = -cross(xyz(:,'head_back',:)-hcom,xyz(:,'head_left',:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3)));
nm = nz.*20+hcom;
xyz.addMarker('htx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

ny = cross(xyz(:,'htx',:)-hcom,xyz(:,'head_back',:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nm = ny.*20+hcom;
xyz.addMarker('hrx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

nx = cross(xyz(:,'hrx',:)-hcom,xyz(:,'htx',:)-hcom);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx.*20+hcom;    
xyz.addMarker('hbx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

xyz.addMarker('hbt',[128,255,128],{{'head_back','head_front',[0,0,1]}},...
                  genRotatedMarker(xyz,'hbx',45,{'hbx','htx'}));
xyz.addMarker('hbr',[128,255,128],{{'head_back','head_front',[0,0,1]}},...
                  genRotatedMarker(xyz,'hbx',45,{'hbx','hrx'}));    
xyz.addMarker('hbrt',[128,255,128],{{'head_back','head_front',[0,0,1]}},...
                  genRotatedMarker(xyz,'hbr',45,{'hbx','htx'}));
xyz.addMarker('hrt',[128,255,128],{{'head_back','head_front',[0,0,1]}},...
                  genRotatedMarker(xyz,'hrx',45,{'hrx','htx'}));
nhm = {'hcom','hbx','hrx','htx','hbt','hbr','hbrt','hrt'};    


if display,
    ind = 10000;
    figure, daspect([1,1,1])
    hold on,plot3(xyz(ind,7,1),xyz(ind,7,2),xyz(ind,7,3),'.b')
    hold on,scatter3(xyz(ind,7,1),xyz(ind,7,2),xyz(ind,7,3),150,'b')
    hold on,plot3(xyz(ind,5,1),xyz(ind,5,2),xyz(ind,5,3),'.c')
    hold on,scatter3(xyz(ind,5,1),xyz(ind,5,2),xyz(ind,5,3),150,'c')
    hold on,plot3(xyz(ind,6,1),xyz(ind,6,2),xyz(ind,6,3),'.b')
    hold on,plot3(xyz(ind,8,1),xyz(ind,8,2),xyz(ind,8,3),'.b')
    hold on,plot3(xyz(ind,'hbx',1),xyz(ind,'hbx',2),xyz(ind,'hbx',3),'*m')
    hold on,scatter3(xyz(ind,'hbx',1),xyz(ind,'hbx',2),xyz(ind,'hbx',3),150,'m')
    hold on,plot3(xyz(ind,'hrx',1),xyz(ind,'hrx',2),xyz(ind,'hrx',3),'*m')
    hold on,scatter3(xyz(ind,'hrx',1),xyz(ind,'hrx',2),xyz(ind,'hrx',3),150,'m')
    hold on,plot3(xyz(ind,'htx',1),xyz(ind,'htx',2),xyz(ind,'htx',3),'*m')
    hold on,scatter3(xyz(ind,'htx',1),xyz(ind,'htx',2),xyz(ind,'htx',3),150,'m')
    hold on,plot3(xyz(ind,'hcom',1),xyz(ind,'hcom',2),xyz(ind,'hcom',3),'+g')

    hold on,plot3(xyz(ind,'hbt',1),xyz(ind,'hbt',2),xyz(ind,'hbt',3),'+k')
    hold on,plot3(xyz(ind,'hbr',1),xyz(ind,'hbr',2),xyz(ind,'hbr',3),'+k')
    hold on,plot3(xyz(ind,'hbrt',1),xyz(ind,'hbrt',2),xyz(ind,'hbrt',3),'+k')
    hold on,plot3(xyz(ind,'hrt',1),xyz(ind,'hrt',2),xyz(ind,'hrt',3),'+k')
    hold on,plot3(sxyz(ind,'hbx',1),sxyz(ind,'hbx',2),sxyz(ind,'hbx',3),'*m')
    hold on,plot3(sxyz(ind,'hrx',1),sxyz(ind,'hrx',2),sxyz(ind,'hrx',3),'*m')
    hold on,plot3(sxyz(ind,'htx',1),sxyz(ind,'htx',2),sxyz(ind,'htx',3),'*m')
    hold on,plot3(sxyz(ind,'hcom',1),sxyz(ind,'hcom',2),sxyz(ind,'hcom',3),'+g')
end 



FileName_coarse = fullfile(Session.spath,[Session.filebase '.xyz-shift.mat']);
if ~exist(FileName_coarse,'file')||overwrite,

    i = [-100:10:100];
    j = [-100:10:100];
    k = [-100:10:0];

    vxy = xyz.vel('head_front',[1,2]);
    ind = vxy.data>2;
    aind = Session.stc{'a'};
    if ~isempty(aind),
        %aind.resample(xyz);
        aind.cast('TimeSeries');
        aind.resample(xyz);
        ind = aind.data&ind;
    end
    

    vxyz = zeros([7,numel(i),numel(j),numel(k)]);

% $$$ x = 51;
% $$$ y = 51;
% $$$ z = 32; 
% $$$ pool = parpool(10);

    txyz = xyz.copy;
    txyz.data = xyz(:,nhm,:);
    txyz.model = xyz.model.rb(nhm);
    for x = 1:numel(i),tic
        for y = 1:numel(j)
            for z = 1:numel(k)

                sxyz = txyz.copy;
                sxyz.data = bsxfun(@plus,nx*i(x)+ny*j(y)+nz*k(z),sxyz.data);


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
                    vxyz(m,x,y,z) = nanvar(ang(bnds(1)<ang(:,m)&ang(:,m)<bnds(2),m));
                    
                end

            end
        end
        toc
    end

    save(FileName_coarse,'i','j','k','vxyz');
else
    load(FileName_coarse);
end




%% Figure of search
mind = [];
for m = 1:7,
    [mind(m,:),mv] = LocalMinimaN(sq(vxyz(m,:,:,:)),100,100);
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
        ind = aind.data&ind;
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
[ijk{:}] = meshgrid(ni,nj,nk);

% COMPUTE surfaces of thresholded region of the estimation domain
mset = [1:7];
for m = mset
    svxyz = log10(sq(nvxyz(m,:,:,:)));
    iopts = [ijk,{svxyz},{prctile(svxyz(nniz(svxyz(:))),2)}];
    isos(m) = isosurface(iopts{:});
end

if display,
    clist = 'rbgymck';    
    figure();hold('on');
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
    quiver3(mpos(:,1),mpos(:,2),mpos(:,3),varLines(:,1),varLines(:,2),varLines(:,3),30)
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
    bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'head_back',:))
nxyz.data(:,nxyz.model.gmi('head_left'),:) =  ...
    bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'head_left',:))
nxyz.data(:,nxyz.model.gmi('head_front'),:) =  ...
    bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'head_front',:))
nxyz.data(:,nxyz.model.gmi('head_right'),:) =  ...
    bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'head_right',:))

% UPDATE MTADxyz object metadata
nxyz.label = 'trb';
nxyz.key  = 't';
nxyz.name = 'tr_corrected_head';
nxyz.updateFilename(Session);      

% SAVE new MTADxyz object
nxyz.save;


