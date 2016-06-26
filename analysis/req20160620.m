function req20160620(Session,varargin);
[display,overwrite] = DefaultArgs(varargin,{false,false});

Session = MTASession.validate(Session);

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



FileName_coarse = fullfile(Session.spath,[Session.filebase,'-req20160620','.xyz-shift.mat']);
if ~exist(FileName_coarse,'file')||overwrite,

    i = [-100:10:100];
    j = [-100:10:100];
    k = [-100:10:0];


    ind = Session.stc{'a'};
    if isempty(ind),
        ind=':';
    end

    % Create binned 
    oang = create(MTADang,Session,xyz);
    eds = linspace(-pi/2,pi/2,40);
    %figure,bar(eds,histc(oang(ind,5,7,2),eds),'histc')
    gHeadPitch = oang(ind,5,7,2);
    [~,ainds] = histc(gHeadPitch,eds);


    vxyz = zeros([7,numel(i),numel(j),numel(k),numel(eds)]);


    txyz = xyz.copy;
    txyz.data = xyz(:,nhm,:);
    txyz.model = xyz.model.rb(nhm);
    for x = 1:numel(i),tic
        for y = 1:numel(j)
            for z = 1:numel(k)
                for b = unique(ainds)'

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

                    minds = b==ainds;
                    for m = 1:7,
                        bnds = prctile(ang(minds,m),[.1,99.1]);
                        vxyz(m,x,y,z,b) = nanvar(ang(bnds(1)<ang(minds,m)&ang(minds,m)<bnds(2),m));
                        
                    end
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


    %% Fine grain search
    ni = [ni(1):2:ni(2)];
    nj = [nj(1):2:nj(2)];
    nk = [nk(1):2:nk(2)];
    ind = Session.stc{'a'}.cast('TimeSeries').resample(xyz);
    if isempty(ind),
        ind = ':';
    else
        ind.data = logical(ind.data);
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

ijk = cell(1,3);
[ijk{:}] = meshgrid(ni,nj,nk);

clist = 'rbgymck';

if display,
    figure
    for m =1:7,
        svxyz = log10(sq(nvxyz(m,:,:,:)));
        iopts = [ijk,{svxyz},{prctile(svxyz(:),2)}];
        p = patch(isosurface(iopts{:}));
        iopts(end) = {p};
        isonormals(iopts{:});
        p.FaceColor = clist(m);
        p.EdgeColor = 'none';
    end
    daspect([1 1 1]); 
    camlight; lighting phong
end


%hbx
% $$$ mind = [];
% $$$ figure,hold on
% $$$ for i = 1:3,
% $$$ [mind,mv] = LocalMinimaN(sq(nvxyz(i,:,:,:)),100,100);
% $$$ plot3(ni(mind(1)),nj(mind(2)),nk(mind(3)),'+g');
% $$$ end
% $$$ for i = 4:7
% $$$ [mind,mv] = LocalMinimaN(sq(nvxyz(i,:,:,:)),100,100);
% $$$ plot3(ni(mind(1)),nj(mind(2)),nk(mind(3)),'+m');
% $$$ end

i =1;
%ni(mind(i,1)),nj(mind(i,2)),nk(mind(i,3))

xyz.addMarker('nhb',[128,255,128],...
              {{'spine_upper','nhb',[0,1,0]},...
               {'nhb'        ,'nhr',[1,0,0]},...
               {'nhb'        ,'nhl',[0,0,1]}},...
              bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_back',:)));
xyz.addMarker('nhl',[128,255,128],...
              {{'nhl','nhf',[0,0,1]}},...
              bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_left',:)));
xyz.addMarker('nhf',[128,255,128],...
              {},...
              bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_front',:)));
xyz.addMarker('nhr',[128,255,128],...
              {{'nhr','nhf',[1,0,0]}},...
              bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_right',:)));

nxyz = Session.load('xyz');
nxyz.data(:,nxyz.model.gmi('head_back'),:) =  bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_back',:));
nxyz.data(:,nxyz.model.gmi('head_left'),:) =  bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_left',:));
nxyz.data(:,nxyz.model.gmi('head_front'),:) =  bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_front',:));
nxyz.data(:,nxyz.model.gmi('head_right'),:) =  bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_right',:));


nxyz.label = 'trs';
nxyz.key  = 's';
nxyz.name = 'trb_state_corrected_head';
nxyz.updateFilename(Session);      
nxyz.save;          



% $$$ xyz = Session.load('xyz');
% $$$ ang = create(MTADang,Session,xyz);
% $$$ nang = create(MTADang,Session,ixyz);
% $$$ 
% $$$ rbn = xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper','nhb','nhl','nhf','nhr'});
% $$$ nxyz = Session.xyz.copy;
% $$$ nxyz.data = xyz(:,rbn.ml,:);
% $$$ nxyz.model = rbn;
% $$$ 
% $$$ rbhn = xyz.model.rb({'nhb','nhl','nhf','nhr'});
% $$$ txyz = Trial.xyz.copy;
% $$$ txyz.addMarker('hcom',[128,255,128],{{'head_back','head_front',[0,0,1]}},nxyz.com(rbhn));
% $$$ txyz.filter('ButFilter',3,2,'low');
% $$$ nxyz.addMarker('fhcom',[128,128,128],...
% $$$                {{'fhcom','nhb',[0,1,0]},...
% $$$                 {'fhcom','nhl',[0,0,1]},...
% $$$                 {'fhcom','nhr',[1,0,0]},...
% $$$                 {'fhcom','nhf',[0,1,0]}},...
% $$$                txyz(:,end,:));
% $$$ 
% $$$ ang =  create(MTADang,Trial,xyz);
% $$$ ang.filter('ButFilter',1,20,'low');
% $$$ nang = create(MTADang,Trial,nxyz);
% $$$ nang.filter('ButFilter',1,20,'low');
% $$$ 
% $$$ txyz = xyz.copy;
% $$$ tnxyz = nxyz.copy;
% $$$ 
% $$$ xyz = txyz.copy;
% $$$ nxyz = tnxyz.copy;
% $$$ %xyz.filter('ButFilter',3,5,'low');
% $$$ %nxyz.filter('ButFilter',3,5,'low');
% $$$ 
% $$$ xyv = xyz.vel([1,7],[1,2]);
% $$$ xyv.filter('ButFilter',3,2.4,'low');
% $$$ xyv.data(xyv.data<.001) = 0.001;
% $$$ xyn = nxyz.vel([1,7],[1,2]);
% $$$ xyn.filter('ButFilter',3,2.4,'low');
% $$$ xyn.data(xyn.data<.001) = 0.001;
% $$$ 
% $$$ ind = Trial.stc{'w'};
% $$$ eds = linspace(-3,2,100);
% $$$ figure,
% $$$ subplot(211)
% $$$ hist2(log10([xyn(ind,1),xyn(ind,2)]),eds,eds);
% $$$ subplot(212)
% $$$ hist2(log10([xyv(ind,1),xyv(ind,2)]),eds,eds);
% $$$ 
% $$$ ind = Trial.stc{'w'};
% $$$ figure,hold on
% $$$ ha = bar(eds,histc(log10(xyn(ind,2)),eds),'histc');
% $$$ ha.FaceColor = 'r';
% $$$ ha.FaceAlpha = .5;
% $$$ hs = bar(eds,histc(log10(xyv(ind,2)),eds),'histc');
% $$$ hs.FaceColor = 'c';
% $$$ hs.FaceAlpha = .5;
% $$$ 
% $$$ figure,
% $$$ for i = 1:1000,    
% $$$     cla
% $$$     xlim([-500,500])
% $$$     ylim([-500,500])
% $$$     zlim([0,350])
% $$$     plotSkeleton(Trial,nxyz,i,'line',nang);
% $$$     pause(.01);
% $$$ end
% $$$ 
% $$$ 
% $$$ [ys,fs,ts] =mtchglong(WhitenSignal(diff([nang(nniz(xyz),'fhcom','nhb',3),ang(nniz(xyz),'fhcom','head_back',3)]),[],1),2^8,ang.sampleRate,2^7,2^7*.875,[],[],[],[1,20]);
% $$$ 
% $$$ sp = [];
% $$$ figure,
% $$$ sp(1)=subplot(2,1,1);imagesc(ts,fs,log10(ys(:,:,1,1))'),axis xy, colormap jet,caxis([-5,-2.3])
% $$$ sp(2)=subplot(2,1,2);imagesc(ts,fs,log10(ys(:,:,2,2))'),axis xy, colormap jet,caxis([-5,-2.3])
% $$$ linkaxes(sp,'xy')


% $$$ ms = nhm(2:end);
% $$$ hfig = figure(3939);clf
% $$$ plot(k,log10(vz));
% $$$ xlabel('z shift from original coordinates (mm)');
% $$$ ylabel('log10 variance of marker to filtered hcom');
% $$$ legend(ms{:});
% $$$ saveas(hfig,fullfile(['/storage/gravio/manuscripts/man2015-jgEd-MoCap/' ...
% $$$                       'p20150708/'],[Trial.filebase '.zshift.eps']),'epsc');
% $$$ 
% $$$ hfig = figure(3939);
% $$$ for m = 1:size(vxy),
% $$$     clf,hold on
% $$$     imagesc(k,k,log10(sq(vxy(m,:,:))));
% $$$     locMin = LocalMinima2(log10(sq(vxy(m,:,:))),10,20);
% $$$     locMin = fliplr(k(locMin));
% $$$     plot(locMin(1),locMin(2),'*w');
% $$$     xlim([k([1,end])]); xlabel('x shift from original coordinates (mm)');
% $$$     ylim([k([1,end])]); ylabel('y shift from original coordinates (mm)');
% $$$     title([nhm{m+1} ' variance of distance to filtered hcom x: ' num2str(locMin(1)) ' y: ' num2str(locMin(2))]);
% $$$ saveas(hfig,fullfile(['/storage/gravio/manuscripts/man2015-jgEd-MoCap/' ...
% $$$                       'p20150708/'],[Trial.filebase '.' nhm{m+1} '.eps']),'epsc');
% $$$ 
% $$$ end
% $$$ 
% $$$ saveas(hfig,fullfile(['/storage/gravio/manuscripts/man2015-jgEd-MoCap/' ...
% $$$                       'p20150708/'],[Trial.filebase '.' nhm{m+1} '.eps']),'epsc');
% $$$ 
% $$$ 
% $$$ %mkdir('/storage/gravio/manuscripts/man2015-jgEd-MoCap/p20150716/');
% $$$ %save(fullfile(Trial.spath,[Trial.filebase '.xyz-shift.mat']),'i','j','k','vxyz');
% $$$ load(fullfile(Trial.spath,[Trial.filebase '.xyz-shift.mat']));
% $$$ 
% $$$ 
% $$$ 
% $$$ rb_h_labels = {'head_back','head_left','head_front','head_right'};
% $$$ rb_h_labels = {'nhb','nhl','nhf','nhr'};
% $$$ hxyz = Trial.xyz.copy;
% $$$ hxyz.model = xyz.model.rb(rb_h_labels);;
% $$$ hxyz.data = xyz(:,rb_h_labels,:);
% $$$ hxyz.addMarker('hcom',[128,255,128],{{'nhb','nhf',[0,0,255]}},xyz.com(hxyz.model));
% $$$ 
% $$$ fhcom = zeros([hxyz.size(1),1,hxyz.size(3)]);
% $$$ fhcom(nniz(xyz),:,:) = ButFilter(hxyz(nniz(hxyz),'hcom',:),3,2/(hxyz.sampleRate/2),'low');
% $$$ hxyz.addMarker('fhcom',[128,255,128],{{'nhb','nhf',[0,0,1]}},fhcom);
% $$$ 
% $$$ 
% $$$ 
% $$$ ang = create(MTADang,Trial,hxyz);
% $$$ 
% $$$ [ys,fs,ts] =mtchglong(WhitenSignal([ang(nniz(xyz),'nhb','fhcom',3),ang(nniz(xyz),'nhr','fhcom',3)],[],1),2^8,ang.sampleRate,2^7,2^7*.875,[],[],[],[1,40]);
% $$$ 
% $$$ sp = [];
% $$$ figure,
% $$$ sp(1)=subplot(2,1,1);
% $$$ imagesc(ts,fs,log10(ys(:,:,1,1))'),axis xy, colormap jet,
% $$$ sp(2)=subplot(2,1,2);
% $$$ imagesc(ts,fs,log10(ys(:,:,2,2))'),axis xy, colormap jet,
% $$$ linkaxes(sp,'xy')
% $$$ 
% $$$ hvel = hxyz.vel({'hcom'},[1,2]);
% $$$ ovel = xyz.vel({'hcom'},[1,2]);
% $$$ 
% $$$ figure, plot(ovel.data),
% $$$ hold on,plot(hvel.data)
% $$$ Lines(Trial.stc{'w'}(:),[],'b');
% $$$ 
% $$$ hvel.filter('ButFilter',3,2.5,'low');
% $$$ ovel.filter('ButFilter',3,2.5,'low');
% $$$ 
% $$$ figure, plot(ovel.data),
% $$$ hold on,plot(hvel.data)
% $$$ Lines(Trial.stc{'w'}(:),[],'b');