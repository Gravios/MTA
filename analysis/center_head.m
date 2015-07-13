Trial = MTATrial('Ed05-20140528');
Trial = MTATrial('Ed10-20140814');
%Trial = MTATrial('jg05-20120317');
%stc_mode = 'qda_filtf1p5';
stc_mode = 'auto_wbhr';
Trial.stc.updateMode(stc_mode);Trial.stc.load;


xyz = Trial.load('xyz');


%% Calculate the Center of Mass of the head for each frame
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});

hcom = MTADxyz('data',xyz.com(rb),'sampleRate',xyz.sampleRate);
xyz.addMarker('hcom',[.7,0,.7],{{'head_back','head_front',[0,0,255]}},hcom.data);

hcom.filter('ButFilter',3,2,'low');
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,255]}},hcom.data);


%% Calculate Rotation Matrix: normal of the head plane as axis of rotation

xyz_hb_b = sq(xyz(:,'head_back',:)-xyz(:,'hcom',:));
xyz_hb_r = sq(xyz(:,'head_right',:)-xyz(:,'hcom',:) );
%xyz_hb_l = sq(xyz(:,'head_right',:)-xyz(:,'hcom',:) );
%xyz_hf_r = sq(xyz(:,'head_right',:)-xyz(:,'hcom',:) );
%xyz_hf_l = sq(xyz(:,'head_right',:)-xyz(:,'hcom',:) );

ax_ord = [2,3,1];
head_norm = cross(xyz_hb_b,xyz_hb_r);
head_norm = multiprod(head_norm,1./sqrt(sum(head_norm.^2,2)),2);
j =1:3;
head_norm = head_norm(:,ax_ord);
head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
j = [ 0,-1, 1;...
      1, 0,-1;...
     -1, 1, 0];
k = [1,3,2;...
     3,1,1;...
     2,1,1];
head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1))...
           .*repmat(j,[1,1,size(head_norm,1)]);
rot_ang = deg2rad(45);
head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
              +sin(rot_ang)*head_cpm...
              +(1-cos(rot_ang))*head_kron;

j =1:3;

% Rotated marker;
nmark = permute(sum(head_rotMat.*permute(reshape(xyz_hb_b(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,1,2]);

xyz.addMarker('head_br45',[.7,0,.7],{{'head_back','head_right',[0,0,255]}},permute(nmark,[1,3,2])+hcom.data);


%% END marker rotation within head reference system

%% New head coordinate refrences with equal norms
  Trial = MTATrial('jg05-20120317');

    
xyz = Trial.load('xyz');
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
                  ButFilter(hcom,3,[2]./(Trial.xyz.sampleRate/2),'low'));
xyz.addMarker('hcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
nz = cross(xyz(:,'head_back',:)-hcom,xyz(:,'head_right',:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3)));
nm = nz.*20+hcom;
xyz.addMarker('htx',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},nm);
ny = cross(xyz(:,'htx',:)-hcom,xyz(:,'head_back',:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nm = ny.*20+hcom;
xyz.addMarker('hrx',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},nm);
nx = cross(xyz(:,'hrx',:)-hcom,xyz(:,'htx',:)-hcom);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx.*20+hcom;    
xyz.addMarker('hbx',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},nm);
xyz.addMarker('hbt',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
                  genRotatedMarker(xyz,'hbx',45,{'hbx','htx'}));
xyz.addMarker('hbr',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
                  genRotatedMarker(xyz,'hbx',45,{'hbx','hrx'}));    
xyz.addMarker('hbrt',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
                  genRotatedMarker(xyz,'hbr',45,{'hbx','htx'}));
xyz.addMarker('hrt',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
                  genRotatedMarker(xyz,'hrx',45,{'hrx','htx'}));
    
    
ind = 1000;

figure, daspect([1,1,1])
hold on,plot3(xyz(ind,7,1),xyz(ind,7,2),xyz(ind,7,3),'.b')
hold on,plot3(xyz(ind,5,1),xyz(ind,5,2),xyz(ind,5,3),'.c')
hold on,plot3(xyz(ind,6,1),xyz(ind,6,2),xyz(ind,6,3),'.b')
hold on,plot3(xyz(ind,8,1),xyz(ind,8,2),xyz(ind,8,3),'.b')
hold on,plot3(xyz(ind,'hbx',1),xyz(ind,'hbx',2),xyz(ind,'hbx',3),'*m')
hold on,plot3(xyz(ind,'hrx',1),xyz(ind,'hrx',2),xyz(ind,'hrx',3),'*m')
hold on,plot3(xyz(ind,'htx',1),xyz(ind,'htx',2),xyz(ind,'htx',3),'*m')
hold on,plot3(xyz(ind,'hcom',1),xyz(ind,'hcom',2),xyz(ind,'hcom',3),'+g')
hold on,plot3(xyz(ind,'hbt',1),xyz(ind,'hbt',2),xyz(ind,'hbt',3),'+k')
hold on,plot3(xyz(ind,'hbr',1),xyz(ind,'hbr',2),xyz(ind,'hbr',3),'+k')
hold on,plot3(xyz(ind,'hbrt',1),xyz(ind,'hbrt',2),xyz(ind,'hbrt',3),'+k')
hold on,plot3(xyz(ind,'hrt',1),xyz(ind,'hrt',2),xyz(ind,'hrt',3),'+k')
hold on,plot3(sxyz(ind,'hbx',1),sxyz(ind,'hbx',2),sxyz(ind,'hbx',3),'*m')
hold on,plot3(sxyz(ind,'hrx',1),sxyz(ind,'hrx',2),sxyz(ind,'hrx',3),'*m')
hold on,plot3(sxyz(ind,'htx',1),sxyz(ind,'htx',2),sxyz(ind,'htx',3),'*m')
hold on,plot3(sxyz(ind,'hcom',1),sxyz(ind,'hcom',2),sxyz(ind,'hcom',3),'+g')


i = [-100:5:100];
j = [-100:5:100];
k = [-100:5:100];
ind = Trial.stc{'a'};
vz = [];
vxy = [];
x = 51;
y = 51;
z = 32; 

%pool = parpool(10);

for x = 95:numel(i)
    for y = 1:numel(j)

        nhm = {'hcom','hbx','hrx','htx','hbt','hbr','hbrt','hrt'};

sxyz = xyz.copy;
sxyz.data = xyz(:,nhm,:);
sxyz.model = xyz.model.rb(nhm);
sxyz.data = bsxfun(@plus,nx*i(x)+ny*j(y)+nz*k(z),sxyz.data);

fhcom = zeros([sxyz.size(1),1,3]);
fhcom(nniz(sxyz),:,:) = ButFilter(sxyz(nniz(sxyz),'hcom',:),3,2/(sxyz.sampleRate/2),'low');
sxyz.addMarker('fhcom',[.7,1,.7],{{'hbx','hcom',[0,0,1]}},fhcom);

ang = create(MTADang,Trial,sxyz);
ang.data(~nniz(ang(:,1,2,1)),:,:,:)=0;

vxy(:,x,y) = [var(ang(ind,'hbx','fhcom',3)),...
           var(ang(ind,'hrx','fhcom',3)),...
           var(ang(ind,'htx','fhcom',3)),...
           ...
           var(ang(ind,'hbt','fhcom',3)),...
           var(ang(ind,'hbr','fhcom',3)),...
           var(ang(ind,'hbrt','fhcom',3)),...
           var(ang(ind,'hrt','fhcom',3))];
end
end


ms = nhm(2:end);
hfig = figure(3939);clf
plot(k,log10(vz));
xlabel('z shift from original coordinates (mm)');
ylabel('log10 variance of marker to filtered hcom');
legend(ms{:});
saveas(hfig,fullfile(['/storage/gravio/manuscripts/man2015-jgEd-MoCap/' ...
                      'p20150708/'],[Trial.filebase '.zshift.eps']),'epsc');

hfig = figure(3939);
for m = 1:size(vxy),
    clf,hold on
    imagesc(k,k,log10(sq(vxy(m,:,:))));
    locMin = LocalMinima2(log10(sq(vxy(m,:,:))),10,20);
    locMin = fliplr(k(locMin));
    plot(locMin(1),locMin(2),'*w');
    xlim([k([1,end])]); xlabel('x shift from original coordinates (mm)');
    ylim([k([1,end])]); ylabel('y shift from original coordinates (mm)');
    title([nhm{m+1} ' variance of distance to filtered hcom x: ' num2str(locMin(1)) ' y: ' num2str(locMin(2))]);
saveas(hfig,fullfile(['/storage/gravio/manuscripts/man2015-jgEd-MoCap/' ...
                      'p20150708/'],[Trial.filebase '.' nhm{m+1} '.eps']),'epsc');

end

saveas(hfig,fullfile(['/storage/gravio/manuscripts/man2015-jgEd-MoCap/' ...
                      'p20150708/'],[Trial.filebase '.' nhm{m+1} '.eps']),'epsc');

save(fullfile(['/storage/gravio/manuscripts/man2015-jgEd-MoCap/' ...
               'p20150708/'],[Trial.filebase '.xyz-shift.mat']),'vz','vxy');
