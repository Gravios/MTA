
configure_default_args();
MjgER2016_load_data();

Trial = Trials{20};

xyz = preproc_xyz(Trial,'trb');


% xyz hcom 
% markersBody xyz in head FOR | yaw, roll, pitch, yaw/dt, pitch/dt, roll/dt

% geometric algebra
% u.v +

xyz.data(~nniz(xyz),:,:) = 0;
xyz.filter('ButFilter',4,30,'low');

ang = create(MTADang,Trial,xyz);


bang = transform_origin(Trial,xyz,'hcom','nose',{'head_left','head_right'});
bang.roll(isnan(bang.roll))=0;
roll = MTADxyz('data',bang.roll,'sampleRate',xyz.sampleRate);


hcom = xyz(:,'hcom',[1,2,3]);

% GENERATE orthogonal basis, origin: head's center of mass
nvec = xyz(:,'nose',[1,2,3])-hcom;
nvec = bsxfun(@rdivide,nvec,sqrt(sum((nvec).^2,3))); 

tvec = cross(nvec,xyz(:,'head_left',[1,2,3])-hcom,3);
tvec = bsxfun(@rdivide,tvec,sqrt(sum((tvec).^2,3))); 

figure,plot(nvec(:,1,3)),hold('on');plot(tvec(:,1,3))

lvec = cross(nvec,tvec,3);
lvec = bsxfun(@rdivide,lvec,sqrt(sum((lvec).^2,3))); 
figure,plot(nvec(:,1,3)),hold('on');plot(tvec(:,1,3));plot(lvec(:,1,3))

nvec = cat(2,nvec,lvec,tvec);
clear('lvec','tvec')

marBu = copy(xyz);
marBu.data = multiprod(sq(xyz(:,'spine_upper',:)-hcom),nvec,[2],[2,3]);

uvec = circshift(hcom,-1)-circshift(hcom,1);
dnvec(:,1) = dot(uvec,nvec(:,1,:),3);
dnvec(:,2) = dot(uvec,nvec(:,2,:),3);
dnvec(:,3) = dot(uvec,nvec(:,3,:),3);

dnvec = MTADxyz('data',dnvec,'sampleRate',xyz.sampleRate);

figure,
hold('on');
ind = Trial.stc{'w'};
plot3(marBu(ind,1),marBu(ind,2),marBu(ind,3),'.')
ind = Trial.stc{'r'};
plot3(marBu(ind,1),marBu(ind,2),marBu(ind,3),'.r')


figure,
hold('on');
ind = Trial.stc{'r'};
scatter3(marBu(ind,1),marBu(ind,2),marBu(ind,3),20,ang(ind,'hcom','nose',2),'Filled');
colormap('jet');
%scatter3(marBu(ind,1),marBu(ind,2),marBu(ind,3),10,dnvec(ind,2),'Filled');

% mean position in HFOR given trajectory
% 
%  min: 2*x(t|dyaw,dpitch,dl,df,du)-x(t-1|dyaw,dpitch,dl,df,du)-x(t+1|dyaw,dpitch,dl,df,du)
%
% distance between head and upper spine as a function of derivative of the movement of the head
%
% dHS = dx/dt + dHS
% 
% minimize the trajectory error as the rat moves through space changing configurations
%
% dz | dy,dx 