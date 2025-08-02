function cag = req20151127(Trial)
%Trial = MTATrial('Ed03-20140624');
% XYZ Positions of Markers
xyz = Trial.load('xyz');


% COM Body Center of Mass
rbb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle'});
xyz.addMarker('hcom',[.7,0,.7],{{'head_back','head_front',[0,0,255]}},xyz.com(rbb));

rba = xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_front'});
xyz.addMarker('acom',[.7,0,.7],{{'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_front',[0,0,255]}},xyz.com(rba));


% FXYZ Filtered Marker Positions {Low Pass 2.5 Hz}
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,.5,'low');


xyz.addMarker('ffh',[.7,0,.7],{},fxyz(:,7,:));
xyz.addMarker('fbh',[.7,0,.7],{},fxyz(:,5,:));

% fxyz = xyz.copy;
% fxyz.filter('ButFilter',3,2.5,'low');


% ANG InterMarker Spherical Coordinates
ang = create(MTADang,Trial,xyz);

% 
% %figure,plot(circ_dist(circshift(ang(:,10,12,1),-2),circshift(ang(:,10,12,1),2)))
% 
% % FANG Filtered Intermarker angles 
% fang = create(MTADang,Trial,fxyz);
% 
% % FVEL Filtered marker speeds in XY plane
% fvel = xyz.vel([],[1,2]);
% fvel.filter('ButFilter',3,2.5,'low');
% fvel.data(fvel.data<0)=.1;
% 
% 
% fac = fvel.copy;
% fac.data = [diff(fac.data);zeros([1,fvel.size(2)])];
% 
% 
% fvel.data = log10(fvel.data);
% 
% % UVEL 
% uvel = xyz.vel([],[3]);
% uvel.filter('ButFilter',3,2.5,'low');
% uvel.data(uvel.data<0)=.1;
% 
% uac = uvel.copy;
% uac.data = [diff(uac.data);zeros([1,uvel.size(2)])];
% 
% uvel.data = log10(uvel.data);

%% End Var Setup

cag = MTADang('data',circ_dist(ang(:,5,7,1),ang(:,11,12,1)),'sampleRate',ang.sampleRate);
cag.filter('ButFilter',3,[1,10],'bandpass')
cag = cag.data;
% 
% msync = Trial.xyz.sync.copy;
% msync.data = msync.sync.data;
% man = MTADfet(Trial.spath,Trial.filebase,...
%                circ_dist(ang(:,5,7,1),ang(:,11,12,1)),...
%                Trial.xyz.sampleRate,...
%                msync,...
%                msync(1),...
%                [],[],[],'casting','hcast','c');
% man.filter('ButFilter',3,[1,10],'bandpass')
% %man.save;

