function fet = fet_all(Trial,varargin)
%function fet = fet_all(Trial)
%
% exhaustive set of features derived from the raw data
%
%
[newSampleRate] = DefaultArgs(varargin,{20},1);

fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'all_Features','fet_all','a');                  


%Testing ARGS
%Trial = MTATrial('jg05-20120317');

% XYZ Positions of Markers
xyz = Trial.load('xyz');
nm = xyz.size(2);

% COM lower Body Center of Mass
xyz.addMarker('bcom',...     Name
              [.7,0,.7],...  Color
              {{'spine_lower', 'acom',[0,0,255]},... Sticks to visually connect
               {'pelvis_root', 'acom',[0,0,255]},... new marker to skeleton
               {'spine_middle','acom',[0,0,255]}},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));

% COM head Center of Mass
xyz.addMarker('hcom',...     Name
              [.7,0,.7],...  Color
              {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
               {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
               {'head_front','hcom',[0,0,255]},...
               {'head_right','hcom',[0,0,255]}},... 
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));

              
% COM subject Center of Mass              
xyz.addMarker('acom',...    Name
              [.7,0,.7],... Color
              {{'spine_lower', 'acom',[0,0,255]},... Sticks to visually connect
               {'pelvis_root', 'acom',[0,0,255]},... new marker to skeleton
               {'spine_middle','acom',[0,0,255]},...
               {'spine_upper', 'acom',[0,0,255]},...
               {'head_back',   'acom',[0,0,255]},...
               {'head_front',  'acom',[0,0,255]}},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_front'})));


xyz.resample(newSampleRate);


% XYZ filtered 
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');


% FVELXY Filtered marker speeds in XY plane
fvelxy = xyz.vel([],[1,2]);
fvelxy.filter('ButFilter',3,2.5,'low');
fvelxy.data(fvelxy.data<0)=.1;

facxy = fvelxy.copy;
facxy.data = [diff(facxy.data);zeros([1,fvelxy.size(2)])];

fvelxy.data = log10(fvelxy.data);



% FVELZ Filtered marker speeds in Z axis
fvelz = fxyz.vel([],[3]);
fvelz.filter('ButFilter',3,2.5,'low');
fvelz.data(fvelz.data<0)=.1;

facz = fvelz.copy;
facz.data = [diff(facz.data);zeros([1,fvelz.size(2)])];

fvelz.data = log10(fvelz.data);

% ANG InterMarker Spherical Coordinates
ang = create(MTADang,Trial,xyz);

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);



%% Derivatives of features

% ANG (theta) between segments 
d = 2;
cang = ang.copy;
cang.data = [circ_dist(ang(:,1,2,d),ang(:,2,3,d)),...
             circ_dist(ang(:,1,3,d),ang(:,3,4,d)),...
             circ_dist(ang(:,1,5,d),ang(:,5,7,d)),...
             circ_dist(ang(:,1,nm+1,d),ang(:,nm+1,4,d)),...
             circ_dist(ang(:,2,3,d),ang(:,3,4,d)),...
             circ_dist(ang(:,2,4,d),ang(:,2,3,d)),...
             circ_dist(ang(:,3,4,d),ang(:,4,7,d)),...
             circ_dist(ang(:,3,5,d),ang(:,5,7,d)),...
             circ_dist(ang(:,4,5,d),ang(:,5,7,d)),...
             circ_dist(ang(:,4,7,d),ang(:,2,3,d)),...
             circ_dist(ang(:,5,7,d),ang(:,2,3,d)),...
             circ_dist(ang(:,nm+1,nm+2,d),ang(:,2,4,d)),...
             circ_dist(ang(:,nm+1,nm+3,d),ang(:,nm+1,nm+2,d))];
cang.filter('ButFilter',3,1.4,'low');
% Add to feature matrix
cang.resample(xyz);
fet.data =  [fet.data,cang.data];


% DRV of pitch
d = 2;
cang = ang.copy;
cang.data = [ang(:,1,2,d),...
             ang(:,1,3,d),...
             ang(:,1,5,d),...
             ang(:,1,nm+1,d),...
             ang(:,2,3,d),...
             ang(:,2,4,d),...
             ang(:,3,4,d),...
             ang(:,3,5,d),...
             ang(:,4,5,d),...
             ang(:,4,7,d),...
             ang(:,5,7,d),...
             ang(:,nm+1,nm+2,d),...
             ang(:,nm+1,nm+3,d),...
             ang(:,nm+2,nm+3,d)];
cang.filter('ButFilter',3,5,'low');
cang.data = diff(cang.data);
% Add to feature matrix
cang.resample(xyz);
fet.data =  [fet.data,cang.data];



%% DRV of Yaw
d = 1;
cang = ang.copy;
cang.data = [ang(:,1,2,d),...
             ang(:,1,3,d),...
             ang(:,1,5,d),...
             ang(:,1,nm+1,d),...
             ang(:,2,3,d),...
             ang(:,2,4,d),...
             ang(:,3,4,d),...
             ang(:,3,5,d),...
             ang(:,4,5,d),...
             ang(:,4,7,d),...
             ang(:,5,7,d),...
             ang(:,nm+1,nm+2,d),...
             ang(:,nm+1,nm+3,d),...
             ang(:,nm+2,nm+3,d)];
cang.filter('ButFilter',3,5,'low');
cang.data = circ_dist(cang.data,circshift(cang.data,1)).*cang.sampleRate;
cang.data = circshift(permute(sum(cang.segs(1:cang.size(1),newSampleRate,0).^2),...
                              [2,3,1]),...
                      [-round(newSampleRate/2),0]);
cang.data(cang.data<1e-6) = 1e-6;
cang.data = log10(cang.data);
% Add to feature matrix
cang.resample(xyz);
fet.data =  [fet.data,cang.data];



% PPC feature
man = Trial.load('fet','lsppc').resample(xyz);


% RHM feature
[rhm,fs] = fet_rhm(Trial,[],'mtchglong',true);
rhm.data = median(rhm(:,fs>6&fs<14),2);
% Add to feature matrix
rhm.resample(xyz);



fet.data =  [fet.data,...
             man.data,...
             rhm.data,...
             fxyz(:,[1:5],3),...
             fvelxy(:,{'spine_lower','spine_upper','head_front','bcom','hcom','acom'}),...
             fvelz(:,{'spine_lower','spine_upper','head_front','bcom','hcom','acom'})];







