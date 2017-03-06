function [fet,featureTitles,featureDesc] = fet_all2(Trial,varargin)
%function fet = fet_all(Trial)
%
% exhaustive set of features derived from the raw data
%
%
[newSampleRate,RefTrial,procOpts] = DefaultArgs(varargin,{20,[],'SPLINE_SPINE_HEAD_EQD'},1);

fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'all_Features','fet_all','a');                  

% Loads preprocessed version of xyz
nm = Trial.xyz.model.N;
[xyz,ss] = preproc_xyz(Trial,procOpts);              
xyz.resample(newSampleRate);
ss.resample(xyz);

% FXYZ filtered 
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');

% FVELXY Filtered marker speeds in XY plane
fvelxy = xyz.vel([],[1,2]);
fvelxy.data(fvelxy.data<1e-4)=1e-4;
fvelxy.data = log10(fvelxy.data);

% FVELZ Filtered marker speeds in Z axis
fvelz = fxyz.vel([],[3]);
fvelz.data(fvelz.data<1e-4)=1e-4;
fvelz.data = log10(fvelz.data);

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);


%% Derivatives of features

% ANG (theta) between segments 
d = 2;
cang = fang.copy;
cang.data = [fang(:,'spine_lower','pelvis_root',2),...          1
             fang(:,'spine_lower','spine_middle',2),...         2
             fang(:,'spine_lower','bcom',2),...                 3
             fang(:,'spine_lower','hcom',2),...                 4
             fang(:,'pelvis_root','spine_middle',2),...         5
             fang(:,'pelvis_root','spine_upper',2),...          6
             fang(:,'spine_middle','spine_upper',2),...         7
             fang(:,'spine_middle','hcom',2),...                8
             fang(:,'spine_upper','hcom',2),...                 9
             fang(:,'spine_upper','hcom',2),...                 10
             fang(:,'bcom','hcom',2),...                        11
             fang(:,'bcom','acom',2)];                         %12
fet.data =  [fet.data,cang.data];


% DRV of pitch
cang = fang.copy;
cang.data = [ang(:,'spine_lower','pelvis_root',2),...          13
             ang(:,'spine_lower','spine_middle',2),...         14
             ang(:,'spine_lower','bcom',2),...                 15
             ang(:,'spine_lower','hcom',2),...                 16
             ang(:,'pelvis_root','spine_middle',2),...         17
             ang(:,'pelvis_root','spine_upper',2),...          18
             ang(:,'spine_middle','spine_upper',2),...         19
             ang(:,'spine_middle','hcom',2),...                20
             ang(:,'spine_upper','hcom',2),...                 21
             ang(:,'bcom','hcom',2),...                        22
             ang(:,'bcom','acom',2),...                        23
             ang(:,'hcom','acom',2)];                         %24
cang.data = [zeros([1,size(cang,2)]);diff(cang.data)];
cang.data = circ_dist(circshift(cang.data,-1),circshift(cang.data,1)).*cang.sampleRate;
cang.data = circshift(permute(sum(cang.segs(1:cang.size(1),round(newSampleRate/2),0).^2),...
                              [2,3,1]),...
                      [round(newSampleRate/4),0]);
cang.data(cang.data<1e-8) = 1e-8;
cang.data = log10(cang.data);
fet.data =  [fet.data,cang.data]; % Add to feature matrix


%% DRV of Yaw
cang = fang.copy;
cang.data = [ang(:,'spine_lower','pelvis_root',1),...          25
             ang(:,'spine_lower','spine_middle',1),...         26
             ang(:,'spine_lower','bcom',1),...                 27
             ang(:,'spine_lower','hcom',1),...                 28
             ang(:,'pelvis_root','spine_middle',1),...         29
             ang(:,'pelvis_root','spine_upper',1),...          30
             ang(:,'spine_middle','spine_upper',1),...         31
             ang(:,'spine_middle','hcom',1),...                32
             ang(:,'spine_upper','hcom',1),...                 33
             ang(:,'spine_upper','hcom',1),...                 34
             ang(:,'head_back','hcom',1),...                   35
             ang(:,'bcom','hcom',1),...                        36
             ang(:,'bcom','acom',1),...                        37
             ang(:,'hcom','acom',1)];                         %38
cang.data = circ_dist(circshift(cang.data,-1),circshift(cang.data,1)).*cang.sampleRate;
cang.data = circshift(permute(sum(cang.segs(1:cang.size(1),round(newSampleRate/2),0).^2),...
                              [2,3,1]),...
                      [round(newSampleRate/4),0]);
cang.data(cang.data<1e-8) = 1e-8;
cang.data = log10(cang.data);
fet.data =  [fet.data,cang.data];


%% PPC feature
try
    man = Trial.load('fet','lsppc');
catch err
    gen_fet_lsppc(Trial);    
    man = Trial.load('fet','lsppc');
end
man.filter('ButFilter',3,2,'low');
man.resample(fxyz);


%% BFET
mxyz = Trial.load('xyz');
rb = mxyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
hcom = mxyz.com(rb);
mxyz.addMarker('fbcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[1]./(mxyz.sampleRate/2),'low'));
mxyz.addMarker('bcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
mxyz.addMarker('fsl',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(mxyz(:,'spine_lower',:),4,[1.5]./(mxyz.sampleRate/2),'low'));
mxyz.addMarker('fsm',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(mxyz(:,'spine_middle',:),4,[1.5]./(mxyz.sampleRate/2),'low'));
mxyz.addMarker('fsu',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(mxyz(:,'spine_upper',:),4,[1.5]./(mxyz.sampleRate/2),'low'));
rb = mxyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = mxyz.com(rb);
mxyz.addMarker('hcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
mxyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[1.5]./(mxyz.sampleRate/2),'low'));
clear('hcom');

% traj NORM body vector projection
shft = 3;
tmar = {'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'};
tvec = [];
for m = 1:numel(tmar),
    tvec(:,m,:) = circshift(mxyz(:,tmar{m},[1,2]),-shft)-circshift(mxyz(:,tmar{m},[1,2]),shft);
end
nind = nniz(tvec);

% body vector
mvec = mxyz(:,'spine_upper',[1,2])-mxyz(:,'fsl',[1,2]);
umvec = bsxfun(@rdivide,bsxfun(@times,permute([1,-1],[1,3,2]),mvec(:,1,[2,1])),sqrt(sum(mvec.^2,3)));

unvec = [];
rotationAngles = deg2rad([0,90]);
for theta = rotationAngles,
    rotMat = repmat(permute([cos(theta),-sin(theta);sin(theta),cos(theta)],[3,1,2]),[size(mvec,1),1,1]);
    unvec(:,end+1,:) = bsxfun(@rdivide,multiprod(mvec,rotMat,[2,3],[2,3]),sqrt(sum(mvec.^2,3)));
end

walkFetRot = [];
for t = rotationAngles;
    for m = 1:numel(tmar),
        walkFetRot(nind,t==rotationAngles,m) = nunity(dot(tvec(nind,m,:),unvec(nind,t==rotationAngles,:),3));
    end
end

nz = nniz(mxyz);
wfet = mxyz.copy;
wfet.data= zeros([size(mxyz,1),size(walkFetRot,2)*size(walkFetRot,3)]);
wfet.data(nz,:) = [reshape(walkFetRot(nz,:),[],size(walkFetRot,2)*size(walkFetRot,3))];
wfs = wfet.segs([],round(mxyz.sampleRate/2));
wfs = MTADxyz('data',reshape(permute(sum(wfs.^2),[2,3,1]),size(wfs,2),[]),'sampleRate',mxyz.sampleRate);
wfs.data = circshift(wfs.data,round(mxyz.sampleRate/4),2);
wfs.data(isnan(wfs.data(:)))=0;
wfs.resample(xyz);

%% SS
sd = sqrt(sum((ss.data-circshift(ss.data,-1,2)).^2,3));
sn = sum(sd(:,2:end-1),2)./sd(:,end);
sv = xyz.copy;
sv.data = sn;


%% AV
sang = [circ_dist(fang(:,'spine_lower','pelvis_root',1),fang(:,'pelvis_root','spine_middle',1)),...
        circ_dist(fang(:,'pelvis_root','spine_middle',1),fang(:,'spine_middle','spine_upper',1)),...
        circ_dist(fang(:,'spine_middle','spine_upper',1),fang(:,'spine_upper','hcom',1))];
av = fang.copy;
av.data = abs(sum(sang,2)-circ_mean(sum(sang,2)));


fet.data =  [fet.data,...
             man.data,...                39
             fxyz(:,{'spine_lower',  ... 40
                    'pelvis_root',   ... 41
                    'spine_middle',  ... 42
                    'spine_upper',   ... 43
                    'hcom'},3),      ... 44
             fvelxy(:,{'spine_lower',... 45
                      'spine_upper', ... 46
                      'head_front',  ... 47
                      'bcom',        ... 48
                      'hcom',        ... 49
                      'acom'}),      ... 50 
             fvelz(:,{'spine_lower', ... 51
                     'spine_upper',  ... 52
                     'head_front',   ... 53
                     'bcom',         ... 54
                     'hcom',         ... 55
                     'acom'}),       ... 56
             zv.data,... 57 58 59 60 61 62 63 64 65 66
             sv.data,... 67
             av.data...  68
            ];


featureTitles = {};
featureDesc = {};



