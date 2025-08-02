function [fet,featureTitles,featureDesc] = fet_all(Trial,varargin)
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


if ~isempty(RefTrial),
    markers = {'spine_lower','pelvis_root','spine_middle','spine_upper',...
               'head_back',  'head_left'  ,'head_front',  'head_right',...
               'bcom',       'hcom',       'acom'};
    mfet = fet.copy;
    mfet.data = xyz(:,markers,3);
    mfet.map_to_reference_session(Trial,RefTrial);
    xyz.data(:,xyz.model.gmi(markers),3) = mfet.data;
end              
              
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
cang.filter('ButFilter',3,5,'low');
cang.data = diff(cang.data);
% Add to feature matrix
cang.resample(xyz);
fet.data =  [fet.data,cang.data];



%% DRV of Yaw
cang = ang.copy;
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
cang.filter('ButFilter',3,5,'low');
cang.data = circ_dist(cang.data,circshift(cang.data,1)).*cang.sampleRate;
cang.data = circshift(permute(sum(cang.segs(1:cang.size(1),newSampleRate,0).^2),...
                              [2,3,1]),...
                      [-round(newSampleRate/2),0]);
cang.resample(xyz);
cang.data(cang.data<1e-6) = 1e-6;
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
mxyz.filter('ButFilter',3,20);
tsh = round(.15*mxyz.sampleRate);
afet = mxyz.copy;
afet.data = circshift(mxyz.data,-tsh)-circshift(mxyz.data,tsh);
afet.data = reshape(afet.data,[],3);
afet.data = permute(bsxfun(@dot,permute(reshape(repmat(mxyz(:,'spine_upper',:)-mxyz(:,'spine_lower',:),[1,mxyz.size(2),1]),[],3),[2,1]),permute(afet.data,[2,1])),[2,1]);
afet.data = reshape(afet.data,[],mxyz.size(2));
afet.resample(fxyz);
zv = afet.copy;
zv.data = log10(abs(afet(:,1))).*sign(afet(:,1));

%% SS
sd = sqrt(sum((ss.data-circshift(ss.data,-1,2)).^2,3));
sn = sum(sd(:,2:end-1),2)./sd(:,end);
sv = Trial.xyz.copy;
sv.data = sn;
sv.resample(fxyz);


%% AV
sang = [circ_dist(fang(:,'spine_lower','pelvis_root',1),fang(:,'pelvis_root','spine_middle',1)),...
        circ_dist(fang(:,'pelvis_root','spine_middle',1),fang(:,'spine_middle','spine_upper',1)),...
        circ_dist(fang(:,'spine_middle','spine_upper',1),fang(:,'spine_upper','hcom',1))];
av = fang.copy;
av.data = abs(sum(sang,2)-circ_mean(sum(sang,2)));


fet.data =  [fet.data,...
             man.data,...        39
             fxyz(:,[1:5],3),... 40,41,42,43,44
             fvelxy(:,{'spine_lower','spine_upper','head_front','bcom','hcom','acom'}),... 45 46 47 48 49 50 
             fvelz(:,{'spine_lower','spine_upper','head_front','bcom','hcom','acom'}),...  51 52 53 54 55 56
             zv.data,... 57
             sv.data,... 58
             av.data...  59
            ];


featureTitles = {};
featureDesc = {};



