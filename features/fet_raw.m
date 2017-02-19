function [fet,featureTitles,featureDesc] = fet_raw(Trial,varargin)
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


% $$$ if ~isempty(RefTrial),
% $$$     markers = {'spine_lower','pelvis_root','spine_middle','spine_upper',...
% $$$                'head_back',  'head_left'  ,'head_front',  'head_right',...
% $$$                'bcom',       'hcom',       'acom'};
% $$$     mfet = fet.copy;
% $$$     mfet.data = xyz(:,markers,3);
% $$$     mfet.map_to_reference_session(Trial,RefTrial);
% $$$     xyz.data(:,xyz.model.gmi(markers),3) = mfet.data;
% $$$ end              
              
xyz.resample(newSampleRate);

% XYZ filtered 
fxyz = xyz.copy;

% FVELXY Filtered marker speeds in XY plane
fvelxy = xyz.vel([],[1,2]);
fvelxy.filter('ButFilter',3,2.5,'low');
fvelxy.data(fvelxy.data<0)=.1;

facxy = fvelxy.copy;
facxy.data = [diff(facxy.data);zeros([1,fvelxy.size(2)])];

fvelxy.data = log10(fvelxy.data);


% FVELZ Filtered marker speeds in Z axis
fvelz = fxyz.vel([],[3]);
fvelz.data(fvelz.data<0)=.1;

facz = fvelz.copy;
facz.data = [diff(facz.data);zeros([1,fvelz.size(2)])];

fvelz.data = log10(fvelz.data);

% ANG InterMarker Spherical Coordinates
ang = create(MTADang,Trial,xyz);

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);


%% Derivatives of features

% DISTANCE between segments 
d = 3;
cang = ang.copy;
cang.data = [fang(:,'spine_lower','pelvis_root',d),...          1
             fang(:,'spine_lower','spine_middle',d),...         2
             fang(:,'spine_lower','bcom',d),...                 3
             fang(:,'spine_lower','hcom',d),...                 4
             fang(:,'pelvis_root','spine_middle',d),...         5
             fang(:,'pelvis_root','spine_upper',d),...          6
             fang(:,'spine_middle','spine_upper',d),...         7
             fang(:,'spine_middle','hcom',d),...                8
             fang(:,'spine_upper','hcom',d),...                 9
             fang(:,'spine_upper','hcom',d),...                 10
             fang(:,'bcom','hcom',d),...                        11
             fang(:,'bcom','acom',d)];                         %12
fet.data =  [fet.data,cang.data];

% PITCH between segments 
d = 2;
cang = ang.copy;
cang.data = [fang(:,'spine_lower','pelvis_root',2),...          13
             fang(:,'spine_lower','spine_middle',2),...         14
             fang(:,'spine_lower','bcom',2),...                 15
             fang(:,'spine_lower','hcom',2),...                 16
             fang(:,'pelvis_root','spine_middle',2),...         17
             fang(:,'pelvis_root','spine_upper',2),...          18
             fang(:,'spine_middle','spine_upper',2),...         19
             fang(:,'spine_middle','hcom',2),...                20
             fang(:,'spine_upper','hcom',2),...                 21
             fang(:,'spine_upper','hcom',2),...                 22
             fang(:,'bcom','hcom',2),...                        23
             fang(:,'bcom','acom',2)];                         %24
fet.data =  [fet.data,cang.data];

% DRV of pitch
cang = fang.copy;
cang.data = [ang(:,'spine_lower','pelvis_root',2),...          25
             ang(:,'spine_lower','spine_middle',2),...         26
             ang(:,'spine_lower','bcom',2),...                 27
             ang(:,'spine_lower','hcom',2),...                 28
             ang(:,'pelvis_root','spine_middle',2),...         29
             ang(:,'pelvis_root','spine_upper',2),...          30
             ang(:,'spine_middle','spine_upper',2),...         31
             ang(:,'spine_middle','hcom',2),...                32
             ang(:,'spine_upper','hcom',2),...                 33
             ang(:,'bcom','hcom',2),...                        34
             ang(:,'bcom','acom',2),...                        35
             ang(:,'hcom','acom',2)];                         %36
%cang.filter('ButFilter',3,5,'low');
cang.data = diff(cang.data);
% Add to feature matrix
cang.resample(xyz);
fet.data =  [fet.data,cang.data];


%% DRV of Yaw
cang = ang.copy;
cang.data = [ang(:,'spine_lower','pelvis_root',1),...          37
             ang(:,'spine_lower','spine_middle',1),...         38
             ang(:,'spine_lower','bcom',1),...                 39
             ang(:,'spine_lower','hcom',1),...                 40
             ang(:,'pelvis_root','spine_middle',1),...         41
             ang(:,'pelvis_root','spine_upper',1),...          42
             ang(:,'spine_middle','spine_upper',1),...         43
             ang(:,'spine_middle','hcom',1),...                44
             ang(:,'spine_upper','hcom',1),...                 45
             ang(:,'spine_upper','hcom',1),...                 46
             ang(:,'head_back','hcom',1),...                   47
             ang(:,'bcom','hcom',1),...                        48
             ang(:,'bcom','acom',1),...                        49
             ang(:,'hcom','acom',1)];                         %50
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
%mxyz.filter('ButFilter',3,20);
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
             man.data,... 51
             fxyz(:,{'spine_lower',  ... 52 Height
                     'pelvis_root',  ... 53
                     'spine_middle', ... 54
                     'spine_upper',  ... 55
                     'hcom'},3),     ... 56
             ...
             fvelxy(:,{'spine_lower',... 57 xy speed
                       'spine_upper',... 58
                       'head_front', ... 59 
                       'bcom',       ... 60 
                       'hcom',       ... 61 
                       'acom'}),     ... 62
             ...
             fvelz(:,{'spine_lower', ... 63 z speed
                      'spine_upper', ... 64
                      'head_front',  ... 65
                      'bcom',        ... 66
                      'hcom',        ... 67
                      'acom'}),      ... 68
             [zeros(1,6);sq(diff(fxyz(:,...
                     {'spine_lower', ... 69 z vel
                      'spine_upper', ... 70
                      'head_front',  ... 71
                      'bcom',        ... 72
                      'hcom',        ... 73
                      'acom'},3)))], ... 74
             zv.data,...  75
             sv.data,...  76
             av.data...   77
            ];


featureTitles = {};
featureDesc = {};



