function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_mis(Trial,varargin)
% $$$ function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_mis(Trial,varargin)
% $$$ defargs = struct('newSampleRate', 12,                       ...
% $$$                  'normalize'    , false,                    ...
% $$$                  'procOpts'     , {'SPLINE_SPINE_HEAD_EQD'});


defargs = struct('newSampleRate', 12,                       ...
                 'normalize'    , false,                    ...
                 'procOpts'     , {'SPLINE_SPINE_HEAD_EQD'});

[newSampleRate,normalize,procOpts] = DefaultArgs(varargin,defargs,'--struct');


% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'req20160310_selected_features','fet_mis','m');                  

% XYZ preprocessed 
[xyz,ss] = preproc_xyz(Trial,newSampleRate,procOpts);
xyz.resample(newSampleRate);

% XYZ filtered 
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');

% FVELXY Filtered marker speeds in XY plane
fvelxy = fxyz.vel({'spine_lower','spine_upper','head_front','acom'},[1,2]);
%fvelxy = xyz.vel({'spine_lower','spine_upper','head_front','acom'},[1,2]);
%fvelxy.filter('ButFilter',3,2.5,'low');
fvelxy.data(fvelxy.data<0)=.1;
fvelxy.data = log10(fvelxy.data);

% FVELZ Filtered marker speeds in Z axis
fvelz = fxyz.vel({'head_front','acom'},[3]);
fvelz.filter('ButFilter',3,2.5,'low');
fvelz.data(fvelz.data<0)=.1;
fvelz.data = log10(fvelz.data);

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);

% SS 
sd = sqrt(sum((ss.data-circshift(ss.data,-1,2)).^2,3));
sn = sum(sd(:,2:end-1),2)./sd(:,end);
sv = Trial.xyz.copy;
sv.data = sn;
sv.resample(fxyz);

% AV 
sang = [circ_dist(fang(:,1,2,1),fang(:,2,3,1)),...
        circ_dist(fang(:,2,3,1),fang(:,3,4,1)),...
        circ_dist(fang(:,3,4,1),fang(:,4,'hcom',1))];
av = fang.copy;
av.data = abs(sum(sang,2)-circ_mean(sum(sang,2)));

% PPC feature
try
    man = Trial.load('fet','lsppc');
catch err
    gen_fet_lsppc(Trial);    
    man = Trial.load('fet','lsppc');
end
man.filter('ButFilter',3,2,'low');
man.resample(fxyz);

% CAT feature
fet.data = [ fang(:,'spine_lower','spine_middle',2),...
             fang(:,'spine_lower','hcom',2),        ...
             fang(:,'pelvis_root','spine_upper',2), ...
             fang(:,'spine_middle','spine_upper',2),...
             fang(:,'spine_middle','hcom',2),       ...
             man.data,                              ...
             fxyz(:,'spine_lower',3),               ...
             fxyz(:,'pelvis_root',3),               ...
             fxyz(:,'spine_middle',3),              ...             
             fxyz(:,'spine_upper',3),               ...             
             fvelxy.data,                           ...
             fvelz.data,                            ...
             sv.data,                               ...
             av.data                                ...
];
