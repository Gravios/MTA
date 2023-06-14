function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_mis_HBREF(Trial,varargin)
% $$$ function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_mis_HBREF(Trial,varargin)
% $$$ defargs = struct('newSampleRate', 12,                       ...
% $$$                  'normalize'    , false,                    ...
% $$$                  'procOpts'     , {'SPLINE_SPINE_HEAD_EQD'});

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate', 12,                       ...
                 'normalize'    , false,                    ...
                 'procOpts'     , {{'SPLINE_SPINE_HEAD_EQD'}});

[newSampleRate,normalize,procOpts] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'req20160310_selected_features','fet_mis_HBREF','m');

% XYZ preprocessed 
[xyz,ss] = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD');
xyz.resample(newSampleRate);
ss.resample(xyz);

% XYZ filtered 
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');

hvfl = fet_href_HXY(Trial,newSampleRate,[],'trb');
bvfl = fet_bref_BXY(Trial,newSampleRate,[],'trb');

% FVELZ Filtered marker speeds in Z axis
fvelz = circshift(fxyz(:,{'hcom','bcom'},[3]),-1)-circshift(fxyz(:,{'hcom','bcom'},[3]),1);

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);

% SS Spine LENGTH
sd = sqrt(sum((ss.data-circshift(ss.data,-1,2)).^2,3));
sn = sum(sd(:,2:end-1),2)./sd(:,end);
sv = copy(Trial.xyz);
sv.data = sn;

% CAT feature
fet.data = [ bvfl.data,                      ... 1-2
             hvfl.data,                      ... 3-4
             fvelz,                          ... 5-6
             fang(:,'spine_lower','hcom',2), ... 7
             fang(:,'spine_middle','hcom',2),... 8
             fxyz(:,'pelvis_root',3),        ... 9
             fxyz(:,'spine_middle',3),       ... 10
             fxyz(:,'hcom',3),               ... 11
             sv.data];                         % 12
fet.data(~nniz(xyz),:)=0;
featureTitles = {};
featureDesc = {};

%---------------------------------------------------------------------------------------------------
 
             
