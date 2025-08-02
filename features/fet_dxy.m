function [fet,featureTitles,featureDesc] = fet_dxy(Trial,varargin)
% function [fet,featureTitles,featureDesc] = fet_xy(Trial,varargin)
%
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%     referenceTrial'  , 'Ed05-20140529.ont.all'
%     referenceFeature', ''
%

Trial = MTATrial.validate(Trial);

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('sampleRate'   , Trial.xyz.sampleRate,                                       ...
                 'normalize'       , false,                                                      ...
                 'procOpts'        , 'trb',                                                  ...
                 'referenceTrial'  , '',                                                         ...
                 'referenceFeature', ''                                                          ...
);
[sampleRate,normalize,procOpts,referenceTrial,referenceFeature] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              sampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'Head XY position','fet_xy','p');



xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
fet.data = [sq(xyz(:,'hcom',[1,2]))];

hxy = sq(xyz(:,'head_front',[1,2]))-fet.data;
hxy = hxy./sqrt(sum(hxy.^2,2));
fet.data =  [atan2(hxy(:,1),hxy(:,2)),fet.data];


fet.data(~nniz(xyz),:)=0;
featureTitles = {};
featureDesc = {};
if nargout>1,
    featureTitles(end+1) = {'Head Yaw  (rad)'};    
    featureDesc(end+1) = {['Head direction']};
    featureTitles(end+1) = {'position x (mm)'};    
    featureDesc(end+1) = {['Position along the X axis']};
    featureTitles(end+1) = {'Position y (mm)'};    
    featureDesc(end+1) = {['Position along the Y axis']};
end


% END MAIN -----------------------------------------------------------------------------------------
