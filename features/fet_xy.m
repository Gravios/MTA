function [fet,featureTitles,featureDesc] = fet_xy(Trial,varargin)
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
              [],'TimeSeries',[],'Head body pitch','fet_HB_pitch','b');                  



xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);

fet.data = [sq(xyz(:,'hcom',[1,2]))];

fet.data(~nniz(xyz),:)=0;
featureTitles = {};
featureDesc = {};
if nargout>1,
    featureTitles(end+1) = {'position x (mm)'};    
    featureDesc(end+1) = {['Position along the X axis']};
    featureTitles(end+1) = {'Position y (mm)'};    
    featureDesc(end+1) = {['Position along the Y axis']};
    featureTitles(end+1) = {'Pitch HCHN-BMBU'};    
    featureDesc(end+1) = {['head pitch relative to xy plane']};    
    featureTitles(end+1) = {'Pitch BMBU'};    
    featureDesc(end+1) = {['upper body pitch relative to xy plane']};
end


% END MAIN -----------------------------------------------------------------------------------------
