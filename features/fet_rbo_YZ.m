function [fet,featureTitles,featureDesc] = fet_rbo_YZ(Trial,varargin)
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
defargs = struct('sampleRate'   , Trial.xyz.sampleRate,                                          ...
                 'normalize'       , false,                                                      ...
                 'procOpts'        , Trial.subject.label,                                        ...
                 'referenceTrial'  , '',                                                         ...
                 'referenceFeature', ''                                                          ...
);
[sampleRate,normalize,procOpts,referenceTrial,referenceFeature] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

procOpts = 'FS04_AC';

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              sampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'Head body pitch','fet_HB_pitch','b');                  


xyz = resample(Trial.load('subject',procOpts),sampleRate);

fet.data = [sq(xyz(:,'Head',[2,3]))];

fet.data(~nniz(xyz),:)=0;
featureTitles = {};
featureDesc = {};
if nargout>1,
    featureTitles(end+1) = {'position Y (mm)'};    
    featureDesc(end+1) = {['Position along the Y axis']};
    featureTitles(end+1) = {'Position Z (mm)'};    
    featureDesc(end+1) = {['Position along the Z axis']};
end


% END MAIN -----------------------------------------------------------------------------------------
