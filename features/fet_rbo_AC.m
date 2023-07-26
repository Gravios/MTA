function [fet,featureTitles,featureDesc] = fet_rbo_AC(Trial,varargin)
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
defargs = struct('sampleRate'   , Trial.subject.sampleRate,                                      ...
                 'normalize'       , false,                                                      ...
                 'procOpts'        , [Trial.subject.label,'_AC'],                                ...
                 'referenceTrial'  , '',                                                         ...
                 'referenceFeature', ''                                                          ...
);
[sampleRate,normalize,procOpts,referenceTrial,referenceFeature] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

%procOpts = 'FS03_AC';
%procOpts = 'FS04_RC';

% INIT Feature
fet = MTADfet(Trial.spath,                                 ... path
              [],                                          ... filename
              [],                                          ... data
              sampleRate,                                  ... samplerate
              Trial.subject.sync.copy,                     ... sync
              Trial.subject.sync.data(1),                  ... origin
              'TimeSeries',                                ... type
              [],                                          ... ext
              'Head Rigid Body Object in Arena Coordinaes',... Name
              'fet_rbo_AC',                                ... Label
              'b');                                          % Key

rbo = resample(Trial.load('subject',procOpts),sampleRate);

fet.data = [sq(rbo(:,'Head',[1,2]))];

fet.data(~nniz(rbo),:)=0;
featureTitles = {};
featureDesc = {};
if nargout>1,
    featureTitles(end+1) = {'position x (mm)'};    
    featureDesc(end+1) = {['Position along the X axis']};
    featureTitles(end+1) = {'Position y (mm)'};    
    featureDesc(end+1) = {['Position along the Y axis']};
end


% END MAIN -----------------------------------------------------------------------------------------
