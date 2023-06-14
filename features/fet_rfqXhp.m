function [fet,featureTitles,featureDesc] = fet_rfqXhp(Trial,varargin)
% function [fet,featureTitles,featureDesc] = fet_head_pitch(Trial,varargin)
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%     referenceTrial'  , 
%     referenceFeature', ''
%   

Trial = MTATrial.validate(Trial);

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate'   , Trial.xyz.sampleRate,                                       ...
                 'normalize'       , false,                                                      ...
                 'procOpts'        , {{}},                                                       ...
                 'referenceTrial'  , 'Ed05-20140529.ont.all',                                    ...
                 'referenceFeature', '',                                                         ...
                 'overwrite'   , false                                                       ...
);
[newSampleRate,normalize,procOpts,referenceTrial,referenceFeature,overwrite] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              Trial.xyz.sampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'respiration frequency and body pitch','fet_rfqXhp','r');

fet.updateFilename(Trial);

if overwrite || ~exist(fet.fpath,'file'),
% XYZ preprocessed 
    xyz = preproc_xyz(Trial,'trb',250);
% LOAD pitches 
    pch = fet_HB_pitch(Trial,250);
    rfq = fet_respiration_freq(Trial,250,65,'',false);
    if ~isempty(referenceTrial),
% MAP to reference trial
        pch.map_to_reference_session(Trial,referenceTrial,referenceFeature);
    end
% CONCATENATE features
    fet.data = [circ_dist(pch(:,3),pch(:,1)),[0;rfq.data;0]];
    fet.data(~nniz(xyz),:)=0;
    fet.save();
else
    load(fet,Trial);
end



fet.resample(newSampleRate);

featureTitles = {};
featureDesc = {};
if nargout>1,
    featureTitles(end+1) = {'Pitch HCHN-BMBU'};    
    featureDesc(end+1) = {['head pitch relative to xy plane']};
    featureTitles(end+1) = {'Head Height'};    
    featureDesc(end+1) = {['z-coordinate of head center']};
end

% END MAIN -----------------------------------------------------------------------------------------