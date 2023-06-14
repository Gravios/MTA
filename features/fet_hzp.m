function [fet,featureTitles,featureDesc] = fet_hzp(Trial,varargin)
% function [fet,featureTitles,featureDesc] = fet_head_pitch(Trial,varargin)
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%     referenceTrial'  , 'Ed05-20140529.ont.all'
%     referenceFeature', ''
%   

Trial = MTATrial.validate(Trial);

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate'   , Trial.xyz.sampleRate,                                       ...
                 'normalize'       , false,                                                      ...
                 'procOpts'        , {{}},                                                       ...
                 'referenceTrial'  , '',                                                         ...
                 'referenceFeature', '',                                                         ...
                 'overwriteFlag'   , false                                                       ...
);
[newSampleRate,normalize,procOpts,referenceTrial,referenceFeature,overwriteFlag] = ...
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
              [],'TimeSeries',[],'Head height and body pitch','fet_hzp','h');

fet.updateFilename(Trial);

if overwriteFlag || ~exist(fet.fpath,'file'),
% XYZ preprocessed 
    xyz = preproc_xyz(Trial,'trb');
% LOAD pitches 
    pch = fet_HB_pitch(Trial);
    if ~isempty(referenceTrial),
% MAP to reference trial
        pch.map_to_reference_session(Trial,referenceTrial,referenceFeature);
    end
% CONCATENATE features
    fet.data = [circ_dist(pch(:,3),pch(:,1)),xyz(:,'hcom',3)];
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