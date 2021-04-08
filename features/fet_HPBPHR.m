function [fet,featureTitles,featureDesc] = fet_HPBPHR(Trial,varargin)
% function [fet,featureTitles,featureDesc] = fet_HPBPHR(Trial,varargin)
% {head-pitch body-pitch head-RHM}
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
                 'procOpts'        , {'trb'},                                                    ...
                 'referenceTrial'  , '',                                                         ...
                 'referenceFeature', '',                                                         ...
                 'overwriteFlag'   , true                                                        ...
);
[newSampleRate,normalize,procOpts,referenceTrial,referenceFeature,overwriteFlag] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'Head Pitch, Body pitch, head speed','fet_HPBPHS','h');                  

fet.updateFilename(Trial);
fet.updatePath(Trial.spath);

if overwriteFlag || ~exist(fullfile(fet.path,fet.filename),'file'),
% XYZ preprocessed 
    xyz = preproc_xyz(Trial,procOpts);
% LOAD pitches 
    pch = fet_HB_pitchB(Trial);
    if ~isempty(referenceTrial),
% MAP to reference trial
        pch.map_to_reference_session(Trial,referenceTrial,referenceFeature);
    end
    pch.resample(newSampleRate);
    [rhm,fs] = fet_rhm(Trial,[],'mtchglong','newSR',newSampleRate);
    rhm.data(rhm.data<1e-10) = 1e-10;
    rhm = mean(log10(rhm(:,fs>5 & fs<14)),2);
    xyz.resample(newSampleRate);
% CONCATENATE features
    fet.data = [pch.data,rhm];
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
    featureTitles(end+1) = {'Pitch BMBU'};    
    featureDesc(end+1) = {['upper body pitch relative to xy plane']};
    featureTitles(end+1) = {'Speed HCOM '};    
    featureDesc(end+1) = {['head speed in xy plane']};
end

% END MAIN -----------------------------------------------------------------------------------------g
