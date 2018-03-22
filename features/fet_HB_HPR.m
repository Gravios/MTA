function [fet,featureTitles,featureDesc] = fet_HB_HPR(Trial,varargin)
% $$$ function [fet,featureTitles,featureDesc] = fet_head_pitch(Trial,varargin)
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
                 'procOpts'        , {{'trb'}},                                                  ...
                 'referenceTrial'  , 'Ed05-20140529.ont.all',                                    ...
                 'referenceFeature', ''                                                          ...
);
[newSampleRate,normalize,procOpts,referenceTrial,referenceFeature] = ...
    DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'Head pitch vs body speed','fet_HB_pitch','b');                  

% XYZ preprocessed 
xyz = preproc_xyz(Trial,procOpts);
[rhm,fs] = fet_rhm(Trial,[],'mtchglong',true);
rhm.data = mean(rhm(:,6<fs&fs<12),2);
rhm.resample(xyz);
rhm.data(rhm.data<1e-9) = 1e-10;
rhm.data = log10(rhm.data);

% PITCHES 
pch = fet_HB_pitch(Trial);
pch.map_to_reference_session(Trial,referenceTrial,referenceFeature);
pch.resample(newSampleRate);

% CONCATENATE features
fet.data = [circ_dist(pch(:,3),pch(:,1)),rhm.data];
fet.data(~nniz(xyz),:)=0;

% SET feature title and descriptions 
featureTitles = {};
featureDesc = {};
if nargout>1,
    % 1. Head Pitch
    featureTitles(end+1) = {'Pitch HCHN-BMBU'};    
    featureDesc(end+1) = {['head pitch relative to xy plane subtracted by body pitch']};
    % 2. RHM
    featureTitles(end+1) = {'RHMp'};
    featureDesc(end+1) = {['log10 6-12Hz Rhythmic Head Motion']};
    
end

% END MAIN -----------------------------------------------------------------------------------------