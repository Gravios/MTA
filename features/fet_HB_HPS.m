function [fet,featureTitles,featureDesc] = fet_HB_HPS(Trial,varargin)
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
xyz.filter('RectFilter');
vxy = xyz.vel({'spine_lower','hcom'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);

% PITCHES 
pch = fet_HB_pitch(Trial);
pch.map_to_reference_session(Trial,referenceTrial,referenceFeature);
pch.resample(newSampleRate);

% CONCATENATE features
fet.data = [circ_dist(pch(:,3),pch(:,1)),pch(:,1),vxy.data];
fet.data(~nniz(xyz),:)=0;

% SET feature title and descriptions 
featureTitles = {};
featureDesc = {};
if nargout>1,
    % 1. Head Pitch
    featureTitles(end+1) = {'Pitch HCHN-BMBU'};    
    featureDesc(end+1) = {['head pitch relative to xy plane subtracted by body pitch']};
    % 2. Body Pitch
    featureTitles(end+1) = {'Pitch BMBU'};    
    featureDesc(end+1) = {['body pitch relative to xy plane']};
    % 3. Body Speed
    featureTitles(end+1) = {'Speed xy BL'};
    featureDesc(end+1) = {['log10 speed of lower body marker']};
    % 4. Head Speed
    featureTitles(end+1) = {'Speed xy HC'};
    featureDesc(end+1) = {['log10 speed of head COM']};
    
end

% END MAIN -----------------------------------------------------------------------------------------