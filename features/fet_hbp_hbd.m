function [fet,featureTitles,featureDesc] = fet_hbp_hbd(Trial,varargin)
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
                 'offset'          , [0,0]                                                       ...
);
[newSampleRate,normalize,procOpts,referenceTrial,referenceFeature,offset] = ...
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
              [],'TimeSeries',[],'Head body pitch and distance','fet_hbp_hbd','b');                  

% XYZ preprocessed 
% NONE
xyz = preproc_xyz(Trial,'trb');

% PITCHES 
pch = fet_HB_pitch(Trial);
if ~isempty(referenceTrial),
% MAP to reference trial
    pch.map_to_reference_session(Trial,referenceTrial,referenceFeature);
end
pch.resample(newSampleRate);
xyz.resample(newSampleRate);


blen = filter(copy(xyz),'ButFilter',4,2,'low');
blen.data = sqrt(sum(diff(blen(:,{'hcom','bcom'},[1,2]),1,2).^2,3));


% CONCATENATE features
fet.data = [circ_dist(pch(:,3),pch(:,1))+offset(1),blen.data+offset(2)];

fet.data(~nniz(xyz),:)=0;
featureTitles = {};
featureDesc = {};
if nargout>1,
    featureTitles(end+1) = {'Pitch HCHN-BMBU'};    
    featureDesc(end+1) = {['head pitch relative to xy plane']};
    featureTitles(end+1) = {'Pitch BMBU'};    
    featureDesc(end+1) = {['upper body pitch relative to xy plane']};
end

% END MAIN -----------------------------------------------------------------------------------------



