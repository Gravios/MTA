function [fet,featureTitles,featureDesc] = fet_hma(Trial,varargin)
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
                 'procOpts'        , 'trb',                                                      ...
                 'xyz'             , []                                                          ...
);
[newSampleRate,normalize,procOpts,xyz] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

% XYZ preprocessed 
if isempty(xyz)
    xyz = preproc_xyz(Trial,procOpts,newSampleRate);
else
    newSampleRate = xyz.sampleRate;
end

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'Head Maze Angle','hma','b');                  

fet.data  = zeros([size(xyz,1),1]);
headAngle = sq(xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2]));
headAngle = atan2(headAngle(:,2),headAngle(:,1));
mazeAngle = sq(xyz(:,'hcom',[1,2]));
mazeAngle = atan2(mazeAngle(:,2),mazeAngle(:,1));
fet.data  = circ_dist(circ_dist( headAngle, mazeAngle),Trial.meta.correction.headYaw);


% CONCATENATE features

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


