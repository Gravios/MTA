function [fet,featureTitles,featureDesc] = fet_hba(Trial,varargin)
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
                 'xyz'             , [],                                                         ...
                 'referenceTrial'  , '',                                                         ...
                 'referenceFeature', '',                                                         ...
                 'offset'          , [0]                                                         ...
);
[newSampleRate,normalize,procOpts,xyz,referenceTrial,referenceFeature,offset] = ...
    DefaultArgs(varargin,defargs,'--struct');
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
              [],'TimeSeries',[],'Head body pitch','fet_HB_pitch','b');                  

hbang = copy(xyz);
xycoor = cat(2,...
             hbang(:,'spine_upper',[1,2])-hbang(:,'bcom',[1,2]),...
             hbang(:,'nose',[1,2])-hbang(:,'hcom',[1,2]));
hbang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
hbang.data = -(circ_dist(hbang.data(:,2),hbang.data(:,1))+Trial.meta.correction.headBody);

% CONCATENATE features
fet.data = [hbang.data+offset];

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


