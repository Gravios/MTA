function [fet,featureTitles,featureDesc] = fet_hbav(Trial,varargin)
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
                 'offset'          , [0],                                                        ...
                 'lowPassCutOff'   , [2.4]                                                       ...                 
);
[newSampleRate,normalize,procOpts,xyz,referenceTrial,referenceFeature,offset,lowPassCutOff] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------





% MAIN ---------------------------------------------------------------------------------------------

fet = fet_hba(Trial,newSampleRate,normalize,procOpts,xyz,xyz,referenceTrial,referenceFeature,offset);
fet.type = 'TimeSeries';
fet.label = 'Head body angular velocity';
fet.name = 'fet_hbav';
fet.key = 'a';

fet.filter('ButFilter', 4, lowPassCutOff, 'low');
fet.data = (circshift(fet.data,-1)-circshift(fet.data,1)) * (fet.sampleRate/2);

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


