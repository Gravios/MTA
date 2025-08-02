function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_hvfl(Trial,varargin)
% function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_href_HXY(Trial,varargin)
% 
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%     filterCutoff:  numeric,  (4) - lowpass filter in Hz
%     theta:         numeric,  (0) - angle to rotate head coordinates
%     markers:       CellARY,  {'hcom','nose'} vector for projection

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate'  , Trial.xyz.sampleRate,                                        ...
                 'normalize'      , false,                                                       ...
                 'procOpts'       , {{'trb'}},                                                   ...
                 'filterCutoff'   , 4,                                                           ...
                 'rotation'       , Trial.meta.correction.headYaw,                               ...
                 'vector'         , {{'hcom','nose'}},                                           ...
                 'marker'         , 'hcom'                                                       ...
);
[newSampleRate,normalize,procOpts,filterCutoff,rotation,vector,marker] =                         ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

if isprop(Trial,'meta')
    if isfield(Trial.meta,'correction');
        theta = Trial.meta.correction.headYaw;
    end
end


% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'head referenced position and motion','fet_href_H','h');

% PREPROC xyz
% FILTER xyz
xyz = preproc_xyz(Trial,procOpts);
xyz.resample(newSampleRate);
xyz.data(~nniz(xyz),:,:) = 0;
xyz.filter('ButFilter',4,filterCutoff,'low');


% COMPUTE hcom projection onto head reference
fet.data = multiprod( sq(circshift(xyz(:,marker,[1,2]),-1)-circshift(xyz(:,marker,[1,2]),1)),   ...
                      transform_vector_to_rotation_matrix( xyz, vector, rotation),              ...
                      2,...
                      [2,3]);
fet.data = fet.data*(xyz.sampleRate/2)/10;

featureTitles = {'speedAP','speed_LAT'};
featureDesc = {['Anteroposterior velocity magnitude within the head frame of reference'],...
               ['Lateral velocity magnitude within the head frame of reference']};
if nargout>1,
end


%---------------------------------------------------------------------------------------------------
