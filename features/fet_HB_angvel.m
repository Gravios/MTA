function [fet,featureTitles,featureDesc] = fet_HB_angvel(Trial,varargin)
% $$$ function [fet,featureTitles,featureDesc] = fet_head_pitch(Trial,varargin)
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%



% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate', Trial.xyz.sampleRate,                                          ...
                 'normalize'    , false,                                                         ...
                 'procOpts'     , {'SPLINE_SPINE_HEAD_EQI'},                                     ...
                 'filtFreq'     , 30                                                             ...
);
[newSampleRate,normalize,procOpts,filtFreq] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'Head body pitch','fet_HB_angvel','b');                  

% XYZ preprocessed 
xyz = preproc_xyz(Trial,procOpts);
xyz.resample(newSampleRate);
%xyz.filter('ButFilter',3,filtFreq,'low');
xyz.filter('RectFilter');

% ANG Filtered Intermarker angles 
ang = create(MTADang,Trial,xyz);

% CAT feature
fet.data = [circ_dist(circshift(ang(:,'pelvis_root','spine_upper',1),-1),...
                      circshift(ang(:,'pelvis_root','spine_upper',1), 1)),...
            circ_dist(circshift(ang(:,'spine_upper','hcom',1)       ,-1),...
                      circshift(ang(:,'spine_upper','hcom',1)       , 1)),...
            circ_dist(circshift(ang(:,'head_back','head_front',1)       ,-1),...
                      circshift(ang(:,'head_back','head_front',1)       , 1))];            

fet.data(~nniz(xyz),:)=0;
featureTitles = {};
featureDesc = {};
if nargout>1,
    featureTitles(end+1) = {'Yaw AngVel BPBU'};    
    featureDesc(end+1) = {['body angular velocity within xy plane']};
    featureTitles(end+1) = {'Yaw AngVel BUHC'};    
    featureDesc(end+1) = {['head body angular velocity within xy plane']};
    featureTitles(end+1) = {'Yaw AngVel HBHF'};    
    featureDesc(end+1) = {['head angular velocity within xy plane']};
end

% END MAIN -----------------------------------------------------------------------------------------