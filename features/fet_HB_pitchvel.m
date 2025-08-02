function [fet,featureTitles,featureDesc] = fet_HB_pitchvel(Trial,varargin)
% $$$ function [fet,featureTitles,featureDesc] = fet_head_pitch(Trial,varargin)
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%



% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate', Trial.xyz.sampleRate,                                          ...
                 'normalize'    , false,                                                         ...
                 'procOpts'     , {'SPLINE_SPINE_HEAD_EQI'});

[newSampleRate,normalize,procOpts] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'Head body pitch','fet_HB_pitchvel','b');                  

% XYZ preprocessed 
xyz = preproc_xyz(Trial,procOpts);
xyz.resample(newSampleRate);
%xyz.filter('ButFilter',3,30,'low');
xyz.filter('RectFilter');

% ANG Filtered Intermarker angles 
ang = create(MTADang,Trial,xyz);

% CAT feature
fet.data = [circ_dist(circshift(ang(:,'pelvis_root','spine_upper',2),-1),...
                      circshift(ang(:,'pelvis_root','spine_upper',2), 1)),...
            circ_dist(circshift(ang(:,'spine_upper','hcom'       ,2),-1),...
                      circshift(ang(:,'spine_upper','hcom'       ,2), 1)),...
            circ_dist(circshift(ang(:,'head_back',  'head_front', 2),-1),...
                      circshift(ang(:,'head_back',  'head_front', 2), 1))];            


fet.data(~nniz(xyz),:)=0;
featureTitles = {};
featureDesc = {};
if nargout>1,

    featureTitles(end+1) = {'Pitch Vel BPBU'};    
    featureDesc(end+1) = {['head pitch relative to xy plane']};
    featureTitles(end+1) = {'Pitch Vel BUHC'};    
    featureDesc(end+1) = {['head body pitch relative to xy plane']};    
    featureTitles(end+1) = {'Pitch Vel HBHF'};    
    featureDesc(end+1) = {['head pitch relative to xy plane']};
end

% END MAIN -----------------------------------------------------------------------------------------