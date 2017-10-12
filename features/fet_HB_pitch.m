function [fet,featureTitles,featureDesc] = fet_HB_pitch(Trial,varargin)
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
              [],'TimeSeries',[],'Head body pitch','fet_HB_pitch','b');                  

% XYZ preprocessed 
xyz = preproc_xyz(Trial,procOpts);
xyz.resample(newSampleRate);
%xyz.filter('RectFilter');

% ANG Filtered Intermarker angles 
ang = create(MTADang,Trial,xyz);

% CAT feature
%fet.data = [ang(:,'pelvis_root','spine_upper',2),ang(:,'spine_upper','hcom',2)];
fet.data = [ang(:,'pelvis_root','spine_upper',2),...
            ang(:,'spine_upper','hcom',2),...
            ang(:,'head_back','head_front',2)];

fet.data(~nniz(xyz),:)=0;
featureTitles = {};
featureDesc = {};
if nargout>1,

    featureTitles(end+1) = {'Pitch BPBU'};    
    featureDesc(end+1) = {['head pitch relative to xy plane']};
    featureTitles(end+1) = {'Pitch BUHC'};    
    featureDesc(end+1) = {['head body pitch relative to xy plane']};
    featureTitles(end+1) = {'Pitch HBHF'};    
    featureDesc(end+1) = {['head pitch relative to xy plane']};
end

% END MAIN -----------------------------------------------------------------------------------------