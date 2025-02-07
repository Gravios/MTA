function [fet,featureTitles,featureDesc] = fet_HB_pitch(Trial,varargin)
% $$$ function [fet,featureTitles,featureDesc] = fet_head_pitch(Trial,varargin)
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%

Trial = MTATrial.validate(Trial);

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate', Trial.xyz.sampleRate,                                          ...
                 'normalize'    , false,                                                         ...
                 'procOpts'     , {{'SPLINE_SPINE_HEAD_EQD'}},                                   ...
                 'referenceTrial'  , 'Ed05-20140529.ont.all',                                    ...
                 'referenceFeature', ''                                                          ...                 
);
[newSampleRate,normalize,procOpts,referenceTrial,referenceFeature] =                             ...
    DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              Trial.xyz.sampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'Head body pitch','fet_HB_pitch','b');                  

% XYZ preprocessed 
try
    xyz = preproc_xyz(Trial,procOpts);
catch err
    disp(err);
    xyz = preproc_xyz(Trial);
end




% ANG Filtered Intermarker angles 
ang = create(MTADang,Trial,xyz);


% REMOVE THIS AT SOME POINT... maybe.. who am I kidding this will
% be here forever.
mnose = 'nose';
if ~ang.model.gmi(mnose)
    mnose = 'head_nose';
end


% CAT feature
%fet.data = [ang(:,'spine_middle','spine_upper',2),...
fet.data = [ang(:,'pelvis_root','spine_upper',2),...
            ang(:,'spine_upper','hcom',2),...
            ang(:,'hcom',mnose,2)];


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


if ~isempty(referenceTrial),
% MAP to reference trial
    fet.map_to_reference_session(Trial,referenceTrial,referenceFeature);
end

fet.resample(newSampleRate);


% END MAIN -----------------------------------------------------------------------------------------