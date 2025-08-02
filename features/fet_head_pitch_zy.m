function [fet,featureTitles,featureDesc] = fet_head_pitch(Trial,varargin)
% $$$ function [fet,featureTitles,featureDesc] = fet_head_pitch(Trial,varargin)
% $$$ defargs = struct('newSampleRate', 12,                       ...
% $$$                  'normalize'    , false);


defargs = struct('newSampleRate', Trial.xyz.sampleRate,    ...
                 'normalize'    , false                    ...
                 );

[newSampleRate,normalize] = DefaultArgs(varargin,defargs,'--struct');


% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'req20160310_selected_features','fet_head_pitch','h');                  

% XYZ preprocessed 
xyz = Trial.load('xyz');
xyz.resample(newSampleRate);

% ANG Filtered Intermarker angles 
ang = create(MTADang,Trial,xyz);

% CAT feature
fet.data = ang(:,'BR','FR',2);

fet.data(~nniz(xyz),:)=0;
featureTitles = {};
featureDesc = {};
if nargout>1,

    featureTitles(end+1) = {'Pitch HBHF'};    
    featureDesc(end+1) = {['head pitch relative to xy plane']};
end