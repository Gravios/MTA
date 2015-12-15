function feature = fet20151007(Trial,varargin)
%function feature = fet20151007(Trial,varargin)
%
%  varargin: [sampleRate] 
%     
%    sampleRate: numeric, the final sampleRate of the feature set
%
% TEST ARGS
%
%Trial = MTATrial('jg05-20120317');
%
% END TEST ARGS


%%
[sampleRate] = DefaultArgs(varargin,{10});

% Constant for now
REF_TRIAL = MTATrial('jg05-20120317');
REF_STATE = REF_TRIAL.stc{'a'};
RFET = fet_tsne(REF_TRIAL,sampleRate);


[feature,fett,fetd] = fet_tsne(Trial,sampleRate);

% Get reference means and standard deviations for the reference features
[~,Rmean,Rstd] = RFET.unity([],[],[],[],REF_STATE);

%NOTE this function only works with the feature set from fet_tsne
feature.normalize_to_reference(RFET);            

% Normalize the feautres using the z-transform relative to 
feature.unity(RFET,[],Rmean,Rstd);
