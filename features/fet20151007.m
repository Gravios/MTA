function feature = fet20151007(Trial,varargin)
%function feature = fet20151007(Trial,varargin)
%
%  varargin: [sampleRate] 
%     
%    sampleRate: numeric, the final sampleRate of the feature set
%
%    normalize:  logical, flag for optional z-transform
%
% TEST ARGS
%
%Trial = MTATrial('jg05-20120317');
%
% END TEST ARGS


%%
[sampleRate,normalize] = DefaultArgs(varargin,{10,false});

% Constant for now
REF_TRIAL = MTATrial('jg05-20120317');
REF_STATE = REF_TRIAL.stc{'a'};
RFET = fet_tsne(REF_TRIAL,sampleRate);


[feature,fett,fetd] = fet_tsne(Trial,sampleRate);

%NOTE this function only works with the feature set from fet_tsne
feature.normalize_to_reference(RFET);            


if normalize,
    % Get reference means and standard deviations for the reference features
    [~,Rmean,Rstd] = unity(RFET,[],[],[],[],REF_STATE);
    % Normalize the feautres using the z-transform relative to 
    feature.unity([],Rmean,Rstd);
end

