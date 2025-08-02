sessionListTag  = 'hand_labeled';%'BHV_S4H5';
sessionList     = get_session_list(sessionListTag);
numSessions     = numel(sessionList);
sampleRate = repmat({119.881035},1,numSessions);

Trials = af(@(Trial) MTATrial.validate(Trial)  , sessionList);
Stc    = cf(@(Trial) Trial.load('stc')         , Trials);
StcNN  = cf(@(Trial) Trial.load('stc','NN0317'), Trials);
xyz    = cf(@(Trial) preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD_NO_TRB'), Trials);
