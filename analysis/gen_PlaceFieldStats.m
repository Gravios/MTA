function gen_PlaceFieldStats(Trial,
[states,stc_mode,pfs_profile,

states = {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};
stc_mode = 'auto_wbhr';

Trial = MTATrial(Trial,trialName,mazeName);
Trial.stc.updateMode(stc_mode);Trial.stc.load;
Trial.load('nq');


