function save(Trial)
stcmode = [];
if ~isempty(Trial.stc)
    stcmode = Trial.stc.mode;
end
trialName = Trial.trialName;
sync = Trial.sync;
save(fullfile(Trial.spath,[Trial.filebase '.trl.mat']),'trialName','sync','stcmode');            
end
