function updateStc2NewVersion(stc);
ds = load(fullfile(stc.path,stc.filename));
stc.states = ds.states;
stc.save(1);