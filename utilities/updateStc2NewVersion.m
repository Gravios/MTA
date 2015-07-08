function stc = updateStc2NewVersion(stc);
ds = load(fullfile(stc.path,stc.filename));
stc.states = ds.states;
if nargout==0,
    stc.save(1);
end