



tlist = 'jg05';
T = get_session_list(tlist);
for t = 1:numel(T)
    Trial = MTATrial.validate(T(t));
    pfs_3d_states(Trial);
end