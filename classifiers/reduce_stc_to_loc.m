function Stc = reduce_stc_to_loc(Stc);

loc = Stc{'w+n'};
loc.label = 'loc';
loc.key   = 'x';

Stc.states{end+1} = loc;

Stc.save(1);