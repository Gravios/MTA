function Stc = reduce_stc_to_loc(Stc);



Stc.states(Stc.gsi('x')) = [];

loc = Stc{'w+n'}.copy;
loc.label = 'loc';
loc.key   = 'x';

Stc.states{end} = loc;
