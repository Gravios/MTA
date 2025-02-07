function Stc = reduce_stc_to_loc(Stc,Trial);

if ischar(Stc)
    Stc = Trial.load('stc',Stc);
end

Stc.states(Stc.gsi('x')) = [];

loc = Stc{'w+n'}.copy;
loc.label = 'loc';
loc.key   = 'x';

Stc.states{end} = loc;

Stc.save(1)