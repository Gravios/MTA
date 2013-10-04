function Pfs = gen_pfs(Trial)
Pfs = {};
for i = 1:length(Trial.Bhv.States),
    Pfs{i} = MTAPlaceField(Trial,[],{Trial.Bhv.States{i}.label},1);
end