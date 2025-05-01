


rind = [pfstats(:,1).patchMFR]>2 | [pfstats(:,2).patchMFR]>2;
rind = rind(1,:,1)';


ppos = nan([size(rind,1),2]);
for ind = find(rind)'
    ppos(ind,:) =sq(pfstats(ind,2).patchCOM(1,1,1,:));
end

figure,
hist(sqrt(sum(ppos.^2,2)))


sum(sqrt(sum(ppos.^2,2))<300)


ds = load('/storage/gravio/data/project/general/jg05-20120312/jg05-20120312.cof.all.trl.mat');

Units = cf(@(T)  T.spk.get_unit_set(T,'placecells'),  Trials); 

% LOAD rate maps -----------------------------------------------------------------------------------
% LOAD place restricted behavior fields
tids = [3:7,17:24,27];
tids = [1,2,8:16,25:30];
for tid = tids
% $$$     s = MTASession.validate(Trials{tid}.filebase);
% $$$     s.spk.create(s)
% $$$     s.save();
% $$$     Trials{tid} = MTATrial.validate(Trials{tid}.filebase);
% $$$     NeuronQuality(Trials{tid}, 'overwrite',true);
    bfs = compute_bhv_ratemaps(Trials{tid},Units{tid},'overwrite',true,'purge',true);
end


for tid = tids
    pfst = pfs_2d_theta(Trials{tid},'overwrite',false);
    pfst.purge_savefile();
    pfst = pfs_2d_theta(Trials{tid},Units{tid},'overwrite',true);
end

for tid = tids
    pfss = pfs_2d_states(Trials{tid},'overwrite',false);
    for sts = 1:numel(pfss)
        pfss{sts}.purge_savefile();
    end
    pfss = pfs_2d_states(Trials{tid},Units{tid},'overwrite',true);
end

for tid = tids
    pfsr = pfs_2d_states( ...
        Trials{tid},...
        [], ...
        [],...
        {'rear&theta','hbhv&theta','lbhv&theta'}, ...
        'overwrite',false);
    for sts = 1:numel(pfsr)
        pfsr{sts}.purge_savefile();
    end
    pfsr = pfs_2d_states( ...
        Trials{tid},...
        Units{tid},...
        [],...
        {'rear&theta','hbhv&theta','lbhv&theta'},...
        'overwrite',true);
end

    


bfs      = cf(@(t,u)  compute_bhv_ratemaps(t,u,'overwrite',true),          Trials, Units);
bfsEx    = bfs{tind};
bfsShuff = cf(@(t,u)  compute_bhv_ratemaps_shuffled(t,u), Trials, Units);
% LOAD all place fields
pfst = cf(@(t,u)  pfs_2d_theta(t,u),    Trials,Units);
pfss = cf(@(t,u)  pfs_2d_states(t,u),   Trials,Units);
pfsr = cf(@(t,u)  pfs_2d_states(t,u,'',{'rear&theta','hbhv&theta','lbhv&theta'}),   Trials,Units);
pfsa = cf(@(s,t)  cat(2,{t},s),       pfss,pfst);
% --------------------------------------------------------------------------------------------------

