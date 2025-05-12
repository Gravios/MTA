


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
tids = [28:30];
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
        Units{tid}, ...
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

    


bfs      = cf(@(t,u)  compute_bhv_ratemaps(t,u),          Trials, Units);
bfsEx    = bfs{tind};
bfsShuff = cf(@(t,u)  compute_bhv_ratemaps_shuffled(t,u), Trials, Units);
% LOAD all place fields
pfst = cf(@(t,u)  pfs_2d_theta(t,u),    Trials,Units);
pfss = cf(@(t,u)  pfs_2d_states(t,u),   Trials,Units);
pfsr = cf(@(t,u)  pfs_2d_states(t,u,'',{'rear&theta','hbhv&theta','lbhv&theta'}),   Trials,Units);
pfsa = cf(@(s,t)  cat(2,{t},s),       pfss,pfst);
% --------------------------------------------------------------------------------------------------



clear('pairs');

pairs = repmat(struct('tid',0,'clm',[0,8],'el',0,'unt',[]),[9,1]);
pairs(1).tid = 17
pairs(1).clm = [0,8];
pairs(1).elc = 6;  pairs(1).unt = [];
pairs(1).elc = 8;  pairs(1).unt = [74, 57, 150, 124];
pairs(1).elc = 8;  pairs(1).unt = [74, 57, 150, 124];
pairs(2).elc = 8;  pairs(2).unt = [62, 50, 141, 105]; 
pairs(3).elc = 8;  pairs(3).unt = [89, 66, 144, 145]; % check 145
pairs(4).elc = 8;  pairs(4).unt = [80, 75, 142, 102]; 
pairs(5).elc = 8;  pairs(5).unt = [59, 49, 133, 128]; 
pairs(6).elc = 8;  pairs(6).unt = [61, 60, 134, 115]; 
pairs(7).elc = 8;  pairs(7).unt = [69, 56, 139, 114];
pairs(8).elc = 8;  pairs(8).unt = [70, 74, 149, 103];

tind = 17;
figure
elc = 6
for unitpair = cmtch{tind}{elc}'
    clf();
    subplot(141);
        plot(pft{tind},unitpair(1),1,'text');
        title(num2str(unitpair(1)));
    subplot(142);
        plot(pft{tind+1},unitpair(2),1,'text');
        title(num2str(unitpair(2)));
    % find next unit
    uind = find(ismember(unitpair(2), cmtch{tind+1}{elc}));
    if ~isempty(uind)
        u3 = cmtch{tind+1}{elc}(uind)
        subplot(143);
            plot(pft{tind+2}, u3, 1, 'text');
            title(num2str(u3));
        uind = find(ismember(u3, cmtch{tind+2}{elc}));
        if ~isempty(uind)
            u4 = cmtch{tind+3}{elc}(uind)
            subplot(144);
                plot(pft{tind+3}, u4, 1, 'text');
                title(num2str(u4));
        end
    end
    waitforbuttonpress();
end



clear('pairs');
pairs(1).tid = 17
pairs(1).clm = [0,8];
pairs(1).chn = 8;  pairs(1).unt = [74, 57, 150, 124];
pairs = repmat(pairs,[9,1]);
pairs(2).chn = 8;  pairs(2).unt = [62, 50, 141, 105]; 
pairs(3).chn = 8;  pairs(3).unt = [89, 66, 144, 145]; % check 145
pairs(4).chn = 8;  pairs(4).unt = [80, 75, 142, 102]; 
pairs(5).chn = 8;  pairs(5).unt = [59, 49, 133, 128]; 
pairs(6).chn = 8;  pairs(6).unt = [61, 60, 134, 115]; 
pairs(7).chn = 8;  pairs(7).unt = [69, 56, 139, 114];
pairs(8).chn = 8;  pairs(8).unt = [70, 74, 149, 103];


clear('pairs');
pairs(1).tid = 21
pairs(1).clm = [0,6];
pairs(1).chn = 6;  pairs(1).unt = [6, 13, 29];
pairs = repmat(pairs,[9,1]);
pairs(2).chn = 7;  pairs(2).unt = [24, 42, 51]; %check unit 24
pairs(3).chn = 7;  pairs(3).unt = [27, 41, 50]; 
pairs(4).chn = 7;  pairs(4).unt = [32, 38, 48]; %check unit 38
pairs(5).chn = 7;  pairs(5).unt = [33, 30, 54]; %check unit 33
pairs(6).chn = 8;  pairs(6).unt = [61, 48, 72]; 
pairs(7).chn = 8;  pairs(7).unt = [63, 65, 69]; %check unit 65
pairs(8).chn = 8;  pairs(8).unt = [73, 51, 67];
pairs(9).chn = 8;  pairs(9).unt = [77, 61, 63]; %check unit 65


figure();
for p = 1:numel(pairs);
    clf();
    nx = numel(pairs(p).unt);
    tid = pairs(p).tid;
    clm = pairs(p).clm;
    for u = 1:nx
        subplot2(2,nx,1,u);
        plot(pft{tid+u-1},      ...
             pairs(p).unt(u),...
             1,              ...
             'text',         ...
             clm);
        title(['Trial: ',num2str(tid+u-1),' Unit:',num2str(pairs(p).unt(u))]);
        subplot2(2,nx,2,u);
        plot(bfs{tid+u-1},...
             pairs(p).unt(u),...
             1,              ...
             'text',         ...
             clm,           ...
             'mazeMask',mask_bhv.mask);
    end
    waitforbuttonpress();
end



