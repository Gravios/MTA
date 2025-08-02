% req20200217
%    Tags: place field head-body distance
%    Status: retired
%    Type: Analysis
%    Author: Justin Graboski
%    Final_Forms: NA
%    Project: MjgER2016: bhvfields
%    Description: place field rate dependence of head-body distance
%    Suplementary_Files: req20200217_args.m
%



%%%<<< LOAD MjgER2016 data
MjgER2016_load_data();
%  MjgER2016_load_data:
%
%  Variables:
%      Trials      units      cluSessionMap      pitchReferenceTrial      FigDir
%      sessionListName        sessionList        states                   numStates      
%      interpParPfsp          interpParDfs
%      
%  Functions:
%      reshape_eigen_vector
%%%>>>

%%%<<< SET MjgER2016 figure 3 global default arguments
req20200217_args();
% MjgER2016_figure3_args:
%
% Default arument overrides:
%     fet_HB_pitchB
%     fet_hbp_hba
%     compute_bhv_ratemaps
%%%>>>


%%%<<< LOAD data
% EXAMPLE Trial : jg05-20120312
tind = 21;
Trial = Trials{tind};
stc = Trials{tind}.stc.copy();
bfs = {};
pft = {};
for tind = 17:23;
    bfs{tind} = compute_bhv_ratemaps(Trials{tind},units{tind});
    pft{tind} = pfs_2d_theta(Trials{tind},units{tind});
end


tind = 22;
hfig = figure
u = units{tind}(1);
while u~=-1,
    subplot(121);
    plot(pft{tind},u,'mean','text',[],true);
    subplot(122);
    plot(bfs{tind},u,1,'colorbar',[],false);
    title(num2str(u));
    u = figure_controls(hfig,u,units{tind});
end

% LOAD place restricted behavior fields
bfs   = cf(@(t,u)   compute_bhv_ratemaps(t,u),  Trials, units);
bfsEx = bfs{tind};
bfsTag = {'HPITCHxBPITCH'};



bfs   = cf(@(t,u)   compute_bhv_ratemaps(t,u),  Trials, units);