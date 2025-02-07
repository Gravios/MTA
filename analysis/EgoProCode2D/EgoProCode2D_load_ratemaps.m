% EgoProCode2D_f1_load_ratemaps.m
% EgoProCode2D_load_data.m


pft = cf(@(t,u)  pfs_2d_theta(t,u),  Trials, units);
placeFieldsNoRear = cf(@(t,u)  pfs_2d_theta(t,u,'theta-groom-sit-rear'),  Trials, units);
rmapNoRear = placeFieldsNoRear;

overwrite = false;
%%%<<< (pfe) ego ratemap
pfe = cf(@(t,u,x,s,p)                        ... Egocentric ratemaps
         compute_ego_ratemap                 ... func name
         (t,u,x,s,p,'overwrite',overwrite),  ... args
             Trials,                         ... MTATrial
             units,                          ... Unit subset, placefields proximal to maze center
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft                                               ... MTAApfs object, theta state placefields 
);
%%%>>>

tind = [1:30];
%%%<<< (pfet) egothp ratemap
pfet = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase.
         compute_egothp_ratemap(t,u,x,s,p,'overwrite',overwrite),...
             Trials(tind),                        ... MTATrial
             units(tind),                         ... Unit subset, placefields away from the maze walls
             xyz(tind),                           ... MTADxyz object, head position
             spk(tind),                           ... MTASpk object, spike time and id collection 
             pft(tind)                            ... MTAApfs object, theta state placefields 
);
%%%>>>



tind = [1:30];
%%%<<< (pfet) egothp ratemap
pfetf = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase.
         compute_egothp_ratemap(t,u,x,s,p,'overwrite',overwrite),...
             Trials(tind),                        ... MTATrial
             units(tind),                         ... Unit subset, placefields away from the maze walls
             xyz(tind),                           ... MTADxyz object, head position
             spk(tind),                           ... MTASpk object, spike time and id collection 
             pft(tind)                            ... MTAApfs object, theta state placefields 
);
%%%>>>


%%%<<< (pfs) egohba ratemap
pfs = cf(@(t,u,x,s,p)                       ... Egocentric ratemap | theta phase , head body angle.
compute_egohba_ratemap(t,u,x,s,p,'overwrite',overwrite),   ...
    Trials,                        ... MTATrial
units,                         ... Unit subset, placefields away from the maze walls
xyz,                           ... MTADxyz object, head position
spk,                           ... MTASpk object, spike time and id collection 
pft                            ... MTAApfs object, theta state placefields 
);
tind = [1:30];
% $$$ pfsh = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase , shuffled head body angle.
% $$$          compute_egohba_ratemap_shuffled(t,u,x,s,p,overwrite),   ...
% $$$              Trials(tind),                        ... MTATrial
% $$$              units(tind),                         ... Unit subset, placefields away from the maze walls
% $$$              xyz(tind),                           ... MTADxyz object, head position
% $$$              spk(tind),                           ... MTASpk object, spike time and id collection 
% $$$              pft(tind)                            ... MTAApfs object, theta state placefields 
% $$$ );
%%%>>>

tind = [1:30];
%tind = [20];
%%%<<< (pfs) egohvf ratemap
pfv = cf(@(t,u,x,s,p)                       ... Egocentric ratemap | theta phase , head body angle.
         compute_egohvf_ratemap(t,u,x,s,p,'overwrite',overwrite),   ...
             Trials(tind),                        ... MTATrial
             units(tind),                         ... Unit subset, placefields away from the maze walls
             xyz(tind),                           ... MTADxyz object, head position
             spk(tind),                           ... MTASpk object, spike time and id collection 
             pft(tind)                            ... MTAApfs object, theta state placefields 
);

% $$$ pfsh = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase , shuffled head body angle.
% $$$          compute_egohba_ratemap_shuffled(t,u,x,s,p,'overwrite',overwrite),   ...
% $$$              Trials,                        ... MTATrial
% $$$              units,                         ... Unit subset, placefields away from the maze walls
% $$$              xyz,                           ... MTADxyz object, head position
% $$$              spk,                           ... MTASpk object, spike time and id collection 
% $$$              pft                            ... MTAApfs object, theta state placefields 
% $$$ );
%%%>>>


tind = [1:30];
%tind = [18];
%%%<<< (pfs) egohvf ratemap
pfl = cf(@(t,u,x,s,p)                       ... Egocentric ratemap | theta phase , head body angle.
         compute_egohvl_ratemap(t,u,x,s,p,'overwrite',overwrite),   ...
             Trials(tind),                        ... MTATrial
             units(tind),                         ... Unit subset, placefields away from the maze walls
             xyz(tind),                           ... MTADxyz object, head position
             spk(tind),                           ... MTASpk object, spike time and id collection 
             pft(tind)                            ... MTAApfs object, theta state placefields 
);

% $$$ pfsh = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase , shuffled head body angle.
% $$$          compute_egohba_ratemap_shuffled(t,u,x,s,p,'overwrite',overwrite),   ...
% $$$              Trials,                        ... MTATrial
% $$$              units,                         ... Unit subset, placefields away from the maze walls
% $$$              xyz,                           ... MTADxyz object, head position
% $$$              spk,                           ... MTASpk object, spike time and id collection 
% $$$              pft                            ... MTAApfs object, theta state placefields 
% $$$ );
%%%>>>



%overwrite = true;
%tind = [1:30];
%%%<<< (pfs) ratemaps allo thp 
ratemapsAlloThp = cf(@(t,u,x,s)                 ... allocentric ratemap | theta phase
         compute_ratemaps_allo_thp(t,u,x,s,'overwrite',overwrite),   ...
             Trials(tind),                        ... MTATrial
             units(tind),                         ... Unit subset, placefields away from the maze walls
             xyz(tind),                           ... MTADxyz object, head position
             spk(tind)                            ... MTASpk object, spike time and id collection 
);
%%%>>>
