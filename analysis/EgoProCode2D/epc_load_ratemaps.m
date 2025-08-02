% epc_load_data.m


tind = ':';
Pfnr = cf(@(t,u)  pfs_2d_theta(t,u,'theta-groom-sit-rear'),  Trials, Units);
rmapNoRear = placeFieldsNoRear;

overwrite = true;
% >>> (pfe) ego ratemap >>> ---------------------------------------------------
pfe = cf(                               ... cellfun
    @(T,U,X,S,P)                        ... Egocentric ratemaps
    compute_ego_ratemap                 ... function 
    (T,U,X,S,P,'overwrite',overwrite),  ... args
         Trials(tind),                  ... { MTATrial ojbects }
         Units(tind),                   ... { array[numeric] }
         Xyz(tind),                     ... { MTADxyz objects }
         Spk(tind),                     ... { MTASpk objects }
         Pfnr(tind)                     ... { MTAApfs objects }
);
% <<< (pfe) ego ratemap <<< ---------------------------------------------------

% >>> (pfet) egothp ratemap >>> -----------------------------------------------
pfet = cf( ...
    @(T,U,X,S,P)                      ... Egocentric ratemap | theta phase.
    compute_egothp_ratemap ...
    (T,U,X,S,P,'overwrite',overwrite),...
         Trials(tind),                        ... MTATrial
         Units(tind),                         ... Unit subset, placefields away from the maze walls
         Xyz(tind),                           ... MTADxyz object, head position
         Spk(tind),                           ... MTASpk object, spike time and id collection 
         Pfnr(tind)                            ... MTAApfs object, theta state placefields 
);
% <<< (pfet) egothp ratemap <<< -----------------------------------------------



% >>> (pfet) egothp ratemap >>> -----------------------------------------------
% $$$ pfetf = cf(@(T,U,X,S,P)                      ... Egocentric ratemap | theta phase.
% $$$          compute_egothp_ratemap(T,U,X,S,P,'overwrite',overwrite),...
% $$$              Trials(tind),                        ... MTATrial 
% $$$             Units(tind),                         ... Unit subset, placefields away from the maze walls
% $$$              Xyz(tind),                           ... MTADxyz object, head position
% $$$              Spk(tind),                           ... MTASpk object, spike time and id collection 
% $$$              Pfnr(tind)                            ... MTAApfs object, theta state placefields 
% $$$ );
% >>> (pfet) egothp ratemap <<< -----------------------------------------------


% >>> (pfs) egohba ratemap >>> ------------------------------------------------
pfs = cf( ...
    @(T,U,X,S,P)          ... Egocentric ratemap | theta phase , head body angle.
    compute_egohba_ratemap(T,U,X,S,P,'overwrite',overwrite),   ...
    Trials(tind),                        ... MTATrial
    Units(tind),                         ... Unit subset, placefields away from the maze walls
    Xyz(tind),                           ... MTADxyz object, head position
    Spk(tind),                           ... MTASpk object, spike time and id collection 
    Pfnr(tind)                            ... MTAApfs object, theta state placefields 
);

pfsh = cf(@(T,U,X,S,P)                      ... Egocentric ratemap | theta phase , shuffled head body angle.
         compute_egohba_ratemap_shuffled(T,U,X,S,P,overwrite),   ...
             Trials(tind),                        ... MTATrial
             Units(tind),                         ... Unit subset, placefields away from the maze walls
             Xyz(tind),                           ... MTADxyz object, head position
             Spk(tind),                           ... MTASpk object, spike time and id collection 
             Pfnr(tind)                            ... MTAApfs object, theta state placefields 
);
% <<< (pfs) egohba ratemap <<< ------------------------------------------------

% >>> (pfs) egohvf ratemap
pfv = cf(@(T,U,X,S,P)                       ... Egocentric ratemap | theta phase , head body angle.
         compute_egohvf_ratemap(T,U,X,S,P,'overwrite',overwrite),   ...
             Trials(tind),                        ... MTATrial
             Units(tind),                         ... Unit subset, placefields away from the maze walls
             Xyz(tind),                           ... MTADxyz object, head position
             Spk(tind),                           ... MTASpk object, spike time and id collection 
             Pfnr(tind)                            ... MTAApfs object, theta state placefields 
);

pfvsh = cf(@(T,U,X,S,P)                      ... Egocentric ratemap | theta phase , shuffled head body angle.
         compute_egohvf_ratemap_shuffled(T,U,X,S,P,'overwrite',overwrite),   ...
             Trials,                        ... MTATrial
             Units,                         ... Unit subset, placefields away from the maze walls
             xyz,                           ... MTADxyz object, head position
             Spk,                           ... MTASpk object, spike time and id collection 
             Pfnr                            ... MTAApfs object, theta state placefields 
);
% <<<

% >>> (pfs) egohvf ratemap
pfl = cf(@(T,U,X,S,P)                       ... Egocentric ratemap | theta phase , head body angle.
         compute_egohvl_ratemap(T,U,X,S,P,'overwrite',overwrite),   ...
             Trials(tind),                        ... MTATrial
             Units(tind),                         ... Unit subset, placefields away from the maze walls
             Xyz(tind),                           ... MTADxyz object, head position
             Spk(tind),                           ... MTASpk object, spike time and id collection 
             Pfnr(tind)                            ... MTAApfs object, theta state placefields 
);

pfsh = cf(@(T,U,X,S,P)                      ... Egocentric ratemap | theta phase , shuffled head body angle.
         compute_egohba_ratemap_shuffled(T,U,X,S,P,'overwrite',overwrite),   ...
             Trials,                        ... MTATrial
             Units,                         ... Unit subset, placefields away from the maze walls
             xyz,                           ... MTADxyz object, head position
             Spk,                           ... MTASpk object, spike time and id collection 
             Pfnr                            ... MTAApfs object, theta state placefields 
);
% <<<

% >>> (pfs) ratemaps allo thp 
ratemapsAlloThp = cf(@(t,u,x,s)                 ... allocentric ratemap | theta phase
         compute_ratemaps_allo_thp(t,u,x,s,'overwrite',overwrite),   ...
             Trials(tind),                        ... MTATrial
             Units(tind),                         ... Unit subset, placefields away from the maze walls
             Xyz(tind),                           ... MTADxyz object, head position
             Spk(tind)                            ... MTASpk object, spike time and id collection 
);
% <<<
