% bpc_load_data - loads the folowing variables and functions
%
%  Variables:
%      Trials
%      Units
%      cluSessionMap
%      pitchReferenceTrial
%      figBasePath
%      sessionListName
%      sessionList
%      stcMode
%      states
%      numStates
%      interpParPfs
%      interpParDfs
%
%  Functions:
%      reshape_eigen_vector


% >>> GLOBAL VARS >>> ---------------------------------------------------------
configure_default_args();
% MjgER2016_figure_BhvPlacefields_args:
%
% Default arument overrides:
%     fet_HB_pitchB
%     compute_bhv_ratemaps
%     compute_bhv_ratemaps_shuffled
%     compute_bhv_ratemap_erpPCA
%
global MTA_FIGURES_PATH
% <<< GLOBAL VARS <<< ---------------------------------------------------------

% >>> LOCAL VARS >>> ----------------------------------------------------------
sessionListName = 'BehaviorPlaceCode';
sessionList = get_session_list_v3(sessionListName);

% >>> LOAD all Trials and cells into cellarray >>> ----------------------------
Trials = af(@(s) MTATrial.validate(s), sessionList);

% LOAD place cells for each Trial into cellarray
Units = cf(@(T)  T.spk.get_unit_set(T,'bhv'),  Trials);
Units = cf(@(T,U) remove_bad_units(T,U), Trials, Units);
UnitsInt = cf(@(T)  T.spk.get_unit_set(T,'inteneurons'),  Trials);
% electrode anatomical location
AnatLocCA1 = cell2mat(cf(@(T) T.meta.anat_loc.CA1, Trials));
AnatLocCA3 = cell2mat(cf(@(T) T.meta.anat_loc.CA3, Trials));

% This only works with this data set 
noUnits = cellfun(@isempty,Units);
AnatLocCA1(noUnits) = [];
AnatLocCA1(noUnits) = [];
Trials(noUnits) = [];
Units(noUnits) = [];
% <<< LOAD all Trials and cells into cellarray <<< ----------------------------

pitchReferenceTrial = 'Ed05-20140529.ont.all';

% SET helper function to reshape eigenvectors
reshape_eigen_vector =                                   ...
    @(V,p)  reshape(V(:,1),p{1}.adata.binSizes')';
 
% 
figBasePath =  ...
    create_directory('/storage/gravio/figures/analysis/BehaviorPlaceCode/');

% >>> Behavioral states >>> ---------------------------------------------------
% SET states to plot
stcMode = 'msnn_ppsvd_raux';
states = {'theta-groom-sit', ...
          'rear&theta',      ...
          'hloc&theta',      ...
          'hpause&theta',    ...
          'lloc&theta',      ...
          'lpause&theta'     ...
};

numStates = numel(states);
% <<< Behavioral states <<< ---------------------------------------------------

% >>> Place Field interpolation definition >>> --------------------------------
interpParPfs = struct(                              ...
    'bins',             {{linspace(-500,500,100)',  ...
                          linspace(-500,500,100)'}},...
    'nanMaskThreshold', 0,                          ...
    'methodNanMap',     'linear',                   ...
    'methodRateMap',    'linear');

interpParDfs = struct(                              ...
    'bins',             {{linspace(-2,2,100)',      ...
                          linspace(-2,2,100)'}},    ...
    'nanMaskThreshold', 0,                          ...
    'methodNanMap',     'linear',                   ...
    'methodRateMap',    'linear');
% <<< Place Field interpolation definition <<< --------------------------------

% <<< LOCAL VARS <<< ----------------------------------------------------------

% >>> CONCATENATE all maps >>> ------------------------------------------------
% ADD global unit id colum
cluSessionMap = [];
for u = 1:numel(Units)
    cluSessionMap = cat( 1,                              ...
                         cluSessionMap,                  ...
                         [u * ones([numel(Units{u}),1]), ...
                          Units{u}(:)]                   ...
    );
end

% ADD global unit id colum
cluSessionMapInt = [];
for u = 1:numel(Units)
    cluSessionMapInt = cat( 1,                              ...
                         cluSessionMapInt,                  ...
                         [u * ones([numel(UnitsInt{u}),1]), ...
                          UnitsInt{u}(:)]                   ...
    );
end
% <<< CONCATENATE all maps <<< ------------------------------------------------



