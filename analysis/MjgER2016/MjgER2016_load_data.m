% MjgER2016_load_data - loads the folowing variables and functions
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


%%%<<< SET default arguments -----------------------------------------------------------------------

configure_default_args();
% MjgER2016_figure_BhvPlacefields_args:
%
% Default arument overrides:
%     fet_HB_pitchB
%     compute_bhv_ratemaps
%     compute_bhv_ratemaps_shuffled
%     compute_bhv_ratemap_erpPCA
%%%>>>


global MTA_FIGURES_PATH

% RETRIVE struct array of meta data for all Trials
sessionListName = 'MjgER2016';
sessionList = get_session_list_v2(sessionListName);

% LOAD all Trials into cellarray
Trials = af(@(s) MTATrial.validate(s), sessionList);

% LOAD place cells for each Trial into cellarray
Units = cf(@(T)  T.spk.get_unit_set(T,'placecells'),  Trials); 
Units = cf(@(T,U) remove_bad_units(T,U), Trials,Units);

% LOAD interneuron ids -> one cell per Trial
UnitsInts = cf(@(T)                                      ...
               T.spk.get_unit_set(T,'interneurons'),     ...
               Trials); 

% CONCATENATE all maps
% ADD global unit id colum
cluSessionMap = [];
for u = 1:numel(Units)
    cluSessionMap = cat( 1,                              ...
                         cluSessionMap,                  ...
                         [u * ones([numel(Units{u}),1]), ...
                          Units{u}(:)]                   ...
    );
end

pitchReferenceTrial = 'Ed05-20140529.ont.all';


% SET helper function to reshape eigenvectors
reshape_eigen_vector =                                   ...
    @(V,p)  reshape(V(:,1),p{1}.adata.binSizes')';
 
% 
figBasePath =  ...
    create_directory('/storage/gravio/figures/analysis/MjgER2016/');

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

