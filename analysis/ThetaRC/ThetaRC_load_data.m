% ThetaRC_load_data - loads the folowing variables and functions
%
%  Variables:
%      Trials
%      units
%      cluSessionMap
%      pitchReferenceTrial
%      figBasePath
%      sessionListName
%      sessionList
%      states
%      numStates
%


project.name = 'ThetaRC';



%% GLOBAL VARIABLES ------------------------------------------------------------------------------
global MTA_FIGURES_PATH
global PROJECT_FIGURE_PATH

PROJECT_FIGURE_PATH = create_directory(getenv(['MTA_PROJECT_' upper(project.name)]));



%% LOCAL PROJECT VARIABLES -----------------------------------------------------------------------
% LOAD Trials and metadata
sessionListName = 'MjgER2016';
sessionList = get_session_list_v2(sessionListName);
Trials = af(@(s) MTATrial.validate(s), sessionList);

% LOAD place cell ids -> one cell per Trial
Units = cf(@(T)  T.spk.get_unit_set(T,'placecells'),  Trials);

% REMOVE place cells with poor sampling or low rates
Units = cf(@(T,U) remove_bad_units(T,U), Trials, Units);

% LOAD interneuron ids -> one cell per Trial
UnitsInts = cf(@(T)  T.spk.get_unit_set(T,'interneurons'),  Trials); 

% LOAD place cell ids for egocentric project -> one cell per Trial
UnitsEgo = cf(@(T)  T.spk.get_unit_set(T,'egocentric'),  Trials); 


% SET the reference trial for computing the inter-subject head pitch correction 
pitchReferenceTrial = 'Ed05-20140529.ont.all';


% SET states
clear('states')
states.name = 'msnn_ppsvd_raux';
states.labels = {{'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta',...
          'lpause&theta'}};

% MAP 
cluSessionMap = [];
for u = 1:numel(Units)
    cluSessionMap = cat(1,cluSessionMap,[u*ones([numel(Units{u}),1]), Units{u}(:)]);
end


sampleRate = 250;
%headCenterCorrection = [-25,-8];
pfsState = 'theta-groom-sit-rear';
hbaBinEdges = -1.5:0.6:1.5;
xyz = cf(@(t) preproc_xyz(t,'trb'),             Trials);
      cf(@(x) x.filter('ButFilter',3,30,'low'), xyz);    
      cf(@(x) x.resample(sampleRate),           xyz);
spk = cf(@(t,u) t.load('spk',sampleRate,'gper',u,'deburst'), Trials, Units);    

