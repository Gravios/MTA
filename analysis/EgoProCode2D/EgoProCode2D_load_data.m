% EgoProCode2D_load_data - loads the folowing variables and functions
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
%      interpParPfs
%      interpParDfs
%
%  Functions:
%      reshape_eigen_vector


project.name = 'EgoProCode2d';



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
units = cf(@(T)  T.spk.get_unit_set(T,'placecells'),  Trials);
% REMOVE place cells with poor sampling or low rates
units = cf(@(T,U) remove_bad_units(T,U), Trials,units);

% LOAD interneuron ids -> one cell per Trial
unitsInts = cf(@(T)  T.spk.get_unit_set(T,'interneurons'),  Trials); 

% LOAD place cell ids for egocentric project -> one cell per Trial
unitsEgo = cf(@(T)  T.spk.get_unit_set(T,'egocentric'),  Trials); 

% SET index arry for each anatomical location in oder of cumulative
% index within cells
% $$$ unitsEgoCA1 = [1:19,44:149];
% $$$ unitsEgoCA3 = [20:41,123:127,150:164];
unitsEgoCA1 = [1:19,44:122,134:149];
unitsEgoCA3 = [20:43,123:133,150:164];


% MAP 
cluSessionMap = [];
for u = 1:numel(units)
    cluSessionMap = cat(1,cluSessionMap,[u*ones([numel(units{u}),1]),units{u}(:)]);
end

egoCluSessionMap = [];
for u = 1:numel(unitsEgo)
    egoCluSessionMap = cat(1,egoCluSessionMap,[u*ones([numel(unitsEgo{u}),1]),unitsEgo{u}(:)]);
end


% SET the reference trial for computing the inter-subject head pitch correction 
pitchReferenceTrial = 'Ed05-20140529.ont.all';


% SET states 
states.name = 'msnn_ppsvd_raux';
states.labels = {{'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta',...
          'lpause&theta'}};

interpParPfs = struct('bins',{{linspace(-500,500,100)',linspace(-500,500,100)'}},...
                   'nanMaskThreshold', 0,...
                   'methodNanMap',     'linear',...
                   'methodRateMap',    'linear');

interpParDfs = struct('bins',{{linspace(-2,2,100)',linspace(-2,2,100)'}},...
                   'nanMaskThreshold', 0,...
                   'methodNanMap',     'linear',...
                   'methodRateMap',    'linear');

sampleRate = 250;
%headCenterCorrection = [-25,-8];
pfsState = 'theta-groom-sit-rear';
hbaBinEdges = -1.5:0.6:1.5;
xyz = cf(@(t) preproc_xyz(t,'trb'),             Trials);
      cf(@(x) x.filter('ButFilter',3,30,'low'), xyz);    
      cf(@(x) x.resample(sampleRate),           xyz);
spk = cf(@(t,u) t.load('spk',sampleRate,'gper',u,'deburst'),Trials,units);    


EgoProCode2D_load_ratemaps();


%% HELPER FUNCTIONS ----------------------------------------------------------------------------
% SET helper function to reshape eigenvectors
reshape_eigen_vector = @(V,p) reshape(V(:,1),p{1}.adata.binSizes')';
