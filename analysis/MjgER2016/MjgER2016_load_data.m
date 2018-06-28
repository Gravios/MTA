% MjgER2016_load_data - loads the folowing variables and functions
%
%  Variables:
%      Trials
%      units
%      cluSessionMap
%      pitchReferenceTrial
%      FigDir
%      sessionListName
%      sessionList
%      states
%      numStates
%      interpParPfsp
%      interpParDfs
%
%  Functions:
%      reshape_eigen_vector

sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);

Trials = af(@(s) MTATrial.validate(s), sessionList);
units = cf(@(T)  select_placefields(T),  Trials); 
units = req20180123_remove_bad_units(units);

cluSessionMap = [];
for u = 1:numel(units)
    cluSessionMap = cat(1,cluSessionMap,[u*ones([numel(units{u}),1]),units{u}(:)]);
end

pitchReferenceTrial = 'Ed05-20140529.ont.all';


% SET helper function to reshape eigenvectors
reshape_eigen_vector = @(V,pfd) reshape(V(:,1),pfd{1}.adata.binSizes')';
 
% 
FigDir = create_directory('/storage/gravio/figures/analysis/parts/MjgER2016/');

% SET states to plot
states = {'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta',...
          'lpause&theta'};
numStates = numel(states);


interpParPfs = struct('bins',{{linspace(-500,500,100),linspace(-500,500,100)}},...
                   'nanMaskThreshold', 0,...
                   'methodNanMap',     'linear',...
                   'methodRateMap',    'linear');

interpParDfs = struct('bins',{{linspace(-2,2,100),linspace(-2,2,100)}},...
                   'nanMaskThreshold', 0,...
                   'methodNanMap',     'linear',...
                   'methodRateMap',    'linear');

