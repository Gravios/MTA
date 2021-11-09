% MjgER2016_load_data - loads the folowing variables and functions
%
%  Variables:
%      Trials
%      units
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

global MTA_FIGURES_PATH

sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);

Trials = af(@(s) MTATrial.validate(s), sessionList);

% $$$ units = cf(@(T)  select_placefields(T),  Trials); 
% $$$ units = req20180123_remove_bad_units(units);
% $$$ cf(@(T,U)  T.spk.set_unit_set(T,'placecells',U),  Trials, units); 


units = cf(@(T)  T.spk.get_unit_set(T,'placecells'),  Trials); 
units = cf(@(T,U) remove_bad_units(T,U), Trials,units);

%units = cf(@(T)  T.spk.get_unit_set('placecells'),  Trials); 


phzCorrection = ...
    [pi*1.25,pi*1.25,                              ... er01 CA3
     pi/2,pi/2,pi/2,                               ... ER06 CA1
     pi/1.25,pi/1.25,                              ... Ed10 CA3
     0,0,0,0,0,                                    ... jg04 CA1
     0,0,0,0,                                      ... jg04 CA3
     pi/4,pi/4,pi/4,pi/4,pi/4,pi/4,pi/4,pi/4,pi/4, ... jg05 CA1
     pi,                                           ... ER06 CA3
     pi/1.25,                                      ... Ed10 CA3
     pi*1.25,                                      ... er01 CA3
     pi/4,                                         ... FS03 CA1
     -pi/8];                                         % jg05 CA3

% GUESS 
% $$$ hbangCorrection =                                  ...
% $$$     [0,0,                                          ... er01
% $$$      0.2,0.2,0.2,                                  ... ER06
% $$$      0,0,                                          ... Ed10
% $$$      0,0,0,0,0,0,0,0,0,                            ... jg04
% $$$      -0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2, ... jg05                   
% $$$      0,                                            ... ER06
% $$$      0,                                            ... Ed10
% $$$      0,                                            ... er01
% $$$      0.1,                                          ... FS03
% $$$      -0.2];                                          % jg05
% $$$ headRotation = ...
% $$$     {0,0,                                           ... er01
% $$$      0,0,0,                                         ... ER06
% $$$      0,0,                                           ... Ed10
% $$$      0,0,0,0,0,0,0,0,0,                             ... jg04
% $$$      0.17, 0.17, 0.17,0.17,0.17, 0.17, 0.17,0.17,0.17,... jg05
% $$$      0,                                             ... ER06
% $$$      0,                                             ... Ed10
% $$$      0,                                             ... er01
% $$$      0,                                             ... FS03
% $$$      0.17};                                           % jg05
% $$$ headRollCorrection =                                     ...
% $$$     {0,0,                                                ... er01
% $$$      -0.25,-0.25,-0.25,                                  ... ER06
% $$$      0,0,                                                ... Ed10
% $$$      0,0,0,0,0,0,0,0,0,                                  ... jg04
% $$$      -0.22, -0.22, -0.22, -0.22,-0.22, -0.22, -0.22,     ... jg05
% $$$      -0.22,-0.22,0,0,0,...% new units - jg05, jg05, ER06, Ed10, er01
% $$$      0,... FS03
% $$$      -0.22};% jg05

% EMPIRICAL 
headBodyCorrection =                                 ...
    {-0.084,-0.084,                                  ... er01
     0.156,0.156,0.156,                              ... ER06
     -0.077,-0.077,                                  ... Ed10
     0.16,0.16,0.16,0.16,0.16,0.16,0.16,0.16,0.16,   ... jg04
     -0.234,-0.234,-0.234,-0.234,-0.234,-0.234,      ... jg05
         -0.234,-0.234,-0.234,                       ... jg05
     0.156,                                          ... ER06
     -0.077,                                         ... Ed10
     -0.084,                                         ... er01
     0.115,                                          ... FS03
     -0.234};                                          % jg05
hbangCorrection = headBodyCorrection;

headRotationCorrection = ...
    {0.225,0.251,                                    ... er01
     -0.155,-0.155,-0.155,                           ... ER06
     0.034,0.034,                                    ... Ed10
     -0.136,-0.136,-0.136,-0.136,-0.136,-0.136,      ... jg04
         -0.136,-0.136,-0.136,                       ... jg04
     0.264,0.264,0.264,0.264,0.264,0.264,0.264,      ... jg05
         0.264,0.264,                                ... jg05
     -0.155,                                         ... ER06
     0.034,                                          ... Ed10
     0.251,                                          ... er01
     0.0,                                            ... FS03
     0.264};                                           % jg05
headRotation = headRotationCorrection;

headRollCorrection =                                 ...
    {-0.012,-0.376,                                  ... er01
     -0.164,-0.164,-0.164,                           ... ER06
     -0.187,-0.187,                                  ... Ed10
     0.194,0.194,0.194,0.194,0.194,0.194,0.194,0.194,... jg04
         0.194,                                      ... jg04
     -0.365,-0.365,-0.365,-0.365,-0.365,-0.365,      ... jg05
         -0.365,-0.365,-0.365,                       ... jg05
     -0.164,                                         ... ER06 
     -0.187,                                         ... Ed10
     -0.376,                                         ... er01
     -0.127,                                         ... FS03
     -0.365};                                          % jg05


unitsInts = cf(@(T)  T.spk.get_unit_set(T,'interneurons'),  Trials); 




cluSessionMap = [];
for u = 1:numel(units)
    cluSessionMap = cat(1,cluSessionMap,[u*ones([numel(units{u}),1]),units{u}(:)]);
end

pitchReferenceTrial = 'Ed05-20140529.ont.all';


% SET helper function to reshape eigenvectors
reshape_eigen_vector = @(V,p) reshape(V(:,1),p{1}.adata.binSizes')';
 
% 
figBasePath = create_directory('/storage/gravio/figures/analysis/MjgER2016/');

% SET states to plot
stcMode = 'msnn_ppsvd_raux';
states = {'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta',...
          'lpause&theta'};
numStates = numel(states);


interpParPfs = struct('bins',{{linspace(-500,500,100)',linspace(-500,500,100)'}},...
                   'nanMaskThreshold', 0,...
                   'methodNanMap',     'linear',...
                   'methodRateMap',    'linear');

interpParDfs = struct('bins',{{linspace(-2,2,100)',linspace(-2,2,100)'}},...
                   'nanMaskThreshold', 0,...
                   'methodNanMap',     'linear',...
                   'methodRateMap',    'linear');

