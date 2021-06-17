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

sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);

Trials = af(@(s) MTATrial.validate(s), sessionList);
units = cf(@(T)  select_placefields(T),  Trials); 
units = req20180123_remove_bad_units(units);


phzCorrection = [pi*1.25,pi*1.25,                             ... er01
                 pi/2,pi/2,pi/2,                    ... ER06
                 pi/1.25,pi/1.25,                             ... Ed10
                 0,0,0,0,0,0,0,0,0,                 ... jg04                 
                 pi/4,pi/4,pi/4,pi/4,pi/4,pi/4,pi/4 ... jg05
                 pi/4,pi/4,pi,pi/1.25,pi*1.25]; % new units - jg05, jg05, ER06, Ed10, er01

headRotation = {0,0,                                         ... er01
                0,0,0,                                       ... ER06
                0,0,                                         ... Ed10
                0,0,0,0,0,0,0,0,0,                           ... jg04
                0.17, 0.17, 0.17, 0.17,0.17, 0.17, 0.17,     ... jg05
                0.17,0.17,0,0,0};% new units - jg05, jg05, ER06, Ed10, er01


unitsInts = {...
    [ 31, 78, 82,125,169,195,203],...                                    er01 20110719
    [ 31, 76,105,147],...                                                er01 20110721
    [ 27, 32, 68, 69,105,124,125,126,128,182,220,221,222,225],...        ER06 20130612
    [ 10, 13, 15, 18, 60, 93,112,113,192,213,215,216,220],...            ER06 20130613
    [  5, 10, 11, 12, 22, 43, 49, 91, 94, 95,110,124,144,145,146,155,... ER06 20130614
       156,168,169,170,171,172,176,177,179,181,182],...                  
    [  9, 11, 12, 29, 30, 43, 45, 46, 51, 52, 54, 55, 56, 87],...        Ed10 20140816
    [ 11, 12, 31, 47, 48, 51, 56, 80, 81, 89],...                        Ed10 20140817
    [  8,  9, 16],...                                                    jg04 20120128
    [ 21],...                                                            jg04 20120129
    [ 24],...                                                            jg04 20120130
    [ 10, 24, 27],...                                                    jg04 20120131
    [ 10],...                                                            jg04 20120201
    [  4,  5],...                                                        jg04 20120210
    [  4],...                                                            jg04 20120211
    [  6],...                                                            jg04 20120212
    [  2],...                                                            jg04 20120213
    [  5, 10, 15, 27, 28, 38, 64, 66,100,114,116,117,121,122],...        jg05 20120309
    [  4,  6,  7,  8, 27, 28, 43, 59, 71, 86, 99,100,101,102],...        jg05 20120310    
    [  2,  7,  9, 25, 41, 44, 46, 70, 71, 72, 98,112,113,118,...         jg05 20120311
     119,135,144,148,188,193,203,204,205],...
    [  3,  7,  8, 15, 16, 43, 45, 50, 76, 77, 92,106,124,184],...        jg05 20120312
    [  1,  5, 34, 60, 69],...                                            jg05 20120315
    [ 17, 28, 49, 55],...                                                jg05 20120316
    [ 15, 16, 17, 52],...                                                jg05 20120317
    [  2,  4, 20, 30, 39, 40, 51],...                                    jg05 20120323
    [  7, 14, 16, 31, 34],...                                            jg05 20120324
    [  5, 13, 27, 66, 71, 72, 75, 91, 94,105,127,157,209,232,233,251,254],...ER06 20130624
    [  4, 17, 20, 22, 23, 25, 29 ,30, 31, 40, 66, 69, 72],...            Ed10 20140815
    [ 93,175] ...                                                        er01 20110722
};


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

