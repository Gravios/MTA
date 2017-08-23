
% TESTING ------------------------------------------
global MTA_PROJECT_PATH
varargin = [];
varargin = {,'nIter',100,'sampleRate',1000};
varargin = {'hand_labeled','fet_box','nIter',100,'sampleRate',1000};
MODEL_TYPE = 'multiSesPatNet';
% DEFARGS ----------------------------------------------------------------------------------------
defargs = struct('sessionList',                 '',                                            ...
                 'featureSet',                  'fet_bref_rev7',                               ...
                 'states',                      {{'walk','rear','turn','pause','groom','sit'}},...
                 'keys',                        {{'w','r','n','p','m','s'}},                   ...
                 'model',                       [],                                            ...
                 'sampleRate',                  10,                                            ...
                 'nNeurons',                    25,                                            ...
                 'nIter',                       10,                                            ...
                 'randomizationMethod',         'WSBNT',                                       ...
                 'map2reference',               true,                                          ...
                 'normalize',                   true,                                          ...
                 'referenceTrial',              'jg05-20120317.cof.all',                       ...
                 'trainingSessionList',         'hand_labeled',                                ...
                 'normalizationSessionList',    'hand_labeled',                                ...
                 'dropIndex',                   [],                                            ...
                 'prctTrain',                   90,                                            ...
                 'postProcessingTag',           ''                                             ...
);
[sessionList,featureSet,states,keys,model,sampleRate,nNeurons,nIter,randomizationMethod,       ...
 map2reference,normalize,referenceTrial,trainingSessionList,normalizationSessionList,          ...
 dropIndex,prctTrain,postProcessingTag] = DefaultArgs(varargin,defargs,'--struct');
% ------------------------------------------------------------------------------------------------    



generate_nn_label_stats_multi_session('train')
generate_nn_label_stats_multi_session('compute')
generate_nn_label_stats_multi_session('optimize')
generate_nn_label_stats_multi_session('display')

% COMPLETE 
% $$$ generate_nn_label_stats_multi_session('train','nIter',100)
% $$$ generate_nn_label_stats_multi_session('compute','nIter',100)
% $$$ generate_nn_label_stats_multi_session('optimize','nIter',100)
% $$$ generate_nn_label_stats_multi_session('display','nIter',100)

% COMPLETE 
% $$$ generate_nn_label_stats_multi_session('train',[],'fet_mis')
% $$$ generate_nn_label_stats_multi_session('compute',[],'fet_mis')
% $$$ generate_nn_label_stats_multi_session('optimize',[],'fet_mis')
% $$$ generate_nn_label_stats_multi_session('display',[],'fet_mis')

% COMPLETE 
% $$$ generate_nn_label_stats_multi_session('train','nIter',100)
% $$$ generate_nn_label_stats_multi_session('compute','nIter',100)
% $$$ generate_nn_label_stats_multi_session('optimize','nIter',100)
% $$$ generate_nn_label_stats_multi_session('display','nIter',100)

% INCOMPLETE 
generate_nn_label_stats_multi_session('train',   [],'fet_mis','nIter',100)
generate_nn_label_stats_multi_session('compute', [],'fet_mis','nIter',100)
generate_nn_label_stats_multi_session('optimize',[],'fet_mis','nIter',100)
generate_nn_label_stats_multi_session('display', [],'fet_mis','nIter',100)

if isempty(sessionList),
    sessionList = trainingSessionList; 
end

% LOAD neural network labeled state collections
statsName = ['MTAC_STATS+TRN+' trainingSessionList                                       ...
             '+LBS+' sessionList '+'                                                      ...
             featureSet                                                                   ...
             '+SR'   num2str(sampleRate)                                                  ...
             'NN'    num2str(nNeurons)                                                    ...
             'NI'    num2str(nIter)                                                       ...
             'M'     num2str(map2reference)                                               ...
             'MREF+' referenceTrial                                                       ...
             '+N'    num2str(normalize)                                                   ...
             'NREF+' normalizationSessionList                                             ...
             '+RND'  randomizationMethod                                                  ...
             '+PRCT' num2str(prctTrain)                                                   ...
             '+STS+' strjoin(keys,'')                                                     ...
             '-'     MODEL_TYPE];
disp(statsName);

generate_nn_label_stats_multi_session('optimize')
generate_nn_label_stats_multi_session('display')

stsRaw = load(fullfile(MTA_PROJECT_PATH,'analysis',[statsName,postProcessingTag,'.mat']));
stsOpt = load(fullfile(MTA_PROJECT_PATH,'analysis',[statsName,postProcessingTag,'_pp.mat']));
for s=1:6,stsRaw.labelingStats{s},end
for s=1:6,stsOpt.labelingStats(s),end
[[stsRaw.labelingStats(:).accuracy]',[stsOpt.labelingStats(:).accuracy]']
mean([[stsRaw.labelingStats(:).accuracy]',[stsOpt.labelingStats(:).accuracy]'])
diff([[stsRaw.labelingStats(:).accuracy]',[stsOpt.labelingStats(:).accuracy]'],1,2)



Trials = af(@(Trial)  MTATrial.validate(Trial), get_session_list(trainingSessionList));
StcHL = cf(@(Trial) Trial.load('stc'),Trials);
        cf(@(stc,states) set(stc,'states',stc(states{:})),...
             StcHL,repmat({states},[1,numel(Trials)]));


figure,
sp = []
sp(end+1) = subplot(211);
plotSTC(Stc{1});
sp(end+1) = subplot(212);
plotSTC(StcCor{1});
linkaxes(sp,'xy');



s = 5
figure,
sp = []
sp(end+1) = subplot(311);
plotSTC(StcHL{s});
sp(end+1) = subplot(312);
plotSTC(stsRaw.stc{s});
sp(end+1) = subplot(313);
plotSTC(stsOpt.stc{s});
linkaxes(sp,'xy');



% Smooth decision boundaries - 200 ms state minimum
bwin = round(.2*xyz.sampleRate)+double(mod(round(.2*xyz.sampleRate),2)==0);
mss = GetSegs(maxState,1:size(maxState,1),bwin,nan);
maxState=circshift(sq(mode(mss))',floor(bwin/2));


sessionList = 'hand_labeled';
trainingSessionList = 'hand_labeled_Ed';
generate_nn_label_stats_multi_session('compute',sessionList,'trainingSessionList',trainingSessionList);
generate_nn_label_stats_multi_session('optimize',sessionList,'trainingSessionList',trainingSessionList);
generate_nn_label_stats_multi_session('display',sessionList,'trainingSessionList',trainingSessionList);


% LOAD neural network labeled state collections
statsName = ['MTAC_STATS+TRN+' trainingSessionList                                       ...
             '+LBS+' sessionList '+'                                                      ...
             featureSet                                                                   ...
             '+SR'   num2str(sampleRate)                                                  ...
             'NN'    num2str(nNeurons)                                                    ...
             'NI'    num2str(nIter)                                                       ...
             'M'     num2str(map2reference)                                               ...
             'MREF+' referenceTrial                                                       ...
             '+N'    num2str(normalize)                                                   ...
             'NREF+' normalizationSessionList                                             ...
             '+RND'  randomizationMethod                                                  ...
             '+PRCT' num2str(prctTrain)                                                   ...
             '+STS+' strjoin(keys,'')                                                     ...
             '-'     MODEL_TYPE];
disp(statsName);

stsRaw = load(fullfile(MTA_PROJECT_PATH,'analysis',[statsName,postProcessingTag,'.mat']));
stsOpt = load(fullfile(MTA_PROJECT_PATH,'analysis',[statsName,postProcessingTag,'_pp.mat']));
for s=1:6,stsRaw.labelingStats(s),end
for s=1:6,stsOpt.labelingStats(s),end
[[cell2mat(cf(@(l) l.accuracy,stsRaw.labelingStats))]',[stsOpt.labelingStats(:).accuracy]']
mean([[stsRaw.labelingStats(:).accuracy]',[stsOpt.labelingStats(:).accuracy]'])
diff([[stsRaw.labelingStats(:).accuracy]',[stsOpt.labelingStats(:).accuracy]'],1,2)
mean(diff([[stsRaw.labelingStats(:).accuracy]',[stsOpt.labelingStats(:).accuracy]'],1,2))
