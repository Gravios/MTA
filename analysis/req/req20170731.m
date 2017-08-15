
generate_nn_label_stats_multi_session('train',[],)
generate_nn_label_stats_multi_session('compute')
generate_nn_label_stats_multi_session('optimize')
generate_nn_label_stats_multi_session('display')

generate_nn_label_stats_multi_session('train',[],'fet_mis')
generate_nn_label_stats_multi_session('compute',[],'fet_mis')
generate_nn_label_stats_multi_session('optimize',[],'fet_mis')
generate_nn_label_stats_multi_session('display',[],'fet_mis')



% TESTING ------------------------------------------

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
stsRaw = load(fullfile(MTA_PROJECT_PATH,'analysis',[statsName,postProcessingTag,'.mat']));




figure,
sp = []
sp(end+1) = subplot(211);
plotSTC(Stc{1});
sp(end+1) = subplot(212);
plotSTC(StcCor{1});
linkaxes(sp,'xy');



s = 3
figure,
sp = []
sp(end+1) = subplot(211);
plotSTC(stsRaw.stc{1});
sp(end+1) = subplot(212);
plotSTC(stc{1});
linkaxes(sp,'xy');



