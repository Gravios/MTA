
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
                 'nIter',                       100,                                            ...
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
% $$$ generate_nn_label_stats_multi_session('train',   [],'fet_bref_rev7','nIter',100)
% $$$ generate_nn_label_stats_multi_session('compute', [],'fet_bref_rev7','nIter',100)
% $$$ generate_nn_label_stats_multi_session('optimize',[],'fet_bref_rev7','nIter',100)
% $$$ generate_nn_label_stats_multi_session('display', [],'fet_bref_rev7','nIter',100)

% COMPLETE 
% $$$ generate_nn_label_stats_multi_session('train',   [],'fet_mis','nIter',100)
% $$$ generate_nn_label_stats_multi_session('compute', [],'fet_mis','nIter',100)
% $$$ generate_nn_label_stats_multi_session('optimize',[],'fet_mis','nIter',100)
% $$$ generate_nn_label_stats_multi_session('display', [],'fet_mis','nIter',100)

% COMPLETE 
% $$$ generate_nn_label_stats_multi_session('train',   [],'fet_bref_SMSU');
% $$$ generate_nn_label_stats_multi_session('compute', [],'fet_bref_SMSU');
% $$$ generate_nn_label_stats_multi_session('optimize',[],'fet_bref_SMSU');
% $$$ generate_nn_label_stats_multi_session('display', [],'fet_bref_SMSU');

% COMPLETE 
% $$$ generate_nn_label_stats_multi_session('train',   [],'fet_bref_TH');
% $$$ generate_nn_label_stats_multi_session('compute', [],'fet_bref_TH');
% $$$ generate_nn_label_stats_multi_session('optimize',[],'fet_bref_TH');
% $$$ generate_nn_label_stats_multi_session('display', [],'fet_bref_TH');

% COMPLETE 
% $$$ generate_nn_label_stats_multi_session('train',   [],'fet_bref_rev8','nIter',10)
% $$$ generate_nn_label_stats_multi_session('compute', [],'fet_bref_rev8','nIter',10)
% $$$ generate_nn_label_stats_multi_session('optimize',[],'fet_bref_rev8','nIter',10)
% $$$ generate_nn_label_stats_multi_session('display', [],'fet_bref_rev8','nIter',10)

% COMPLETE 
% $$$ generate_nn_label_stats_multi_session('train',   [],'fet_bref_rev8','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('compute', [],'fet_bref_rev8','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('optimize',[],'fet_bref_rev8','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('display', [],'fet_bref_rev8','nIter',10,'randomizationMethod','WSBT')

% COMPLETE 
% $$$ generate_nn_label_stats_multi_session('train',   [],'fet_bref_rev9','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('compute', [],'fet_bref_rev9','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('optimize',[],'fet_bref_rev9','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('display', [],'fet_bref_rev9','nIter',10,'randomizationMethod','WSBT')

% COMPLETE 
% $$$ generate_nn_label_stats_multi_session('train',   [],'fet_bref_rev10','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('compute', [],'fet_bref_rev10','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('optimize',[],'fet_bref_rev10','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('display', [],'fet_bref_rev10','nIter',10,'randomizationMethod','WSBT')

% COMPLETE 
% $$$ generate_nn_label_stats_multi_session('train',   [],'fet_bref_rev11','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('compute', [],'fet_bref_rev11','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('optimize',[],'fet_bref_rev11','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('display', [],'fet_bref_rev11','nIter',10,'randomizationMethod','WSBT')

% COMPLETE 
% $$$ generate_nn_label_stats_multi_session('train',   [],'fet_bref_rev12','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('compute', [],'fet_bref_rev12','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('optimize',[],'fet_bref_rev12','nIter',10,'randomizationMethod','WSBT')
% $$$ generate_nn_label_stats_multi_session('display', [],'fet_bref_rev12','nIter',10,'randomizationMethod','WSBT')

% COMPLETE 
generate_nn_label_stats_multi_session('train',   [],'fet_mis','nIter',10,'randomizationMethod','WSBT')
generate_nn_label_stats_multi_session('compute', [],'fet_mis','nIter',10,'randomizationMethod','WSBT')
generate_nn_label_stats_multi_session('optimize',[],'fet_mis','nIter',10,'randomizationMethod','WSBT')
generate_nn_label_stats_multi_session('display', [],'fet_mis','nIter',10,'randomizationMethod','WSBT')

% INCOMPLETE 
generate_nn_label_stats_multi_session('train',   [],'fet_mis','nIter',100,'randomizationMethod','WSBT')
generate_nn_label_stats_multi_session('compute', [],'fet_mis','nIter',100,'randomizationMethod','WSBT')
generate_nn_label_stats_multi_session('optimize',[],'fet_mis','nIter',100,'randomizationMethod','WSBT')
generate_nn_label_stats_multi_session('display', [],'fet_mis','nIter',100,'randomizationMethod','WSBT')


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


% WEIGHTED network output


% COMPLETE 
% $$$ generate_nn_label_stats_multi_session_weighted('train',   [],'fet_bref_rev7','nIter',100)
% $$$ generate_nn_label_stats_multi_session_weighted('compute', [],'fet_bref_rev7','nIter',100)
% $$$ generate_nn_label_stats_multi_session_weighted('optimize',[],'fet_bref_rev7','nIter',100)
generate_nn_label_stats_multi_session_weighted('display', [],'fet_bref_rev7','nIter',100)


% $$$ generate_nn_label_stats_multi_session_weighted('train',   [],'fet_mis','nIter',100,'randomizationMethod','WSBT')
generate_nn_label_stats_multi_session_weighted('compute', [],'fet_mis','nIter',100,'randomizationMethod','WSBT')
generate_nn_label_stats_multi_session_weighted('optimize',[],'fet_mis','nIter',100,'randomizationMethod','WSBT')
generate_nn_label_stats_multi_session_weighted('display', [],'fet_mis','nIter',100,'randomizationMethod','WSBT')



% Score fragmentation

% number of fragments within hand labeled periods measured by counting the number of neural-network
% label centers within each hand label.

Trial = MTATrial.validate('jg05-20120317.cof.all');
states = {'walk','rear','turn','pause','groom','sit'};
stch = Trial.load('stc','hand_labeled_rev3_jg');
stcw = Trial.load('stc','MTAC_BATCH+hand_labeled+1+fet_mis+SR10NN25NI100M1MREF+jg05-20120317.cof.all+N1NREF+hand_labeled+RNDWSBT+PRCT90+STS+wrnpms+multiSesPatNet_weighted');
stcp = Trial.load('stc','MTAC_BATCH+hand_labeled+1+fet_mis+SR10NN25NI100M1MREF+jg05-20120317.cof.all+N1NREF+hand_labeled+RNDWSBT+PRCT90+STS+wrnpms+multiSesPatNet_weighted_ppsvd');
stcnp = Trial.load('stc','MTAC_BATCH+hand_labeled+1+fet_mis+SR10NN25NI100M1MREF+jg05-20120317.cof.all+N1NREF+hand_labeled+RNDWSBT+PRCT90+STS+wrnpms+multiSesPatNet_ppsvd');

stcn = Trial.load('stc','msnn_ppsvd');


figure(); plot_stcs(stch,stcnp,stcp);



sres = [];
sbhv = [];
for s = 1:numel(states)
    sper = stcp{states{s},119.881035};
    sbhv = cat(1,sbhv,repmat(s,[size(sper,1),1]));
    sres = cat(1,sres,round(mean(sper.data,2)));
end

stateCountH = cell2mat(cf(@(s)  size(s,1),  stch(states)))
stateCountN = cell2mat(cf(@(s)  size(s,1),  stcw(states)))


fragCnt = {};
stateFragMismatch = zeros([numel(states),numel(states)]);
for s = 1:numel(states),
    sper = stch{states{s},119.881035};
    for p = 1:size(sper,1),
        sind = WithinRanges(sres,sper(p,:));
        fragCnt{s}(p) = sum(sind);
        bhvs = sbhv(sind);
        stateFragMismatch(:,s) = stateFragMismatch(:,s)+histcounts(bhvs,1:numel(states)+1)';
    end
end


figure,hist(fragCnt{1},10)


sts = 'n';
figure,
subplot(311);hist(log10(diff(stch{sts}.data,1,2)),100)
subplot(312);hist(log10(diff(stcw{sts}.data,1,2)),100)
subplot(313);hist(log10(diff(stcp{sts}.data,1,2)),100)
af(@(ax) xlim(ax,[0,3.5]), findobj(gcf,'Type','axes'));
% CCG of state centers



