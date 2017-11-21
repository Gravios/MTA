function [varargout] = bhv_nn_multi_session_weighted_patternnet(varargin)
%function [stc,varargout] = bhv_nn_multi_patternnet(Trial,varargin)
% 
%
% varargin:
%    states,
%    stcMode,
%    featureSet,
%    sampleRate,
%    model,
%    nNeurons,
%    nIter,
%    randomizationMode,
%
% varargout:
%     Stc
%     labelingStats
%     cumulativeNetworkOutput
%     model

%    Stc,
%    d_state,
%    labelingStats,
%    labelingStatsMulti,

global MTA_PROJECT_PATH
%global MTA_VERBOSITY

MODEL_TYPE = 'multiSesPatNet';

varargout = cell([1,nargout-1]);

%% DEFARGS ----------------------------------------------------------------------------------------
defArgs = struct('sessionList',                 '',                                            ...
                 'featureSet',                  'fet_bref_rev7',                               ...
                 'states',                      {{'walk','rear','turn','pause','groom','sit'}},...
                 'keys',                        {{'w','r','n','p','m','s'}},                   ...
                 'model',                       '',                                            ...
                 'sampleRate',                  10,                                            ...
                 'nNeurons',                    25,                                            ...
                 'nIter',                       100,                                           ...
                 'randomizationMethod',         'WSBNT',                                       ...
                 'map2reference',               true,                                          ...
                 'normalize',                   true,                                          ...
                 'referenceTrial',              'jg05-20120317.cof.all',                       ...
                 'trainingSessionList',         'hand_labeled',                                ...
                 'normalizationSessionList',    'hand_labeled',                                ...
                 'dropIndex',                   [],                                            ...
                 'prctTrain',                   90                                             ...
);
[sessionList,featureSet,states,keys,model,sampleRate,nNeurons,nIter,randomizationMethod,       ...
 map2reference,normalize,referenceTrial,trainingSessionList,normalizationSessionList,          ...
 dropIndex,prctTrain] = DefaultArgs(varargin,defArgs,'--struct');
% ------------------------------------------------------------------------------------------------



%% MAIN -------------------------------------------------------------------------------------------

trainModel = false;
if isempty(sessionList),
% LOAD Trials for training neural network    
    Trials = af(@(Trial)  MTATrial.validate(Trial), get_session_list(trainingSessionList));
    trainModel = true;
elseif iscell(sessionList) & isa(sessionList{1},'MTASession'),
% ASSIGN cell array of MTATrials for behavioral labeling    
    Trials = sessionList;
elseif isa(sessionList,'MTATrial'),
    Trials = {sessionList};
else
% LOAD Trials for behavioral labeling
    Trials = af(@(Trial)  MTATrial.validate(Trial), get_session_list(sessionList));

end

% LOAD xyz as a reference for synchronization
xyz    = cf(@(Trial)  Trial.load('xyz')       , Trials);

% REMOVE unwanted Trials from list
if isempty(sessionList),
    Trials(dropIndex) = [];
    xyz(dropIndex) = [];
end

% LOAD the state collections
StcHL = cf(@(Trial) Trial.load('stc'),Trials);
if ~any(cell2mat(cf(@(s) strcmp(s.mode,'default'),StcHL))),
    cf(@(stc,states) set(stc,'states',stc(states{:})),...
         StcHL,repmat({states},[1,numel(Trials)]));
end
        
% CREATE model name from parameters if no model name is provided
if isempty(model),
    dropIndexTag = '';
    if ~isempty(dropIndex), dropIndexTag = num2str(dropIndex);
        dropIndexTag(dropIndexTag==[' '])='_';
    end
    model = ['MTAC_BATCH+' trainingSessionList '+'                        ...
             dropIndexTag '+'                                             ...
             featureSet                                                   ...
             '+SR'   num2str(sampleRate)                                  ...
             'NN'    num2str(nNeurons)                                    ...
             'NI'    num2str(nIter)                                       ...
             'M'     num2str(map2reference)                               ...
             'MREF+' referenceTrial                                       ...                          
             '+N'    num2str(normalize)                                   ...             
             'NREF+' normalizationSessionList                             ...
             '+RND'  randomizationMethod                                  ...
             '+PRCT' num2str(prctTrain)                                   ...
             '+STS+' strjoin(keys,'')                                     ...
             '+'     MODEL_TYPE];
end

% CREATE model path from model name 
modelPath = fullfile(fileparts(mfilename('fullpath')),model);
if ~exist(modelPath,'dir'),mkdir(modelPath);end
disp(model);


% ENCAPSULATE basic vars in cell arrays, replicated for cellfun operations
numTrials  = numel(Trials);
numStates  = repmat({numel(states)},[1,numTrials]);
numIter    = repmat({nIter}        ,[1,numTrials]);
states     = repmat({states}       ,[1,numTrials]);
sampleRate = repmat({sampleRate}   ,[1,numTrials]);


% LOAD the feature set used for mapping behavior
features = cf(@(Trial,fetSet,sr) feval(fetSet,Trial,sr,false),...
               Trials,...
              repmat({featureSet},[1,numTrials]),...
              sampleRate);


% MAP features to reference session
if map2reference,
    xyzls  = cf(@(Trial)  Trial.load('xyz'),        Trials);
             cf(@(f,x)    x.resample(f),            features,xyzls);
    refTrial = MTATrial.validate(referenceTrial);
    cf(@(f,t,r) f.map_to_reference_session(t,r),...
       features, Trials, repmat({refTrial},[1,numTrials]));
    for s = 1:numTrials, features{s}.data(~nniz(xyzls{s}),:,:) = 0;end

    if normalize,
% CONDITIONAL normalization, use multiple sessions for normalization        
        [refMean,refStd] = load_normalization_parameters_unity(featureSet,...
                                                          refTrial.filebase,...
                                                          normalizationSessionList);
        cf(@(f,m,s) f.unity(@nan,m,s), ...
           features,...
           repmat({refMean},[1,numTrials]),...
           repmat({refStd}, [1,numTrials]));
    end

elseif normalize,
% CONDITIONAL normalization, use multiple sessions for normalization    
    zfrCat = cf(@(f) get(f,'data'),    features); 
    zfrCat = cat(1,zfrCat{:});
    zfrMean = nanmean(zfrCat(nniz(zfrCat),:,:));
    zfrStd  = nanstd( zfrCat(nniz(zfrCat),:,:));
    cf(@(w,m,s) set(w,'data',nunity(w,[],m,s)),...
       features,...
       repmat({zfrMean},[1,numTrials]),...
       repmat({zfrStd}, [1,numTrials]));
    clear('zfrCat','zfrMean','zfrStd');
end


% INITIALIZE outputs
p_state = cf(@(x,n)  zeros([x.size(1),n]),    xyz, numStates);
d_state = cf(@(x,n)  zeros([x.size(1),n]),    xyz, numStates);

% INITIALIZE network output collection variable
cumulativeNetworkOutput = cf(@(f) f.copy('empty'),    features);


% INITIALIZE labeling timeperiods
labelingEpochs = cf(@(t,x)  resample(t.stc{'a'}.cast('TimeSeries'),x), Trials,xyz);


% LOAD mapminmax parameters
psa = load_normalization_parameters_mapminmax(features{1}.label,...
                                              [],... need feature.treatmentRecord
                                              'hand_labeled');
% LOAD hand labeled state matrix
shl = {};
if ~any(cell2mat(cf(@(s) strcmp(s.mode,'default'),StcHL))),
    shl = cf(@(s,x,sts) MTADxyz('data',double(0<stc2mat(s,x,sts)),...
                                'sampleRate',x.sampleRate),...
             StcHL, xyz, states);
end


labelingStatsMulti = repmat(struct(...
    'precision',     zeros([nIter,numStates{1}]),...
    'sensitivity',   zeros([nIter,numStates{1}]),...
    'accuracy',      zeros([nIter,1]))   ,[numTrials,1]);

% TRAIN or LABEL Trials
for iter = 1:nIter,
    try,        
% $$$         if trainModel,
% $$$             disp([mfilename,': training iteration: ',num2str(iter)]);
% $$$ % RESAMPLE feature matrix
% $$$             switch randomizationMethod
% $$$               case 'ERS' % equal_restructured_sampling
% $$$                 % Can be found in MTA:classifiers:bhv_nn_multi_patternnet.m
% $$$               case 'WSB'   % whole state bootstrap
% $$$                 [StcRnd,~,trainingFeatures] = ...
% $$$                     cf(@(s,f,sts) resample_whole_state_bootstrap(s,f,sts),...
% $$$                        StcHL,features,states);
% $$$                 trainingEpochs = [];
% $$$               case 'WSBN'  % whole state bootstrap noisy
% $$$                 [StcRnd,~,trainingFeatures] = ...
% $$$                     cf(@(s,f,sts) resample_whole_state_bootstrap_noisy(s,f,sts),...
% $$$                        StcHL,features,states);
% $$$                 trainingEpochs = [];
% $$$               case 'WSBNT' % whole state bootstrap noisy with trimmed boundaries
% $$$                 [StcRnd,~,trainingFeatures] = ...
% $$$                     cf(@(s,f,sts) resample_whole_state_bootstrap_noisy_trim(s,f,sts),...
% $$$                        StcHL,features,states);
% $$$                 trainingEpochs = [];
% $$$               case 'WSBT' % whole state bootstrap noisy with trimmed boundaries
% $$$                 [StcRnd,~,trainingFeatures] = ...
% $$$                     cf(@(s,f,sts) resample_whole_state_bootstrap_trim(s,f,sts),...
% $$$                        StcHL,features,states);
% $$$                 trainingEpochs = [];
% $$$               case 'rndsamp'
% $$$                 % Can be found in MTA:classifiers:bhv_nn_multi_patternnet.m                
% $$$             end
% $$$             
% $$$ 
% $$$ % CAST state collection object into state label matrix
% $$$             [smat] = cf(@(stc,fet,states) stc2mat(stc,fet,states), ...
% $$$                         StcRnd,trainingFeatures,states);
% $$$             smat = cat(1,smat{:});
% $$$             
% $$$ % COCATENATE feature matrix
% $$$             trainingFeatures = cf(@(f) f.data, trainingFeatures);
% $$$             trainingFeatures = cat(1,trainingFeatures{:});
% $$$ 
% $$$ % SELECT indecies which are not zero, nan or inf 
% $$$             ind    = any(smat,2)&nniz(trainingFeatures);
% $$$ 
% $$$ % CREATE network object for training
% $$$             net = patternnet(nNeurons);
% $$$             %net.trainParam.showWindow = true;
% $$$             net.trainParam.showWindow = false;
% $$$             net.inputs{1}.processFcns(2) = [];            
% $$$             [net,tr] = train(net,mapminmax('apply',trainingFeatures(ind,:)',psa),~~smat(ind,:)');
% $$$             save(fullfile(modelPath,[model,'-',num2str(iter),'.mat']),'net','tr');
% $$$ 
% $$$             disp([mfilename,': training iteration: ',num2str(iter),': complete']);
% $$$             
% $$$         end

% LOAD network        
        load(fullfile(modelPath,[model,'-',num2str(iter),'.mat']));

        networkOutput = cf(@(f,p)  net(mapminmax('apply',f.data',p))',  features,repmat({psa},[1,numTrials]));
        networkOutput = cf(@(d,s)  MTADxyz('data',d,'sampleRate',s),    networkOutput,sampleRate);
                        cf(@(d,x)  d.resample(x),                       networkOutput,xyz);
% ACCUMULATE network output
        cf(@(c,n) set(c,'data',cat(3,c.data,n.data)),   cumulativeNetworkOutput,networkOutput);

% COMPUTE labeling statistics of individual patternnets
        if ~isempty(shl),
% SELECT maximum network output as state, winner takes all
            [~,maxState] = cf(@(c)    max(mean(c.data,3),[],2),   networkOutput);
% UNLABEL regions where labels cannot exist 
            for s = 1:numTrials, maxState{s}(~nniz(xyz{s}),:) = 0; end
% CREATE state collection with simple mode name
            Stc = cf(@(s) s.copy,StcHL);
            cf(@(s) set(s,'states',{}), Stc);
% APPEND network output labeles to state collection
            for i = 1:numStates{1},
                sts = cf(@(m,i) ThreshCross(m==i,0.5,1),    maxState,repmat({i},[1,numTrials]));
                try
                    sts = cf(@(m) bsxfun(@plus,m,[1,0]),    sts);
                end    
                cf(@(sts,s,t,x,label,key) s.addState(t.spath,t.filebase,sts,x.sampleRate,...
                                   x.sync.copy,x.origin,label,key,'TimePeriods'),...
                   sts,Stc,Trials,xyz,repmat(states{1}(i),[1,numTrials]),repmat(keys(i),[1,numTrials]))
            end
            
            ysm = cf(@(s,x) MTADxyz('data',double(0<stc2mat(s,x)),'sampleRate',x.sampleRate),...
                     Stc,xyz);
            for s = 1:numTrials, ysm{s}.data(~nniz(xyz{s}),:) = 0;  end

            ind = cf(@(s,y,l) any(s.data,2)&any(y.data,2)&l.data==1, shl,ysm,labelingEpochs);
            tcm = cf(@(s,y,i) confmat(s(i,:),y(i,:)),      shl,ysm,ind); 
            % #DEP: netlab
            for s = 1:numTrials,
                labelingStatsMulti(s).precision(iter,:) = round(diag(tcm{s})'./sum(tcm{s},2)',4);
                labelingStatsMulti(s).sensitivity(iter,:) = round(diag(tcm{s})'./sum(tcm{s}),4);
                labelingStatsMulti(s).accuracy(iter,:) = sum(diag(tcm{s}))/sum(tcm{s}(:));
            end
        end
        
    catch err,
        disp(err);
        for e = 1:numel(err.stack), disp(err.stack(e)); end        
        keyboard
    end
    
end



% COMPUTE scores for each state 
if trainModel || ~exist(fullfile(modelPath,[model,'-weights.mat']),'file'),
    ranks = [1:nIter]./sum(1:nIter);
    for i = 1:numStates{1},
        stateScore = af(@(l,i)  l.sensitivity(:,i)+l.precision(:,i),  ...
                        labelingStatsMulti,repmat([i],[numTrials,1]));
        stateScore = cat(2,stateScore{:});
        [~,stateRankInd] = sort(mean(stateScore,2));
        stateRanks(:,i) = ranks(stateRankInd);
    end
    save(fullfile(modelPath,[model,'-weights.mat']),'stateRanks');
else
    load(fullfile(modelPath,[model,'-weights.mat']));
end


if nargout==0,
    return,
end


% FILTER network output to reduce 
%cf(@(c) c.filter('ButFilter',3,1,'low'),cumulativeNetworkOutput);
disp([mfilename,': labeling states']);

%[labelingStatsMulti(:,:).accuracy]

% SELECT maximum network output as state, winner takes all
cf(@(c,w)  set(c,'data',bsxfun(@times,c.data,permute(w,[3,2,1]))),...
                         cumulativeNetworkOutput,repmat({stateRanks},[1,numTrials]));
mcno = cf(@(c)    sum(c.data,3),   cumulativeNetworkOutput);
[~,maxState] = cf(@(c)    max(sum(c,3),[],2),   mcno);
% UNLABEL regions where labels cannot exist 
for s = 1:numTrials, maxState{s}(~nniz(xyz{s}),:) = 0; end


% CREATE state collection with simple mode name
Stc = cf(@(s) s.copy,StcHL);
cf(@(s) s.updateMode(['mswnnN' num2str(trainModel) '+' trainingSessionList]), Stc);
cf(@(s) set(s,'states',{}), Stc);


% APPEND network output labeles to state collection
for i = 1:numStates{1},
    sts = cf(@(m,i) ThreshCross(m==i,0.5,1),    maxState,repmat({i},[1,numTrials]));
    try
        sts = cf(@(m) bsxfun(@plus,m,[1,0]),    sts);
    end    
    cf(@(sts,s,t,x,label,key) s.addState(t.spath,t.filebase,sts,x.sampleRate,...
                       x.sync.copy,x.origin,label,key,'TimePeriods'),...
       sts,Stc,Trials,xyz,repmat(states{1}(i),[1,numTrials]),repmat(keys(i),[1,numTrials]))
end

% SAVE state collection
cf(@(s) s.save(true),    Stc);
% UPDATE state collection mode with model name
cf(@(s,m) s.updateMode(m),Stc,repmat({[model,'_weighted']},[1,numTrials]));
% SAVE state collection
cf(@(s) s.save(true),    Stc);


% if MTA_VERBOSITY > 1, 
disp([mfilename,': labeling states: complete']);
disp([mfilename,': computing labeling statistics']);
% end


% CONVERT state collection into state matrix for inter labeler comparison
ysm = cf(@(s,x) MTADxyz('data',double(0<stc2mat(s,x)),'sampleRate',x.sampleRate),...
         Stc,xyz);
for s = 1:numTrials, 
    ysm{s}.data(~nniz(xyz{s}),:) = 0; 
end

% INITIALIZE inter labeler stats
labelingStats = repmat(struct('confusionMatrix',zeros([numStates{1},numStates{1}]),...
                               'precision',     zeros([1,numStates{1}]),...
                               'sensitivity',   zeros([1,numStates{1}]),...
                               'accuracy',      0 ),[1,numTrials]);


% COMPUTE inter labeler stats if valid comparative Stc was loaded above
if ~isempty(shl),
    ind = cf(@(s,y,l) any(s.data,2)&any(y.data,2)&l.data==1, shl,ysm,labelingEpochs);
    tcm = cf(@(s,y,i) confmat(s(i,:),y(i,:)),      shl,ysm,ind); % #DEP: netlab
    af(@(l,t,s) setfield(l,'confusionMatrix',round(t{1}./s{1},2)),labelingStats,tcm,sampleRate);
    for s = 1:numTrials,
        labelingStats(s).confusionMatrix = round(tcm{s}./sampleRate{s},2);
        labelingStats(s).precision = round(diag(tcm{s})'./sum(tcm{s},2)',4).*100;
        labelingStats(s).sensitivity = round(diag(tcm{s})'./sum(tcm{s}),4).*100;
        labelingStats(s).accuracy = sum(diag(tcm{s}))/sum(tcm{s}(:));
    end
end

disp([mfilename,': computing labeling statistics: complete']);

disp([mfilename,': complete']);

if nargout>=1, varargout{1} = Stc;                end
if nargout>=2, varargout{2} = labelingStats;      end
if nargout>=3, varargout{3} = cumulativeNetworkOutput; end
if nargout>=4, varargout{4} = model; end
