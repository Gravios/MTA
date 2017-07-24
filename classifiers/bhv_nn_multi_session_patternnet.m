function [varargout] = bhv_nn_multi_session_patternnet(varargin)
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
%    Stc,
%    d_state,
%    labelingStats,
%    labelingStatsMulti,

global MTA_PROJECT_PATH

MODEL_TYPE = 'multiSesPatNet';

varargout = cell([1,nargout-1]);

% Default Arguments for uninitiated variables
defArgs = struct('sessionList',                 '',                                            ...
                 'featureSet',                  'fet_bref_rev7',                               ...
                 'sampleRate',                  10,                                            ...
                 'states',                      {{'walk','rear','turn','pause','groom','sit'}},...
                 'model',                       [],                                            ...
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
                 'stcTag',                      'msnn'                                         ...
);


[sessionList,featureSet,sampleRate,states,model,nNeurons,nIter,randomizationMethod,...
 map2reference,normalize,referenceTrial,trainingSessionList,normalizationSessionList,...
 dropIndex,prctTrain,stcTag] = DefaultArgs(varargin,defArgs,'--struct');

if isempty(sessionList),
    Trials = af(@(Trial)  MTATrial.validate(Trial), get_session_list(trainingSessionList));    
    trainModel = true;
else
    Trials = af(@(Trial)  MTATrial.validate(Trial), get_session_list(sessionList));
    trainModel = false;
end

xyz    = cf(@(Trial)  Trial.load('xyz')       , Trials);



Trials(dropIndex) = [];
xyz(dropIndex) = [];

% CREATE model name and directory from parameters if none is provided
if isempty(model),
    dropTag = '';
    if ~isempty(dropIndex), dropTag = num2str(dropIndex);
        dropTag(dropTag==[' '])='_';
    end
    model = ['MTAC_BATCH-' trainingSessionList '-'            ...
             dropTag '+'                                      ...
             featureSet                                       ...
             '+SR'  num2str(sampleRate)                       ...
             'NN'  num2str(nNeurons)                         ...
             'NI'   num2str(nIter)                            ...
             'M'    num2str(map2reference)                    ...
             'MREF+' referenceTrial                          ...                          
             '+N'    num2str(normalize)                        ...             
             'NREF+' referenceTrial                          ...
             '+RND' randomizationMethod                        ...
             'PRCT' num2str(prctTrain)                        ...
             '-'    MODEL_TYPE];
end
modelPath = fullfile(fileparts(mfilename('fullpath')),model);
if ~exist(modelPath,'dir'),mkdir(modelPath);end

numTrials  = numel(Trials);
numStates  = repmat({numel(states)},[1,numTrials]);
numIter    = repmat({nIter}        ,[1,numTrials]);
states     = repmat({states}       ,[1,numTrials]);
sampleRate = repmat({sampleRate}   ,[1,numTrials]);

% LOAD the feature set
features = cf(@(Trial,fetSet,sr) feval(fetSet,Trial,sr,false),...
               Trials,...
              repmat({featureSet},[1,numTrials]),...
              sampleRate);

% LOAD the state collections
StcHL = cf(@(Trial) Trial.load('stc'),Trials);
        cf(@(stc,states) set(stc,'states',stc(states{:})),...
             StcHL,states);
keys   = cf(@(s) s.key,  StcHL{1}(states{1}));
labels = cf(@(s) s.label,StcHL{1}(states{1}));

% MAP features to reference session
if map2reference,
    refTrial = MTATrial.validate(referenceTrial);
    cf(@(f,t,r) f.map_to_reference_session(t,r),...
       features, Trials, repmat({refTrial},[1,numTrials]));
    for s = 1:numTrials, features{s}.data(~nniz(xyz{s}),:,:) = 0;end

% CONDITIONAL normalization        
    if normalize,
        [refMean,refStd] = load_normalization_parameters_unity(featureSet,...
                                                          refTrial.filebase,...
                                                          normalizationSessionList);
        cf(@(f,m,s) f.unity(@nan,m,s), ...
           features,...
           repmat({refMean},[1,numTrials]),...
           repmat({refStd}, [1,numTrials]));
    end
% CONDITIONAL normalization    
elseif normalize,
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
labelingStatsMulti.confussionMatrix = cf(@(i,n) zeros([i,n,n]),numIter,numStates);
labelingStatsMulti.precision        = cf(@(i,n) zeros([i,n])  ,numIter,numStates);
labelingStatsMulti.sensitivity      = cf(@(i,n) zeros([i,n])  ,numIter,numStates);
labelingStatsMulti.accuracy         = cf(@(i,n) zeros([i,1])  ,numIter);

% INITIALIZE network output collection variable
cumulativeNetworkOutput = cf(@(f) f.copy('empty'),    features);

% INITIALIZE labeling timeperiods
labelingEpochs = cf(@(Trial)    Trial.stc{'a'}.cast('TimeSeries'),    Trials);

% LOAD mapminmax parameters
psa = load_normalization_parameters_mapminmax(features{1}.label,...
                                              [],... need feature.treatmentRecord
                                              'hand_labeled');
% LOAD hand labeled state matrix
shl = cf(@(s,x,sts) MTADxyz('data',double(0<stc2mat(s,x,sts)),...
                            'sampleRate',x.sampleRate),...
         StcHL, xyz, states);



for iter = 1:nIter,
    try,        
        if trainModel,
            
            switch randomizationMethod
              case 'ERS' % equal_restructured_sampling
                % Can be found in MTA:classifiers:bhv_nn_multi_patternnet.m
              case 'WSB'   % whole state bootstrap
                [StcRnd,labelingEpochs,trainingFeatures] = ...
                    cf(@(s,f,sts) resample_whole_state_bootstrap(s,f,sts),...
                       StcHL,features,states);
                trainingEpochs = [];
              case 'WSBN'  % whole state bootstrap noisy
                [StcRnd,labelingEpochs,trainingFeatures] = ...
                    cf(@(s,f,sts) resample_whole_state_bootstrap_noisy(s,f,sts),...
                       StcHL,features,states);
                trainingEpochs = [];
              case 'WSBNT' % whole state bootstrap noisy with trimmed boundaries
                [StcRnd,labelingEpochs,trainingFeatures] = ...
                    cf(@(s,f,sts) resample_whole_state_bootstrap_noisy_trim(s,f,sts),...
                       StcHL,features,states);
                trainingEpochs = [];
              case 'rndsamp'
                % Can be found in MTA:classifiers:bhv_nn_multi_patternnet.m                
            end
            

% CAST Stc object into state timeseries matrix
            [smat] = cf(@(stc,fet,states) stc2mat(stc,fet,states), ...
                        StcRnd,trainingFeatures,states);
            smat = cat(1,smat{:});
% COCATENATE feature matrix
            trainingFeatures = cf(@(f) f.data, trainingFeatures);
            trainingFeatures = cat(1,trainingFeatures{:});
% SELECT indecies which are not zero, nan or inf 
            ind    = any(smat,2)&nniz(trainingFeatures);

% CREATE network object for training
            net = patternnet(nNeurons);
            net.trainParam.showWindow = true;
            %net.trainParam.showWindow = false;
            net.inputs{1}.processFcns(2) = [];            
            [net,tr] = train(net,mapminmax('apply',trainingFeatures(ind,:)',psa),~~smat(ind,:)');
            save(fullfile(modelPath,[model,'-',num2str(iter),'.mat']),'net','tr');

        end

% LOAD network        
        load(fullfile(modelPath,[model,'-',num2str(iter),'.mat']));

        networkOutput = cf(@(f,p)  net(mapminmax('apply',f.data',p))',  features,repmat({psa},[1,numTrials]));
        networkOutput = cf(@(d,s)  MTADxyz('data',d,'sampleRate',s),    networkOutput,sampleRate);
                        cf(@(d,x)  d.resample(x),                       networkOutput,xyz);
% ACCUMULATE network output
        cf(@(c,n) set(c,'data',cat(3,c.data,n.data)),   cumulativeNetworkOutput,networkOutput);
        
        
    catch err,
        disp(err);
        for e = 1:numel(err.stack), disp(err.stack(e)); end        
        keyboard
    end
    
end

if nargout==0,
    return,
end



%cf(@(c) c.filter('ButFilter',3,1,'low'),cumulativeNetworkOutput);



%  the winners from the losers
[~,maxState] = cf(@(c)    max(mean(c.data,3),[],2),   cumulativeNetworkOutput);
%maxState     = cf(@(d,s)  MTADxyz('data',d,'sampleRate',s.sampleRate),    maxState,xyz);
for s = 1:numTrials, 
    maxState{s}(~nniz(xyz{s}),:) = 0; 
end


Stc = cf(@(s) s.copy,StcHL);
cf(@(s) s.updateMode(['msnnN' num2str(trainModel) '+' trainingSessionList]), Stc);
cf(@(s) set(s,'states',{}), Stc);

for i = 1:numel(states),
    sts = cf(@(m,i) ThreshCross(m==i,0.5,1),    maxState,repmat({i},[1,numTrials]));
    try
        sts = cf(@(m) bsxfun(@plus,m,[1,0]),    sts);
    end

    
    cf(@(sts,s,t,x,label,key) s.addState(t.spath,t.filebase,sts,x.sampleRate,...
                       x.sync.copy,x.origin,label,key,'TimePeriods'),...
       sts,Stc,Trials,xyz,repmat(labels(i),[1,numTrials]),repmat(keys(i),[1,numTrials]))
end



        % if an stc was provided get comparison stats
ysm = cf(@(s,x) MTADxyz('data',double(0<stc2mat(s,x)),'sampleRate',x.sampleRate),...
         Stc,xyz);


% Stats in comparision to the collection of labels specified in the stcMode
labelingStats = repmat(struct('confusionMatrix',zeros([numStates{1},numStates{1}]),...
                               'precision',     zeros([1,numStates{1}]),...
                               'sensitivity',   zeros([1,numStates{1}]),...
                               'accuracy',      0 ),[1,numTrials]);

    ind = cf(@(s,y,l) any(s.data,2)&any(y.data,2), shl,ysm);
    tcm = cf(@(s,y,i) confmat(s(i,:),y(i,:)),      shl,ysm,ind); % #DEP: netlab
    af(@(l,t,s) setfield(l,'confusionMatrix',round(t{1}./s{1},2)),labelingStats,tcm,sampleRate);
for s = 1:numTrials,
    labelingStats(s).confusionMatrix = round(tcm{s}./sampleRate{s},2);
    labelingStats(s).precision = round(diag(tcm{s})'./sum(tcm{s},2)',4).*100;
    labelingStats(s).sensitivity = round(diag(tcm{s})'./sum(tcm{s}),4).*100;
    labelingStats(s).accuracy = sum(diag(tcm{s}))/sum(tcm{s}(:));
end

% Copy State Collection object to store new labeled periods
% $$$ Stc = Trial.stc.copy;
% $$$ Stc.updateMode([MODEL_TYPE '-' Model_Information.StcMode...
% $$$                 '-' cell2mat(Model_Information.state_keys)]);
% $$$ Stc.states = {};



% Create new StateCollection ... well copy
% $$$ Stc = Trial.stc.copy;
% $$$ Stc.updateMode([MODEL_TYPE '-' Model_Information.Trial '-'...
% $$$                 Model_Information.StcMode...
% $$$                 '-' cell2mat(Model_Information.state_keys)]);
% $$$ Stc.states = {};

% Create new StateCollection ... well copy
Stc = Trial.stc.copy;
Stc.updateMode([model_out '-' cell2mat(Model_Information.state_keys)]);
Stc.states = {};


% Populate Stc object with the new states
for i = 1:numel(Model_Information.state_labels),

    sts = ThreshCross(maxState==i,0.5,1);
    if ~isempty(sts),
        sts = bsxfun(@plus,sts,[1,0]);
    end
    
    Stc.addState(Trial.spath,...
             Trial.filebase,...
             sts,...
             xyz.sampleRate,...
             features.sync.copy,...
             features.origin,...
             Model_Information.state_labels{i},...
             Model_Information.state_keys{i},...
             'TimePeriods');
end

if nargout>=1, varargout{1} = Stc;                end
if nargout>=2, varargout{2} = d_state.data;       end
if nargout>=6, varargout{6} = p_state;            end
if nargout>=3, varargout{3} = labelingStats;      end
if nargout>=4, varargout{4} = labelingStatsMulti; end
if nargout>=5, varargout{5} = model_out;          end
