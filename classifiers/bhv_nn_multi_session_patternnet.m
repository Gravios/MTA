function [varargout] = bhv_nn_multi_session_patternnet(sessionList,varargin)
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

MODEL_TYPE = 'NNmultiPN';



Trials = af(@(Trial)  MTATrial.validate(Trial), get_session_list(sessionList));
xyz    = cf(@(Trial)  Trial.load('xyz')       , Trials);


varargout = cell([1,nargout-1]);

% Default Arguments for uninitiated variables
defArgs = struct('featureSet',                  'fet_bref_rev7',                               ...
                 'sampleRate',                  10,                                            ...
                 'states',                      {{'walk','rear','turn','pause','groom','sit'}},...
                 'model',                       [],                                            ...
                 'nNeurons',                    25,                                            ...
                 'nIter',                       10,                                            ...
                 'randomizationMethod',         'WSBNT',                                       ...
                 'map2reference',               true,                                          ... 
                 'normalize',                   true,                                          ...
                 'tag',                         '',                                            ...
                 'prctTrain',                   [],                                            ...
                 'normalizationSessionList',    'hand_labeled',                                ...
                 'referenceTrial',              'jg05-20120317',                               ...
                 'dropIndex',                   [],                                            ...
                 'trainModel',                  false                                          ...
);


[featureSet,sampleRate,states,model,nNeurons, nIter,randomizationMethod,...
 map2reference,normalize,tag,prctTrain,normalizationSessionList,referenceTrial,dropIndex] = ...
    DefaultArgs(varargin,defArgs,'--struct');

Trials(dropIndex) = [];
xyz(dropIndex) = [];


% CREATE model name from parameters if none is provided
if isempty(model),
    if ~isempty(prctTrain), prctTrainTag = ['PRT',prctTrain]; else, prctTrainTag = ''; end    
    model = ['MTAC_BATCH-' tag sessionList '+'         ...
             featureSet                                ...
             '+SR'  num2str(sampleRate)                ...
             'N'    num2str(normalize)                 ...             
             'REF+' referenceTrial                     ...
             '+NN'  num2str(nNeurons)                  ...
             'NI'   num2str(nIter)                     ...
             prctTrainTag                              ...
             '-'    MODEL_TYPE];
end
modelPath = fullfile(fileparts(mfilename('fullpath')),model);



%numTrials  = repmat({numel(Trials)},1,numel(Trials));
numTrials  = numel(Trials);
numStates  = repmat({numel(states)},1,numel(Trials));
numIter  = repmat({nIter},1,numel(Trials));
states     = repmat({states}       ,1,numel(Trials));
sampleRate = repmat({sampleRate}   ,1,numel(Trials));

% LOAD the feature set
features = cf(@(Trial,fetSet,sr) feval(fetSet,Trial,sr,false),...
               Trials,...
              repmat({featureSet},1,numel(Trials)),...
              sampleRate);

% LOAD the state collections
StcHL = cf(@(Trial) Trial.load('stc'),Trials);
        cf(@(stc,states) set(stc,'states',stc(states{:})),...
             StcHL,states);
        
%keys = cellfun(@subsref,StcHL(states{:}),repmat({substruct('.','key')},size(states)));


% MAP features to reference session
if map2reference,
    refTrial = MTATrial.validate(referenceTrial);
    

    % Map features via linear or circular shift to the training Session
    % of the Neural Network.
    cf(@(f,t,r) f.map_to_reference_session(t,r),...
       features, Trials, repmat({refTrial},1,numel(Trials)));

    for s = 1:numel(Trials), features{s}.data(~nniz(xyz{s}),:,:) = 0;end


    if normalize,
        [refMean,refStd] = load_normalization_parameters_unity(featureSet,...
                                                          refTrial.filebase,...
                                                          normalizationSessionList);
        cf(@(f,m,s) f.unity(@nan,m,s), ...
           features,...
           repmat({refMean},1,numel(Trials)),...
           repmat({refStd},1,numel(Trials)));
    end
elseif normalize,
    zfrCat = cf(@(f) get(f,'data'),features);
    zfrCat = cat(1,zfrCat{:});
    zfrMean = nanmean(zfrCat(nniz(zfrCat),:,:));
    zfrStd = nanstd(zfrCat(nniz(zfrCat),:,:));
    cf(@(w,m,s) set(w,'data',nunity(w,[],m,s)),...
       features,...
       repmat({zfrMean},1,numel(Trials)),...
       repmat({zfrStd},1,numel(Trials)));
    clear('zfrCat','zfrMean','zfrStd')
end



% Use xyz as base for labeling output ... quit staring at this part.

% INITIALIZE outputs
p_state = cf(@(x,n)  zeros([x.size(1),n]),    xyz, numStates);
d_state = cf(@(x,n)  zeros([x.size(1),n]),    xyz, numStates);
labelingStatsMulti.confussionMatrix = cf(@(i,n) zeros([i,n,n]),numIter,numStates);
labelingStatsMulti.precision        = cf(@(i,n) zeros([i,n])  ,numIter,numStates);
labelingStatsMulti.sensitivity      = cf(@(i,n) zeros([i,n])  ,numIter,numStates);
labelingStatsMulti.accuracy         = cf(@(i,n) zeros([i,1])  ,numIter);


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
        % 00:30 FTS
        if iter==1&&trainModel
            model = [model '_RAND_' randomizationMethod];
        end
        if iter==1,
            model_out = model;
            model_path = fileparts(mfilename('fullpath'));
            mkdir(fullfile(model_path,model));
            model = [model '/' model];
        end


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
            
            % Concatenate stc's and features
            %StcRnd = 
            

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
            %net.trainParam.showWindow = true;
            net.trainParam.showWindow = false;
            net.inputs{1}.processFcns(2) = [];            
            [net,tr] = train(net,mapminmax('apply',trainingFeatures(ind,:)',psa),~~smat(ind,:)');
            save(fullfile(modelPath,[modelPath,'-',num2str(iter),'.mat']),'net','tr','Model_Information');
        end

% LOAD network        
        load(fullfile(modelPath,[modelPath,'-',num2str(iter),'.mat']));

        d_state = cf(@(f,p)  net(mapminmax('apply',f.data',p))',  features,repmat({psa},[1,numTrials]));
        d_state = cf(@(d,s)  MTADxyz('data',d,'sampleRate',s),    d_state,sampleRate);
                  cf(@(d,x)  d.resample(x),                       d_state,xyz);

        d_state = cf(@(d) d.data, d_state);;
        n_state = d_state;

%  the winners from the losers
        [~,maxState] = cf(@(d)    max(d,[],2),   d_state);
        %maxState     = cf(@(d,s)  MTADxyz('data',d,'sampleRate',s),    maxState,sampleRate);        
        %               cf(@(m,f)  set(m,'data',m(~nniz(f),:) = 0;
        for s = 1:numTrials, maxState{s}(~nniz(features{s}),:) = 0; end
        

% Smooth decision boundaries - 200 ms state minimum
% $$$         bwin = round(.2*sampleRate)+double(mod(round(.2*sampleRate),2)==0);
% $$$         mss = GetSegs(maxState,1:size(maxState,1),bwin,nan);
% $$$         maxState=circshift(sq(mode(mss))',floor(bwin/2));

% $$$ % Populate Stc object with the new states
% $$$ for i = 1:numel(Model_Information.state_labels),
% $$$ 
% $$$     sts = ThreshCross(maxState==i,0.5,1);
% $$$     if ~isempty(sts),
% $$$         sts = bsxfun(@plus,sts,[1,0]);
% $$$     end
% $$$     
% $$$     Stc.addState(Trial.spath,...
% $$$              Trial.filebase,...
% $$$              sts,...
% $$$              xyz.sampleRate,...
% $$$              feature.sync.copy,...
% $$$              feature.origin,...
% $$$              Model_Information.state_labels{i},...
% $$$              Model_Information.state_keys{i},...
% $$$              'TimePeriods');
% $$$ end
        
        
        
        % Label States
        [Stc,ps,Model_Information,ns] = bhv_nn (Trial,         ... Trial
                                             false,         ... ifTrain
                                             trainingStates,... States
                                             Trial.stc.copy,... StateCollection
                                             features,      ... feature_set
                                             [model '_' num2str(iter)]); % model name

        % if an stc was provided get comparison stats
        ysm = MTADxyz('data',double(0<stc2mat(Stc,xyz)),'sampleRate',xyz.sampleRate); 
        d_state = nansum(cat(3,ysm.data,d_state),3);
        p_state = p_state +ns;
        
        if ~isempty(stcMode)            
            if nargout>=4,
                labelingEpochs.resample(xyz);
    
                ind = any(shl.data,2)&any(ysm.data,2)&labelingEpochs.data;
    
                tcm = confmat(shl(ind&labelingEpochs,:),ysm(ind&labelingEpochs,:)); % DEP: netlab
                labelingStatsMulti.confusionMatrix(iter,:,:) = round(tcm./xyz.sampleRate,2);
                labelingStatsMulti.precision(iter,:) = round(diag(tcm)./sum(tcm,2),4).*100;
                labelingStatsMulti.sensitivity(iter,:) = round(diag(tcm)'./sum(tcm),4).*100;
                labelingStatsMulti.accuracy(iter) = sum(diag(tcm))/sum(tcm(:));
            end
        end
        
    catch err,
        keyboard
        warning([err.identifier,', ',err.message]);
        continue
    end
    
end

if nargout==0,
    return,
end


% $$$ 
d_state = MTADxyz('data',p_state,'sampleRate',xyz.sampleRate);
% $$$ 
% $$$ 
% $$$ figure,
% $$$ sp    = subplot(311);imagesc(d_state.data');caxis([20,100]);
% $$$ fds = d_state.copy;fds.filter('ButFilter',5,1,'low');
% $$$ sp(2) = subplot(312);imagesc(fds.data');caxis([20,100]);
% $$$ sp(3) = subplot(313);imagesc(shl.data');caxis([0,1]);
% $$$ linkaxes(sp,'xy')
% Determine winning states based on the the labels of nurmerous
% neural networks.
d_state.filter('ButFilter',3,2,'low');

[~,maxState] = max(d_state.data,[],2);
maxState(~nniz(xyz)) = 0;
maxState(~any(d_state.data,2)) = 0;
% Smooth decision boundaries - 200 ms state minimum
% $$$ bwin = round(.2*xyz.sampleRate)+double(mod(round(.2*xyz.sampleRate),2)==0);
% $$$ mss = GetSegs(maxState,1:size(maxState,1),bwin,nan);
% $$$ maxState=circshift(sq(mode(mss))',floor(bwin/2));


% Stats in comparision to the collection of labels specified in the stcMode
labelingStats = struct();
if ~isempty(stcMode),
    ysm = MTADxyz('data',zeros([shl.size]),'sampleRate',xyz.sampleRate); 
    ysm.data = ysm.data';
    ysm.data([1:size(ysm,1):size(ysm,2).*size(ysm,1)]+maxState'-1) = 1;
    ysm.data = ysm.data';
    ind = any(shl.data,2)&any(ysm.data,2)&labelingEpochs.data;
    tcm = confmat(shl(ind&labelingEpochs,:),ysm(ind&labelingEpochs,:)); % #DEP: netlab
    labelingStats.confusionMatrix = round(tcm./xyz.sampleRate,2);
    labelingStats.precision = round(diag(tcm)./sum(tcm,2),4).*100;
    labelingStats.sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
    labelingStats.accuracy = sum(diag(tcm))/sum(tcm(:));
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
if nargout>=2, varargout{2} = d_state.data;            end
if nargout>=3, varargout{3} = labelingStats;      end
if nargout>=4, varargout{4} = labelingStatsMulti; end
if nargout>=5, varargout{5} = model_out; end
if nargout>=6, varargout{6} = p_state; end