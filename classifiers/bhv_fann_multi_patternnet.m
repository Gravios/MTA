function [varargout] = bhv_fann_multi_patternnet(Trial,varargin)
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


MODEL_TYPE = 'NN_multiFANN';


% If Trial is not a MTASession try loading it.
Trial = MTATrial.validate(Trial);

varargout = cell([1,nargout-1]);

% Default Arguments for uninitiated variables
defArgs = {...
           ... states
               {'walk','rear','turn','pause','groom','sit'}             ...
           ...
           ... stcMode
               '',                                                      ...
           ...
           ... featureSet 
               'fet_tsne_rev3',                                         ...
           ...
           ... sampleRate
               12,                                                      ...
           ...
           ... model
               [],                                                      ...
           ...
           ... nNeurons          
               100,                                                     ...
           ...
           ... nIter
               100,                                                     ...
           ...
           ... randomizationMethod
               '',                                                      ...
           ...
           ... map2reference
               false,                                                   ... 
           ...
           ... normalize
               false                                                    ... 
           ...
           ... tag
               ''                                                       ...
           ...
           ... targetState
               ''                                                       ...
          };

[states,stcMode,featureSet,sampleRate,model,nNeurons,...
 nIter,randomizationMethod,map2reference,normalize,tag,...
 targetState] = DefaultArgs(varargin,defArgs);


% Load the feature set
if isa(featureSet,'MTADfet'),
    features = featureSet.copy;
    featureSet = features.label;
else,
    features = feval(featureSet,Trial,sampleRate,false);
end


% Load the hand labels
if ischar(stcMode),
    Trial.load('stc',stcMode);
    StcHL = Trial.stc.copy;
elseif isa(stcMode,'MTAStateCollection'),
    StcHL = stcMode.copy;
    stcMode = stcMode.mode;
end
StcHL.states = StcHL(states{:});
keys = cellfun(@subsref,StcHL(states{:}),repmat({substruct('.','key')},size(states)));


% Default mode is labeling
train = false;

% If the model name is empty then create a composite name and 
% toggle train to create new neural network models
if isempty(model),
    model = ['MTAC_BATCH-' tag targetState featureSet ...
             '_SR_'  num2str(sampleRate) ...
             '_NORM_' num2str(normalize) ...             
             '_REF_' Trial.filebase ...
             '_STC_' stcMode ...
             '_NN_'  num2str(nNeurons)...
             '_NI_'  num2str(nIter)...
             '_' MODEL_TYPE];
    train = true;
end


if map2reference,
    % if model exists and the feautures should be mapped to a reference
    % session, parse said reference session from the model
    % NOTE: Replace this with a hash reference
    filebasePattern = '[a-zA-Z]{1,2}[0-9]{2}[-][0-9]{8,8}(\.[a-zA-Z0-9]+){2,2}';
    refSession = regexp(model,filebasePattern,'match');
    refSession = strsplit(refSession{1},'.');
    refSession = refSession([1,3,2]);
    refSession = MTATrial(refSession{:});
    
    % parse and load the stcMode used in nn training
    stcModePattern = ['STC_(.+)_NN_' num2str(nNeurons)];
    refStcMode = regexp(model,stcModePattern,'tokens');
    if ~isempty(refStcMode),
        refSession.load('stc',refStcMode{1}{1});
    end
    
    % Map features via linear or circular shift to the training Session
    % of the Neural Network.
    features.map_to_reference_session(Trial,refSession);
    if normalize,
        rfet = feval(featureSet,refSession,sampleRate,false);
        [~,refMean,refStd] = nunity(rfet(refSession.stc{'a'},:));
        features.unity([],refMean,refStd);
    end
elseif normalize,
   [~,fetMean,fetStd] = nunity(features(Trial.stc{'a'},:));
    features.unity([],fetMean,fetStd);
end



% Use xyz as base for labeling output ... quit staring at this part.
xyz = Trial.load('xyz');

nStates = numel(states);
if ~isempty(targetState),
    nStates = 2;
end

% Initialize outputs
p_state = zeros([xyz.size(1),nStates]);
d_state = zeros([xyz.size(1),nStates]);
labelingStatsMulti.confussionMatrix = zeros([nIter,nStates,nStates]);
labelingStatsMulti.precision =   zeros([nIter,nStates]);
labelingStatsMulti.sensitivity = zeros([nIter,nStates]);
labelingStatsMulti.accuracy =    zeros([nIter,1]);

% Initialize labeling timeperiods
labelingEpochs = Trial.stc{'a'}.cast('TimeSeries');

if ~isempty(targetState),
    compState = states(cellfun(@isempty,regexp(states,['^',targetState,'$'])));
    trainingStates = {[strjoin(compState,'+')],targetState};
else
    trainingStates = states;
end


% load hand labeled state matrix
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,trainingStates)),'sampleRate',xyz.sampleRate);


for iter = 1:nIter,
    try,        
        % 00:30 FTS
        if iter==1&&train
            model = [model '_RAND_' randomizationMethod];
        end
        if iter==1,
            model_out = model;
            model_path = fileparts(mfilename('fullpath'));
            mkdir(fullfile(model_path,model));
            model = [model '/' model];
        end


        if train,

            switch randomizationMethod
              case 'ERS' % equal_restructured_sampling
                nRndPeriods = 100; 
                prctTrain = 70;
                trainingEpochs = [];
                % Create MTAStateColletion based on the blocks of
                % data assymbled during the resampling.
                StcRnd = StcHL.copy;
                StcRnd.updateMode(['RAND_ERS' num2str(nRndPeriods) '_' StcHL.mode]);
                StcRnd.states = {};
                StcLab = StcHL.copy;
                StcLab.states = {};

                trainingFeatures = features.copy;
                trainingFeatures.clear;
                
                stateBlockSize = 0;

                trainingPerInds = {};
                labelingPerInds = {};
                for s = StcHL(trainingStates{:}),
                    s = s{1};
                    rprs = randperm(s.size(1));
                    % select random half of the periods
                    trainingPerInds(end+1) = {rprs(1:round(s.size(1).*prctTrain/100))};
                    labelingPerInds(end+1) = {rprs((round(s.size(1).*prctTrain/100)+1):s.size(1))};
                    l = s.copy;
                    l.data = s.data(labelingPerInds{end},:);
                    s.data = s.data(trainingPerInds{end}(...
                                          randi(numel(trainingPerInds{end}),nRndPeriods,1)),:);
                    s.resample(features);
                    
                    
                    % Select first 0.1 seconds of state
                    s_start = s.copy;
                    s_start.data = [s_start.data(:,1),...
                                    bsxfun(@plus,s_start.data(:,1),round(0.1*s.sampleRate))];
                    s_start = s_start&s;
                    for per = s_start.data',
                        trainingFeatures.data = cat(1,...
                                                    trainingFeatures.data,...
                                                    features(per(1):per(2),:)...
                                                );
                    end
                    
                    % Randomly sample 0.5 seconds of middle section
                    % (oversample if state is less than 0.7 seconds)
                    s_middle = s.copy;
                    %s_middle = s+[.1,-.1];
                    for per = s_middle.data',
                        try
                            rsamp = randi(per(2)-per(1),...
                                          min([per(2)-per(1),round(0.5*sampleRate)]),...
                                          1);
                        catch
                            continue;
                        end
                        trainingFeatures.data = cat(1,...
                                                    trainingFeatures.data,...
                                                    features((per(1)+rsamp-1):(per(1)+rsamp),:)...
                                                );
                    end
                    
                    % Select last 0.1 seconds of state
                    s_end = s.copy;
                    s_end.data = [bsxfun(@minus,s_end.data(:,2),round(0.1*s.sampleRate)),...
                                  s_end.data(:,2)];
                    s_end = s_end&s;
                    for per = s_end.data',
                        trainingFeatures.data = cat(1,...
                                                    trainingFeatures.data,...
                                                    features(per(1):per(2),:)...
                                                );
                    end

                    % At this point I should have a synthetic
                    % feature matrix which consists of big blocks
                    % of state data
                    
                    % keep track of which blocks contain the states
                    s.data = [stateBlockSize+1,trainingFeatures.size(1)];
                    StcRnd.states{end+1} = s.copy;
                    StcLab.states{end+1} = l.copy;
                    stateBlockSize = trainingFeatures.size(1);
                end

                labelingEpochs = MTADepoch([],[],...
                                           any(stc2mat(StcLab,features,trainingStates),2),...
                                           features.sampleRate,... 
                                           features.sync.copy, ... 
                                           features.origin,    ... 
                                           'TimeSeries',       ...
                                           [],[],              ...
                                           'labeling',         ... label
                                           'l'                 ... key
                );

                        
                   
              case 'WSB'   % whole state bootstrap
                [StcRnd,labelingEpochs,trainingFeatures] = ...
                    resample_whole_state_bootstrap(StcHL,features,states);
                trainingEpochs = [];
              case 'WSBN'  % whole state bootstrap noisy
                [StcRnd,labelingEpochs,trainingFeatures] = ...
                    resample_whole_state_bootstrap_noisy(StcHL,features,states);
                trainingEpochs = [];
              case 'WSBNT' % whole state bootstrap noisy with trimmed boundaries
                [StcRnd,labelingEpochs,trainingFeatures] = ...
                    resample_whole_state_bootstrap_noisy_trim(StcHL,features,states);
                trainingEpochs = [];
              case 'rndsamp'
                rndInd = randperm(features.size(1))';
                rndInd = rndInd(1:floor(features.size(1)/2));
                trndInd = false([features.size(1),1]);
                trndInd(rndInd) = true;
                rndInd = trndInd;
                trainingFeatures = features.copy;
                StcRnd = trainingStates;


                trainingEpochs = MTADepoch([],[],              ...
                                           rndInd,             ...
                                           features.sampleRate,... 
                                           features.sync.copy, ... 
                                           features.origin,    ... 
                                           'TimeSeries',       ...
                                           [],[],              ...
                                           'training',         ... label
                't'                ... key
                );

                labelingEpochs = MTADepoch([],[],...
                                           ~rndInd,...
                                           features.sampleRate,... 
                                           features.sync.copy, ... 
                                           features.origin,    ... 
                                           'TimeSeries',       ...
                                           [],[],              ...
                                           'labeling',         ... label
                'l'                ... key
                );
            end

            
            % Train Model
            bhv_fann (Trial,                              ... Trial
                    true,                               ... ifTrain
                    trainingStates,                     ... States
                    StcRnd,                             ... StateCollection
                    trainingFeatures,                   ... feature set
                    [model '_' num2str(iter)],          ... model name
                    'subset',trainingEpochs);                                        
            

        end
        
        % Label States
        [Stc,ps,Model_Information] = bhv_fann (Trial,     ... Trial
                                             false,     ... ifTrain
                                             trainingStates,    ... States
                                             StcHL,     ... StateCollection
                                             features,  ... feature set
                                             [model '_' num2str(iter)]); % model name

        


        ysm = MTADxyz('data',double(0<stc2mat(Stc,xyz)),'sampleRate',xyz.sampleRate); 
        d_state = ysm.data+d_state;
        p_state = p_state +ps;
        if nargout>=4,
            

            labelingEpochs.resample(xyz);

            ind = any(shl.data,2)&any(ysm.data,2)&labelingEpochs.data;

            tcm = confmat(shl(ind&labelingEpochs,:),ysm(ind&labelingEpochs,:)); % DEP: netlab
            labelingStatsMulti.confusionMatrix(iter,:,:) = round(tcm./xyz.sampleRate,2);
            labelingStatsMulti.precision(iter,:) = round(diag(tcm)./sum(tcm,2),4).*100;
            labelingStatsMulti.sensitivity(iter,:) = round(diag(tcm)'./sum(tcm),4).*100;
            labelingStatsMulti.accuracy(iter) = sum(diag(tcm))/sum(tcm(:));
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
d_state = MTADxyz('data',d_state,'sampleRate',xyz.sampleRate);
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
[~,maxState] = max(d_state.data,[],2);
maxState(~nniz(xyz),:) = 0;

% Smooth decision boundaries - 200 ms state minimum
% $$$ bwin = round(.2*xyz.sampleRate)+double(mod(round(.2*xyz.sampleRate),2)==0);
% $$$ mss = GetSegs(maxState,1:size(maxState,1),bwin,nan);
% $$$ maxState=circshift(sq(mode(mss))',floor(bwin/2));


% Stats in comparision to the collection of labels specified in the stcMode
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
