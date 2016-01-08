function [varargout] = bhv_nn_multi_patternnet(Trial,varargin)
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


MODEL_TYPE = 'NN_multiPN';


% If Trial is not a MTASession try loading it.
if ischar(Trial),
    Trial = MTATrial(Trial);
elseif iscell(Trial),
    Trial = MTATrial(Trial{:});
end

varargout = cell([1,nargout-1]);

% Default Arguments for uninitiated variables
defArgs = {...
           ... states
               {'walk','rear','turn','pause','groom','sit'}             ...
           ...
           ... stcMode
               'hand_labeled_rev2',                                     ...
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
           ... randomizationMode
               ''                                                       ...
};

[states,stcMode,featureSet,sampleRate,model,nNeurons,nIter,randomizationMode] = DefaultArgs(varargin,defArgs);


% Load the feature set
features = feval(featureSet,Trial,sampleRate,false);

% Default mode is labeling
train = false;

% If the model name is empty then create a composite name and 
% toggle train to create new neural network models
if isempty(model),
    model = ['MTAC_BATCH-' featureSet ...
             '_SR_'  num2str(sampleRate) ...
             '_REF_' Trial.filebase ...
             '_NN_'  num2str(nNeurons) ];
    %model = 'fet_tsne_REFjg0520120317_NN';
    train = true;
end


% Use xyz as base for labeling output ... quit staring at this part.
xyz = Trial.load('xyz');

% Initialize outputs
d_state = zeros([xyz.size(1),numel(states)]);
labelingStatsMulti.confussionMatrix = zeros([nIter,numel(states),numel(states)]);
labelingStatsMulti.precision =   zeros([nIter,numel(states)]);
labelingStatsMulti.sensitivity = zeros([nIter,numel(states)]);
labelingStatsMulti.accuracy =    zeros([nIter,1]);

% Initialize labeling timeperiods
labelingEpochs = Trial.stc{'a'}.cast('TimeSeries');

% Load the hand label
if ~isempty(stcMode),
    Trial.load('stc',stcMode);
    StcHL = Trial.stc.copy;
end
keys = subsref(Trial.stc.list_state_attrib('key'),...
               substruct('()',{Trial.stc.gsi(states)}));


for iter = 1:nIter,
    try,
        if train,
            % Create two partitions for validation
            %blocs  = reshape(1mod(features.size(1), 4.*20 )
            rndInd = randperm(features.size(1))';
            rndInd = rndInd(1:floor(features.size(1)/2));
            trndInd = false([features.size(1),1]);
            trndInd(rndInd) = true;
            rndInd = trndInd;


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



            % Train Model
            bhv_nn (Trial,                                           ... Trial
            true,                                            ... ifTrain
            states,                                          ... States
            features,                                        ... feature set
            [model '_' num2str(iter)],                                           ... model name
            'subset',trainingEpochs);                                        
        

        end
        % Label States
        [Stc,~,Model_Information] = bhv_nn (Trial,    ... Trial
                                             false,    ... ifTrain
                                             states,   ... States
                                             features, ... feature set
                                             [model '_' num2str(iter)]); % model name

        


        ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,states)),'sampleRate',xyz.sampleRate); 
        d_state = ysm.data+d_state;
        
        if nargout>=4,
            shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);

            labelingEpochs.resample(xyz);

            ind = any(shl.data,2)&any(ysm.data,2)&labelingEpochs.data;

            tcm = confmat(shl(ind&labelingEpochs,:),ysm(ind&labelingEpochs,:)); % DEP: netlab
            labelingStatsMulti.confusionMatrix(iter,:,:) = round(tcm./xyz.sampleRate,2);
            labelingStatsMulti.precision(iter,:) = round(diag(tcm)./sum(tcm,2),4).*100;
            labelingStatsMulti.sensitivity(iter,:) = round(diag(tcm)'./sum(tcm),4).*100;
            labelingStatsMulti.accuracy(iter) = sum(diag(tcm))/sum(tcm(:));
        end
    catch
        continue
    end
    
end



% Determine winning states based on the the labels of nurmerous
% neural networks.
[~,maxState] = max(d_state,[],2);
maxState(~nniz(xyz),:) = 0;

% Smooth decision boundaries - 200 ms state minimum
% $$$ bwin = round(.2*xyz.sampleRate)+double(mod(round(.2*xyz.sampleRate),2)==0);
% $$$ mss = GetSegs(maxState,1:size(maxState,1),bwin,nan);
% $$$ maxState=circshift(sq(mode(mss))',floor(bwin/2));


% Stats in comparision to the collection of labels specified in the stcMode
ysm = MTADxyz('data',zeros([shl.size]),'sampleRate',xyz.sampleRate); 
ysm.data = ysm.data';
ysm.data([1:6:size(ysm,2).*6]+maxState'-1) = 1;
ysm.data = ysm.data';
ind = any(shl.data,2)&any(ysm.data,2)&labelingEpochs.data;
tcm = confmat(shl(ind&labelingEpochs,:),ysm(ind&labelingEpochs,:)); % #DEP: netlab
labelingStats.confusionMatrix = round(tcm./xyz.sampleRate,2);
labelingStats.precision = round(diag(tcm)./sum(tcm,2),4).*100;
labelingStats.sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
labelingStats.accuracy = sum(diag(tcm))/sum(tcm(:));

% Copy State Collection object to store new labeled periods
Stc = Trial.stc.copy;
Stc.updateMode([MODEL_TYPE '-' Model_Information.StcMode...
                '-' cell2mat(Model_Information.state_keys)]);
Stc.states = {};



% Create new StateCollection ... well copy
Stc = Trial.stc.copy;
Stc.updateMode([MODEL_TYPE '-' Model_Information.Trial '-'...
                Model_Information.StcMode...
                '-' cell2mat(Model_Information.state_keys)]);
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
if nargout>=2, varargout{2} = d_state;            end
if nargout>=3, varargout{3} = labelingStats;      end
if nargout>=4, varargout{4} = labelingStatsMulti; end
    