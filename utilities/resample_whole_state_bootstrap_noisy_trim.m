function [StcRnd,labelingEpochs,trainingFeatures] = resample_whole_state_bootstrap_noisy_trim(StcHL,features,states,varargin)

if ~isempty(varargin),
    prctTrain = varargin{1};
else
    prctTrain = 90;
end

trainingEpochs = [];

% Create MTAStateColletion based on the blocks of
% data assymbled during the resampling.
StcRnd = StcHL.copy;
StcRnd.updateMode(['RAND_WSB_' StcHL.mode]);
StcRnd.states = {};
StcLab = StcHL.copy;
StcLab.states = {};

trainingFeatures = features.copy;
trainingFeatures.clear;
tmpFeatures = features.copy;
tmpFeatures.clear;

stateBlockSize = 7500;

trainingPerInds = {};
labelingPerInds = {};

for s = StcHL(states{:}),
    s = s{1};
    rprs = randperm(s.size(1));
    % select random sets of the periods
    trainingPerInds(end+1) = {rprs(1:round(s.size(1).*prctTrain/100))};
    labelingPerInds(end+1) = {rprs((round(s.size(1).*prctTrain/100)+1):s.size(1))};
    l = s.copy;
    l.data = s.data(labelingPerInds{end},:);
    s.data = s.data(trainingPerInds{end},:);
    s = s+[0.1,-0.1];
    s.resample(features);

    tmpFeatures.clear;
    tmpFeatures.data = features(s,:);

    s.data = [trainingFeatures.size(1)+1,trainingFeatures.size(1)+stateBlockSize];
    StcRnd.states{end+1} = s.copy;
    StcLab.states{end+1} = l.copy;
    
    trainingFeatures.data = [trainingFeatures.data;
                        tmpFeatures(randi(tmpFeatures.size(1),...
                                          stateBlockSize,...
                                          1),...
                                    :)];
end

trainingFeatures.data = trainingFeatures.data+randn(trainingFeatures.size)/5;

labelingEpochs = MTADepoch([],[],...
                           any(stc2mat(StcLab,features,states),2),...
                           features.sampleRate,... 
                           features.sync.copy, ... 
                           features.origin,    ... 
                           'TimeSeries',       ...
                           [],[],              ...
                           'labeling',         ... label
                           'l'                ... key
);


