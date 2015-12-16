% req20151216
% Primary goal - randomized train/label validation (50/50) bhv_nn



% Net feature set fet_all
% Train a logistic regresion model using features of jg05-20120317
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');

states = {'walk','rear','turn','pause','groom','sit'};

%features = fet_tsne(Trial,20,false);
features = fet20151007(Trial,20,false);
model = 'fet20151007_REFjg0520120317_NN';

nNeurons = 100;
nIter = 100;

cm = zeros([100,numel(states),numel(states)]);
precision =   zeros([100,numel(states)]);
sensitivity = zeros([100,numel(states)]);
for iter = 1:nIter,
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
    model = 'fet_tsne_REFjg0520120317_NN';
    bhv_nn (Trial,                                           ... Trial
            true,                                            ... ifTrain
            states,                                          ... States
            features,                                        ... feature set
            model,                                           ... model name
            'subset',trainingEpochs);                                        

    % Label States
    Stc = bhv_nn (Trial,                                         ... Trial
                  false,                                         ... ifTrain
                  states,                                        ... States
                  features,                                      ... feature set
                  model);                                          % model name

    
    %
    %StcHL = Trial.load('stc','hand_labeled_rev1');
    StcHL = Trial.load('stc','hand_labeled_rev2');
    xyz = Trial.load('xyz');
    shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
    ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,states)),'sampleRate',xyz.sampleRate); 

    labelingEpochs.resample(xyz);

    ind = any(shl.data,2)&any(ysm.data,2)&labelingEpochs.data;

    tcm = confmat(shl(ind&labelingEpochs,:),ysm(ind&labelingEpochs,:)); % DEP: netlab
    precision(iter,:) = round(diag(tcm)'./sum(tcm),4).*100;
    sensitivity(iter,:) = round(diag(tcm)./sum(tcm,2),4).*100;
    cm(iter,:,:) = round(tcm./xyz.sampleRate,2);

end





bhv_nn (Trial,                                               ... Trial
        true,                                                ... ifTrain
        states,                                              ... States
        {'fet20151007',Trial,20,false},                          ... feature set
        'fet20151007_REFjg0520120317_NN');                                % model name

        
bhv_nn (Trial,                                               ... Trial
        true,                                                ... ifTrain
        states,                                              ... States
        {'fet_tsne',Trial,20,true},                          ... feature set
        'Ufet_tsne_REFjg0520120317_NN');                                % model name

        
bhv_nn (Trial,                                               ... Trial
        true,                                                ... ifTrain
        states,                                              ... States
        {'fet20151007',Trial,20,true},                          ... feature set
        'Ufet20151007_REFjg0520120317_NN');                                % model name
        
        
        

Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
states = {'walk','rear','turn','pause','groom','sit'};
features = fet_tsne(Trial,20,false);


StcHL = Trial.load('stc','hand_labeled_rev2');
xyz = Trial.load('xyz');
shl = MTADxyz('data',double(0<stc2mat(StcHL,features,states)),'sampleRate',features.sampleRate);



rndInd = randperm(features.size(1))';
rndInd = rndInd(1:floor(features.size(1)/2));
trndInd = false([features.size(1),1]);
trndInd(rndInd) = true;
rndInd = trndInd;

smat = stc2mat(Trial.stc,features,states);
ind = any(smat,2)&rndInd;

net = patternnet(nNeurons);
[net,tr] = train(net,features(ind,:)',~~smat(ind,:)');

d_state = net(features.data')';

% Separate the winners from the losers
[~,maxState] = max(d_state,[],2);

maxState(~nniz(features),:) = 0;
msi = 0:numel(states):numel(states).*(features.size(1)-1);
ysl = zeros([features.size(1),numel(states)])';
ysl(maxState+msi') = 1;
ysl = ysl';

cm = confmat(shl(~ind,:),ysl(~ind,:)); % DEP: netlab                
precision = round(diag(cm)'./sum(cm),4).*100;
sensitivity = round(diag(cm)./sum(cm,2),4).*100;
cm = round(cm./xyz.sampleRate,2);




