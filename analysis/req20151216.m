% req20151216
% Primary goal - randomized train/label validation (50/50) bhv_nn



% Net feature set fet_all
% Train a logistic regresion model using features of jg05-20120317
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');

states = {'walk','rear','turn','pause','groom','sit'};
featureSet = 'fet_tsne_rev1';
fetSampleRate = 20;%Hz
nNeurons = 150;
nIter = 200;


features = feval(featureSet,Trial,fetSampleRate,false);
model = ['MTCA_BATCH-' featureSet '-REF-' Trial.name '-NN'];

%features = fet20151007(Trial,fetSampleRate,false);
%model = 'fet20151007_REFjg0520120317_NN';


cm = zeros([nIter,numel(states),numel(states)]);
precision =   zeros([nIter,numel(states)]);
sensitivity = zeros([nIter,numel(states)]);
accuracy =    zeros([nIter,1]);
for iter = 1:nIter,
    try,
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
            [model '_' num2str(iter)],                                           ... model name
            'subset',trainingEpochs);                                        

    % Label States
    Stc = bhv_nn (Trial,                                         ... Trial
                  false,                                         ... ifTrain
                  states,                                        ... States
                  features,                                      ... feature set
                  [model '_' num2str(iter)]);                                          % model name

    
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
    accuracy(iter) = sum(diag(tcm))/sum(tcm(:));
    catch
        continue
    end
    
end

%save(['/storage/gravio/data/project/general/analysis/req20151216-' ...
%     model '-' Trial.filebase '.mat']);
load(['/storage/gravio/data/project/general/analysis/req20151216-' ...
     model '-' Trial.filebase '.mat']);


% Plot the distribution of model precision
plot_precision(precision,states,'fet_tsne')


%% Inter subject labeling using jg05 train nn's

Trial = MTATrial('Ed01-20140707');
featureName = 'fet_tsne';
%features = fet_tsne(Trial,fetSampleRate,false);

featureName = 'fet20151007';

features = feval(featureName,Trial,fetSampleRate,false);


cm = zeros([nIter,numel(states),numel(states)]);
precision =   zeros([nIter,numel(states)]);
sensitivity = zeros([nIter,numel(states)]);

xyz = Trial.load('xyz');
StcHL = Trial.load('stc','hand_labeled_rev1');
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);

for iter = 1:nIter,
    %Label states
     model = 'fet_tsne_REFjg0520120317_NN';
    Stc = bhv_nn (Trial,                                         ... Trial
                  false,                                         ... ifTrain
                  states,                                        ... States
                  features,                                      ... feature set
                  [model '_' num2str(iter)]);                      % model name


    ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,states)),'sampleRate',xyz.sampleRate); 


    ind = any(shl.data,2)&any(ysm.data,2);

    tcm = confmat(shl(ind,:),ysm(ind,:)); % DEP: netlab
    sensitivity(iter,:) = round(diag(tcm)'./sum(tcm),4).*100;
    precision(iter,:) = round(diag(tcm)./sum(tcm,2),4).*100;
    cm(iter,:,:) = round(tcm./xyz.sampleRate,2);
end

% Plot the distribution of model precision
plot_precision(precision,states,featureName)




bhv_nn (Trial,                                               ... Trial
        true,                                                ... ifTrain
        states,                                              ... States
        {'fet20151007',Trial,fetSampleRate,false},                          ... feature set
        'fet20151007_REFjg0520120317_NN');                                % model name

        
bhv_nn (Trial,                                               ... Trial
        true,                                                ... ifTrain
        states,                                              ... States
        {'fet_tsne',Trial,fetSampleRate,true},                          ... feature set
        'Ufet_tsne_REFjg0520120317_NN');                                % model name

        
bhv_nn (Trial,                                               ... Trial
        true,                                                ... ifTrain
        states,                                              ... States
        {'fet20151007',Trial,fetSampleRate,true},                          ... feature set
        'Ufet20151007_REFjg0520120317_NN');                                % model name
        
        
        

Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
states = {'walk','rear','turn','pause','groom','sit'};
features = fet_tsne(Trial,fetSampleRate,false);


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




