% req20151216
% Primary goal - randomized train/label validation (50/50) bhv_nn

Trial = MTATrial('jg05-20120317');
mod.states     = {'walk','rear','turn','pause','groom','sit'};
mod.stcMode    = 'hand_labeled_rev2';
mod.featureSet = 'fet_tsne_rev5';
mod.sampleRate = 12;
mod.nNeurons   = 100;
mod.randomizationMethod = 'equal_restructured_sampling';

argin = struct2varargin(mod);
[stc,d_state,ls,lsm,model] = bhv_nn_multi_patternnet(Trial,argin{:});



Trial = MTATrial('jg05-20120317');
clear('mod');
mod.states     = {'walk','rear','turn','pause','groom','sit'};
mod.stcMode    = 'hand_labeled_rev2';
mod.featureSet = 'fet_tsne_rev5';
mod.sampleRate = 12;
mod.nNeurons   = 100;
mod.randomizationMethod = 'whole_state_bootstrap';

argin = struct2varargin(mod);
[stc,d_state,ls,lsm,model] = bhv_nn_multi_patternnet(Trial,argin{:});




Trial = MTATrial('Ed03-20140624');
mod.states     = {'walk','rear','turn','pause','groom','sit'};
mod.stcMode    = 'hand_labeled_rev2_alt';
mod.featureSet = 'fet_tsne_rev5';
mod.model      = model;
mod.sampleRate = 12;
mod.nNeurons   = 100;
mod.randomizationMethod = 'equal_restructured_sampling';

argin = struct2varargin(mod);
[stc,d_state,ls,lsm] = bhv_nn_multi_patternnet(Trial,argin{:});



Trial = MTATrial('Ed03-20140624');
mod.states     = {'walk','rear','turn','pause','groom','sit'};
mod.stcMode    = 'hand_labeled_rev1';
mod.featureSet = 'fet_tsne_rev5';
mod.model      = model;
mod.sampleRate = 12;
mod.nNeurons   = 100;
mod.randomizationMethod = 'whole_state_bootstrap';

argin = struct2varargin(mod);
[stc,d_state,ls,lsm] = bhv_nn_multi_patternnet(Trial,argin{:});




Trial = MTATrial('Ed03-20140624');
mod.states     = {'walk','rear','turn','pause','groom','sit'};
mod.stcMode    = 'hand_labeled_rev1';
mod.featureSet = 'fet_tsne_rev3';
mod.model      = 'MTAC_BATCH-fet_tsne_rev3_SR_12_REF_jg05-20120317.cof.all_NN_100';
mod.sampleRate = 12;
mod.nNeurons   = 100;

argin = struct2varargin(mod);
[stc,d_state,ls,lsm] = bhv_nn_multi_patternnet(Trial,argin{:});


Trial = MTATrial('Ed03-20140625');
mod.states     = {'walk','rear','turn','pause','groom','sit'};
mod.stcMode    = 'hand_labeled_rev1';
mod.featureSet = 'fet_tsne_rev3';
mod.model      = 'MTAC_BATCH-fet_tsne_rev3_SR_12_REF_jg05-20120317.cof.all_NN_100';
mod.sampleRate = 12;
mod.nNeurons   = 100;
argin = struct2varargin(mod);

[stc,d_state,ls,lsm] = bhv_nn_multi_patternnet(Trial,argin{:});


Trial = MTATrial('Ed03-20140625');
mod.states     = {'walk','rear','turn','pause','groom','sit'};
mod.stcMode    = 'hand_labeled_rev1';
mod.featureSet = 'fet_tsne_rev3';
mod.model      = 'MTAC_BATCH-fet_tsne_rev3_SR_12_REF_Ed03-20140624.cof.all_NN_100';
mod.sampleRate = 12;
mod.nNeurons   = 100;
argin = struct2varargin(mod);

[stc,d_state,ls,lsm] = bhv_nn_multi_patternnet(Trial,argin{:});


Trial = MTATrial('Ed03-20140624');
mod.states     = {'walk','rear','turn','pause','groom','sit'};
mod.stcMode    = 'NN_multiPN-jg05-20120317.cof.all-hand_labeled_rev2-wrnpms';
mod.featureSet = 'fet_tsne_rev3';
mod.model      = 'MTAC_BATCH-fet_tsne_rev3_SR_12_REF_Ed03-20140625.cof.all_NN_100';
mod.sampleRate = 12;
mod.nNeurons   = 100;
argin = struct2varargin(mod);

[stc,d_state,ls,lsm] = bhv_nn_multi_patternnet(Trial,argin{:});







%% STuff below has been moved into the classifier function named
%% bhv_nn_multi_patternnet

% $$$ 
% $$$ % Net feature set fet_all
% $$$ % Train a logistic regresion model using features of jg05-20120317
% $$$ Trial = MTATrial('jg05-20120317');
% $$$ Trial.load('stc','hand_labeled_rev2');
% $$$ 
% $$$ states = {'walk','rear','turn','pause','groom','sit'};
% $$$ featureSet = 'fet_tsne_rev2';
% $$$ 
% $$$ nIter = 200;
% $$$ nNeurons = 150;
% $$$ 
% $$$  
% $$$ sampleRate = 10;%Hz
% $$$ features = feval(featureSet,Trial,sampleRate,false);
% $$$ 
% $$$ model = ['MTAC_BATCH-' featureSet ...
% $$$          '_SR_'  num2str(sampleRate) ...
% $$$          '_REF_' Trial.name ...
% $$$          '_NN_'  num2str(nNeurons) ];
% $$$ 
% $$$ 
% $$$ %features = fet20151007(Trial,sampleRate,false);
% $$$ %model = 'fet20151007_REFjg0520120317_NN';
% $$$ 
% $$$ 
% $$$ cm = zeros([nIter,numel(states),numel(states)]);
% $$$ precision =   zeros([nIter,numel(states)]);
% $$$ sensitivity = zeros([nIter,numel(states)]);
% $$$ accuracy =    zeros([nIter,1]);
% $$$ for iter = 1:nIter,
% $$$     try,
% $$$     % Create two partitions for validation
% $$$     %blocs  = reshape(1mod(features.size(1), 4.*20 )
% $$$     rndInd = randperm(features.size(1))';
% $$$     rndInd = rndInd(1:floor(features.size(1)/2));
% $$$     trndInd = false([features.size(1),1]);
% $$$     trndInd(rndInd) = true;
% $$$     rndInd = trndInd;
% $$$ 
% $$$ 
% $$$     trainingEpochs = MTADepoch([],[],              ...
% $$$                                rndInd,             ...
% $$$                                features.sampleRate,... 
% $$$                                features.sync.copy, ... 
% $$$                                features.origin,    ... 
% $$$                                'TimeSeries',       ...
% $$$                                [],[],              ...
% $$$                                'training',         ... label
% $$$                                't'                ... key
% $$$     );
% $$$ 
% $$$     labelingEpochs = MTADepoch([],[],...
% $$$                                ~rndInd,...
% $$$                                features.sampleRate,... 
% $$$                                features.sync.copy, ... 
% $$$                                features.origin,    ... 
% $$$                                'TimeSeries',       ...
% $$$                                [],[],              ...
% $$$                                'labeling',         ... label
% $$$                                'l'                ... key
% $$$     );
% $$$ 
% $$$ 
% $$$ 
% $$$     % Train Model
% $$$     model = 'fet_tsne_REFjg0520120317_NN';
% $$$     bhv_nn (Trial,                                           ... Trial
% $$$             true,                                            ... ifTrain
% $$$             states,                                          ... States
% $$$             features,                                        ... feature set
% $$$             [model '_' num2str(iter)],                                           ... model name
% $$$             'subset',trainingEpochs);                                        
% $$$ 
% $$$     % Label States
% $$$     Stc = bhv_nn (Trial,                                         ... Trial
% $$$                   false,                                         ... ifTrain
% $$$                   states,                                        ... States
% $$$                   features,                                      ... feature set
% $$$                   [model '_' num2str(iter)]);                                          % model name
% $$$ 
% $$$     
% $$$     %
% $$$     %StcHL = Trial.load('stc','hand_labeled_rev1');
% $$$     StcHL = Trial.load('stc','hand_labeled_rev2');
% $$$     xyz = Trial.load('xyz');
% $$$     shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
% $$$     ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,states)),'sampleRate',xyz.sampleRate); 
% $$$ 
% $$$     labelingEpochs.resample(xyz);
% $$$ 
% $$$     ind = any(shl.data,2)&any(ysm.data,2)&labelingEpochs.data;
% $$$ 
% $$$     tcm = confmat(shl(ind&labelingEpochs,:),ysm(ind&labelingEpochs,:)); % DEP: netlab
% $$$     precision(iter,:) = round(diag(tcm)'./sum(tcm),4).*100;
% $$$     sensitivity(iter,:) = round(diag(tcm)./sum(tcm,2),4).*100;
% $$$     cm(iter,:,:) = round(tcm./xyz.sampleRate,2);
% $$$     accuracy(iter) = sum(diag(tcm))/sum(tcm(:));
% $$$     catch
% $$$         continue
% $$$     end
% $$$     
% $$$ end
% $$$ 
% $$$ %save(['/storage/gravio/data/project/general/analysis/req20151216-' ...
% $$$ %     model '-' Trial.filebase '.mat']);
% $$$ load(['/storage/gravio/data/project/general/analysis/req20151216-' ...
% $$$      model '-' Trial.filebase '.mat']);
% $$$ 
% $$$ 
% $$$ % Plot the distribution of model precision
% $$$ plot_precision(precision,states,'fet_tsne')
% $$$ 
% $$$ 
% $$$ %% Inter subject labeling using jg05 train nn's
% $$$ 
% $$$ Trial = MTATrial('Ed01-20140707');
% $$$ featureName = 'fet_tsne';
% $$$ %features = fet_tsne(Trial,sampleRate,false);
% $$$ 
% $$$ featureName = 'fet20151007';
% $$$ 
% $$$ features = feval(featureName,Trial,sampleRate,false);
% $$$ 
% $$$ 
% $$$ cm = zeros([nIter,numel(states),numel(states)]);
% $$$ precision =   zeros([nIter,numel(states)]);
% $$$ sensitivity = zeros([nIter,numel(states)]);
% $$$ 
% $$$ xyz = Trial.load('xyz');
% $$$ StcHL = Trial.load('stc','hand_labeled_rev1');
% $$$ shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
% $$$ 
% $$$ for iter = 1:nIter,
% $$$     %Label states
% $$$      model = 'fet_tsne_REFjg0520120317_NN';
% $$$     Stc = bhv_nn (Trial,                                         ... Trial
% $$$                   false,                                         ... ifTrain
% $$$                   states,                                        ... States
% $$$                   features,                                      ... feature set
% $$$                   [model '_' num2str(iter)]);                      % model name
% $$$ 
% $$$ 
% $$$     ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,states)),'sampleRate',xyz.sampleRate); 
% $$$ 
% $$$ 
% $$$     ind = any(shl.data,2)&any(ysm.data,2);
% $$$ 
% $$$     tcm = confmat(shl(ind,:),ysm(ind,:)); % DEP: netlab
% $$$     sensitivity(iter,:) = round(diag(tcm)'./sum(tcm),4).*100;
% $$$     precision(iter,:) = round(diag(tcm)./sum(tcm,2),4).*100;
% $$$     cm(iter,:,:) = round(tcm./xyz.sampleRate,2);
% $$$ end
% $$$ 
% $$$ % Plot the distribution of model precision
% $$$ plot_precision(precision,states,featureName)
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ bhv_nn (Trial,                                               ... Trial
% $$$         true,                                                ... ifTrain
% $$$         states,                                              ... States
% $$$         {'fet20151007',Trial,sampleRate,false},                          ... feature set
% $$$         'fet20151007_REFjg0520120317_NN');                                % model name
% $$$ 
% $$$         
% $$$ bhv_nn (Trial,                                               ... Trial
% $$$         true,                                                ... ifTrain
% $$$         states,                                              ... States
% $$$         {'fet_tsne',Trial,sampleRate,true},                          ... feature set
% $$$         'Ufet_tsne_REFjg0520120317_NN');                                % model name
% $$$ 
% $$$         
% $$$ bhv_nn (Trial,                                               ... Trial
% $$$         true,                                                ... ifTrain
% $$$         states,                                              ... States
% $$$         {'fet20151007',Trial,sampleRate,true},                          ... feature set
% $$$         'Ufet20151007_REFjg0520120317_NN');                                % model name
% $$$         
% $$$         
% $$$         
% $$$ 
% $$$ Trial = MTATrial('jg05-20120317');
% $$$ Trial.load('stc','hand_labeled_rev2');
% $$$ states = {'walk','rear','turn','pause','groom','sit'};
% $$$ features = fet_tsne(Trial,sampleRate,false);
% $$$ 
% $$$ 
% $$$ StcHL = Trial.load('stc','hand_labeled_rev2');
% $$$ xyz = Trial.load('xyz');
% $$$ shl = MTADxyz('data',double(0<stc2mat(StcHL,features,states)),'sampleRate',features.sampleRate);
% $$$ 
% $$$ 
% $$$ 
% $$$ rndInd = randperm(features.size(1))';
% $$$ rndInd = rndInd(1:floor(features.size(1)/2));
% $$$ trndInd = false([features.size(1),1]);
% $$$ trndInd(rndInd) = true;
% $$$ rndInd = trndInd;
% $$$ 
% $$$ smat = stc2mat(Trial.stc,features,states);
% $$$ ind = any(smat,2)&rndInd;
% $$$ 
% $$$ net = patternnet(nNeurons);
% $$$ [net,tr] = train(net,features(ind,:)',~~smat(ind,:)');
% $$$ 
% $$$ d_state = net(features.data')';
% $$$ 
% $$$ % Separate the winners from the losers
% $$$ [~,maxState] = max(d_state,[],2);
% $$$ 
% $$$ maxState(~nniz(features),:) = 0;
% $$$ msi = 0:numel(states):numel(states).*(features.size(1)-1);
% $$$ ysl = zeros([features.size(1),numel(states)])';
% $$$ ysl(maxState+msi') = 1;
% $$$ ysl = ysl';
% $$$ 
% $$$ cm = confmat(shl(~ind,:),ysl(~ind,:)); % DEP: netlab                
% $$$ precision = round(diag(cm)'./sum(cm),4).*100;
% $$$ sensitivity = round(diag(cm)./sum(cm,2),4).*100;
% $$$ cm = round(cm./xyz.sampleRate,2);
% $$$ 
% $$$ 
% $$$ 
% $$$ 


%Batch stuff
out = {};
nNeurons = 60:10:150;
for nn = nNeurons,
    Trial = MTATrial('jg05-20120317');
    mod.states     = {'walk','rear','turn','pause','groom','sit'};
    mod.stcMode    = 'hand_labeled_rev2';
    mod.featureSet = 'fet_tsne_rev3';
    mod.model      = '';
    mod.sampleRate = 10;
    mod.nNeurons   = nn;
    argin = struct2varargin(mod);
    [out{end+1}.stc,out{end+1}.d_state,out{end+1}.ls,out{end+1}.lsm] =...
        bhv_nn_multi_patternnet(Trial,argin{:});
end


cellfun(@subsref,out,repmat({substruct('.','ls','.','accuracy')},size(out)))

figure,plot(cellfun(@subsref,out,repmat({substruct('.','ls','.','accuracy')},size(out))))



% Maybe add this at some point for fun :) 1000 nn models with different num of neurons.
Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
StcHL = Trial.load('stc','hand_labeled_rev2');
labelingEpochs = Trial.stc{'a'}.cast('TimeSeries');

figure,
imagesc(sum(reshape(cell2mat(cellfun(@subsref,out,repmat({substruct('.','d_state')},size(out)),'UniformOutput',false)),[],6,10),3)')


ysm = MTADxyz('data',sum(reshape(cell2mat(...
    cellfun(@subsref,out,repmat({substruct('.','d_state')},size(out)),'UniformOutput',false)),...
                  [],numel(mod.states),numel(out)),3),...
              'sampleRate',Trial.xyz.sampleRate);
ysm.filter('ButFilter',3,3,'low');
ysm.data = bsxfun(@eq,mind,1:6);
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,mod.states)),'sampleRate',xyz.sampleRate);

ind = any(shl.data,2)&any(ysm.data,2)&labelingEpochs.data;

tcm = confmat(shl(ind&labelingEpochs,:),ysm(ind&labelingEpochs,:)); % DEP: netlab
con = round(tcm./xyz.sampleRate,2);
pre = round(diag(tcm)./sum(tcm,2),4).*100;
sen = round(diag(tcm)'./sum(tcm),4).*100;
acc = sum(diag(tcm))/sum(tcm(:));

Trial = MTATrial('jg05-20120317');
Trial = MTATrial('Ed03-20140624');
xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz = fxyz.filter('ButFilter',3,2.4,'low');
vxy = fxyz.vel([],[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);
ang = create(MTADang,Trial,xyz);
 


name = 'pbj'; label = 'pb'; key = 'b';
zv = MTADfet.encapsulate(Trial,...
                         ang(:,10,12,3),...
                         xyz.sampleRate,...
                         name,label,key);

dspec = struct('nFFT',2^8,'Fs',zv.sampleRate,...
               'WinLength',2^7,'nOverlap',2^7*.875,...
               'FreqRange',[1,30]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',false,[],dspec,false);

name = 'zv and wag pw'; label = 'zvwp'; key = 'w';
zfet = MTADfet.encapsulate(Trial,...
                         log10(nanmedian(ys(:,round(fs)==10,1,1),2)),...
                         ys.sampleRate,...
                         name,label,key);
zfet.data(~nniz(zfet),:)=-20;





% $$$ figure,plot(mang.data)
% $$$ Lines(Trial.stc{'w'}(:),[],'b');
%eds = linspace(.2,1.5,100);

arv = zfet.copy;eds = linspace(-12,-3,100);

arv = ys.copy;eds = linspace(-10,-2,100);
arv.data = log10(arv.data);
arv = MTADang('data',ang(:,10,12,3),'sampleRate',ang.sampleRate);eds = linspace(30,70,100);

f = 1;
figure,hold on
ind = Trial.stc{'p'};
hs = bar(eds,histc(arv(ind,f),eds),'histc');hs.FaceAlpha=0.5;hs.FaceColor='c';
ind = Trial.stc{'w'};
hs = bar(eds,histc(arv(ind,f),eds),'histc');hs.FaceAlpha=0.5;hs.FaceColor='b';
ind = Trial.stc{'m'};
hs = bar(eds,histc(arv(ind,f),eds),'histc');hs.FaceAlpha=0.5;hs.FaceColor='m';
ind = Trial.stc{'n'};
hs = bar(eds,histc(arv(ind,f),eds),'histc');hs.FaceAlpha=0.5;hs.FaceColor='g';
ind = Trial.stc{'s'};
hs = bar(eds,histc(arv(ind,f),eds),'histc');hs.FaceAlpha=0.5;hs.FaceColor='y';
ind = Trial.stc{'r'};
hs = bar(eds,histc(arv(ind,f),eds),'histc');hs.FaceAlpha=0.5;hs.FaceColor='r';




arv2 = vxy.copy;a2i = 1;  eds2 = linspace(-3,2,100);
arv = ys.copy;  a1i = 19;  eds = linspace(-6,2,100);
arv.resample(arv2);
arv.data = log10(arv.data);

arv2 = zv.copy;a2i = 1;  eds2 = linspace(.2,1.4,100);
arv = MTADang('data',ang(:,10,12,3),'sampleRate',ang.sampleRate);
a1i=1;eds = linspace(30,60,100);


figure,
subplot(2,3,1);
ind = Trial.stc{'a'};
hist2([arv(ind,a1i),arv2(ind,a2i)],eds,eds2);
subplot(2,3,2);
ind = Trial.stc{'a-w'};
hist2([arv(ind,a1i),arv2(ind,a2i)],eds,eds2);
subplot(2,3,3);
ind = Trial.stc{'w'};
hist2([arv(ind,a1i),arv2(ind,a2i)],eds,eds2);
subplot(2,3,4);
ind = Trial.stc{'a-m'};
hist2([arv(ind,a1i),arv2(ind,a2i)],eds,eds2);
subplot(2,3,5);
ind = Trial.stc{'m'};
hist2([arv(ind,a1i),arv2(ind,a2i)],eds,eds2);
subplot(2,3,6);
ind = Trial.stc{'p'};
hist2([arv(ind,a1i),arv2(ind,a2i)],eds,eds2);
ForAllSubplots('caxis([0,30]),grid(''on'')');
    
ind = Trial.stc{'a&w'};
s = arv(ind,:);
ind = Trial.stc{'a&p'};
n = arv(ind,:);
dprime = (mean(s)-mean(n))./(.5*sqrt(var(s)+var(n)));


