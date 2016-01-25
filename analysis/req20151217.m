
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');

Trial = MTATrial('Ed03-20140625');
Trial.load('stc','hand_labeled_rev1_Ed');
train = true;
states = {'walk','rear','turn','pause','groom','sit'};
featureSet = 'fet_tsne_rev10';
sampleRate = 12;
ifNormalize = false;
nNeurons = 100;
model = ['MTAC_' featureSet ...
         '_SR_'  num2str(sampleRate) ...
         '_REF_' Trial.name ...
         '_NN_'  num2str(nNeurons) ];

features = feval(featureSet,Trial,sampleRate,ifNormalize);
if ifNormalize, [~,rm,rs] = unity(features); end
 
 % Train Modela
bhv_nn (Trial,train,states,features,model,'nNeurons',nNeurons);


% $$$ Trial = MTATrial('er01-20110719');
% $$$ StcHL = Trial.load('stc','hand_labeled_rev1');
% $$$ 
% $$$ Trial = MTATrial('Ed01-20140707');
% $$$ StcHL = Trial.load('stc','hand_labeled_rev1');

Trial = MTATrial('Ed03-20140624');
%StcHL = Trial.load('stc','hand_labeled_rev1');
Trial.load('stc','hand_labeled_rev1_Ed');
%Trial.load('stc','hand_labeled_rev1_jg');
StcHL = Trial.stc.copy;

Trial = MTATrial('Ed05-20140529');
Trial.load('stc','hand_labeled_rev1_Ed');
StcHL = Trial.stc.copy;

Trial = MTATrial('Ed03-20140625');
Trial.load('stc','hand_labeled_rev1_Ed');
StcHL = Trial.stc.copy;

Trial = MTATrial('Ed05-20140529','all','ont');
Trial.load('stc','hand_labeled_rev1_Ed');
StcHL = Trial.stc.copy;

Trial = MTATrial('Ed01-20140707');
Trial.load('stc','hand_labeled_rev1_jg');
StcHL = Trial.stc.copy;

% It's needed ... I'm not telling why.
xyz = Trial.load('xyz');

% Load and Correct for inter Trial difference
features = feval(featureSet,Trial,sampleRate,ifNormalize);
features.map_to_reference_session(Trial,{'jg05-20120317','all','cof'});
%features.map_to_reference_session(Trial,{'jg05-20120317','all','cof'},true);
if ifNormalize, features.unity([],rm,rs); end

% Create state matrix (N x k) N=samples, k=states, Domain:Boolean
Stc = bhv_nn (Trial,false,states,features,model);
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,states)),'sampleRate',xyz.sampleRate); 

% Select clean periods
aper = Trial.stc{'a'}.cast('TimeSeries');
ind = any(shl.data,2)&any(ysm.data,2)&logical(aper.data);

% Compute confusion matrix, precision and sensitivity
tcm = confmat(shl(ind,:),ysm(ind,:)); % DEP: netlabe
sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
precision = round(diag(tcm)./sum(tcm,2),4).*100;
cm = round(tcm./xyz.sampleRate,2);
acc = round(sum(diag(tcm))/sum(tcm(:)),4);

acm = cat(2,cm,precision);
acm = cat(1,acm,[sensitivity,acc]);


sprintf(strjoin(cat(2,states,{'precision'},{'\n'},...                                              header
                    repmat(cat(2,{'%s'},repmat({'%2.2f'},1,numel(states)+1),{'\n'}),1,numel(states)+1)... cm
                    ),'\t'),acm',states{:},states{:})


f = 4
figure,hold on
eds = linspace(0,170,100);
hs = bar(eds,histc(tfet(:,f),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';
hs = bar(eds,histc(trfet(:,f),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';

Trial = MTATrial('jg05-20120317');
tfet = fet_tsne(Trial,20,false);
ofet = fet_tsne(Trial,20,false);
features = fet20151007(Trial,20,false);

figure,hold on
plot(tfet(:,1))
plot([ofet(:,1),features(:,1)])



%% Stop here this hnn is shit
%% Heirarchical nn labeling 
Trial = MTATrial('jg05-20120317');

Trial.load('stc','hand_labeled_rev2');
features = fet_tsne(Trial,20,false);

nNeurons = 10;
states = {'rear'};
model = 'MTAC_fet_tsne_REFjg0520120317_NN_rear';



% Train Model
      bhv_nn (Trial,true,states,features,model,false,true);

Trial = MTATrial('Ed01-20140707');
StcHL = Trial.load('stc','hand_labeled_rev1');
%features = fet20151007(Trial,20,false);
features = fet_tsne(Trial,20,false);

% Label Trial
Stc = bhv_nn (Trial,false,states,features,model,false,true);


xyz = Trial.load('xyz');
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,{states{1},['gper-' states{1}]})),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,{states{1},'other'})),'sampleRate',xyz.sampleRate); 

ind = any(shl.data,2)&any(ysm.data,2);

tcm = confmat(shl(ind,:),ysm(ind,:)); % DEP: netlab
sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
precision = round(diag(tcm)./sum(tcm,2),4).*100;
cm = round(tcm./xyz.sampleRate,2);



xyz = Trial.load('xyz');

hvar = xyz.copy;
hvar.filter('ButFilter',3,2.4,'low');
hvar = hvar.vel(1,[1,2]);

hvar = features;
hind = 6;

%hind = 1;
%%
eds = linspace(-3,2,100);
figure, hold on
ind = Trial.stc{'a-w-n'};
hs = bar(eds,histc(log10(hvar(ind,hind)),eds),'histc');
hs.FaceColor = 'c';
hs.FaceAlpha = .5;
ind = Trial.stc{'w'};
hs = bar(eds,histc(log10(hvar(ind,hind)),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .5;




f = 2; s ='n' ;
figure,
eds = 50:.5:120;
eds = -3:.05:2;
subplot(311),
bar(eds,hist(rfet(rStc{s},f),eds),'histc'),Lines(mean(rfet(rStc{s},f)),[],'r')
subplot(312),
bar(eds,hist(features(StcHL{s},f),eds),'histc'),Lines(mean(features(StcHL{s},f)),[],'r')
subplot(313),
bar(eds,hist(nfet(StcHL{s},f),eds),'histc'),Lines(mean(nfet(StcHL{s},f)),[],'r')
