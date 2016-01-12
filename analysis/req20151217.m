Trial = MTATrial('jg05-20120317');

Trial.load('stc','hand_labeled_rev2');
states = {'walk','rear','turn','pause','groom','sit'};

%featureSet = 'fet_trial20160107';%'fet_tsne_rev3';
featureSet = 'fet_tsne_rev5';
%featureSet = 'fet_tsne_rev4';

sampleRate = 12;
ifNormalize = false;
features = feval(featureSet,Trial,sampleRate,ifNormalize);
%[~,rm,rs] = unity(features);
nNeurons = 100;

model = ['MTAC_' featureSet ...
         '_SR_'  num2str(sampleRate) ...
         '_REF_' Trial.name ...
         '_NN_'  num2str(nNeurons) ];
 

 % Train Model
bhv_nn (Trial,true,states,features,model,'nNeurons',nNeurons);


Trial = MTATrial('er01-20110719');
StcHL = Trial.load('stc','hand_labeled_rev1');

Trial = MTATrial('Ed01-20140707');
StcHL = Trial.load('stc','hand_labeled_rev1');

Trial = MTATrial('Ed03-20140624');
%StcHL = Trial.load('stc','hand_labeled_rev1');
StcHL = Trial.load('stc','hand_labeled_rev2_alt');


%features = fet20151007(Trial,sampleRate,ifNormalize,featureSet,false);
%features = feval(featureSet,Trial,sampleRate,ifNormalize);
features = feval(featureSet,Trial,Trial.xyz.sampleRate,false);
features.normalize(Trial,{'jg05-20120317','all','cof'});
features.resample(sampleRate);
features.unity([],rm,rs);

Stc = bhv_nn (Trial,false,states,features,model);


xyz = Trial.load('xyz');
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,states)),'sampleRate',xyz.sampleRate); 

aper = Trial.stc{'a'}.cast('TimeSeries');

ind = any(shl.data,2)&any(ysm.data,2)&logical(aper.data);

tcm = confmat(shl(ind,:),ysm(ind,:)); % DEP: netlabe
sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
precision = round(diag(tcm)./sum(tcm,2),4).*100;
cm = round(tcm./xyz.sampleRate,2);


%(truePositives + trueNegatives) / totalPopulation
acc = sum(diag(tcm))/sum(tcm(:));

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