
Trial = MTATrial('jg05-20120317');
stcMode = 'hand_labeled_rev2_jg';

% $$$ Trial = MTATrial('Ed03-20140625');
% $$$ stcMode = 'hand_labeled_rev1_Ed';

train = true;
Trial.load('stc',stcMode);
states = {'walk','rear','turn','pause','groom','sit'};
featureSet = 'fet_tsne_rev13';
sampleRate = 12;
ifNormalize = true;
nNeurons = 100;
model = ['MTAC_' featureSet ...
         '_NORM_'  num2str(ifNormalize)...
         '_SR_'  num2str(sampleRate) ...
         '_REF_' Trial.name ...
         '_STC_' stcMode ...
         '_NN_'  num2str(nNeurons) ];

features = feval(featureSet,Trial,sampleRate);
if ifNormalize, [~,rm,rs] = unity(features); end
%features.data =  features.data./5;


 % Train Model
bhv_nn (Trial,train,states,features,model,'nNeurons',nNeurons);

bhv_lgr (Trial,train,states,features);
[StcLGR,dst] = bhv_lgr (Trial,false,states,features);

 
% $$$ Trial = MTATrial('er01-20110719');
% $$$ StcHL = Trial.load('stc','hand_labeled_rev1');
% $$$ 
% $$$ Trial = MTATrial('Ed01-20140707');
% $$$ StcHL = Trial.load('stc','hand_labeled_rev1');

Trial = MTATrial('Ed03-20140624');
%StcHL = Trial.load('stc','hand_labeled_rev1');
Trial.load('stc','hand_labeled_rev1_Ed');
Trial.load('stc','hand_labeled_rev1_jg');
StcHL = Trial.stc.copy;

Trial = MTATrial('Ed03-20140625');
Trial.load('stc','hand_labeled_rev1_Ed');
StcHL = Trial.stc.copy;

% $$$ Trial = MTATrial('er01-20110719');
% $$$ Trial.load('stc','hand_labeled_rev1_jg');
% $$$ StcHL = Trial.stc.copy;
% $$$ 
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2_jg');
StcHL = Trial.stc.copy;
% $$$ 
Trial = MTATrial('Ed05-20140529','all','ont');
Trial.load('stc','hand_labeled_rev1_Ed');
StcHL = Trial.stc.copy;

Trial = MTATrial('Ed01-20140707');
Trial.load('stc','hand_labeled_rev1_jg');
Trial.load('stc','hand_labeled_rev1_Ed');
StcHL = Trial.stc.copy;

% It's needed ... I'm not telling why.
xyz = Trial.load('xyz');

% Load and Correct for inter Trial difference
%features = feval(featureSet,Trial,sampleRate,ifNormalize);
features = feval(featureSet,Trial,sampleRate);
features.map_to_reference_session(Trial,{'jg05-20120317','all','cof'});
%features.map_to_reference_session(Trial,{'jg05-20120317','all','cof'},true);
if ifNormalize, features.unity([],rm,rs); end
%features.data =  features.data./5;


% Create state matrix (N x k) N=samples, k=states, Domain:Boolean
Stc = bhv_nn (Trial,false,states,features,model);
%Stc = bhv_lgr (Trial,false,states,features);
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,states)),'sampleRate',xyz.sampleRate); 
%ysm = MTADxyz('data',double(0<stc2mat(StcLGR,  xyz,states)),'sampleRate',xyz.sampleRate); 

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


                    
                    
%% New heuristics
td = [];
for s = Stc{'w'}.data',
    td(end+1) = sqrt(sum((xyz(s(2),1,:)-xyz(s(1),1,:)).^2,3));
end                    



td = [];
for s = Stc{'r'}.data',
    td(end+1) = median(xyz(s(1):s(2),5,3));
end                    
th = [];
for s = StcHL{'r'}.data',
    th(end+1) = median(xyz(s(1):s(2),5,3));
end                    
figure,hold on,
plot(diff(StcHL{'r'}.data,1,2)/120,th','.b')
plot(diff(Stc{'r'}.data,1,2)/120,td','.r')


z = Trial.xyz.copy;
z.data = xyz(:,5,3);
rl =  z.segs(Stc{'r'}(:,1)-60,120,nan);

z = Trial.xyz.copy;
z.data = xyz(:,5,3);
rh =  z.segs(StcHL{'r'}(:,1)-60,120,nan);

z = Trial.xyz.copy;
z.data = xyz(:,5,3);
rh =  z.segs(StcHL{'r'}(:,1)-120,240,nan);


c = convn(mean(rh(90:150,:),2)',rh(:,2)','full');
figure,plot(c)
figure,plot(diff(c))

figure,sp = [];
sp(end+1) = subplot(211);
imagesc(shl.data')
sp(end+1) = subplot(212);
imagesc(ysm.data')
linkaxes(sp,'xy');

s = 6;
xlim(Stc{'r'}(s,:)+[-60,60])
StcHL{'r'}(s,:)
     


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
