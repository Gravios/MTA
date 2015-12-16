Trial = MTATrial('jg05-20120317');


Trial.load('stc','hand_labeled_rev2');


states = {'walk','rear','turn','pause','groom','sit'};
features = fet_tsne(Trial,20,false);
model = 'MTAC_fet_tsne_REFjg0520120317_NN';



nNeurons = 100;
nIter = 100;

  % Train Model

bhv_nn (Trial,true,states,features,model);




Trial = MTATrial('Ed01-20140707');
StcHL = Trial.load('stc','hand_labeled_rev1');
features = fet20151007(Trial,20,false);
Stc = bhv_nn (Trial,false,states,features,model);


xyz = Trial.load('xyz');
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,states)),'sampleRate',xyz.sampleRate); 

ind = any(shl.data,2)&any(ysm.data,2);

tcm = confmat(shl(ind,:),ysm(ind,:)); % DEP: netlab
precision = round(diag(tcm)'./sum(tcm),4).*100;
sensitivity = round(diag(tcm)./sum(tcm,2),4).*100;
cm = round(tcm./xyz.sampleRate,2);