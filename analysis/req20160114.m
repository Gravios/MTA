

%Trial = MTATrial('jg05-20120317');
Trial = MTATrial('Ed03-20140624');

StcHL = Trial.load('stc','hand_labeled_rev2_alt');
StcML = Trial.load('stc','NN_multiPN-jg05-20120317.cof.all-RAND_wsb_hand_labeled_rev2-wrnpms');

Stc = bhv_nn (Trial,false,states,features,model);
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,states)),'sampleRate',xyz.sampleRate); 


