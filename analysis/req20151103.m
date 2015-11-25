Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
states = {'walk','rear','turn','groom','pause','sit'};


x = fet_tsne(Trial,15,true);
t = MTADxyz('data',stc2mat(Trial.stc,x,states),'sampleRate',x.sampleRate);
ind = Trial.stc{'a'};


net = patternnet(40);
view(net);
[net,tr] = train(net,x(ind,:)',~~t(ind,:)');