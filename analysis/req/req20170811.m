
Trial = MTATrial.validate('jg05-20120317.cof.all');
features = fet_bref(Trial);
Stc = Trial.load('stc','hand_labeled_rev3_jg');

[smat] = stc2mat(Stc,features,states);
ind = any(smat,2)&nniz(features);



net = layrecnet(1:2:60,10);
[Xs,Xi,Ai,Ts] = preparets(net,num2cell(features(ind,:)'),num2cell(smat(ind,:)'));
net = train(net,Xs,Ts,Xi,Ai);
view(net)
Y = net(Xs,Xi,Ai);
perf = perform(net,Y,Ts)


[X,T] = simpleseries_dataset;
net = layrecnet(1:2,10);
[Xs,Xi,Ai,Ts] = preparets(net,X,T);
net = train(net,Xs,Ts,Xi,Ai);
view(net)
Y = net(Xs,Xi,Ai);
perf = perform(net,Y,Ts)