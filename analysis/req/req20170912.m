% NAME :    req20170912
% PURPOSE : recursive neural network
%           
%


% Training set
Trial = MTATrial.validate('jg05-20120317.cof.all');
StcHL = Trial.load('stc','hand_labeled_rev3_jg');

% Testing set
Trial = MTATrial.validate('Ed05-20140529.ont.all');
StcHL = Trial.load('stc','hand_labeled_rev1_Ed');



features = fet_bref_rev14(Trial);
features.map_to_reference_session(Trial,'jg05-20120317.cof.all');

xyz = Trial.load('xyz');



states = {'walk','rear','turn','pause','groom','sit'};
keys   = {'w','r','n','p','m','s'};
numStates = numel(states);

shl = MTADxyz('data',double(0<stc2mat(StcHL,features,states)),'sampleRate',features.sampleRate);

aper = StcHL{'a'};
aper.cast('TimeSeries');
aper.resample(features);
nind = nniz(features) & logical(aper.data);

%net = layrecnet([1:10],[10,10,10]);
[Xs,Xi,Ai,Ts] = preparets(net,con2seq(features(nind,:)'),con2seq(shl(nind,:)'));

% $$$ net = train(net,Xs,Ts,Xi,Ai);
% $$$ 
% $$$ save(fullfile(Trial.spath,[Trial.filebase,'-recnet.mat']),'net');


%view(net)
Y = zeros(size(shl));
g = net(Xs,Xi,Ai);
gind = find(nind);
Y(gind(1:end-10),:) = circshift(reshape(cell2mat(g),6,[])',10);
% $$$ perf = perform(net,Y,Ts)

labs = features.copy;
labs.data = Y;
labs.resample(xyz);

figure();
subplot(211);plot(Y);
subplot(212);plotSTC(StcHL,features.sampleRate);
linkaxes(findobj(gcf(),'Type','axes'),'x');





% CHECK Labeling stats



[~,maxState] = max(labs.data,[],2);
% UNLABEL regions where labels cannot exist 
maxState(~nniz(features),:) = 0;


% CREATE state collection with simple mode name
StcNN = StcHL.copy();
StcNN.states = {};

% APPEND network output labeles to state collection
for i = 1:numel(states)
    sts = ThreshCross(maxState==i,0.5,1);
    try sts = bsxfun(@plus,sts,[1,0]); end
    StcNN.addState(Trial.spath,Trial.filebase,sts,xyz.sampleRate,...
                       xyz.sync.copy,xyz.origin,states{i},keys{i},'TimePeriods');
end




shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);

% CONVERT state collection into state matrix for inter labeler comparison
ysm = MTADxyz('data',double(0<stc2mat(StcNN,xyz)),'sampleRate',xyz.sampleRate);
ysm.data(~nniz(xyz),:) = 0; 


% INITIALIZE inter labeler stats
labelingStats = struct('confusionMatrix',zeros([numStates,numStates]),...
                       'precision',     zeros([1,numStates]),...
                       'sensitivity',   zeros([1,numStates]),...
                       'accuracy',      0 );


% COMPUTE inter labeler stats if valid comparative Stc was loaded above

ind = any(shl.data,2)&any(ysm.data,2);
tcm = confmat(shl(ind,:),ysm(ind,:));
labelingStats.confusionMatrix = round(tcm./xyz.sampleRate,2);
labelingStats.precision = round(diag(tcm)'./sum(tcm,2)',4).*100;
labelingStats.sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
labelingStats.accuracy = sum(diag(tcm))/sum(tcm(:));

