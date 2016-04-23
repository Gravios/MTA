
% Testing Fast Artificial Neural Network (FANN) c library


% dbstop in bhv_nn.m at 146 just before the patternnet

fetMat = feature(ind,:);
stsMat = ~~smat(ind,:);

fid = fopen('/storage/gravio/data/fann/train_data.data','w+');
% print header values
fprintf(fid,'%i %i %i \n',size(fetMat,1),size(fetMat,2),size(stsMat,2));
% print over loop input answer pairs
for ind = 1:size(fetMat,1),
    fprintf(fid,[repmat('%12.12f ',[1,size(fetMat,2)]),'\n'],fetMat(ind,:));
    fprintf(fid,[repmat('%i ',[1,size(stsMat,2)]),'\n'],stsMat(ind,:));
end

fclose(fid);



% dbstop at 145 in bhv_nn_multi_patternnet.m

fid = fopen('/storage/gravio/data/fann/test_data.data','w+');

% print over loop input
for ind = 1:size(features,1),
    fprintf(fid,[repmat('%12.12f ',[1,size(features,2)]),'\n'],features(ind,:));
end

fclose(fid);


mstc = load('/storage/gravio/data/fann/results.data','-ASCII');

Trial = MTATrial.validate(s);
xyz = Trial.load('xyz');
xyz.resample(12);
labelingEpochs.resample(xyz);
StcHL = Trial.stc.copy;
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);

d_state = MTADxyz('data',mstc,'sampleRate',xyz.sampleRate);
[~,maxState] = max(d_state.data,[],2);
maxState(~nniz(xyz),:) = 0;

ysm = MTADxyz('data',zeros([shl.size]),'sampleRate',xyz.sampleRate); 
ysm.data = ysm.data';
ysm.data([1:size(ysm,1):size(ysm,2).*size(ysm,1)]+maxState'-1) = 1;
ysm.data = ysm.data';
ind = any(shl.data,2)&any(ysm.data,2)&labelingEpochs.data;
tcm = confmat(shl(ind&labelingEpochs,:),ysm(ind&labelingEpochs,:)) % #DEP: netlab

precision = round(diag(tcm)./sum(tcm,2),4).*100
sensitivity = round(diag(tcm)'./sum(tcm),4).*100
accuracy = sum(diag(tcm))/sum(tcm(:)).*100



rlist = SessionList('training_hand_labeled');
slist = {'hand_labeled_jg';'hand_labeled_Ed'};
fetSet  = 'fet_tsne_rev15';
sampleRate = 12;
nNeurons = 100;
nIter = 2;
states = {'walk','rear','turn','pause','groom','sit'};
rndMethod = 'WSBNT';
norm = true;
mref = true;


%for s = rlist
s = rlist(1);
clear('mod');
mod.states     = states;
mod.stcMode    = s.stcMode;
mod.featureSet = fetSet;
mod.sampleRate = sampleRate;
mod.nNeurons   = nNeurons;
mod.nIter      = nIter;
mod.randomizationMethod = rndMethod;
mod.normalize = norm;
argin = struct2varargin(mod);
bfo = cell([1,6]);
[bfo{:}] = bhv_fann_multi_patternnet(MTATrial.validate(s),argin{:});
%end


s = rlist(2);
clear('mod');
mod.states     = states;
mod.stcMode    = s.stcMode;
mod.featureSet = fetSet;
mod.sampleRate = sampleRate;
mod.nNeurons   = nNeurons;
mod.nIter      = nIter;
mod.randomizationMethod = rndMethod;
mod.normalize = norm;
mod.model = bfo{5};
argin = struct2varargin(mod);
tfo = cell([1,6]);
[tfo{:}] = bhv_fann_multi_patternnet(MTATrial.validate(s),argin{:});

