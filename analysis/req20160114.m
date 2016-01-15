

%Trial = MTATrial('jg05-20120317');
Trial = MTATrial('Ed03-20140624');
states = {'walk','rear','turn','pause','groom','sit'};

% Load the state Collections
StcHL = Trial.load('stc','hand_labeled_rev2_alt');
StcML = Trial.load('stc','NN_multiPN-jg05-20120317.cof.all-RAND_wsb_hand_labeled_rev2-wrnpms');

% Create state matrix (N x k) N=samples, k=states, Domain:Boolean
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(0<stc2mat(StcML,  xyz,states)),'sampleRate',xyz.sampleRate); 




% Select clean periods
aper = Trial.stc{'a'}.cast('TimeSeries');
ind = any(shl.data,2)&any(ysm.data,2)&logical(aper.data);

tau = .1; % second symetric mask at state change points
tind = zeros(size(ind));

[~,mind] = max(shl.data,[],2);
mask = MTADepoch('data',abs(diff(maskS)),'sampleRate',shl.sampleRate,'type','TimeSeries');
mask.cast('TimePeriods');



% Compute confusion matrix, precision and sensitivity
tcm = confmat(shl(ind,:),ysm(ind,:)); % DEP: netlabe
sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
precision = round(diag(tcm)./sum(tcm,2),4).*100;
cm = round(tcm./xyz.sampleRate,2);
acc = round(sum(diag(tcm))/sum(tcm(:)),4);



