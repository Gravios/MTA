

Trial = MTATrial('jg05-20120317');
StcHL = Trial.load('stc','hand_labeled_rev1');
%StcML = Trial.load('stc','hand_labeled_rev2');


%Trial = MTATrial('Ed03-20140624');
%StcHL = Trial.load('stc','hand_labeled_rev2_alt');
%StcML = Trial.load('stc','NN_multiPN-jg05-20120317.cof.all-RAND_wsb_hand_labeled_rev2-wrnpms');



xyz = Trial.load('xyz');
states = {'walk','rear','turn','pause','groom','sit'};



% Create state matrix (N x k) N=samples, k=states, Domain:Boolean
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(0<stc2mat(StcML,  xyz,states)),'sampleRate',xyz.sampleRate); 





tau = .1; % second symetric mask at state change points
[~,mind] = max(shl.data,[],2);
mask = MTADepoch('data',        abs(diff(mind)),...
                 'sampleRate',  shl.sampleRate,...
                 'syncPeriods', Trial.sync.copy,...
                 'syncOrigin',  Trial.sync.data(1),...
                 'type',       'TimeSeries');
mask.cast('TimePeriods');
mask = mask+[-tau,tau];
mask.cast('TimeSeries');


% Select clean periods
aper = Trial.stc{'a'}.cast('TimeSeries');
ind = any(shl.data,2)&any(ysm.data,2)&logical(aper.data)&logical(mask.data);


% Compute confusion matrix, precision and sensitivity
tcm = confmat(shl(ind,:),ysm(ind,:)); % DEP: netlabe
sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
precision = round(diag(tcm)./sum(tcm,2),4).*100;
cm = round(tcm./xyz.sampleRate,2);
acc = round(sum(diag(tcm))/sum(tcm(:)),4);



