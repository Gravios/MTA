
% 

Trial = MTATrial.validate('ER06-20130613.cof.all');

Trial = MTATrial.validate('Ed03-20140624.cof.all');

Trial = MTATrial.validate('jg05-20120317.cof.all');

%compute_session_rhm_distribution(Trial,[],'loc','NN0317R')




parspec = empty_spec;
xyz = Trial.load('xyz');
varargin = {};
[sampleRate,mode,wsig,defspec,overwrite,newSR,type] = DefaultArgs(varargin,{'xyz','mtchglong',1,def_spec_parm(xyz),false,[],'mta'});

fs = []; ts = [];
% create a ridgid body model
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
% find the center of mass of the model
hcom = xyz.com(rb);
% add coordinates of the model's center of mass to the xyz object
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(hcom,3,[2]./(Trial.xyz.sampleRate/2),'low'));
% if xyz sampling rat e is greater than 120 Hz then resample it to 120 Hz
if ischar(sampleRate), sampleRate = Trial.(sampleRate).sampleRate;end
if isa(sampleRate,'MTAData'),
    xyz.resample(sampleRate);    
elseif sampleRate > 120, 
    xyz.resample(120); 
end
xyz.filter('ButFilter',3,55,'low');
% convert all iter-marker coordinates from cartesian to spherical
ang = create(MTADang,Trial,xyz);
% Create a feature object based on the xyz object
fet = Trial.xyz.copy;
bang = ButFilter(ang(:,'head_back','fhcom',3),3,[.5,50]./(Trial.ang.sampleRate/2),'bandpass');
fet.data = [0;ButFilter(diff(bang),3,[.5,50]/(ang.sampleRate/2),'bandpass')];
% update fields for modified options passed through function handle 
defspec.FreqRange = [1,30];
dsf = fieldnames(defspec);
for i = 1:length(dsf),parspec.(dsf{i}) = defspec.(dsf{i});end
% create placeholder for data 
data = zeros(fet.size);
% whiten signal with a 10 second window
[data(nniz(fet.data),:),ARmodel] = WhitenSignal(fet.data(nniz(fet.data),:),[],true);
[data(nniz(fet.data),:),ARmodel] = WhitenSignal(fet.data(nniz(fet.data),:),round(fet.sampleRate.*30),true,[],5);
% compute spectrogram of feature
[ys,fs,ts] = spec(str2func(mode),data,parspec);
% Modify time stamps and spec; add padding (0's)
ts = ts+(parspec.WinLength/2)/fet.sampleRate;
ssr = 1/diff(ts(1:2));
pad = round([ts(1),mod(fet.size(1)-round(parspec.WinLength/2),parspec.WinLength)/fet.sampleRate].*ssr)-[1,0];
szy = size(ys);
rhm = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
ts = cat(2,([1:pad(1)]./ssr),ts',([pad(1)+size(ts,1)]+[1:pad(2)])./ssr)';


figure,
sp(1) = subplot(211);imagesc(ts,fs,log10(rhm.data)');
axis('xy');
colormap('jet');
caxis([-7,-3])
sp(2) = subplot(212);imagesc(ts,fs,log10(rhm.data)');
axis('xy');
colormap('jet');
caxis([-7,-3])
linkaxes(sp,'xy');

% log10 and unity
% $$$ rhm.data  = log10(rhm.data);
% $$$ rhm.data(rhm<-9) = nan;
% $$$ rhm.data(nniz(rhm.data))=nan;


figure,imagesc(ts,fs,bsxfun(@ldivide,log10(rhm.data)', ...
                         sum(log10(rhm.data)'))),axis xy,colormap jet

rpra = log10(1./bsxfun(@rdivide,log10(rhm.data)',sum(log10(rhm.data)')))';
figure,imagesc(1:size(rhm,1),fs,rpra),axis xy,colormap jet
caxis([1.8,1.9])

mrpra = rhm.copy;
mrpra.data = median(rpra(:,6<fs&fs<12),2);
mrp = rhm.copy;
mrp.data = median(rhm.data(:,6<fs&fs<12),2);



ind = Trial.stc{'x+p'};
ind.cast('TimeSeries',mrp);
ind.data = ind.data&nniz(mrp.data)
[State, hmm, decode] = gausshmm(mrp(ind.data),3);

headPitchState = zeros([size(xyz,1),2]);
[~,stsInd] = max(cat(1,hmm.state.Mu));
headPitchState(ind.data==1,1) = State'==stsInd;


figure,imagesc(1:size(rhm,1),fs,log10(rhm.data)'),
axis xy,
colormap jet
caxis([-7,-3])
hold on
plot(headPitchState(:,1))



sesList = get_session_list('jg05');

trialList = arrayfun(@(x) MTATrial.validate(x),sesList,'UniformOutput',false);

xyzList = cellfun(@(x) x.load('xyz'),trialList,'UniformOutput',false);


rhmList = cellfun(@(x) fet_rhm(x,[],'mtchglong',true),trialList,'UniformOutput',false);


% load stc NN0317R 
cellfun(@(x) x.load('stc','NN0317R'),trialList,'UniformOutput',false);


ind = Trial.stc{'x+p'};
ind.cast('TimeSeries',mrp);
ind.data = ind.data&nniz(mrp.data)
[State, hmm, decode] = gausshmm(mrp(ind.data),3);

headPitchState = zeros([size(xyz,1),2]);
[~,stsInd] = max(cat(1,hmm.state.Mu));
headPitchState(ind.data==1,1) = State'==stsInd;


figure,imagesc(1:size(rhm,1),fs,log10(rhm.data)'),
axis xy,
colormap jet
caxis([-7,-3])
hold on
plot(headPitchState(:,1))
