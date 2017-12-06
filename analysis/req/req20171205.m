% req20171205
% Description: Determine the respiration frequency distribution during
%              immobility 

Trial = MTATrial.validate('Ed05-20140528.cof.all');
stc   = Trial.load('stc','msnn_ppsvd');
xyz   = preproc_xyz(Trial,xyzProcOpts);
xyz.filter('ButFilter',3,2.4,'low');

% SPEED 
vh = xyz.vel({'spine_lower','head_front'},[1,2]);
vh.resample(ys);
vh.data = log10(abs(vh.data));
vh.data(~nniz(vh(:))) = nan;

% RHYTHMIC HEAD MOTION 
rhm = fet_rhm(Trial);
rhm.resample(xyz);
rhm.data(~nniz(rhm.data(:))) = 1;

% RESPIRATION 
ncp = fet_ncp(Trial,rhm,'mta',ncpChannel);
ncp.filter('RectFilter',5,5);
ncp.filter('ButFilter',3,[0.5,14],'bandpass');

% FIND timepoints of inhalation
ncpPeaks = LocalMinima(ncp.data,8,-500);
ncpPeaks(ncpPeaks>size(ncp,1)) = [];
% COMPUTE instantaneous respiration rate (Hz)
ncpFreq = 1./((ncpPeaks-circshift(ncpPeaks,1))/ncp.sampleRate);
% SMOOTH respiration rate (Hz)
ncpFreq = median(GetSegs(ncpFreq,circshift([1:size(ncpFreq,1)]',5),11));


hfig = figure();
tind = stc{'t',ncp.sampleRate};
pind = stc{'p',ncp.sampleRate};

subplot(131);% Pause
ind = WithinRanges(ncpPeaks,pind.data);
plot(vh(ncpPeaks(ind),2),ncpFreq(ind)+randn([1,sum(ind)]),'.')
xlim([-3,3]);  ylim([0,16]);  grid('on');

subplot(132);% Pause & Theta
ind = pind&tind;
ind = WithinRanges(ncpPeaks,ind.data);
plot(vh(ncpPeaks(ind),2),ncpFreq(ind)+randn([1,sum(ind)]),'.')
xlim([-3,3]);  ylim([0,16]);  grid('on');

subplot(133);% Pause - Theta
ind = pind-tind;
ind = WithinRanges(ncpPeaks,ind.data);
plot(vh(ncpPeaks(ind),2),ncpFreq(ind)+randn([1,sum(ind)]),'.')
xlim([-3,3]);  ylim([0,16]);  grid('on');
