function out = DetectStructuredOscilations(Trial,varargin)

[lfp,FreqRange,MinCycles,MaxDuration,Thr,ModelTemplate] = DefaultArgs(varargin,{Trial.lfp.copy,[160,220],5,.1,5,'default'});

%%Diagnostic Vars
% $$$ Trial = MTATrial('jg05-20120310');
% $$$ lfp = Trial.lfp.copy;
% $$$ lfp.load(Trial,[65:72]);
% $$$ FreqRange = [160,220];
% $$$ MinCycles = 5;
% $$$ %MaxDuration = [];
% $$$ MaxDuration = 0.1;
% $$$ Thr = 5;
% $$$ ModelTemplate = 'default';
%%%%%%%%%%%%%%%%%%%%%%%%%




xyz = Trial.load('xyz');
fx = ButFilter(lfp.data,2,FreqRange/lfp.sampleRate*2,'bandpass');

MinDuration = xyz.sampleRate/mean(FreqRange)*MinCycles;

WinLen = round(lfp.sampleRate/FreqRange(1)*5); 
Window = gausswin(WinLen, 5);
Window = Window/sum(Window);
amp = sqrt(filtfilt(Window,1,fx.^2));
amp = MTADlfp('data',amp,'sampleRate',lfp.sampleRate);
amp.resample(xyz);


switch ModelTemplate
  case 'default'
    ModelTemplate = load(fullfile(Trial.path.MTAPath,'DetectStructedOscilations.default.mdl'));

  case 'ripples'
     
    % Make interative gui for making models
    %ripper = [81196,81208];
    %msave(fullfile(Trial.path.MTAPath,'DetectStructedOscilations.default.mdl'),urpexp-bsh)
    
    %rpexp = amp(ripper(1):ripper(2),:);
    rpexp = load(fullfile(Trial.path.MTAPath,'DetectStructedOscilations.default.mdl'));

    for i= 1:lfp.size(2),rpmean(i) = mean(amp(amp(:,i)>2,i));end
    for i= 1:lfp.size(2),rpstd(i) = std(amp(amp(:,i)>2,i));end
    urpow = (amp.data-rpmean(1))./rpstd(1);
    urpexp = rpexp;
    %urpexp = urpexp/sum(urpexp(:));

    [~,mchan] = max(urpexp,[],2)
    mchan = round(mean(mchan));
    
    turpexp = urpexp;
    fsh = .0002;
    bshift = [];
    for i = 1:80,
        turpexp = turpexp-fsh;
        crp = conv2(urpow,turpexp,'same');
        bshift(end+1) = skewness(zscore(crp(:,mchan).*urpow(:,mchan)));    
    end
    %figure,plot([1:80]*fsh,bshift,'.')
    [~,bsh] = max(bshift);
    bsh = bsh*fsh;
    crp = conv2(urpow,urpexp-bsh,'same');
    amp.data = crp(:,mchan).*urpow(:,mchan);
  case 'STDMASK'

  otherwise % Load custom model
    
end

%figure,plot(unity(ButFilter(crp(:,4),3,[.1,30]/(amp.sampleRate/2),'bandpass').*max(urpow,[],2)),'b'),
%hold on,plot(max(urpow,[],2),'r');


zamp = amp.copy;
zamp.data = zscore(amp.data);


if any(Thr>20)
    Thresh = prctile(amp.data,Thr);
else
    Thresh = Thr*std(amp.data)+mean(amp.data);
end
[PowerPeaksInd,PowerPeaks] = LocalMinima(-amp.data, ceil(1.5*MinDuration), -Thresh(1));
PowerPeaks = -PowerPeaks;


if length(Thresh)==1
    Thresh(2) = 0.5*median(abs(PowerPeaks)); 
end

OscPeriods = ThreshCross(amp, Thresh(2), 1*MinDuration);

% filter our too long oscillation periods 
if ~isempty(MaxDuration)
    OscPeriods = OscPeriods(diff(OscPeriods,1,2)<MaxDuration*xyz.sampleRate,:);
end

%select Peaks within OscPeriods
WithinPeriods = WithinRanges(PowerPeaksInd, OscPeriods);
PowerPeaks =  PowerPeaks(WithinPeriods);
PowerPeaksInd =  PowerPeaksInd(WithinPeriods);
out.t = PowerPeaksInd;
out.pow = amp(PowerPeaksInd);
out.zpow = zamp(PowerPeaksInd);
out.sampleRate = xyz.sampleRate;

%pp = WithinRangesLongA(PowerPeaksInd, OscPeriods,[1:size(OscPeriods,1)]',[],[],1);
%[~, OscPerInd] = find(pp>0);
out.per = OscPeriods;%(OscPerInd,:);
out.len = diff(out.per,1,2)/xyz.sampleRate*1000;

bootInd =nan;
while ~isempty(bootInd)
    bootInd = find(out.t(1:size(out.per,1)) -out.per(:,1)<0,1,'first');
    out.t(bootInd)=[];
    out.pow(bootInd)=[];
    out.zpow(bootInd)=[];
end

%lm = LocalMinima(fx, 1, 0);
    
%out.troughs = lm(WithinRanges(lm, out.per));


%select Oscillaiton Periods that have Power Peak inside - WILL HAVE TO WAIT
% $$$ if 0
% $$$     %[xxl xil yil] = NearestNeighbour(OscPeriods(:,1),PowerPeaks, 'left', MinDuration); 
% $$$     %[xxr xir yir] = NearestNeighbour(OscPeriods(:,2),PowerPeaks, 'right', MinDuration);
% $$$     %[WithPeaks il ir] = intersect(xil, xir);
% $$$     
% $$$     pp = WithinRangesLongA(PowerPeaks, OscPeriods,[1:size(OscPeriods,1)]',[],[],1);
% $$$     PeaksInPeriod = full(sum(pp,1))';
% $$$     PeriodsWithPeaks  = full(sum(pp,2));
% $$$     
% $$$     OscPeriods1 = OscPeriods(PeriodsWithPeaks,:);
% $$$     PowerPeaks1 = PowerPeaks
% $$$    
% $$$     
% $$$     PeriodLen = diff(OscPeriods,1,2);
% $$$     
% $$$     out.OscPeriods = OscPeriods;
% $$$ end
% $$$ if 0
% $$$ clf
% $$$ plot([1:length(amp)]/SampleRate*1000,([fx amp]));
% $$$ hold on
% $$$ Lines([],Thresh, 'r');
% $$$ Lines(PowerPeaks/SampleRate*1000,[],'k');
% $$$ Lines(OscPeriods(:)/SampleRate*1000,[],'c');
% $$$ % Lines(PowerPeaks(yil)/SampleRate*1000,[],'r');
% $$$ % Lines(PowerPeaks(yir)/SampleRate*1000,[],'r');
% $$$ % Lines(OscPeriods(:,1)/SampleRate*1000,[],'g');
% $$$ % Lines(OscPeriods(:,2)/SampleRate*1000,[],'g');
% $$$ end





