function out = DetectOscillations(x, FreqRange,varargin)
% function out = DetectOscillations(x, FreqRange, MinCycles, SampleRate, Thr, MaxDuration)
% detects oscillatory bursts in FreqRange of MinCycles above Threshold
% by default MinCycles is 5 cycles of the center of FreqRange
% Threshold is computed as Thr percentile of the RMS distribution
% or from Thr * std of RMS if Thr is below 20 (e.g. 3 or 5 would be treated
% as std)
% if Thr is a scalar, this is used to detect peaks of the RMS and 1/2
% of the median peaks amplitude will be used for duration determination.
% if Thr is a vector, then first element is for the peak detection,
% second - for the duration.
% MaxDuration (msec)- maximal duration above threshold 

[ MinCycles, SampleRate, Thr, MaxDuration ] = DefaultArgs(varargin,{5, 1250, 90, []});

fx = ButFilter(x,2,FreqRange/SampleRate*2,'bandpass');

%tic; amp = abs(hilbert(fx));toc % too slow and not optimal, needs
%smoothing

MinDuration = SampleRate/mean(FreqRange)*MinCycles;

WinLen = round(SampleRate/FreqRange(1)*5); 
Window = gausswin(WinLen, 5);
Window = Window/sum(Window);

amp = sqrt(filtfilt(Window,1,fx.^2));
zamp = zscore(amp);

if any(Thr>20)
    Thresh = prctile(amp,Thr);
else
    Thresh = Thr*std(amp)+mean(amp);
end

PowerPeaks = LocalMinima(-amp, 1.5*MinDuration, -Thresh(1));

if length(Thresh)==1
    Tresh(2) = 0.5*median(abs(PowerPeaks)); 
end

OscPeriods = ThreshCross(amp, Thresh(2), 1*MinDuration);

% filter our too long oscillation periods 
if ~isempty(MaxDuration)
    OscPeriods = OscPeriods(diff(OscPeriods,1,2)<MaxDuration*1000/SampleRate,:);
end

%select Peaks within OscPeriods
WithinPeriods = WithinRanges(PowerPeaks, OscPeriods);
PowerPeaks =  PowerPeaks(WithinPeriods);

out.t = PowerPeaks;
out.pow = amp(PowerPeaks);
out.zpow = zamp(PowerPeaks);

pp = WithinRangesLongA(PowerPeaks, OscPeriods,[1:size(OscPeriods,1)]',[],[],1);
[~, OscPerInd] = find(pp>0);
out.per = OscPeriods(OscPerInd,:);
out.len = diff(out.per,1,2)/SampleRate*1000;

lm = LocalMinima(fx, 1, 0);
    
out.troughs = lm(WithinRanges(lm, out.per));


%select Oscillaiton Periods that have Power Peak inside - WILL HAVE TO WAIT
if 0
    %[xxl xil yil] = NearestNeighbour(OscPeriods(:,1),PowerPeaks, 'left', MinDuration);
    %[xxr xir yir] = NearestNeighbour(OscPeriods(:,2),PowerPeaks, 'right', MinDuration);
    %[WithPeaks il ir] = intersect(xil, xir);
    
    pp = WithinRangesLongA(PowerPeaks, OscPeriods,[1:size(OscPeriods,1)]',[],[],1);
    PeaksInPeriod = full(sum(pp,1))';
    PeriodsWithPeaks  = full(sum(pp,2));
    
    OscPeriods1 = OscPeriods(PeriodsWithPeaks,:);
    PowerPeaks1 = PowerPeaks
   
    
    PeriodLen = diff(OscPeriods,1,2);
    
    out.OscPeriods = OscPeriods;
end
if 0
clf
plot([1:length(amp)]/SampleRate*1000,([fx amp]));
hold on
Lines([],Thresh, 'r');
Lines(PowerPeaks/SampleRate*1000,[],'k');
Lines(OscPeriods(:)/SampleRate*1000,[],'c');
% Lines(PowerPeaks(yil)/SampleRate*1000,[],'r');
% Lines(PowerPeaks(yir)/SampleRate*1000,[],'r');
% Lines(OscPeriods(:,1)/SampleRate*1000,[],'g');
% Lines(OscPeriods(:,2)/SampleRate*1000,[],'g');

end