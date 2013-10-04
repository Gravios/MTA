% [ ThetaPhase ThetaAmp TotPhase] = ThetaPhase(Eeg, FreqRange, FilterOrd, Ripple, SampleRate)
%
% takes a 1-channel Eeg file (assumes 1250 Hz) and produces
% instantaneous theta phase and amplitude.  TotPhase is
% unwrapped phase.
%
% Theta is filtered with filtfilt and a cheby2 filter with parameters
% FreqRange, FilterOrd, Ripple defautls [4 10], 4, 20.
% good values for gamma are [40 100], 8, 20.
%
% if no args are provided, will plot some diagnostics
function [ThetaPhase, ThetaAmp, TotPhase, Eegf] = myThetaPhase(Eeg, varargin)

[FreqRange, FilterOrd, Ripple, SampleRate]=DefaultArgs(varargin,{[4 10], 4, 20, 1250});


if min(size(Eeg)>1)
	error('Eeg should be 1 channel only!');
elseif size(Eeg,1)==1
    Eeg = Eeg(:);
end
NFreq = SampleRate/2;
%[b a] = Scheby2(FilterOrd, Ripple, FreqRange/NFreq);
%Eegf = Sfiltfilt(b,a,Eeg);
Eegf = ButFilter(Eeg,2,FreqRange/NFreq,'bandpass');
% remove constant term to avoid bias
Eegf = Eegf - mean(Eegf);
if nargout>0, clear Eeg; end;
Hilb = Shilbert(Eegf);
if nargout>0, clear Eegf; end;
ThetaPhase = angle(Hilb);

if nargout>=2
	ThetaAmp = abs(Hilb);
end
if nargout>=3
	TotPhase = unwrap(ThetaPhase);
end

if nargout==0
    subplot(3,1,1);
    [h w s] = freqz(b, a, 2048, 1250);
    plot(w,abs(h));
    grid on
    title('frequncy response of filter');
    
    subplot(3,1,2)
    xr = (1:length(Eeg))*1000/1250;
    plot(xr, [Eeg, Eegf]);
    title('eeg');
    legend('raw', 'filtered');
    
    subplot(3,1,3);
    plot(xr, [Eegf, ThetaPhase*std(Eegf)]);
    clear ThetaPhase
    title('extracted phase');
    legend('filtered wave', 'phase');
end

