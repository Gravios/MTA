% [Co f] = Comodugram(x, nFFT, SampleRate, FreqRange, WinLength, NW, Detrend)
%
% Takes an input sequence x and does a multi-pane
% plot showing correlated changes in power in frequency
% bands 
%
% nFFT and SampleRate are just like for specgram() function
%
% FreqRange = [fLow fHigh] allows you to view only a certain 
% frequency range.  Specify it in Hz.If 2 rows - freq.range is different 
% for each signal
% NB x is of the form x(Time, Channel)   !!!IT CHANGED!!!
%
% NW is an argument for the multitaper method, as is Detrend
%
% optional output argument nTimeBins gives number of points in the regression
% - so you can do significance testing.
%
% Example: [Co nTimeBins] = Comodugram(eeg([1,5],:)',1024,1250,[1 100]);
%
%         Where: eeg is the filename loaded with bload, 1 and 5 are the channels
%                you've selected to compare (the : here is where you would put
%                a time restriction, such as a REM period). 1024 is nFFT, 1250
%                is the sampling rate, and [1 100] is the frequency range you
%                wish to compare.
%
%          NOTE: you MUST transpose the first input statement (x)

function [Co, f] = Comodugram(Trial,x,varargin)

parspec = empty_spec;

defargs.states   = Trial.stc{'a'};
defargs.mode    = 'mtcsglong';
defargs.defspec = def_spec_parm(x);

[states,mode, defspec] = DefaultArgs(varargin,defargs,'--struct');


dsf = fieldnames(defspec);
for i = 1:length(dsf),parspec.(dsf{i}) = defspec.(dsf{i});end


Clip = 0;
% if Clip is 1, correlations below 0 will be replaced by 0.

nChannels = size(x,2);
nSamples = size(x,1);

% calculate spectrograms
[spex,f,ts] = spec(str2func(mode),x.data,parspec,false);
spex = log(2*(parspec.NW*2-1)*spex)  - log(repmat(mean(spex), [size(spex,1),1, 1]));
spex = abs(spex);
% Modify time stamps and spec; add padding (0's)
ts = ts+(parspec.WinLength/2)/x.sampleRate;
ssr = 1/diff(ts(1:2));
pad = round([ts(1),mod(x.size(1)-round(parspec.WinLength/2),parspec.WinLength)/x.sampleRate].*ssr)-[1,0];
szy = size(spex);
spex = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),spex,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
ts = cat(2,([1:pad(1)]./ssr),ts',([pad(1)+size(ts,1)]+[1:pad(2)])./ssr)';

s = 1;
for state = states,
    state = state{1};
    if ischar(state),
        state = Trial.stc{state};
    end
    state.cast('TimeSeries');
    state.resample(spex);

    tspex = spex(state.data==1,:,:,:);

    % find frequency bins to consider

    nTimeBins = size(tspex,1);
    nFreqBins = size(tspex,2);

    % calculate correlation coefficients
    DataMat = reshape(tspex, [nTimeBins, nFreqBins*nChannels]);
    
    CorrMat = corrcoef(DataMat);
    if (Clip) CorrMat = clip(CorrMat, 0, 1); end;

    % produce output array and plot(if required)
    C = zeros(nFreqBins,nFreqBins);

    for i=1:nChannels
	for j=1:nChannels
            
            C = CorrMat((i-1)*nFreqBins + (1:nFreqBins), (j-1)*nFreqBins + (1:nFreqBins));
            
            if (nargout<1)
                subplot(nChannels, nChannels, j + (i-1) * nChannels);
                imagesc(f, f, C(:,:));
                set(gca,'ydir','norm');
                drawnow;
            else
                Co(:,:,i,j,s) = C(:,:);
            end
        end
    end
    
    s = s + 1;

end

% 
% if nargout >=1
% 	Co = C(:,:);
% end