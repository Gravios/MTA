function [Co,Po,fs] = Comodugram(Trial,x,varargin)
% function [Co, fs] = Comodugram(Trial,x,varargin)
%
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



% DEFARGS ------------------------------------------------------------------------------------------
parspec = empty_spec;
defargs = struct('states',Trial.stc{'a'},...
                 'mode','mtcsglong',...
                 'defspec',def_spec_parm(x));
[states, mode, defspec] = DefaultArgs(varargin,defargs,'--struct');

dsf = fieldnames(defspec);
for i = 1:length(dsf),parspec.(dsf{i}) = defspec.(dsf{i});end
% if Clip is 1, correlations below 0 will be replaced by 0.
Clip = 0;
%---------------------------------------------------------------------------------------------------


nChannels = size(x,2);
nSamples = size(x,1);

%x.data = nunity(x.data);

% calculate spectrograms
[spex,fs,ts ] = fet_spec(Trial,x,mode,false,[],parspec,[],true);
%[spex,f,ts] = spec(str2func(mode),x.data,parspec,false);
% ??? Normalization ???
thetaState = logical(get(resample(cast([Trial.stc{'theta-groom-sit&gper'}],'TimeSeries'),spex),'data'));
gthind = thetaState;
gthind(thetaState) = nniz(spex(thetaState));

spex.data = bsxfun(@minus,...
                   log(2*(parspec.NW*2-1)*spex.data),...
                   log(mean(spex.data(gthind,:,:,:))));
%spex.data = abs(spex.data);

s = 1;
for state = states,
    state = state{1};
    if ischar(state),
        state = [Trial.stc{state}];
    end
    tspex = spex(logical(get(resample(cast(state,'TimeSeries'),spex),'data')),:,:,:);

    % FIND frequency bins to consider
    nTimeBins = size(tspex,1);
    nFreqBins = size(tspex,2);

    % CALCULATE correlation coefficients
    dataMat = reshape(tspex, [nTimeBins, nFreqBins*nChannels]);
    dataMat = dataMat(nniz(dataMat),:);
    nTimeBins = size(dataMat,1);
    [corrMat,pvalMat]  = corrcoef(dataMat);
    
    if (Clip) corrMat = clip(corrMat, 0, 1); end;

    % produce output array and plot(if required)
    for i=1:nChannels
	for j=i+1:nChannels
            C = corrMat((i-1)*nFreqBins + (1:nFreqBins), (j-1)*nFreqBins + (1:nFreqBins));
            Co(:,:,i,j,s) = C(:,:);
            P = pvalMat((i-1)*nFreqBins + (1:nFreqBins), (j-1)*nFreqBins + (1:nFreqBins));
            Po(:,:,i,j,s) = P(:,:);
        end
    end
    s = s + 1;
end

