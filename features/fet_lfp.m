function [flfp,fs,ts] = fet_lfp(Trial,varargin)
%function [fncp,fs,ts] = fet_lfp(Trial,varargin)
%
%rgs(varargin,{'xyz','',2,false,10^11,2^9,[],[],[],[],[1,30]});
%  varargin:
%       sampleRate
%       mode
%       type : MTADlfp or matrix
%       chans
%       overwrite
%       nFFT
%       WinLength
%       nOverlap
%       NW
%       Detrend
%       nTapers
%       FreqRange
%
% Note: The final sampleRate will not be the given sampleRate unless nOverlap is
%       set to -1
[sampleRate,mode,type,chans,overwrite,nFFT,WinLength,nOverlap,NW,Detrend,nTapers,FreqRange] = DefaultArgs(varargin,{'xyz','wcsd','MTADlfp',2,false,10^11,2^9,[],[],'linear',[],[1,30]});

if isempty(nOverlap),
    nOverlap = WinLength*.875;
elseif nOverlap == -1
    nOverlap = WinLength-1;
end


fs = []; ts = [];

flfp = Trial.lfp.copy;

if ~isempty(chans),
    flfp.load(Trial,chans);
end

%% resample lfp
if ischar(sampleRate), 
    sampleRate = Trial.load(sampleRate);
    sampleSize = sampleRate.size(1);
    sampleRate = sampleRate.sampleRate;
    flfp.resample(sampleRate);
elseif isa(sampleRate,'MTAData')
    if sampleRate.isempty,sampleRate.load(Trial);end
    sampleSize = sampleRate.size(1);
    sampleRate = sampleRate.sampleRate;
    flfp.resample(sampleRate);
else
    sampleSize = round(flfp.size(1)/flfp.sampleRate*sampleRate);
    flfp.resample(sampleRate);
end



switch mode
  case 'wcsd'
    % Load ARmodel and Whiten Signal
    try load(fullfile(Trial.path.arm,'fet_rhm.arm.mat')); end
    if exist('ARmodel','var')||overwrite,
        flfp = WhitenSignal(flfp.data,[],[],ARmodel);
    else
        [flfp,ARmodel] = WhitenSignal(flfp.data,[],true);
        save(fullfile(Trial.path.arm,[mfilename '.arm.mat']),'ARmodel');
    end
    
    % Calculate the Cross Spectral Density and Power Spectral Density
    [ys,fs,ts] = mtcsdglong(flfp,nFFT,sampleRate,WinLength,nOverlap,NW,Detrend,nTapers,FreqRange);
    
end

ts = ts+WinLength*.5/sampleRate;
ssr = 1/diff(ts(1:2));
pad = round([ts(1),size(fet,1)./flfp.sampleRate-ts(end)].*ssr)-[1,0];
%pad = round([ts(1),mod(sampleSize-WinLength*.5,WinLength)/sampleRate].*ssr)-[1,0];
szy = size(ys);
ts = cat(1,zeros([pad(1),1]),ts,zeros([pad(2),1]));

switch type
  case 'MTADlfp'
    flfp = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
  case 'matrix'
    flfp = cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)]));
end




