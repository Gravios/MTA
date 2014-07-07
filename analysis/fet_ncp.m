function fncp = fet_ncp(Trial,varargin)
[sampleRate,mode,chans] = DefaultArgs(varargin,{'lfp','',2});

if ischar(sampleRate), sampleRate = Trial.(sampleRate).sampleRate;end

fncp = Trial.lfp.copy;

if ~isempty(chans),
    fncp.load(Trial,chans);
end
Trial.load('xyz');
fncp.resample(Trial.xyz);

switch mode
case 'spectral'
[ys,fs,ts] = mtcsdglong(nang,2^8,fncp.sampleRate,2^7,2^7-1,[],'linear',[],[1,30]);
fncp = zeros(xyz.size(1),size(ys,2));
fncp((2^6+1):(size(ys)+2^6),:) = ys;
fncp = MTADlfp('data',fncp,'sampleRate',xyz.sampleRate);
case 'wspectral'
load(fullfile(Trial.path.MTAPath,'fet_rhm.arm'),'-mat');
[nang] = WhitenSignal(fncp.data,[],[],ARmodel);
[ys,fs,ts] = mtcsdglong(nang,2^8,fncp.sampleRate,2^7,2^7-1,[],'linear',[],[1,30]);
fncp = zeros(xyz.size(1),size(ys,2));
fncp((2^6+1):(size(ys)+2^6),:) = ys;
fncp = MTADlfp('data',fncp,'sampleRate',xyz.sampleRate);
case 'default'
fncp = MTADlfp('data',fncp.data,'sampleRate',xyz.sampleRate);
otherwise
fncp = fncp.data;
end




