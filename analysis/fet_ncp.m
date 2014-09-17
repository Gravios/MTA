function [fncp,fs,ts] = fet_ncp(Trial,varargin)
%function [fncp,fs,ts] = fet_ncp(Trial,varargin)
%[sampleRate,mode,chans] = DefaultArgs(varargin,{'xyz','',2});
%SampleRate is always xyz sampleRate at the moment
[sampleRate,mode,chans,overwrite] = DefaultArgs(varargin,{'xyz','',2,false});

if ischar(sampleRate), sampleRate = Trial.(sampleRate).sampleRate;end

fs = [];
ts = [];

fncp = Trial.lfp.copy;

if ~isempty(chans),
    fncp.load(Trial,chans);
end

%SampleRate is always xyz sampleRate at the moment
Trial.load('xyz');
fncp.resample(Trial.xyz);

switch mode
  case 'wcsd'
% $$$     try,load(fullfile(Trial.path.MTAPath,'fet_ncp.arm.mat'));end
% $$$     if exist('ARmodel','var')||overwrite,
% $$$         bang = WhitenSignal(fncp.data,[],true,ARmodel);
% $$$     else
% $$$         [bang,ARmodel] = WhitenSignal(fncp.data,[],true);
% $$$         save(fullfile(Trial.path.MTAPath,'fet_rhm.arm.mat'),'ARmodel');
% $$$     end
    fncp = WhitenSignal(fncp.data,[],1);

    [ys,fs,ts] = mtcsdglong(fncp,2^9,Trial.xyz.sampleRate,2^7,2^7*.875,[],'linear',[],[1,20]);
    ts = ts+(2^6)/Trial.xyz.sampleRate;
    ssr = 1/diff(ts(1:2));
    pad = round([ts(1),mod(Trial.xyz.size(1)-2^6,2^7)/Trial.xyz.sampleRate].*ssr)-[1,0];
    szy = size(ys);
    fncp = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
    ts = cat(1,zeros([pad(1),1]),ts,zeros([pad(2),1]));

  case 'csd'
    
    [ys,fs,ts] = mtcsdglong(fncp.data,2^9,Trial.xyz.sampleRate,2^7,2^7*.875,[],'linear',[],[1,20]);
    ts = ts+(2^6)/Trial.xyz.sampleRate;
    ssr = 1/diff(ts(1:2));
    pad = round([ts(1),mod(Trial.xyz.size(1)-2^6,2^7)/Trial.xyz.sampleRate].*ssr)-[1,0];
    szy = size(ys);
    fncp = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
    ts = cat(1,zeros([pad(1),1]),ts,zeros([pad(2),1]));

% $$$   case 'spectral'
% $$$     [ys,fs,ts] = mtcsdglong(nang,2^8,Trial.xyz.sampleRate,2^7,2^7-1,[],'linear',[],[1,30]);
% $$$     fncp = zeros(Trial.xyz.size(1),size(ys,2));
% $$$     fncp((2^6+1):(size(ys)+2^6),:) = ys;
% $$$     fncp = MTADlfp('data',fncp,'sampleRate',Trial.xyz.sampleRate);
% $$$   case 'wspectral'
% $$$     load(fullfile(Trial.path.MTAPath,'fet_rhm.arm'),'-mat');
% $$$     [nang] = WhitenSignal(fncp.data,[],[],ARmodel);
% $$$     [ys,fs,ts] = mtcsdglong(nang,2^8,fncp.sampleRate,2^7,2^7-1,[],'linear',[],[1,30]);
% $$$     ts = ts+2^6/Trial.xyz.sampleRate;
% $$$     fncp = zeros(Trial.xyz.size(1),size(ys,2));
% $$$     fncp((2^6+1):(size(ys)+2^6),:) = ys;
% $$$     fncp = MTADlfp('data',fncp,'sampleRate',Trial.xyz.sampleRate);
% $$$     ts = cat(1,ts(1)-1/fncp.sampleRate*flipud(cumsum(padding(:,1)+1)),ts,ts(end)+1/fncp.sampleRate*cumsum(padding(:,1)+1));
% $$$     
% $$$   case 'Swspectral'
% $$$     try,load(fullfile(Trial.path.MTAPath,'fet_ncp.arm.mat'));end
% $$$     if exist('ARmodel','var')||overwrite,
% $$$         bang = WhitenSignal(bang,[],true,ARmodel);
% $$$     else
% $$$         [bang,ARmodel] = WhitenSignal(bang,[],true);
% $$$         save(fullfile(Trial.path.MTAPath,'fet_rhm.arm.mat'),'ARmodel');
% $$$     end
% $$$     fncp = WhitenSignal(fncp.data);
% $$$ 
% $$$     [ys,fs,ts] = mtcsdglong(fncp,2^9,Trial.xyz.sampleRate,2^7,2^7*.875,[],'linear',[],[1,20]);
% $$$     ts = ts+(2^6)/Trial.xyz.sampleRate;
% $$$     ssr = 1/diff(ts(1:2));
% $$$     pad = round([ts(1),mod(Trial.xyz.size(1)-2^6,2^7)/Trial.xyz.sampleRate].*ssr)-[1,0];
% $$$     szy = size(ys);
% $$$     fncp = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
% $$$     ts = cat(1,zeros([pad(1),1]),ts,zeros([pad(2),1]));

    
  case 'default'
    fncp = MTADlfp('data',fncp.data,'sampleRate',Trial.xyz.sampleRate);
  otherwise
    fncp = fncp.data;
end




