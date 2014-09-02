function [rhm,fs,ts] = fet_rhm(Trial,varargin)
% [rhm,fs,ts] = fet_rhm(Trial,varargin)
% [sampleRate,mode,windowSize] = DefaultArgs(varargin,{Trial.xyz.sampleRate,'spectral',1});
% Need to update the spectral window size to adapt to the xyz.sampleRate

[sampleRate,mode,windowSize,overwrite] = DefaultArgs(varargin,{Trial.xyz.sampleRate,'',1,false});

% Search for computationally optimal window size 
% given windowSize in seconds Default is 1 second
% swins = 2.^[1:12];
% [~,swi]=min(abs(swins-sampleRate)); 
% windowSize = swins(swi);

fs = [];
ts = [];

xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gausswin(5)./sum(gausswin(5)));

rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(hcom,3,[2]./(Trial.ang.sampleRate/2),'low'));


if xyz.sampleRate > 120, 
    xyz.resample(120); 
end

ang = Trial.ang.copy;
ang.create(Trial,xyz);
bang = ButFilter(ang(:,'head_back','fhcom',3),3,[1,30]./(Trial.ang.sampleRate/2),'bandpass');
bang = [0;diff(ButFilter(diff(bang),3,[1,30]/(ang.sampleRate/2),'bandpass'));0];
switch mode
  case 'spectral'
    [ys,fs,ts] = mtcsdglong(bang,2^8,Trial.ang.sampleRate,2^7,2^7-1,[],'linear',[],[1,30]);
    rhm = zeros(xyz.size(1),size(ys,2));
    rhm((2^6+1):(size(ys)+2^6),:) = ys;
    rhm = MTADlfp('data',rhm,'sampleRate',xyz.sampleRate);

  case 'wspectral'
    try,load(fullfile(Trial.path.MTAPath,'fet_rhm.arm.mat'));end
    if exist('ARmodel','var'),
    bang = WhitenSignal(bang,[],[],ARmodel);
    else
        [bang,ARmodel] = WhitenSignal(bang);
        save(fullfile(Trial.path.MTAPath,'fet_rhm.arm.mat'),'ARmodel');
    end
    [ys,fs,ts] = mtcsdglong(bang,2^8,Trial.ang.sampleRate,2^7,2^7-1,[],'linear',[],[1,30]);
    ts = ts+(2^6)/xyz.sampleRate;
    ssr = 1/diff(ts(1:2));
    pad = round([ts(1),mod(xyz.size(1)-2^6,2^7)/xyz.sampleRate].*ssr)-[1,0];
    szy = size(ys);
    rhm = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
    ts = cat(1,zeros([pad(1),1]),ts,zeros([pad(2),1]));
  
  case 'Swspectral'
    try,load(fullfile(Trial.path.MTAPath,'fet_rhm.arm.mat'));end
    if exist('ARmodel','var')||overwrite,
        bang = WhitenSignal(bang,[],[],ARmodel);
    else
        [bang,ARmodel] = WhitenSignal(bang);
        save(fullfile(Trial.path.MTAPath,'fet_rhm.arm.mat'),'ARmodel');
    end
    [ys,fs,ts] = mtcsdglong(bang,2^9,ang.sampleRate,2^7,2^7*.875,[],'linear',[],[1,20]);
    ts = ts+(2^6)/xyz.sampleRate;
    ssr = 1/diff(ts(1:2));
    pad = round([ts(1),mod(xyz.size(1)-2^6,2^7)/xyz.sampleRate].*ssr)-[1,0];
    szy = size(ys);
    rhm = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
    ts = cat(1,zeros([pad(1),1]),ts,zeros([pad(2),1]));

  case 'Sspectral'
    [ys,fs,ts] = mtcsdglong(bang,2^9,ang.sampleRate,2^7,2^7*.875,[],'linear',[],[1,20]);
    ts = ts+(2^6)/xyz.sampleRate;
    ssr = 1/diff(ts(1:2));
    pad = round([ts(1),mod(xyz.size(1)-2^6,2^7)/xyz.sampleRate].*ssr)-[1,0];
    szy = size(ys);
    rhm = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
    ts = cat(1,zeros([pad(1),1]),ts,zeros([pad(2),1]));

  case 'default'
    rhm = MTADlfp('data',bang,'sampleRate',xyz.sampleRate);
  otherwise
    rhm = bang;
end






