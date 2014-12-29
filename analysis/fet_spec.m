function [rhm,fs,ts] = fet_spec(Trial,fet,varargin)
[sampleRate,mode,windowSize,overwrite] = DefaultArgs(varargin, ...
                                                  {Trial.xyz.sampleRate,'',1,false});


fs = []; ts = [];
xyz = Trial.load('xyz');


% if xyz sampling rate is greater than 120 Hz then resample it to 120 Hz
if fet.sampleRate > 120, 
    fet.resample(120); 
end
fet.filter(gausswin(5)./sum(gausswin(5)));


switch mode

  case 'csd'
    [ys,fs,ts] = mtcsdglong(fet.data,2^9,fet.sampleRate,2^7,2^7*.875,[],'linear',[],[1,20]);
    ts = ts+(2^6)/fet.sampleRate;
    ssr = 1/diff(ts(1:2));
    pad = round([ts(1),mod(fet.size(1)-2^6,2^7)/fet.sampleRate].*ssr)-[1,0];
    szy = size(ys);
    rhm = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
    ts = cat(1,zeros([pad(1),1]),ts,zeros([pad(2),1]));

  case 'wcsd'
    try,load(fullfile(Trial.path.MTAPath,[mfilename,'.arm.mat']));end
    if exist('ARmodel','var')||overwrite,
        fet.data(nniz(fet.data),:) = WhitenSignal(fet.data(nniz(fet.data),:),[],true,ARmodel);
    else
        [fet.data(nniz(fet.data),:),ARmodel] = WhitenSignal(fet.data(nniz(fet.data),:),[],true);
        save(fullfile(Trial.path.MTAPath,[mfilename,'.arm.mat']),'ARmodel');
    end
    [ys,fs,ts] = mtcsdglong(fet.data,2^9,fet.sampleRate,2^7,2^7*.875,[],'linear',[],[1,20]);
    ts = ts+(2^6)/fet.sampleRate;
    ssr = 1/diff(ts(1:2));
    pad = round([ts(1),mod(fet.size(1)-2^6,2^7)/fet.sampleRate].*ssr)-[1,0];
    szy = size(ys);
    rhm = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
    ts = cat(1,zeros([pad(1),1]),ts,zeros([pad(2),1]));

  case 'default'
    rhm = MTADlfp('data',fet.data,'sampleRate',fet.sampleRate);
  otherwise
    rhm = fet.data;
end






