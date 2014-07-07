function [rhm,fs,ts] = fet_rhm(Trial,varargin)
% [rhm,fs,ts] = fet_rhm(Trial,varargin)
% [sampleRate,mode,windowSize] = DefaultArgs(varargin,{Trial.xyz.sampleRate,'spectral',1});
% Need to update the spectral window size to adapt to the xyz.sampleRate

[sampleRate,mode,windowSize] = DefaultArgs(varargin,{Trial.xyz.sampleRate,'',1});

% Search for computationally optimal window size 
% given windowSize in seconds Default is 1 second
% swins = 2.^[1:12];
% [~,swi]=min(abs(swins-sampleRate)); 
% windowSize = swins(swi);
%Trial = MTATrial('Ed05-20140528');
%Trial = MTATrial('jg05-20120317');

fwin = gtwin(1.25,Trial.xyz.sampleRate);
swin = gtwin(.1,Trial.xyz.sampleRate);

xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gausswin(5)./sum(gausswin(5)));

rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
rbb = Trial.xyz.model.rb({'spine_lower','pelvis_root'});

hcom = xyz.com(rb);
xyz.addMarker('hcom',[.7,0,.7],{{'head_back','head_front',[0,0,1]}},hcom);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},permute(Filter0(fwin,hcom),[1,3,2]));


ang = Trial.ang.copy;
ang.create(Trial,xyz);
bang = ButFilter(ang(:,'head_back','fhcom',3),3,[1,30]./(Trial.ang.sampleRate/2),'bandpass');

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
    [ys,fs,ts] = mtcsdlong(bang,2^8,Trial.ang.sampleRate,2^7,2^7-1,[],'linear',[],[1,30]);
    for sigs = 1:size(bang,2),
        ys(:,:,sigs,sigs) = mtcsglong(bang(:,sigs),2^8,Trial.ang.sampleRate,2^7,2^7-1,[],'linear',[],[1,30]);
    end
    %[ys,fs,ts] = mtcsdglong(bang,2^8,Trial.ang.sampleRate,2^7,2^7-1,[],'linear',[],[1,30]);
    rhm = zeros(xyz.size(1),size(ys,2));
    rhm((2^6+1):(size(ys)+2^6),:) = ys;
    rhm = MTADlfp('data',rhm,'sampleRate',xyz.sampleRate);
  case 'default'
    rhm = MTADlfp('data',bang,'sampleRate',xyz.sampleRate);
  otherwise
    rhm = bang;
end






