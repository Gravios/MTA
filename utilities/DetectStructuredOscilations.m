function DetectStructuredOscilations(Trial,varargin)

[lfp,FreqRange,MinCycles,Thr,ModelTemplate] = DefaultArgs(varagin,{Trial.lfp.copy,[160,220],5,5,'default'});

%%Diagnostic Vars
Trial = MTATrial('jg05-20120310');
lfp = Trial.lfp.copy;
lfp.load(Trial,[65:72]);
FreqRange = [160,220];
MinCycles = 5;
Thr = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%




xyz = Trial.load('xyz');
fx = ButFilter(lfp.data,2,FreqRange/lfp.sampleRate*2,'bandpass');

MinDuration = xyz.sampleRate/mean(FreqRange)*MinCycles;

WinLen = round(lfp.sampleRate/FreqRange(1)*5); 
Window = gausswin(WinLen, 5);
Window = Window/sum(Window);
amp = sqrt(filtfilt(Window,1,fx.^2));
amp = MTADlfp('data',amp,'sampleRate',lfp.sampleRate);
amp.resample(xyz);


switch ModelTemplate
  case 'default'
    ModelTemplate = load(fullfile(Trial.path.MTAPath,'DetectStructedOscilations.default.mdl.mat'));

  case 'ripples'
     
    % Make interative gui for making models
    %ripper = [81196,81208];
    %msave(fullfile(Trial.path.MTAPath,'DetectStructedOscilations.default.mdl'),urpexp-bsh)
    
    %rpexp = amp(ripper(1):ripper(2),:);
    rpexp = load(fullfile(Trial.path.MTAPath,'DetectStructedOscilations.default.mdl'));

    for i= 1:lfp.size(2),rpmean(i) = mean(amp(amp(:,i)>2,i));end
    for i= 1:lfp.size(2),rpstd(i) = std(amp(amp(:,i)>2,i));end
    urpow = (amp.data-rpmean(1))./rpstd(1);
    urpexp = rpexp;
    %urpexp = urpexp/sum(urpexp(:));

    [~,mchan] = max(urpexp,[],2)
    mchan = round(mean(mchan));
    
    turpexp = urpexp;
    fsh = .0002;
    bshift = [];
    for i = 1:80,
        turpexp = turpexp-fsh;
        crp = conv2(urpow,turpexp,'same');
        bshift(end+1) = skewness(zscore(crp(:,mchan).*urpow(:,mchan)));    
    end
    %figure,plot([1:80]*fsh,bshift,'.')
    [~,bsh] = max(bshift);
    bsh = bsh*fsh;
    crp = conv2(urpow,urpexp-bsh,'same');
    amp.data = crp(:,mchan).*urpow(:,mchan);
  case 'STDMASK'

  otherwise % Load custom model
    
end

%figure,plot(unity(ButFilter(crp(:,4),3,[.1,30]/(amp.sampleRate/2),'bandpass').*max(urpow,[],2)),'b'),
%hold on,plot(max(urpow,[],2),'r');


zamp = amp.copy;
zamp.data = zscore(amp.data);


if any(Thr>20)
    Thresh = prctile(amp.data,Thr);
else
    Thresh = Thr*std(amp.data)+mean(amp.data);
end
PowerPeaks = LocalMinima(-amp.data, ceil(1.5*MinDuration), -Thresh(1));







