function Stc = label_ripples(Trial,varargin)
%function label_ripples(Trial,varargin)
%
% Inputs: 
%    Trial:      (MTASession)             
%    Stc:        (MTAStateCollection)     Trial.stc.copy()
%    chans:      (numericArray)           68:69
%    freqRange:  (numericArray)           [150,250]
%    specWindow: (numeric)                50/Trial.lfp.sampleRate
%
% NOTE - Requires theta periods to have been defined
% NOTE - Expects lfp sample rate to be between 1000 and 1500 Hz

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('Stc',               Trial.stc.copy(),                                         ...
                 'chans',             [40,48,56,64,68],                                         ...
                 'freqRange',         [140,220]                                                ...
);
[Stc,chans,freqRange] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

Trial = MTATrial.validate(Trial);
Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',chans);

specArgsGamma = struct('nFFT',2^7,...
                       'Fs',  lfp.sampleRate,...
                       'WinLength',2^6,...
                       'nOverlap',2^6*0.875,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',freqRange);

[ysg,fsg,tsg] = fet_spec(Trial,lfp,[],[],[],specArgsGamma,[],true);


% COMPUTE mean inter-channel coherence
sp = [];
coherenceWeights = zeros([size(ysg,1),1]);
interChanCount = 0;
for s = 1:numel(chans-1),
    for c = s+1:numel(chans),
        coherenceWeights = coherenceWeights + ysg(:,:,s,c);
        interChanCount = interChanCount+1;
    end
end
coherenceWeights = coherenceWeights./interChanCount;

% GET normalization coefficients 
ind = get(resample(cast(Stc{'t'},'TimeSeries'),ysg),'data')&nniz(ysg);
ripMean = mean(mean(mean(log10(ysg(ind,fsg>freqRange(1)&fsg<freqRange(2),:)),2),3));
ripStd  = std(mean(mean(log10(ysg(ind,fsg>freqRange(1)&fsg<freqRange(2),:)),2),3));

% NORMALIZE ripple power
ripplePower = nunity(mean(mean(log10(ysg(:,fsg>freqRange(1)&fsg<freqRange(2),:)),2),3),...
                     @nan,...
                     ripMean,...
                     ripStd);
% WEIGHT ripple power by inter-channel coherence 
ripplePowerWeighted = ysg.copy('empty');
ripplePowerWeighted.data = ripplePower.*mean(coherenceWeights(:,fsg>freqRange(1)&fsg<freqRange(2)),2);
ripplePowerWeighted.filter('ButFilter',3,10,'low');


% FIND ripple periods
rper = ThreshCross(ripplePowerWeighted.data,0.75,5);
for p = 1:size(rper,1),
    if max(ripplePowerWeighted(rper(p,1):rper(p,2)))<1,
        rper(p,:) = nan;
    end
end
rper(~nniz(rper),:) = [];
rmerge = 5 > (circshift(rper(:,1),-1)-rper(:,2));
rper = [rper(~rmerge,1),rper(~circshift(rmerge,-1),2)];


% ADD state to collection
Stc.addState(Trial.spath,...
             Trial.filebase,...
             rper,...
             ysg.sampleRate,...
             Trial.sync.copy,...
             Trial.sync.data(1),...
             'ripple','R');

% SAVE state to trial project folder
Stc{'R'}.save(1);

% END MAIN -----------------------------------------------------------------------------------------


% $$$ figure
% $$$ hold('on');
% $$$ plot([ripplePowerWeighted.data]);
% $$$ plot((1:size(lfp,1))./lfp.sampleRate.*ripplePowerWeighted.sampleRate,nunity(lfp(:,[end-1,end]))./5-2);
% $$$ plot((1:size(lfp,1))./lfp.sampleRate.*ripplePowerWeighted.sampleRate,uic/10-5);
% $$$ plot(rper(:),ripplePowerWeighted(rper(:)),'*m');
% $$$ Lines([],0.5,'k');
