function Trial = label_thetarc(Trial,varargin)
%function Trial = label_thetarc(Trial,varargin)
%[Stc,channels,overwrite] = DefaultArgs(varargin,{[],1,false});
[Stc,channels,overwrite] = DefaultArgs(varargin,{[],1,false});

% LOAD metadata and lfp
par = LoadPar(fullfile(Trial.spath, [Trial.name '.xml'])); 
lfp = LoadBinary(                               ...
    fullfile(Trial.spath, [Trial.name '.lfp']), ...
    channels,                                   ...
    par.nChannels)';


% RC theta spectrum
[ys,fs,ts] = mtchglong(diff(lfp,1,2), ...
                          2^11, ...
                          par.lfpSampleRate,...
                          2^10, ...
                          2^9, ...
                          3, ...
                          [],[],[1,35]);

ts = ts+((2^9)/2)/par.lfpSampleRate;
ssr = 1/diff(ts(1:2));
gpad = round([ts(1),size(lfp,1)./par.lfpSampleRate-ts(end)].*ssr)-[1,0];
szg = size(ys);
gpad(gpad<0) = 0;

ys = cat(1,                           ...
         zeros([gpad(1),szg(2)]), ...
         ys,                         ...
         zeros([gpad(2),szg(2)]));

ts = cat(1,                                    ...
         ([1:gpad(1)]./ssr)',                  ...
         ts,                                   ...
         ([gpad(1)+size(ts,1)]+[1:gpad(2)])'./ssr);

%nys = ys;
nys = bsxfun( ...
    @rdivide, ...
    ys(:,:)', ...
    sum(ys(:,:),2)')';

%nys = RectFilter(RectFilter(RectFilter(nys(:,:),3,1)',3,1)',3,1);

% THETA DELTA ratio 
fqin  =  (5 < fs & fs < 8);
fqout =  fs < 2;
tdRatio =   log10(mean(nys(:,fqin ),2)) ...
          -log10( mean(nys(:,fqout),2));

% GENERATE theta periods via GHMM
hmmStates = zeros(size(tdRatio));
nind = ~isnan(tdRatio) & ~isinf(tdRatio) & tdRatio~=0;
[hmmStates(nind),thhmm,thdec] = gausshmm(tdRatio(nind),2,1,0);
for ii = 1:2
    tdRatioStates(ii) = mean(tdRatio(find(hmmStates==ii)));
end

[~,tdStateInd] = max(tdRatioStates);
thetaState = zeros([numel(ts),1]);
thetaState = hmmStates==tdStateInd;
thetaPeriods = ThreshCross(thetaState,0.5,1);

% RESAMPLE to lfp sample rate
tper = ts(thetaPeriods).*par.lfpSampleRate;

% SAVE sts file
msave(fullfile(Trial.spath,[Trial.name,'.sts.thetarc']),tper);

Stc.states(Stc.gsi('c')) = [];
data = load(fullfile(Trial.spath,[Trial.name '.sts.thetarc']));
sync = Trial.lfp.sync.copy;
lsync = sync.sync.copy;
lsync.resample(Trial.lfp.sampleRate);
data = IntersectRanges(lsync.data,data)-lsync.data(1)+1;
Stc.addState(Trial.spath,...
             Trial.filebase,...
             data,...
             Trial.lfp.sampleRate,...
             Trial.sync.copy,...
             Trial.sync.data(1),...
             'thetarc','c');

Stc{'c'}.save(1);

Stc.save(1);
Trial.stc = Stc;
Trial.save;


