function Trial = label_thetarc(Trial,varargin)
%function Trial = label_thetarc(Trial,varargin)
%[Stc,channels,overwrite] = DefaultArgs(varargin,{[],1,false});
[Stc,channels,overwrite] = DefaultArgs(varargin,{[],1,false});

% LOAD metadata and lfp
par = LoadPar(fullfile(Trial.spath, [Trial.name '.xml'])); 
lfp = LoadBinary(fullfile(Trial.spath, [Trial.name '.lfp']),channels,par.nChannels)';

% RC theta spectrum
[gys,gfs,gts] = mtchglong(diff(lfp,1,2),2^11,par.lfpSampleRate,2^10,2^10*0.5,3,[],[],[1,25]);
gts = gts+((2^10*0.5)/2)/par.lfpSampleRate;
gssr = 1/diff(gts(1:2));
gpad = round([gts(1),size(lfp,1)./par.lfpSampleRate-gts(end)].*gssr)-[1,0];
szg = size(gys);
gys = cat(1,zeros([gpad(1),szg(2:end)]),gys,zeros([gpad(2),szg(2:end)]));
gts = cat(2,([1:gpad(1)]./gssr),gts',([gpad(1)+size(gts,1)]+[1:gpad(2)])./gssr)';

% THETA DELTA ratio kinda
tdRatio = mean(log(gys(:, 6<fs & fs<9)),2)-mean(log(gys(:,3>fs|(fs>12&fs<15))),2);

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
tper = gts(thetaPeriods).*par.lfpSampleRate;
% SAVE sts file
msave(fullfile(Trial.spath,[Trial.name,'.sts.thetarc']),tper);

% $$$ figure
% $$$ subplot(211);
% $$$ imagesc(gts,gfs,log(gys)');
% $$$ axis('xy');
% $$$ colormap('parula');
% $$$ caxis([15,18]);
% $$$ Lines(gts(thetaPeriods(:,1)),[],'g');
% $$$ Lines(gts(thetaPeriods(:,2)),[],'m');
% $$$ subplot(212);
% $$$ imagesc(ts,fs,log(ys(:,:,1))');
% $$$ axis('xy');
% $$$ colormap('parula');
% $$$ %caxis([15,18])
% $$$ caxis([18,21])
% $$$ Lines(tper(:,1),[],'g');
% $$$ Lines(tper(:,2),[],'m');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'xy');

if isempty(Stc),
    Stc = Trial.stc.copy;
end

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

