
wp = Trial.stc{'w'}.copy;
wp.resample(1);


units = [69,72,75,78,81,84,87];

Trial.lfp.load(Trial,units);
theta = ButFilter(Trial.lfp.data,3,[6,12]./(Trial.lfp.sampleRate/2),'bandpass');

g = hilbert(theta);
ph = phase(g);

gamma =  ButFilter(Trial.lfp.data,3,[20,210]./(Trial.lfp.sampleRate/2),'bandpass');

[wgamma,arm]= WhitenSignal(gamma,[],1);

gs = [];
for i= 1:size(wgamma,2),
    [gs(:,:,i),gfs,gts] = mtchglong(wgamma(1:400000,i),2^8,Trial.lfp.sampleRate,2^7,2^7-1,[],[],[],[25,200]);
end

sp = [];
for i= 1:size(wgamma,2),
    sp(i) = subplot(size(wgamma,2),1,i);
    imagesc(t,f,log10(gs(:,:,i)')),axis xy
end
linkaxes(sp,'xy');
Lines(wp.data(:),[],'k');
ForAllSubplots('caxis([2.7,4])')


fgs = Filter0(gausswin(7)./sum(gausswin(7)),Filter0(gausswin(151)./sum(gausswin(151)),gs(:,:,2))');
lfgs = log10(fgs);


gb = LocalMinima2(-lfgs,-2.4,5);
figure,hold on,imagesc(lfgs),axis xy,caxis([2,3]),scatter(gb(:,2),gb(:,1),6);

% FastICA
[icasig] = fastica(lfgs);
figure,imagesc(gts,gfs,icasig),axis xy

zlfgs = unity(lfgs')';
gb = LocalMinima2(-zlfgs,-1,5);

gb = LocalMinima2(-zlfgs,-.5,[5,20]);
figure,hold on,imagesc(gts,gfs,zlfgs),axis xy,
caxis([.5,3]),scatter(gts(gb(:,2)),gfs(gb(:,1)),6);

tzlfgs = zlfgs;
tzlfgs(zlfgs<.5)=0;
[ticasig] = fastica(tzlfgs);

