
Trial = MTATrial('Ed03-20140624');

[ys,fs,ts] = fet_rhm(Trial,[],'mtchglong');
[ys,fs,ts] = fet_ncp(Trial,[],'mtchglong');


figure,
imagesc(ts,fs,log10(ys.data)');
axis('xy');
caxis([-7,-3]);

rhm = fet_rhm(Trial,[],'mta');


rhmpow = ys.copy;
rhmpow.data = log10(median(ys(:,6<fs&fs<13),2));
~nniz(rhmpow))

[State, hmm, decode] = gausshmm(rhmpow.data,3);

Stateall = zeros(Trial.xyz.size(1),1);
Stateall = State;


figure
sp(1) = subplot(211);
imagesc(ts,fs,log10(ys.data)');
axis('xy');
caxis([-7,-3]);
sp(2) = subplot(212);
plot(ts,rhmpow.data);
Lines([],-4,'k');
% $$$ plot(Stateall*10+1,'c')
% $$$ Lines(Trial.stc{'r',ys.sampleRate}.data(:),[],'r');
% $$$ Lines(Trial.stc{'w',ys.sampleRate}.data(:),[],'k');
linkaxes(sp,'x');




