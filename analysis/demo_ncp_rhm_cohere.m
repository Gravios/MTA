sname = 'Ed03-20140625';
tname = 'all';
mode = {'height','hangle'};
marker = 'spine_lower';
stc_mode = 'auto_wbhr';
Trial = MTATrial(sname,tname);

Trial.stc.updateMode(stc_mode);
Trial.stc.load;


xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.2,xyz.sampleRate));

[rhm,fsr,tsr] = fet_rhm(Trial,[],'Swspectral');
[ncp,fsn,tsn] = fet_ncp(Trial,[],'Swspectral',66);

[rhm_fet] = fet_rhm(Trial,[],'default');
[ncp_fet] = fet_ncp(Trial,[],'default');

wfet = WhitenSignal([rhm_fet.data,ncp_fet.data]);
[ys,fs,ts,phi,fstat] = mtchglong(wfet,2^9,ang.sampleRate,2^7,2^7*.875,[],'linear',[],[1,20]);
szy = size(ys);
padding = zeros([round(2^6/ang.sampleRate/diff(ts(1:2))),szy(2:end)]);
ys = MTADlfp('data',cat(1,padding,ys,padding),'sampleRate',1/diff(ts(1:2)));
ts = ts+2^6/ang.sampleRate;
ts = cat(1,ts(1)-(1/ys.sampleRate)*flipud(cumsum(padding(:,1)+1)),ts,ts(end)+(1/ys.sampleRate)*cumsum(padding(:,1)+1));

ftime =(1:rhm_fet.size(1))/rhm_fet.sampleRate;
ang = Trial.ang.copy;
ang.create(Trial,xyz);

ny = 4;
figure,s = 1;
sp(s)=subplot(ny,1,s);s=s+1;
hold on
%plot(ftime,(ncp_fet.data-nanmean(ncp_fet.data))./nanstd(ncp_fet.data),'g');
plot(ftime,ButFilter((ncp_fet.data-nanmean(ncp_fet.data))./nanstd(ncp_fet.data),3,[5,15]./(Trial.ang.sampleRate/2),'bandpass'),'g');
plot(ftime,4*ButFilter(unity(rhm_fet.data),3,[5,15]./(Trial.ang.sampleRate/2),'bandpass'),'r');
plot(ftime,ang(:,5,7,2).*5,'m');
Lines([],-5*.8,'k');
sp(s)=subplot(ny,1,s);s=s+1;
imagesc(ts,fsr,log10(rhm.data')),axis xy,caxis([-7,-2])
%imagesc(ts,fsr,log10(ys(:,:,1,1))'),axis xy,caxis([-7,-2])
sp(s)=subplot(ny,1,s);s=s+1;
imagesc(ts,fs,ys(:,:,1,2)');axis xy,caxis([.4,.8])
sp(s)=subplot(ny,1,s);s=s+1;
imagesc(ts,fsn,log10(ncp.data')),axis xy,caxis([1,6])
%imagesc(ts,fsn,log10(ys(:,:,2,2))'),axis xy,caxis([1,6])
linkaxes(sp,'x');


