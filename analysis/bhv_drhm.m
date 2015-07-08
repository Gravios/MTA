
Trial = MTATrial('Ed10-20140814');
Trial = MTATrial('Ed10-20140814');

xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter;

rhm = fet_rhm(Trial,'mode','default');
%rhm.filter(gtwin(.,rhm.sampleRate));
ncp = ButFilter(fet_ncp(Trial,'chans',66),3,[1,20]/(rhm.sampleRate/2),'bandpass');
%ncp = fet_ncp(Trial,'chans',66);

%drhm = [0;diff(diff(rhm.data));0];
drhm = [0;diff(diff(ButFilter(rhm.data,3,[3,20]/(rhm.sampleRate/2),'bandpass')));0];
drhm(~nniz(drhm))=0;
frhm = ButFilter([0;diff(ButFilter(diff(ButFilter(rhm.data,3,[3,20]/(rhm.sampleRate/2),'bandpass')),3,[3,20]/(rhm.sampleRate/2),'bandpass'));0],3,[3,20]/(rhm.sampleRate/2),'bandpass');
frhm(~nniz(frhm))=0;


ang = Trial.ang.copy;
ang.load(Trial);
dang = [0;diff(diff(ButFilter(ang(:,5,7,2),3,[3,20]/(ang.sampleRate/2),'bandpass')));0];
bang = [0;diff(diff(ButFilter(ang(:,4,7,3),3,[3,20]/(ang.sampleRate/2),'bandpass')));0];
%vel = xyz.vel(7,3);
%vel = [diff(vel.data);0];
figure,plot(drhm),hold on,plot(unity(ncp)/4,'r'),plot(dang*50,'m'),plot(bang,'g')


[yv,fv,tv,phiv,fstv] = mtchglong([bang,dang*50,unity(ncp)./4],2^8,rhm.sampleRate,2^7,2^7*.875,3,'linear',[],[1,30]);
tv = tv+(2^6)/rhm.sampleRate;
ssr = 1/diff(tv(1:2));
pad = round([tv(1),mod(rhm.size(1)-2^6,2^7)/rhm.sampleRate].*ssr)-[1,0];
szy = size(yv);
yv = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),yv,zeros([pad(2),szy(2:end)])),...
             'sampleRate',ssr);
tv = cat(1,zeros([pad(1),1]),tv,zeros([pad(2),1]));

ny = 5;
figure,sp=[];
sp(end+1)=subplot(ny,1,numel(sp)+1);imagesc(tv,fv,log10(yv(:,:,1,1))'),axis xy,caxis([-5,-2]),
%set(sp(end),'outerposition',get(sp(end),'outerposition').*[1,1,1,1.05])
sp(end+1)=subplot(ny,1,numel(sp)+1);imagesc(tv,fv,yv(:,:,1,3)'),axis xy,caxis([.7,1]),
%set(sp(end),'outerposition',get(sp(end),'outerposition').*[1,1,1,1.05])
sp(end+1)=subplot(ny,1,numel(sp)+1);imagesc(tv,fv,log10(yv(:,:,2,2)')),axis xy,caxis([-5,-2]),
%set(sp(end),'outerposition',get(sp(end),'outerposition').*[1,1,1,1.05])
sp(end+1)=subplot(ny,1,numel(sp)+1);imagesc(tv,fv,yv(:,:,2,3)'),axis xy,caxis([.7,1]),
%set(sp(end),'outerposition',get(sp(end),'outerposition').*[1,1,1,1.05])
sp(end+1)=subplot(ny,1,numel(sp)+1);imagesc(tv,fv,log10(yv(:,:,3,3)')),axis xy,caxis([-7,-2]),
%set(sp(end),'outerposition',get(sp(end),'outerposition').*[1,1,1,1.05])
%sp(end+1)=subplot(ny,1,numel(sp)+1);plot(drhm),hold on,plot(unity(ncp)/4,'r')
linkaxes(sp,'xy');
