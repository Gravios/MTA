function bhv_body_head(Trial)


%% Load Session and State Collection
sname = 'jg05-20120311';
tname = 'all';
marker = 'spine_lower';
stc_mode = 'auto_wbhr';

Trial = MTATrial(sname,tname);
Trial.stc.updateMode(stc_mode);
Trial.stc.load;



xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.1,xyz.sampleRate));
vel = xyz.vel([1:5,7],[1,2]);

wvel=WhitenSignal(vel.data,[],1);

[yv,fv,tv,phiv,fstv] = mtchglong(wvel,2^8,vel.sampleRate,2^7,2^7*.875,3,'linear',[],[1,30]);
tv = tv+(2^6)/xyz.sampleRate;
ssr = 1/diff(tv(1:2));
pad = round([tv(1),mod(xyz.size(1)-2^6,2^7)/xyz.sampleRate].*ssr)-[1,0];
szy = size(yv);
yv = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),yv,zeros([pad(2),szy(2:end)])),...
             'sampleRate',ssr);
tv = cat(1,zeros([pad(1),1]),tv,zeros([pad(2),1]));


myv = yv.copy;


vl = [];
for i= 1:6,
vl(:,i) = yv(:,2,i,i);
end
hist2([mean(vl,2),var(mean(vl,2),100,100);

hist2([1./(clip(mean(vl(nniz(vl),:),2),0,6).*clip(var(vl(nniz(vl),:),[],2),0,15)),vl(nniz(vl),1)],100,100);

figure,sp=[];
sp(end+1)=subplot(411);imagesc(yv(:,:,1,4)'),axis xy,caxis([0.5,1]),Lines(Trial.stc{'w',myv.sampleRate}(:)-2^6/diff(tv(1:2)),[],'w');
sp(end+1)=subplot(412);imagesc(yv(:,:,1,6)'),axis xy,caxis([0.5,1])
sp(end+1)=subplot(413);imagesc(sum(yv(:,:,1,1:6),4)'),axis xy,caxis([0.5,6]),Lines(Trial.stc{'w',myv.sampleRate}(:)-2^6/diff(tv(1:2)),[],'w');
sp(end+1)=subplot(414);imagesc(yv(:,:,1,4)'.*sum(yv(:,:,1,1:6),4)'),axis xy
linkaxes(sp,'xy');

figure,hist(unity(log10(myv(nniz(myv),2,1,1))),1000)



figure,hist2([unity(log10(myv(nniz(myv),2,1,1))),unity(log10(myv(nniz(myv),2,6,6)))],100,100),caxis([0,75])

umyv = myv.copy;
umyv.data(nniz(myv),:,:,:) = unity(log10(myv(nniz(myv),:,:,:)));

figure,plot(unity(log10(myv(nniz(myv),1,1,1))),unity(log10(myv(nniz(myv),1,7,7))),'.')
figure,plot(umyv(nniz(myv),1,1,1),umyv(nniz(myv),1,7,7),'.')
hold on
sind = Trial.stc{'w'};
plot(umyv(sind,1,1,1),umyv(sind,1,7,7),'.g')
sind = Trial.stc{'r'};
plot(umyv(sind,1,1,1),umyv(sind,1,7,7),'.r')
sind = Trial.stc{'h'};
plot(umyv(sind,1,1,1),umyv(sind,1,7,7),'.m')


ang = Trial.ang.copy;
ang.create(Trial,xyz);

rhm = fet_rhm(Trial);

bang = [rhm,ButFilter(ang(:,2,4,3),3,[1,30]./(Trial.ang.sampleRate/2),'bandpass'),ButFilter(ang(:,1,3,3),3,[1,30]./(Trial.ang.sampleRate/2),'bandpass')];

[bang] = WhitenSignal(bang,[],1);


[ys,fs,ts] = mtcsdglong(bang,2^8,Trial.ang.sampleRate,2^7,2^7*.875,[],'linear',[],[1,30]);

figure,
sp(1)=subplot(311);
imagescnan({ts,fs,log10(ys(:,:,1,1)')},[-8,-2.5],0,1);axis xy
sp(2)=subplot(312);
imagescnan({ts,fs,angle(ys(:,:,1,3)')},[],1,1);axis xy
sp(3)=subplot(313);
imagescnan({ts,fs,log10(ys(:,:,3,3)')},[-8,-4],0,1);axis xy
linkaxes(sp,'xy')


[yc,fc,tc,phic,fstc] = mtchglong(bang,2^8,Trial.ang.sampleRate,2^7,2^7*.875,[],'linear',[],[1,30]);

figure,
sp(1)=subplot(411);
imagescnan({tc,fc,log10(yc(:,:,1,1)')},[-8,-2.5],0,1);axis xy
sp(2)=subplot(412);
imagescnan({tc,fc,yc(:,:,1,3)'},[0.4,1],0,1);axis xy
sp(3)=subplot(413);
imagescnan({tc,fc,fstc(:,:,3)'},[0,10],0,1);axis xy
sp(4)=subplot(414);
imagescnan({tc,fc,log10(yc(:,:,3,3)')},[-8,-4],0,1);axis xy
linkaxes(sp,'xy')
][