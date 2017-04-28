%bhv_psd

%% Figure 14 - Spectral power for each behavior
%sname = 'jg05-20120317';
sname = 'jg05-20120310';
%sname = 'jg05-20120309';
disp = 0;
%chans = [1:2:8];
chans = [65:1:96];
marker = 'spine_lower';

Trial = MTATrial(sname,'all');
ang = Trial.load('ang');
xyz = Trial.load('xyz');
lfp = Trial.load('lfp',chans);
Stc = Trial.load('stc','NN0317');

wlfp = WhitenSignal([lfp.data(:,6),lfp.data],[],1);

matlabpool open 12

tl=[];
fl=[];
yl=[];
spec.nfft = 2^11;
spec.window = 2^10;
parfor i = 1:size(wlfp,2),
    [yl(:,:,i),fl(:,i),tl(:,i)] = mtchglong(wlfp(:,i),...
                                            spec.nfft,...
                                            lfp.sampleRate,...
                                            spec.window,...
                                            spec.window*0.875,...
                                            [],'linear',[],[1,40]);
end

fl = fl(:,1);
tl = tl(:,1);
yldSampleRate = 1/diff(tl(1:2,1));
tshift = round(spectral.window/2/lfp.sampleRate*yldSampleRate);
yld = MTADlfp('data',yl,'sampleRate',yldSampleRate);
yld.data = cat(1,zeros([tshift,yld.size(2:end)]),yld.data);
yld.data(yld.data==0)=nan;
yld.data = log10(yld.data);
yld.data = (yld.data-repmat(nanmean(yld.data),[yld.size(1),1,1]))./repmat(nanstd(yld.data),[yld.size(1),1,1]);


% $$$ 
% $$$ 
% $$$ 
% $$$ sts='w'
% $$$ evt=1;
% $$$ figure,cnt=1;
% $$$ for i=linspace(-2,2,40).*yld.sampleRate
% $$$ imagesc(fl,1:13,sq(nanmean(yld(round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)+i-tshift),:,:),1))'),
% $$$ caxis([-1,1])
% $$$ text( 35,13,num2str((i)./yld.sampleRate),'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7])
% $$$ pause(.2)
% $$$ frms(cnt) = getframe;
% $$$ cnt=cnt+1;
% $$$ end
% $$$ 
% $$$ prm.fps = 5;
% $$$ prm.loop = inf;
% $$$ makeGif(frms,'/gpfs01/sirota/bach/homes/gravio/figures/spect_rear_on_theta.gif',prm);
% $$$ 
% $$$ 
% $$$ th=[];
% $$$ fh=[];
% $$$ yh=[];
% $$$ spectral.nfft = 2^9;
% $$$ spectral.window = 2^7;
% $$$ spectral.freq = [40,120];
% $$$ parfor i = 1:Trial.lfp.size(2),
% $$$     [yh(:,:,i),fh(:,i),th(:,i)] = mtchglong(wlfp(:,i),spectral.nfft,Trial.lfp.sampleRate,spectral.window,spectral.window*0.875,[],[],[],spectral.freq);
% $$$ end
% $$$ fh = fh(:,1);
% $$$ th = th(:,1);
% $$$ yhd = MTADlfp('data',yh,'sampleRate',1/diff(th(1:2,1)));
% $$$ yhd.data(yhd.data==0)=nan;
% $$$ yhd.data = log10(yhd.data);
% $$$ yhd.data = (yhd.data-repmat(nanmedian(yhd.data),[yhd.size(1),1,1]))./repmat(nanstd(yhd.data),[yhd.size(1),1,1]);
% $$$ tshift = round(spectral.window/2/Trial.lfp.sampleRate*yhd.sampleRate);
% $$$ 
% $$$ sts='r'
% $$$ evt=2;
% $$$ figure,cnt=1;
% $$$ for i=linspace(-2,2,200).*yhd.sampleRate
% $$$ imagesc(fh,1:yhd.size(3),sq(nanmean(yhd(round(Trial.stc{sts,yhd.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)+i-tshift),:,:),1))'),
% $$$ caxis([-.51,.51])
% $$$ text( 80,30,num2str(i./yhd.sampleRate),'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7])
% $$$ pause(.2)
% $$$ frms(cnt) = getframe;
% $$$ cnt=cnt+1;
% $$$ end
% $$$ 
% $$$ prm.fps = 10;
% $$$ prm.loop = inf;
% $$$ makeGif(frms,'/gpfs01/sirota/bach/homes/gravio/figures/spect_rear_off_gamma.gif',prm);
% $$$ 
% $$$ 
% $$$ figure,
% $$$ s = 'r'
% $$$ sp1 = subplot(121);hold on
% $$$ sp2 = subplot(122);hold on
% $$$ colors= jet(numel(chans))';
% $$$ 
% $$$ for c = colors
% $$$ plot(sp1,fl,-.5*find(ismember(colors',c','rows'))+nanmean(log10(yld(Trial.stc{s,yld.sampleRate}.data(2:end-1,:),:,ismember(colors',c','rows')))),'color',c)
% $$$ plot(sp2,fh,-.21*find(ismember(colors',c','rows'))+nanmean(log10(yhd(Trial.stc{s,yhd.sampleRate},:,ismember(colors',c','rows')))),'color',c)
% $$$ end
% $$$ 
% $$$ 
% $$$ figure,
% $$$ %subplot(121),
% $$$ chan = 15;
% $$$ st1 = 'l'; st2 = 'g'; 
% $$$ st1 = 'r'; st2 = 'w';
% $$$ boundedline(fl,mean(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),2*std(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),'-r',fl,mean(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),2*std(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),'-b','alpha')
% $$$ 
% $$$ figure,
% $$$ %subplot(121),
% $$$ st1 = 'l'; st2 = 'g';
% $$$ boundedline(fl,mean(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),2*std(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),'-r',fl,mean(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),2*std(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),'-b','alpha')
% $$$ 
% $$$ 
% $$$ subplot(122),
% $$$ 
% $$$ boundedline(fh,mean(log10(yhd(Trial.stc{'r',yhd.sampleRate}.data(2:end-1,:),:,3))),std(log10(yhd(Trial.stc{'r',yhd.sampleRate}.data(2:end-1,:),:,3))),'-r',fh,mean(log10(yhd(Trial.stc{'w',yhd.sampleRate}.data(2:end-1,:),:,3))),std(log10(yhd(Trial.stc{'w',yhd.sampleRate}.data(2:end-1,:),:,3))),'-b','alpha')
% $$$ 
% $$$ 
% $$$ evt = 1;
% $$$ sts = 'r';
% $$$ figure,imagesc(sq(nanmean(GetSegs(sq(nanmean(yld.data(:,fl<12&fl>6,:),2)),round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yld.sampleRate),round(4*yld.sampleRate),0),2))')
% $$$ 
% $$$ evt = 1;
% $$$ sts = 'r';
% $$$ figure,imagesc(sq(nanmean(GetSegs(sq(nanmean(yld.data(:,fl<25&fl>18,:),2)),round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yld.sampleRate),round(4*yld.sampleRate),0),2))')
% $$$ 
% $$$ 
% $$$ evt = 1;
% $$$ sts = 'r';
% $$$ figure,imagesc(sq(nanmean(GetSegs(sq(nanmean(yhd.data(:,fh<120&fh>80,:),2)),round(Trial.stc{sts,yhd.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yhd.sampleRate),round(4*yhd.sampleRate),0),2))')
% $$$ 
% $$$ tdr = sq(mean(yl(:,fl>6&fl<12,:),2)./(mean(yl(:,fl<5,:),2)+mean(yl(:,fl<18&fl>13,:),2))/2);
% $$$ tdr = MTADlfp('data',tdr,'sampleRate',1/diff(tl(1:2,1)));
% $$$ 
% $$$ 
% $$$ tdr.data = (tdr.data-repmat(nanmedian(tdr.data),[tdr.size(1),1,1]))./repmat(nanstd(tdr.data),[tdr.size(1),1,1]);
% $$$ 
% $$$ 
% $$$ evt = 1;
% $$$ sts = 'r';
% $$$ figure,imagesc(sq(nanmean(GetSegs(tdr.data,round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yld.sampleRate-tshift/yld.sampleRate),round(4*yld.sampleRate),0),2))')
% $$$ 
% $$$ sts = 'w';
% $$$ [U,V,D] = svd(cov(yld(round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-tshift/yld.sampleRate),:,21)));
% $$$ figure,imagesc(fl,fl,U),axis xy

s = 'rwgl';

figure,
for i = 1:numel(s)
subplot(1,numel(s),i)
imagesc(fl,1:32,sq(nanmean(yld(Trial.stc{s(i),yld.sampleRate}.data+tshift,:,:)))');
%caxis([-1,1.2])
%subplot(1,numel(s),i),imagesc(fl,1:32,sq(nanstd(yld(Trial.stc{s(i),yld.sampleRate}.data+tshift,:,:)))');
%caxis([1.6,3])
end

Trial.xyz.filter(gausswin(31)./sum(gausswin(31)));
spd = MTADxyz('data',Trial.vel,'sampleRate',Trial.xyz.sampleRate);
spd.resample(yld);

figure,
for i = 1:numel(s)
    subplot(1,numel(s),i)
    plot(log10(spd(Trial.stc{s(i)},7)),nanmean(yld(Trial.stc{s(i)},fl<12&fl>6,1),2),'.')
    ylim([-4,4])
    xlim([-2,2])
end


for i = 1:numel(s)
%[rho,pval] =corr(log10(spd(Trial.stc{s(i)},7)),sq(nanmean(yld(Trial.stc{s(i)},fl<12&fl>6,:),2)));
[tvcov] = nancov([log10(spd(Trial.stc{s(i)},7)),sq(nanmean(yld(Trial.stc{s(i)},fl<12&fl>6,:),2))]);
tvcorr(:,i) = tvcov(:,1)./(sqrt(diag(tvcov))*sqrt(tvcov(1)));
end
figure,
imagesc(tvcorr(2:end,:))


