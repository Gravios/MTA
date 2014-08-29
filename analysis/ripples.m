function ripples(Session)

Session = MTASession('jg05-20120310',...
                   {{'ripples',64:71},...
                     'ufr','CluRes',...
                    {'lfp',64:71},'nq'});

%Session = Session.load_ripples(64:71,1);
Session = Session.load_lfp(64:71);

rlfp = ButFilter(Session.lfp,11,[150,250]./(Session.lfpSampleRate/2),'bandpass');
%rlfp = ButFilter(Session.lfp,11,[150,250]./(Session.lfpSampleRate/2),'bandpass');

wrlfp = WhitenSignal(rlfp,[],1);

frange=[50,350];
nffts = 2^7;
winlen = 2^5;
y=[];t=[];f=[];

[y,f,t] = mtchglong(wrlfp(1:100000,:),...
                    nffts,...
                    Session.lfpSampleRate,...
                    winlen,...
                    winlen*0.875,[],[],[],frange);

figure
for i = 1:8
sp(i) = subplot(810+i);
imagesc(t,f,log10(y(:,:,i,i))'),axis xy
Lines(Session.ripples(1:150,1)./Session.lfpSampleRate,[],'k',[],2);
end
linkaxes(sp,'xy');

wlfp = WhitenSignal(Session.lfp,[],1);

frange = [100,300];
nffts = 2^6;
windows = 2^5;
y = [];
t = [];
f = [];

rip_lfp = GetSegs(Session.lfp,Session.ripples(:,1),250,0);
rip_wlfp = GetSegs(wlfp,Session.ripples(:,1),250,0);
rip_rlfp = GetSegs(rlfp,Session.ripples(:,1),250,0);
rip_flfp = GetSegs(flfp,Session.ripples(:,1),250,0);

tl = 0:1/Session.lfpSampleRate:249/Session.lfpSampleRate;

swr = cat(1,Session.nq.SpkWidthR);
punits = find(swr>.56);
iunits = find(swr<.56);

pufr = sum(Session.ufr(:,punits),2);
iufr = sum(Session.ufr(:,iunits),2);

pufr = sum((Session.ufr(:,punits)>0),2)/size(Session.map,1);
iufr = sum((Session.ufr(:,iunits)>0),2)/size(Session.map,1);

rip_pufr = GetSegs(pufr,Session.ripples(:,1),250,0);
rip_iufr = GetSegs(iufr,Session.ripples(:,1),250,0);

figure(301),
ii = find(Session.ripples(:,2)>1300);
ri = 1;
set(gcf,'CurrentCharacter','l');
while 1
    i = ii(ri);
    clf
    [ys,fs,ts,phis,FStatss] = mtchglong(sq(rip_wlfp(:,i,:)),...
                                   nffts,...
                                   Session.lfpSampleRate,...
                                   windows,...
                                   windows*0.875,[],[],[],frange);
    sp6=subplot(611),imagesc(ts+diff(ts(1:2))/2,fs,log10(ys(:,:,5,5))'),axis xy
    sp7=subplot(612),plot(ts+diff(ts(1:2))/2,median(ys(:,:,5,5),2));
    sp1=subplot(613),plot(tl,sq(rip_flfp(:,i,:))+repmat([0:1000:7000],size(rip_flfp,1),1));
    ylim(sp1,[-1000,8000]);
    sp2=subplot(614),plot(tl,sq(rip_rlfp(:,i,:))+repmat([0:1000:7000],size(rip_lfp,1),1));
    ylim(sp2,[-1000,8000]);
    sp4=subplot(615),plot(ts+diff(t(1:2))/2,Filter0(gausswin(7)/sum(gausswin(7)),1./sum(sq(circ_var(phis(:,:,1,:),[],[],2)),2)));
    ylim(sp4,[0,500]);
    sp5=subplot(616),plot(tl,Filter0(gausswin(7)/sum(gausswin(7)),rip_pufr(:,i).*4),tl,rip_iufr(:,i));
    ylim(sp5,[0,.8]);
    linkaxes([sp1,sp2,sp4,sp5,sp6,sp7],'x')
    waitforbuttonpress
    whatkey = get(gcf,'CurrentCharacter');
    switch double(whatkey)
      case double('i')
        ri = input('Enter ri #: ');
      case double('n')
        ri = ri+1;
      case double('p')
        ri=ri-1;
      case double('q')
        return
    end
end


[y,f,t,phi,FStats] = mtchglong(wlfp,...
                               nffts,...
                               Session.lfpSampleRate,...
                               windows,...
                               windows*0.875,[],[],[],frange);
t = t+diff(t(1:2))/2;




pufrt = d2t(pufr,Session.lfpSampleRate,0)';
[~,utyt,~] = NearestNeighbour(pufrt,t);

pufrm = Filter0(gausswin(7)/sum(gausswin(7)),(pufr.*mean(iufr(iufr>0.01))/mean(pufr(pufr>0.01))+std(iufr(iufr>0.01))/std(pufr(pufr>0.01))));
iufrm = Filter0(gausswin(7)/sum(gausswin(7)),iufr);

%prt = pufr(utyt);
%irt = iufr(utyt);


prt = pufrm(utyt);
irt = iufrm(utyt);

%rrt = 

fsphi = Filter0(gausswin(7)/sum(gausswin(7)),1./sum(sq(circ_var(phi(:,:,1,:),[],[],2)),2));



ind = 1;
for shift = 1:4:80,
subplotfit(ind,length(1:4:80))
fsphis = diff(fsphi(1+shift:(end)));
prts = diff(prt(1:(end-shift)));
[out,xbin,ybin,pos] = hist2([prts,fsphis],16,100);
imagesc(xbin,ybin,log10(out)'),axis xy,colorbar
%plot(prts(fsphis>50&prts>0.05),log10(fsphis(fsphis>50&prts>0.05)),'.');
ind  = ind+1;
end


ind = 1;
for shift = 1:4:80,
subplotfit(ind,length(1:4:80))
prts = prt(1:(end-shift));
fsphis = fsphi(1+shift:(end));
ymed = sq(median(y(shift+1:end,:,4,4),2));
[out,xbin,ybin,pos] = hist2([prts(fsphis>50&prts>0.05),log10(ymed(fsphis>50&prts>0.05))],20,100);
imagesc(xbin,ybin,log10(out)'),axis xy,colorbar
ind  = ind+1;
end

ind = 1;
for shift = 1:3:30,
subplotfit(ind,length(1:5:50))
fsphis = fsphi(1:end-shift);
ymed = sq(median(y(shift+1:end,:,4,4),2));
[out,xbin,ybin,pos] = hist2([log10(fsphis(fsphis>50)),log10(ymed(fsphis>50))],100,100);
imagesc(xbin,ybin,log10(out)'),axis xy,colorbar
ind  = ind+1;
end





rfp = LocalMinima(-fsphi,10,-50)';
rpr = LocalMinima(-prt,10,-10^0.002)';
rpr = LocalMinima(-prt,10,-0.002)';

rpr = rpr(2:end-1);
rfp = rfp(2:end-1);



ymed = [];
for i=1:8,
ymed(:,i) = sq(median(y(:,:,i,i),2));
end
intymed = sum(ymed,2);


% $$$ figure,hist(intymed(LocalMinima(-fsphi,20,-50)),2000)
% $$$ figure,hist(intymed(LocalMinima(-prt,20,-0.001)),2000)
% $$$ 
% $$$ binSize = 8;
% $$$ halfBins = 80;
% $$$ normalization = 'count';
% $$$ 
% $$$ [prccg,prtbin,prpairs] = Trains2CCG({rpr,rfp},{1,2},binSize,halfBins,1/diff(t(1:2)),normalization);
% $$$ 
% $$$ figure,bar(prtbin,prccg(:,1,2))



miy = mean(log10(intymed));
siy =  std(log10(intymed));
dte = 2;
iyt = 10^(miy+dte*siy);

ind = 1;
figure
for shift = 1:2:40,
subplotfit(ind,length(1:2:40))
frfet = prt(rfp+shift)+prt(rfp+shift)./irt(rfp+shift);
%frfet(~isfinite(frfet)) = max(frfet(isfinite(frfet)));
frfet(frfet==inf)=0;
frfet = clip(frfet,0,5);
frfet(isnan(frfet)) = min(frfet(~isnan(frfet)));
phifet = log10(fsphi(rfp));
%[out,xbin,ybin,pos] = hist2([frfet(intymed(rfp-5)>iyt),phifet(intymed(rfp-5)>iyt)],50,50);
%imagesc(xbin,ybin,out'),axis xy,colorbar
plot(frfet(intymed(rfp-5)>iyt),10.^phifet(intymed(rfp-5)>iyt),'.')
ind  = ind+1;
end

figure,plot(frfet(intymed(rfp-5)>iyt),phifet(intymed(rfp-5)>iyt),'.')


shift=13
frfet = prt(rfp+shift)+prt(rfp+shift)./irt(rfp+shift);
%frfet(~isfinite(frfet)) = max(frfet(isfinite(frfet)));
frfet(frfet==inf)=0;
frfet = clip(frfet,0,5);
frfet(isnan(frfet)) = min(frfet(~isnan(frfet)));
pripind = rfp(find(frfet>1.3&intymed(rfp-1)>iyt));


%Trains2CCG(frfet>1.3,intymed(rfp-1)>iyt

figure(401),
ri = 1;
set(gcf,'CurrentCharacter','l');
while 1
i = ri
    clf
    rip_wlfp = GetSegs(wlfp,round(t(pripind(i))*Session.lfpSampleRate)-125,250,0);
    rip_flfp = GetSegs(flfp,round(t(pripind(i))*Session.lfpSampleRate)-125,250,0);
    rip_rlfp = GetSegs(rlfp,round(t(pripind(i))*Session.lfpSampleRate)-125,250,0);
    [ys,fs,ts,phis,FStatss] = mtchglong(sq(rip_wlfp),...
                                   nffts,...
                                   Session.lfpSampleRate,...
                                   windows,...
                                   windows*0.875,[],[],[],frange);


   rip_pufr = GetSegs(pufr,round(t(pripind(i))*Session.lfpSampleRate)-125,250,0);
   rip_iufr = GetSegs(iufr,round(t(pripind(i))*Session.lfpSampleRate)-125,250,0);


    sp6=subplot(611),imagesc(ts+diff(ts(1:2))/2,fs,log10(ys(:,:,5,5))'),axis xy
    sp7=subplot(612),plot(ts+diff(ts(1:2))/2,median(ys(:,:,5,5),2));
    sp1=subplot(613),plot(tl,sq(rip_flfp)+repmat([0:1000:7000],size(rip_flfp,1),1));
    ylim(sp1,[-1000,8000]);
    sp2=subplot(614),plot(tl,sq(rip_rlfp)+repmat([0:1000:7000],size(rip_lfp,1),1));
    ylim(sp2,[-1000,8000]);
    sp4=subplot(615),plot(ts+diff(t(1:2))/2,Filter0(gausswin(7)/sum(gausswin(7)),1./sum(sq(circ_var(phis(:,:,1,:),[],[],2)),2)));
    ylim(sp4,[0,500]);
    sp5=subplot(616),plot(tl,Filter0(gausswin(7)/sum(gausswin(7)),rip_pufr.*4),tl,rip_iufr);
    ylim(sp5,[0,.8]);
    linkaxes([sp1,sp2,sp4,sp5,sp6,sp7],'x')
    waitforbuttonpress
    whatkey = get(gcf,'CurrentCharacter');
    switch double(whatkey)
      case double('i')
        ri = input('Enter ri #: ');
      case double('n')
        ri = ri+1;
      case double('p')
        ri=ri-1;
      case double('q')
        return
    end
end


% $$$ %% Detect Ripples
% $$$ if ~exist([Session.spath.nlx Session.name '.spw'],'file'),
% $$$     DetectRipples([Session.spath.nlx Session.name],59:62,[],3);
% $$$ end
% $$$ Rips = load([Session.spath.nlx Session.name '.spw'],'file');
% $$$ [~,rpind] = SelectPeriods(Rips(:,1),[Session.syncPeriods(1),Session.syncPeriods(end)],'d',1,1);
% $$$ Rips = Rips(rpind);
% $$$ 

% $$$ binSize = 2;%ms
% $$$ halfBins = 32;
% $$$ normalization = 'count';
% $$$ prip = [rt-64,rt+64];
% $$$ [ripRes ripind] = SelectPeriods(Res,prip,'d',1,0);
% $$$ ripClu = Clu(ripind);
% $$$ [rip_ccg,rip_tbin,rip_pairs] = Trains2CCG({ripRes},{ripClu},binSize,halfBins,Session.lfpSampleRate,normalization);
% $$$ %[rip_ccg,rip_tbin,pairs] = CCG(ripRes,ripClu,binSize,halfBins,Session.lfpSampleRate,[],normalization,[]);


%% 20140812
Trial = MTATrial('jg05-20120310');
Trial.lfp.load(Trial,[57:72]);


[yl,fl,tl,phil,fstl] = mtchglong(WhitenSignal(Trial.lfp(1:1000000,9:16),[],1),2^12,Trial.lfp.sampleRate,2^10,2^10-8,5,'linear',[],[1,20]);
tl = tl+2^9/1250;

[ym,fm,tm,phim,fstm] = mtchglong(WhitenSignal(Trial.lfp(1:1000000,9:16),[],1),2^10,Trial.lfp.sampleRate,2^8,2^8*.875,5,'linear',[],[60,240]);
tm = tm+2^7/1250;


[yh,fh,th,phih,fsth] = mtchglong(WhitenSignal(Trial.lfp(1:1000000,9:16),[],1),2^7,Trial.lfp.sampleRate,2^6,2^6*.875,5,'linear',[],[150,250]);
th = th+2^5/1250;

figure,imagesc(tl+2^9/1250,fl,log10(yl(:,:,5,5))'),axis xy
figure,imagesc(th+2^3/1250,fh,log10(yh(:,:,5,5))'),axis xy

tpow = zeros([size(yl,1),8]);
for i=1:8,
tpow(:,i) = median(log10(yl(:,fl>5&fl<12,i,i)),2);
end

rpow = zeros([size(yh,1),8]);
for i=1:8,
rpow(:,i) = median(log10(yh(:,:,i,i)),2);
end

rpow = zeros([size(ym,1),8]);
for i=1:8,
rpow(:,i) = median(log10(ym(:,fm<200&fm>160,i,i)),2);
end


figure,imagesc(th,1:8,rpow'),axis xy,caxis([0,5])
tsf = [1:1000000]'/Trial.lfp.sampleRate;
for i=1:8
hold on,plot(tsf,bsxfun(@plus,unity(ButFilter(Trial.lfp(1:1000000,i+8),3,[4,13]/(Trial.lfp.sampleRate/2),'bandpass'))/4,i))
hold on,plot(tsf,bsxfun(@plus,unity(ButFilter(Trial.lfp(1:1000000,i+8),3,[150,250]/(Trial.lfp.sampleRate/2),'bandpass'))/4,i))
end





figure
tsf = [1:1000000]'/Trial.lfp.sampleRate;
hold on,plot(repmat(tsf,[1,8]),bsxfun(@plus,unity(ButFilter(Trial.lfp(1:1000000,9:16),3,[4,13]/(Trial.lfp.sampleRate/2),'bandpass'))/4,[1:8]*2))
hold on,plot(bsxfun(@plus,unity(ButFilter(Trial.lfp(1:1000000,9:16),3,[150,250]/(Trial.lfp.sampleRate/2),'bandpass'))/4,[1:8]*2))


rpexp = rpow(th>216.10&th<216.15,:);
%rpexp = rpow(tm>216.10&tm<216.15,:);

for i= 1:8,rpmean(i) = mean(rpow(rpow(:,i)>2,i));end
for i= 1:8,rpstd(i) = std(rpow(rpow(:,i)>2,i));end

%urpow = bsxfun(@ldivide,bsxfun(@minus,rpow,rpmean(1)),rpstd(1));
urpow = (rpow-rpmean(1))./rpstd(1);
%urpow  = rpow;
%urpexp = bsxfun(@ldivide,bsxfun(@minus,rpexp,rpmean),rpstd);
urpexp = rpexp;
urpexp = urpexp/sum(urpexp(:));
crp = conv2(urpow,urpexp,'same');

urpexpds = urpexp([3:6,2,1]);
crpd = conv2(urpow./sum(urpexp(:)),urpexpds,'same');
urpexpus = urpexp([8,7,1:6]);
crpu = conv2(urpow./sum(urpexp(:)),urpexpus,'same');
crp = max([mean(crp(:,3:4),2),mean(crpd(:,3:4),2),mean(crpu(:,3:4),2)],[],2);


figure,imagesc(th,1:8,crp'),caxis([0,4])

figure,plot(ButFilter(crp(:,1),3,[.1,30]/((1/diff(th(1:2)))/2),'bandpass'))
figure,hist(ButFilter(crp(1:10:end,1),3,[.1,30]/((1/diff(th(1:2)))/2),'bandpass'),1000)
figure,hist(log10(ButFilter(crp(:,4),3,[.1,30]/((1/diff(th(1:2)))/2),'bandpass')),1000)

figure,imagesc(th(th>12&th<750),1:16,unity([rpow(th>12&th<750,:),tpow(tl>12&tl<750,:)])'),axis xy,caxis([-1,4])
hold on,plot(th(th>12&th<750),crp(th>12&th<750,4)*7-15,'m')

figure,plot(ButFilter(crp(:,4),3,[.1,30]/((1/diff(th(1:2)))/2),'bandpass').*max(urpow,[],2))

ripfet = crp(:,4).*max(urpow,[],2);
ripfet = ButFilter(ripfet,3,[5]/((1/diff(th(1:2)))/2),'high');
srfet = GetSegs(ripfet,LocalMinima(-ripfet)-3,5);
srfet = sum(srfet)';
% $$$ srfet(srfet<0)=.001;
% $$$ srfet = log10(srfet);
figure,hist(srfet(1:end),100)

figure,plot(th,ButFilter(crp(:,1),3,[.1,30]/((1/diff(th(1:2)))/2),'bandpass'))