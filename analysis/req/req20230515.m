

MjgER2016_load_data();



Trial = Trials{20};

xyz = preproc_xyz(Trial,'trb');
vxy = vel(filter(copy(xyz),'ButFilter',4,2,'low'),'hcom',[1,2]);
lvxy = vel(filter(copy(xyz),'ButFilter',4,0.1,'low'),'hcom',[1,2]);


lfp = Trial.load('lfp',Trial.meta.channelGroup.theta);
phz = load_theta_phase(Trial,lfp.sampleRate);

stc = Trial.stc.copy();
int = select_units(Trial,'int');
spk = Trial.load('spk',lfp.sampleRate,'',int);


thetaTroughs = abs(phz.data-pi)<0.1&[0;diff(phz.data-pi)]>0;

thetaTroughs = LocalMinima(-convn(thetaTroughs,ones([1,21]),'same'),0,21);

rlfp = Trial.load('lfp',Trial.meta.channelGroup.thetarc);
rlfp.data = diff(rlfp.data,1,2);


flfp = filter(lfp.copy(),'ButFilter',4,[5,12],'bandpass');

frlfp = filter(copy(rlfp),'ButFilter',4,[5,12],'bandpass');
rmins = LocalMinima(frlfp.data-16000,10,-2e4);

glfp = Trial.load('lfp',[69,73]);
glfp.data = diff(glfp.data,1,2);

mlfp = Trial.load('lfp',[41,48]);
mlfp.data = diff(mlfp.data,1,2);

zlfp = Trial.load('lfp',[74,77]);
zlfp.filter('ButFilter',4,[1,150],'bandpass');
zlfp.data = diff(zlfp.data,1,2);
zlfp.filter('ButFilter',4,1,'high');


dlfp = resample(Trial.load('lfp',[65:2:96]),250);
dlfp.data = diff(dlfp.data,1,2);
ddlfp = copy(dlfp);
ddlfp.data = diff(dlfp.data(:,[11,12]),1,2);

figure,
hold('on');
plot([1:size(lfp,1)]./lfp.sampleRate,lfp.data, 'b', 'LineWidth',2);
plot([1:size(lfp,1)]./lfp.sampleRate,rlfp.data,'r', 'LineWidth',2);
plot([1:size(lfp,1)]./lfp.sampleRate,dlfp.data,'m', 'LineWidth',2);


figure,
hold('on');
plot([1:size(lfp,1)]./lfp.sampleRate,lfp.data,'b', 'LineWidth',2);
plot([1:size(lfp,1)]./lfp.sampleRate,rlfp.data,'r', 'LineWidth',2);
%plot([1:size(lfp,1)]./lfp.sampleRate,frlfp.data,'g', 'LineWidth',2);
plot([1:size(lfp,1)-1]./lfp.sampleRate,diff(frlfp.data)*12,'m', 'LineWidth',2);
plot([1:size(lfp,1)-1]./lfp.sampleRate,diff(flfp.data)*12,'k',  'LineWidth',2);
plot([1:size(lfp,1)-1]./lfp.sampleRate,diff(flfp.data)*12-diff(frlfp.data)*12,'c', 'LineWidth',2);
%plot([1:size(lfp,1)]./lfp.sampleRate,zlfp.data,'k', 'LineWidth',2);

drlfp = lfp.copy();
drlfp.data = [frlfp.data];
drphz = drlfp.phase([5,13]);
drthetaTroughs = abs(drphz.data-pi)<0.1&[0;diff(drphz.data-pi)]>0;
drthetaTroughs = LocalMinima(-convn(drthetaTroughs,ones([1,21]),'same'),0,21);


defspec = struct('nFFT',2^8,'Fs',zlfp.sampleRate,...
                            'WinLength',2^7,'nOverlap',2^7*.875,...
                            'FreqRange',[4,150]);
[rhm,fs,ts] = fet_spec(Trial,zlfp,'mtchglong',true,[],defspec);

figure,
hold('on');
imagesc(ts,fs,log10(rhm.data)');
axis('xy');
colormap('jet');
plot([1:size(lfp,1)]./lfp.sampleRate,zlfp.data./2000+30,'k', 'LineWidth',2);

figure,
hold('on');
plot(ts,mean(log10(rhm(:,fs<120&fs>80)),2))
plot(ts,mean(log10(rhm(:,fs<60&fs>40)),2))
plot([1:size(lfp,1)]./lfp.sampleRate,zlfp.data./2000,'k', 'LineWidth',2);

rhm.data = log10(rhm.data);
rhm.resample(dlfp);




figure,
hold('on');
for c = 1:size(dlfp,2),
plot([1:size(dlfp,1)]./dlfp.sampleRate,dlfp(:,c)-c*5000,'b', 'LineWidth',2);
end
plot([1:size(lfp,1)]./lfp.sampleRate,rlfp.data,'m', 'LineWidth',2);
plot([1:size(dlfp,1)]./dlfp.sampleRate,diff(dlfp.data(:,[11,12]),1,2)+4000,'r', 'LineWidth',2);

figure();
hold('on');
plot([1:size(dlfp,1)]./dlfp.sampleRate,dlfp.data(:,[11]),'b', 'LineWidth',2);
plot([1:size(dlfp,1)]./dlfp.sampleRate,dlfp.data(:,[12]),'b', 'LineWidth',2);
plot([1:size(dlfp,1)]./dlfp.sampleRate,diff(dlfp.data(:,[11,12]),1,2)+16000,'r', 'LineWidth',2);
plot([1:size(lfp,1)]./lfp.sampleRate,rlfp.data-6000,'m', 'LineWidth',2);
plot([1:size(lfp,1)]./lfp.sampleRate,frlfp.data-16000,'m', 'LineWidth',2);
Lines(rmins./lfp.sampleRate,[],'k');


fdlfp = filter(copy(ddlfp),'ButFilter',4,2,'low');

tpowCT = sum(sqrt(segs(filter(copy(lfp),'ButFilter',4,[5,11],'bandpass'),thetaTroughs-256,512).^2))';
tpowRC = sum(sqrt(segs(filter(copy(rlfp),'ButFilter',4,[5,11],'bandpass'),thetaTroughs-256,512).^2))';


dpow = fdlfp(round(thetaTroughs/1250*250));

tvel = log10(vxy(round(thetaTroughs/lfp.sampleRate.*xyz.sampleRate)));

sper = [stc{'x&t-s-m',lfp.sampleRate}];
sper = [stc{'x&t',lfp.sampleRate}];
sper = [stc{'x+p&t',lfp.sampleRate}];
sper = [stc{'t&s',lfp.sampleRate}];
% $$$ sper = [stc{'s-t',lfp.sampleRate}];
figure,
ind = WithinRanges(thetaTroughs,sper.data) & nniz(tvel) & nniz(phzOffset) & nniz(sqrt(phzOffsetVar))<30;
subplot(211);
plot(tvel(ind),...
     dpow(ind),...
     '.');
subplot(212);
plot(tvel(ind),...
     intOffset(ind),...
     '.'); 


tfreq = 1./(-(thetaTroughs-circshift(thetaTroughs,-1))./lfp.sampleRate);
figure
plot(tfreq(ind),tvel(ind),'.');

figure();
sper = [stc{'s&t',lfp.sampleRate}];
ind = WithinRanges(thetaTroughs,sper.data) & nniz(tvel) & nniz(phzOffset) & nniz(sqrt(phzOffsetVar))<50;
hold('on');
plot(phzOffset(ind),tvel(ind),'.');
mean(phzOffset(ind))

figure
sper = [stc{'x&t',lfp.sampleRate}];
ind = WithinRanges(thetaTroughs,sper.data) & nniz(tvel) & nniz(phzOffset) & nniz(sqrt(phzOffsetVar))<10;
plot(phzOffset(ind),tvel(ind),'.r');
mean(phzOffset(ind))


figure
sper = [stc{'x+p&t',lfp.sampleRate}];
ind = WithinRanges(thetaTroughs,sper.data) & nniz(tvel) & nniz(phzOffset) & nniz(sqrt(phzOffsetVar))<50;
plot(phzOffset(ind),log10(tpowRC(ind)),'.r');
mean(phzOffset(ind))

figure();
hold('on');
sper = [stc{'x+p&t',lfp.sampleRate}];
ind = WithinRanges(thetaTroughs,sper.data) & nniz(tvel) & nniz(phzOffset) & nniz(sqrt(phzOffsetVar))<20;
plot(log10(tpowRC(ind))-0.2,log10(tpowCT(ind)),'.b');
sper = [stc{'s&t',lfp.sampleRate}];
ind = WithinRanges(thetaTroughs,sper.data) & nniz(tvel) & nniz(phzOffset) & nniz(sqrt(phzOffsetVar))<20;
plot(log10(tpowRC(ind))-0.2,log10(tpowCT(ind)),'.r');
line([5,6.4],[5,6.4]);


%%% THIS
figure,
subplot(211);
sper = [stc{'x+p&t',lfp.sampleRate}];
ind = WithinRanges(thetaTroughs,sper.data) & nniz(tvel) & nniz(phzOffset) & nniz(sqrt(phzOffsetVar))<20;
hist2([log10(tpowRC(ind))-0.2,log10(tpowCT(ind))],linspace([5.2,6.2,20]),linspace([5.2,6.2,20]),'yprob');
line([5.2,6.2],[5.35,6.1],'Color','m');
%line([5.45,5.95],[5.2,6.2],'Color','m');
caxis([0,0.3])
subplot(212);
sper = [stc{'s&t',lfp.sampleRate}];
ind = WithinRanges(thetaTroughs,sper.data) & nniz(tvel) & nniz(phzOffset) & nniz(sqrt(phzOffsetVar))<20;
hist2([log10(tpowRC(ind))-0.2,log10(tpowCT(ind))],linspace([5.2,6.2,20]),linspace([5.2,6.2,20]),'yprob');
line([5.2,6.2],[5.35,6.1],'Color','m');
%line([5.45,5.95],[5.2,6.2],'Color','m');
caxis([0,0.3])
colormap('jet');



figure
%plot(phzOffset(ind),tvel(ind),'.');
hist2([phzOffset(ind),tvel(ind)],linspace(-10,25,20),linspace(0.5,1.8,20),'xprob')
%hist2([phzOffset(ind),log10(tpowRC(ind))./log10(tpowCT(ind))],linspace(-10,25,20),linspace(0.98,1.1,20),'xprob')
colormap('jet')
caxis([0,0.1])


figure
%plot(phzOffset(ind),tvel(ind),'.');
hist2([intOffset(ind),tvel(ind)],linspace(-20,10,20),linspace(0.5,1.8,20),'yprob')
colormap('jet')
caxis([0,0.1])

figure
plot(tfreq(ind),phzOffset(ind),'.');



[rho,pval] = corr(tvel(ind), phzOffset(ind))


[rho,pval] = corr(tpowRC(ind), phzOffset(ind))
[rho,pval] = corr(tpowCT(ind), phzOffset(ind))

[rho,pval] = corr(tpowCT(ind), tpowRC(ind))

[rho,pval] = corr(tvel(ind), intOffset(ind))

[rho,pval] = corr(tvel(ind), tfreq(ind))




figure,
plot(intOffset,...
     phzOffset,...
     '.');


nfft = 128;
FreqRange = [0,1000/diff(bins(1:2))];
fo = ([1:nfft/2+1]-1)'*lfp.sampleRate/nfft;
select = find( fo > FreqRange(1,1) & fo < FreqRange(1,2));
fo = fo(select);

figure,plot(lfp.data),Lines(ttr,[],'k');    
Lines(ttr(1:2:end),[],'r');    
    ttr = thetaTroughs(9975:10975);

u =6;
fftOutAll = nan([255,1000]);
mspk1 = spk(units(u));
mspk2 = spk(units(11));
phzOffset = nan([numel(thetaTroughs),1]);
intOffset = nan([numel(thetaTroughs),1]);
phzOffsetVar = nan([numel(thetaTroughs),1]);
intOffsetVar = nan([numel(thetaTroughs),1]);
for tt = 3:numel(thetaTroughs)-3,
    ttr = thetaTroughs(tt-2:tt+2);
    tspk1 = mspk1(WithinRanges(mspk1,ttr([1,end])'));    
    tspk2 = mspk2(WithinRanges(mspk2,ttr([1,end])'));        
    [mccg, bins] = CCG([tspk1;tspk2;ttr],...
                       [ones(size(tspk1));2*ones(size(tspk2));3.*ones(size(ttr))],...
                       1,...
                       64,...
                       spk.sampleRate,...
                       [1,2,3],...
                   'count');
% $$$     figure,bar(bins,mccg(:,1,3));
phzOffset(tt) = sum(bins'.*mccg(:,1,3)./sum(mccg(:,1,3)));
intOffset(tt) = sum(bins'.*mccg(:,1,2)./sum(mccg(:,1,2)));

phzOffsetVar(tt) = sum(((bins'-phzOffset(tt)).*mccg(:,1,2)./(sum(mccg(:,1,2))-1)).^2);
intOffsetVar(tt) = sum(((bins'-intOffset(tt)).*mccg(:,1,2)./(sum(mccg(:,1,2))-1)).^2);

% $$$ fftOut = fft(mccg(:,1,2)',64).*sqrt(2);
% $$$ fftOutAll(:,tt) = fftOut(select).*conj(fftOut(select));
end

sper = [stc{'x&t-s-m',lfp.sampleRate}];
sper = [stc{'x&t',lfp.sampleRate}];
sper = [stc{'x+p&t',lfp.sampleRate}];
sper = [stc{'t&s',lfp.sampleRate}];

u = 15
int(u)
mspk = spk(int(u));
mspk = mspk(WithinRanges(mspk,sper.data));
mtroughs= thetaTroughs(WithinRanges(thetaTroughs,sper.data));
figure
subplot(121);
[mccg, bins] = CCG([mspk;mtroughs],...
                       [ones(size(mspk));2.*ones(size(mtroughs))],...
                       4,...
                       128,...
                       spk.sampleRate,...
                       [1,2],...
                       'hz');
bar(bins,mccg(:,1,2));
Lines([],15,'k');
Lines([],2,'k');
subplot(122);
mtroughs= drthetaTroughs(WithinRanges(drthetaTroughs,sper.data));    
[mccg, bins] = CCG([mspk;mtroughs],...
                       [ones(size(mspk));2.*ones(size(mtroughs))],...
                       4,...
                       128,...
                       spk.sampleRate,...
                       [1,2],...
                       'hz');
bar(bins,mccg(:,1,2));
Lines([],15,'k');
Lines([],2,'k');
linkaxes(findobj(gcf,'Type','Axes'),'y');





figure,
hold('on');
plot(lfp.data/3000);
plot(thetaTroughs,phzOffset);
plot(thetaTroughs,sqrt(phzOffsetVar));



figure,
hold('on');
plot(lfp.data/1000);
plot(thetaTroughs,intOffset);
plot(thetaTroughs,sqrt(intOffsetVar));

figure,plot(fo,fftOutAll(:,tt))


figure,imagesc(1:1000,fo,log10(fftOutAll));
axis('xy');
caxis([-4,-1.7])
colormap jet


figure,imagesc(1:1000,fo,...
               bsxfun(@rdivide,...
                      log10(imgaussfilt(fftOutAll*10000,1.5)),...
                      log10(sum(imgaussfilt(fftOutAll*10000,1.5)))...
                      ) ...
               );
axis('xy');
colormap jet





trainNextIndex = 1;
for ind = 129:32:size(xyz,1)-129,
    twin = [-128,128]+ind;
% $$$     trainStartIndex = trainNextIndex;
% $$$     trainEndIndex = trainStartIndex;
% $$$     while true
% $$$         if mspk(trainEndIndex) < twin(2)
% $$$             trainEndIndex = trainEndIndex + 1;
% $$$         else
% $$$             trainNextIndex = trainEndIndex;
% $$$             break
% $$$         end
% $$$     end
% $$$     if (trainEndIndex - trainStartIndex) < 2 
% $$$         continue
% $$$     end
% $$$     tspk = mspk(trainStartIndex:trainEndIndex);
    tspk = mspk(WithinRanges(mspk,twin));
    tphz = thetaTroughs
    [mccg, bins] = CCG([tspk],...
                       [ones(size(tspk))],...
                       1,...
                       96,...
                       spk.sampleRate,...
                       [1],...
                   'hz');
    accg(ind,:,u) = mccg(:,1,1)';
end
end

figure,imagesc(accg(1:32:end,:)')
    

% $$$ FreqRange = [0,xyz.sampleRate/2];
% $$$ fo = ([1:256/2+1]-1)'*xyz.sampleRate/256;
% $$$ select = find( fo > FreqRange(1,1) & fo < FreqRange(1,2));
% $$$ fo = fo(select);


fftOutAll = nan([size(xyz,1),127,numel(units)]);
for u = 1:numel(units)
for i = 129:32:size(xyz,1)-129
    fftOut = fft(accg(i,:,u)',256).*sqrt(2);
    fftOutAll(i,:,u) = fftOut(select).*conj(fftOut(select));
end
end

figure,plot(imgaussfilt(log10(fftOutAll(1:32:end,1,u))))

u = 6
figure,imagesc((1:32:size(fftOutAll,1))./xyz.sampleRate,fo,...
               bsxfun(@rdivide,...
                      imgaussfilt(log10(fftOutAll(1:32:end,:,u))',1.5),...
                      sum(imgaussfilt(log10(fftOutAll(1:32:end,:,u))',1.5))...
                      ) ...
               );

figure,imagesc((1:32:size(fftOutAll,1))./xyz.sampleRate,fo,...
               bsxfun(@rdivide,...
                      imgaussfilt(log10(mean(fftOutAll(1:32:end,:,[5,6,11,13]),3,'omitnan'))',1.5),...
                      sum(imgaussfilt(log10(mean(fftOutAll(1:32:end,:,[5,6,11,13]),3,'omitnan'))',1.5))...
                      ) ...
               );
axis('xy');
colormap('jet');
title(num2str(units(u)))
% $$$ figure,imagesc((1:64:size(fftOutAll,1))./xyz.sampleRate,fo,bsxfun(@rdivide,imgaussfilt(log10(fftOutAll(1:64:end,:))',2),sum(imgaussfilt(log10(fftOutAll(1:64:end,:))',2))));
%caxis([3,8]);   
caxis([0.007,0.011])
colormap('jet');
tper = [stc{'t',1}]
Lines(tper.data(:,1),[],'g');
Lines(tper.data(:,2),[],'m');

accg15 = accg;

figure,imagesc((1:60:size(fftOutAll,1))./xyz.sampleRate,...
               fo,...
               bsxfun(@rdivide,log10(fftOutAll(1:60:end,:))'
axis('xy');
%caxis([3,8]);   
caxis([0.006,0.013])
colormap('jet');




size(CCG([mspk],...
                       [ones(size(mspk))],...
                       1,...
                       80,...
                       spk.sampleRate,...
                       [1],...
                   'hz'))



u = 13;
figure,
subplot (221) ,hist2([log10(vxy(1:32:end)),log10(fftOutAll(1:32:end,18,u))],linspace([-3,2,50]),linspace([-1,8,50]),'xprob'),caxis([0,0.1,]),colormap('jet');
subplot (222) ,hist2([log10(vxy(1:32:end)),log10(fftOutAll(1:32:end,1,u))],linspace([-3,2,50]),linspace([2,8,50]),'yprob'),caxis([0,0.1,]),colormap('jet');
subplot (223) ,hist2([log10(fftOutAll(1:32:end,1,u)),log10(fftOutAll(1:32:end,14,u))],linspace([3,8,50]),linspace([-1,8,50]),'yprob'),caxis([0,0.1,]),colormap('jet');



figure,hist2([log10(vxy(1:32:end)),log10(fftOutAll(1:32:end,18,u))],linspace([-3,2,50]),linspace([-1,8,50]))



figure,plot(log10(vxy(1:32:end)),log10(fftOutAll(1:32:end,1,u)),'.');

figure,plot(log10(vxy(1:32:end)),sum(log10(fftOutAll(1:32:end,30:58,u)),2),'.');


figure,plot(log10(fftOutAll(1:32:end,12,u)),log10(fftOutAll(1:32:end,18,u)),'.');



units = select_units(Trial,'pyr');
%dc = accumulate_decoding_vars_simple(Trial,units);
dc = accumulate_decoding_vars(Trial,units);



dcrmins = round(rmins./lfp.sampleRate.*dc.sampleRate);

refind = false([size(xyz,1),1]);
refind(dc.ind) = true;

rmnind = false([size(xyz,1),1]);
rmnind(dcrmins) = true;

ind = find( refind & rmnind );

ind = ismember(dc.ind,swr) & dc.stcm(:,1)~=1  & (dc.stcm(:,11)==11);
ind = ismember(dc.ind,dcrmins(log10(-frlfp(rmins))>3.8)) & dc.stcm(:,1)~=1  & (dc.stcm(:,11)==11);
ind = ismember(dc.ind,dcrmins) & dc.stcm(:,1)~=1  & (dc.stcm(:,5)==5);

ind = dc.stcm(:,1)==1 & dc.stcm(:,3)==3 & ~(dc.stcm(:,11)==11 |dc.stcm(:,10)==10);



swrPer = [Trial.stc{'R&s',dc.sampleRate}];
swr = mean(swrPer.data,2);

figure,plot(dc.sax(ind,1),mean(rhm(dc.ind(ind),fs>80&fs<120),2)./mean(rhm(dc.ind(ind),fs>40&fs<60),2),'.')
figure,plot(dc.sax(ind,2),mean(rhm(dc.ind(ind),fs>80&fs<120),2)./mean(rhm(dc.ind(ind),fs>40&fs<60),2),'.')

figure,plot(dc.sax(ind,2),mean(rhm(dc.ind(ind),fs>80&fs<120),2),'.')
figure,plot(dc.sax(ind,2),mean(rhm(dc.ind(ind),fs>40&fs<60),2),'.')

[mccg,bins] = CCG([swr;dcrmins], ...
                  [ones(size(swr)); 2*ones(size(dcrmins))], ...
                  1,...
                  50,...
                  dc.sampleRate,...
                  [1,2],...
                  'count');
figure,
bar(bins,mccg(:,1,2));

Lines(bins([35,65]),[],'k');



ii = 10175;
ind = ismember(dc.ind,dcrmins(ii)-30:dcrmins(ii)+30);
figure;
subplot(211);
plot(dc.sax(ind,1),dc.sax(ind,2),'-.');
xlim([-450,450]);
ylim([-450,450]);
subplot(212);
plot(dc.esax(ind,2),dc.esax(ind,1),'-.');
xlim([-650,650]);
ylim([-650,650]);
figure
bar(sum(dc.uinc(ind,:)))

units(find(sum(dc.uinc(ind,:))>0))


figure,
hist2([dc.esax(ind,2),dc.esax(ind,1)],...
      linspace(-800,800,50),...
      linspace(-800,800,50));
colormap('jet');

figure,
hist2([dc.sax(ind,1),dc.sax(ind,2)],...
      linspace(-450,450,50),...
      linspace(-450,450,50));
colormap('jet');


figure,
hist2([dc.com(ind,1),dc.com(ind,2)],...
      linspace(-450,450,50),...
      linspace(-450,450,50));
colormap('jet');




rsegs = sum(sqrt(ButFilter(GetSegs(lfp,rmins-64-300,128,0),4,[150,250]./(lfp.sampleRate/2),'bandpass').^2))';




figure,hist2([log10(-frlfp(rmins)),log10(rsegs)],50,50)


