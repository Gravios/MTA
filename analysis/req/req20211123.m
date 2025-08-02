% Identification of lfp and csd components which correlate with differential modulation of hippocampal
% input pathways.
configure_default_args();

MjgER2016_load_data();

trialId = 18;
sampleRate = 250;
sampleRate = 30;
nBin = 20;

% RASTER STUFF
Trial = Trials{trialId};    % jg05-20120312.cof.all

xyz = preproc_xyz(Trial,'trb',sampleRate);
fxyz = filter(xyz.copy(),'ButFilter',3,14,'low');
vxy = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
stc = Trial.stc.copy();

lvxy = filter(copy(vxy),'ButFilter',4,2,'low');;
lvxy.data(lvxy.data<=0.0001) = 0.0001;
lvxy.data = log10(lvxy.data);
lvxyLim = [-1.5,2];
lvxyBin = linspace([lvxyLim,nBin]);
lvxyCtr = mean([lvxyBin(1:end-1);lvxyBin(2:end)]);
lvxyInd = discretize(lvxy.data,lvxyBin);

Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',[65:96]);
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
[lys,lfs,lts] = fet_spec(Trial,lfp,[],[],[],specArgsTheta);
lys.model = specArgsTheta;
lys.model.ts = lts;
lys.model.fs = lfs;
lys.ext = 'spectra';
lys.label = 'LIN32_LOW';
lys.name = 'LIN32-HPC-CA1-DG';
lys.key = 't';
lys.sync = copy(lfp.sync);
lys.orign = lfp.origin;
lys.update_filename(Trial);
lys.path = lfp.path;
lys.save();

[tys,tfs,tts] = fet_spec(Trial,elfp,[],[],[],specArgsTheta);

elfp = copy(lfp);
elfp.data = lfp(:,[7,14,19]);
specArgsTheta = struct('nFFT',2^10,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^9,...
                  'nOverlap',2^9*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
[tys,tfs,tts] = fet_spec(Trial,elfp,[],[],[],specArgsTheta);



Trial.lfp.filename = [Trial.name,'.lfp'];
mlfp = Trial.load('lfp',[57,60,64]);
mlfp.filter('ButFilter',4,30,'low');
specArgsTheta = struct('nFFT',2^10,...
                  'Fs',  mlfp.sampleRate,...
                  'WinLength',2^9,...
                  'nOverlap',2^9*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
[mys,mfs,mts] = fet_spec(Trial,mlfp,[],false,[],specArgsTheta);



% $$$ Trial.lfp.filename = [Trial.name,'.lfp'];
% $$$ clfp = Trial.load('lfp',[41,48,49,56,57,64]);
% $$$ clfp.filter('ButFilter',4,30,'low');

elfp = lfp.copy();
elfp.data = clfp(:,5)-clfp(:,6);
specArgsTheta = struct('nFFT',2^10,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^9,...
                  'nOverlap',2^9*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
[cys,cfs,cts] = fet_spec(Trial,elfp,[],false,[],specArgsTheta);



plfp = copy(clfp);
plfp.data = diff(clfp(:,[5,6]),1,2);
phzProx = plfp.phase([4,12]);
phzLim = [-pi,pi];
phzBin = linspace([phzLim,37]);
phzCtr = mean([phzBin(1:end-1);phzBin(2:end)]);
phzInd = discretize(phzProx.data,phzBin);



dc = accumulate_decoding_vars(Trial,                               ...
                              units{trialId},                      ...
                              sessionList(trialId).thetaRefGeneral,...
                              phzCorrection(trialId),              ...
                              headRotation{trialId},               ...
                              hbangCorrection{trialId});

hvfl = nan([size(xyz,1),1]);
%hvfl(dc.ind,:) = dc.hvfl(:,2); hvflLim = [-60,60];
hvfl(dc.ind,:) = dc.hvfl(:,1); hvflLim = [-20,80];
hvflBin = linspace([hvflLim,nBin]);
hvflCtr = mean([hvflBin(1:end-1);hvflBin(2:end)]);
hvflInd = discretize(hvfl(:,1),hvflBin);


% $$$ dc = accumulate_decoding_vars(Trial,                               ...
% $$$                               units{trialId},                      ...
% $$$                               sessionList(trialId).thetaRefGeneral,...
% $$$                               phzCorrection(trialId),              ...
% $$$                               headRotation{trialId},               ...
% $$$                               hbangCorrection{trialId});


states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};
stateColors = 'krggbbmy';

figure();
sax = gobjects([0,1]);
sax(end+1) = subplot(211);
hold('on');
plot([1:size(clfp,1)]./clfp.sampleRate,clfp(:,1)-clfp(:,2),'k');
plot([1:size(clfp,1)]./clfp.sampleRate,clfp(:,3)-clfp(:,4),'r');
plot([1:size(clfp,1)]./clfp.sampleRate,clfp(:,5)-clfp(:,6),'b');
plot([1:size(clfp,1)]./clfp.sampleRate,((clfp(:,3)-clfp(:,4))-(clfp(:,5)-clfp(:,6)))*2.5 -6000,'m');
Lines([],0,'k')
Lines([],-6000,'k')
sax(end+1) = subplot(212);
hold(sax(end),'on');
plotSTC(Trial.stc,1,'text',states,stateColors);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
ylim([1,9]);
xlabel('Time (s)');
linkaxes(sax,'x');


elfp = lfp.copy();
elfp.data = ((clfp(:,3)-clfp(:,4))-(clfp(:,5)-clfp(:,6)));
specArgsTheta = struct('nFFT',2^9,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^8,...
                  'nOverlap',2^8*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,80]);
[dys,dfs,dts] = fet_spec(Trial,elfp,[],[],[],specArgsTheta);



figure();
sax = gobjects([0,1]);
sax(end+1) = subplot(411);
hold('on');
imagesc(cts,cfs,log10(lys(:,:,19))'); 
axis('xy'); colormap('jet');
sax(end+1) = subplot(412);
hold('on');
imagesc(cts,cfs,log10(cys(:,:,1))'); 
axis('xy'); colormap('jet');
sax(end+1) = subplot(413);
hold('on');
imagesc(dts,dfs,log10(dys(:,:,1))'); 
axis('xy'); colormap('jet');
sax(end+1) = subplot(414);
hold(sax(end),'on');
plotSTC(Trial.stc,1,'text',states,stateColors);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
ylim([1,9]);
xlabel('Time (s)');
linkaxes(sax,'x');



tpow = copy(tys);
tpow.data = sq(mean(log10(tys(:,tfs>6&tfs<10,:)),2));
tpow.resample(xyz);
tpowLim = [3.5,6.5];
tpowBin = linspace([tpowLim,nBin]);
tpowCtr = mean([tpowBin(1:end-1);tpowBin(2:end)]);
tpowInd = discretize(tpow.data,tpowBin);

mpow = copy(mys);
mpow.data = sq(mean(log10(mys(:,mfs>6&mfs<10,:)),2));
mpow.resample(xyz);
mpowLim = [7.75,9.25];
mpowBin = linspace([mpowLim,nBin]);
mpowCtr = mean([mpowBin(1:end-1);mpowBin(2:end)]);
mpowInd = discretize(mpow.data,mpowBin);

cpow = copy(cys);
cpow.data = sq(mean(log10(cys(:,cfs>6&cfs<10,:)),2));
cpow.resample(xyz);
cpowLim = [7.75,9.25];
cpowBin = linspace([cpowLim,nBin]);
cpowCtr = mean([cpowBin(1:end-1);cpowBin(2:end)]);
cpowInd = discretize(cpow.data,cpowBin);

cSpow = copy(clfp);
cSpow.data = clfp(:,5)-clfp(:,6);
cSpow.filter('ButFilter',4,[5,12],'bandpass');
cSpow.data = sq(sqrt(sum(cSpow.segs(1:size(cSpow,1),2^8).^2)))';
cSpow.data = circshift(cSpow.data,-2^7);
cSpow.resample(xyz);
cSpow.data = log10(cSpow.data);
cSpowLim = [3.8,5];
cSpowBin = linspace([cSpowLim,nBin]);
cSpowCtr = mean([cSpowBin(1:end-1);cSpowBin(2:end)]);
cSpowInd = discretize(cSpow.data,cSpowBin);

dpow = copy(dys);
dpow.data = sq(mean(log10(dys(:,10:13,:)),2));
dpow.resample(xyz);
dpowLim = [0,3];
dpowBin = linspace([dpowLim,nBin]);
dpowCtr = mean([dpowBin(1:end-1);dpowBin(2:end)]);
dpowInd = discretize(dpow.data,dpowBin);

lpow = copy(mys);
lpow.data = sq(mean(log10(mys(:,mfs<15,2)),2));
lpow.resample(xyz);
lpowLim = [7.4,8.8];
lpowBin = linspace([lpowLim,nBin]);
lpowCtr = mean([lpowBin(1:end-1);lpowBin(2:end)]);
lpowInd = discretize(lpow.data,lpowBin);



ppow = copy(mpow);
ppow.data = mpow(:,1)./(mpow(:,2)); ppowLim = [0.95,1.07];
%ppow.data = mpow(:,3)./(mpow(:,2)); ppowLim = [0.95,1.07];
%ppow.data = mpow(:,1)./(mpow(:,3)); ppowLim = [0.95,1.07];
ppowBin = linspace([ppowLim,nBin]);
ppowCtr = mean([ppowBin(1:end-1);ppowBin(2:end)]);
ppowInd = discretize(ppow.data,ppowBin);

rpow = copy(tpow);
%rpow.data = log2(cpow(:)./(tpow(:,1))); rpowLim = [0.55,0.95];
%rpow.data = abs(log2(cSpow(:)./(tpow(:,1)))); rpowLim = [0,0.5];
rpow.data = cpow(:)./(mpow(:,2)); rpowLim = [0.975,1.1];
%rpow.data = cpow(:)./(mpow(:,1)); rpowLim = [0.95,1.1];
%rpow.data = cpow(:)./(mpow(:,3)); rpowLim = [0.95,1.1];
%rpow.data = (cSpow(:)./mpow(:,2)); rpowLim = [0.45,0.6];
%rpow.data = (cSpow(:)./mpow(:,2)); rpowLim = [0.45,0.6];
%rpow.data = mpow(:,1)./(mpow(:,2)); rpowLim = [0.95,1.07];
%rpow.data = mpow(:,3)./(mpow(:,2)); rpowLim = [0.95,1.07];
%rpow.data = mpow(:,1)./(mpow(:,3)); rpowLim = [0.95,1.07];
%rpow.data = mean([mpow(:,1)./mpow(:,3),mpow(:,1)./mpow(:,2)],2); rpowLim = [0.95,1.07];
%rpow.data = unity(tpow(:,7))-unity(cpow(:));
%rpowLim = [0.4,1.6];
%rpowLim = [0.8,1.1];
%rpowLim = [1,1.1];
%rpowLim = [0,3];
%rpowLim = [-3,3];
%rpowLim = [0.5,1.4];
rpowBin = linspace([rpowLim,nBin]);
rpowCtr = mean([rpowBin(1:end-1);rpowBin(2:end)]);
rpowInd = discretize(rpow.data,rpowBin);
 


tper = [Trial.stc{'t-s-m'}];
tper.cast('TimeSeries');
tper.resample(xyz);

pyrCSD = copy(clfp);
pyrCSD.data = ((clfp(:,3)-clfp(:,4))-(clfp(:,5)-clfp(:,6)));
%pyrCSD.data = ((clfp(:,3)-clfp(:,4))-(clfp(:,5)-clfp(:,6)));

xmin = LocalMinima(abs(circ_dist(phzProx.data,phzCtr(18))),70,0.1);
vq = interp1(xmin./pyrCSD.sampleRate,pyrCSD.data(xmin),[1:size(tpow)]./tpow.sampleRate)';
nind = logical(tper.data) & nniz(tpow(:,19)) & nniz(vq);
figure,plot(tpow(nind,14),vq(nind)','.')

ccc = nan([32,numel(phzCtr)]);
for p = 1:numel(phzCtr)
    xmin = LocalMinima(abs(circ_dist(phzProx.data,phzCtr(p))),60,pi/2);
    vq = interp1(xmin./pyrCSD.sampleRate,pyrCSD.data(xmin),[1:size(tpow)]./tpow.sampleRate)';
    nind = logical(tper.data) & nniz(tpow(:,19)) & nniz(vq);
    for c = 1:32
        ccc(c,p) = corr(tpow(nind,c)./tpow(nind,19),vq(nind));
    end
end
figure,imagesc(phzCtr,1:32,ccc)
colormap('jet');



mazeCntrDist = sqrt(sum(xyz(:,'hcom',[1,2]).^2,3));


ind(dc.ind(dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>2&dc.ucnt>2)) = true;
ind(dc.ind(dc.stcm(:,2)==2&dc.stcm(:,1)==1)) = true;

ind = false([size(xyz,1),1]);
ind(dc.ind((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1 ...
           & dc.ucnt>2 & mazeCntrDist(dc.ind)<350 dc.iphz==18 )) = true;
ind(dc.ind((dc.stcm(:,4)==4|dc.stcm(:,6)==6)&dc.stcm(:,1)==1 & dc.ucnt>2 &  mazeCntrDist(dc.ind)<350)) = true;
ind(dc.ind((dc.stcm(:,6)==6)&dc.stcm(:,1)==1 & dc.ucnt>2 & mazeCntrDist(dc.ind)<350)) = true;
ind(dc.ind((dc.stcm(:,4)==4)&dc.stcm(:,1)==1)) = true;
ind(dc.ind((dc.stcm(:,6)==6)&dc.stcm(:,1)==1)) = true;
%ind(dc.ind(dc.stcm(:,3)==3&dc.stcm(:,1)==1)) = true;
mcdlim = [0,350];

    
    

efet = 'ecom';
erd = nan([size(xyz,1),1]);    
erd(dc.ind) = sqrt(sum(dc.(efet).^2,2));
erdLim = [0,600];
erdBin = linspace([erdLim,40]);
erdInd = discretize(erd,erdBin);

erf = nan([size(xyz,1),1]);    
erf(dc.ind) = dc.(efet)(:,1);
erfLim = [-300,400];
erfBin = linspace([erfLim,40]);
erfInd = discretize(erf,erfBin);

erl = nan([size(xyz,1),1]);    
erl(dc.ind) = dc.(efet)(:,2);
erlLim = [-200,200];
erlBin = linspace([erlLim,40]);
erlInd = discretize(erl,erlBin);

ind = false([size(xyz,1),1]);
ind(dc.ind((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1 ...
           & dc.ucnt>1 & mazeCntrDist(dc.ind)<350 &dc.phz>pi )) = true;


ind = false([size(xyz,1),1]);
ind(dc.ind(dc.stcm(:,8)==8 & dc.stcm(:,1)~=1 ...
           & dc.ucnt>1 & mazeCntrDist(dc.ind)<350  )) = true;

% PAUSE ALL
ind = false([size(xyz,1),1]);
ind(dc.ind((dc.stcm(:,6)==6) ...
           & dc.stcm(:,1)==1 ... %& lvxy(dc.ind,2)<0.5 ...
           & dc.post>0.051 ...
           & dc.ucnt>2 & mazeCntrDist(dc.ind)<325  )) = true;

% PAUSE ALL
ind = false([size(xyz,1),1]);
ind(dc.ind((dc.stcm(:,4)==4) ...
           & dc.stcm(:,1)==1 ... %& lvxy(dc.ind,2)<0.5 ...
           & dc.post>0.05 ...
           & dc.ucnt>3 & mazeCntrDist(dc.ind)<350  )) = true;

% PAUSE ASC
ind = false([size(xyz,1),1]);
ind(dc.ind((dc.stcm(:,6)==6)&dc.stcm(:,1)==1 & abs(dc.hvfl(:,2))<20 ... ...
           & dc.ucnt>1 & mazeCntrDist(dc.ind)<350 &dc.phz>pi &dc.phz<1.9*pi )) = true;

% PAUSE DEC
ind = false([size(xyz,1),1]);
ind(dc.ind((dc.stcm(:,6)==6)&dc.stcm(:,1)==1 & abs(dc.hvfl(:,2))<20 ... ...
           & dc.ucnt>1 & mazeCntrDist(dc.ind)<350 & dc.phz<pi & dc.phz>0.2 )) = true;

% LOC ALL
ind = false([size(xyz,1),1]);
ind(dc.ind((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1 ...
           & dc.ucnt>1 & mazeCntrDist(dc.ind)<350 )) = true;

% LOC ASC
ind = false([size(xyz,1),1]);
ind(dc.ind((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1 ...
           & dc.ucnt>1 & mazeCntrDist(dc.ind)<350&dc.phz>pi &dc.phz<1.9*pi )) = true;

% LOC DEC
ind = false([size(xyz,1),1]);
ind(dc.ind((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1 ...
           & dc.ucnt>1 & mazeCntrDist(dc.ind)<350 & dc.phz<pi & dc.phz>0.2 )) = true;


% ALL ALL
ind = false([size(xyz,1),1]);
ind(dc.ind((dc.stcm(:,3)==3|dc.stcm(:,5)==5|dc.stcm(:,6)==6) ...
           &dc.stcm(:,1)==1  ...
           & dc.post>0.05 ...
           & dc.ucnt>3 & mazeCntrDist(dc.ind)<300)) = true;

% ALL ASC
ind = false([size(xyz,1),1]);
ind(dc.ind((dc.stcm(:,3)==3|dc.stcm(:,5)==5|dc.stcm(:,6)==6)&dc.stcm(:,1)==1 ...
           & dc.ucnt>1 & mazeCntrDist(dc.ind)<350&dc.phz>pi &dc.phz<1.9*pi )) = true;

% ALL DEC
ind = false([size(xyz,1),1]);
ind(dc.ind((dc.stcm(:,3)==3|dc.stcm(:,5)==5|dc.stcm(:,6)==6)&dc.stcm(:,1)==1 ...
           & dc.ucnt>1 & mazeCntrDist(dc.ind)<350 & dc.phz<pi & dc.phz>0.2 )) = true;


iphz = nan([size(xyz,1),1]);
iphz(dc.ind) = dc.phz;
phzBins = linspace(0,2*pi,25);

figure();
hist2([erf(ind),iphz(ind)],erfBin,phzBins,'xprob');

figure();
hist2([lpow(ind),erd(ind)],lpowBin,erdBin,'yprob');
caxis([0,0.1])
corr(cpow(ind),erd(ind))




xVal = rpow(:,1); xInd = rpowInd(:,1); xBin = rpowBin; xCtr = rpowCtr;
%yVal = cpow(:,1); yInd = cpowInd(:,1); yBin = cpowBin; yCtr = cpowCtr;
%yVal = mpow(:,3); yInd = mpowInd(:,1); yBin = mpowBin; yCtr = mpowCtr;
yVal = lpow(:,1); yInd = lpowInd(:,1); yBin = lpowBin; yCtr = lpowCtr;
%yVal = ppow(:,1); yInd = ppowInd(:,1); yBin = ppowBin; yCtr = ppowCtr;
%yVal = lvxy(:,1); yInd = lvxyInd(:,1); yBin = lvxyBin; yCtr = lvxyCtr;
%yVal = hvfl(:,1); yInd = hvflInd(:,1); yBin = hvflBin; yCtr = hvflCtr;
%yVal = hvfl(:,1); yInd = hvflInd(:,1); yBin = hvflBin; yCtr = hvflCtr;
aVal = erf;
figure,
subplot(141);
nind =  ind & nniz(xInd) & nniz(yInd);
outc = hist2([xVal(nind,1),yVal(nind)],xBin,yBin);
imagesc(xCtr,yCtr,log10(outc)');
Lines(xCtr(round(numel(xCtr)./2)),[],'k');Lines([],yCtr(round(numel(yCtr)./2)),'k');
axis('xy');
colormap('jet');
subplot(142);
nind =  ind & nniz(xInd) & nniz(yInd);
out = accumarray([xInd(nind,1),yInd(nind)],aVal(nind),[numel(xCtr),numel(yCtr)],@mean);
out(outc(:)<20) = nan;
imagesc(xCtr,yCtr,out');
Lines(xCtr(round(numel(xCtr)./2)),[],'k');Lines([],yCtr(round(numel(yCtr)./2)),'k');
axis('xy');
colormap('jet');
caxis([-60,100])
subplot(143);
out = accumarray([xInd(nind,1),yInd(nind)],log10(vxy(nind,2)),[numel(xCtr),numel(yCtr)],@mean);
out(outc(:)<20) = nan;
imagesc(xCtr,yCtr,out');
Lines(xCtr(round(numel(xCtr)./2)),[],'k');Lines([],yCtr(round(numel(yCtr)./2)),'k');
axis('xy');
colormap('jet');
subplot(144);
out = accumarray([xInd(nind,1),yInd(nind)],log10(vxy(nind,2)),[numel(xCtr),numel(yCtr)],@std);
out(outc(:)<20) = nan;
imagesc(xCtr,yCtr,out');
Lines(xCtr(round(numel(xCtr)./2)),[],'k');Lines([],yCtr(round(numel(yCtr)./2)),'k');
axis('xy');
colormap('jet');
corr(xVal(nind,1),erf(nind,1))
corr(xVal(ind,1),erf(ind,1))



figure,plot(xVal(nind),erf(nind),'.')
figure,plot(xVal(ind),erf(ind),'.')

figure,plot(log10(lpow(ind)./rpow(ind)),erf(ind),'.')
corr((lpow(ind)-mean(lpow(ind)))./(rpow(ind)-mean(rpow(ind))),erf(ind))

figure();
plot(unity(lpow(ind)),unity(rpow(ind)),'.')

rots = -pi:0.1:pi,
nfc =[];
for r = 1:numel(rots);
nf = multiprod([unity(lpow(ind)),unity(rpow(ind))],[cos(rots(r)),-sin(rots(r));sin(rots(r)),cos(rots(r))],[2],[1,2]);
% $$$ figure,
% $$$ subplot(121);
% $$$ plot(nf(:,1),nf(:,2),'.');
% $$$ subplot(122);
% $$$ plot(nf(:,2),erf(ind),'.');
nfc(r,1) = corr(nf(:,1),erf(ind));
nfc(r,2) = corr(nf(:,2),erf(ind));
end
figure,plot(rots,nfc(:,1))




rots = -pi:0.1:pi,
nfc =[];
for r = 1:numel(rots);
nf = multiprod([unity(lpow(ind,1)),unity(rpow(ind))],[cos(rots(r)),-sin(rots(r));sin(rots(r)),cos(rots(r))],[2],[1,2]);
nfc(r,1) = corr(nf(:,1),lvxy(ind,2));
nfc(r,2) = corr(nf(:,2),lvxy(ind,2));
end
figure,plot(rots,nfc(:,1))

figure,
subplot(121);
plot(nf(:,1),lvxy(ind,2),'.');
subplot(122);
plot(nf(:,2),lvxy(ind,2),'.');

figure,scatter(nf(:,1),nf(:,2),10,lvxy(ind,2),'filled');

corr(mpow(ind,2),lvxy(ind,2))

corr(mpow(ind,2),lvxy(ind,2))

figure,
subplot(121)
plot(nf(:,1),lvxy(ind,2),'.')
subplot(122)
plot(nf(:,2),lvxy(ind,2),'.r')


rot = 0.8;
nf = multiprod([unity(ppow(ind)),unity(lpow(ind))],[cos(rot),-sin(rot);sin(rot),cos(rot)],[2],[1,2]);
figure,plot(nf(:,1),nf(:,2),'.');
% $$$ figure,plot(nf(:,2),erf(ind),'.');
corr(nf(:,1),erf(ind))
corr(nf(:,2),erf(ind))

figure();
hist2([cpow(ind),erf(ind)],cpowBin,erfBin,'yprob');
caxis([0,0.1])
corr(cpow(ind),erf(ind))
figure();
hist2([cSpow(ind),erf(ind)],cSpowBin,erfBin,'yprob');
caxis([0,0.1])
corr(cSpow(ind),erf(ind))

figure();
subplot(121);
hist2([rpow(ind,1),erf(ind,1)],rpowBin,erfBin,'yprob');
caxis([0,0.1])
corr(rpow(ind,1),erf(ind,1))
subplot(122);
hist2([rpow(nind,1),erf(nind,1)],rpowBin,erfBin,'yprob');
caxis([0,0.1])
hold('on');plot(rpow(nind,1),erf(nind,1),'.')
corr(rpow(nind,1),erf(nind,1))

figure();
hist2([dpow(ind),erf(ind)],dpowBin,erfBin,'yprob');
caxis([0,0.1])
corr(dpow(ind),erf(ind))
figure();
hist2([dpow(ind),erd(ind)],dpowBin,erdBin,'yprob');
caxis([0,0.15])
corr(dpow(ind),erf(ind))
figure,
plot(dpow(ind),log10(vxy(ind)),'.');

% cpow is NOT indicative rad/lm influence of decoded postion at small time windows
% RECOMPUTING dc with 88 ms timewindow instead of 44ms
% 300 ms shows relationship during pause. Could head direction changes and 
% 44 ms window dc was computed with 2d rate maps over 4d - could have an effect
% 
% could phase dependent firing rate also be speed dependent? most likely.
%
