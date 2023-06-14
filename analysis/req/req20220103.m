% req20220103
%     Tags: realtime state segmentation lm rad ratio thetarc
%     Status: Active
%     Type: Analysis
%     Author: Justin Graboski
%     Final_Forms: NA
%     Project: General
%     Description: theta return current feature used for realtime ephys state segmentation
%

% Notes:
%     quarter second time window does not capture the feature of low <20Hz vs high >40 and <200 Hz as well
%     as the filtered means of the eigth of a second spectra.

MjgER2016_load_data();
sessionList = get_session_list_v2('MjgER2016');
sampleRate = 30;
trialId = 21;


Trial = Trials{trialId};
unitSubset = units{trialId};
meta = sessionList(trialId);

xyz = preproc_xyz(Trial,'trb',sampleRate);
fxyz = filter(xyz.copy(),'ButFilter',3,14,'low');
vxy = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
vxyz = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2,3]);
lvxy = copy(vxy);
lvxy.data(lvxy.data<=0.0001) = 0.0001;
lvxy.data = log10(lvxy.data);

figure();
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxyz)]./vxy.sampleRate,vxyz(:,2)/5+1,'r');
plot([1:size(vxyz)]./vxy.sampleRate,vxyz(:,1)/5+1,'g');
xlabel('Time (s)');


Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',meta.subject.channelGroup.theta);


flfp = diff(get(Trial.load('lfp',[meta.subject.channelGroup.thetarc]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2);
%flfp = diff(get(Trial.load('lfp',[1,8]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2);

% $$$ flfp = [diff(get(Trial.load('lfp',[33,40]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2),...
% $$$         diff(get(Trial.load('lfp',[41,48]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2),...
% $$$         diff(get(Trial.load('lfp',[49,56]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2),...
% $$$         diff(get(Trial.load('lfp',[57,64]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2)];

rlfp = Trial.load('lfp',[61]);

%rlfp & flfp
% $$$ 
% $$$ 
% $$$ figure,
% $$$ subplot(211);
% $$$ imagesc([1:size(fslfp,1)]./lfp.sampleRate,1:8,nunity(fslfp(:,:))')
% $$$ hold('on');
% $$$ plot([1:size(fslfp,1)]./lfp.sampleRate,nunity(diff(fslfp(:,[1,end]),1,2))+4.5,'m','LineWidth',2);
% $$$ plot([1:size(fslfp,1)]./lfp.sampleRate,sfslfp(:,1)+4.5,'r','LineWidth',2);
% $$$ colormap('jet');
% $$$ axis('xy');
% $$$ caxis([-3,3]);
% $$$ Lines([],4.5,'k');
% $$$ subplot(212);
% $$$ hold('on');
% $$$ plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
% $$$ ylim([1,9]);
% $$$ plot([1:size(vxyz)]./vxy.sampleRate,vxyz(:,2)/5+1,'r');
% $$$ plot([1:size(vxyz)]./vxy.sampleRate,vxyz(:,1)/5+1,'g');
% $$$ xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ 
% $$$ figure,plot(slfp(1:1e5,:))
% $$$ hold('on');
% $$$ plot(flfp(1:1e5,4),'m','LineWidth',2)








unitsInt = Trial.spk.get_unit_set(Trial,'interneurons');

int = Trial.load('spk', Trial.lfp.sampleRate, '', unitsInt, '');

% $$$ % LOAD theta phase
% $$$ pchan = [1,8,9,16,17,24,25,32,33,40,41,48,49,56,57,64];
% $$$ for c = 1:numel(pchan)
% $$$ phz = load_theta_phase(Trial,...
% $$$                        Trial.lfp.sampleRate,...
% $$$                        pchan(c),...
% $$$                        meta.subject.correction.thetaPhase);
% $$$ for ii = 1:numel(unitsInt),
% $$$ phzMean(c,ii) = circ_mean(phz(int(unitsInt(ii))));
% $$$ phzR(c,ii) = circ_r(phz(int(unitsInt(ii))));
% $$$ end
% $$$ end


%sort_interneurons_by_tpp(Trial,meta)

phz = load_theta_phase(Trial,...
                       Trial.lfp.sampleRate,...
                       meta.subject.channelGroup.theta,...
                       meta.subject.correction.thetaPhase);

% ORDER interneurons based on theta phase preference
intPhzPref = zeros([numel(unitsInt),1]);
for ii = 1:numel(unitsInt),
    intPhzPref(ii) = circ_mean(phz(int(unitsInt(ii))));
    intPhzR(ii) = circ_r(phz(int(unitsInt(ii))));
end
intPhzPref(intPhzPref<0) = intPhzPref(intPhzPref<0)+2*pi;
[mpv,mpi] = sort(intPhzPref,'descend');
mrv = intPhzR(mpi);
unitsInt = unitsInt(mpi);


%     5.22099827436091       9                     4.84224310697765    11
%     4.43219199687354      10                     4.79422023722339     7
%     3.93546757990167       7                     4.64561500480617    13
%     3.83702698028012      14                     3.65005269713263     1
%     3.78686563321746      13                     2.93272011688267    14
%     3.54861455639414       6                     2.69081170479908    10
%     3.44208801573508       8                     2.54160077488562    12
%     3.36913474655993      12                      2.3985019436955     6
%     3.02889248200655       5                     2.38542610007183     4
%     2.85983023151636       3                     2.29636991494833     9
%     2.68403989509002      11                     2.23340627953952     3
%     2.59567869954665       2                     2.19082410252469     8
%     2.59076386975001       4                     1.98523738734222     2
%     2.19652065581479       1                     1.54319131641419     5

ufr = Trial.load('ufr', lvxy, [], unitsInt, 0.12, 'boxcar');
% $$$ fufr = Trial.load('ufr', lvxy, [], unitsInt, 0.5, 'gauss');
fufr = Trial.load('ufr', lvxy, [], unitsInt, 1.5, 'gauss');
ufrp = Trial.load('ufr', lvxy, [], unitsSubset, 0.12, 'boxcar');
%ufrl = Trial.load('ufr', lvxy, [], unitsInt, 0.24, 'boxcar');
% $$$ 
% $$$ phzThresh = 3;
% $$$ figure
% $$$ subplot(211);hold('on');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,RectFilter(sum(ufr.data,2)./numel(unitsInt),9,3));
% $$$ % $$$ plot([1:size(vxy)]./vxy.sampleRate,RectFilter(sum(ufr(:,mpv>phzThresh),2)./sum(mpv>phzThresh),9,3));
% $$$ % $$$ plot([1:size(vxy)]./vxy.sampleRate,RectFilter(sum(ufr(:,mpv<3),2)./sum(mpv<3),9,3));
% $$$ plot([1:size(vxy)]./vxy.sampleRate,log10(RectFilter(sum(ufr(:,mpv>phzThresh),2)./sum(mpv>phzThresh),9,3)./RectFilter(sum(ufr(:,mpv<2.5),2)./sum(mpv<2.5),9,3)));
% $$$ plot([1:size(vxy)]./vxy.sampleRate,RectFilter(sum(ufr.data,2)./numel(unitsInt),9,3)- ...
% $$$      log10(RectFilter(sum(ufr(:,mpv>phzThresh),2)./sum(mpv>phzThresh),9,3)./RectFilter(sum(ufr(:,mpv<3),2)./sum(mpv<3),9,3)));
% $$$ Lines([],0,'k');
% $$$ Lines([],0.5,'k');
% $$$ subplot(212);
% $$$ hold('on');
% $$$ plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
% $$$ ylim([1,9]);
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
% $$$ xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ 
% $$$ 
% $$$ 

figure,
udcor = [];
for u = 1:numel(unitsInt),
    subplot(2,7,u);
% $$$     sufr = fufr(dc.ind(((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>1)),u);
% $$$     dcom = sqrt(sum(dc.ecom((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>1,:).^2,2));
    sufr = ufr(dc.ind((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1)),u);
    dcom = sqrt(sum(dc.ecom((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1),:).^2,2));
    nind = nniz(sufr) & nniz(dcom) & sufr>0.0001;
    plot(sufr(nind),dcom(nind),'.');
    udcor(u) = corr(sufr(nind),dcom(nind));
end


figure,
uvcor = [];
for u = 1:numel(unitsInt),
    subplot(2,7,u);
    sufr = fufr(dc.ind(((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2)),u);
    dcom = sqrt(sum(dc.ecom((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2,:).^2,2));
% $$$     sufr = ufr(dc.ind((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1)),u);
% $$$     dcom = lvxy(dc.ind((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1)),2);
    nind = nniz(sufr) & nniz(dcom) & sufr>0.0001;
    plot(sufr(nind),dcom(nind),'.');
    uvcor(u) = corr(sufr(nind),dcom(nind));
end




stcm = stc2mat(Trial.stc,xyz,{'theta','rear','loc','pause','sit','groom'});
gufr = Trial.load('ufr',xyz,int,unitsInt,0.125,'gauss');
%ind = (stcm(:,3)==3|stcm(:,4)==4);
%ind = (stcm(:,5)==5);
ind = (stcm(:,4)==4) & stcm(:,1)~=1;
suI = gufr(ind,:);
svI = lvxy(ind,:);
nClust = 2;
[idx,C,SUMD,D] = kmeans(suI,nClust,'Distance','correlation');

kCov = zeros(size(gufr,2),size(gufr,2),nClust);
for vind = 1:nClust,
    kCov(:,:,vind) = (bsxfun(@minus,suI(idx==vind,:),median(suI(idx==vind,:)))' ...
                       *bsxfun(@minus,suI(idx==vind,:),median(suI(idx==vind,:)))) ...
                      ./sum(idx==vind);
end

velBins = linspace(-2.5,1.8,20);
velInds = discretize(lvxy(:,2),velBins);
figure,
for vind = 1:nClust,
    subplot2(nClust,2,vind,1);
        imagesc(kCov(:,:,vind));
    subplot2(nClust,2,vind,2);
    histogram(svI(idx==vind,2),velBins);
end

out = [];
for vind =  1:nClust
    out(:,vind) = -.5*log(det(kCov(:,:,vind)))-0.5*(multiprod(bsxfun(@minus,gufr(:,:),median(suI(idx==vind,:))),...
                  multiprod(inv(kCov(:,:,vind)),bsxfun(@minus,gufr(:,:),median(suI(idx==vind,:))),[1,2],[2]),2,2));
end
figure,plot(out)


% $$$ figure();
% $$$ ind = (stcm(:,5)==5);
% $$$ plot(sum(gufr(ind,:),2)./14,lvxy(ind,2),'.')
% $$$ hold('on');
% $$$ ind = (stcm(:,3)==3|stcm(:,4)==4);
% $$$ plot(sum(gufr(ind,:),2)./14,lvxy(ind,2),'.r')
% $$$ 

udThresh = 0.1;
udThresh = 0.0;
irat = copy(vxy);
irat.data = log10(  sum(fufr(:,udcor>0&abs(udcor)>udThresh),2)./sum(udcor>0&abs(udcor)>udThresh) ...
                  ./sum(fufr(:,udcor<0&abs(udcor)>udThresh),2)./sum(udcor<0&abs(udcor)>udThresh));
% $$$ 
% $$$ frat = copy(vxy);
% $$$ frat.data = log10(  RectFilter(sum(ufr(:,udcor>0&abs(udcor)>udThresh),2)./sum(udcor>0&abs(udcor)>udThresh),21,3) ...
% $$$                   ./RectFilter(sum(ufr(:,udcor<0&abs(udcor)>udThresh),2)./sum(udcor<0&abs(udcor)>udThresh),21,3));
% $$$ 
% $$$ uvThresh = 0.0;
% $$$ vrat = copy(vxy);
% $$$ vrat.data = log10(  RectFilter(sum(ufr(:,udcor>0&abs(udcor)>uvThresh),2)./sum(udcor>0&abs(udcor)>uvThresh),9,3) ...
% $$$                   ./RectFilter(sum(ufr(:,udcor<0&abs(udcor)>uvThresh),2)./sum(udcor<0&abs(udcor)>uvThresh),9,3));
% $$$ 
% $$$ rat = copy(vxy);
% $$$ nrat.data = log(RectFilter(p1,21,3)./RectFilter(p2,21,3));
% $$$ 
% $$$ prat = copy(vxy);
% $$$ prat.data = log(RectFilter(p1,31,3)+RectFilter(p2,31,3));
% $$$ 
% $$$ 
% $$$ rdThresh = randperm(14);
% $$$ irat.data = log10(  RectFilter(sum(ufr(:,rdThresh(1:7)),2)./7,9,3) ...
% $$$                   ./RectFilter(sum(ufr(:,rdThresh(8:14)),2)./7,9,3));
% $$$ 
phzThresh = 2.5;
%phzThresh = 3;
trat = copy(vxy);
trat.data = log10(  (sum(fufr(:,mpv>phzThresh),2)./sum(mpv>phzThresh))  ...
                  ./(sum(fufr(:,mpv<phzThresh),2)./sum(mpv<phzThresh)));

figure
plot(unity(hfet.data),unity(sum(fufr.data,2)),'.');

figure
subplot(211);hold('on');
%plot([1:size(vxy)]./vxy.sampleRate,unity(sum(fufr.data,2)));
plot([1:size(vxy)]./vxy.sampleRate,trat.data);
Lines([],0,'k');
Lines([],0.5,'k');
subplot(212);
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

% DECODING vars -----------------------------------------------------------------------------

dc = accumulate_decoding_vars(Trial,                               ...
                              unitSubset,                          ...
                              meta.subject.channelGroup.theta,     ...
                              meta.subject.correction.thetaPhase,  ...
                              meta.subject.correction.headYaw,     ...
                              meta.subject.correction.headBody);


% $$$ figure
% $$$ plot(glfp(1:1e6,:))

% $$$ figure
% $$$ plot(diff(glfp(1:1e6,[1,2]),1,2))
% $$$ hold('on')
% $$$ plot(glfp(1:1e6,1))

% $$$ elfp = copy(lfp);
% $$$ elfp.data = diff(glfp(1:1e6,[1,3]),1,2);
% $$$ elfp = copy(lfp);
% $$$ elfp.data = RectFilter(glfp,3,3);
% $$$ specArgsTheta = struct('nFFT',2^8,...
% $$$                   'Fs',  elfp.sampleRate,...
% $$$                   'WinLength',2^7,...
% $$$                   'nOverlap',2^7*0.875,...
% $$$                   'NW',3,...
% $$$                   'Detrend',[],...
% $$$                   'nTapers',[],...
% $$$                   'FreqRange',[1,250]);
% $$$ [mys,mfs,mts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);
% $$$ 
% $$$ mxyz = preproc_xyz(Trial,'trb');
% $$$ mxyz.resample(mys);
% $$$ mfxyz = filter(mxyz.copy(),'ButFilter',3,14,'low');
% $$$ mvxy = vel(filter(mxyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
% $$$ mlvxy = copy(mvxy);
% $$$ mlvxy.data(mlvxy.data<=0.0001) = 0.0001;
% $$$ mlvxy.data = log10(mlvxy.data);
% $$$ 
% $$$ 
% $$$ mNBins = 31;
% $$$ mlvxyBinEdgs = linspace(-3,2,mNBins);
% $$$ mlvxyBinCntr = mean([mlvxyBinEdgs(1:end-1);mlvxyBinEdgs(2:end)]);
% $$$ mlvxyBinInds = discretize(mlvxy(:,2),mlvxyBinEdgs);
% $$$ 
% $$$ mspc = [];
% $$$ for v = 1:mNBins-1
% $$$     tspc = log10(mys(v == mlvxyBinInds,:,4));
% $$$     mspc(:,v) = mean(tspc(nniz(tspc),:));
% $$$ end
% $$$ 
% $$$ figure;imagesc(mfs,mlvxyBinCntr,bsxfun(@rdivide,mspc,sum(mspc,2))');axis('xy');colormap('jet');colorbar();
% $$$ figure;imagesc(mfs,mlvxyBinCntr,mspc');axis('xy');colormap('jet');colorbar();



fwin = gausswin(2^10);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = rlfp(:,1);
elfp.filter('ButFilter', 4, [20], 'low');
lpfet = copy(elfp);
lpfet.data = log10(conv(lpfet.data.^2,fwin,'same'));
lpfet.resample(xyz);

% $$$ fwin = gausswin(128);
% $$$ elfp = copy(lfp);
% $$$ elfp.data = flfp(:,1);
% $$$ elfp.filter('ButFilter', 4, [1,20], 'bandpass');
% $$$ lfet = copy(elfp);
% $$$ lfet.data = log10(conv(lfet.data.^2,fwin,'same'));
% $$$ lfet.resample(xyz);
% $$$ lfet.filter('ButFilter', 4, [1], 'low');

fwin = gausswin(64);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = flfp(:,1);
%elfp.filter('ButFilter', 4, [50,200], 'bandpass');
%elfp.filter('ButFilter', 4, [50,100], 'bandpass');
%elfp.filter('ButFilter', 4, [125,200], 'bandpass');
elfp.filter('ButFilter', 4, [75,125], 'bandpass');
hfet = copy(elfp);
hfet.data = log10(conv(hfet.data.^2,fwin,'same'));
hfet.filter('ButFilter', 4, [1], 'low');
hfet.resample(xyz);

% ifet 
ifet = copy(fufr);
ifet.data = log10(sum(fufr.data,2)./size(fufr,2));

% rfet
fwin = gausswin(2^8);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = rlfp(:,1);
elfp.filter('ButFilter', 4, [140,220], 'bandpass');
rfet = copy(elfp);
rfet.data = log10(conv(rfet.data.^2,fwin,'same'));
rfet.resample(xyz);

%dfet
fwin = gausswin(2^12);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = rlfp(:,1);
elfp.filter('ButFilter', 4, [4], 'low');
dfet = copy(elfp);
dfet.data = log10(conv(dfet.data.^2,fwin,'same'));
ttdfet = copy(dfet);
dfet.resample(xyz);
%tfet
elfp = copy(lfp);
elfp.data = rlfp(:,1);
elfp.filter('ButFilter', 4, [5,10], 'bandpass');
tfet = copy(elfp);
tfet.data = log10(conv(tfet.data.^2,fwin,'same'));
ttfet = copy(tfet);
tfet.resample(xyz);

tdfet = copy(ttfet);
tdfet.data = log10(ttfet.data./ttdfet.data);
tdfet.resample(xyz);


% drfet
fwin = gausswin(2^12);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = flfp(:,1);
elfp.filter('ButFilter', 4, [4], 'low');
drfet = copy(elfp);
drfet.data = log10(conv(drfet.data.^2,fwin,'same'));
ttdrfet = copy(drfet);
drfet.resample(xyz);
% trfet
elfp = copy(lfp);
elfp.data = flfp(:,1);
elfp.filter('ButFilter', 4, [5,10], 'bandpass');
trfet = copy(elfp);
trfet.data = log10(conv(trfet.data.^2,fwin,'same'));
ttrfet = copy(trfet);
trfet.resample(xyz);

tdrfet = copy(ttdrfet);
tdrfet.data = log10(ttrfet.data./ttdrfet.data);
tdrfet.resample(xyz);

xts = [1:size(trfet,1)]./sampleRate;
figure,
subplot(211);hold('on');
plot(xts,tfet.data)
plot(xts,trfet.data)
Lines([],-0.05,'k')
subplot(212);
    hold('on');
    plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


figure,
ind = stcm(:,1)==1&(stcm(:,3)==3|stcm(:,4)==4|stcm(:,2)==2);
%ind = [Trial.stc{'w+p+r+n+s'}];
subplot(221);
hist2([trfet(ind),tfet(ind)],linspace(5,8,50),linspace(5,7,50));
grid('on');set(gca,'GridColor','w');
subplot(222);
hist2([dfet(ind),tfet(ind)],linspace(4,7,50),linspace(4.5,7,50));
grid('on');set(gca,'GridColor','w');
ind = stcm(:,1)==1 & stcm(:,5)==5;
subplot(223);
hist2([trfet(ind),tfet(ind)],linspace(5,8,50),linspace(5,7,50));
grid('on');set(gca,'GridColor','w');
subplot(224);
hist2([dfet(ind),tfet(ind)],linspace(4,7,50),linspace(4.5,7,50));
grid('on');set(gca,'GridColor','w');

% $$$ drfet = rfet.copy();
% $$$ drfet.data = diff([0;rfet.data]);
% $$$ ind = [stc{'s'}];
% $$$ ind = [stc{'w+n+p+r'}];
% $$$ %plot(rfet(ind),lpfet(ind),'.');
% $$$ sum(rfet(ind)<5.5&rfet(ind)>4.25 &abs(drfet(ind))<0.2)
% $$$ sum((rfet(ind)>5.5|rfet(ind)<4.25)& abs(drfet(ind))<0.2)
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ ind = [stc{'s'}];
% $$$ %plot(rfet(ind),lpfet(ind),'.');
% $$$ plot(rfet(ind),diff([0;rfet(ind)]),'.');
% $$$ ind = [stc{'w+n+p+r'}];
% $$$ %plot(rfet(ind),lpfet(ind),'.');
% $$$ plot(rfet(ind),diff([0;rfet(ind)]),'.');
% $$$ 
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ ind = [stc{'s'}];
% $$$ %plot(lpfet(ind),lpfet(ind),'.');
% $$$ plot(lpfet(ind),diff([0;lpfet(ind)]),'.');
% $$$ ind = [stc{'w+n+p+r'}];
% $$$ %plot(lpfet(ind),lpfet(ind),'.');
% $$$ plot(lpfet(ind),diff([0;lpfet(ind)]),'.');
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ ind = [stc{'s'}];
% $$$ plot(rfet(ind),lpfet(ind),'.');
% $$$ %plot(rfet(ind),diff([0;rfet(ind)]),'.');
% $$$ ind = [stc{'w+n+r'}];
% $$$ plot(rfet(ind),lpfet(ind),'.');



% $$$ fwin = gausswin(2^10);
% $$$ fwin = fwin./sum(fwin);
% $$$ elfp = copy(lfp);
% $$$ elfp.data = flfp(:,1);
% $$$ elfp.filter('ButFilter', 4, [5,12], 'bandpass');
% $$$ dfet = copy(elfp);
% $$$ dfet.data = log10(conv(dfet.data.^2,fwin,'same'));
% $$$ dfet.resample(xyz);

% $$$ sfet = copy(sfslfp);
% $$$ sfet.resample(xyz);
% $$$ dsfet = sfet.copy;
% $$$ dsfet.data = nunity(circshift(dsfet.data,1)-dsfet.data);


% $$$ 
% $$$ figure();
% $$$ subplot(211);
% $$$     hold('on');
% $$$     plot([1:size(hfet)]./lfet.sampleRate,hfet.data)
% $$$     plot([1:size(lfet)]./hfet.sampleRate,lfet.data)
% $$$ subplot(212);
% $$$     hold('on');
% $$$     plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
% $$$     ylim([1,9]);
% $$$     plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
% $$$     plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
% $$$     xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ 
% $$$ ind = [Trial.stc{'s+w+n+r+p'}];
% $$$ ind = [Trial.stc{'w+n+p'}];
% $$$ [Uf,Sf,Vf] = svd([lfet(ind)-mean(lfet(ind)),hfet(ind)-mean(hfet(ind))],0);
% $$$ Vf
% $$$ % $$$ Vf = [0.771753449181619,        -0.635921861297655;
% $$$ % $$$      0.635921861297655,         0.771753449181619];
% $$$ Vf = [-0.83534705515978,         0.549722927879022; ...
% $$$         -0.549722927879022,         -0.83534705515978];
% $$$ Vf = [0.916079244726519,        -0.400997278521052;...
% $$$          0.400997278521052,         0.916079244726519];
% $$$ 
% $$$ aind = [Trial.stc{'w+p'}];
% $$$ aind.cast('TimeSeries');
% $$$ aind.resample(xyz);
% $$$ Vv = [];
% $$$ for v = 1:19,
% $$$     ind = v==velInds & logical(aind.data);
% $$$     ml(v) = mean(lfet(ind));
% $$$     hl(v) = mean(hfet(ind));
% $$$     [Uv,Sv,Vv(:,:,v)] = svd([lfet(ind)-ml(v),hfet(ind)-hl(v)],0);
% $$$ end    
% $$$ figure,plot(abs(sq(Vv(:,1,:)))')
% $$$ figure,plot(ml,hl,'-+');
% $$$ 
% $$$ mlfet = lfet.copy();  mlfet.data = mlfet.data - mean(mlfet(ind));
% $$$ mhfet = hfet.copy();  mhfet.data = mhfet.data - mean(mhfet(ind));
% $$$ 
% $$$ mlhfet = lfet.copy();
% $$$ mlhfet.data = multiprod([mlfet.data,mhfet.data],Vf,2,[1,2]);
% $$$ 
% $$$ inds = {[Trial.stc{'s+w+r+n+p'}],[Trial.stc{'s'}],[Trial.stc{'w'}],[Trial.stc{'p'}]};
% $$$ figure,
% $$$ normf = '';      clim = [0,300];
% $$$ for s = 1:numel(inds),
% $$$     subplot(2,2,s);
% $$$     ind = inds{s};
% $$$     hist2(multiprod([mlfet(ind),mhfet(ind)],Vf,2,[1,2]),linspace(-2,2,50),linspace(-1,1,50),normf);
% $$$     Lines([],0,'m');
% $$$     Lines(0,[],'m');
% $$$     caxis(clim);
% $$$ end




elfp = copy(lfp);
elfp.data = flfp(:,1);
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  elfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[0.5,25]);
[lys,lfs,lts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);


elfp = copy(lfp);
elfp.data = rlfp(:,1);
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  elfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[0.5,25]);
[rys,rfs,rts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);

figure
subplot(311);
imagesc(lts,lfs,log10(lys.data)');
axis('xy');
colormap('jet');
caxis([3,5.5])
subplot(312);
imagesc(rts,rfs,log10(rys.data)');
axis('xy');
colormap('jet');
caxis([4,6.5])
subplot(313);
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

% $$$ % $$$ elfp = copy(lfp);
% $$$ % $$$ elfp.data = diff(glfp(1:1e6,[1,3]),1,2);
% $$$ elfp = copy(lfp);
% $$$ specArgsTheta = struct('nFFT',2^11,...
% $$$                   'Fs',  elfp.sampleRate,...
% $$$                   'WinLength',2^1s0,...
% $$$                   'nOverlap',2^10*0.875,...
% $$$                   'NW',3,...
% $$$                   'Detrend',[],...
% $$$                   'nTapers',[],...
% $$$                   'FreqRange',[1,30]);
% $$$ [rys,rfs,rts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);


% $$$ elfp = copy(lfp);
% $$$ specArgsTheta = struct('nFFT',2^8,...
% $$$                   'Fs',  elfp.sampleRate,...
% $$$                   'WinLength',2^7,...
% $$$                   'nOverlap',2^7*0.875,...
% $$$                   'NW',3,...
% $$$                   'Detrend',[],...
% $$$                   'nTapers',[],...
% $$$                   'FreqRange',[50,280]);
% $$$ [tys,tfs,tts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure
% $$$ subplot(411);
% $$$ imagesc(mts,mfs,log10(mean(mys.data,3))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(412);
% $$$ imagesc(mts,mfs,log10(std(mys.data,[],3))');
% $$$ %imagesc(tts,tfs,log10(tys.data)');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(413);
% $$$ imagesc(mts,mfs,log10(mean(mys.data,3))'.^2./log10(std(mys.data,[],3))');
% $$$ %imagesc(tts,tfs,log10(tys.data)');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(414);
% $$$ hold('on');
% $$$ plotSTC(Trial.stc,1,'text',states,'krggbb');
% $$$ ylim([1,9]);
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
% $$$ xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ 
% $$$ figure
% $$$ subplot(511);
% $$$ imagesc(mts,mfs,nunity(log10(mys(:,:,1)))');
% $$$ hold('on');
% $$$ plot((1:size(elfp,1))./elfp.sampleRate,nunity(elfp.data(:,3))*10+30,'w');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(512);
% $$$ imagesc(mts,mfs,nunity(log10(mys(:,:,2)))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(513);
% $$$ imagesc(mts,mfs,nunity(log10(mys(:,:,3)))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(514);
% $$$ imagesc(mts,mfs,nunity(log10(mys(:,:,4)))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(515);
% $$$ hold('on');
% $$$ plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
% $$$ ylim([1,9]);
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
% $$$ xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ 
% $$$ figure
% $$$ subplot(511);
% $$$ imagesc(mts,mfs,(log10(mys(:,:,2)))');
% $$$ hold('on');
% $$$ plot((1:size(elfp,1))./elfp.sampleRate,nunity(elfp.data(:,3))*10+30,'w');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ grid('on')
% $$$ caxis([3.75,5.25])
% $$$ subplot(512);
% $$$ imagesc(mts,mfs,(log10(mys(:,:,3)))');
% $$$ %imagesc(mts,mfs,(log10(mys(:,:,2)))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ grid('on')
% $$$ caxis([3.75,5.25])
% $$$ subplot(513);
% $$$ imagesc(mts,mfs,(log10(mys(:,:,4)))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ grid('on')
% $$$ caxis([4,5.25])
% $$$ subplot(514);
% $$$ hold('on');
% $$$ imagesc(mts,mfs,bsxfun(@rdivide,RectFilter(log10(mys(:,:,4)),3,1),sum(RectFilter(log10(mys(:,:,4)),3,1),2))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ grid('on');
% $$$ caxis([0.019,0.022]);
% $$$ plot((1:size(elfp,1))./elfp.sampleRate,nunity(elfp.data(:,3))*10+30,'w');
% $$$ subplot(515);
% $$$ hold('on');
% $$$ plotSTC(Trial.stc,1,'text',states,'krggbb');
% $$$ ylim([1,9]);
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,(tdRatio(:,1)+1)*5,'r');
% $$$ xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ figure
% $$$ subplot(211);
% $$$ %plot(mts,mean((log10(mys(:,:,3))),2)');
% $$$ hold('on');
% $$$ plot(mts,mean((log10(mys(:,mfs<200&mfs>40,3))),2));
% $$$ plot(mts,mean(RectFilter(log10(mys(:,mfs<200&mfs>40,3)),11,3),2));
% $$$ %plot(mts,mean(RectFilter(log10(mys(:,mfs<100,4)),7,3),2));
% $$$ plot(mts,mean(log10(mys(:,mfs<20,3)),2),'r');
% $$$ plot(mts,mean(RectFilter(log10(mys(:,mfs<20,3)),11,3),2),'k');
% $$$ plot(mts,mean(RectFilter(log10(mys(:,mfs<20,3)),21,3),2),'c');
% $$$ Lines([],4,'k');
% $$$ Lines([],4.5,'k');
% $$$ Lines([],5,'k');
% $$$ subplot(212);
% $$$ hold('on');
% $$$ plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
% $$$ ylim([1,9]);
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
% $$$ xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ 
% $$$ fmys = copy(mys);
% $$$ fmys.data =  RectFilter(log10(mys(:,:,:)),21,3);
% $$$ 
% $$$ rcInd = 1;
% $$$ tpm = copy(mys);
% $$$ tpm.data = mean(fmys(:,mfs<20,rcInd),2);
% $$$ tpm.resample(vxy);
% $$$ 
% $$$ mpm = copy(mys);
% $$$ mpm.data = mean(fmys(:,mfs>140 & mfs<200,rcInd),2);
% $$$ mpm.resample(vxy);
% $$$ 
% $$$ ipm = copy(mys);
% $$$ ipm.data = mean(fmys(:,mfs>40 & mfs<200,rcInd),2);
% $$$ ipm.resample(vxy);
% $$$ 
% $$$ gpm = copy(mys);
% $$$ gpm.data = mean(fmys(:,mfs>50 & mfs<100,rcInd),2);
% $$$ gpm.resample(vxy);
% $$$ 
% $$$ hpm = copy(mys);
% $$$ hpm.data = mean(fmys(:,mfs>200,rcInd),2);
% $$$ hpm.resample(vxy);
% $$$ 
% $$$ tdpow = copy(rys);
% $$$ tdpow.resample(vxy);
% $$$ tdpow.data = log(abs(log(mean(tdpow(:,rfs>5 & rfs<12),2))-log(mean(tdpow(:,rfs<5 | (rfs>12 & rfs<15)),2))));
% $$$ 
% $$$ tpow = copy(rys);
% $$$ tpow.resample(vxy);
% $$$ tpow.data = log(mean(tpow(:,rfs>5 & rfs<12),2));
% $$$ 
% $$$ 
% $$$ % ERROR in stc s+t is not correct
% $$$ figure
% $$$ normf = '';      clim = [0,500];
% $$$ % $$$ normf = '';      clim = [0,100];
% $$$ % $$$ normf = 'xprob'; clim = [0,0.15];
% $$$ ind = [Trial.stc{'s+w+p+r+n'}];
% $$$ %ind = [Trial.stc{'w+p+r+n'}];
% $$$ %ind = [Trial.stc{'r'}];
% $$$ %ind = [Trial.stc{'s'}];
% $$$ subplot(251);
% $$$ hist2([tpm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(252);
% $$$ hist2([gpm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(253);
% $$$ hist2([mpm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(254);
% $$$ hist2([ipm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(255);
% $$$ hist2([hpm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(256);
% $$$ hist2([tpm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(257);
% $$$ hist2([gpm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(258);
% $$$ hist2([mpm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(259);
% $$$ hist2([ipm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(2,5,10);
% $$$ hist2([hpm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ 
% $$$ 
% $$$ 
% $$$ ind = [Trial.stc{'s+w+r+n+p'}];
% $$$ tpmm = copy(tpm); tpmm.data = tpmm.data-mean(tpmm(ind));
% $$$ ipmm = copy(ipm); ipmm.data = ipmm.data-mean(ipmm(ind));
% $$$ %mpmm = copy(mpm); mpmm.data = mpmm.data-mean(mpmm(ind));
% $$$ 
% $$$ 
% $$$ %ind = [Trial.stc{'s'}];
% $$$ figure
% $$$ hist2([tpmm(ind),ipmm(ind)],linspace(-1.25,1.25,50),linspace(-1.25,1.25,50),normf);
% $$$ caxis(clim)
% $$$ 
% $$$ ind = [Trial.stc{'s+w+r+n+p'}];
% $$$ [U,S,V] = svd([tpm(ind)-mean(tpm(ind)),ipm(ind)-mean(ipm(ind))],0);
% $$$ V = [0.771753449181619,        -0.635921861297655;
% $$$      0.635921861297655,         0.771753449181619];
% $$$ 
% $$$ 
% $$$ %[U,S,V] = svd([tpm(ind)-mean(tpm(ind)),mpm(ind)-mean(mpm(ind))],0);
% $$$ % $$$ figure,
% $$$ % $$$ plot(U(:,1),U(:,2),'.')
% $$$ % $$$ 
% $$$ % $$$ figure
% $$$ % $$$ hist2([U(:,1),U(:,2)],linspace(-3.2e-3,-2.2e-3,50),linspace(-8e-3,8e-3,50))
% $$$ 
% $$$ 
% $$$ figure,
% $$$ normf = '';
% $$$ subplot(221);
% $$$ ind = [Trial.stc{'w+r+n+p+s'}];
% $$$ hist2(multiprod([tpmm(ind),ipmm(ind)],V,2,[1,2]),linspace(-1.25,1.25,50),linspace(-.5,.5,50),normf);
% $$$ %line([-6.85,-5.5],[-0.4,0.4],'Color','r');
% $$$ line([-.7,0.8],[0.5,-0.5],'Color','r');
% $$$ subplot(222);
% $$$ ind = [Trial.stc{'p'}];
% $$$ hist2(multiprod([tpmm(ind),ipmm(ind)],V,2,[1,2]),linspace(-1.25,1.25,50),linspace(-.5,.5,50),normf);
% $$$ %line([-6.85,-5.5],[-0.4,0.4],'Color','r');
% $$$ line([-.7,0.8],[0.5,-0.5],'Color','r');
% $$$ subplot(223);
% $$$ ind = [Trial.stc{'w+r+n'}];
% $$$ hist2(multiprod([tpmm(ind),ipmm(ind)],V,2,[1,2]),linspace(-1.25,1.25,50),linspace(-.5,.5,50),normf);
% $$$ %line([-6.85,-5.5],[-0.4,0.4],'Color','r');
% $$$ line([-.7,0.8],[0.5,-0.5],'Color','r');
% $$$ subplot(224);
% $$$ ind = [Trial.stc{'s'}];
% $$$ hist2(multiprod([tpmm(ind),ipmm(ind)],V,2,[1,2]),linspace(-1.25,1.25,50),linspace(-.5,.5,50),normf);
% $$$ %line([-6.85,-5.5],[-0.4,0.4],'Color','r');
% $$$ line([-.7,0.8],[0.5,-0.5],'Color','r');
% $$$ 
% $$$ 
% $$$ dc = accumulate_decoding_vars(Trial,                               ...
% $$$                               units{trialId},                      ...
% $$$                               sessionList(trialId).thetaRefGeneral,...
% $$$                               phzCorrection(trialId),              ...
% $$$                               headRotation{trialId},               ...
% $$$                               hbangCorrection{trialId});
% $$$ 
% $$$ figure
% $$$ subplot(141);hold('on');
% $$$     ind = [Trial.stc{'s'}]; 
% $$$     set(histogram(lpfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'); 
% $$$     ind = [Trial.stc{'w'}]; 
% $$$     set(set(histogram(lpfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     xlim([5.5,7.5]);
% $$$ subplot(142);hold('on');
% $$$     ind = [Trial.stc{'s'}]; 
% $$$     set(histogram(tdfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'); 
% $$$     ind = [Trial.stc{'w'}]; 
% $$$     set(set(histogram(tdfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     xlim([-0.1,0.2]);
% $$$ subplot(143);hold('on');
% $$$     ind = [Trial.stc{'s'}]; 
% $$$     set(histogram(rfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'); 
% $$$     ind = [Trial.stc{'w'}]; 
% $$$     set(set(histogram(rfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     xlim([4,6]);
% $$$ subplot(144);hold('on');
% $$$     ind = [Trial.stc{'s'}]; 
% $$$     set(histogram(dfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'); 
% $$$     ind = [Trial.stc{'w'}]; 
% $$$     set(set(histogram(dfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     xlim([5,7.5]);

    
% $$$ 
% $$$ mfufr = copy(fufr);
% $$$ mfufr.data = sum(fufr(ind,:),2)./size(fufr,2);
% $$$ figure
% $$$ fufrSteps = 0:4:16;
% $$$ for f = 1:size(afet,2)
% $$$ subplot(1,size(afet,2),f);hold('on');
% $$$     for u = 1:numel(fufrSteps)-1
% $$$         ind = mfufr(:,1)>fufrSteps(u) & mfufr(:,1)<fufrSteps(u+1);
% $$$         set(set(histogram(afet(ind,f)-afetMeanW(f),linspace([afetLims(f,:),100]),'Normalization', ...
% $$$                           'probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     end
% $$$     ind = [Trial.stc{'w'}]; 
% $$$     set(set(histogram(afet(ind,f)-afetMeanW(f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     xlim(afetLims(f,:));
% $$$ end    


% $$$ figure,
% $$$ ind = ':';
% $$$ hist2([afet(ind,3),sum(fufr(ind,:),2)./size(ufr,2)],25,25);


afet = copy(lpfet);    
%afet.data = [lpfet(:,1),tdfet(:,1),rfet(:,1)];
%afet.data = [lpfet(:,1),tdfet(:,1),rfet(:,1),dfet(:,1),mlhfet(:,1),mlhfet(:,2)];
%afet.data = [lpfet(:,1),tdfet(:,1),rfet(:,1),dfet(:,1),lfet(:,1),hfet(:,1)];
%afet.data = [lpfet(:,1),tdfet(:,1),rfet(:,1),lfet(:,1),hfet(:,1)];
%afet.data = [lpfet(:,1),tdfet(:,1),tdrfet(:,1),rfet(:,1),lfet(:,1),hfet(:,1)];
%afet.data = [lpfet(:,1),tdfet(:,1),tdrfet(:,1),rfet(:,1),tfet(:,1),trfet(:,1),lfet(:,1),hfet(:,1)];
%afet.data = [dfet(:,1), lpfet(:,1), rfet(:,1), trfet(:,1), hfet(:,1), trat(:,1), ifet(:,1)];
afet.data = [dfet(:,1), lpfet(:,1), rfet(:,1), trfet(:,1), hfet(:,1)];
afetMeanW = median(afet(Trial.stc{'w'},:));
afet.data = bsxfun(@minus,afet.data,afetMeanW);
afetLims = [-1,1.5; -1,2;  -0.2,0.2; -1.5,1.25; -1.5,1.5; -1.5,1.5;-1.5,0.5; -1,0.8; -1.5,0.25  ];




figure();
for a = 1:size(afet,2)
subplot(2,4,a);
hold('on');
ind = stcm(:,1)==1&(stcm(:,3)==3|stcm(:,4)==4|stcm(:,2)==2);
set(histogram(afet(ind,a),linspace([afetLims(a,:),50]),'Normalization','probability'),'EdgeColor','none'); 
ind = stcm(:,1)==1 & stcm(:,5)==5;
set(histogram(afet(ind,a),linspace([afetLims(a,:),50]),'Normalization','probability'),'EdgeColor','none'); 
ind = stcm(:,1)~=1 & stcm(:,5)==5;
set(histogram(afet(ind,a),linspace([afetLims(a,:),50]),'Normalization','probability'),'EdgeColor','none'); 
end

%afetLims = [-1,1.5;  -0.175,0.1; -1,2; -1.5,1.25; -1.5,1.5; -1.5,1.5,];
% $$$ 
% $$$ 
% $$$ figure
% $$$ for f = 1:size(afet,2)
% $$$ subplot(1,size(afet,2),f);hold('on');
% $$$     ind = [Trial.stc{'s'}]; 
% $$$     set(histogram(afet(ind,f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'); 
% $$$     ind = [Trial.stc{'w'}]; 
% $$$     set(set(histogram(afet(ind,f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     ind = [Trial.stc{'n+p+r'}]; 
% $$$     set(set(histogram(afet(ind,f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     %ind = abs(nunity(afet(:,3)))<2 & sum(ufr(:,:),2)./size(ufr,2)>0.8;
% $$$     %ind = sum(fufr(ind,:),2)./size(fufr,2)>10;
% $$$     set(set(histogram(afet(ind,f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     xlim(afetLims(f,:));
% $$$ end    
% $$$ 
% $$$ nfet = copy(lpfet);
% $$$ nfet.data = [lpfet(:,1),tdrfet(:,1),rfet(:,1),lfet(:,1),hfet(:,1)];
% $$$ figure
% $$$ nind = nniz(nfet);
% $$$ for a = 1:size(nfet,2)
% $$$     for b = a+1:size(nfet,2)
% $$$         subplot2(size(nfet,2),size(nfet,2),a,b);
% $$$         hist2([nfet(nind,a),nfet(nind,b)],50,50);
% $$$         caxis([0,1000]);
% $$$     end
% $$$ end
% $$$ 
% $$$ 

figure
%nind = nniz(afet);
%nind = stcm(:,1)==1 & stcm(:,5)==5;
nind = stcm(:,1)~=1 & stcm(:,5)==5;
nind = stcm(:,1)==1 & (stcm(:,3)==3|stcm(:,4)==4);
for a = 1:size(afet,2)
    for b = a+1:size(afet,2)
        subplot2(size(afet,2),size(afet,2),a,b);
        hist2([afet(nind,a),afet(nind,b)],linspace([afetLims(a,:),30]),linspace([afetLims(b,:),30]));
        caxis([0,100]);        
    end
end

nind = stcm(:,1)==1 & (stcm(:,3)==3|stcm(:,4)==4);
[U,S,V] = svd(afet(nind,:),0);

[LU, LR, FSr, VT] = erpPCA( afet(nind,:), 5);

figure,
subplot(121);
imagesc(LR')
subplot(122);
plot(VT(:,4))


clear('hmm');
updateOM = 1;
hmm.K = 6;
hmm = hmminit(afet.data,hmm,'full');
hmm.train.cyc = 100;
hmm.obsmodel='Gauss';
hmm.train.obsupdate=ones([1,hmm.K])*updateOM;
hmm.train.init = 1;
hmm = hmmtrain(afet.data,size(afet,1),hmm);

% trcpow hrcpow tpow

save(fullfile(Trial.path.project,'analysis',['test_hmm_model.mat']),'hmm');
load(fullfile(Trial.path.project,'analysis',['test_hmm_model.mat']),'hmm');
 
diag(hmm.P)

% COMPUTE hmm states
[decode] = hmmdecode(afet.data,size(afet,1),hmm);
decode.q_star = decode.q_star';


figure
subplot(311);
ind = stcm(:,1)~=1&stcm(:,5)==5;
bar(accumarray(decode.q_star(ind),ones([sum(ind),1]),[hmm.K,1],@sum))
subplot(312);
ind = stcm(:,1)==1&stcm(:,5)==5;
bar(accumarray(decode.q_star(ind),ones([sum(ind),1]),[hmm.K,1],@sum))
subplot(313);
ind = stcm(:,1)==1&(stcm(:,2)==2|stcm(:,3)==3|stcm(:,4)==4);
bar(accumarray(decode.q_star(ind),ones([sum(ind),1]),[hmm.K,1],@sum))


afetGmean = [];
for grp = 1:hmm.K
    for aft = 1:size(afet,2)
        afetGmean(grp,aft) = mean(afet(decode.q_star==grp,aft));
    end
end
figure,
imagesc(afetGmean');

grpSymbol = '*o^s+<d';
figure();hold('on');
    ind = [Trial.stc{'w+r+n+p+s'}]; cast(ind,'TimeSeries'); resample(ind,vxy); ind = logical(ind.data);
    hist2([xcomp.data(ind),ycomp.data(ind)],xcomp.edgs,ycomp.edgs);
    grid('on'); Lines([],0,'w'); Lines(0,[],'w'); colormap('parula'); colorbar();
    for grp = 1:hmm.K
        plot(mean(xcomp.data(decode.q_star==grp)),mean(ycomp.data(decode.q_star==grp)),['m',grpSymbol(grp)]);
    end
    xlim(xcomp.edgs([1,end]));
    ylim(ycomp.edgs([1,end]));    

% ASSIGN new group order (automate later)
grpNewOrder = [5,1,4,2,6,3];


grpLabels = {'activeTheta','ImmobileTheta','LowTheta','REMTheta','ImmobileStuff','SWR'};

% REORDER grps
hmmState = decode.q_star+hmm.K;
for grp = 1:hmm.K
    hmmState((grpNewOrder(grp)+hmm.K)==hmmState) = grp;
end
thmmState = hmmState;

thmmState(thmmState==2) = 1;
thmmState(thmmState==3) = 1;
thmmState(thmmState==4) = 1;

thmmState(~xor(thmmState==1,thmmState==4)) = 3;
%figure,plot(thmmState)
thmmState(thmmState==4) = 2;

thmmStateSegs = circshift(GetSegs(thmmState,1:size(thmmState,1),301,0),151,2);
thmmStateSegs = circshift(GetSegs(thmmState,1:size(thmmState,1),301,0),151,2);

mhmmState = sq(mode(thmmStateSegs))';

fhmmState = hmmState;
fhmmState(thmmState==2&mhmmState==1) = 3;
fhmmState(thmmState==1&mhmmState==2) = 4;
fhmmState(hmmState==5&mhmmState==3) = 5;
fhmmState(hmmState==5&mhmmState==1) = 5;
fhmmState(hmmState==6&(mhmmState==1|mhmmState==2)) = 6;


figure,
subplot(311);
    hold('on');
    plot([1:size(hmmState,1)]./sampleRate,hmmState)
    plot([1:size(hmmState,1)]./sampleRate,thmmState)
    %plot([1:size(hmmState,1)]./sampleRate,sq(mode(thmmStateSegs)))
subplot(312);
    hold('on');
    plot([1:size(hmmState,1)]./sampleRate,hmmState)
subplot(313);
    hold('on');
    plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


derr = nan([size(xyz,1),1]);
derr(dc.ind(dc.ucnt>2)) = sqrt(sum(dc.esax(dc.ucnt>2,:).^2,2));
derf = nan([size(xyz,1),1]);
derf(dc.ind(dc.ucnt>2)) = dc.ecom(dc.ucnt>2,1);
figure();
subplot(511);
    imagesc(rts,rfs,log10(rys.data)');
    axis('xy');
    colormap('jet');
    caxis([4,6.5])
subplot(512);
    imagesc(lts,lfs,log10(lys.data)');
    axis('xy');
    colormap('jet');
    caxis([2,6])
subplot(513);
    hold('on');
    plot([1:size(derr,1)]./sampleRate,derf)
    plot([1:size(derr,1)]./sampleRate,derr)
    Lines([],0,'k')
    Lines([],100,'k')
subplot(514);
    hold('on');
    plot([1:size(hmmState,1)]./sampleRate,hmmState)
    plot([1:size(hmmState,1)]./sampleRate,thmmState)    
    plot([1:size(hmmState,1)]./sampleRate,rfet.data)
    %plot([1:size(hmmState,1)]./sampleRate,fhmmState)
    ylim([0.5,hmm.K+0.5]);
    grid('on');
subplot(515);
    hold('on');
    plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

nfet = copy(afet);
nfet.data = [lpfet(:,1),hfet(:,1),trat(:,1) ifet(:,1)];
nfetMeanW = median(nfet(Trial.stc{'w'},:));
nfet.data = bsxfun(@minus,nfet.data,nfetMeanW);
nfetLims = [-1,1.5; -1.5,1.5; -1.5,1.25; -1.5,1.5; -1.5,1.5];


figure
%nind = thmmState==1;
%nind = thmmState==1 & stcm(:,3)==3;
nind = thmmState==1 & stcm(:,5)==5;
for a = 1:size(nfet,2)
    for b = a+1:size(nfet,2)
        subplot2(size(nfet,2),size(nfet,2),a,b);
        hist2([nfet(nind,a),nfet(nind,b)],linspace([nfetLims(a,:),30]),linspace([nfetLims(b,:),30]));
        caxis([0,2000]);        
    end
end

clear('hmm2');
updateOM = 1;
hmm2.K = 3;
hmm2 = hmminit(nfet(thmmState==1,:),hmm2,'full');
hmm2.train.cyc = 100;
hmm2.obsmodel='Gauss';
hmm2.train.obsupdate=ones([1,hmm2.K])*updateOM;
hmm2.train.init = 1;
hmm2 = hmmtrain(nfet.data,size(nfet,1),hmm2);

diag(hmm2.P)

% COMPUTE hmm states
[decode2] = hmmdecode(nfet.data,size(nfet,1),hmm2);
decode2.q_star = decode2.q_star';

%hmmState2 = decode2.q_star+hmm2.K;
hmmState2 = decode2.q_star;
hmmState2(thmmState~=1) = 0;




figure
subplot(311);
ind = hmmState==1 & stcm(:,5)==5;
bar(accumarray(decode2.q_star(ind),ones([sum(ind),1]),[hmm.K,1],@sum))
subplot(312);
ind = hmmState==1 &(stcm(:,2)==2|stcm(:,3)==3|stcm(:,4)==4);
bar(accumarray(decode2.q_star(ind),ones([sum(ind),1]),[hmm.K,1],@sum))
subplot(313);
ind = hmmState==1 &(stcm(:,6)==6);
bar(accumarray(decode2.q_star(ind),ones([sum(ind),1]),[hmm.K,1],@sum))


derr = nan([size(xyz,1),1]);
derr(dc.ind(dc.ucnt>2)) = sqrt(sum(dc.esax(dc.ucnt>2,:).^2,2));
derf = nan([size(xyz,1),1]);
derf(dc.ind(dc.ucnt>2)) = dc.ecom(dc.ucnt>2,1);
figure();
subplot(511);
    imagesc(rts,rfs,log10(rys.data)');
    axis('xy');
    colormap('jet');
    caxis([4.25,6.25])
subplot(512);
    imagesc(lts,lfs,log10(lys.data)');
    axis('xy');
    colormap('jet');
    caxis([3,5])
subplot(513);
    hold('on');
    plot([1:size(derr,1)]./sampleRate,derf)
    plot([1:size(derr,1)]./sampleRate,derr)
    Lines([],0,'k')
    Lines([],100,'k')
subplot(514);
    hold('on');
    plot([1:size(hmmState,1)]./sampleRate,hmmState)
    plot([1:size(hmmState,1)]./sampleRate,hmmState2)    
    plot([1:size(hmmState,1)]./sampleRate,rfet.data)
    %plot([1:size(hmmState,1)]./sampleRate,fhmmState)
    ylim([0.5,hmm.K+0.5]);
    grid('on');
subplot(515);
    hold('on');
    plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

derr = zeros([size(xyz,1),1]);
derr(dc.ind(dc.ucnt>2)) = sqrt(sum(dc.esax(dc.ucnt>2,:).^2,2));

figure();
subplot(411)
    ind = thmmState==1&hmmState2==3 & (stcm(:,3)==3|stcm(:,4)==4);
    histogram(derr(ind),linspace(0,600,100),'Normalization','probability');
subplot(412)
    ind = thmmState==1&hmmState2==2 & (stcm(:,3)==3|stcm(:,4)==4);
    histogram(derr(ind),linspace(0,600,100),'Normalization','probability');
subplot(413)
    ind = thmmState==1&hmmState2==1 & (stcm(:,3)==3|stcm(:,4)==4);
    histogram(derr(ind),linspace(0,600,100),'Normalization','probability');
subplot(414)
    ind = thmmState==1&hmmState2==2 & stcm(:,5)==5;
    histogram(derr(ind),linspace(0,600,100),'Normalization','probability');
ForAllSubplots('ylim([0,0.4])');


figure,
subplot(141);
ind = stcm(:,1)==1&stcm(:,3)==3;
hist2([tdfet(ind),tdrfet(ind)],linspace(-0.05,0.2,40),linspace(-0.05,0.2,40));
subplot(142);
ind = stcm(:,1)==1&stcm(:,4)==4;
hist2([tdfet(ind),tdrfet(ind)],linspace(-0.05,0.2,40),linspace(-0.05,0.2,40));
subplot(143);
ind = stcm(:,1)==1&stcm(:,5)==5;
hist2([tdfet(ind),tdrfet(ind)],linspace(-0.05,0.2,40),linspace(-0.05,0.2,40));
subplot(144);
ind = stcm(:,1)~=1&stcm(:,5)==5&hmmState==3;
hist2([tdfet(ind),tdrfet(ind)],linspace(-0.05,0.2,40),linspace(-0.05,0.2,40));


figure,
subplot(141);
ind = stcm(:,1)==1&stcm(:,3)==3;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);
subplot(142);
ind = stcm(:,1)==1&stcm(:,4)==4;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);
subplot(143);
ind = stcm(:,1)==1&stcm(:,5)==5;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);
subplot(144);
ind = stcm(:,1)~=1&stcm(:,5)==5&hmmState==3;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);

figure,
subplot(141);
ind = stcm(:,1)==1&stcm(:,3)==3&hmmState==4;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);
subplot(142);
ind = stcm(:,1)==1&stcm(:,4)==4&hmmState==4;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);
subplot(143);
ind = stcm(:,1)==1&stcm(:,5)==5&hmmState==4;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);
subplot(144);
ind = stcm(:,1)~=1&stcm(:,5)==5&hmmState==4;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);






clear('xcomp','ycomp');
% $$$ xcomp.data = rfet(:,1);
% $$$ xcomp.edgs = linspace(3.5,6.5,30);
% $$$ xcomp.data = dfet(:,1);
% $$$ xcomp.edgs = linspace(5,8,40);
xcomp.data = lpfet(:,1);
xcomp.edgs = linspace(5.5,7.5,40);
ycomp.data = tdfet(:,1);
ycomp.edgs = linspace(-0.1,0.2,40);
% $$$ ycomp.data = dfet(:,1);
% $$$ ycomp.edgs = linspace(5,8,40);

figure,
for grp = 1:hmm.K
% $$$     mind = dc.ucnt>=1 ...
% $$$     & (  dc.stcm(:,3)==3 ...
% $$$        | dc.stcm(:,5)==5 ...
% $$$        | dc.stcm(:,2)==2 ...
% $$$        | dc.stcm(:,4)==4 ...
% $$$        | dc.stcm(:,6)==6) ...
% $$$     & decode.q_star(dc.ind)==grp;
% $$$     mind = dc.ucnt>1 ...
% $$$     & (  dc.stcm(:,3)==3 ...
% $$$        | dc.stcm(:,5)==5) ...
% $$$     & decode.q_star(dc.ind)==grp;
% $$$     mind = dc.ucnt>=1 ...
% $$$            & dc.stcm(:,1)~=1 ...
% $$$     & (  dc.stcm(:,4)==4 ...
% $$$        | dc.stcm(:,6)==6) ...
% $$$     & decode.q_star(dc.ind)==grp;
% $$$ mind = dc.ucnt>=1 & (dc.stcm(:,1)~=1&dc.stcm(:,8)==8 ) & decode.q_star(dc.ind)==grp;
% $$$ mind = dc.ucnt>1 & (dc.stcm(:,1)==1&dc.stcm(:,8)==8 ) & decode.q_star(dc.ind)==grp;


mind = dc.ucnt>1 & (dc.stcm(:,1)~=1&dc.stcm(:,7)==7 ) & decode.q_star(dc.ind)==grp;
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,derr);
climM = [0,600];climS = [0,250];
zmean = zcomp.mean;
zmean(zcomp.count<20) = nan;
zstd  = zcomp.std;
zstd(zcomp.count<20) = nan;
zcount = zcomp.count;
zcount(zcomp.count<20) = nan;
subplot2(4,hmm.K,1,grp);
    imagesc(xcomp.ctrs,ycomp.ctrs,zmean');
    Lines([],0,'w');Lines(0,[],'w');
    axis('xy'); colormap('jet'); colorbar();caxis([climM]);
subplot2(4,hmm.K,2,grp);
    imagesc(xcomp.ctrs,ycomp.ctrs,zstd');
    Lines([],0,'w'); Lines(0,[],'w');
    axis('xy'); colormap('jet'); colorbar();caxis([climS]);
subplot2(4,hmm.K,3,grp);
    imagesc(xcomp.ctrs,ycomp.ctrs,zcount');
    Lines([],0,'w'); Lines(0,[],'w');
    axis('xy'); colormap('jet'); colorbar();
subplot2(4,hmm.K,4,1);
    ind = [Trial.stc{'w+r+n+p+s'}]; cast(ind,'TimeSeries'); resample(ind,vxy); ind = logical(ind.data);
    hist2([xcomp.data(ind),ycomp.data(ind)],xcomp.edgs,ycomp.edgs);
    grid('on'); Lines([],0,'w'); Lines(0,[],'w'); colormap('jet'); colorbar();
subplot2(4,hmm.K,4,2);
    ind = [Trial.stc{'w+r+n'}]; cast(ind,'TimeSeries'); resample(ind,vxy); ind = logical(ind.data);
    hist2([xcomp.data(ind),ycomp.data(ind)],xcomp.edgs,ycomp.edgs);
    grid('on'); Lines([],0,'w'); Lines(0,[],'w'); colormap('jet'); colorbar();
subplot2(4,hmm.K,4,3);
    ind = [Trial.stc{'s'}]; cast(ind,'TimeSeries'); resample(ind,vxy); ind = logical(ind.data);
    hist2([xcomp.data(ind),ycomp.data(ind)],xcomp.edgs,ycomp.edgs);
    grid('on'); Lines([],0,'w'); Lines(0,[],'w'); colormap('jet'); colorbar();
end



stsInds = {};
for grp = 1:hmm.K
stsInds{grp} = { dc.ucnt>1 & ( dc.stcm(:,2)==2 ) & decode.q_star(dc.ind)==grp,...
            ...
            dc.ucnt>1 & ( dc.stcm(:,3)==3  ...
                        | dc.stcm(:,5)==5) ...
           & decode.q_star(dc.ind)==grp,...
            ...
            dc.ucnt>1 & ( dc.stcm(:,4)==4  ...
                        | dc.stcm(:,6)==6) ...
           & decode.q_star(dc.ind)==grp,...
           ...
           dc.ucnt>1 & ( dc.stcm(:,1)==1 & dc.stcm(:,7)==7 ) & decode.q_star(dc.ind)==grp,...
           dc.ucnt>1 & ( dc.stcm(:,1)~=1 & dc.stcm(:,7)==7 ) & decode.q_star(dc.ind)==grp,...            
           dc.ucnt>1 & ( dc.stcm(:,1)==1 & dc.stcm(:,8)==8 ) & decode.q_star(dc.ind)==grp,...
           dc.ucnt>1 & ( dc.stcm(:,1)~=1 & dc.stcm(:,8)==8 ) & decode.q_star(dc.ind)==grp ...
};
end
stsLabels = {'Rear','Loc','Pause','GroomT','GroomN','SitT','SitN'};

figure,
for sts = 1:numel(stsInds{1})
for grp = 1:hmm.K
    mind = stsInds{grp}{sts};
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
    subplot2(numel(stsInds{1}),hmm.K,sts,grp);
    set(histogram(derr(nniz(derr)),linspace(0,600,50)),'EdgeColor','none');
    if grp==1,ylabel(stsLabels{sts});end
    if sts==1,title(['HMM Grp: ',num2str(grp)]);end
end    
end

ss = 6
sind = stcm(:,ss)==ss;
%sind = true([size(stcm,1),1]);
cmat = [sum(stcm(:,1)==1 & sind & decode.q_star==1),...
        sum(stcm(:,1)==1 & sind & decode.q_star==2),...
        sum(stcm(:,1)==1 & sind & decode.q_star==3),...
        sum(stcm(:,1)==1 & sind & decode.q_star==4),...
        sum(stcm(:,1)==1 & sind & decode.q_star==5);...
        ...
        sum(stcm(:,1)~=1 & sind & decode.q_star==1),...
        sum(stcm(:,1)~=1 & sind & decode.q_star==2),...
        sum(stcm(:,1)~=1 & sind & decode.q_star==3),...
        sum(stcm(:,1)~=1 & sind & decode.q_star==4),...
        sum(stcm(:,1)~=1 & sind & decode.q_star==5)]
bsxfun(@rdivide,cmat,sum(cmat))

sum(cmat)./sum(sum(cmat))


clear('xcomp','ycomp','zcomp','wcomp');
nbins = 15;
xcomp.data = tipmmv(:,1);
xcomp.edgs = linspace(-1,1,nbins);
ycomp.data = tipmmv(:,2);
ycomp.edgs = linspace(-.4,.4,nbins);
% $$$ zcomp.data = tdpow(:,1);
% $$$ zcomp.edgs = linspace(-5,1,nbins);
%zcomp.data = lvxy(:,2);
%zcomp.edgs = linspace(-3,2,nbins);
zcomp.data = irat(:,1);
zcomp.edgs = linspace(-1,2,nbins);

cmatY = [];
cmatN = [];
for s = 1:size(stcm,2)-1,
    for g = 1:hmm.K
        cmatY(g,s) = sum(stcm(:,1)==1 & stcm(:,s+1)==s+1 & decode.q_star==g);
        cmatN(g,s) = sum(stcm(:,1)~=1 & stcm(:,s+1)==s+1 & decode.q_star==g);
    end
end





bsxfun(@rdivide,cmat,sum(cmat))





figure
hold('on');
mind = dc.ucnt>2 & dc.stcm(:,1)==1;
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
pos = cell([1,3]);
[pos{:}] = ndgrid(xcomp.ctrs,ycomp.ctrs,zcomp.ctrs);
phl = patch(isosurface(pos{:},wcomp.count,10));
phl.FaceColor = [1,0,0];
%phl.FaceAlpha = 0.3;

mind = dc.ucnt>2 & dc.stcm(:,1)~=1;
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,10));
phl.FaceColor = [0,0,1];
phl.FaceAlpha = 0.3;


figure,imagesc(wcomp.mean(:,:,15)');axis('xy');colorbar();colormap('jet');

figure
hold('on');
% $$$ mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5 | dc.stcm(:,4)==4 |dc.stcm(:,6)==6);
% $$$ mind = dc.ucnt>3 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5 );
mind = dc.ucnt>3 & (dc.stcm(:,4)==4 |dc.stcm(:,6)==6);
derr = zeros(size(xcomp.data));
%derr(dc.ind(mind)) = dc.ucnt(mind);
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
derr(dc.ind(mind)) = dc.ecom(mind,1);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
pos = cell([1,3]);
[pos{:}] = ndgrid(xcomp.ctrs,ycomp.ctrs,zcomp.ctrs);
phl = patch(isosurface(pos{:},wcomp.count,30));
phl.FaceColor = [1,0,0];

mind = dc.ucnt>3 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5 );
derr = zeros(size(xcomp.data));
%derr(dc.ind(mind)) = dc.ucnt(mind);
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
derr(dc.ind(mind)) = dc.ecom(mind,1);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
pos = cell([1,3]);
[pos{:}] = ndgrid(xcomp.ctrs,ycomp.ctrs,zcomp.ctrs);
phl = patch(isosurface(pos{:},wcomp.count,30));
phl.FaceColor = [1,0,1];

mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)==1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomvp.count,30));
phl.FaceColor = [0,0,1];
phl.FaceAlpha = 0.3;

mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,30));
phl.FaceColor = [0,1,0];
phl.FaceAlpha = 0.3;
% $$$ phl = patch(isosurface(pos{:},wcomp.count,50));
% $$$ phl.FaceColor = [0,1,0];
% $$$ phl.FaceAlpha = 0.3;



clear('xcomp','ycomp','zcomp','wcomp');
nbins = 15;
xcomp.data = lpfet(:,1);
xcomp.edgs = linspace(3.5,7.5,nbins);
% $$$ xcomp.data = mlhfet(:,1);
% $$$ xcomp.edgs = linspace(-2,2,nbins);
% $$$ xcomp.data = mlhfet(:,2);
% $$$ xcomp.edgs = linspace(-1.5,1.5,nbins);
% $$$ ycomp.data = mlhfet(:,2);
% $$$ ycomp.edgs = linspace(-1.5,1.5,nbins);
ycomp.data = rfet(:,1);
ycomp.edgs = linspace(2.5,7,nbins);
% $$$ zcomp.data = lvxy(:,2);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
% $$$ zcomp.data = log10(sqrt(sfet(:,1).^2)).*cos(pi/6)+log10(sqrt(dsfet(:,1).^2)).*sin(pi/6);
% $$$ zcomp.edgs = linspace(-5,2,nbins);
% $$$ zcomp.data = rfet(:,1);
% $$$ zcomp.edgs = linspace(2.5,6,nbins);
zcomp.data = tdfet(:,1);
zcomp.edgs = linspace(-0.1,0.2,nbins);
% $$$ zcomp.data = irat(:,1);
% $$$ zcomp.edgs = linspace(-3,3,nbins);
% $$$ zcomp.data = tfet(:,1);
% $$$ zcomp.edgs = linspace(-5.5,8,nbins);
% $$$ zcomp.data = frat(:,1);
% $$$ zcomp.edgs = linspace(-3,3,nbins);
% $$$ zcomp.data = vrat(:,1);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
% $$$ zcomp.data = trat(:,1);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
% $$$ zcomp.data = nrat(:,1);
% $$$ zcomp.edgs = linspace(-2,2,nbins);
% $$$ zcomp.data = prat(:,1);
% $$$ zcomp.edgs = linspace(0.5,6,nbins);

figure();
hold('on');
ind = stc{'s'};
plot3(sdfet(ind),rfet(ind,1),tdfet(ind,1),'.')
ind = stc{'w+n+r'};
plot3(sdfet(ind),rfet(ind,1),tdfet(ind,1),'.')

figure();
hold('on');
ind = stc{'s'};
plot(rfet(ind),mlhfet(ind,2),'.')
ind = stc{'w+n+r'};
plot(rfet(ind),mlhfet(ind,2),'.')


isothresh = [500];
isothresh = [300];
isothresh = [250];
isothresh = [200];
isothresh = [100];
isothresh = [50];
isothresh = [20];

figure();
hold('on');
mvec = [];

mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
[pos{:}] = ndgrid(xcomp.ctrs,ycomp.ctrs,zcomp.ctrs);
%phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.std));
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
phl.FaceAlpha = 0.5;
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
caxis([0,400])
%caxis([0,16])
mvec(end+1,:) = [mean(xcomp.data(nniz(derr))),...
                 mean(ycomp.data(nniz(derr))),...
                 mean(zcomp.data(nniz(derr)))];
scatter3(mvec(end,1),...
         mvec(end,2),...
         mvec(end,3),...
         100,...
         'k',...
         'Filled');



mind = dc.ucnt>1 ...
       & dc.stcm(:,1)==1 ...
       & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2 );
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
phl.FaceAlpha = 0.5;
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
mvec(end+1,:) = [mean(xcomp.data(nniz(derr))),...
                 mean(ycomp.data(nniz(derr))),...
                 mean(zcomp.data(nniz(derr)))];
scatter3(mvec(end,1),...
         mvec(end,2),...
         mvec(end,3),...
         100,...
         'k',...
         'Filled');

bvec = diff(mvec(1:2,:));
bvec = bvec./sqrt(sum(bvec.^2,2));

tfet = sum(bsxfun(@times,bsxfun(@minus,[xcomp.data,ycomp.data,zcomp.data],mvec(2,:)),bvec),2);
t2fet = tfet;

figure,
subplot(211);
    mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2 );
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
    hist2([tfet(nniz(derr)),t2fet(nniz(derr))],linspace(-2,3,30),linspace(-1,0.6,30));
    Lines(-0.5,[],'k');Lines([],-0.25,'k');
    caxis([0,1000])
subplot(212);
    mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
    hist2([tfet(nniz(derr)),t2fet(nniz(derr))],linspace(-2,3,30),linspace(-1,0.6,30));
    Lines(-0.5,[],'k');Lines([],-0.25,'k');
    caxis([0,1000])

mind = dc.ucnt>1 ...
       & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2|dc.stcm(:,8)==8);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));

cmat = [sum(stcm(:,1)==1 & (tfet(:)>-0.3&t2fet(:)>-0.75)),...
        sum(stcm(:,1)~=1 & (tfet(:)>-0.3&t2fet(:)>-0.75));...
        sum(stcm(:,1)==1 & (tfet(:)<-0.3|t2fet(:)<-0.75)),...
        sum(stcm(:,1)~=1 & (tfet(:)<-0.3|t2fet(:)<-0.75))]
bsxfun(@rdivide,cmat,sum(cmat))
    
figure();
    mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2|dc.stcm(:,8)==8);
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
    hist2([tfet(nniz(derr)),t2fet(nniz(derr))],linspace(-2,3,30),linspace(-1,0.6,30));

    

figure();
mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2 );
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
hhdl = histogram(tfet(nniz(derr)),linspace(-2.5,4,100));
hhdl.EdgeColor = 'none';
hold('on');
mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
hhdl = histogram(tfet(nniz(derr)),linspace(-2.5,4,100));
hhdl.EdgeColor = 'none';
hhdl.FaceAlpha = 0.5;

figure();
subplot(211);
plot([1:size(tfet,1)]./sampleRate,tfet);
subplot(212);
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


clear('xcomp','ycomp','zcomp','wcomp');
nbins = 15;

sdfet = sfet.copy();
sdfet.data = log10(sqrt(sfet(:,1).^2)).*cos(pi/6)+log10(sqrt(dsfet(:,1).^2)).*sin(pi/6);

afet = [tdfet(:,1),lpfet(:,1)];
%afet = [mlhfet(:,1),mlhfet(:,2),tdfet(:,1),irat(:,1)];
%afet = [tdfet(:,1),sdfet(:,1),irat(:,1)];


mvec = [];
mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
mvec(end+1,:)= mean(afet(nniz(derr),:));
mind = dc.ucnt>1 ...
       & dc.stcm(:,1)==1 ...
       & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2 );
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
mvec(end+1,:)= mean(afet(nniz(derr),:));


wfet = bsxfun(@minus,afet,mvec(2,:));
ind = stcm(:,1)==1 & (stcm(:,2)==2 | stcm(:,3)==3 | stcm(:,4)==4);
wsigma = cov(wfet(ind,:));
wscore = copy(xyz);
wscore.data = exp(-1/2*(multiprod(wfet,multiprod((wsigma.*40).^-1,wfet,[1,2],2),2,2)))...
          ./sqrt((2*pi).^size(wfet,2).*det(wsigma));

ssfet = bsxfun(@minus,afet,mvec(1,:));
ind = stcm(:,1)~=1 & (stcm(:,5)==5);
ssigma = cov(wfet(ind,:));
sscore = copy(xyz);
sscore.data = exp(-1/2*(multiprod(ssfet,multiprod((ssigma.*40).^-1,ssfet,[1,2],2),2,2)))...
          ./sqrt((2*pi).^size(ssfet,2).*det(ssigma));


figure,
hold('on');
ind = stcm(:,1)~=1 & (stcm(:,5)==5);
plot(log10(sscore(ind)),log10(wscore(ind)),'.')
ind = stcm(:,1)==1 & (stcm(:,2)==2 | stcm(:,3)==3 | stcm(:,4)==4);
plot(log10(sscore(ind)),log10(wscore(ind)),'.')

figure,
subplot(211);
    ind = stcm(:,1)~=1 & (stcm(:,5)==5);
    hist2([log10(sscore(ind)),log10(wscore(ind))],linspace(-20,10,30),linspace(-100,150,30));
    caxis([0,1000]);
    Lines([],25,'k');Lines(0,[],'k');    
subplot(212);
    ind = stcm(:,1)==1 & (stcm(:,2)==2 | stcm(:,3)==3 | stcm(:,4)==4);
    hist2([log10(sscore(ind)),log10(wscore(ind))],linspace(-20,10,30),linspace(-100,150,30));
    caxis([0,1000]);
    Lines([],25,'k');Lines(0,[],'k');

arot = 0.05;
nsscore = log10(sscore(:,1)).*cos(arot)+log10(wscore(:,1)).*sin(arot);   
nwscore = log10(sscore(:,1)).*-sin(arot)+log10(wscore(:,1)).*cos(arot);   


figure,
hold('on');
ind = stcm(:,1)~=1 & (stcm(:,5)==5);
plot(nsscore(ind),nwscore(ind),'.')
ind = stcm(:,1)==1 & (stcm(:,5)==5);
plot(nsscore(ind),nwscore(ind),'.g')
ind = stcm(:,1)==1 & (stcm(:,2)==2 | stcm(:,3)==3 | stcm(:,4)==4);
plot(nsscore(ind),nwscore(ind),'.r')
ind = stcm(:,1)~=1 & (stcm(:,2)==2 | stcm(:,3)==3 | stcm(:,4)==4);
plot(nsscore(ind),nwscore(ind),'.m')

wthresh = -50;
sthresh = 1;
nind = nniz(xyz);
cmat = [sum(stcm(nind,1)==1&nsscore(nind)<sthresh&nwscore(nind)>wthresh),...
        sum(stcm(nind,1)~=1&nsscore(nind)<sthresh&nwscore(nind)>wthresh);...
        sum(stcm(nind,1)==1&(nsscore(nind)>sthresh|nwscore(nind)<wthresh)),...
        sum(stcm(nind,1)~=1&(nsscore(nind)>sthresh|nwscore(nind)<wthresh))]
bsxfun(@rdivide,cmat,sum(cmat))



    
wthresh = 1e4;
nind = nniz(xyz);
cmat = [sum(stcm(nind,1)==1&wscore(nind)<wthresh),sum(stcm(nind,1)~=1&wscore(nind)<wthresh);...
        sum(stcm(nind,1)==1&wscore(nind)>wthresh),sum(stcm(nind,1)~=1&wscore(nind)>wthresh)]
bsxfun(@rdivide,cmat,sum(cmat))


derr = zeros(size(xcomp.data));
derr(dc.ind) = sqrt(sum(dc.ecom(:,:).^2,2));

figure();
subplot(211);
hold('on');
plot([1:size(wscore,1)]./sampleRate,double(nsscore(:)<sthresh&nwscore(:)>wthresh).*600);
plot([1:size(wscore,1)]./sampleRate,derr,'r');
%plot([1:size(wscore,1)]./sampleRate,log10(wscore.data));
subplot(212);
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


nsscore(:)<sthresh&nwscore(:)>wthresh

figure();
mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2 );
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
hhdl = histogram(log10(wscore(nniz(derr))),linspace(-200,200,100));
hhdl.EdgeColor = 'none';
hold('on');
mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
hhdl = histogram(log10(wscore(nniz(derr))),linspace(-200,200,100));
hhdl.EdgeColor = 'none';
hhdl.FaceAlpha = 0.5;



mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)==1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
phl.FaceAlpha = 0.5;

mind = dc.ucnt>1 & (dc.stcm(:,6)==6);
derr = zeros(size(xcomp.data));
%derr(dc.ind(mind)) = dc.ecom(mind,1);
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
%phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.std));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
%phl.FaceAlpha = 0.5;


mind = dc.ucnt>2 & (dc.stcm(:,4)==4);
derr = zeros(size(xcomp.data));
%derr(dc.ind(mind)) = dc.ecom(mind,1);
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
%phl.FaceAlpha = 0.5;

mind = dc.ucnt>2 & (dc.stcm(:,2)==2);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
%phl.FaceAlpha = 0.5;

%mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5 );
mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2 );
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
%phl.FaceAlpha = 0.5;



%norm = 'yprob';
norm = '';
figure
subplot(311);
ind = stcm(:,1)==1|stcm(:,5)==5;
hist2([log10(sfet(ind,1).^2),log10(dsfet(ind,1).^2)],linspace(-5,1.5,30),linspace(-5,1.5,30),norm);
caxis([0,200])
%caxis([0,50])
Lines([],0,'k');
Lines([],-1,'k');
Lines(0,[],'k');
subplot(312);
ind = stcm(:,1)~=1&stcm(:,5)==5;
%ind(ind) = rand([sum(ind),1])>0.86;
hist2([log10(sfet(ind,1).^2),log10(dsfet(ind,1).^2)],linspace(-5,1.5,30),linspace(-5,1.5,30),norm);
caxis([0,200])
%caxis([0,50])
Lines([],0,'k');
Lines([],-1,'k');
Lines(0,[],'k');
subplot(313);
ind = stcm(:,1)==1&stcm(:,2)==2;
hist2([log10(sfet(ind,1).^2),log10(dsfet(ind,1).^2)],linspace(-5,1.5,30),linspace(-5,1.5,30),norm);
caxis([0,200])
Lines([],0,'k');
Lines([],-1,'k');
Lines(0,[],'k');


figure
mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5 );
derr = zeros(size(sfet.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ecom(mind,1);
%hist2([log10(sfet(nniz(derr),1).^2),derr(nniz(derr))],linspace(-5,1.5,30),linspace(-150,150,30),'xprob');
hist2([log10(sfet(nniz(derr),1).^2),derr(nniz(derr))],linspace(-5,1.5,30),linspace(0,200,30),'xprob');
corr(log10(sfet(nniz(derr),1).^2),derr(nniz(derr)))


xlim(xcomp.edgs([1,end]));ylim(ycomp.edgs([1,end]));zlim(zcomp.edgs([1,end]));

xlim([-1.4,1.4]);ylim([-0.5,0.5]);zlim([-1,2]);

mind = dc.ucnt>2 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
mind = dc.ucnt>2 & (dc.stcm(:,8)==8&dc.stcm(:,1)==1);
mind = dc.ucnt>2 & (dc.stcm(:,3)==3&dc.stcm(:,1)==1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.esax(mind,:).^2,2));
figure,plot3(trat(dc.ind(mind)),irat(dc.ind(mind)),derr(dc.ind(mind)),'.r')
xlim([-2,1]);ylim([-1,1]);zlim([0,800]);


mind = dc.ucnt>2 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
mind = dc.ucnt>2 & (dc.stcm(:,8)==8&dc.stcm(:,1)==1);
mind = dc.ucnt>2 & (dc.stcm(:,3)==3&dc.stcm(:,1)==1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.esax(mind,:).^2,2));
figure,plot3(log10(abs(out(dc.ind(mind),1))),log10(abs(out(dc.ind(mind),2))),derr(dc.ind(mind)),'.r')
xlim([-2,1]);ylim([-1,1]);zlim([0,800]);

view(-169.499998031745,       2.40000023759164);

 
 
 
 



% $$$ phl = patch(isosurface(pos{:},wcomp.count,50));
% $$$ phl.FaceColor = [0,1,0];
% $$$ phl.FaceAlpha = 0.3;




%mind = dc.ucnt>2 & (dc.stcm(:,4)==4 |dc.stcm(:,6)==6);
%mind = dc.ucnt>3 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5);
mind = dc.ucnt>3 & (dc.stcm(:,1)~=1 & dc.stcm(:,8)==8);
derr = zeros(size(xcomp.data));
%derr(dc.ind(mind)) = dc.ecom(mind,1);
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);

figure,
wmean = wcomp.mean;
wmean(wcomp.count(:)<20) = nan;
for s = 1:nbins-1
    subplot(2,7,s);
    set(pcolor(xcomp.ctrs,ycomp.ctrs,wmean(:,:,s)'),'EdgeColor','none');
    axis('xy')
    colormap('jet');
    caxis([0,400]);
    %caxis([-40,80]);
end



swrper = [Trial.stc{'R'}];
swrper.cast('TimeSeries');
swrper.resample(xyz);

wper = [Trial.stc{'w'}];
wper.cast('TimeSeries');
wper.resample(xyz);


figure();
hold('on');
plot(mpv,mean(ufr(logical(swrper.data),:))')
plot(mpv,mean(ufr(logical(wper.data),:))')
xlim([0,2*pi]);



clear('xcomp','ycomp','zcomp','wcomp');
nbins = 15;
xcomp.data = mlhfet(:,1);
xcomp.edgs = linspace(-2,2,nbins);
ycomp.data = mlhfet(:,2);
ycomp.edgs = linspace(-1.5,1.5,nbins);
% $$$ zcomp.data = lvxy(:,2);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
zcomp.data = irat(:,1);
zcomp.edgs = linspace(-3,3,nbins);
% $$$ zcomp.data = frat(:,1);
% $$$ zcomp.edgs = linspace(-3,3,nbins);
% $$$ zcomp.data = vrat(:,1);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
% $$$ zcomp.data = trat(:,1);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
%fet = dc.ecom(:,1);
fet = sqrt(sum(dc.ecom(:,:).^2,2)); climM = [0,400]; climS = [0,200];
%fet = dc.ucnt; climM = [0,16]; climS = [0,8];


isothresh = [100];
isothresh = [50];
figure();
hold('on');
sp = gobjects([0,0]);

sp(end+1) = subplot(221);
    mind = dc.ucnt>2 & (dc.stcm(:,4)==4|dc.stcm(:,6)==6)&dc.stcm(:,1)==1;
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = fet(mind,1);
    [xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
    phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
    phl.FaceColor = 'interp';
    phl.EdgeColor = 'none';
    colormap('jet');
    sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
    caxis(climM)
    colorbar();    
    xlim(xcomp.edgs([1,end]));ylim(ycomp.edgs([1,end]));zlim(zcomp.edgs([1,end]));    
sp(end+1) = subplot(223);
    mind = dc.ucnt>2 & (dc.stcm(:,4)==4|dc.stcm(:,6)==6)&dc.stcm(:,1)==1;
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = fet(mind,1);    
    [xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
    phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.std));
    phl.FaceColor = 'interp';
    phl.EdgeColor = 'none';
    colormap('jet');
    sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
    caxis(climS)
    colorbar();    
sp(end+1) = subplot(222);
    mind = dc.ucnt>2 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = fet(mind,1);        
    [xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
    phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
    phl.FaceColor = 'interp';
    phl.EdgeColor = 'none';
    colormap('jet');
    sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
    caxis(climM)
    colorbar();
sp(end+1) = subplot(224);
    mind = dc.ucnt>2 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = fet(mind,1);    
    [xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
    phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.std));
    phl.FaceColor = 'interp';
    phl.EdgeColor = 'none';
    colormap('jet');
    sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
    caxis(climS)
    colorbar();    
    %linkprop(sp,{'XLim','YLim','ZLim','CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle'});
linkprop(sp,{'XLim','YLim','ZLim','View'});
    %linkprop(sp,{'View'});    
