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
trialId = 7;


Trial = Trials{trialId};
unitSubset = units{trialId};
meta = sessionList(trialId);


Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',meta.subject.channelGroup.theta);

% LOAD 
flfp = [diff(get(Trial.load('lfp',[33,40]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2),...
        diff(get(Trial.load('lfp',[41,48]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2),...
        diff(get(Trial.load('lfp',[49,56]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2),...
        diff(get(Trial.load('lfp',[57,64]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2)];
rlfp = Trial.load('lfp',[61]);

% LOAD rat positions
xyz = preproc_xyz(Trial,'trb',sampleRate);
fxyz = filter(xyz.copy(),'ButFilter',3,14,'low');
vxy = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
vxyz = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2,3]);
lvxy = copy(vxy);
lvxy.data(lvxy.data<=0.0001) = 0.0001;
lvxy.data = log10(lvxy.data);

% LOAD behavioral state labels
stcm = stc2mat(Trial.stc,xyz,{'theta','rear','loc','pause','sit','groom'});

unitsInt = Trial.spk.get_unit_set(Trial,'interneurons');
int = Trial.load('spk', Trial.lfp.sampleRate, '', unitsInt, '');

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

ufr = Trial.load('ufr', lvxy, [], unitsInt, 0.12, 'boxcar');
fufr = Trial.load('ufr', lvxy, [], unitsInt, 0.5, 'gauss');
% $$$ ufrp = Trial.load('ufr', lvxy, [], unitsSubset, 0.12, 'boxcar');

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
% $$$ figure,
% $$$ udcor = [];
% $$$ for u = 1:numel(unitsInt),
% $$$     subplot(2,7,u);
% $$$ % $$$     sufr = ufr(dc.ind(((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2)),u);
% $$$ % $$$     dcom = sqrt(sum(dc.ecom((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2,:).^2,2));
% $$$     sufr = ufr(dc.ind((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1)),u);
% $$$     dcom = sqrt(sum(dc.ecom((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1),:).^2,2));
% $$$     nind = nniz(sufr) & nniz(dcom) & sufr>0.0001;
% $$$     plot(sufr(nind),dcom(nind),'.');
% $$$     udcor(u) = corr(sufr(nind),dcom(nind));
% $$$ end
% $$$ 
% $$$ 
% $$$ figure,
% $$$ uvcor = [];
% $$$ for u = 1:numel(unitsInt),
% $$$     subplot(2,7,u);
% $$$     sufr = ufr(dc.ind(((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2)),u);
% $$$     dcom = sqrt(sum(dc.ecom((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2,:).^2,2));
% $$$ % $$$     sufr = ufr(dc.ind((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1)),u);
% $$$ % $$$     dcom = lvxy(dc.ind((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1)),2);
% $$$     nind = nniz(sufr) & nniz(dcom) & sufr>0.0001;
% $$$     plot(sufr(nind),dcom(nind),'.');
% $$$     uvcor(u) = corr(sufr(nind),dcom(nind));
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ gufr = Trial.load('ufr',xyz,int,unitsInt,0.125,'gauss');
% $$$ %ind = (stcm(:,2)==2|stcm(:,3)==3|stcm(:,4)==4|stcm(:,5)==5);
% $$$ ind = (stcm(:,5)==5);
% $$$ %ind = (stcm(:,4)==4|stcm(:,5)==5) & stcm(:,1)~=1;
% $$$ suI = gufr(ind,:);
% $$$ svI = lvxy(ind,:);
% $$$ nClust = 2;
% $$$ [idx,C,SUMD,D] = kmeans(suI,nClust,'Distance','correlation');
% $$$ 
% $$$ kCov = zeros(size(gufr,2),size(gufr,2),nClust);
% $$$ for vind = 1:nClust,
% $$$     kCov(:,:,vind) = (bsxfun(@minus,suI(idx==vind,:),median(suI(idx==vind,:)))' ...
% $$$                        *bsxfun(@minus,suI(idx==vind,:),median(suI(idx==vind,:)))) ...
% $$$                       ./sum(idx==vind);
% $$$ end
% $$$ 
% $$$ [W,H] = nnmf(suI,2);
% $$$ 
% $$$ p1 = multiprod(gufr(:,:),H(1,:)',2,1);
% $$$ p2 = multiprod(gufr(:,:),H(2,:)',2,1);
% $$$ figure,plot(p1);hold('on');plot(p2);
% $$$ figure,plot(log(p1./p2))
% $$$ 
% $$$ 
% $$$ velBins = linspace(-2.5,1.8,20);
% $$$ velInds = discretize(lvxy(:,2),velBins);
% $$$ figure,
% $$$ for vind = 1:nClust,
% $$$     subplot2(nClust,2,vind,1);
% $$$         imagesc(kCov(:,:,vind));
% $$$     subplot2(nClust,2,vind,2);
% $$$     histogram(svI(idx==vind,2),velBins);
% $$$ end
% $$$ 
% $$$ out = [];
% $$$ for i =  1:nClust
% $$$     out(:,i) = -.5*log(det(kCov(:,:,i)))...
% $$$         -.5*dot((bsxfun(@minus,gufr(:,:),median(suI(idx==vind,:)))/kCov(:,:,i))',...
% $$$                  bsxfun(@minus,gufr(:,:),median(suI(idx==vind,:)))')';
% $$$ end
% $$$ figure,plot(out)
% $$$ 
% $$$ 
% $$$ figure();
% $$$ ind = (stcm(:,5)==5);
% $$$ plot(sum(gufr(ind,:),2)./14,lvxy(ind,2),'.')
% $$$ hold('on');
% $$$ ind = (stcm(:,3)==3|stcm(:,4)==4);
% $$$ plot(sum(gufr(ind,:),2)./14,lvxy(ind,2),'.r')
% $$$ 
% $$$ udThresh = 0.1;
% $$$ udThresh = 0.0;
% $$$ irat = copy(vxy);
% $$$ irat.data = log10(  RectFilter(sum(ufr(:,udcor>0&abs(udcor)>udThresh),2)./sum(udcor>0&abs(udcor)>udThresh),9,3) ...
% $$$                   ./RectFilter(sum(ufr(:,udcor<0&abs(udcor)>udThresh),2)./sum(udcor<0&abs(udcor)>udThresh),9,3));
% $$$ 
% $$$ uvThresh = 0.0;
% $$$ vrat = copy(vxy);
% $$$ vrat.data = log10(  RectFilter(sum(ufr(:,udcor>0&abs(udcor)>uvThresh),2)./sum(udcor>0&abs(udcor)>uvThresh),9,3) ...
% $$$                   ./RectFilter(sum(ufr(:,udcor<0&abs(udcor)>uvThresh),2)./sum(udcor<0&abs(udcor)>uvThresh),9,3));

% DECODING vars -----------------------------------------------------------------------------

dc = accumulate_decoding_vars(Trial,                               ...
                              unitSubset,                          ...
                              meta.subject.channelGroup.theta,     ...
                              meta.subject.correction.thetaPhase,  ...
                              meta.subject.correction.headYaw,     ...
                              meta.subject.correction.headBody);



%fwin = 128;
rcChan = 4;
%rcChan = 1;
fwin = gausswin(64);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = RectFilter(flfp(:,rcChan),3,3);
elfp.filter('ButFilter', 4, [75,125], 'bandpass');
hfet = copy(elfp);
hfet.data = log10(conv(hfet.data.^2,fwin,'same'));
hfet.filter('ButFilter', 4, [1], 'low');
hfet.resample(xyz);


fwin = gausswin(256);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = rlfp(:,1);
elfp.filter('ButFilter', 4, [125,250], 'bandpass');
rfet = copy(elfp);
rfet.data = log10(conv(rfet.data.^2,fwin,'same'));
rfet.resample(xyz);

fwin = gausswin(2^10);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = rlfp(:,1);
elfp.filter('ButFilter', 4, [20], 'low');
lpfet = copy(elfp);
lpfet.data = log10(conv(lpfet.data.^2,fwin,'same'));
lpfet.resample(xyz);

% $$$ drfet = rfet.copy();
% $$$ drfet.data = diff([0;rfet.data]);
% $$$ ind = [stc{'s'}];
% $$$ ind = [stc{'w+n+p+r'}];


fwin = gausswin(2^10);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = RectFilter(rlfp(:,1),3,3);
elfp.filter('ButFilter', 4, [5,12], 'bandpass');
dfet = copy(elfp);
dfet.data = log10(conv(dfet.data.^2,fwin,'same'));
elfp = copy(lfp);
elfp.data = RectFilter(rlfp(:,1),3,3);
elfp.filter('ButFilter', 4, [4], 'low');
tdfet = copy(elfp);
tdfet.data = log10(conv(tdfet.data.^2,fwin,'same'));
tdfet.data = log10(dfet.data)-log10(tdfet.data);
tdfet.resample(xyz);



fwin = gausswin(2^10);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = flfp(:,4);
elfp.filter('ButFilter', 4, [5,12], 'bandpass');
drfet = copy(elfp);
drfet.data = log10(conv(drfet.data.^2,fwin,'same'));
elfp = copy(lfp);
elfp.data = flfp(:,4);
elfp.filter('ButFilter', 4, [4], 'low');
tdrfet = copy(elfp);
tdrfet.data = log10(conv(tdrfet.data.^2,fwin,'same'));
tdrfet.data = log10(drfet.data)-log10(tdrfet.data);
%rfet.filter('ButFilter', 4, [1], 'low');
tdrfet.resample(xyz);
%figure;plot(unity(ipm.data));hold('on');plot(unity(hfet.data));

fwin = gausswin(2^11);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = rlfp(:,1);
elfp.filter('ButFilter', 4, [5,12], 'bandpass');
tfet = copy(elfp);
tfet.data = log10(conv(tfet.data.^2,fwin,'same'));
tfet.resample(xyz);


fwin = gausswin(2^11);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = flfp(:,4);
elfp.filter('ButFilter', 4, [5,12], 'bandpass');
trfet = copy(elfp);
trfet.data = log10(conv(trfet.data.^2,fwin,'same'));
trfet.resample(xyz);


fwin = gausswin(2^10);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = flfp(:,4);
elfp.filter('ButFilter', 4, [5,12], 'bandpass');
dfet = copy(elfp);
dfet.data = log10(conv(dfet.data.^2,fwin,'same'));
dfet.resample(xyz);

fwin = gausswin(128);
elfp = copy(lfp);
elfp.data = flfp(:,rcChan);
elfp.filter('ButFilter', 4, [1,20], 'bandpass');
lfet = copy(elfp);
lfet.data = log10(conv(lfet.data.^2,fwin,'same'));
lfet.resample(xyz);
lfet.filter('ButFilter', 4, [1], 'low');


% $$$ sfet = copy(sfslfp);
% $$$ sfet.resample(xyz);
% $$$ dsfet = sfet.copy;
% $$$ dsfet.data = nunity(circshift(dsfet.data,1)-dsfet.data);



figure();
subplot(211);
    hold('on');
    plot([1:size(hfet)]./lfet.sampleRate,hfet.data)
    plot([1:size(lfet)]./hfet.sampleRate,lfet.data)
subplot(212);
    hold('on');
    plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


ind = [Trial.stc{'s+w+n+r+p'}];
ind = [Trial.stc{'w+n+p'}];
[Uf,Sf,Vf] = svd([lfet(ind)-mean(lfet(ind)),hfet(ind)-mean(hfet(ind))],0);
Vf
% $$$ Vf = [0.771753449181619,        -0.635921861297655;
% $$$      0.635921861297655,         0.771753449181619];
Vf = [-0.83534705515978,         0.549722927879022; ...
        -0.549722927879022,         -0.83534705515978];
Vf = [0.916079244726519,        -0.400997278521052;...
         0.400997278521052,         0.916079244726519];

aind = [Trial.stc{'w+p'}];
aind.cast('TimeSeries');
aind.resample(xyz);
Vv = [];
for v = 1:19,
    ind = v==velInds & logical(aind.data);
    ml(v) = mean(lfet(ind));
    hl(v) = mean(hfet(ind));
    [Uv,Sv,Vv(:,:,v)] = svd([lfet(ind)-ml(v),hfet(ind)-hl(v)],0);
end    
figure,plot(abs(sq(Vv(:,1,:)))')
figure,plot(ml,hl,'-+');

mlfet = lfet.copy();  mlfet.data = mlfet.data - mean(mlfet(ind));
mhfet = hfet.copy();  mhfet.data = mhfet.data - mean(mhfet(ind));

mlhfet = lfet.copy();
mlhfet.data = multiprod([mlfet.data,mhfet.data],Vf,2,[1,2]);

inds = {[Trial.stc{'s+w+r+n+p'}],[Trial.stc{'s'}],[Trial.stc{'w'}],[Trial.stc{'p'}]};
figure,
normf = '';      clim = [0,300];
for s = 1:numel(inds),
    subplot(2,2,s);
    ind = inds{s};
    hist2(multiprod([mlfet(ind),mhfet(ind)],Vf,2,[1,2]),linspace(-2,2,50),linspace(-1,1,50),normf);
    Lines([],0,'m');
    Lines(0,[],'m');
    caxis(clim);
end


% $$$ 
% $$$ % $$$ elfp = copy(lfp);
% $$$ % $$$ elfp.data = diff(glfp(1:1e6,[1,3]),1,2);
% $$$ elfp = copy(lfp);
% $$$ elfp.data = RectFilter(glfp,3,3);
% $$$ specArgsTheta = struct('nFFT',2^11,...
% $$$                   'Fs',  elfp.sampleRate,...
% $$$                   'WinLength',2^10,...
% $$$                   'nOverlap',2^10*0.875,...
% $$$                   'NW',3,...
% $$$                   'Detrend',[],...
% $$$                   'nTapers',[],...
% $$$                   'FreqRange',[1,250]);
% $$$ [lys,lfs,lts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);
% $$$ 
% $$$ 
% $$$ % $$$ elfp = copy(lfp);
% $$$ % $$$ elfp.data = diff(glfp(1:1e6,[1,3]),1,2);
% $$$ elfp = copy(lfp);
% $$$ specArgsTheta = struct('nFFT',2^11,...
% $$$                   'Fs',  elfp.sampleRate,...
% $$$                   'WinLength',2^10,...
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





figure
subplot(411);
imagesc(mts,mfs,log10(mean(mys.data,3))');
axis('xy');
colormap('jet');
subplot(412);
imagesc(mts,mfs,log10(std(mys.data,[],3))');
%imagesc(tts,tfs,log10(tys.data)');
axis('xy');
colormap('jet');
subplot(413);
imagesc(mts,mfs,log10(mean(mys.data,3))'.^2./log10(std(mys.data,[],3))');
%imagesc(tts,tfs,log10(tys.data)');
axis('xy');
colormap('jet');
subplot(414);
hold('on');
plotSTC(Trial.stc,1,'text',states,'krggbb');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


figure
subplot(511);
imagesc(mts,mfs,nunity(log10(mys(:,:,1)))');
hold('on');
plot((1:size(elfp,1))./elfp.sampleRate,nunity(elfp.data(:,3))*10+30,'w');
axis('xy');
colormap('jet');
subplot(512);
imagesc(mts,mfs,nunity(log10(mys(:,:,2)))');
axis('xy');
colormap('jet');
subplot(513);
imagesc(mts,mfs,nunity(log10(mys(:,:,3)))');
axis('xy');
colormap('jet');
subplot(514);
imagesc(mts,mfs,nunity(log10(mys(:,:,4)))');
axis('xy');
colormap('jet');
subplot(515);
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


figure
subplot(511);
imagesc(mts,mfs,(log10(mys(:,:,2)))');
hold('on');
plot((1:size(elfp,1))./elfp.sampleRate,nunity(elfp.data(:,3))*10+30,'w');
axis('xy');
colormap('jet');
grid('on')
caxis([3.75,5.25])
subplot(512);
imagesc(mts,mfs,(log10(mys(:,:,3)))');
%imagesc(mts,mfs,(log10(mys(:,:,2)))');
axis('xy');
colormap('jet');
grid('on')
caxis([3.75,5.25])
subplot(513);
imagesc(mts,mfs,(log10(mys(:,:,4)))');
axis('xy');
colormap('jet');
grid('on')
caxis([4,5.25])
subplot(514);
hold('on');
imagesc(mts,mfs,bsxfun(@rdivide,RectFilter(log10(mys(:,:,4)),3,1),sum(RectFilter(log10(mys(:,:,4)),3,1),2))');
axis('xy');
colormap('jet');
grid('on');
caxis([0.019,0.022]);
plot((1:size(elfp,1))./elfp.sampleRate,nunity(elfp.data(:,3))*10+30,'w');
subplot(515);
hold('on');
plotSTC(Trial.stc,1,'text',states,'krggbb');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
plot([1:size(vxy)]./vxy.sampleRate,(tdRatio(:,1)+1)*5,'r');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

figure
subplot(211);
%plot(mts,mean((log10(mys(:,:,3))),2)');
hold('on');
plot(mts,mean((log10(mys(:,mfs<200&mfs>40,3))),2));
plot(mts,mean(RectFilter(log10(mys(:,mfs<200&mfs>40,3)),11,3),2));
%plot(mts,mean(RectFilter(log10(mys(:,mfs<100,4)),7,3),2));
plot(mts,mean(log10(mys(:,mfs<20,3)),2),'r');
plot(mts,mean(RectFilter(log10(mys(:,mfs<20,3)),11,3),2),'k');
plot(mts,mean(RectFilter(log10(mys(:,mfs<20,3)),21,3),2),'c');
Lines([],4,'k');
Lines([],4.5,'k');
Lines([],5,'k');
subplot(212);
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


fmys = copy(mys);
fmys.data =  RectFilter(log10(mys(:,:,:)),21,3);

rcInd = 1;
tpm = copy(mys);
tpm.data = mean(fmys(:,mfs<20,rcInd),2);
tpm.resample(vxy);

mpm = copy(mys);
mpm.data = mean(fmys(:,mfs>140 & mfs<200,rcInd),2);
mpm.resample(vxy);

ipm = copy(mys);
ipm.data = mean(fmys(:,mfs>40 & mfs<200,rcInd),2);
ipm.resample(vxy);

gpm = copy(mys);
gpm.data = mean(fmys(:,mfs>50 & mfs<100,rcInd),2);
gpm.resample(vxy);

hpm = copy(mys);
hpm.data = mean(fmys(:,mfs>200,rcInd),2);
hpm.resample(vxy);

tdpow = copy(rys);
tdpow.resample(vxy);
tdpow.data = log(abs(log(mean(tdpow(:,rfs>5 & rfs<12),2))-log(mean(tdpow(:,rfs<5 | (rfs>12 & rfs<15)),2))));

tpow = copy(rys);
tpow.resample(vxy);
tpow.data = log(mean(tpow(:,rfs>5 & rfs<12),2));


% ERROR in stc s+t is not correct
figure
normf = '';      clim = [0,500];
% $$$ normf = '';      clim = [0,100];
% $$$ normf = 'xprob'; clim = [0,0.15];
ind = [Trial.stc{'s+w+p+r+n'}];
%ind = [Trial.stc{'w+p+r+n'}];
%ind = [Trial.stc{'r'}];
%ind = [Trial.stc{'s'}];
subplot(251);
hist2([tpm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
caxis(clim)
subplot(252);
hist2([gpm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
caxis(clim)
subplot(253);
hist2([mpm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
caxis(clim)
subplot(254);
hist2([ipm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
caxis(clim)
subplot(255);
hist2([hpm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
caxis(clim)
subplot(256);
hist2([tpm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
caxis(clim)
subplot(257);
hist2([gpm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
caxis(clim)
subplot(258);
hist2([mpm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
caxis(clim)
subplot(259);
hist2([ipm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
caxis(clim)
subplot(2,5,10);
hist2([hpm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
caxis(clim)



ind = [Trial.stc{'s+w+r+n+p'}];
tpmm = copy(tpm); tpmm.data = tpmm.data-mean(tpmm(ind));
ipmm = copy(ipm); ipmm.data = ipmm.data-mean(ipmm(ind));
%mpmm = copy(mpm); mpmm.data = mpmm.data-mean(mpmm(ind));


%ind = [Trial.stc{'s'}];
figure
hist2([tpmm(ind),ipmm(ind)],linspace(-1.25,1.25,50),linspace(-1.25,1.25,50),normf);
caxis(clim)

ind = [Trial.stc{'s+w+r+n+p'}];
[U,S,V] = svd([tpm(ind)-mean(tpm(ind)),ipm(ind)-mean(ipm(ind))],0);
V = [0.771753449181619,        -0.635921861297655;
     0.635921861297655,         0.771753449181619];


%[U,S,V] = svd([tpm(ind)-mean(tpm(ind)),mpm(ind)-mean(mpm(ind))],0);
% $$$ figure,
% $$$ plot(U(:,1),U(:,2),'.')
% $$$ 
% $$$ figure
% $$$ hist2([U(:,1),U(:,2)],linspace(-3.2e-3,-2.2e-3,50),linspace(-8e-3,8e-3,50))


figure,
normf = '';
subplot(221);
ind = [Trial.stc{'w+r+n+p+s'}];
hist2(multiprod([tpmm(ind),ipmm(ind)],V,2,[1,2]),linspace(-1.25,1.25,50),linspace(-.5,.5,50),normf);
%line([-6.85,-5.5],[-0.4,0.4],'Color','r');
line([-.7,0.8],[0.5,-0.5],'Color','r');
subplot(222);
ind = [Trial.stc{'p'}];
hist2(multiprod([tpmm(ind),ipmm(ind)],V,2,[1,2]),linspace(-1.25,1.25,50),linspace(-.5,.5,50),normf);
%line([-6.85,-5.5],[-0.4,0.4],'Color','r');
line([-.7,0.8],[0.5,-0.5],'Color','r');
subplot(223);
ind = [Trial.stc{'w+r+n'}];
hist2(multiprod([tpmm(ind),ipmm(ind)],V,2,[1,2]),linspace(-1.25,1.25,50),linspace(-.5,.5,50),normf);
%line([-6.85,-5.5],[-0.4,0.4],'Color','r');
line([-.7,0.8],[0.5,-0.5],'Color','r');
subplot(224);
ind = [Trial.stc{'s'}];
hist2(multiprod([tpmm(ind),ipmm(ind)],V,2,[1,2]),linspace(-1.25,1.25,50),linspace(-.5,.5,50),normf);
%line([-6.85,-5.5],[-0.4,0.4],'Color','r');
line([-.7,0.8],[0.5,-0.5],'Color','r');


% $$$ dc = accumulate_decoding_vars(Trial,                               ...
% $$$                               units{trialId},                      ...
% $$$                               sessionList(trialId).thetaRefGeneral,...
% $$$                               phzCorrection(trialId),              ...
% $$$                               headRotation{trialId},               ...
% $$$                               hbangCorrection{trialId});

figure
subplot(141);hold('on');
    ind = [Trial.stc{'s'}]; 
    set(histogram(lpfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'); 
    ind = [Trial.stc{'w'}]; 
    set(set(histogram(lpfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
    xlim([5.5,7.5]);
subplot(142);hold('on');
    ind = [Trial.stc{'s'}]; 
    set(histogram(tdfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'); 
    ind = [Trial.stc{'w'}]; 
    set(set(histogram(tdfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
    xlim([-0.1,0.2]);
subplot(143);hold('on');
    ind = [Trial.stc{'s'}]; 
    set(histogram(rfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'); 
    ind = [Trial.stc{'w'}]; 
    set(set(histogram(rfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
    xlim([4,6]);
subplot(144);hold('on');
    ind = [Trial.stc{'s'}]; 
    set(histogram(dfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'); 
    ind = [Trial.stc{'w'}]; 
    set(set(histogram(dfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
    xlim([5,7.5]);

    

    
afetLims = [-1,1.5; -0.175,0.1; -1,2; -1.5,1.25];
figure
for f = 1:size(afet,2)
subplot(1,size(afet,2),f);hold('on');
    ind = [Trial.stc{'s'}]; 
    set(histogram(afet(ind,f)-afetMeanW(f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'); 
    ind = [Trial.stc{'w'}]; 
    set(set(histogram(afet(ind,f)-afetMeanW(f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
    ind = [Trial.stc{'n+p+r'}]; 
    set(set(histogram(afet(ind,f)-afetMeanW(f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
    %ind = abs(nunity(afet(:,3)))<2 & sum(ufr(:,:),2)./size(ufr,2)>0.8;
    ind = sum(fufr(ind,:),2)./size(fufr,2)>10;
    set(set(histogram(afet(ind,f)-afetMeanW(f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
    xlim(afetLims(f,:));
end    

mfufr = copy(fufr);
mfufr.data = sum(fufr(ind,:),2)./size(fufr,2);
figure
fufrSteps = 0:4:16;
for f = 1:size(afet,2)
subplot(1,size(afet,2),f);hold('on');
    for u = 1:numel(fufrSteps)-1
        ind = mfufr(:,1)>fufrSteps(u) & mfufr(:,1)<fufrSteps(u+1);
        set(set(histogram(afet(ind,f)-afetMeanW(f),linspace([afetLims(f,:),100]),'Normalization', ...
                          'probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
    end
    xlim(afetLims(f,:));
end    


figure,
%ind = abs(nunity(afet(:,3)))<;
ind = ':';
hist2([afet(ind,3),sum(fufr(ind,:),2)./size(ufr,2)],25,25);

figure,
ind = sum(fufr(ind,:),2)./size(fufr,2)>10;
hist2([afet(ind,1),afet(ind,3)],25,25);



% $$$ xcomp.data = ;
% $$$ xcomp.edgs = linspace(-2,2,nbins);
% $$$ xcomp.data = mlhfet(:,2);
% $$$ xcomp.edgs = linspace(-1.5,1.5,nbins);


afet = copy(lpfet);    
%afet.data = [lpfet(:,1),tdfet(:,1),rfet(:,1)];
%afet.data = [lpfet(:,1),tdfet(:,1),rfet(:,1),dfet(:,1),mlhfet(:,1),mlhfet(:,2)];
%afet.data = [lpfet(:,1),tdfet(:,1),rfet(:,1),lfet(:,1),hfet(:,1)];
afet.data = [lpfet(:,1),tdfet(:,1),tdrfet(:,1),rfet(:,1),tfet(:,1),trfet(:,1),lfet(:,1),hfet(:,1)];
afetMeanW = median(afet(Trial.stc{'w'},:),'omitnan');
afet.data = bsxfun(@minus,afet.data,afetMeanW);

afetLims = [-1,1.5; -0.175,0.1; -1,2; -1.5,1.25; -1.5,1.5; -1.5,1.5];
figure
for f = 1:size(afet,2)
subplot(1,size(afet,2),f);hold('on');
    ind = [Trial.stc{'s'}]; 
    set(histogram(afet(ind,f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'); 
    ind = [Trial.stc{'w'}]; 
    set(set(histogram(afet(ind,f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
    ind = [Trial.stc{'n+p+r'}]; 
    set(set(histogram(afet(ind,f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
    %ind = abs(nunity(afet(:,3)))<2 & sum(ufr(:,:),2)./size(ufr,2)>0.8;
    ind = sum(fufr(ind,:),2)./size(fufr,2)>10;
    set(set(histogram(afet(ind,f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
    xlim(afetLims(f,:));
end    



figure();
nind = nniz(afet);
for a = 1:size(afet,2)
    for b = a+1:size(afet,2)
        subplot2(size(afet,2),size(afet,2),a,b);
        hist2([afet(nind,a),afet(nind,b)],50,50);
    end
end

clear('hmm');
updateOM = 1;
hmm.K = 5;
hmm = hmminit(afet.data,hmm,'full');
hmm.train.cyc = 100;
hmm.obsmodel='Gauss';
hmm.train.obsupdate=ones([1,hmm.K])*updateOM;
hmm.train.init = 1;
hmm = hmmtrain(afet.data,size(afet,1),);
% trcpow hrcpow tpow

load(fullfile(Trial.path.project,'analysis',['test_hmm_model.mat']),'hmm');

hmm.P

[decode] = hmmdecode(afet.data,size(afet,1),hmm);
decode.q_star = decode.q_star';


figure();
subplot(211);
plot([1:size(decode.q_star,1)]./sampleRate,decode.q_star)
ylim([0,hmm.K]+0.5);
subplot(212);
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');



clear('xcomp','ycomp');
% $$$ xcomp.data = rfet(:,1);
% $$$ xcomp.edgs = linspace(3.5,6.5,30);
xcomp.data = lpfet(:,1);
xcomp.edgs = linspace(5.5,7.5,40);
ycomp.data = tdfet(:,1);
ycomp.edgs = linspace(-0.1,0.2,40);
% $$$ ycomp.data = dfet(:,1);
% $$$ ycomp.edgs = linspace(5,8,40);

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
    mind = dc.ucnt>=1 ...
           & dc.stcm(:,1)~=1 ...
    & (  dc.stcm(:,4)==4 ...
       | dc.stcm(:,6)==6) ...
    & decode.q_star(dc.ind)==grp;
% $$$ mind = dc.ucnt>=1 & (dc.stcm(:,1)~=1&dc.stcm(:,8)==8 ) & decode.q_star(dc.ind)==grp;
% $$$ mind = dc.ucnt>1 & (dc.stcm(:,1)==1&dc.stcm(:,8)==8 ) & decode.q_star(dc.ind)==grp;
% $$$ mind = dc.ucnt>1 & (dc.stcm(:,1)~=1&dc.stcm(:,7)==7 ) & decode.q_star(dc.ind)==grp;
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



figure,
for grp = 1:hmm.K
    mind = dc.ucnt>1 & (dc.stcm(:,1)~=1&dc.stcm(:,7)==7 ) & decode.q_star(dc.ind)==grp;
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
subplot2(1,hmm.K,1,grp);
    set(histogram(derr(nniz(derr)),linspace(0,400,50)),'EdgeColor','none');
end    

cmat = [sum(stcm(:,1)==1 & decode.q_star==1),sum(stcm(:,1)==1 & decode.q_star==2),sum(stcm(:,1)==1 & decode.q_star==3),sum(stcm(:,1)==1 & decode.q_star==4),sum(stcm(:,1)==1 & decode.q_star==5);...
sum(stcm(:,1)~=1 & decode.q_star==1),sum(stcm(:,1)~=1 & decode.q_star==2),sum(stcm(:,1)~=1 & decode.q_star==3),sum(stcm(:,1)~=1 & decode.q_star==4),sum(stcm(:,1)~=1 & decode.q_star==5)]
bsxfun(@rdivide,cmat,sum(cmat))




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
