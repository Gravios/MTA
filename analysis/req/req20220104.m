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
trialId = 18;

Trial = Trials{trialId};
unitSubset = units{trialId};
meta = sessionList(trialId);


Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',meta.subject.channelGroup.theta);
glfp = [diff(Trial.load('lfp',[33,40]).data,1,2),...
        diff(Trial.load('lfp',[41,48]).data,1,2),...
        diff(Trial.load('lfp',[49,56]).data,1,2),...
        diff(Trial.load('lfp',[57,64]).data,1,2)];
flfp = [diff(get(Trial.load('lfp',[33,40]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2),...
        diff(get(Trial.load('lfp',[41,48]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2),...
        diff(get(Trial.load('lfp',[49,56]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2),...
        diff(get(Trial.load('lfp',[57,64]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2)];
% $$$ glfp = [diff(Trial.load('lfp',[ 1, 8]).data,1,2),...
% $$$         diff(Trial.load('lfp',[ 9,16]).data,1,2),...
% $$$         diff(Trial.load('lfp',[17,24]).data,1,2),...
% $$$         diff(Trial.load('lfp',[25,32]).data,1,2)];
figure();
plot(glfp(:,[3,4]));
hold('on');
plot(flfp(:,[3,4]));

xyz = preproc_xyz(Trial,'trb',sampleRate);
fxyz = filter(xyz.copy(),'ButFilter',3,14,'low');
vxy = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
lvxy = copy(vxy);
lvxy.data(lvxy.data<=0.0001) = 0.0001;
lvxy.data = log10(lvxy.data);


unitsInt = Trial.spk.get_unit_set(Trial,'interneurons');

int = Trial.load('spk', Trial.lfp.sampleRate, '', unitsInt, '');

% LOAD theta phase
pchan = [1,8,9,16,17,24,25,32,33,40,41,48,49,56,57,64];
for c = 1:numel(pchan)
    phz = load_theta_phase(Trial,...
                           Trial.lfp.sampleRate,...
                           pchan(c),...
                           meta.subject.correction.thetaPhase);
    for ii = 1:numel(unitsInt),
        phzMean(c,ii) = circ_mean(phz(int(unitsInt(ii))));
        phzR(c,ii) = circ_r(phz(int(unitsInt(ii))));
    end
end

figure,plot(phzMean(:,5),phzR(:,5),'-+')


phz = load_theta_phase(Trial,...
                       Trial.lfp.sampleRate,...
                       meta.subject.channelGroup.theta,...
                       meta.subject.correction.thetaPhase);


figure,
for ii = 1:numel(unitsInt),
    subplot(2,7,ii);
    %intPhzPref(ii) = circ_mean(phz(int(unitsInt(ii))));
    rose(phz(int(unitsInt(ii))),37);
end

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


figure,
for ii = 1:numel(unitsInt),
    subplot(2,7,ii);
    %intPhzPref(ii) = circ_mean(phz(int(unitsInt(ii))));
    rose(phz(int(unitsInt(ii))),37);
end


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
%ufrl = Trial.load('ufr', lvxy, [], unitsInt, 0.24, 'boxcar');

phzThresh = 3;
figure
subplot(211);hold('on');
plot([1:size(vxy)]./vxy.sampleRate,RectFilter(sum(ufr.data,2)./numel(unitsInt),9,3));
% $$$ plot([1:size(vxy)]./vxy.sampleRate,RectFilter(sum(ufr(:,mpv>phzThresh),2)./sum(mpv>phzThresh),9,3));
% $$$ plot([1:size(vxy)]./vxy.sampleRate,RectFilter(sum(ufr(:,mpv<3),2)./sum(mpv<3),9,3));
plot([1:size(vxy)]./vxy.sampleRate,log10(RectFilter(sum(ufr(:,mpv>phzThresh),2)./sum(mpv>phzThresh),9,3)./RectFilter(sum(ufr(:,mpv<2.5),2)./sum(mpv<2.5),9,3)));
plot([1:size(vxy)]./vxy.sampleRate,RectFilter(sum(ufr.data,2)./numel(unitsInt),9,3)- ...
     log10(RectFilter(sum(ufr(:,mpv>phzThresh),2)./sum(mpv>phzThresh),9,3)./RectFilter(sum(ufr(:,mpv<3),2)./sum(mpv<3),9,3)));
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



figure,
for u = 1:14,
    subplot(2,7,u);
% $$$     sufr = ufr(dc.ind(((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2)),u);
% $$$     dcom = sqrt(sum(dc.ecom((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2,:).^2,2));
    sufr = ufr(dc.ind((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1)),u);
    dcom = sqrt(sum(dc.ecom((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1),:).^2,2));
    nind = nniz(sufr) & nniz(dcom) & sufr>0.0001;
    plot(sufr(nind),dcom(nind),'.');
    udcor(u) = corr(sufr(nind),dcom(nind));
end


figure,
for u = 1:14,
    subplot(2,7,u);
    sufr = ufr(dc.ind(((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2)),u);
    dcom = sqrt(sum(dc.ecom((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2,:).^2,2));
% $$$     sufr = ufr(dc.ind((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1)),u);
% $$$     dcom = lvxy(dc.ind((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1)),2);
    nind = nniz(sufr) & nniz(dcom) & sufr>0.0001;
    plot(sufr(nind),dcom(nind),'.');
    uvcor(u) = corr(sufr(nind),dcom(nind));
end



stcm = stc2mat(Trial.stc,xyz,{'theta','rear','loc','pause','sit'});
gufr = Trial.load('ufr',xyz,int,unitsInt,0.125,'gauss');
%ind = (stcm(:,2)==2|stcm(:,3)==3|stcm(:,4)==4|stcm(:,5)==5);
ind = (stcm(:,5)==5);
%ind = (stcm(:,4)==4|stcm(:,5)==5) & stcm(:,1)~=1;
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

[W,H] = nnmf(suI,2);

p1 = multiprod(gufr(:,:),H(1,:)',2,1);
p2 = multiprod(gufr(:,:),H(2,:)',2,1);
figure,plot(p1);hold('on');plot(p2);
figure,plot(log(p1./p2))


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
for i =  1:nClust
    out(:,i) = -.5*log(det(kCov(:,:,i)))...
        -.5*dot((bsxfun(@minus,gufr(:,:),median(suI(idx==vind,:)))/kCov(:,:,i))',...
                 bsxfun(@minus,gufr(:,:),median(suI(idx==vind,:)))')';
end
figure,plot(out)


figure();
ind = (stcm(:,5)==5);
plot(sum(gufr(ind,:),2)./14,lvxy(ind,2),'.')
hold('on');
ind = (stcm(:,3)==3|stcm(:,4)==4);
plot(sum(gufr(ind,:),2)./14,lvxy(ind,2),'.r')

udThresh = 0.1;
udThresh = 0.0;
irat = copy(vxy);
irat.data = log10(  RectFilter(sum(ufr(:,udcor>0&abs(udcor)>udThresh),2)./sum(udcor>0&abs(udcor)>udThresh),9,3) ...
                  ./RectFilter(sum(ufr(:,udcor<0&abs(udcor)>udThresh),2)./sum(udcor<0&abs(udcor)>udThresh),9,3));

frat = copy(vxy);
frat.data = log10(  RectFilter(sum(ufr(:,udcor>0&abs(udcor)>udThresh),2)./sum(udcor>0&abs(udcor)>udThresh),21,3) ...
                  ./RectFilter(sum(ufr(:,udcor<0&abs(udcor)>udThresh),2)./sum(udcor<0&abs(udcor)>udThresh),21,3));

uvThresh = 0.0;
vrat = copy(vxy);
vrat.data = log10(  RectFilter(sum(ufr(:,udcor>0&abs(udcor)>uvThresh),2)./sum(udcor>0&abs(udcor)>uvThresh),9,3) ...
                  ./RectFilter(sum(ufr(:,udcor<0&abs(udcor)>uvThresh),2)./sum(udcor<0&abs(udcor)>uvThresh),9,3));

rat = copy(vxy);
nrat.data = log(RectFilter(p1,21,3)./RectFilter(p2,21,3));

prat = copy(vxy);
prat.data = log(RectFilter(p1,31,3)+RectFilter(p2,31,3));


% $$$ rdThresh = randperm(14);
% $$$ irat.data = log10(  RectFilter(sum(ufr(:,rdThresh(1:7)),2)./7,9,3) ...
% $$$                   ./RectFilter(sum(ufr(:,rdThresh(8:14)),2)./7,9,3));

phzThresh = 2.5;
phzThresh = 3;
trat = copy(vxy);
trat.data = log10(  RectFilter(sum(ufr(:,mpv>phzThresh),2)./sum(mpv>phzThresh),9,3)  ...
                  ./RectFilter(sum(ufr(:,mpv<phzThresh),2)./sum(mpv<phzThresh),9,3));



figure
subplot(211);hold('on');
%plot([1:size(vxy)]./vxy.sampleRate,vrat.data);
plot([1:size(vxy)]./vxy.sampleRate,p1);
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


glvxy = vel(filter(preproc_xyz(Trial,'trb',sampleRate),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
glvxy.resample(lfp);
glvxy.data(glvxy.data<=0.0001) = 0.0001;
glvxy.data = log10(glvxy.data);
elfp = copy(lfp);
elfp.data = RectFilter(flfp(:,rcChan),3,3);
gphz = elfp.phase();

%fwin = 128;
rcChan = 3;
%rcChan = 1;
fwin = gausswin(64);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = RectFilter(flfp(:,rcChan),3,3);
%elfp.filter('ButFilter', 4, [30,50], 'bandpass');
elfp.filter('ButFilter', 4, [50,75], 'bandpass');
%elfp.filter('ButFilter', 4, [75,125], 'bandpass');
%elfp.filter('ButFilter', 4, [125,175], 'bandpass');
gfet = copy(elfp);
gfet.data = log10(conv(gfet.data.^2,fwin,'same'));
% $$$ gfet.filter('ButFilter', 4, [1], 'low');
% $$$ gfet.resample(xyz);
% $$$ figure;plot(unity(ipm.data));hold('on');plot(unity(hfet.data));
% $$$ figure;plot(unity(gfet.data));hold('on');plot(unity(hfet.data));

[gmin,gval] = LocalMinima(-gfet.data,100,0);

tper = Trial.stc{'t'};
gmin = gmin(WithinRanges(gmin,tper.data));

figure,
%shift = [-1,-1,0];
shift = [0,0,0];
gvmin = gmin(glvxy(gmin,1)>0&glvxy(gmin,1)<0.5);
gvval = gval(glvxy(gmin,1)>0&glvxy(gmin,1)<0.5);
subplot(231);
hist2([gphz(gvmin),-gvval],linspace(-pi,pi,13),linspace([4,6,13]+shift));
subplot(234);
hist2([gphz(gvmin),-gvval],linspace(-pi,pi,13),linspace([4,6,13]+shift),'xprob');
caxis([0,0.15])
gvmin = gmin(glvxy(gmin,1)>0.5&glvxy(gmin,1)<1);
gvval = gval(glvxy(gmin,1)>0.5&glvxy(gmin,1)<1);
subplot(232);
hist2([gphz(gvmin),-gvval],linspace(-pi,pi,13),linspace([4,6,13]+shift));
subplot(235);
hist2([gphz(gvmin),-gvval],linspace(-pi,pi,13),linspace([4,6,13]+shift),'xprob');
caxis([0,0.15])
gvmin = gmin(glvxy(gmin,1)>1&glvxy(gmin,1)<2);
gvval = gval(glvxy(gmin,1)>1&glvxy(gmin,1)<2);
subplot(233);
hist2([gphz(gvmin),-gvval],linspace(-pi,pi,13),linspace([4,6,13]+shift));
subplot(236);
hist2([gphz(gvmin),-gvval],linspace(-pi,pi,13),linspace([4,6,13]+shift),'xprob');
caxis([0,0.15])
colormap('jet');


sgmins = round(gmin./1250*30);

gvmin = sgmins(lvxy(sgmins,1)>0.5&lvxy(sgmins,1)<2);
gvval = gval(lvxy(sgmins,1)>0.5&lvxy(sgmins,1)<2);
figure,plot(sqrt(sum(dc.ecom(dc.ind(ismember(gvmin,dc.ind)),:).^2,2)),gvval(ismember(gvmin,dc.ind)),'.')
figure,hist2([sqrt(sum(dc.ecom(dc.ind(ismember(gvmin,dc.ind)),:).^2,2)),...
              -gvval(ismember(gvmin,dc.ind))],...
             linspace(0,800,30),...
             linspace(4,6,20),'xprob');

%fwin = 128;
rcChan = 3;
%rcChan = 1;
fwin = gausswin(64);
%fwin = gausswin(128);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = RectFilter(glfp(:,rcChan),3,3);
%elfp.filter('ButFilter', 4, [50,200], 'bandpass');
%elfp.filter('ButFilter', 4, [50,100], 'bandpass');
elfp.filter('ButFilter', 4, [75,125], 'bandpass');
hfet = copy(elfp);
hfet.data = log10(conv(hfet.data.^2,fwin,'same'));
hfet.filter('ButFilter', 4, [1.5], 'low');
hfet.resample(xyz);
% $$$ figure;plot(unity(ipm.data));hold('on');plot(unity(hfet.data));


fwin = gausswin(128);
elfp = copy(lfp);
elfp.data = glfp(:,rcChan);
elfp.filter('ButFilter', 4, [1,20], 'bandpass');
lfet = copy(elfp);
lfet.data = log10(conv(lfet.data.^2,fwin,'same'));
lfet.resample(xyz);
lfet.filter('ButFilter', 4, [1.5], 'low');
% $$$ figure;plot(nunity(tpm.data));hold('on');plot(nunity(lfet.data));


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


%ind = [Trial.stc{'s+w+n+r+p'}];
ind = [Trial.stc{'w+n+p'}];
[Uf,Sf,Vf] = svd([lfet(ind)-mean(lfet(ind)),hfet(ind)-mean(hfet(ind))],0);
% $$$ Vf = [0.771753449181619,        -0.635921861297655;
% $$$      0.635921861297655,         0.771753449181619];
% $$$ Vf = [-0.83534705515978,         0.549722927879022; ...
% $$$         -0.549722927879022,         -0.83534705515978];


mlfet = lfet.copy();  mlfet.data = mlfet.data - mean(mlfet(ind));
mhfet = hfet.copy();  mhfet.data = mhfet.data - mean(mhfet(ind));
mlhfet = lfet.copy();
mlhfet.data = multiprod([mlfet.data,mhfet.data],Vf,2,[1,2]);


velBins = linspace(-2.5,1.8,20);
velInds = discretize(lvxy(:,2),velBins);
figure,
for vind = 1:nClust,
    subplot2(nClust,2,vind,1);
        imagesc(kCov(:,:,vind));
    subplot2(nClust,2,vind,2);
    histogram(svI(idx==vind,2),velBins);
end


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

% $$$ mlhfet = lfet.copy();
% $$$ mlhfet.data = [lfet.data,hfet.data];

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



% $$$ elfp = copy(lfp);
% $$$ elfp.data = diff(glfp(1:1e6,[1,3]),1,2);
elfp = copy(lfp);
elfp.data = RectFilter(glfp,3,3);
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  elfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,250]);
[lys,lfs,lts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);


% $$$ elfp = copy(lfp);
% $$$ elfp.data = diff(glfp(1:1e6,[1,3]),1,2);
elfp = copy(lfp);
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  elfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,30]);
[rys,rfs,rts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);


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


tipmmv = multiprod([tpmm.data,ipmm.data],V,2,[1,2]);
%tipmmv = multiprod([tpmm.data,mpmm.data],V,2,[1,2]);


xcomp.data = tipmmv(:,1);
xcomp.edgs = linspace(-1,1,30);
ycomp.data = tipmmv(:,2);
ycomp.edgs = linspace(-.4,.4,30)

% trcpow hrcpow tpow

mind = dc.ucnt>2 & (dc.stcm(:,3)==3 |dc.stcm(:,4)==4);
mind = dc.ucnt>2 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5 | dc.stcm(:,4)==4 |dc.stcm(:,6)==6);
mind = dc.ucnt>2 & (dc.stcm(:,4)==4 |dc.stcm(:,6)==6);
mind = dc.ucnt>1 & (dc.stcm(:,8)==8);
mind = dc.ucnt>1 & (dc.stcm(:,1)~=1&dc.stcm(:,8)==8 );
mind = dc.ucnt>1 & (dc.stcm(:,1)==1&dc.stcm(:,8)==8 );
mind = dc.ucnt>1 & dc.stcm(:,1)==1;
mind = dc.ucnt>1 & dc.stcm(:,2)==2;
mind = dc.ucnt>1;

derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,derr);

derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = dc.ecom(mind,1);
[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,derr);

derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,derr);

derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = ipm(dc.ind(mind));
[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,derr);

derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = tpm(dc.ind(mind));
[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,derr);

derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = tpow(dc.ind(mind));
[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,derr);

derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = tdpow(dc.ind(mind));
[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,derr);

zmean = zcomp.mean;
zmean(zcomp.count<20) = nan;
zstd  = zcomp.std;
zstd(zcomp.count<20) = nan;

figure,
subplot(231);
imagesc(xcomp.ctrs,ycomp.ctrs,zmean');
Lines([],0,'w');Lines(0,[],'w');
axis('xy'); colormap('jet'); colorbar();
subplot(232);
imagesc(xcomp.ctrs,ycomp.ctrs,zstd');
Lines([],0,'w'); Lines(0,[],'w');
axis('xy'); colormap('jet'); colorbar();
subplot(234);
ind = [Trial.stc{'w+r+n+p+s'}]; cast(ind,'TimeSeries'); resample(ind,vxy); ind = logical(ind.data);
hist2([xcomp.data(ind),ycomp.data(ind)],xcomp.edgs,ycomp.edgs);
grid('on'); Lines([],0,'w'); Lines(0,[],'w'); colormap('jet'); colorbar();
subplot(235);
ind = [Trial.stc{'w+r+n'}]; cast(ind,'TimeSeries'); resample(ind,vxy); ind = logical(ind.data);
hist2([xcomp.data(ind),ycomp.data(ind)],xcomp.edgs,ycomp.edgs);
grid('on'); Lines([],0,'w'); Lines(0,[],'w'); colormap('jet'); colorbar();
subplot(236);
ind = [Trial.stc{'s'}]; cast(ind,'TimeSeries'); resample(ind,vxy); ind = logical(ind.data);
hist2([xcomp.data(ind),ycomp.data(ind)],xcomp.edgs,ycomp.edgs);
grid('on'); Lines([],0,'w'); Lines(0,[],'w'); colormap('jet'); colorbar();



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


%3dfig
clear('xcomp','ycomp','zcomp','wcomp');
nbins = 15;
% $$$ xcomp.data = tipmmv(:,1);
% $$$ xcomp.edgs = linspace(-1.4,1.4,nbins);
% $$$ ycomp.data = tipmmv(:,2);
% $$$ ycomp.edgs = linspace(-.5,.5,nbins);
% $$$ xcomp.data = mlhfet(:,1);
% $$$ xcomp.edgs = linspace(-7,11,nbins);
% $$$ ycomp.data = mlhfet(:,2);
% $$$ ycomp.edgs = linspace(3.5,5,nbins);
xcomp.data = mlhfet(:,1);
xcomp.edgs = linspace(-2,2,nbins);
ycomp.data = mlhfet(:,2);
ycomp.edgs = linspace(-1.5,1.5,nbins);
% $$$ zcomp.data = tdpow(:,1);
% $$$ zcomp.edgs = linspace(-5,2,nbins);
% $$$ zcomp.data = lvxy(:,2);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
zcomp.data = irat(:,1);
zcomp.edgs = linspace(-1.5,1.5,nbins);%zcomp.edgs = linspace(-3,3,nbins);
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



isothresh = [100];
isothresh = [50];
figure();
hold('on');

mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
pos = cell([1,3]);
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

mind = dc.ucnt>2 & (dc.stcm(:,4)==4 |dc.stcm(:,6)==6);
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

mind = dc.ucnt>1 & (dc.stcm(:,2)==2);
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

%loco
mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5 );
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
