
% The distance of the actual position of the rat's head to the
% decoded positon is normally taken as the decoding error. We
% projected the error onto the frame of reference of the rat's head
% to decompose the lateral and anteroposterior components to find
% if the error varies with behavioral variables normally associated
% with navigation.
%
% Egocentric 
% Movement direction has been demonstrated to be correlated with
% the ascending phase of theta, where EC3 input drives the CA1
% place cells (Huxter&csicsvari 2008). 
%
% If the 


% plot - maze <- posterior & subject at position with axes over head frame of reference
%

% TODO - generate lat(hba) subplots for 


configure_default_args();
EgoProCode2D_load_data();


rat = load_patch_model('rat');
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false);
subject = update_subject_patch(subject,'body',[],false);
              

% CA1
tind = [3,4,5,17,18,19,20,21,22,23,29];
%tind = [6,7,26,27,30];
sampleRate = 250;
halfSpikeWindow = 0.015;

global AP
% compute_ratemaps ---------------------------------------------------------------------------------
AP.compute_ratemaps =                                                                            ...
    struct('get_featureSet',            @fet_xy,                                                 ...
           'sampleRate',                16,                                                      ...
           'pfsArgs',                   struct('states',           'theta-groom-sit-rear',       ...
                                               'binDims',          [50,50],                      ...
                                               'SmoothingWeights', [2.4,2.4],                    ...
                                               'numIter',          1,                            ...
                                               'boundaryLimits',   [-500,500;-500,500],          ...
                                               'halfsample',       false)                        ...
           );
%---------------------------------------------------------------------------------------------------                                                                                                                                                                                                                          

dca = cf(@(T,U) accumulate_decoding_vars( T, U, sampleRate, ...
                                           halfSpikeWindow), Trials(tind),units(tind));
%dcao = dca;
%dca = dcan


figure,
hold('on');
plot(dca{7}.ind./250,dca{7}.ecom(:,2))
plot(dca{7}.ind./250,dca{7}.phz*10,'r')
plot(dca{7}.ind./250,-dca{7}.hbang.data*100,'m')


lfp = Trials{20}.load('lfp',Trials{20}.meta.channelGroup.theta);
lfp.resample(sampleRate);


phz = load_theta_phase(Trials{20},sampleRate);
thetaTroughs = abs(phz.data-pi)<0.1&[0;diff(phz.data-pi)]>0;
thetaTroughs = LocalMinima(-convn(thetaTroughs,ones([1,5]),'same'),0,5);


figure,
subplot(311);
hold('on'),
ind = ismember(dca{7}.ind,thetaTroughs);
plot(dca{7}.ind(ind)./sampleRate,atan2(dca{7}.ecom(ind,2),dca{7}.ecom(ind,1))*10,'.b')
plot(dca{7}.ind(ind)./sampleRate,zeros([sum(ind),1]),'+m');
ind = ismember(dca{7}.ind,thetaTroughs+80);
plot(dca{7}.ind(ind)./sampleRate,atan2(dca{7}.ecom(ind,2),dca{7}.ecom(ind,1))*10,'.r')
plot(dca{7}.ind(ind)./sampleRate,dca{7}.hdist(ind,1)/10)
plot([1:size(lfp,1)]./sampleRate,lfp.data./500,'b')
Lines([],0,'k');
Lines([],300,'k');
subplot(312);
hold('on'),
ind = ismember(dca{7}.ind,thetaTroughs);
plot(dca{7}.ind(ind)./sampleRate,-dca{7}.hbang(ind,1),'.b')
plot(dca{7}.ind(ind)./sampleRate,zeros([sum(ind),1]),'+m');
ind = ismember(dca{7}.ind,thetaTroughs);
plot(dca{7}.ind(ind)./sampleRate,-dca{7}.hbang(ind,1),'.r')
Lines([],0,'k');
subplot(313);
plotSTC(Trials{20}.stc,1)
linkax('x')


figure
hold('on');
ind = ismember(dca{7}.ind,thetaTroughs)& dca{7}.stcm(:,1)==1 & dca{7}.stcm(:,3)==3;
histogram(atan2(dca{7}.ecom(ind,2),dca{7}.ecom(ind,1)),linspace([-pi,pi,50]))
ind = ismember(dca{7}.ind,thetaTroughs+100)& dca{7}.stcm(:,1)==1& dca{7}.stcm(:,3)==3;
histogram(atan2(dca{7}.ecom(ind,2),dca{7}.ecom(ind,1)),linspace([-pi,pi,50]))


figure
hold('on');
ind = ismember(dca{7}.ind,thetaTroughs)& dca{7}.stcm(:,1)==1 & ...
      dca{7}.stcm(:,3)==3 &dca{7}.hdist.data<300& dca{7}.hbang.data<-0.2 ;;
histogram(dca{7}.ecom(ind,2),linspace([-pi*100,pi*100,50]))
ind = ismember(dca{7}.ind,thetaTroughs+80) & dca{7}.stcm(:,1)==1& ...
      dca{7}.stcm(:,3)==3 &dca{7}.hdist.data<300& dca{7}.hbang.data<-0.2 ;;
histogram(dca{7}.ecom(ind,2),linspace([-pi*100,pi*100,50]))


figure
hold('on');
ind = ismember(dca{7}.ind,thetaTroughs)& dca{7}.stcm(:,1)==1 & ...
      dca{7}.stcm(:,4)==4 &dca{7}.hdist.data<300 & dca{7}.hbang.data<-0.2 ;

figure
ind = ismember(dca{7}.ind,thetaTroughs) & dca{7}.stcm(:,1)==1& ...
      dca{7}.stcm(:,3)==3 &dca{7}.hdist.data<300;
histogram2(dca{7}.ecom(ind,2),dca{7}.ecom(ind,1),...
           linspace([-pi*100,pi*100,20]),...
           linspace([-pi*100,pi*100,20]),'DisplayStyle','tile')



t = 5;
figure,
subplot(4,1,1:3)
hold('on');

ind = WithinRanges(dca{t}.phz,[4,6]) & dca{t}.ucnt>=2 & dca{t}.ucnt<=8;
cnv = dca{t}.ecom(ind,2)-conv(dca{t}.ecom(ind,2),ones([121,1])./121,'same');
plot(dca{t}.ind(ind)./250,dca{t}.ecom(ind,2))
plot(dca{t}.ind(ind)./250,dca{t}.ecom(ind,2),'.k')
plot(dca{t}.ind./250,-dca{t}.hbang(:,1)*100,'m')
plot(dca{t}.ind./250,dca{t}.phz*10,'r')
plot(dca{t}.ind(ind)./250,cnv,'k');
%plot(dca{t}.ind./250,dca{t}.bvfl(:,2)*10+1000,'k')
plot(dca{t}.ind./250,dca{t}.hvfl(:,1)*10,'b')
plot(dca{t}.ind./250,dca{t}.hvfl(:,2)*10,'g')
subplot(414)
plotSTC(Trials{tind(t)}.stc,1)
ylim([-1,7])
Lines([],0,'k');
linkx();


t = 7;
figure,
subplot(4,1,[1:3])
hold('on')
plot(dca{t}.ind./250,dca{t}.ucnt(:,1),'k')
subplot(414)
plotSTC(Trials{tind(t)}.stc,1);
ylim(gca(),[0,7]);
linkx();



t = 7;
figure,
subplot(4,1,[1:3])
%ind = WithinRanges(dca{t}.phz,[4,6]) & dca{t}.ucnt>=1 & dca{t}.ucnt<=8;
ind = dca{t}.ucnt>=1 & dca{t}.ucnt<=8;
hold('on')
plot(dca{t}.ind(ind)./250,dca{t}.ecom(ind,2),'k')
plot(dca{t}.ind(ind)./250,dca{t}.ecom(ind,2),'.r')
% $$$ plot(dca{t}.ind(ind)./250,atan2(dca{t}.ecom(ind,2),dca{t}.ecom(ind,1)),'k')
% $$$ plot(dca{t}.ind(ind)./250,atan2(dca{t}.ecom(ind,2),dca{t}.ecom(ind,1)),'.r')
plot(dca{t}.ind./250,dca{t}.phz*10,'r')
plot(dca{t}.ind./250,-dca{t}.hbang(:,1)*100,'m')
plot(dca{t}.ind./250,dca{t}.hdist(:,1),'c')
plot(dca{t}.ind./250,dca{t}.hvfl(:,2)*10,'g')
Lines([],0,'k');
subplot(414)
plotSTC(Trials{tind(t)}.stc,1);
ylim([0,7]);
linkx();

id = WithinRanges(dca{7}.phz,[4,6]) & dca{7}.ucnt>=2 & dca{7}.ucnt<=8 ...
     &nniz(dca{7}.ecom);

[IMF, RESIDUAL] = emd(dca{7}.ecom(id,2));

figure,
subplot(5,1,[1:3]);
ts = dca{7}.ind(id)./sampleRate;
offs =  0;
for j = 1:10
    plot(ts,IMF+offs);
    offs=offs+100;
end    
subplot(5,1,[4]);
plot(ts,dca{7}.phz(id));
subplot(5,1,[5]);
plotSTC(Trials{20}.stc,1)
ylim([-1,7])
linkx();




ind = WithinRanges(dca{7}.phz,[4,6]) & dca{7}.ucnt>=2 & dca{7}.ucnt<=4 ...
      & (dca{7}.stcm(:,3)==3)...
      ...      & (dca{7}.stcm(:,3)==3|dca{7}.stcm(:,4)==4|dca{7}.stcm(:,5)==5)...
      & dca{7}.hdist.data<300;
figure,plot(dca{7}.ecom(ind,2),dca{7}.bvfl(ind,2),'.')

figure,hist2([dca{7}.ecom(ind,2),dca{7}.bvfl(ind,2)],50,50,'log','mud')

figure,hist2([dca{7}.esax(ind,2),dca{7}.hbang(ind,1)],50,50,'log','mud')
figure,hist2([dca{7}.ecom(ind,2),dca{7}.hvang(ind,1)],50,50,'log','mud')

figure,hist2([dca{7}.esax(ind,2),dca{7}.hbang(ind,1)],50,50,'log','mud')


figure()
for t = 1:numel(dca)
ind = WithinRanges(dca{t}.phz,[4,6]) & dca{t}.ucnt>=2 & dca{t}.ucnt<=4 ...
      & (dca{t}.stcm(:,3)==3|dca{t}.stcm(:,4)==4|dca{t}.stcm(:,5)==5)...
      & dca{t}.hdist.data<300;
subplot(3,5,t)
hist2([dca{t}.esax(ind,2),-dca{t}.hbang(ind,1)],11,11,'log','mud')
colormap(gca,'jet')
end


s = 3
figure()
for t = 1:numel(dca)
ind = WithinRanges(dca{t}.phz,[4,6]) & dca{t}.ucnt>=2 & dca{t}.ucnt<=4 ...
      ...& (dca{t}.stcm(:,3)==3|dca{t}.stcm(:,4)==4|dca{t}.stcm(:,5)==5)...
       & (dca{t}.stcm(:,s)==s)...
      & dca{t}.hdist.data<300;
subplot2(4,11,1,t);
hist2([dca{t}.esax(ind,2),-dca{t}.hbang(ind,1)],7,7,'log','mud');
colormap(gca,'jet');
subplot2(4,11,2,t);
hist2([dca{t}.esax(ind,2),dca{t}.hvflP(ind,2)],7,7,'log','mud');
colormap(gca,'jet');
subplot2(4,11,3,t);
hist2([dca{t}.esax(ind,2),dca{t}.hvfl(ind,2)],7,7,'log','mud');
colormap(gca,'jet')
subplot2(4,11,4,t);
hist2([dca{t}.esax(ind,2),dca{t}.hvflF(ind,2)],7,7,'log','mud');
colormap(gca,'jet');
end



s = 5;
decoded = struct('fwd',[],...
                 'lat',[],...
                 'tfwd',[],...
                 'tlat',[],...
                 'hvf',[],...
                 'hvl',[],...
                 'bvl',[],...
                 'hvlF',[],...
                 'hvlP',[],...                 
                 'hav',[],...
                 'hba',[],...
                 'phz',[]);
for t = [1:3,5:8,11],
%for t = [1:11],
%for t = [5:8],
    %for t = [1:3,5:8,11],
mind =   dca{t}.stcm(:,1)==1                                              ...
...        & (dca{t}.stcm(:,3)==3|dca{t}.stcm(:,5)==5)   ...
        & (dca{t}.stcm(:,3)==3|dca{t}.stcm(:,4)==4|dca{t}.stcm(:,5)==5)   ...
...        & (dca{t}.stcm(:,s)==s)   ...                           
        &  dca{t}.hvfl(:,1)>-5                                            ...
        & dca{t}.ucnt>=1 & dca{t}.ucnt<=5                                 ...
        & sqrt(sum(dca{t}.xyz(:,'hcom',[1,2]).^2,3))>0                    ...
        & sqrt(sum(dca{t}.xyz(:,'hcom',[1,2]).^2,3))<300;
    decoded.fwd = cat(1,decoded.fwd,dca{t}.ecom(mind,1));
    decoded.lat = cat(1,decoded.lat,dca{t}.ecom(mind,2));%-10*double(t>4)+10*double(t<4));
    decoded.tfwd = cat(1,decoded.tfwd,dca{t}.tcom(mind,1));
    decoded.tlat = cat(1,decoded.tlat,dca{t}.tcom(mind,2));%+20*double(t>4));
    decoded.hvf = cat(1,decoded.hvf,dca{t}.hvfl(mind,1));
    decoded.hvl = cat(1,decoded.hvl,dca{t}.hvfl(mind,2));
    decoded.bvl = cat(1,decoded.bvl,dca{t}.bvfl(mind,2));    
    decoded.hvlF = cat(1,decoded.hvlF,dca{t}.hvflF(mind,2));
    decoded.hvlP = cat(1,decoded.hvlP,dca{t}.hvflP(mind,2));        
    decoded.hav = cat(1,decoded.hav,dca{t}.hvang(mind,1));
    decoded.hba = cat(1,decoded.hba,dca{t}.hbang(mind,1));
    decoded.phz = cat(1,decoded.phz,dca{t}.phz(mind,1));
end




ind = WithinRanges(decoded.phz,phzBin.edges(3:4)) ...
      & randn(size(decoded.hba))>1 ...
      & abs(decoded.hba)<0.4;
ind = WithinRanges(decoded.phz,[5,6]) ... 
     & randn(size(decoded.hba))>1 ...
      & abs(decoded.hba)<0.6;


ind = WithinRanges(decoded.phz,[1,2.5]) ... 
     & randn(size(decoded.hba))>1 ...
      & abs(decoded.hba)<0.4;

figure,
hist2([decoded.fwd(ind),decoded.hvf(ind)],linspace(-400,400,30),-5: ...
      10:55,'','mud');


ind = WithinRanges(decoded.phz,[5,6]) ... 
     & randn(size(decoded.hba))>0.5; ...
figure,
hist2([decoded.lat(ind),decoded.hba(ind)],linspace(-400,400,30),-1.2:0.12:1.2,'xprob','');


figure,
hist2([decoded.lat(ind),decoded.hav(ind)],linspace(-400,400,30),-.5:0.1:.5,'xprob','');

ind = WithinRanges(decoded.phz,[5,6]) ... 
      & randn(size(decoded.hba))>0.5 ...
      & decoded.hba<-0.2;      
figure,
hist2([decoded.lat(ind),decoded.hvl(ind)],linspace(-400,400,30),-40:10:40,'xprob','');
colormap('jet');

figure,
hold('on');
for hbaInd = 1:hbaBin.count
ind = WithinRanges(decoded.phz,[5,6]) ... 
      & randn(size(decoded.hba))>0.5 ...
      & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1]));      
%set(histogram(decoded.lat(ind),linspace(-300,300,40)),'EdgeColor','none','FaceAlpha',.3); 
[F,xi] = ksdensity(decoded.lat(ind)/10)
plot(xi,F,'-','Color',hbaBin.color(hbaInd,:));
end

figure();
hold('on');
histogram(decoded.fwd(ind&decoded.hvf<5 & decoded.hvf>-5 & abs(decoded.hvl<15)),linspace(-400,400,30),'Normalization','probability')
histogram(decoded.fwd(ind&decoded.hvf<50& decoded.hvf>40 & abs(decoded.hvl<15)),linspace(-400,400,30),'Normalization','probability')


medf = [...
median(decoded.fwd(ind&decoded.hvf<5 & decoded.hvf>-5)),...
median(decoded.fwd(ind&decoded.hvf<15 & decoded.hvf>5)),...
median(decoded.fwd(ind&decoded.hvf<25 & decoded.hvf>15)),...
median(decoded.fwd(ind&decoded.hvf<35 & decoded.hvf>25)),...
median(decoded.fwd(ind&decoded.hvf<45 & decoded.hvf>35)),...
median(decoded.fwd(ind&decoded.hvf<55& decoded.hvf>45))];
figure,plot(medf)

[h,p] = ttest2(decoded.fwd(ind&decoded.hvf<50&decoded.hvf>40),decoded.fwd(ind&decoded.hvf<20&decoded.hvf>10))


dhvf = discretize(decoded.hvf(ind),-5:10:55);
dfwd = decoded.fwd(ind);
nind = nniz([dfwd,dhvf]);
out = accumarray(dhvf(nind),dfwd(nind)-20,[6,1],@mean);
figure,plot(out)


ind = WithinRanges(decoded.phz,phzBin.edges(3:4)) ...
      & randn(size(decoded.hba))>1 ...
      & abs(decoded.hba)<1.2;
stats =[];
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind)]);
stats = cat(1,stats,STATS);
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hvlP(ind)]);
stats = cat(1,stats,STATS);
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hvl(ind)]);
stats = cat(1,stats,STATS);
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hvlF(ind)]);
stats = cat(1,stats,STATS);
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind),decoded.hvlP(ind)]);
stats = cat(1,stats,STATS);
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind),decoded.hvl(ind)]);
stats = cat(1,stats,STATS);
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind),decoded.hvlF(ind)]);
stats = cat(1,stats,STATS);
tstats =[];
[B,BINT,R,RINT,STATS] = regress(decoded.tlat(ind),[ones([sum(ind),1]),decoded.hba(ind)]);
tstats = cat(1,tstats,STATS);
[B,BINT,R,RINT,STATS] = regress(decoded.tlat(ind),[ones([sum(ind),1]),decoded.hvlP(ind)]);
tstats = cat(1,tstats,STATS);
[B,BINT,R,RINT,STATS] = regress(decoded.tlat(ind),[ones([sum(ind),1]),decoded.hvl(ind)]);
tstats = cat(1,tstats,STATS);
[B,BINT,R,RINT,STATS] = regress(decoded.tlat(ind),[ones([sum(ind),1]),decoded.hvlF(ind)]);
tstats = cat(1,tstats,STATS);
[B,BINT,R,RINT,STATS] = regress(decoded.tlat(ind),[ones([sum(ind),1]),decoded.hba(ind),decoded.hvlP(ind)]);
tstats = cat(1,tstats,STATS);
[B,BINT,R,RINT,STATS] = regress(decoded.tlat(ind),[ones([sum(ind),1]),decoded.hba(ind),decoded.hvl(ind)]);
tstats = cat(1,tstats,STATS);
[B,BINT,R,RINT,STATS] = regress(decoded.tlat(ind),[ones([sum(ind),1]),decoded.hba(ind),decoded.hvlF(ind)]);
tstats = cat(1,tstats,STATS);

round(stats,5)
round(tstats,5)
[std(decoded.lat(ind)),std(decoded.tlat(ind))]


[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hvl(ind)]);
STATS

[B,BINT,R,RINT,STATS] = regress(decoded.fwd(ind),[ones([sum(ind),1]),decoded.hvf(ind)]);
STATS



figure,hist2([decoded.hvl,decoded.hba],linspace([-80,80,24]),linspace([-1.5,1.5,24]));

ind = decoded.hba>0.2& randn(size(decoded.hba))>0.5; [mean(decoded.lat(ind)),mean(decoded.tlat(ind))]
ind = abs(decoded.hba)<0.2& randn(size(decoded.hba))>0.5;[mean(decoded.lat(ind)),mean(decoded.tlat(ind))]
ind = decoded.hba<-0.2& randn(size(decoded.hba))>0.5;[mean(decoded.lat(ind)),mean(decoded.tlat(ind))]
figure,hist2([decoded.tlat(ind),decoded.lat(ind)],linspace([-400,400,24]),linspace([-400,400,24]));
set(gca(),'ColorScale','log');
colormap(gca(),'jet');
Lines([],0,'w');
Lines(0,[],'w');
[mean(decoded.lat(ind)),mean(decoded.tlat(ind))]

norm = 'xprob';
%norm = '';
nlim = [0,0.2];
figure,
for p = 1:phzBin.count
ind = WithinRanges(decoded.phz,phzBin.edges([p,p+1])) ...
      & randn(size(decoded.hba))>0;

% $$$ [B,BINT,RL,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind)]);
% $$$ STATS
% $$$ figure,cdfplot(RL)
% $$$ [B,BINT,RT,RINT,STATS] = regress(decoded.tlat(ind),[ones([sum(ind),1]),decoded.hba(ind)]);
% $$$ STATS
% $$$ hold('on');
% $$$ cdfplot(RT);
% $$$ xlim([-300,300])

subplot2(6,4,phzBin.count+1-p,1);
    hist2([decoded.lat(ind),decoded.hba(ind)],linspace(-300,300,18),linspace(-1.2,1.2,18),norm)
    colormap('jet'); Lines([],0,'w'); Lines(0,[],'w'); caxis(nlim);
subplot2(6,4,phzBin.count+1-p,2);
    hist2([decoded.tlat(ind),decoded.hba(ind)],linspace(-300,300,18),linspace(-1.2,1.2,18),norm)
    colormap('jet'); Lines([],0,'w'); Lines(0,[],'w'); caxis(nlim);
subplot2(6,4,phzBin.count+1-p,3);
    hist2([decoded.lat(ind),decoded.hvl(ind)],linspace(-300,300,18),linspace(-60,60,18),norm)
    colormap('jet'); Lines([],0,'w'); Lines(0,[],'w'); caxis(nlim);
subplot2(6,4,phzBin.count+1-p,4);
    hist2([decoded.tlat(ind),decoded.hvl(ind)],linspace(-300,300,18),linspace(-60,60,18),norm)
    colormap('jet'); Lines([],0,'w'); Lines(0,[],'w'); caxis(nlim);

subplot2(6,4,2*phzBin.count+1-p,1);
    hist2([decoded.fwd(ind),decoded.hba(ind)],linspace(-300,300,18),linspace(-1.2,1.2,18),norm)
    colormap('jet'); Lines([],0,'w'); Lines(0,[],'w'); caxis(nlim);
subplot2(6,4,2*phzBin.count+1-p,2);
    hist2([decoded.tfwd(ind),decoded.hba(ind)],linspace(-300,300,18),linspace(-1.2,1.2,18),norm)
    colormap('jet'); Lines([],0,'w'); Lines(0,[],'w'); caxis(nlim);
subplot2(6,4,2*phzBin.count+1-p,3);
    hist2([decoded.fwd(ind),decoded.hvl(ind)],linspace(-300,300,18),linspace(-60,60,18),norm)
    colormap('jet'); Lines([],0,'w'); Lines(0,[],'w'); caxis(nlim);
subplot2(6,4,2*phzBin.count+1-p,4);
    hist2([decoded.tfwd(ind),decoded.hvl(ind)],linspace(-300,300,18),linspace(-60,60,18),norm)
    colormap('jet'); Lines([],0,'w'); Lines(0,[],'w'); caxis(nlim);
end




%mBinHvl.edges = linspace(-35,35,8);
%mBinHvl.edges = linspace(-35,35,8);
mBinHvl.edges = linspace(-0.5,0.5,8);
mBinHvl.centers = mean([mBinHvl.edges(1:end-1);mBinHvl.edges(2:end)]);
mBinHvl.count = numel(mBinHvl.edges)-1;
mBinHba.edges = linspace(-1.2,1.2,8);
mBinHba.centers = mean([mBinHba.edges(1:end-1);mBinHba.edges(2:end)]);
mBinHba.count = numel(mBinHba.edges)-1;
mout = zeros([mBinHba.count,mBinHvl.count]);
sout = zeros([mBinHba.count,mBinHvl.count]);
tout = zeros([mBinHba.count,mBinHvl.count]);
cout = zeros([mBinHba.count,mBinHvl.count]);
for a = 1:mBinHba.count
    for v = 1:mBinHvl.count
        ind =   WithinRanges(decoded.phz,phzBin.edges([3,4])) ...
              & WithinRanges(decoded.hav,mBinHvl.edges([v,v+1])) ...
              & WithinRanges(decoded.hba,mBinHba.edges([a,a+1])) ...              
              & randn(size(decoded.hba))>0;
        mout(a,v) = mean(decoded.lat(ind),'omitnan');
        sout(a,v) = std(decoded.lat(ind),'omitnan');
        [~,tout(a,v)] = ttest(decoded.lat(ind));
        cout(a,v) = sum(ind);
    end
end

mask = ones([mBinHba.count,mBinHvl.count]);
mask(cout<50) = nan;
%mask(tout>abs(norminv(1-(1-0.05)^(1/sum(~isnan(mask(:))))))) = nan;
mask(tout>(1-(1-0.05)^(1/sum(~isnan(mask(:)))))) = nan;
%mask(cout<100) = nan;

figure();
%imagescnan({mBinHba.centers,mBinHvl.centers,(mout.*mask)'},[-80,80],'linear',true,'colorMap',@jet);
subplot(211);
imagescnan({mBinHba.centers,mBinHvl.centers,(mout.*mask)'},[-60,60],'linear',true,'colorMap',@jet);
axis('xy');
subplot(212);
imagescnan({mBinHba.centers,mBinHvl.centers,(sout.*mask)'},[50,150],'linear',true,'colorMap',@jet);
axis('xy');



figure,
hist2([decoded.fwd,decoded.hba],linspace(-300,300,24),linspace(-1.2,1.2,24),'xprob')
colormap('jet');
Lines([],0,'w');
Lines(0,[],'w');


[hfig,fig,fax,sax] = set_figure_layout(figure(666003),'A4','portrait',[],2,2,0.2,0.2);
xBlockOffset = 0;
yBlockOffset = 0;

clims = [100,2000];
[yind, yOffSet, xind, xOffSet] = deal( 1, 0, 1, 0);        
sax(end+1) = axes('Units','centimeters',                                                        ...
                  'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,                         ...
                              fig.page.ypos(yind+yBlockOffset)+yOffSet,                         ...
                              fig.subplot.width,                                                ...
                              fig.subplot.height],                                              ...
                  'FontSize', 8,                                                                ...
                  'LineWidth',1);
hold(sax(end),'on');


nx = 3;
ny = 2;
clims = [100,2000];
clims = [10,500];
figure
sax = gobjects([0,1]);
sax(end+1) = subplot2(ny,nx,1,1);hold('on');
ind = decoded.hba>0.2& randn(size(decoded.hba))>0;
%ind = decoded.hba>0.2&WithinRanges(decoded.phz,phzBin.edges([3,4]))& randn(size(decoded.hba))>0;
hist2([decoded.lat(ind),decoded.fwd(ind)],linspace(-250,250,20),linspace(-200,300,20));
plot(mean(decoded.lat(ind)),mean(decoded.fwd(ind)),'*k');
sax(end).ColorScale = 'log';
axis     (sax(end),'tight');
colormap (sax(end),'jet');
caxis    (sax(end), clims);
Lines(0,[],'w');
Lines([],0,'w');
sax(end).XTickLabel=[];
sax(end).YTickLabel=[];        
d = 1;
subject = struct(rat);
subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
subject = update_subject_patch(subject,'body', hbaBin.count+1-d,  true,hbaBin.edges,hbaBin.centers);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);    

sax(end+1) = subplot2(ny,nx,1,2);hold('on');
hold('on');
%ind = abs(decoded.hba)<0.2&WithinRanges(decoded.phz,phzBin.edges([3,4]))& randn(size(decoded.hba))>0;
ind = abs(decoded.hba)<0.2&randn(size(decoded.hba))>0;
hist2([decoded.lat(ind),decoded.fwd(ind)],linspace(-250,250,20),linspace(-200,300,20));
plot(mean(decoded.lat(ind)),mean(decoded.fwd(ind)),'*k');
set(gca(),'ColorScale','log');     axis('tight');
colormap(gca(),'jet');             caxis(clims);
Lines(0,[],'w');                   Lines([],0,'w');
set(gca(),'XTickLabel',[]);        set(gca(),'YTickLabel',[]);
d = 2;
subject = struct(rat);
subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
subject = update_subject_patch(subject,'body', hbaBin.count+1-d,  true,hbaBin.edges,hbaBin.centers);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);    

sax(end+1) = subplot2(ny,nx,1,3);hold('on');
%ind = decoded.hba<-0.2 &WithinRanges(decoded.phz,phzBin.edges([3,4]))& randn(size(decoded.hba))>0;
ind = decoded.hba<-0.2 &randn(size(decoded.hba))>0;
hist2([decoded.lat(ind),decoded.fwd(ind)],linspace(-250,250,20),linspace(-200,300,20));
plot(mean(decoded.lat(ind)),mean(decoded.fwd(ind)),'*k');
set(gca(),'ColorScale','log');     axis('tight');
colormap(gca(),'jet');             caxis(clims);
Lines(0,[],'w');                   Lines([],0,'w');
set(gca(),'XTickLabel',[]);        set(gca(),'YTickLabel',[]);
d = 3;
subject = struct(rat);
subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
subject = update_subject_patch(subject,'body', hbaBin.count+1-d,  true,hbaBin.edges,hbaBin.centers);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);    


sax = gobjects([0,1]);
sax(end+1) = subplot2(ny,nx,2,1);hold('on');
ind = decoded.hba>0.2& randn(size(decoded.hba))>0;
%ind = decoded.hba>0.2&WithinRanges(decoded.phz,phzBin.edges([3,4]))& randn(size(decoded.hba))>0;
hist2([decoded.tlat(ind),decoded.tfwd(ind)],linspace(-250,250,20),linspace(-200,300,20));
plot(mean(decoded.tlat(ind)),mean(decoded.tfwd(ind)),'*k');
sax(end).ColorScale = 'log';
axis     (sax(end),'tight');
colormap (sax(end),'jet');
caxis    (sax(end), clims);
Lines(0,[],'w');
Lines([],0,'w');
sax(end).XTickLabel=[];
sax(end).YTickLabel=[];        
d = 1;
subject = struct(rat);
subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
subject = update_subject_patch(subject,'body', hbaBin.count+1-d,  true,hbaBin.edges,hbaBin.centers);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);    

sax(end+1) = subplot2(ny,nx,2,2);hold('on');
hold('on');
%ind = abs(decoded.hba)<0.2&WithinRanges(decoded.phz,phzBin.edges([3,4]))& randn(size(decoded.hba))>0;
ind = abs(decoded.hba)<0.2&randn(size(decoded.hba))>0;
hist2([decoded.tlat(ind),decoded.tfwd(ind)],linspace(-250,250,20),linspace(-200,300,20));
plot(mean(decoded.tlat(ind)),mean(decoded.tfwd(ind)),'*k');
set(gca(),'ColorScale','log');     axis('tight');
colormap(gca(),'jet');             caxis(clims);
Lines(0,[],'w');                   Lines([],0,'w');
set(gca(),'XTickLabel',[]);        set(gca(),'YTickLabel',[]);
d = 2;
subject = struct(rat);
subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
subject = update_subject_patch(subject,'body', hbaBin.count+1-d,  true,hbaBin.edges,hbaBin.centers);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);    

sax(end+1) = subplot2(ny,nx,2,3);hold('on');
%ind = decoded.hba<-0.2 &WithinRanges(decoded.phz,phzBin.edges([3,4]))& randn(size(decoded.hba))>0;
ind = decoded.hba<-0.2 &randn(size(decoded.hba))>0;
hist2([decoded.tlat(ind),decoded.tfwd(ind)],linspace(-250,250,20),linspace(-200,300,20));
plot(mean(decoded.tlat(ind)),mean(decoded.tfwd(ind)),'*k');
set(gca(),'ColorScale','log');     axis('tight');
colormap(gca(),'jet');             caxis(clims);
Lines(0,[],'w');                   Lines([],0,'w');
set(gca(),'XTickLabel',[]);        set(gca(),'YTickLabel',[]);
d = 3;
subject = struct(rat);
subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
subject = update_subject_patch(subject,'body', hbaBin.count+1-d,  true,hbaBin.edges,hbaBin.centers);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);    




hbaBnds = {[0.2,1.2],[-0.2,0.2],[-1.2,-0.2]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
havBnds = {[-0.3,-0.018],[-0.018,-0.009],[-0.009,0.009],[0.009,0.018],[0.018,0.3]};
hvlBnds = {[-50,-5],[-5,5],[5,50]};
figure
norm = 'xprob';
xBinEdg = linspace(-300,300,8);
xBinCtr = mean([xBinEdg(1:end-1);xBinEdg(2:end)]);
yBinEdg = linspace(0.5,2*pi-0.5,4);
yBinCtr = mean([yBinEdg(1:end-1);yBinEdg(2:end)]);
clims = [0,0.4];
%clims = 'auto';
out = zeros([7,3,3,3]);
for h = 1:numel(hbaBnds)
    for v = 1:numel(hvlBnds)
        subplot2(numel(hvlBnds),numel(hbaBnds),v,h);
        ind = WithinRanges(decoded.hba,hbaBnds{h}) & ...
              WithinRanges(decoded.hvl,hvlBnds{v}) & ...
              randn(size(decoded.hvl))>0.5;
        out(:,:,h,v) = hist2([decoded.tlat(ind), decoded.phz(ind)],xBinEdg,yBinEdg,norm);
        imagesc(xBinCtr, yBinCtr, out(:,:,h,v)');
        colormap('jet');caxis(clims);
        axis('xy');
    end
end





%%

%% HEAD BODY ANGLE ---------------------------------------------------------------------------------
medD = [];
skwD = [];
stdD = [];
medR = [];
skwR = [];
stdR = [];
medL = [];
skwL = [];
stdL = [];
medC = [];
skwC = [];
stdC = [];

rdists = 0:5:30;
for r = 1:numel(rdists)
    clear('xcomp','ycomp','zcomp','ccomp');
    xcomp.label = 'hba';            xcomp.data = [];    xcomp.edgs = [-1.25,-0.2,0.2,1.25];
    ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5,6 );
    ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-80,80];
    fcomp.data = [];
    mcomp.data = [];
    %for t = 1:10
    %for t = [2,3,4,5,7,8,9]
    for t = [1:3,5:8,11],
        %for t = [6]    
        dc = dca{t};
        
        mang = sq(dca{t}.xyz(:,'nose',[1,2])-dca{t}.xyz(:,'hcom',[1,2]));
        mbang = atan2(mang(:,2),mang(:,1));
        bang =  sq(dca{t}.xyz(:,'hcom',[1,2]));
        bbang = atan2(bang(:,2),bang(:,1));
        mbbang = circ_dist(bbang,mbang);
        
        mind =  dc.stcm(:,1)==1                                      ...
                & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5)  ...               
                & dc.hvfl(:,1)>-2 ...
                & dc.ucnt>=2 & dc.ucnt<10 ...
                & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)<rdists(r)+10 ...
                & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)>rdists(r);
% $$$                        & sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))>200;

        %mind(mind==true) = randn([sum(mind),1])>0.5;
        mind(mind==true) = randn([sum(mind),1])>0;        
        xcomp.data = cat(1, xcomp.data, -dc.hbang(mind,1));
        ycomp.data = cat(1, ycomp.data, dc.phz(mind));
        %ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+0*double(t<5)-12.5*double(t>=5));    
        %ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+25*double(t<4)-12.5*double(t>=4));
        ccomp.data = cat(1, ccomp.data, dc.esax(mind,2));
        fcomp.data = cat(1, fcomp.data, dc.esax(mind,2));
        mcomp.data = cat(1, mcomp.data, mbbang(mind));
    end

    [xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,ccomp);
    zmean   = zcomp.mean;    zmean(zcomp.count<5)  = nan;
    zmedian = zcomp.median;  zmean(zcomp.count<5)  = nan;
    zstd    = zcomp.std;     zstd(zcomp.count<5)   = nan;
    zcount  = zcomp.count;   zcount(zcomp.count<5) = nan;
    sectors = linspace(-pi,pi,11);
    for s = 1:numel(sectors)-1
        bind = mcomp.data<sectors(s+1) & mcomp.data>sectors(s) & ycomp.data>4& ycomp.data<6;
        ind = bind & xcomp.data>0.2;
        medR(r,s) = median(ccomp.data(ind,1)./10);%-rmodel(r,s);
        skwR(r,s) = skew(ccomp.data(ind,1)./10);
        stdR(r,s) = std(ccomp.data(ind,1)./10);    
        ind = bind & xcomp.data<0.2& xcomp.data>-0.2;
        medC(r,s) = median(ccomp.data(ind,1)./10);%-rmodel(r,s);
        skwC(r,s) = skew(ccomp.data(ind,1)./10);
        stdC(r,s) = std(ccomp.data(ind,1)./10);    
        ind = bind & xcomp.data<-0.2;
        medL(r,s) = median(ccomp.data(ind,1)./10);%-rmodel(r,s);
        skwL(r,s) = skew(ccomp.data(ind,1)./10);
        stdL(r,s) = std(ccomp.data(ind,1)./10);    
        medD(r,s) = medR(r,s)-medL(r,s);
        skwD(r,s) = skwR(r,s)-skwL(r,s);
        stdD(r,s) = stdR(r,s)-stdL(r,s);
    end
end

sectorc = mean([sectors(2:end);sectors(1:end-1)]);

% $$$ figure,imagesc(medR')


[THETA,RR] = meshgrid(sectors,[rdists,35]);
% $$$ 
% $$$ THETA = cat(2,THETA(:,end),THETA);
% $$$ THETA = cat(1,THETA(end,:),THETA);
% $$$ 
% $$$ RR = cat(2,RR(:,end),RR);
% $$$ RR = cat(1,RR(end,:),RR);


[A,B] = pol2cart(circ_dist(THETA,diff(sectors([1,2]))),RR);
%[A,B] = pol2cart(THETA,RR);

% SUPFIG - Non-Uniformity of decoded lateral position: subject-arena orientation and distance to maze center
figure
cmedD = medL; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,1);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-7.5,7.5]);
ylabel(cax,'cm');
title({'Leftward','Median'});
cmedD = stdL; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,1);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([0,20]);
ylabel(cax,'Std')
cmedD = skwL; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,1);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-2,2]);
ylabel(cax,'Skew')

cmedD = medC; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,2);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-7.5,7.5]);
ylabel(cax,'cm');
title({'Centered','Median'})
cmedD = stdC; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,2);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([0,20]);
ylabel(cax,'Std');
cmedD = skwC; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,2);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-2,2]);
title('Skew')

cmedD = medR; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,3);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-7.5,7.5])
ylabel(cax,'cm');
title({'Rightward','Median'});
cmedD = stdR; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,3);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([0,20]);
ylabel(cax,'Std');
cmedD = skwR; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,3);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-2,2]);
ylabel(cax,'Skew');
       
cmedD = medD; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,4);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-7.5,7.5]);
title({'Rightward-Leftward','\DeltaMedian'})
ylabel(cax,'cm');
cmedD = stdD; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,4);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([0,5]);
ylabel(cax,'\DeltaStd');
cmedD = skwD; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,4);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-1,1]);
ylabel(cax,'\DeltaSkew');

ForAllSubplots('daspect([1,1,1]);');
ForAllSubplots('ylim([-40,40]);');
ForAllSubplots('xlim([-40,40]);');
fax = axes('Position',[0,0,1,1],'Visible','off');
text(fax,0.1,0.9,{'Non-Uniformity of decoded lateral position ',...
                  'depedent on subject-arena orientation and distance to maze center'});
set(gcf(),'PaperOrientation','landscape');





%% Forward


%% HEAD BODY ANGLE ---------------------------------------------------------------------------------
medD = [];
skwD = [];
stdD = [];
medR = [];
skwR = [];
stdR = [];
medL = [];
skwL = [];
stdL = [];
medC = [];
skwC = [];
stdC = [];

rdists = 0:5:30;
for r = 1:numel(rdists)
    clear('xcomp','ycomp','zcomp','ccomp');
    xcomp.label = 'hba';            xcomp.data = [];    xcomp.edgs = [-1.25,-0.2,0.2,1.25];
    ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5,6 );
    ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-80,80];
    fcomp.data = [];
    mcomp.data = [];
    %for t = 1:10
    %for t = [2,3,4,5,7,8,9]
    for t = [1:3,5:8,11],
        %for t = [6]    
        dc = dca{t};
        
        mang = sq(dca{t}.xyz(:,'nose',[1,2])-dca{t}.xyz(:,'hcom',[1,2]));
        mbang = atan2(mang(:,2),mang(:,1));
        bang =  sq(dca{t}.xyz(:,'hcom',[1,2]));
        bbang = atan2(bang(:,2),bang(:,1));
        mbbang = circ_dist(bbang,mbang);
        
        mind =  dc.stcm(:,1)==1                                      ...
                & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5)  ...               
                & dc.hvfl(:,1)>-2 ...
                & dc.ucnt>=2 & dc.ucnt<10 ...
                & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)<rdists(r)+10 ...
                & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)>rdists(r);
% $$$                        & sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))>200;

        %mind(mind==true) = randn([sum(mind),1])>0.5;
        mind(mind==true) = randn([sum(mind),1])>0;        
        xcomp.data = cat(1, xcomp.data, -dc.hbang(mind,1));
        ycomp.data = cat(1, ycomp.data, dc.phz(mind));
        %ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+0*double(t<5)-12.5*double(t>=5));    
        %ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+25*double(t<4)-12.5*double(t>=4));
        ccomp.data = cat(1, ccomp.data, dc.esax(mind,1));
        fcomp.data = cat(1, fcomp.data, dc.esax(mind,1));
        mcomp.data = cat(1, mcomp.data, mbbang(mind));
    end

    [xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,ccomp);
    zmean   = zcomp.mean;    zmean(zcomp.count<5)  = nan;
    zmedian = zcomp.median;  zmean(zcomp.count<5)  = nan;
    zstd    = zcomp.std;     zstd(zcomp.count<5)   = nan;
    zcount  = zcomp.count;   zcount(zcomp.count<5) = nan;
    sectors = linspace(-pi,pi,11);
    for s = 1:numel(sectors)-1
        bind = mcomp.data<sectors(s+1) & mcomp.data>sectors(s) & ycomp.data>4& ycomp.data<6;
        ind = bind & xcomp.data>0.2;
        medR(r,s) = median(ccomp.data(ind,1)./10);%-rmodel(r,s);
        skwR(r,s) = skew(ccomp.data(ind,1)./10);
        stdR(r,s) = std(ccomp.data(ind,1)./10);    
        ind = bind & xcomp.data<0.2& xcomp.data>-0.2;
        medC(r,s) = median(ccomp.data(ind,1)./10);%-rmodel(r,s);
        skwC(r,s) = skew(ccomp.data(ind,1)./10);
        stdC(r,s) = std(ccomp.data(ind,1)./10);    
        ind = bind & xcomp.data<-0.2;
        medL(r,s) = median(ccomp.data(ind,1)./10);%-rmodel(r,s);
        skwL(r,s) = skew(ccomp.data(ind,1)./10);
        stdL(r,s) = std(ccomp.data(ind,1)./10);    
        medD(r,s) = medR(r,s)-medL(r,s);
        skwD(r,s) = skwR(r,s)-skwL(r,s);
        stdD(r,s) = stdR(r,s)-stdL(r,s);
    end
end

sectorc = mean([sectors(2:end);sectors(1:end-1)]);

% $$$ figure,imagesc(medR')


[THETA,RR] = meshgrid(sectors,[rdists,35]);
% $$$ 
% $$$ THETA = cat(2,THETA(:,end),THETA);
% $$$ THETA = cat(1,THETA(end,:),THETA);
% $$$ 
% $$$ RR = cat(2,RR(:,end),RR);
% $$$ RR = cat(1,RR(end,:),RR);


[A,B] = pol2cart(circ_dist(THETA,diff(sectors([1,2]))),RR);
%[A,B] = pol2cart(THETA,RR);

% SUPFIG - Non-Uniformity of decoded anterioposterior position:
% subject-arena orientation and distance to maze center
opt.view = [0,90];
opt.clim = [-5,15];
figure
cmedD = medL; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,1);surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis(opt.clim);
ylabel(cax,'cm');
title({'Leftward','Median'});
cmedD = stdL; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,1);surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis([0,20]);
ylabel(cax,'Std')
cmedD = skwL; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,1);surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis([-2,2]);
ylabel(cax,'Skew')

cmedD = medC; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,2);surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis(opt.clim);
ylabel(cax,'cm');
title({'Centered','Median'})
cmedD = stdC; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,2);surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis([0,20]);
ylabel(cax,'Std');
cmedD = skwC; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,2);surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis([-2,2]);
title('Skew')

cmedD = medR; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,3);surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis(opt.clim)
ylabel(cax,'cm');
title({'Rightward','Median'});
cmedD = stdR; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,3);surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis([0,20]);
ylabel(cax,'Std');
cmedD = skwR; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,3);surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis([-2,2]);
ylabel(cax,'Skew');
       
cmedD = medD; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,4);surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis([-5,5]);
title({'Rightward-Leftward','\DeltaMedian'})
ylabel(cax,'cm');
cmedD = stdD; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,4);surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis([0,5]);
ylabel(cax,'\DeltaStd');
cmedD = skwD; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,4);surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis([-1,1]);
ylabel(cax,'\DeltaSkew');

ForAllSubplots('daspect([1,1,1]);');
ForAllSubplots('ylim([-40,40]);');
ForAllSubplots('xlim([-40,40]);');
fax = axes('Position',[0,0,1,1],'Visible','off');
text(fax,0.1,0.9,{'Non-Uniformity of decoded anterioposterior position ',...
                  'depedent on subject-arena orientation and distance to maze center'});
set(gcf(),'PaperOrientation','landscape');

