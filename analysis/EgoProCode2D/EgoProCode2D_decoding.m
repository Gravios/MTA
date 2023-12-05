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
halfSpikeWindow = 0.020;

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


[hfig,fig,fax,sax] = set_figure_layout(figure(667001),'A4','portrait',[],2.5,2.5,0.6,0.6);
% GENERATE Copula between LAT and HBA for states W+P+N 
for t = 1:numel(dca)
    ind = WithinRanges(dca{t}.phz,[4,6])                                  ...
          & dca{t}.stcm(:,1)==1                                           ...
          ... & dca{t}.ucnt>=2 & dca{t}.ucnt<=5 ...
          & (dca{t}.stcm(:,3)==3|dca{t}.stcm(:,4)==4|dca{t}.stcm(:,5)==5) ...
          & dca{t}.hdist.data<00;
    sax(end+1) = subplot2(5,3,ceil(t/3)+1,mod(t,3)+1);
    hist2([round(dca{t}.ecom(ind,2)./10,2),-dca{t}.hbang(ind,1)],11,11,'','mud');
    colormap(gca(),'jet');
    title({Trials{tind(t)}.name});
    xlabel('cm');
    ylabel('rad');
end
pause(0.5);
fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
xlim([0,hfig.Position(3)]);
ylim([0,hfig.Position(4)]);
pause(0.5);
text(0.5 * fig.page.width,                                                 ...
     0.9 * fig.page.height,                                                ...
     {'Copula Between Lateral Decoded Position in the Head Frame of Reference VS the Head Body Angle',...
      'for Individual Sessions'},...
     'HorizontalAlignment','Center'                                        ...
);





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


% ACCUMULATE vars for specific subset of decoding
s = 5;
decoded = struct('fwd',[],...
                 'lat',[],...
                 'tfwd',[],...
                 'tlat',[],...
                 'xyz',[],...
                 'dst',[],...
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
        &  dca{t}.hvfl(:,1)>0                                            ...
        & dca{t}.ucnt>=4 & dca{t}.ucnt<=10                                 ...
        & sqrt(sum(dca{t}.xyz(:,'hcom',[1,2]).^2,3))>0                    ...
        & sqrt(sum(dca{t}.xyz(:,'hcom',[1,2]).^2,3))<500;
    decoded.fwd = cat(1,decoded.fwd,dca{t}.ecom(mind,1));
    decoded.lat = cat(1,decoded.lat,dca{t}.ecom(mind,2));%-10*double(t>4)+10*double(t<4));
    decoded.tfwd = cat(1,decoded.tfwd,dca{t}.tcom(mind,1));
    decoded.tlat = cat(1,decoded.tlat,dca{t}.tcom(mind,2));%+20*double(t>4));
    decoded.xyz =  cat(1,decoded.xyz, dca{t}.xyz(mind,{'hcom','nose'},:));
    decoded.dst =  cat(1,decoded.dst, dca{t}.hdist(mind));
    decoded.hvf = cat(1,decoded.hvf,dca{t}.hvfl(mind,1));
    decoded.hvl = cat(1,decoded.hvl,dca{t}.hvfl(mind,2));
    decoded.bvl = cat(1,decoded.bvl,dca{t}.bvfl(mind,2));    
    decoded.hvlF = cat(1,decoded.hvlF,dca{t}.hvflF(mind,2));
    decoded.hvlP = cat(1,decoded.hvlP,dca{t}.hvflP(mind,2));        
    decoded.hav = cat(1,decoded.hav,dca{t}.hvang(mind,1));
    decoded.hba = cat(1,decoded.hba,dca{t}.hbang(mind,1));
    decoded.phz = cat(1,decoded.phz,dca{t}.phz(mind,1));
end
headAngle = sq(decoded.xyz(:,2,[1,2])-decoded.xyz(:,1,[1,2]));
headAngle = atan2(headAngle(:,2),headAngle(:,1));
mazeAngle =  sq(decoded.xyz(:,1,[1,2]));
mazeAngle = atan2(mazeAngle(:,2),mazeAngle(:,1));
decoded.hma = circ_dist(headAngle,mazeAngle);
 
[Xe,Ye] = pol2cart(decoded.hma,decoded.dst);
decoded.correction.lat = mfun(beta,[Xe(:),Ye(:)]);
decoded.clat = (decoded.lat-decoded.correction.lat)-8;
decoded.clt = (decoded.lat-decoded.correction.lat)-8;


ind = WithinRanges(decoded.phz,[4.5,5.5]) ... 
     & randn(size(decoded.hba))>1 ...
      & abs(decoded.hba)<1.2;
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind)]);STATS
[B,BINT,R,RINT,STATS] = regress(decoded.clt(ind),[ones([sum(ind),1]),decoded.hba(ind)]);STATS

ind = WithinRanges(decoded.phz,[1.5,2.5]) ... 
     & randn(size(decoded.hba))>0 ...
      & abs(decoded.hba)<0.4;

figure,
hist2([decoded.fwd(ind),decoded.hvf(ind)],linspace(-400,400,30),0: ...
      10:60,'','mud');


figure,
hist2([decoded.fwd(ind),decoded.hvf(ind)],linspace(-400,400,30),-5: ...
      15:65,'','mud');


ind = WithinRanges(decoded.phz,[4.5,5.5]) ... 
      & randn(size(decoded.hba))>0.5 ...
      & abs(decoded.hba)<1.2;
figure,
hist2([decoded.clt(ind),decoded.hba(ind)],linspace(-400,400,8),-1.2:0.2:1.2,'','mud');

figure,
hist2([decoded.clt(ind),decoded.hba(ind)],linspace(-400,400,8),-1.2:0.2:1.2,'xprob','');
figure,plot(decoded.clt(ind),decoded.hba(ind),'.')


figure,
hist2([decoded.lat(ind),decoded.hav(ind)],linspace(-400,400,8),-.5:0.1:.5,'','mud');

figure
hist2([decoded.lat(ind),decoded.hav(ind)],linspace(-400,400,8),-.5:0.1:.5,'xprob','');
colormap('jet');


figure
hist2([decoded.clt(ind),decoded.hav(ind)],linspace(-400,400,8),-.5:0.1:.5,'xprob','');
colormap('jet');

ind = WithinRanges(decoded.phz,[5,6]) ... 
      & randn(size(decoded.hba))>0.5 ...
      & decoded.hba<-0.2;      
figure,
hist2([decoded.lat(ind),decoded.hvl(ind)],linspace(-400,400,30),-40:10:40,'xprob','');
colormap('jet');

figure,
hold('on');
for hbaInd = 1:hbaBin.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) ... 
      & randn(size(decoded.hba))>0 ...
      & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1]));      
%set(histogram(decoded.lat(ind),linspace(-300,300,40)),'EdgeColor','none','FaceAlpha',.3); 
[F,xi] = ksdensity(decoded.lat(ind)/10);
plot(xi,F,'-','Color',hbaBin.color(hbaInd,:));
end
xlim([-40,40]);



figure();
hold('on');
histogram(decoded.fwd(ind&decoded.hvf<5 & decoded.hvf>-5 & abs(decoded.hvl<15)),linspace(-400,400,30),'Normalization','probability')
histogram(decoded.fwd(ind&decoded.hvf<50& decoded.hvf>40 & abs(decoded.hvl<15)),linspace(-400,400,30),'Normalization','probability')


ind = WithinRanges(decoded.phz,[1,2.5]) ... 
     & randn(size(decoded.hba))>0 ...
      & abs(decoded.hba)<0.4;

medf = [...
median(decoded.fwd(ind&decoded.hvf<5 & decoded.hvf>-5)),...
median(decoded.fwd(ind&decoded.hvf<15 & decoded.hvf>5)),...
median(decoded.fwd(ind&decoded.hvf<25 & decoded.hvf>15)),...
median(decoded.fwd(ind&decoded.hvf<35 & decoded.hvf>25)),...
median(decoded.fwd(ind&decoded.hvf<45 & decoded.hvf>35)),...
median(decoded.fwd(ind&decoded.hvf<55& decoded.hvf>45))];
figure,plot(medf)

[h,p] = ttest2(decoded.fwd(ind&decoded.hvf<50&decoded.hvf>40),decoded.fwd(ind&decoded.hvf<20&decoded.hvf>10))


figure,
hold('on');
for hvfInd = 1:hvfBin.count-1
ind = WithinRanges(decoded.phz,[1.5,2.5]) ... 
      & randn(size(decoded.hba))>0 ...
      & WithinRanges(decoded.hvf,hvfBin.edges(hvfInd+[1,2]));      
[F,xi] = ksdensity(decoded.fwd(ind)/10);
plot(xi,F,'-','Color',hbaBin.color(hvfInd,:));
end
xlim(gca(),[-30,30]);

dhva = discretize(decoded.hva(ind),-5:10:55);
dfwd = decoded.fwd(ind);
nind = nniz([dfwd,dhvf]);
out = accumarray(dhvf(nind),dfwd(nind)-20,[6,1],@mean);
figure,plot(out)


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
mBinHvl.edges = linspace(-45,45,8);
%mBinHvl.edges = linspace(-0.5,0.5,8);
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
              & WithinRanges(decoded.hvl,mBinHvl.edges([v,v+1])) ...
              & WithinRanges(decoded.hba,mBinHba.edges([a,a+1])) ...              
              & randn(size(decoded.hba))>0;
        mout(a,v) = mean(decoded.clat(ind),'omitnan');
        sout(a,v) = std(decoded.clat(ind),'omitnan');
        [~,tout(a,v)] = ttest(decoded.clat(ind));
        cout(a,v) = sum(ind);
    end
end
mask = ones([mBinHba.count,mBinHvl.count]);
mask(cout<100) = nan;
%mask(tout>abs(norminv(1-(1-0.05)^(1/sum(~isnan(mask(:))))))) = nan;
%mask(tout>(1-(1-0.1)^(1/sum(~isnan(mask(:)))))) = nan;
%mask(cout<100) = nan;
figure();
%imagescnan({mBinHba.centers,mBinHvl.centers,(mout.*mask)'},[-80,80],'linear',true,'colorMap',@jet);
subplot(211);
imagescnan({mBinHba.centers,mBinHvl.centers,(mout.*mask)'},[-6,6],'linear',true,'colorMap',@jet);
axis('xy');
subplot(212);
imagescnan({mBinHba.centers,mBinHvl.centers,(sout.*mask)'},[5,15],'linear',true,'colorMap',@jet);
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

rdists = 0:5:35;
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


[THETA,RR] = meshgrid(sectors,[rdists,rdists(end)]);


%[A,B] = pol2cart(circ_dist(THETA,diff(sectors([1,2]))),RR);
[X,Y] = pol2cart(THETA,RR);

figure()
hold('on');
for t = 1:numel(sectors)-1
    for r = 1:numel(rdists)
        patch([X(r,t),X(r,t+1),X(r+1,t+1),X(r+1,t)],...
              [Y(r,t),Y(r,t+1),Y(r+1,t+1),Y(r+1,t)],...
              medC(r,t));
    end
end
colormap(gca,'jet')
caxis([-10,10])




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





%% Ego Forward
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

rdists = 5:5:45;
for r = 1:numel(rdists)
    clear('xcomp','ycomp','zcomp','ccomp');
    xcomp.label = 'hba';            xcomp.data = [];    xcomp.edgs = [-1.25,-0.2,0.2,1.25];
    ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5,6 );
    ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-80,80];
    fcomp.data = [];
    mcomp.data = [];
    for t = [1:3,5:8,11],
        %for t = [6]    
        dc = dca{t};
        
        headAngle = sq(dca{t}.xyz(:,'nose',[1,2])-dca{t}.xyz(:,'hcom',[1,2]));
        headAngle = atan2(headAngle(:,2),headAngle(:,1));
        
        mazeAngle =  sq(dca{t}.xyz(:,'hcom',[1,2]));
        mazeAngle = atan2(mazeAngle(:,2),mazeAngle(:,1));
        
        headMazeAngle = circ_dist(headAngle, mazeAngle);
        
        mind =  dc.stcm(:,1)==1                                             ... Theta
                & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5)         ... Walk,Turn,Pause              
                & dc.hvfl(:,1)>-2                                           ... FwdHeadSpeed
                & dc.ucnt>=4 & dc.ucnt<=10                                    ... UnitCount
                & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)<rdists(r)+5 ... DistanceFromMazeCenter
                & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)>rdists(r)-5;

        %mind(mind==true) = randn([sum(mind),1])>0; % Subsample
        
        xcomp.data = cat(1, xcomp.data, dc.hbang(mind,1));
        ycomp.data = cat(1, ycomp.data, dc.phz(mind));
        
        ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)./10);
        %ccomp.data = cat(1, ccomp.data, dc.esax(mind,1)./10);        
        mcomp.data = cat(1, mcomp.data, headMazeAngle(mind));
    end
    
    [xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,ccomp);
    zmean   = zcomp.mean;    zmean(zcomp.count<5)  = nan;
    zmedian = zcomp.median;  zmean(zcomp.count<5)  = nan;
    zstd    = zcomp.std;     zstd(zcomp.count<5)   = nan;
    zcount  = zcomp.count;   zcount(zcomp.count<5) = nan;
    sectors = linspace(-pi,pi,11);
    
    sectorsL = linspace(-pi,pi,11);
    for s = 1:numel(sectors)-1
        if r>1
            bind = mcomp.data<sectors(s+1) & mcomp.data>sectors(s) ...
                   & ycomp.data>4.5& ycomp.data<5;
        else
            bind = ycomp.data>4.5& ycomp.data<5;
        end
        

        ind  = bind & xcomp.data<-0.2;        
        medR(r,s) = median(ccomp.data(ind,1));%-rmodel(r,s);
        skwR(r,s) = skew(ccomp.data(ind,1));
        stdR(r,s) = std(ccomp.data(ind,1));
        
        ind = bind & xcomp.data<0.2& xcomp.data>-0.2;
        medC(r,s) = median(ccomp.data(ind,1));%-rmodel(r,s);
        skwC(r,s) = skew(ccomp.data(ind,1));
        stdC(r,s) = std(ccomp.data(ind,1));
        
        ind = bind & xcomp.data>0.2;
        medL(r,s) = median(ccomp.data(ind,1));%-rmodel(r,s);
        skwL(r,s) = skew(ccomp.data(ind,1));
        stdL(r,s) = std(ccomp.data(ind,1));
        
        medD(r,s) = medR(r,s)-medL(r,s);
        skwD(r,s) = skwR(r,s)-skwL(r,s);
        stdD(r,s) = stdR(r,s)-stdL(r,s);
    end
end

sectorc = mean([sectors(2:end);sectors(1:end-1)]);
rdistc = mean([rdists(2:end);rdists(1:end-1)]);

rdiste = [0,rdists+2.5];

%[THETA,RR] = meshgrid(sectors,[rdists,rdists(end)]);
[THETA,RR] = meshgrid(sectors,[rdiste]);
[THETAC,RRC] = meshgrid(sectorc,rdists);

%[A,B] = pol2cart(circ_dist(THETA,diff(sectors([1,2]))),RR);
[X,Y] = pol2cart(THETA,RR);
[Xc,Yc] = pol2cart(THETAC,RRC);

figure();
subplot(131);
hold('on');
for t = 1:numel(sectors)-1 
    for r = 1:numel(rdists)-1
        patch([X(r,t),X(r,t+1),X(r+1,t+1),X(r+1,t)],...
              [Y(r,t),Y(r,t+1),Y(r+1,t+1),Y(r+1,t)],...
              medL(r,t));
    end
end
colormap(gca,'jet')
caxis([-15,15])

subplot(132);
hold('on');
for t = 1:numel(sectors)-1
    for r = 1:numel(rdists)-1
        patch([X(r,t),X(r,t+1),X(r+1,t+1),X(r+1,t)],...
              [Y(r,t),Y(r,t+1),Y(r+1,t+1),Y(r+1,t)],...
              medC(r,t));
    end
end
colormap(gca,'jet')
caxis([-15,15])

subplot(133);
hold('on');
for t = 1:numel(sectors)-1
    for r = 1:numel(rdists)-1
        patch([X(r,t),X(r,t+1),X(r+1,t+1),X(r+1,t)],...
              [Y(r,t),Y(r,t+1),Y(r+1,t+1),Y(r+1,t)],...
              medR(r,t));
    end
end
colormap(gca,'jet')
caxis([-15,15])


subplot(133);
hold('on');
for t = 1:numel(sectors)-1
    for r = 1:numel(rdists)-1
        patch([X(r,t),X(r,t+1),X(r+1,t+1),X(r+1,t)],...
              [Y(r,t),Y(r,t+1),Y(r+1,t+1),Y(r+1,t)],...
              medR(r,t)-medL(r,t));
    end
end
colormap(gca,'jet')
caxis([-10,10])




figure();
subplot(131);
hold('on');
patch(bsxfun(@plus,repmat(X(:),[1,4]),[-dx,dx,dx,-dx]),...
      bsxfun(@plus,repmat(Y(:),[1,4]),[-dy,dy,dy,-dy]),...
      medL(:));
colormap(gca,'jet')
caxis([-15,15])

nx =3;
ny =1;
sax = gobjects([0,1]);
figure,
sax(end+1) = subplot2(ny,nx,1,1);
hold('on');
for t = 1:numel(sectors)-1
    for r = 1:numel(rdists)-1
        plot3(mean(mean(X(r:r+1,t:t+1))),...
              mean(mean(Y(r:r+1,t:t+1))),medC(r,t),'o')
    end
end
daspect(gca(), [1,1,0.3])
grid   (gca(), 'on');
xlim   (gca(), [-50,50])
ylim   (gca(), [-50,50])
% trajectory values of headMazeAngle and headMazeDistance -> pol2cart




sax(end+1) = subplot2(ny,nx,1,2);
hold('on');
for t = 1:numel(sectors)-1
    for r = 1:numel(rdists)-1
        plot3(mean(mean(X(r:r+1,t:t+1))),...
              mean(mean(Y(r:r+1,t:t+1))),...
              stdC(r,t),'o')
    end
end
daspect(gca(), [1,1,0.1])
grid   (gca(), 'on');
xlim   (gca(), [-50,50])
ylim   (gca(), [-50,50])
zlim   (gca(), [11,16])
% trajectory values of headMazeAngle and headMazeDistance -> pol2cart


sax(end+1) = subplot2(ny,nx,1,3);
hold('on');
for t = 1:numel(sectors)-1
    for r = 1:numel(rdists)-1
        plot3(mean(mean(X(r:r+1,t:t+1))),...
              mean(mean(Y(r:r+1,t:t+1))),...
              skwC(r,t),'o')
    end
end
daspect(gca(), [1,1,0.1])
grid   (gca(), 'on');
xlim   (gca(), [-50,50])
ylim   (gca(), [-50,50])
% trajectory values of headMazeAngle and headMazeDistance -> pol2cart
linkprop(sax,'View');



figure()
hold('on');
plot(Yc(:),medC(:),'ob')
plot(Yc(:),medR(:),'or')
plot(Yc(:),medL(:),'og')




figure()
hold('on');
plot(Yc(:),medC(:)-medR(:),'or')
plot(Yc(:),medC(:)-medL(:),'og')

figure,
histogram(medR(:)-medL(:),linspace(-10,10,20));









var = 'hba';
clear('decLat');
decLat.vars = {bins.hba};
Ya = repmat(Yc(:),[3,1]);
medA = cat(1,medL(:),medC(:),medR(:));
ind =~isnan(Ya(:))&~isnan(medA(:));
decLat.fit.p = polyfit(Ya(ind),medA(ind),1);

figure()
hold('on');
plot(Yc(:),medC(:),'ob');
plot(Yc(:),medR(:),'or');
plot(Yc(:),medL(:),'og');
xlim([-50,50]);
ylim([-30,30]);
line([min(Yc(:)),max(Yc(:))],polyval(decLat.fit.p,xlim),'Color','k');


[~,decoded.correction.lat] = pol2cart(decoded.hma,decoded.dst);
decoded.correction.lat = polyval(p,decoded.correction.lat);





mfun = @(beta,x) beta(2).*(sum(x.^2,2)).*sin(atan2(x(:,2),x(:,1)))+beta(1);
[beta,R,J,COVB,MSE,EMI] = nlinfit([Xc(:),Yc(:)].*10,medC(:).*10,mfun,[0,0.0001]);

%figure,plot3(Xc(:),Yc(:),mfun(beta,[Xc(:),Yc(:)]),'o');

figure,
hold('on');
plot3(Xc(:),Yc(:),medC(:),'ok');
plot3(Xc(:),Yc(:),R(:),'or');

[~,decoded.correction.lat] = pol2cart(decoded.hma,decoded.dst);
decoded.correction.lat = polyval(decLat.fit.p,decoded.correction.lat);

[Xe,Ye] = pol2cart(decoded.hma,decoded.dst);
decoded.correction.lat = mfun(beta,[Xe(:),Ye(:)]);

decoded.clat = (decoded.lat-decoded.correction.lat)-8;
decoded.clt = (decoded.lat-decoded.correction.lat)-8;


decLat.hba.raw.cntr.mean        = [];
decLat.hba.raw.cntr.median      = [];
decLat.hba.raw.cntr.std         = [];
decLat.hba.raw.cntr.ttest.h     = [];
decLat.hba.raw.cntr.ttest.pval  = [];
decLat.hba.raw.cntr.ttest.ci    = [];
decLat.hba.raw.cntr.ttest.stats = [];

decLat.hba.raw.peri.mean        = [];
decLat.hba.raw.peri.median      = [];
decLat.hba.raw.peri.std         = [];
decLat.hba.raw.peri.ttest.h     = [];
decLat.hba.raw.peri.ttest.pval  = [];
decLat.hba.raw.peri.ttest.ci    = [];
decLat.hba.raw.peri.ttest.stats = [];

decLat.hba.res.cntr.mean        = [];
decLat.hba.res.cntr.median      = [];
decLat.hba.res.cntr.std         = [];
decLat.hba.res.cntr.ttest.h     = [];
decLat.hba.res.cntr.ttest.pval  = [];
decLat.hba.res.cntr.ttest.ci    = [];
decLat.hba.res.cntr.ttest.stats = [];

decLat.hba.res.peri.mean        = [];
decLat.hba.res.peri.median      = [];
decLat.hba.res.peri.std         = [];
decLat.hba.res.peri.ttest.h     = [];
decLat.hba.res.peri.ttest.pval  = [];
decLat.hba.res.peri.ttest.ci    = [];
decLat.hba.res.peri.ttest.stats = [];


decLat.hba.res.all.mean        = [];
decLat.hba.res.all.median      = [];
decLat.hba.res.all.std         = [];
decLat.hba.res.all.ttest.h     = [];
decLat.hba.res.all.ttest.pval  = [];
decLat.hba.res.all.ttest.ci    = [];
decLat.hba.res.all.ttest.stats = [];



figure,
subplot(222);
hold('on');
for hbaInd = 1:bins.hba.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) ... 
      & randn(size(decoded.hba))>0.0 ...
      & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1])) ...
      & WithinRanges(decoded.dst,[0,200]);      
[F,xi] = ksdensity((decoded.lat(ind)-decoded.correction.lat(ind))/10);
plot(xi,F,'-','Color',hbaBin.color(hbaInd,:));
% STATS 
decLat.hba(hbaInd).res.cntr.mean = mean((decoded.lat(ind))/10);
decLat.hba(hbaInd).res.cntr.median = median((decoded.lat(ind))/10);
decLat.hba(hbaInd).res.cntr.std  = std( (decoded.lat(ind))/10);
[decLat.hba(hbaInd).res.cntr.ttest.h,...
 decLat.hba(hbaInd).res.cntr.ttest.pval,...
 decLat.hba(hbaInd).res.cntr.ttest.ci,...
 decLat.hba(hbaInd).res.cntr.ttest.stats] = ttest((decoded.lat(ind)-decoded.correction.lat(ind))/10);

% $$$ [~,xi] = NearestNeighbour(dxi,decLat.hba(hbaInd).res.cntr.median);
% $$$ line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',bins.hba.color(hbaInd,:));

end
xlim([-40,40]);
title('Corrected: head-maze-distance < 20cm')

subplot(224)
hold('on');
for hbaInd = 1:bins.hba.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) ... 
      & randn(size(decoded.hba))>0.5 ...
      & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1])) ...
      & WithinRanges(decoded.dst,[200,500]);      
[F,xi] = ksdensity((decoded.lat(ind)-decoded.correction.lat(ind))/10);
plot(xi,F,'-','Color',hbaBin.color(hbaInd,:));
% STATS 
decLat.hba(hbaInd).res.peri.mean = mean((decoded.lat(ind))/10);
decLat.hba(hbaInd).res.peri.std  = std( (decoded.lat(ind))/10);
[decLat.hba(hbaInd).res.peri.ttest.h,...
 decLat.hba(hbaInd).res.peri.ttest.pval,...
 decLat.hba(hbaInd).res.peri.ttest.ci,...
 decLat.hba(hbaInd).res.peri.ttest.stats] = ttest((decoded.lat(ind)-decoded.correction.lat(ind))/10);
end
xlim([-40,40]);
title('Corrected: head-maze-distance > 20cm')

subplot(221);
hold('on');
for hbaInd = 1:bins.hba.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) ... 
      & randn(size(decoded.hba))>0.5 ...
      & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1])) ...
      & WithinRanges(decoded.dst,[0,200]);      
%set(histogram(decoded.lat(ind),linspace(-300,300,40)),'EdgeColor','none','FaceAlpha',.3); 
[F,xi] = ksdensity((decoded.lat(ind))/10);
plot(xi,F,'-','Color',hbaBin.color(hbaInd,:));
% STATS 
decLat.hba(hbaInd).raw.cntr.mean = mean((decoded.lat(ind))/10);
decLat.hba(hbaInd).raw.cntr.std  = std( (decoded.lat(ind))/10);
[decLat.hba(hbaInd).raw.cntr.ttest.h,...
 decLat.hba(hbaInd).raw.cntr.ttest.pval,...
 decLat.hba(hbaInd).raw.cntr.ttest.ci,...
 decLat.hba(hbaInd).raw.cntr.ttest.stats] = ttest((decoded.lat(ind)-decoded.correction.lat(ind))/10);
end
xlim([-40,40]);
title('Uncorrected: head-maze-distance < 20cm')

subplot(223);
hold('on');
for hbaInd = 1:bins.hba.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) ... 
      & randn(size(decoded.hba))>.5 ...
      & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1])) ...
      & WithinRanges(decoded.dst,[200,500]);      
%set(histogram(decoded.lat(ind),linspace(-300,300,40)),'EdgeColor','none','FaceAlpha',.3); 
[F,xi] = ksdensity((decoded.lat(ind))/10);
plot(xi,F,'-','Color',hbaBin.color(hbaInd,:));
% STATS 
decLat.hba(hbaInd).raw.peri.mean = mean((decoded.lat(ind))/10);
decLat.hba(hbaInd).raw.peri.std  = std( (decoded.lat(ind))/10);
[decLat.hba(hbaInd).raw.peri.ttest.h,...
 decLat.hba(hbaInd).raw.peri.ttest.pval,...
 decLat.hba(hbaInd).raw.peri.ttest.ci,...
 decLat.hba(hbaInd).raw.peri.ttest.stats] = ttest((decoded.lat(ind))/10);
end
xlim([-40,40]);
title('Uncorrected: head-maze-distance > 20cm')




figure()
subplot(221);hold('on');
for hbaInd = 1:bins.hba.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) & randn(size(decoded.hba))>.5 & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1])) & WithinRanges(decoded.dst,[0,200]);      
[F] = cdfplot((decoded.lat(ind))/10);  F.Color = hbaBin.color(hbaInd,:);
end, xlim([-40,40]); title('Uncorrected: head-maze-distance > 20cm');

subplot(222); hold('on');
for hbaInd = 1:bins.hba.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) & randn(size(decoded.hba))>.5 & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1])) & WithinRanges(decoded.dst,[0,200]);      
[F] = cdfplot((decoded.lat(ind)-decoded.correction.lat(ind))/10);  F.Color = hbaBin.color(hbaInd,:);
end, xlim([-20,20]); title('corrected: head-maze-distance > 20cm');

subplot(223); hold('on');
for hbaInd = 1:bins.hba.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) & randn(size(decoded.hba))>.5 & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1])) & WithinRanges(decoded.dst,[200,500]);      
[F] = cdfplot((decoded.lat(ind))/10);  F.Color = hbaBin.color(hbaInd,:);
end, xlim([-40,40]); title('Uncorrected: head-maze-distance > 20cm');

subplot(224); hold('on');
for hbaInd = 1:bins.hba.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) & randn(size(decoded.hba))>.5 & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1])) & WithinRanges(decoded.dst,[200,500]);      
[F] = cdfplot((decoded.lat(ind)-decoded.correction.lat(ind))/10);  F.Color = hbaBin.color(hbaInd,:);
end, xlim([-20,20]); title('corrected: head-maze-distance > 20cm');





figure(); hold('on');
for hbaInd = 1:bins.hba.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) & randn(size(decoded.hba))>.5 & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1]));      
[F,stats] = cdfplot((decoded.lat(ind)-decoded.correction.lat(ind))/10-1.8); F.Color = hbaBin.color(hbaInd,:);
end, xlim([-30,30]); 



[hfig,fig,fax,sax] = set_figure_layout(figure(666003),'A4','landscape',[],1.8,1.8,0.2,0.4);

hfig = open('/storage/share/Projects/EgoProCode2D/fig/untitled.fig');


sax = 

    [yind, yOffSet, xind, xOffSet] = deal(1, 0, 1, 0);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                                  fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    title(sax(end),{'Place Field'});                                                                                                                                                                                                                                                                                                                                                           
                                        

figure(); hold('on');
for hbaInd = 1:bins.hba.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) & randn(size(decoded.hba))>0.5 & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1]));      
[F,xi] = ksdensity((decoded.lat(ind)-decoded.correction.lat(ind))/10-0.8);
decLat.hba(hbaInd).res.all.mean = mean((decoded.lat(ind)-decoded.correction.lat(ind))/10-0.8);
decLat.hba(hbaInd).res.all.std  = std( (decoded.lat(ind))/10);
[decLat.hba(hbaInd).res.all.ttest.h,...
 decLat.hba(hbaInd).res.all.ttest.pval,...
 decLat.hba(hbaInd).res.all.ttest.ci,...
 decLat.hba(hbaInd).res.all.ttest.stats] = ttest((decoded.lat(ind)-decoded.correction.lat(ind))/10-0.8);
plot(xi,F,'-','Color',hbaBin.color(hbaInd,:));
[~,xid] = NearestNeighbour(xi,decLat.hba(hbaInd).res.all.mean);
line(xi(xid)*[1,1],[0,F(xid)],'color',bins.hba.color(hbaInd,:));
end, xlim([-30,30]); 
title('Decoded Lateral Position');

hbaInd = 1; indL = WithinRanges( decoded.phz, [4.5,5.5]) & randn(size(decoded.hba))>1 & WithinRanges( decoded.hba, bins.hba.edges(hbaInd+[0,1]));
hbaInd = 2; indC = WithinRanges( decoded.phz, [4.5,5.5]) & randn(size(decoded.hba))>1 & WithinRanges( decoded.hba, bins.hba.edges(hbaInd+[0,1]));
hbaInd = 3; indR = WithinRanges( decoded.phz, [4.5,5.5]) & randn(size(decoded.hba))>1 & WithinRanges( decoded.hba, bins.hba.edges(hbaInd+[0,1]));

[h,p,r,s] =ttest2((decoded.lat(indC)-decoded.correction.lat(indC))/10,(decoded.lat(indR)-decoded.correction.lat(indR))/10);

[h,p,r,s] =ttest2((decoded.lat(indC)-decoded.correction.lat(indC))/10,(decoded.lat(indL)-decoded.correction.lat(indL))/10);

decoded.clat = (decoded.lat-decoded.correction.lat)-8;
decoded.clt = (decoded.lat-decoded.correction.lat)-8;




[A,B] = pol2cart(circ_dist(THETA,diff(sectors([1,2]))),RR);
%[A,B] = pol2cart(THETA,RR);

% SUPFIG - Non-Uniformity of decoded anterioposterior position:
% subject-arena orientation and distance to maze center
opt.view = [0,90];
opt.clim = [-10,10];
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
ForAllSubplots('ylim([-45,45]);');
ForAllSubplots('xlim([-45,45]);');
fax = axes('Position',[0,0,1,1],'Visible','off');
text(fax,0.1,0.9,{'Non-Uniformity of decoded anterioposterior position ',...
                  'depedent on subject-arena orientation and distance to maze center'});
set(gcf(),'PaperOrientation','landscape');





% Fig 2 decoding subplots


figure(); hold('on');
for hbaInd = 1:bins.hba.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) & randn(size(decoded.hba))>.5 & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1]));      
[F,stats] = cdfplot((decoded.lat(ind)-decoded.correction.lat(ind))/10-1.8); F.Color = hbaBin.color(hbaInd,:);
end, xlim([-30,30]); 



[hfig,fig,fax,sax] = set_figure_layout(figure(666003),'A4','landscape',[],1.8,1.8,0.2,0.4);

hfig = open('/storage/share/Projects/EgoProCode2D/fig/untitled.fig');


[yind, yOffSet, xind, xOffSet] = deal(1, 0, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                    fig.page.ypos(yind)+yOffSet+globalYOffset,...
                    fig.subplot.width,                        ...
                    fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold(sax(end),'on');
title(sax(end),{'Place Field'});
                                        


figure(); hold('on');
for hbaInd = 1:bins.hba.count
ind = WithinRanges(decoded.phz,[4.5,5.5]) & randn(size(decoded.hba))>0.5 & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1]));      
[F,xi] = ksdensity((decoded.lat(ind)-decoded.correction.lat(ind))/10-0.8);
decLat.hba(hbaInd).res.all.mean = mean((decoded.lat(ind)-decoded.correction.lat(ind))/10-0.8);
decLat.hba(hbaInd).res.all.std  = std( (decoded.lat(ind))/10);
[decLat.hba(hbaInd).res.all.ttest.h,...
 decLat.hba(hbaInd).res.all.ttest.pval,...
 decLat.hba(hbaInd).res.all.ttest.ci,...
 decLat.hba(hbaInd).res.all.ttest.stats] = ttest((decoded.lat(ind)-decoded.correction.lat(ind))/10-0.8);
plot(xi,F,'-','Color',hbaBin.color(hbaInd,:));
[~,xid] = NearestNeighbour(xi,decLat.hba(hbaInd).res.all.mean);
line(xi(xid)*[1,1],[0,F(xid)],'color',bins.hba.color(hbaInd,:));
end, xlim([-30,30]); 
title('Decoded Lateral Position');



mask = ones([mBinHba.count,mBinHvl.count]);
mask(cout<100) = nan;
figure();
imagescnan({mBinHba.centers,mBinHvl.centers,(mout.*mask)'},[-6,6],'linear',true,'colorMap',@jet);
axis('xy');
