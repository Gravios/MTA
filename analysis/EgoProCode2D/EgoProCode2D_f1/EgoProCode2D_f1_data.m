
% interphase shift -> pairwise 
% add to sup fig center of mass for each placefield

% add lines for each phase pref group fig1 last subplots

% EgoProCode2D_f1_data.m

configure_default_args();
EgoProCode2D_load_data();

exampleUnit.trialIndex = 18;
exampleUnit.close.Xlims = [-200,400];
exampleUnit.close.Ylims = [-400,200];
exampleUnit.id = 11;
exampleUnit.maxRate = 25;
exampleUnit.index = find(units{exampleUnit.trialIndex}==exampleUnit.id);
exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];


exampleUnit.allo.xlims = [-300,300];
exampleUnit.allo.ylims = [-300,300];
exampleUnit.ego.xlims = [-300,300];
exampleUnit.ego.ylims = [-250,350];

exampleLFP.trialIndex = 18;
exampleLFP.unitId = 11;


uidsCA1 = unitsEgoCA1;
uidsCA3 = unitsEgoCA3;




% $$$ tid=19; unts = [10,27,50,53,66,172]
% $$$ tid=20; unts = [25];
% $$$ 
% $$$ 
% $$$ 
% $$$ usub = 0+[1:10];
% $$$ tid = 20;
% $$$ figure,
% $$$ disp(numel(units{tid}))
% $$$ for p = 1:3,
% $$$     for u = usub,
% $$$         if u > numel(units{tid})
% $$$             break;
% $$$         end
% $$$         subplot2(numel(usub),3,u-usub(1)+1,p);
% $$$         plot(ratemapsAlloThp{tid}{p},units{tid}(u),'text',1,[0,20],'colorMap',@jet);
% $$$         title(units{tid}(u))
% $$$     end
% $$$ end
% $$$ 
% $$$ usub = 5;
% $$$ tid = 20;
% $$$ figure,
% $$$ disp(numel(units{tid}))
% $$$ for p = 1:3,
% $$$     for u = usub,
% $$$         if u > numel(units{tid})
% $$$             break;
% $$$         end
% $$$         subplot2(numel(usub),3,u-usub(1)+1,p);
% $$$         plot(ratemapsAlloThp{tid}{p},units{tid}(u),'text',1,[0,20],'colorMap',@jet);
% $$$         title(units{tid}(u))
% $$$     end
% $$$ end
% $$$ 
% $$$ u = 2
% $$$ figure,
% $$$ for p = 1:3
% $$$     subplot(1,3,p)
% $$$     rmap = plot(ratemapsAlloThp{tid}{p},units{tid}(u));
% $$$     imagesc(log10(rmap'))
% $$$     caxis([-1.3,1.3])
% $$$     colormap(jet);
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ pfsBinGrid = cell([1,2]);
% $$$ [pfsBinGrid{:}] = ndgrid(ratemapsAlloThp{t}{p}.adata.bins{:});
% $$$ 
% $$$ %u = units{tid}==25;
%u = units{tid}==31;
%t = 1;tid = 18; u = units{tid}==61;
% $$$ t = 18;tid = 18; u = units{tid}==11;
% $$$ t = 20;tid = 20; u = units{tid}==25;
% $$$ t =  3;tid =  3; u = units{tid}==173;
% $$$ t = 30;tid = 30; u = units{tid}==29;
% $$$ 
% $$$ figure,
% $$$     [~,mpx] = pft{t}.maxRate(units{tid}(u));
% $$$ for p = 1:3
% $$$     subplot(3,3,p)
% $$$     rmape = plot(pfet{t}{p},units{tid}(u));
% $$$     imagescnan({pfet{1}{p}.adata.bins{:},log2(rmape)},...
% $$$                'colorLimits', [0,4],...
% $$$                'nanRGB',[0.75,0.75,0.75],...
% $$$                'colorMap',@jet);
% $$$     hold('on');
% $$$     Lines([],0,'w');
% $$$     Lines(0,[],'w');
% $$$     circle(0,0,250,'r-');    
% $$$     axis('xy')    
% $$$     caxis([0,4])
% $$$     xlim([-500,500]);
% $$$     ylim([-500,500]);
% $$$     colormap(jet);
% $$$ 
% $$$     subplot(3,3,p+3)
% $$$     rmapa = plot(ratemapsAlloThp{t}{p},units{tid}(u));
% $$$     imagescnan({ratemapsAlloThp{t}{p}.adata.bins{:},log2(rmapa)'},...
% $$$                'colorLimits', [0,4],...
% $$$                'nanRGB',[0.75,0.75,0.75],...
% $$$                'colorMap',@jet);
% $$$     hold('on');
% $$$     circle(mpx(1),mpx(2),250,'r-');    
% $$$     axis('xy')
% $$$     caxis([0,4])
% $$$     colormap(jet);
% $$$     
% $$$     subplot(3,3,p+6);
% $$$     gtinda = ~isnan(rmapa);
% $$$     gtinde = ~isnan(rmape);
% $$$     rmapa(isnan(rmapa)|rmapa==0) = eps;
% $$$     rmape(isnan(rmape)|rmape==0) = eps;
% $$$     plot(sort(rmapa(:)),sort(rmape(:)),'.');
% $$$     title(num2str(sum(rmape(:)>2)./(sum(rmapa(:)>2)+sum(rmape(:)>2))));
% $$$     line([0,14],[0,14]);
% $$$     
% $$$     si = [0,0];
% $$$     sp = [0,0];
% $$$     
% $$$     MRate = mean(rmapa(gtinda),'omitnan');
% $$$     si(1) = nansum(1/sum(double(gtinda(:))).*(rmapa(gtinda)./MRate) .* log2(rmapa(gtinda)./MRate));
% $$$     sp(1) = nansum(1/sum(double(gtinda(:))).*rmapa(gtinda)).^2./nansum(1/sum(double(gtinda(:))).*rmapa(gtinda).^2);
% $$$     MRate = mean(rmape(gtinde),'omitnan');
% $$$     si(2) = nansum(1/sum(double(gtinde(:))).*(rmape(gtinde)./MRate) .* log2(rmape(gtinde)./MRate));
% $$$     sp(2) = nansum(1/sum(double(gtinde(:))).*rmape(gtinde)).^2./nansum(1/sum(double(gtinde(:))).*rmape(gtinde).^2);
% $$$     round([si,sp],2)'
% $$$ end
% $$$ 
% $$$ 
% $$$ figure,
% $$$ histogram(log10(rmapa(:)+0.001),linspace(-3,1.3,50));
% $$$ hold('on');
% $$$ histogram(log10(rmape(:)+0.001),linspace(-3,1.3,50));

% $$$ tid = 19; u = units{tid}==146;
% $$$ figure
% $$$ for p = 1:3
% $$$     for a = 1:3
% $$$         subplot2(3,3,4-p,a);     
% $$$ plot(pfs{tid}{p,a},...
% $$$      units{tid}(u),...
% $$$      1, ...
% $$$      '',...
% $$$      [0,16],...
% $$$      false,...
% $$$      'colorMap',@jet,...
% $$$      'flipAxesFlag',true);
% $$$         hold('on');
% $$$         Lines(0,[],'w');/
% $$$         Lines([],0,'w');
% $$$     
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ tid = 21; u = units{tid}==6;
% $$$ figure
% $$$ for p = 1:3
% $$$     for a = 1:4
% $$$         subplot2(3,4,4-p,a);
% $$$         try
% $$$         plot(pfv{tid}{p,a},...
% $$$              units{tid}(u),...
% $$$              1, ...
% $$$              '',...
% $$$              [0,8],...
% $$$              false,...
% $$$              'colorMap',@jet,...
% $$$              'nanColor',[.7,.7,.7],...
% $$$              'flipAxesFlag',true);
% $$$         hold('on');
% $$$         Lines(0,[],'w');
% $$$         Lines([],0,'w');
% $$$         end
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ 
k = 1;
for tid = 1:numel(Trials)
for u = 1:numel(unitsEgo{tid})
    [mxr,mpx] = placeFieldsNoRear{tid}.maxRate(unitsEgo{tid}(u));
    for p = 1:3
    rmape = plot(pfet{tid}{p},unitsEgo{tid}(u));
    rmapa = plot(ratemapsAlloThp{tid}{p},unitsEgo{tid}(u));
    gtinda = ~isnan(rmapa);
    gtinde = ~isnan(rmape);
    rmapa(isnan(rmapa)|rmapa==0) = eps;
    rmape(isnan(rmape)|rmape==0) = eps;
    MRate = mean(rmapa(gtinda),'omitnan');
    si(1) = nansum(1/sum(double(gtinda(:))).*(rmapa(gtinda)./MRate) .* log2(rmapa(gtinda)./MRate));
    sp(1) = nansum(1/sum(double(gtinda(:))).*rmapa(gtinda)).^2./nansum(1/sum(double(gtinda(:))).*rmapa(gtinda).^2);
    MRate = mean(rmape(gtinde),'omitnan');
    si(2) = nansum(1/sum(double(gtinde(:))).*(rmape(gtinde)./MRate) .* log2(rmape(gtinde)./MRate));
    sp(2) = nansum(1/sum(double(gtinde(:))).*rmape(gtinde)).^2./nansum(1/sum(double(gtinde(:))).*rmape(gtinde).^2);
    sia(k,p,:) = si;
    spa(k,p,:) = sp;
    end
    mxa(k,:) = mxr;
    k = k+1;
end
end


figure,
trlI = 20
unit = 21;
[mxr,mxp] = ratemapsAlloThp{trlI}{2}.maxRate(unit);
for phzI = 1:pc
    subplot2(3,2,pc+1-phzI,1);
    plot(ratemapsAlloThp{trlI}{phzI}, unit, [], 'colorbar', [0,mxr*1.5], true, [], false,[],@jet);
        ylim([-200,200]+mxp(2));
        xlim([-200,200]+mxp(1));
        daspect([1,1,1]);
    subplot2(3,2,pc+1-phzI,2);
        plot(pfet{trlI}{phzI}, unit, [], 'colorbar', [0,mxr*1.5], false, [], true,[],@jet);
        ylim([-150,250]);
        xlim([-200,200]);    
        daspect([1,1,1]);
end


trlI = 20;
phz = load_theta_phase(Trials{trlI}, sampleRate);
headYawCorrection = Trials{trlI}.meta.correction.headYaw;
headCenterCorrection = Trials{trlI}.meta.correction.headCenter;
% COMPUTE head basis
hvec = xyz{trlI}(:,'nose',[1,2])-xyz{trlI}(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(headYawCorrection),-sin(headYawCorrection);...
                  sin(headYawCorrection), cos(headYawCorrection)],...
                 [2,3],...
                 [1,2]);
hvfl = fet_href_HXY(Trials{trlI}, sampleRate, [], 'trb', 2.4);
hafl = circshift(hvfl.data,1)-hvfl.data;
%bvfl = fet_bref_BXY(Trials{trlI}, sampleRate, [], 'trb', 5);
%bafl = circshift(bvfl(:,1),1)-bvfl(:,1);

sts = [Trials{trlI}.stc{'walk&theta'}];
sts = [Trials{trlI}.stc{'walk+pause&theta'}];
%sts = [Trials{trlI}.stc{'pause&theta'}];
sts.resample(xyz{trlI});

unit = 103;
ures = spk{trlI}(unit,sts);
%ures = spk{trlI}(unit);
ures = ures(within_ranges(phz(ures),[4,5.5]));
[mxr,mxp] = ratemapsAlloThp{trlI}{2}.maxRate(unit);
%[mxr,mpx] = pft{trlI}.maxRate(unit);
exy = [bsxfun(@plus,                                            ...
              multiprod(bsxfun(@minus,                          ...
                               mxp,                             ...
                               sq(xyz{trlI}(:,'hcom',[1,2]))),  ...
                        hvec,2,[2,3]),                          ...
              headCenterCorrection)];
figure,
ind =  ures(abs(exy(ures,2))<200 ...
       & within_ranges(exy(ures,1),[-300,300]));
subplot(2,3,1);
scatter(hvfl(ind,1) , exy(ind,1), 10, hafl(ind,1)*125,'Filled');
grid('on'); hold('on'); colormap('jet'); colorbar(); caxis([-100,100]); xlim([-20,100]); ylim([-300,300]);
plot(0,median(exy(ind,1)),'*m');
subplot(2,3,2);
scatter(hafl(ind,1)*125, exy(ind,1), 10, hvfl(ind,1),'Filled');
grid('on'); hold('on'); colormap('jet'); colorbar(); caxis([-30,100]); xlim([-100,100]); ylim([-300,300]);
plot(0,median(exy(ind,1)),'*m');
subplot(2,3,3);
histogram(hafl(ind,1), linspace(-2,2,24));
R = corrcoef(hafl(ind,1), exy(ind,1))
R = corrcoef(hvfl(ind,1), exy(ind,1))
median(hafl(ind,1))
ures = spk{trlI}(unit,sts);
ures = ures(within_ranges(phz(ures),[0.5,2]));
ind =  ures(abs(exy(ures,2))<200 ...
       & within_ranges(exy(ures,1),[-300,300]));
subplot(2,3,4);
scatter(hvfl(ind,1) , exy(ind,1), 10, hafl(ind,1)*125,'Filled');
grid('on'); hold('on'); colormap('jet'); colorbar(); caxis([-100,100]); xlim([-20,100]); ylim([-300,300]);
plot(0,median(exy(ind,1)),'*m');
subplot(2,3,5);
scatter(hafl(ind,1)*125, exy(ind,1), 10, hvfl(ind,1),'Filled');
grid('on'); hold('on'); colormap('jet'); colorbar(); caxis([-30,100]); xlim([-100,100]); ylim([-300,300]);
plot(0,median(exy(ind,1)),'*m');
subplot(2,3,6);
histogram(hafl(ind,1), linspace(-2,2,24));
median(hafl(ind,1))
% $$$ 
% $$$ ind =  abs(exy(:,2))<200 ...
% $$$        & within_ranges(exy(:,1),[-100,50]) ...
% $$$        & within_ranges(1:numel(exy(:,1)),sts);
% $$$ figure,
% $$$ histogram(hafl(ind,1),linspace(-2,2,30))
% $$$ median(hafl(ind,1))
% $$$ mean(hafl(ind,1))
% $$$ figure,
% $$$ plot(hvfl(ind,1),hafl(ind,1),'.')
% $$$ hold('on')
% $$$ ures = spk{trlI}(unit,sts);
% $$$ ures = ures(within_ranges(phz(ures),[0.5,2.5]));
% $$$ %ures = ures(within_ranges(phz(ures),[4,6]));
% $$$ ind =  ures(abs(exy(ures,2))<200 ...
% $$$        & within_ranges(exy(ures,1),[-300,300]));
% $$$ plot(hvfl(ind,1),hafl(ind,1),'.r','MarkerSize',10)
% $$$ Lines([],0,'k');
% $$$ figure,hist(hafl(ind,1),linspace(-2,2,20))
% $$$ median(hafl(ind,1))
% $$$ mean(hafl(ind,1))


normalization = 'count';
ures = spk{trlI}(unit,sts);
ures = ures(within_ranges(phz(ures),[0.5,1.5]));
dres =  ures(abs(exy(ures,2))<150 ...
       & within_ranges(hvfl(ures,1),[-10,50]) ...
       & within_ranges(exy(ures,1),[-200,200]));
ures = spk{trlI}(unit,sts);
ures = ures(within_ranges(phz(ures),[3,6]));
ares =  ures(abs(exy(ures,1))<250 & within_ranges(hvfl(ures,1),[-10,50]));
abins = linspace(-250,250,14);
%abins = linspace(-20,100,6);
abinc = mean([abins(1:end-1);abins(2:end)]);
accg = zeros(201,5);
for a = 1:numel(abinc)
    adres = dres;
    aares = ares;
    %aares = aares(within_ranges(exy(aares,1),abins([a,a+1])));
    adres = adres(within_ranges(exy(adres,1),abins([a,a+1])));
    [mccg,t,pairs] = CCG([adres;aares],[ones(size(adres));2*ones(size(aares))],2,100,250,[1,2],normalization);
    accg(:,a) = mccg(:,1,2);
end

%STOP THAT
figure
imagesc(t,abinc,accg')
axis('xy')
colormap('jet')

normalization = 'count';
ures = spk{trlI}(unit,sts);
ures = ures(within_ranges(phz(ures),[0.5,2.5]));
dres =  ures(abs(exy(ures,2))<200 & within_ranges(hvfl(ures,1),[-10,80]) ...
       & within_ranges(exy(ures,1),[-200,200]));
ures = spk{trlI}(unit,sts);
ures = ures(within_ranges(phz(ures),[3,6]));
ares =  ures(abs(exy(ures,2))<200 & within_ranges(hvfl(ures,1),[-10,80]));
figure
% $$$ hist2([ (NearestNeighbour(ares,dres,'both')-dres)./250,exy(dres,1)],linspace(-0.5,0.5,36),linspace(-200,200,14));
[outl,xb,yb] = hist2([ (NearestNeighbour(ares,dres,'left')-dres)./250,exy(dres,1)],linspace(-0.8,0.8,44),linspace(-200,200,14));
[outr     ] = hist2([ (NearestNeighbour(ares,dres,'right')-dres)./250,exy(dres,1)],linspace(-0.8,0.8,44),linspace(-200,200,14));
imagesc(xb,yb,[outl+outr]');
%imagesc(xb,yb,bsxfun(@rdivide,[outl+outr],sum([outl+outr],1))');
axis('xy');
colormap('cool');
Lines(0,[],'k');

% dec before ace +
% dec after ace -


figure
hist2([ (NearestNeighbour(ares(2:2:end),ares(1:2:end-1))-ares(1:2:end-1))./250,exy(ares(1:2:end-1),1)],linspace(-0.2,0.2,18),linspace(-200,200,20))
colormap('cool');
Lines(0,[],'k');

nu = abs((NearestNeighbour(ares(2:2:end),ares(1:2:end-1))-ares(1:2:end-1))./250)<0.1;


figure,
for hvfI = 1:bins.hvf.count
    subplot(1, bins.hvf.count, hvfI);
    ind = ures(within_ranges(hvfl(ures,1),bins.hvf.edges(hvfI:hvfI+1)) ...
               & abs(exy(ures,2))<200 ...
               & within_ranges(exy(ures,1),[-200,300]));
scatter( hvfl(ind,2), exy(ind,1), 10, hvfl(ind,1),'Filled');
hold('on');
plot(0,median(exy(ind,1)),'*m');
colormap('jet')
colorbar();
caxis([-20,60]);
xlim([-200,200]);
ylim([-200,300]);
end

figure,
ind =  ures(abs(exy(ures,2))<200 ...
       & within_ranges(exy(ures,1),[-200,300]));
scatter( hvfl(ind,1), exy(ind,1), 10, hvfl(ind,1),'Filled');
hold('on');
plot(0,median(exy(ind,1)),'*m');
colormap('jet')
colorbar();
caxis([-20,60]);
xlim([-20,100]);
ylim([-200,300]);


%caxis([0,80]);


18,42



% $$$ 
% $$$ figure,
% $$$ ind = mxa<10;
% $$$ for p = 1:3
% $$$     subplot2(3,1,4-p,1);
% $$$     histogram(sq(sia(ind,p,2)-sia(ind,p,1)),linspace(-0.6,0.6,20));
% $$$ end
% $$$ 
% $$$ 
% $$$ figure,
% $$$ ind = mxa<10;
% $$$ for p = 1:3
% $$$     subplot2(3,2,4-p,1);
% $$$     histogram(sq(sia(ind,p,2)-sia(ind,p,1)),linspace(-0.7,0.7,20));
% $$$ end
% $$$ ind = ~ind;
% $$$ for p = 1:3
% $$$     subplot2(3,2,4-p,2);
% $$$     histogram(sq(sia(ind,p,2)-sia(ind,p,1)),linspace(-0.7,0.7,20));
% $$$ end
% $$$ 
% $$$ 
% $$$ figure,
% $$$ ind = uidsCA1
% $$$ for p = 1:3
% $$$     subplot2(3,2,4-p,1);
% $$$     histogram(sq(sia(ind,p,2)-sia(ind,p,1)),linspace(-0.7,0.7,20));
% $$$ end
% $$$ ind = uidsCA3
% $$$ for p = 1:3
% $$$     subplot2(3,2,4-p,2);
% $$$     histogram(sq(sia(ind,p,2)-sia(ind,p,1)),linspace(-0.7,0.7,20));
% $$$ end
% $$$ 
% $$$ 
% $$$ figure,
% $$$ ind = uidsCA1
% $$$ for p = 1:3
% $$$     subplot2(3,2,4-p,1);
% $$$     histogram(sq(sia(ind,p,2)),linspace(-0,5,20));
% $$$ end
% $$$ ind = uidsCA3
% $$$ for p = 1:3
% $$$     subplot2(3,2,4-p,2);
% $$$     histogram(sq(sia(ind,p,2)),linspace(0,5,20));
% $$$ end
% $$$ 
% $$$ 
% $$$ figure,
% $$$ ind = uidsCA1;
% $$$ for p = 1:3
% $$$     subplot2(3,2,4-p,1);
% $$$     plot(log2(mxa(ind)),sia(ind,p,2)-sia(ind,p,1),'.','MarkerSize',20);
% $$$     xlabel('Allo SI (bits)');
% $$$     ylabel('Ego SI (bits)');
% $$$ end
% $$$ ind = uidsCA3;
% $$$ for p = 1:3
% $$$     subplot2(3,2,4-p,2);
% $$$     plot(log2(mxa(ind)),sia(ind,p,2)-sia(ind,p,1),'.','MarkerSize',20);
% $$$     xlabel('Allo SI (bits)');
% $$$     ylabel('Ego SI (bits)');
% $$$ end
% $$$ 
% $$$ plot(mxa,sia(:,3,2)-sia(:,3,1),'.','MarkerSize',20);
% $$$ 
% $$$ 
% $$$ phaseLabels = {'descend','trough','ascend'}
% $$$ figure,
% $$$ ind = uidsCA1;
% $$$ for p = 1:3
% $$$     subplot2(3,2,4-p,1);
% $$$     scatter(sia(ind,p,1),sia(ind,p,2),40,log2(mxa(ind)),'filled');
% $$$     line([0,4],[0,4],'Color','k');
% $$$     xlabel('Allo SI (bits)');
% $$$     ylabel('Ego SI (bits)');
% $$$     title(phaseLabels{p})
% $$$     if p ==3 
% $$$             title(['CA1 ',phaseLabels{p}])
% $$$     end
% $$$     
% $$$ end
% $$$ ind = uidsCA3;
% $$$ for p = 1:3
% $$$     subplot2(3,2,4-p,2);
% $$$     scatter(sia(ind,p,1),sia(ind,p,2),40,log2(mxa(ind)),'filled');    
% $$$     line([0,4],[0,4],'Color','k');    
% $$$     xlabel('Allo SI (bits)');
% $$$     ylabel('Ego SI (bits)');    
% $$$     title(phaseLabels{p})
% $$$     if p ==3 
% $$$             title(['CA3 ',phaseLabels{p}])
% $$$     end
% $$$     
% $$$ end
% $$$ colormap('jet');
% $$$ 
% $$$ 

% ANALYSIS 
% Place field size comparison between ego and allo fields
% PROBLEM 
% estimation of field size is difficult in the allo due to deformation away from the gaussian placefield
% RECOMENDATION 
% note primary field position from theta state allo map and then compute the spatial information of the
% cell for each phase bin within a specified radius of the field position in theta-allo-pfs
% $$$ 
% $$$ % $$$ 
% $$$ % $$$ 
% $$$ figure,
% $$$ tid = 20;
% $$$ for p = 1:3,
% $$$     for u = 6
% $$$         subplot2(1,3,1,p);
% $$$         plot(ratemapsAlloThp{tid}{p},units{tid}(46),'text',1,[0,8],'colorMap',@jet);
% $$$     end
% $$$ end


%%% LOAD behavioral data at spike times
spkv = MjgER2016_load_spikeVars(Trials,units,sessionList,[],[],[],[],overwrite); 


%%% LOAD patch modelm 
rat = load_patch_model('rat');


%%% 3 Part theta phase color map
pclr = cool(3);

%%% state periods
sper = Trials{exampleUnit.trialIndex}.stc{'x+p&t&a'}; % State Periods

%%% head direction rotation matrix
hvec = transform_vector_to_rotation_matrix(xyz{exampleUnit.trialIndex},{'hcom','nose'}, Trial.meta.correction.headYaw);


[mxr,mxp] = pft{exampleUnit.trialIndex}.maxRate(units{exampleUnit.trialIndex}(exampleUnit.index));        
pfstrj = MTADfet.encapsulate( ...
    Trials{exampleUnit.trialIndex},...
    multiprod(bsxfun(@minus,...
                     mxp,...
                     sq(xyz{exampleUnit.trialIndex}(:,'hcom',[1,2]))),...
              hvec,2,[2,3]),...
    sampleRate,...
    'pfstrj','ptrj','t');


%%% ???
spkPhzBin = [0, 2*pi/3, 4*pi/3, 2*pi];
spktemp = Trials{exampleUnit.trialIndex}.load('spk',250,sper,exampleUnit.id,'');
phztemp = load_theta_phase(Trials{exampleUnit.trialIndex},250);
spkPhzInd = discretize(phztemp(spktemp.res), spkPhzBin);        


% rat model
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false);
subject = update_subject_patch(subject,'body',[],false);


%%%<<< COMPUTATION - theta phase resolved egocentric field center of mass --------------
ucounter = 1;
egoMeanRmapPos = [];
egoMaxRmapPos = [];
egoSize = [];
egoMeanRmapRate = [];
egoMaxRmapRate = [];

t = 20;p=1;
binSubsetX = abs(pfet{t}{p}.adata.bins{1})<300;
binSubsetY = abs(pfet{t}{p}.adata.bins{2})<300;
mapPosition = cell([1,2]);
[mapPosition{:}] = ndgrid(pfet{t}{p}.adata.bins{1}(binSubsetX),...
                          pfet{t}{p}.adata.bins{2}(binSubsetY));
mapPosition = cat(numel(mapPosition)+1,mapPosition{:});


for t = 1:numel(Trials)
    for u = 1:numel(unitsEgo{t})
        unit = unitsEgo{t}(u);
        for p = 1:3;
            rmap = plot(pfet{t}{p},unit,[],[],[],false);
            rmap = rmap(binSubsetX,binSubsetY);
            nanmap = double(~isnan(rmap));
            nanmap(nanmap==0) = nan;
            rmap = rmap.*fliplr(nanmap);
            rmap(isnan(rmap)) = 0;
            rmap(rmap<2) = 0;
            nrmap =rmap./sum(rmap(:),'omitnan');
            rmapCenter = sq(sum(sum(bsxfun(@times,nrmap,mapPosition),'omitnan'),'omitnan'))';
            egoMeanRmapPos(ucounter,p,:) = rmapCenter./10 + [-2,1.6]; % Convert from mm to cm and correction
            egoSize(ucounter,p) = sum(nniz(nrmap(:)));
            egoMeanRmapRate(ucounter,p,:) = mean(nonzeros(rmap));
            [~,maxPos] = max(nrmap(:));
            egoMaxRmapRate(ucounter,p) = rmap(maxPos);
            if ~isempty(maxPos)
                [maxX,maxY] = ind2sub(size(nrmap),maxPos);
                egoMaxRmapPos(ucounter,p,:) = mapPosition(maxX,maxY,:);
            else
                egoMaxRmapPos(ucounter,p,:) = nan([1,1,1,2]);
            end
        end
        ucounter = ucounter+1;
    end
end


figure,
uinds = uidsCA1(ismember(egoCluSessionMap(uidsCA1,1),[3,4,5]));
uinds = uidsCA1(ismember(egoCluSessionMap(uidsCA1,1),[18:25]));
subplot(311);
histogram(egoMeanRmapPos(uinds,2,1)-egoMeanRmapPos(uinds,1,1),linspace(-5,20,20))
title('trg-dsc');
subplot(312);
histogram(egoMeanRmapPos(uinds,3,1)-egoMeanRmapPos(uinds,2,1),linspace(-5,20,20))
title('asc-trg');
subplot(313);
histogram(egoMeanRmapPos(uinds,3,1)-egoMeanRmapPos(uinds,1,1),linspace(-5,20,20))
title('asc-dsc');

% STATS
egoMeanPosPvalForward = nan([6,6]);
egoMeanPosMeanForward = nan([6,6]);
tMeanPos = repmat(egoMeanRmapPos(:,:,1),[1,2]);
for x = 1:6
    for p = 1:6
        if x==p, continue, end
        if p <= 3
            uind1 = uidsCA1;
        else
            uind1 = uidsCA3;
        end
        if x <= 3
            uind2 = uidsCA1;
        else
            uind2 = uidsCA3;
        end
        egoMeanPosPvalForward(x,p) = log10(ranksum(tMeanPos(uind2,x), tMeanPos(uind1,p)));
        egoMeanPosMeanForward(x,p) = median(tMeanPos(uind2,p)) - median(tMeanPos(uind1,x));
% $$$         [~,tpval] = ttest(tMeanPos(uind2,x)-tMeanPos(uind1,p));
% $$$         egoMeanPosPvalForward(x,p) = log10(tpval);        
% $$$         egoMeanPosMeanForward(x,p) = median(tMeanPos(uind2,p) - tMeanPos(uind1,x));        
    end
end
egoMeanPosPvalForward(egoMeanPosPvalForward > -1.301) = nan;
egoMeanPosPvalForward




egoMeanPosPairwiseLabels = {'trgh - dsc',...
                            'asc - trgh',...
                            'asc - dsc'};

%                                                  units phz dim                     units phz dim

egoMeanPosPairwiseDiffFwdCA1 = zeros([numel(uidsCA1),3]);
egoMeanPosPairwiseDiffFwdCA3 = zeros([numel(uidsCA3),3]);
egoMeanPosPairwiseDiffLatCA1 = zeros([numel(uidsCA1),3]);
egoMeanPosPairwiseDiffLatCA3 = zeros([numel(uidsCA3),3]);
egoMeanPosPairwiseDiffFwdCA1(:,1) = egoMeanRmapPos(uidsCA1,  2,  1) - egoMeanRmapPos(uidsCA1,  1,  1);
egoMeanPosPairwiseDiffFwdCA1(:,2) = egoMeanRmapPos(uidsCA1,  3,  1) - egoMeanRmapPos(uidsCA1,  2,  1);
egoMeanPosPairwiseDiffFwdCA1(:,3) = egoMeanRmapPos(uidsCA1,  3,  1) - egoMeanRmapPos(uidsCA1,  1,  1);
egoMeanPosPairwiseDiffFwdCA3(:,1) = egoMeanRmapPos(uidsCA3,  2,  1) - egoMeanRmapPos(uidsCA3,  1,  1);
egoMeanPosPairwiseDiffFwdCA3(:,2) = egoMeanRmapPos(uidsCA3,  3,  1) - egoMeanRmapPos(uidsCA3,  2,  1);
egoMeanPosPairwiseDiffFwdCA3(:,3) = egoMeanRmapPos(uidsCA3,  3,  1) - egoMeanRmapPos(uidsCA3,  1,  1);
egoMeanPosPairwiseDiffLatCA1(:,1) = egoMeanRmapPos(uidsCA1,  2,  2) - egoMeanRmapPos(uidsCA1,  1,  2);
egoMeanPosPairwiseDiffLatCA1(:,2) = egoMeanRmapPos(uidsCA1,  3,  2) - egoMeanRmapPos(uidsCA1,  2,  2);
egoMeanPosPairwiseDiffLatCA1(:,3) = egoMeanRmapPos(uidsCA1,  3,  2) - egoMeanRmapPos(uidsCA1,  1,  2);
egoMeanPosPairwiseDiffLatCA3(:,1) = egoMeanRmapPos(uidsCA3,  2,  2) - egoMeanRmapPos(uidsCA3,  1,  2);
egoMeanPosPairwiseDiffLatCA3(:,2) = egoMeanRmapPos(uidsCA3,  3,  2) - egoMeanRmapPos(uidsCA3,  2,  2);
egoMeanPosPairwiseDiffLatCA3(:,3) = egoMeanRmapPos(uidsCA3,  3,  2) - egoMeanRmapPos(uidsCA3,  1,  2);




egoMeanPosPvalLateral = nan([6,6]);
egoMeanPosMeanLateral = nan([6,6]);
tMeanPos = repmat(egoMeanRmapPos(:,:,2),[1,2]);
for x = 1:6
    for p = 1:6
        if x==p, continue, end
        if p <= 3
            uind1 = uidsCA1;
        else
            uind1 = uidsCA3;
        end
        if x <= 3
            uind2 = uidsCA1;
        else
            uind2 = uidsCA3;
        end
        egoMeanPosPvalLateral(x,p) = log10(ranksum(tMeanPos(uind2,x), tMeanPos(uind1,p)));
        egoMeanPosMeanLateral(x,p) = median(tMeanPos(uind2,p)) - median(tMeanPos(uind1,x));
    end
end
egoMeanPosPvalLateral(egoMeanPosPvalLateral>-1.301) = nan;



%%%>>>
%---------------------------------------------------------------------------------------





% $$$ uids = [1:164];


% DETECTING tp-ego-ratemap centers

% ACCUMULATE place fields
t = 20;
[pfstats] = PlaceFieldStats(Trials{t},pfet{t}{2},20,2,true,false);
for t = 1:30,
    disp(['[INFO] accumulating placefield stats: ' Trials{t}.filebase]);
    for u = unitsEgo{t},
        [pfstats(end+1,1)] = PlaceFieldStats(Trials{t},pfet{t}{1},u,2,false,false);
        [pfstats(end  ,2)] = PlaceFieldStats(Trials{t},pfet{t}{2},u,2,false,false);
        [pfstats(end  ,3)] = PlaceFieldStats(Trials{t},pfet{t}{3},u,2,false,false);
    end
end
pfstats(1,:) = [];




apfstats = CatStruct(pfstats(1,:),[],2);
for field = fieldnames(pfstats)'
    field = field{1};
    if numel(apfstats.(field)) > 3
        apfstats.(field) = permute(apfstats.(field),[5,1,2,3,4]);
    end
    for uid = 2:size(egoCluSessionMap,1)
        apfstats.(field) = cat(1,apfstats.(field),cat(2,pfstats(uid,:).(field)));
    end
end

% COMPUTE center point of all good patches
apfstats.patchCenters = nan([size(egoCluSessionMap,1),3,2,2]);
for uid = 1:size(egoCluSessionMap,1),
    t = egoCluSessionMap(uid,1);
    for s = 1:3
        for p = 1:2
            if apfstats.patchPFR(uid,s,p)>1
                rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
                rinds = rinds(nniz(rinds),:);
                apfstats.patchCenters(uid,s,p,:) = [mean(pfet{t}{s}.adata.bins{1}(rinds(:,1))),...
                                    mean(pfet{t}{s}.adata.bins{2}(rinds(:,2)))];
            end
        end
    end
end


% $$$ % DIAGNOSTIC figure
% $$$ t = 20
% $$$ figure,
% $$$ for u = units{t};
% $$$     uid = find(ismember(egoCluSessionMap,[t,u],'rows'));
% $$$     mrate = max(apfstats.peakFR(uid,:));
% $$$     clf();
% $$$     for s = 1:3
% $$$         subplot(1,3,s);
% $$$         plot(pfet{t}{s},u,1,'colorbar',[0,mrate]);
% $$$         hold('on');
% $$$         for p = 1:2
% $$$             if apfstats.patchPFR(uid,s,p)>1
% $$$                 rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
% $$$                 rinds = rinds(nniz(rinds),:);
% $$$                 if size(rinds,1)>40
% $$$                     plot(pfet{t}{s}.adata.bins{1}(rinds(:,1)),...
% $$$                          pfet{t}{s}.adata.bins{2}(rinds(:,2)),...
% $$$                          'm*');
% $$$                     circle(apfstats.patchCenters(uid,s,p,1),apfstats.patchCenters(uid,s,p,2),150);
% $$$                 end        
% $$$             end
% $$$         end
% $$$         title(num2str(u))
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end


% ACCUMULATE independent patches 
%    IF two patchCenters are closer than 150 mm then count patch with the 
%       final patch -> mean center and radius = 150 + dist(center1,center2)


% COMPUTE center point of all good patches
apfstats.patchCenters = nan([size(egoCluSessionMap,1),3,2,2]);
for uid = 1:size(egoCluSessionMap,1),
    t = egoCluSessionMap(uid,1);
    for s = 1:3
        for p = 1:2
            if apfstats.patchPFR(uid,s,p)>1
                rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
                rinds = rinds(nniz(rinds),:);
                apfstats.patchCenters(uid,s,p,:) = [mean(pfet{t}{s}.adata.bins{1}(rinds(:,1))),...
                                    mean(pfet{t}{s}.adata.bins{2}(rinds(:,2)))];
            end
        end
    end
end

cdist = nan([size(egoCluSessionMap,1),3,3,2,2]);
for uid = 1:size(egoCluSessionMap,1),
    t = egoCluSessionMap(uid,1);
    for s1 = 1:3
        for s2 = 1:3
            for p1  = 1:2
                for p2 = 1:2
                    cdist(uid,s1,s2,p1,p2) = sqrt(sum(sq(apfstats.patchCenters(uid,s1,p1,:) - ...
                                                         apfstats.patchCenters(uid,s2,p2,:)).^2));
                end
            end
        end
    end
end
tempPatchCenter = nan([size(egoCluSessionMap,1),3,3,2,2,2]);
tempPatchRadius = nan([size(egoCluSessionMap,1),3,3,2,2]);
% MERGE close patches
for uid = 1:size(egoCluSessionMap,1),
    for s1 = 1:3
        for s2 = 1:3
            for p1 = 1:2
                for p2 = 1:2                
                    if cdist(uid,s1,s2,p1,p2) < 150           ...
                            && apfstats.patchPFR(uid,s1,p1)>1 ...
                            && apfstats.patchPFR(uid,s2,p2)>1 ...
                            && apfstats.patchPFR(uid,s1,p2)>1 ...
                            && apfstats.patchPFR(uid,s2,p1)>1 ...
                            && sum(nniz(apfstats.patchRateMap(uid,s1,p1)))>40 ...
                            && sum(nniz(apfstats.patchRateMap(uid,s1,p2)))>40 ...
                            && sum(nniz(apfstats.patchRateMap(uid,s2,p1)))>40 ...
                            && sum(nniz(apfstats.patchRateMap(uid,s2,p2)))>40 
                        % 18 comparisons
                        % don't look at same state patch centers 
                        tempPatchCenter(uid,s1,s2,p1,p2,:) = (apfstats.patchCenters(uid,s1,p1,:)+apfstats.patchCenters(uid,s2,p2,:))./2;
                        tempPatchRadius(uid,s1,s2,p1,p2) = 150 + cdist(uid,s1,s2,p1,p2);
                        cdist(uid,s1,s2,p1,p2) = 0;
                        apfstats.patchCenters(uid,s1,p1,:) = tempPatchCenter(uid,s1,s2,p1,p2,:);
                        apfstats.patchCenters(uid,s2,p2,:) = tempPatchCenter(uid,s1,s2,p1,p2,:);
                    end
                end
            end
        end
    end
end

% RECOMPUTE patch distance after merger
cdist = nan([size(egoCluSessionMap,1),3,3,2,2]);
for uid = 1:size(egoCluSessionMap,1),
    t = egoCluSessionMap(uid,1);
    for s1 = 1:3
        for s2 = 1:3
            for p1  = 1:2
                for p2 = 1:2
                    cdist(uid,s1,s2,p1,p2) = sqrt(sum(sq(apfstats.patchCenters(uid,s1,p1,:) - ...
                                                         apfstats.patchCenters(uid,s2,p2,:)).^2));
                end
            end
        end
    end
end
% MERGE close patches again
tempPatchCenter = nan([size(egoCluSessionMap,1),3,3,2,2,2]);
tempPatchRadius = nan([size(egoCluSessionMap,1),3,3,2,2]);
for uid = 1:size(egoCluSessionMap,1),
    for s1 = 1:3
        for s2 = 1:3
            for p1 = 1:2
                for p2 = 1:2                
                    if cdist(uid,s1,s2,p1,p2) < 150           
                        % 18 comparisons
                        % don't look at same state patch centers 
                        tempPatchCenter(uid,s1,s2,p1,p2,:) = (apfstats.patchCenters(uid,s1,p1,:)+apfstats.patchCenters(uid,s2,p2,:))./2;
                        tempPatchRadius(uid,s1,s2,p1,p2) = 150 + cdist(uid,s1,s2,p1,p2);
                        cdist(uid,s1,s2,p1,p2) = 0;
                        apfstats.patchCenters(uid,s1,p1,:) = tempPatchCenter(uid,s1,s2,p1,p2,:);
                        apfstats.patchCenters(uid,s2,p2,:) = tempPatchCenter(uid,s1,s2,p1,p2,:);
                    end
                end
            end
        end
    end
end

% RECOMPUTE patch distance after merger
cdist = nan([size(egoCluSessionMap,1),3,3,2,2]);
for uid = 1:size(egoCluSessionMap,1),
    t = egoCluSessionMap(uid,1);
    for s1 = 1:3
        for s2 = 1:3
            for p1  = 1:2
                for p2 = 1:2
                    cdist(uid,s1,s2,p1,p2) = sqrt(sum(sq(apfstats.patchCenters(uid,s1,p1,:) - ...
                                                         apfstats.patchCenters(uid,s2,p2,:)).^2));
                end
            end
        end
    end
end
% SELECT final patches
tempPatchCenter = nan([size(egoCluSessionMap,1),3,3,2,2,2]);
tempPatchRadius = nan([size(egoCluSessionMap,1),3,3,2,2]);
for uid = 1:size(egoCluSessionMap,1),
    for s1 = 1:3
        for s2 = 1:3
            for p1 = 1:2
                for p2 = 1:2                
                    if cdist(uid,s1,s2,p1,p2) < 150 
                        % 18 comparisons
                        % don't look at same state patch centers 
                        tempPatchCenter(uid,s1,s2,p1,p2,:) = (apfstats.patchCenters(uid,s1,p1,:)+apfstats.patchCenters(uid,s2,p2,:))./2;
                        tempPatchRadius(uid,s1,s2,p1,p2) = 150 + cdist(uid,s1,s2,p1,p2);
                    else
                        tempPatchCenter(uid,s1,s2,p1,1,:) = apfstats.patchCenters(uid,s1,p1,:);
                        tempPatchCenter(uid,s1,s2,p1,2,:) = apfstats.patchCenters(uid,s2,p2,:);
                        tempPatchRadius(uid,s1,s2,p1,p2) = 150 + cdist(uid,s1,s2,p1,p2);
                    end
                end
            end
        end
    end
end

apfstats.finalPatchCenters = nan([size(egoCluSessionMap,1),10,2]);
for uid = 1:size(egoCluSessionMap,1)
    t = egoCluSessionMap(uid,1);
    u = egoCluSessionMap(uid,2);
    [ patches, patchInd ] = unique(reshape(sq(tempPatchCenter(uid,:,:,:,:)),[],2),'rows');
    patches = patches(nniz(patches),:);
    for p1 = 1:size(patches,1)
        for p2 = p1+1:size(patches,1)
            patchDist = sqrt(sum((patches(p1,:)-patches(p2,:)).^2));
            if patchDist < 150
                patches(p1,:) = (patches(p1,:)+patches(p2,:))/2;
                patches(p2,:) = nan([1,2]);
            end
        end
    end
    
    for p =  1:size(patches,1)
        if  isnan(patches(p,1))
            continue
        end
        [~,mapindx] = NearestNeighbour(pfet{t}{s}.adata.bins{1},patches(p,1));
        [~,mapindy] = NearestNeighbour(pfet{t}{s}.adata.bins{2},patches(p,2));        
        rmap = cf(@(pp) plot(pp,u,[],[],[],false), pfet{t});
        if ~any(cellfun(@(rr) rr(mapindx,mapindy), rmap)>1)
            patches(p,:) = nan;
        end
    end
    patches(~nniz(patches),:) = [];
    apfstats.finalPatchCenters(uid,1:size(patches,1),:) = patches;
end


% $$$ 
% $$$ t = 4
% $$$ figure,
% $$$ for u = units{t};
% $$$     uid = find(ismember(egoCluSessionMap,[t,u],'rows'));
% $$$     mrate = max(apfstats.peakFR(uid,:));
% $$$     clf();
% $$$     for s = 1:3
% $$$         subplot(1,3,s);
% $$$         plot(pfet{t}{s},u,1,'colorbar',[0,mrate]);
% $$$         hold('on');
% $$$         for p = 1:2
% $$$             if apfstats.patchPFR(uid,s,p)>1
% $$$                 rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
% $$$                 rinds = rinds(nniz(rinds),:);
% $$$                 if size(rinds,1)>40
% $$$                     plot(pfet{t}{s}.adata.bins{1}(rinds(:,1)),...
% $$$                          pfet{t}{s}.adata.bins{2}(rinds(:,2)),...
% $$$                          'm*');                    
% $$$                 end        
% $$$             end
% $$$         end
% $$$         for p = 1:sum(nniz(apfstats.finalPatchCenters(uid,:,1)'))
% $$$             circle(apfstats.finalPatchCenters(uid,p,1),apfstats.finalPatchCenters(uid,p,2),175);
% $$$         end
% $$$         title(num2str(u))
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end

%% SECOND attempt
numUnits = size(egoCluSessionMap,1);
numStates = numel(pfet{1});
maxNumPatches = 2;
rthresh = 1.5;
dthresh = 100;
athresh = 9600;%16000
numPatches =numel(pfet{1})*maxNumPatches;
patchArea = reshape(apfstats.patchArea,[],numPatches);
patchPFR  = reshape(apfstats.patchPFR,[],numPatches);
patchInds = reshape(apfstats.patchRateInd,[],numPatches,2,370);
patchCntr = nan([numUnits,numPatches,2]);
for uid = 1:numUnits
    for p = 1:numPatches
        if patchArea(uid,p) < athresh  |  patchPFR(uid,p) < rthresh
            patchInds(uid,p,:,:) = nan;
            patchArea(uid,p) = nan;
            patchPFR(uid,p) = nan;
        end
    end
end

for uid = 1:numUnits
    t = egoCluSessionMap(uid,1);
    for p = 1:numPatches
        if ~isnan( patchArea(uid,p) )
            rinds = sq(patchInds(uid,p,:,:))';
            rinds = rinds(nniz(rinds),:);
            patchCntr(uid,p,:) = [mean(pfet{t}{1}.adata.bins{1}(rinds(:,1))),...
                                  mean(pfet{t}{1}.adata.bins{2}(rinds(:,2)))];
        end
    end
end

patchDist = nan([numUnits,numPatches,numPatches]);
for uid = 1:numUnits
    for p1 = 1:numPatches    
        for p2 = 1:numPatches
            patchDist(uid,p1,p2) = sqrt(sum((patchCntr(uid,p1,:)-patchCntr(uid,p2,:)).^2));
        end
    end
end



patchCntrF = nan([size(patchCntr)]);
patchCntrI = nan([size(patchCntr)]);
for uid = 1:numUnits
    [cnt,sind] = sort(sum(sq(patchDist(uid,:,:))<dthresh  &  sq(patchDist(uid,:,:))~=0),'descend');
    for p = 1:numel(cnt),
        pweights = permute(repmat(sq(patchPFR(uid,sq(patchDist(uid,sind(p),:))<dthresh)./sum(patchPFR(uid,sq(patchDist(uid,sind(p),:))<dthresh))),...
                          [2,1]),...
                           [3,2,1]);
        patchCntrI(uid,p,:) = sum(patchCntr(uid,sq(patchDist(uid,sind(p),:))<dthresh,:).*pweights,2);
    end
    patchCntrUnique = unique(sq(patchCntrI(uid,:,:)),'rows');
    patchCntrUnique(~nniz(patchCntrUnique),:) = [];
    patchCntrF(uid,1:size(patchCntrUnique,1),:) = patchCntrUnique;
end


patchCntr = patchCntrF;
patchCntrF = nan([size(patchCntr)]);
patchCntrI = nan([size(patchCntr)]);
patchDist = nan([numUnits,numPatches,numPatches]);
for uid = 1:numUnits
    for p1 = 1:numPatches    
        for p2 = 1:numPatches
            patchDist(uid,p1,p2) = sqrt(sum((patchCntr(uid,p1,:)-patchCntr(uid,p2,:)).^2));
        end
    end
end
for uid = 1:numUnits
    [cnt,sind] = sort(sum(sq(patchDist(uid,:,:))<dthresh  &  sq(patchDist(uid,:,:))~=0),'descend');
    for p = 1:numel(cnt),
        patchCntrI(uid,p,:) = mean(patchCntr(uid,sq(patchDist(uid,sind(p),:))<dthresh,:),2);
    end
    patchCntrUnique = unique(sq(patchCntrI(uid,:,:)),'rows');
    patchCntrUnique(~nniz(patchCntrUnique),:) = [];
    patchCntrF(uid,1:size(patchCntrUnique,1),:) = patchCntrUnique;
end

patchDistF = nan([numUnits,numPatches,numPatches]);
for uid = 1:numUnits
    for p1 = 1:numPatches    
        for p2 = 1:numPatches
            patchDistF(uid,p1,p2) = sqrt(sum((patchCntrF(uid,p1,:)-patchCntrF(uid,p2,:)).^2));
        end
    end
end



% DIAGNOSTIC figure
t = 20
figure,
for u = unitsEgo{t};
    uid = find(ismember(egoCluSessionMap,[t,u],'rows'));
    mrate = max(apfstats.peakFR(uid,:));
    clf();
    for s = 1:3
        subplot(1,3,s);
        plot(pfet{t}{s},u,1,'colorbar',[0,mrate],false);
        hold('on');
        for p = 1:2
            if apfstats.patchPFR(uid,s,p)>1
                rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
                rinds = rinds(nniz(rinds),:);
                if size(rinds,1)>40
                    plot(pfet{t}{s}.adata.bins{1}(rinds(:,1)),...
                         pfet{t}{s}.adata.bins{2}(rinds(:,2)),...
                         'm*');
                    circle(apfstats.patchCenters(uid,s,p,1),apfstats.patchCenters(uid,s,p,2),150);
                end        
            end
        end
        title(num2str(u))
    end
    waitforbuttonpress();
end



%%%<<< EgoProCode2D_f1_ego_maxrate_stats

rateRows = log10(egoMaxRmapRate(uidsCA1,:));
[mx,mind] = max(log10(egoMaxRmapRate(uidsCA1,:)),[],2);
EgoProCode2D_f1_ego_maxrate_stats_ca1_rows = [];
EgoProCode2D_f1_ego_maxrate_stats_ca1_rowlengths = [];
for m = 1:3 
    [~,sind] = sort(mx(mind==m));
    myrows = rateRows(mind==m,:);
    EgoProCode2D_f1_ego_maxrate_stats_ca1_rows = cat(1, EgoProCode2D_f1_ego_maxrate_stats_ca1_rows,myrows(sind,:));
    EgoProCode2D_f1_ego_maxrate_stats_ca1_rowlengths(m) = size(myrows,1);
end

rateRows = log10(egoMaxRmapRate(uidsCA3,:));
[mx,mind] = max(log10(egoMaxRmapRate(uidsCA3,:)),[],2);
EgoProCode2D_f1_ego_maxrate_stats_ca3_rows = [];
EgoProCode2D_f1_ego_maxrate_stats_ca3_rowlengths(m) = [];
for m = 1:3 
    [~,sind] = sort(mx(mind==m));
    myrows = rateRows(mind==m,:);
    EgoProCode2D_f1_ego_maxrate_stats_CA3_rows = cat(1, EgoProCode2D_f1_ego_maxrate_stats_CA3_rows,myrows(sind,:));
    EgoProCode2D_f1_ego_maxrate_stats_ca3_rowlengths(m) = size(myrows,1);
end

%%%>>>


%%%<<< EgoProCode2D_f1_ego_field_size_stats_ca1
rateRows = log10(egoSize(uidsCA1,:)*4/10000);
[mx,mind] = max(log10(egoSize(uidsCA1,:)*0.02^2),[],2);
EgoProCode2D_f1_ego_field_size_stats_ca1_rows = [];
EgoProCode2D_f1_ego_field_size_stats_ca1_rowlengths(m) = [];
for m = 1:3 
    [~,sind] = sort(mx(mind==m));
    myrows = rateRows(mind==m,:);
    EgoProCode2D_f1_ego_field_size_stats_ca1_rows = cat(1, EgoProCode2D_f1_ego_field_size_stats_ca1_rows,myrows(sind,:));
    EgoProCode2D_f1_ego_field_size_stats_ca1_rowlengths(m) = size(myrows,1);
end
