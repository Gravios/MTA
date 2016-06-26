
dbstop at 107 in map_to_reference_session.m

distributions_features

dbcont

dbcont

f = 10;
vbins = linspace(-3,2,100);
% Default Font Size set to 8
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)


% Figure Preformating
hfig = figure(2016061402);clf
set(hfig,'PaperPositionMode','auto');
set(hfig,'units','centimeters');
set(hfig,'Position',[0,0,18,16])



% Axes Preformating

axes('Units','centimeters',...
     'Position',[2,6,3,3]);
imagesc(vbins,vbins,tarMean{f}');
axis xy
colormap jet
caxis([40,100]);
title('mean height target');
ylabel({'head speed','log10(cm/s)'});

axes('Units','centimeters',...
     'Position',[2,2,3,3]);
imagesc(vbins,vbins,tarStd{f}');
axis xy
colormap jet
caxis([1,15]);
title('std height target');
xlabel({'body speed','log10(cm/s)'});
ylabel({'head speed','log10(cm/s)'});

axes('Units','centimeters',...
     'Position',[6,6,3,3]);
imagesc(vbins,vbins,refMean{f}');
axis xy
colormap jet
caxis([40,100])
hcb = colorbar('manual');
set(hcb,'Units','centimeters');
set(hcb,'Position',[9.5,6,0.5,3]);
title('mean height reference');



axes('Units','centimeters',...
     'Position',[6,2,3,3]);
imagesc(vbins,vbins,refStd{f}');
axis xy
colormap jet
caxis([0.1,15])
hcb = colorbar('manual');
set(hcb,'Units','centimeters');
set(hcb,'Position',[9.5,2,0.5,3]);
title('std height reference');
xlabel({'body speed','log10(cm/s)'});

inSync = [1;features.size(1)];
nnz = nniz([tarMean{f}(:),refMean{f}(:)])...
      & tarStd{f}(:)<stdThresh{f==fetInds}...
      & tarStd{f}(:)<stdThresh{f==fetInds};
nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarStd{f}(:)<15 & refStd{f}(:)<15;
mzd = diffFun{f==fetInds}(tarMean{f}(:),refMean{f}(:));
mshift = nanmean(mzd(nnz));
dm = tarMean{f}-refMean{f};
dm(~nnz) = nan;



axes('Units','centimeters')
[~,hcb] = imagescnan({vbins,vbins,dm'},[-10,20],0,1,[0,0,0]);
axis xy
colormap jet
set(hcb,'location','EastOutside');
set(hcb,'Units','centimeters')
set(hcb,'Position',[15.5,6,0.5,3]);
set(gca,'Position',[12,6,3,3]);
title('mean height difference')


axes('Units','centimeters',...
     'Position',[12,2,3,3]);
bar(linspace(-30,40,200),histc(mzd(nnz),linspace(-30,40,200)),'histc');
xlim([-30,40])
Lines(mean(mzd(nnz)),[],'c');
title('distribution of differences')
ylabel('count')
xlabel({'height difference (mm)','target - reference'})


dbstop at 47 in mean_embeded_feature_vbvh.m

matInd = [20,30];
axes('Units','centimeters',...
     'Position',[1.5,11,3,1]);
hist(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)),200)
subDist_m = mean(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
subDist_s = std(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
Lines(subDist_m,[],'c');
Lines([subDist_m-subDist_s,subDist_m+subDist_s],[],'r');
ylim([0,30])
xlim([40,70])
xlabel('Height mm')


matInd = [40,60];
axes('Units','centimeters',...
     'Position',[1.5,12.5,3,1]);
hist(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)),200)
ylim([0,10])
xlim([40,100])
subDist_m = mean(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
subDist_s = std(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
Lines(subDist_m,[],'c');
Lines([subDist_m-subDist_s,subDist_m+subDist_s],[],'r');


matInd = [70,80];
axes('Units','centimeters',...
     'Position',[1.5,14,3,1]);
hist(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)),200)
ylim([0,10])
xlim([50,100])
subDist_m = mean(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
subDist_s = std(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
Lines(subDist_m,[],'c');
Lines([subDist_m-subDist_s,subDist_m+subDist_s],[],'r');
title('height distribution')


dbcont


matInd = [20,30];
axes('Units','centimeters',...
     'Position',[6,11,3,1]);
hist(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)),200)
ylim([0,10])
xlim([45,55])
subDist_m = mean(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
subDist_s = std(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
Lines(subDist_m,[],'c');
Lines([subDist_m-subDist_s,subDist_m+subDist_s],[],'r');
xlabel('Height mm')

matInd = [40,60];
axes('Units','centimeters',...
     'Position',[6,12.5,3,1]);
hist(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)),200)
ylim([0,10])
xlim([40,100])
subDist_m = mean(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
subDist_s = std(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
Lines(subDist_m,[],'c');
Lines([subDist_m-subDist_s,subDist_m+subDist_s],[],'r');


matInd = [70,80];
axes('Units','centimeters',...
     'Position',[6,14,3,1]);
hist(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)),200)
ylim([0,10])
xlim([40,100])
subDist_m = mean(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
subDist_s = std(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
Lines(subDist_m,[],'c');
Lines([subDist_m-subDist_s,subDist_m+subDist_s],[],'r');
title('height distribution')


print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_1',...
                     ['fig_feature_alignment.eps']))

