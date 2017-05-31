% req20170524 ---------------------------------------------------------------
%
%  Status: active
%  Type: Optimization
%  Description: Variability of marker placement between sessions
%               necessary for c
%               Req 1.0: stats between neural network classifier (NNC) labels
%               Req 1.1: stats between expert labeler classifier (ELC) labels 
%               Req 1.2: stats between expert labeler classifier (ELC) labels 
%                        and neural network classifier (NNC) labels 
%  Current Bugs: 
%  Final_Forms: NA

 



%diagnostic_MTADfet_map_to_reference_session_20170524.m


global bTar;
global bRef;
global eds;
global binnedFeatureHistTar;
global binnedFeatureHistRef;

eds = 0:100;



dbstop at 160 in map_to_reference_session.m


distributions_features
% FIRST dbstop within map_to_reference_session.m
global bTar;
global bRef;
global eds;
global binnedFeatureHistTar;
global binnedFeatureHistRef;


f = 8;
isCirc = false;
vbins = linspace(-3,2,100);

% Figure Preformating
hfig = figure(2016061402);clf
set(hfig,'PaperPositionMode','auto');
set(hfig,'units','centimeters');
set(hfig,'Position',[0,0,18,16])

% MEAN Target
axes('Units','centimeters',...
     'Position',[2,6,3,3]);
figure,
imagesc(reshape(tarMean{f},25^2,25^2)');
%imagesc(vbins,vbins,tarMean{f}');
axis xy
colormap jet
if isCirc,caxis([-pi/2,pi/2]),else,caxis([40,100]);end
title('mean height target');
ylabel({'head speed','log10(cm/s)'});

% STD Target
axes('Units','centimeters',...
     'Position',[2,2,3,3]);
imagesc(reshape(tarStd{f},25^2,25^2)');
%imagesc(vbins,vbins,tarStd{f}');
axis xy
colormap jet
if isCirc,caxis([0,0.2]),else,caxis([0,15]);end
title('std height target');
xlabel({'body speed','log10(cm/s)'});
ylabel({'head speed','log10(cm/s)'});

% MEAN Reference
axes('Units','centimeters',...
     'Position',[6,6,3,3]);
imagesc(vbins,vbins,refMean{f}');
axis xy
colormap jet
if isCirc,caxis([-pi/2,pi/2]),else,caxis([40,100]);end
hcb = colorbar('manual');
set(hcb,'Units','centimeters');
set(hcb,'Position',[9.5,6,0.5,3]);
title('mean height reference');

% MEAN Reference
axes('Units','centimeters',...
     'Position',[6,2,3,3]);
imagesc(vbins,vbins,refStd{f}');
axis xy
colormap jet
if isCirc,caxis([0,0.2]),else,caxis([0,15]);end
hcb = colorbar('manual');
set(hcb,'Units','centimeters');
set(hcb,'Position',[9.5,2,0.5,3]);
title('std height reference');
xlabel({'body speed','log10(cm/s)'});


inSync = [1;features.size(1)];
% $$$ nnz = nniz([tarMean{f}(:),refMean{f}(:)])...
% $$$       & tarKur{f}(:)<stdThresh{f==fetInds}...
% $$$       & refKur{f}(:)<stdThresh{f==fetInds};
%nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarKur{f}(:)<5 & refKur{f}(:)<5;
%nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarKur{f}(:)<1 & refKur{f}(:)<1;
nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarStd{f}(:)<10 & refStd{f}(:)<10;
nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarKur{f}(:)<10 & refKur{f}(:)<10;
mzd = diffFun{f==fetInds}(tarMean{f}(:),refMean{f}(:));
mshift = nanmean(mzd(nnz));
dm = tarMean{f}-refMean{f};
dm(~nnz) = nan;

figure,
% VALID points
axes('Units','centimeters')
[~,hcb] = imagescnan({vbins,vbins,dm'},[-10,20],0,1,[0,0,0]);
%[~,hcb] = imagescnan({vbins,vbins,dm'},[-0.2,.20],0,1,[0,0,0]);
axis xy
colormap jet
set(hcb,'location','EastOutside');
set(hcb,'Units','centimeters')
set(hcb,'Position',[15.5,6,0.5,3]);
set(gca,'Position',[12,6,3,3]);
title('mean height difference')
 
% 
axes('Units','centimeters',...
     'Position',[12,2,3,3]);
bar(linspace(-30,40,200),histc(mzd(nnz),linspace(-30,40,200)),'histc');
%bar(linspace(-.4,.4,200),histc(mzd(nnz),linspace(-.4,.4,200)),'histc');
%xlim([-.3,.4])
xlim([-30,40])
Lines(mean(mzd(nnz)),[],'c');
title('distribution of differences')
ylabel('count')
xlabel({'height difference (mm)','target - reference'})

