global bTar;
global bRef;
global eds;
global binnedFeatureHistTar;
global binnedFeatureHistRef;

eds = 0:100;



dbstop at 166 in map_to_reference_session.m


%distributions_features
distributions_features('featureSet','fet_head_pitch','state','gper');
distributions_features('featureSet','fet_head_pitch','state','walk','mapToReference' , true);
distributions_features('featureSet','fet_head_pitch','state', 'gper',...
                       'normalize',false, 'mapToReference' , true);
distributions_features('featureSet','fet_bref','state','walk');
% FIRST dbstop within map_to_reference_session.m
global bTar;
global bRef;
global eds;
global binnedFeatureHistTar;
global binnedFeatureHistRef;


f = 1;
isCirc = false;
vbins = linspace(-3,2,100);

clm = prctile([tarMean{f}(nniz(tarMean{f}(:)));refMean{f}(nniz(refMean{f}(:)))],[5,95]);

% Figure Preformating
hfig = figure(2016061402);clf
set(hfig,'PaperPositionMode','auto');
set(hfig,'units','centimeters');
set(hfig,'Position',[0,0,18,16])

% MEAN Target
axes('Units','centimeters',...
     'Position',[2,6,3,3]);
imagesc(vbins,vbins,tarMean{f}(:,:,10,10)');
axis xy
colormap jet
caxis(clm)
title('mean height target');
ylabel({'head speed','log10(cm/s)'});

% STD Target
axes('Units','centimeters',...
     'Position',[2,2,3,3]);
imagesc(vbins,vbins,tarStd{f}');
axis xy
colormap jet
caxis([0,abs(diff(clm))*2])
title('std height target');
xlabel({'body speed','log10(cm/s)'});
ylabel({'head speed','log10(cm/s)'});

% MEAN Reference
axes('Units','centimeters',...
     'Position',[6,6,3,3]);
imagesc(vbins,vbins,refMean{f}');
axis xy
colormap jet
caxis(clm)
hcb = colorbar('manual');
set(hcb,'Units','centimeters');
set(hcb,'Position',[9.5,6,0.5,3]);
title('mean height reference');

% STD Reference
axes('Units','centimeters',...
     'Position',[6,2,3,3]);
imagesc(vbins,vbins,refStd{f}');
axis xy
colormap jet
caxis([0,abs(diff(clm))]*2)
hcb = colorbar('manual');
set(hcb,'Units','centimeters');
set(hcb,'Position',[9.5,2,0.5,3]);
title('std height reference');
xlabel({'body speed','log10(cm/s)'});

f = 1;
inSync = [1;features.size(1)];
% $$$ nnz = nniz([tarMean{f}(:),refMean{f}(:)])...
% $$$       & tarKur{f}(:)<stdThresh{f==fetInds}...
% $$$       & refKur{f}(:)<stdThresh{f==fetInds};
%nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarKur{f}(:)<5 & refKur{f}(:)<5;
%nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarKur{f}(:)<1 & refKur{f}(:)<1;
%nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarStd{f}(:)<20 & refStd{f}(:)<10;
%nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarKur{f}(:)<20 & refKur{f}(:)<10;
%nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarStd{f}(:)<20 & refStd{f}(:)<20;
%nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarKur{f}(:)<5 & refKur{f}(:)<5& tarStd{f}(:)<.5 & refStd{f}(:)<.5;
nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarStd{f}(:)<2 & refStd{f}(:)<2& tarKur{f}(:)<10 & refKur{f}(:)<10;
nnz = nniz([tarMean{f}(:),refMean{f}(:)]) & tarStd{f}(:)<30 & refStd{f}(:)<30& tarKur{f}(:)<20 & refKur{f}(:)<20;
nnz = nniz([tarMean{f}(:),refMean{f}(:)])&tarCnt{f}(:)>minimumOccupancy*features.sampleRate&refCnt{f}(:)>minimumOccupancy*features.sampleRate;%,...

mzd = diffFun{f==fetInds}(tarMean{f}(:),refMean{f}(:));
mshift = nanmean(mzd(nnz));
dm = tarMean{f}-refMean{f};
dm(~nnz) = nan;

% $$$ 
% $$$ % VALID points
% $$$ axes('Units','centimeters')
% $$$ [~,hcb] = imagescnan({vbins,vbins,dm'},[clm(1),clm(2)],0,1,[0,0,0]);
% $$$ %[~,hcb] = imagescnan({vbins,vbins,dm'},[-0.2,.20],0,1,[0,0,0]);
% $$$ axis xy
% $$$ colormap jet
% $$$ set(hcb,'location','EastOutside');
% $$$ set(hcb,'Units','centimeters')
% $$$ set(hcb,'Position',[15.5,6,0.5,3]);
% $$$ set(gca,'Position',[12,6,3,3]);
% $$$ title('mean height difference')
% $$$  
clf
axes('Units','centimeters',...
     'Position',[5,2,5,5]);
%bar(linspace(clm(1),clm(2),200),histc(mzd(nnz),linspace(clm(1),clm(2),200)),'histc');
bar(linspace(-120,120,100),histc(dm(nnz),linspace(-120,120,100)),'histc');
%xlim([-.3,.4])
xlim([-120,120])
Lines(mean(mzd(nnz)),[],'c');
Lines(median(mzd(nnz)),[],'g');
title('distribution of differences')
ylabel('count')
xlabel({'height difference (mm)','target - reference'})


dbstop at 61 in mean_embeded_feature_vbvh.m


[tarMean,tarStd,tarKur] = mean_embeded_feature_vbvh(features,   Trial,fetInds,[],verbose);


dind = 8;
acvar = Data(mind,dind);    
    
    
matInd = [15,20];
axes('Units','centimeters',...
     'Position',[1.5,11,3,1]);
hist(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)),200)
subDist_m = mean(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
subDist_s = std(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
Lines(subDist_m,[],'c');
Lines([subDist_m-subDist_s,subDist_m+subDist_s],[],'r');
% $$$ ylim([0,30])
% $$$ xlim([40,70])
xlabel('Height mm')


matInd = [12,12];
axes('Units','centimeters',...
     'Position',[1.5,12.5,3,1]);
hist(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)),200)
%ylim([0,10])
%xlim([40,100])
subDist_m = mean(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
subDist_s = std(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
Lines(subDist_m,[],'c');
Lines([subDist_m-subDist_s,subDist_m+subDist_s],[],'r');


matInd = [5,8];
axes('Units','centimeters',...
     'Position',[1.5,14,3,1]);
hist(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)),200)
%ylim([0,10])
%xlim([50,100])
subDist_m = mean(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
subDist_s = std(acvar(all(bsxfun(@eq,manifoldIndex,matInd),2)));
Lines(subDist_m,[],'c');
Lines([subDist_m-subDist_s,subDist_m+subDist_s],[],'r');
title('height distribution')



dbstop at 61 in mean_embeded_feature_vbvh.m
[tarMean,tarStd,tarKur] = mean_embeded_feature_vbvh(features,   Trial,fetInds,[],verbose);
[refMean,refStd,refKur] = mean_embeded_feature_vbvh(rfet,    RefTrial,fetInds,[],verbose);
dind = 8;
acvar = Data(mind,dind);    


bTar = [];
for i = 1:NBINS
    for j = 1:NBINS
        bTar(:,i,j) = histc(acvar(all(bsxfun(@eq,manifoldIndex,[i,j]),2)),eds);
    end
end
binnedFeatureHistTar = reshape(bTar,size(bTar,1),[]);


bRef = [];
for i = 1:NBINS
    for j = 1:NBINS
        bRef(:,i,j) = histc(acvar(all(bsxfun(@eq,manifoldIndex,[i,j]),2)),eds);
    end
end
binnedFeatureHistRef = reshape(bRef,size(bRef,1),[]);

figure,
subplot(211)
imagesc(binnedFeatureHistTar'),caxis([0,100])
subplot(212)
imagesc(binnedFeatureHistRef'),caxis([0,100])



figure,
subplot(221);
imagesc(sq(sum(bTar,2))')
subplot(222);
imagesc(sq(sum(bTar,3))')

subplot(223);
imagesc(sq(sum(bRef,2))')
subplot(224);
imagesc(sq(sum(bRef,3))')

[U,S,V] = svd(binnedFeatureHistTar);
erpPCA.m


nind = sum(binnedFeatureHistTar,2);
[LU,LR,FSr,VT] = erpPCA(binnedFeatureHistTar(nind~=0,:)',3);
nind = sum(binnedFeatureHistRef,2);
[LUr,LRr,FSrr,VTr] = erpPCA(binnedFeatureHistRef(nind~=0,:)',3);
figure,hold on,
for i = 1:3, subplot2(2,3,1,i);plot(LR(:,i)),axis tight,end
for i = 1:3, subplot2(2,3,2,i);plot(LRr(:,i)),axis tight,end
figure,hold on,
for i = 1:3, subplot2(2,3,1,i);imagesc(reshape(FSr(:,i),25,25)');end
for i = 1:3, subplot2(2,3,2,i);imagesc(reshape(FSrr(:,i),25,25)');end




[Ur,Sr,Vr] = svd(binnedFeatureHistRef);
dbcont

figure,hold on,
for i = 1:3, subplot2(2,3,1,i);imagesc(reshape(V(:,i),25,25)');end
for i = 1:3, subplot2(2,3,2,i);imagesc(reshape(Vr(:,i),25,25)');end

figure,hold on,
for i = 1:3, subplot2(2,3,1,i);plot(U(:,i)),axis tight,end
for i = 1:3, subplot2(2,3,2,i);plot(Ur(:,i)),axis tight,end

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

