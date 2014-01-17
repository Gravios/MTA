Trial = MTATrial('jg05-20120310','all');
pfw = MTAPlaceField(Trial,[],'walk',0);
%pfw_ps =MTAPlaceField(Trial,[],'walk',0,'all','xy',[],100000,1000);
%pfw_ps =MTAPlaceField(Trial,[],'walk',0,'all','xy',[],100000,10000);
%pfr = MTAPlaceField(Trial,[],'rear',0);
%pfrs = MTAPlaceField(Trial,[],'theta',0,'all','pfcrz');
pfw =MTAApfs(Trial,[],'walk',0,[],[30,30],[1.2,1.2]);
pfwl =MTAApfs(Trial,[],'lwalk',0,[],[30,30],[1.2,1.2]);
pfwh =MTAApfs(Trial,[],'hwalk',0,[],[30,30],[1.2,1.2]);

pfw_ps =MTAApfs(Trial,[7,9,6,5,20,89,25,71],'walk',0,[],[30,30],[1.2,1.2],'posShuffle',2e4,'numIter',1e4);
pfwl_ps =MTAApfs(Trial,[7,9,6,5,20,89,25,71],'lwalk',0,[],[30,30],[1.2,1.2],'posShuffle',2e4,'numIter',1e4);
pfwh_ps =MTAApfs(Trial,[7,9,6,5,20,89,25,71],'hwalk',0,[],[30,30],[1.2,1.2],'posShuffle',2e4,'numIter',1e4);

pfw_ss =MTAApfs(Trial,[7,9,6,5,20,89,25,71],'walk',0,[],[30,30],[1.2,1.2],'spkShuffle','r','numIter',1e4);
pfwl_ss =MTAApfs(Trial,[7,9,6,5,20,89,25,71],'lwalk',0,[],[30,30],[1.2,1.2],'spkShuffle','r','numIter',1e4);
pfwh_ss =MTAApfs(Trial,[7,9,6,5,20,89,25,71],'hwalk',0,[],[30,30],[1.2,1.2],'spkShuffle','r','numIter',1e4);

% $$$ 
% $$$ figure,
% $$$ unit = 1;
% $$$ while unit~=-1
% $$$ subplot2(2,1,1,1)
% $$$ pfrs.plot(unit,1)
% $$$ subplot2(2,1,1,2)
% $$$ pfw.plot(unit,1)
% $$$ title(num2str(unit))
% $$$ unit = figure_controls(gcf,unit)
% $$$ end
u = 9;
figure,
subplot2(5,3,1,1),pfw.plot(u),         subplot2(5,3,1,2),pfwl.plot(u),         subplot2(5,3,1,3),pfwh.plot(u)
subplot2(5,3,2,1),pfw_ss.plot(u,'std'),subplot2(5,3,2,2),pfwl_ss.plot(u,'std'),subplot2(5,3,2,3),pfwh_ss.plot(u,'std')
subplot2(5,3,3,1),pfw_ss.plot(u,'mean'),      subplot2(5,3,3,2),pfwl_ss.plot(u,'mean'),      subplot2(5,3,3,3),pfwh_ss.plot(u,'mean')

subplot2(5,3,4,1)
imagesc(reshape(sq(1./sum((repmat(max(pfw_ps.data.rateMap(:,pfw_ps.data.clu==u,:)),size(pfwh.data.rateMap,1),1,1)-repmat(pfw.data.rateMap(:,u,1),1,1,pfw_ps.parameters.numIter))<0,3)),33,33)'),axis xy
caxis([0,.05]),subplot2(5,3,4,2)
imagesc(reshape(sq(1./sum((repmat(max(pfwl_ps.data.rateMap(:,pfwl_ps.data.clu==u,:)),size(pfwh.data.rateMap,1),1,1)-repmat(pfwl.data.rateMap(:,u,1),1,1,pfwl_ps.parameters.numIter))<0,3)),33,33)'),axis xy
caxis([0,.05]),subplot2(5,3,4,3)
imagesc(reshape(sq(1./sum((repmat(max(pfwh_ps.data.rateMap(:,pfwh_ps.data.clu==u,:)),size(pfwh.data.rateMap,1),1,1)-repmat(pfwh.data.rateMap(:,u,1),1,1,pfwh_ps.parameters.numIter))<0,3)),33,33)'),axis xy

caxis([0,.05]),subplot2(5,3,5,1)
imagesc(reshape(sq(1./sum((repmat(max(pfw_ss.data.rateMap(:,pfw_ss.data.clu==u,:)),size(pfwh.data.rateMap,1),1,1)-repmat(pfw.data.rateMap(:,u,1),1,1,pfw_ss.parameters.numIter))<0,3)),33,33)'),axis xy
caxis([0,.05]),subplot2(5,3,5,2)
imagesc(reshape(sq(1./sum((repmat(max(pfwl_ss.data.rateMap(:,pfwl_ss.data.clu==u,:)),size(pfwh.data.rateMap,1),1,1)-repmat(pfwl.data.rateMap(:,u,1),1,1,pfwl_ss.parameters.numIter))<0,3)),33,33)'),axis xy
caxis([0,.05]),subplot2(5,3,5,3)
imagesc(reshape(sq(1./sum((repmat(max(pfwh_ss.data.rateMap(:,pfwh_ss.data.clu==u,:)),size(pfwh.data.rateMap,1),1,1)-repmat(pfwh.data.rateMap(:,u,1),1,1,pfwh_ss.parameters.numIter))<0,3)),33,33)'),axis xy
caxis([0,.05])

figure,
imagesc(reshape(sq(mean(pfw_ps.data.rateMap(:,pfw_ps.data.clu==u,:),3)./std(pfw_ps.data.rateMap(:,pfw_ps.data.clu==u,:),[],3)),33,33)'),axis xy



% pval straight up
for unit = 1:numel(pfw_ps.data.clu)
    if size(pfw_ps.rateMap{unit},3) > 1,
        pscoreMap{unit} = 1./sum(pfw_ps.rateMap{unit}<permute(reshape(repmat(pfw.rateMap{unit},size(pfw_ps.rateMap{unit},3),1),length(pfw_ps.ybin),size(pfw_ps.rateMap{unit},3),length(pfw_ps.xbin)),[1,3,2]),3);
    end
end

figure
imagesc(sq(1./sum((repmat(max(pfwl_ps.data.rateMap(:,3,:)),size(pfwl.data.rateMap,1),1,1)-repmat(pfwl.data.rateMap(:,25,1),1,1,pfwl_ps.parameters.numIter))>0,3)))

figure
imagesc(reshape(sq(1./sum((repmat(max(pfwl_ps.data.rateMap(:,1,:)),size(pfwl.data.rateMap,1),1,1)-repmat(pfwl.data.rateMap(:,20,1),1,1,pfwl_ps.parameters.numIter))<0,3)),33,33)'),axis xy
% pval multiple comp
u = 71;
imagesc(reshape(sq(1./sum((repmat(max(pfwh_ps.data.rateMap(:,pfwh_ps.data.clu==u,:)),size(pfwh.data.rateMap,1),1,1)-repmat(pfwh.data.rateMap(:,u,1),1,1,pfwh_ps.parameters.numIter))<0,3)),33,33)'),axis xy

for unit = 1:sum(pfw_ps.calculation_completion_map),
    pscoreMapMC{unit}=zeros(length(pfw_ps.ybin),length(pfw_ps.xbin));
    if size(pfw_ps.rateMap{unit},3) > 1,
        pscoreMapMC{unit} = 1./sum(repmat(max(max(pfw_ps.rateMap{unit})),length(pfw_ps.xbin),length(pfw_ps.xbin))<permute(reshape(repmat(pfw.rateMap{unit},size(pfw_ps.rateMap{unit},3),1),length(pfw_ps.ybin),size(pfw_ps.rateMap{unit},3),length(pfw_ps.xbin)),[1,3,2]),3);
    end
end



figure,
unit = 1;
while unit~=-1,
subplot2(2,1,1,1)
pfw.plot(unit,1)
subplot2(2,1,1,2)
imagesc(pscoreMap{unit}'),axis xy,colorbar,caxis([0,0.001])
title(num2str(unit))
unit = figure_controls(gcf,unit)
end


figure,
unit = 1;
while unit~=-1,
subplot2(2,1,1,1)
pfw.plot(unit,1)
subplot2(2,1,1,2)
imagesc(pscoreMapMC{unit}'),axis xy,colorbar,caxis([0,0.001])
title(num2str(unit))
unit = figure_controls(gcf,unit)
end