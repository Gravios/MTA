Trial = MTATrial('jg05-20120310',[],'all');
pfw = MTAPlaceField(Trial,[],'walk',0);
%pfw_ps =MTAPlaceField(Trial,[],'walk',0,'all','xy',[],100000,1000);
pfw_ps =MTAPlaceField(Trial,[],'walk',0,'all','xy',[],100000,10000);
pfr = MTAPlaceField(Trial,[],'rear',0);
pfrs = MTAPlaceField(Trial,[],'theta',0,'all','pfcrz');
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



% pval straight up
for unit = 1:sum(pfw_ps.calculation_completion_map),
    pscoreMap{unit}=zeros(length(pfw_ps.ybin),length(pfw_ps.xbin));
    if size(pfw_ps.rateMap{unit},3) > 1,
        pscoreMap{unit} = 1./sum(pfw_ps.rateMap{unit}<permute(reshape(repmat(pfw.rateMap{unit},size(pfw_ps.rateMap{unit},3),1),length(pfw_ps.ybin),size(pfw_ps.rateMap{unit},3),length(pfw_ps.xbin)),[1,3,2]),3);
    end
end


% pval multiple comp

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