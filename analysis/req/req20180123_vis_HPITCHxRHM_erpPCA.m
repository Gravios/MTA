
bins = pfd{1,pfindex}.adata.bins;

% helper function to reshape eigenvectors
reshape_eigen_vector = @(V,pfd) reshape(V(:,1),pfd{1}.adata.binSizes')';

% LOAD Behavioral state contours
[stateContourMaps,stateContourHandles] =                           ...
    bhv_contours(sessionListName,                                  ... sessionListName
                 'fet_HB_HPR',                                     ... featureSet
                 [1,2],                                            ... featureInd
                 {linspace(-2,2,50),linspace(-9,-3,50)},           ... featureBin
                 'Ed05-20140529.ont.all',                          ... referenceTrial
                 {'lloc&theta','lpause&theta','hloc&theta','hpause&theta'},... states
                 'bcrm'                                            ... stateColors
);

hfig = figure(666008);clf();
hfig.Units = 'centimeters';
hfig.Position = [0.5,0.5,24,6];
hfig.PaperPositionMode = 'auto';


hax  = gobjects([1,numComp+1]);
hcax = zeros([numComp,2]);
fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan(size(validDims{pfindex}));
    fpc{i}(validDims{pfindex}) = eigVec{pfindex}(:,i);
end

fpcMinMax = [min(cellfun(@min,fpc)),max(cellfun(@max,fpc))];

for i = 1:numComp,
    hax(i) = subplot(1,numComp+1,i);
    imagescnan({bins{:},abs(reshape_eigen_vector(fpc{i},pfd(1,pfindex)))},...
               fpcMinMax,'linear',false,[0,0,0],1,1);                % PRINT eigenvectors
    colorbar();                                           
    axis('xy');
    axis('tight');
    hold('on');
    for s = 1:numel(stateContourHandles),                            % OVERLAY state Contours
        copyobj(stateContourHandles{s},hax(i));
    end
    xlim([-2,pi/2]);
    ylim([-9,-3]);
end
hax(numComp+1) = subplot(1,numComp+1,numComp+1);
plot(eigVar{pfindex}(1:numComp,4),'-+');
xlim([0,numComp+1]);

af(@(h) set(h,'Units','centimeters'),                hax);    
af(@(h) set(h,'Position',[h.Position(1:2),1.5,1.5]), hax);


figName = ['erpPCA_HPITCHxRHM_v',version,'_',sessionListName];
print(hfig,'-depsc2',fullfile(figDir,analDir,[figName,'.eps']));        
print(hfig,'-dpng',  fullfile(figDir,analDir,[figName,'.png']));
