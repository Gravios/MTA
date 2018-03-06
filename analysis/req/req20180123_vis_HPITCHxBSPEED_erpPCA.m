
bins = pfd{1,pfindex}.adata.bins;

% helper function to reshape eigenvectors
reshape_eigen_vector = @(V,pfd) reshape(V(:,1),pfd{1}.adata.binSizes')';


% LOAD Behavioral state contours
[stateContourMaps,stateContourHandles] =                           ...
    bhv_contours(sessionListName,                                  ... sessionListName
                 'fet_HB_HPS',                                     ... featureSet
                 [1,2],                                            ... featureInd
                 {{linspace(-2,2,50),linspace(-2,2,50)}},          ... featureBin
                 'Ed05-20140529.ont.all',                          ... referenceTrial
                 {{'lloc+lpause&theta','hloc+hpause&theta',        ... states
                   'rear&theta'}},                                 ...
                 'wcr'                                             ... stateColors
);


hfig = figure(666003);clf();
hfig.Units = 'centimeters';
hfig.Position = [0.5,0.5,20,6];
hfig.PaperPositionMode = 'auto';


nV   = 5;
hax  = gobjects([1,nV]);
hcax = zeros([nV,2]);
fpc  = cell([1,nV]);
for i = 1:nV,
    fpc{i} = nan([zdims(1),1]);
    fpc{i}(validDimsInds) = V(:,i);
end
fpcMinMax = [min(cellfun(@min,fpc)),max(cellfun(@max,fpc))];

for  i = 1:nV,
    hax(i) = subplot(1,nV,i);
    imagescnan({bins{:},abs(reshape_eigen_vector(fpc{i},pfd(:,pfindex)))},... % PRINT eigenvectors
               fpcMinMax,'linear',false,[0,0,0],1,1);       
    axis('xy');
    axis('tight');
    hold('on');
    xlim([-pi/2,2]);
    ylim([-2,2]);
    for s = 1:numel(stateContourHandles),                              % OVERLAY state Contours
        copyobj(stateContourHandles{s},hax(i));
    end
end

af(@(h) set(h,'Units','centimeters'),            hax);    
af(@(h) set(h,'Position',[h.Position(1:2),1.5,1.5]), hax);

figName = ['erpPCA_HPITCHxBSPEED_v',version,'_',sessionListName];
print(hfig,'-depsc2',fullfile(figDir,[figName,'.eps']));        
print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));
