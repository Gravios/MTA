
Trial = MTATrial('jg05-20120317',[],'all')
srange = round(sum(diff(Trial.Bhv.getState('walk').state,1,2))/2);
niter = 10000;

pfwps = MTAPlaceField(Trial,[],'walk',1,'all','xy',[],srange, ...
                      niter);






load /data/homes/gravio/data/analysis/jg05-20120310/jg05-20120310.cof.all.Pfs.xy.head_front.walk.n100000bs10000sm3bn50.mat
Trial = MTATrial('jg05-20120310',[],'all')

for unit=1:size(pfwps.cluMap,1),
    if size(pfwps.rateMap{unit},3)~=pfwps.numBSiterations, 
        tpval=nan(pfwps.nbins,pfwps.nbins);,
    else
        tpval = 1./sum(repmat(max(max(pfwps.rateMap{unit})),pfwps.nbins,pfwps.nbins)...
                       <permute(repmat(permute(pfw.rateMap{unit}, ...
                                               circshift([1:ndims(pfw.rateMap{unit})+1]',1)'),pfwps.numBSiterations,1),circshift([1:ndims(pfw.rateMap{unit})+1]',-1)'),3);
        tpval(isinf(tpval))=nan;
    end
    pval(:,:,unit) = tpval;
end


figure



unit = 1;
while unit~=-1
subplot(131)
pfw.plot(unit,1)
subplot(132)
imagescnan({pfw.xbin,pfw.ybin,pval(:,:,unit)'},[],[],1,[0,0,0]);axis xy
subplot(133)
pfws = pfw.rateMap{unit};
pfws(pval(:,:,unit)>0.05|isnan(pval(:,:,unit))&~(isnan(pval(:,:,unit))&isnan(pfws(:,:)))) = 0;
imagescnan({pfw.xbin,pfw.ybin,pfws'},[],[],1,[0,0,0]);axis xy
unit = figure_controls(gcf,unit);
end