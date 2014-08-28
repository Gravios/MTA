function [errorPeriods,hbflr,hrlbf,etrig] = FindErrorPeriods(Trial)
xyz = Trial.xyz.copy;
xyz.load(Trial);

hbflr = Trial.transformOrigin(xyz,'head_back','head_front',{'head_left','head_right'});
hrlbf = Trial.transformOrigin(xyz,'head_right','head_left',{'head_back','head_front'});

efet = [hbflr.transVec(:,:,2),hrlbf.transVec(:,:,2)];
efmean = zeros(size(efet,2),1);
for i = 1:size(efet,2),
    efmean(i) = mean(efet(~isnan(efet(:,i)),i));
end
etrig =sum(abs(efet-repmat(efmean',size(efet,1),1)),2);

ets = std(etrig(~isnan(etrig)));
etm = mean(etrig(~isnan(etrig)));
errorInd = find(abs((etrig-etm)./ets)>1);

errorState = zeros(size(etrig));
errorState(errorInd)=1;

errorPeriods = ThreshCross(errorState,0.5,0);

