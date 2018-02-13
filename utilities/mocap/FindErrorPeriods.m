function [errorPeriods,hfcl,etrig] = FindErrorPeriods(Trial,varargin)
%function [errorPeriods,hbflr,hrlbf,etrig] = FindErrorPeriods(Trial,varargin)
[xyz,markers] = DefaultArgs(varargin,{'pos',{'head_back','head_left','head_front','head_right'}},1);

if ischar(xyz),
    xyz = Trial.load('xyz',xyz);
end


hfcl{1} = transform_origin(Trial,xyz,markers{1},markers{3},{markers{2},markers{4}});
hfcl{2} = transform_origin(Trial,xyz,markers{3},markers{1},{markers{4},markers{2}});
hfcl{3} = transform_origin(Trial,xyz,markers{4},markers{2},{markers{1},markers{3}});
hfcl{4} = transform_origin(Trial,xyz,markers{2},markers{4},{markers{3},markers{1}});

efet = [hfcl{1}.transVec(:,1,2),hfcl{2}.transVec(:,1,2)];
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

