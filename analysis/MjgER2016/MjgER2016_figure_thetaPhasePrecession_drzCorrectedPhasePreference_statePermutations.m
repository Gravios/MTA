% MjgER2016_figure_thetaPhasePrecession_correctedPhasePreferenceDecomposition.m
% DRZ := Directed Rate Zone (gaussian)
% TPP := Theta Phase Preference
%
% MjgER2016_load_data(); { MjgER2016_figure_thetaPhasePrecession.m }
% Trials;                { MjgER2016_figure_thetaPhasePrecession.m }
% units;                 { MjgER2016_figure_thetaPhasePrecession.m }
% cluSessionMapSubset;   { MjgER2016_figure_thetaPhasePrecession.m }
% unitSubset;            { MjgER2016_figure_thetaPhasePrecession.m }
% stsColor = 'rgb';      { MjgER2016_figure_thetaPhasePrecession_drzCorrectedPhasePreference.m }
% statesIndsGPE;         { MjgER2016_figure_thetaPhasePrecession_drzCorrectedPhasePreference.m }


%%%<<< DEFINE default vars
perm.states = stateLabels(statesIndsGPE);
perm.sigma = sigma;
perm.statePairs = {};
perm.numIter = 100;
nSts = numel(perm.states);
perm.statePairs = cat(1,...
                      perm.states(nonzeros(repmat(1:nSts,   nSts,1).*double(~eye(nSts)))),...
                      perm.states(nonzeros(repmat([1:nSts]',1,nSts).*double(~eye(nSts)))))';
perm.stateComparisons = {{'rear','high'},{'rear','low'},{'low','high'}};
%%%>>>

%%%<<< COLLECT all pair-wise permutation rate maps
% LOAD drzphz fields
perm.ce = cf(@(s1,s2) ...
             cf(@(T,u) ...
                MTAApfs(T,u,'tag',['ddtp-','s',num2str(perm.sigma),'-',s1,'-',s2,'-1']), ...
                Trials, units), ...
             perm.statePairs(:,1),perm.statePairs(:,2));
for n = 2:perm.numIter,
    ceTemp = cf(@(s1,s2) ...
                cf(@(T,u) ...
                   MTAApfs(T,u,'tag',['ddtp-','s',num2str(perm.sigma),'-',s1,'-',s2,'-',num2str(n)]), ...
                   Trials, units), ...
                perm.statePairs(:,1),perm.statePairs(:,2));    
    for p = 1:size(perm.statePairs,1),
        for t = 1:numel(Trials),
            perm.ce{p}{t}.data.rateMap = cat(3,perm.ce{p}{t}.data.rateMap,ceTemp{p}{t}.data.rateMap);
        end
    end
end
%%%>>>

%%%<<< DECAPSULATE drzphz fields
perm.ceMat = cf(@(pfs) ...
                cf(@(p,u) ...
                   p.data.rateMap(:,ismember(p.data.clu,u),:),  ...
                   pfs,units), ...
                perm.ce);
% COLLECT cluster Ids
perm.clu     = cf(@(pfs) cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),pfs,units), perm.ce);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
perm.ceMat = cf(@(pfs) cat(2,pfs{:}), perm.ceMat);
perm.clu     = cf(@(c)   cat(2,c{:}),perm.clu);
assert(all(cell2mat(cf(@(c) all(c), cf(@(c,o) c==o, perm.clu,circshift(perm.clu,1,2))))),'Clusters do not match');
perm.ceMat = cat(4,perm.ceMat{:});
perm.clu = perm.clu{1};
tlu=cat(2,tlu{:});
perm.clu=[tlu',perm.clu'];
% SORT drzphz field maps by the unitSubset
[perm.clu,rind] = sortrows(perm.clu);
perm.clu   = perm.clu(unitSubset,:);
perm.ceMat = perm.ceMat(:,rind,:,:);
perm.ceMat = perm.ceMat(:,unitSubset,:,:);
perm.ceMat = reshape(permute(perm.ceMat,[1,5,2,3,4]),[perm.ce{1}{1}.adata.binSizes',...
                    numel(unitSubset),perm.numIter,size(perm.statePairs,1)]);
%%%>>>

%%%<<< COLLECT drz corrected rates
% FOR each unit, for each state, for each phase, find the maximum rate and position along the gdz

% GET bins
[perm.drzBins,perm.phzBins] = deal(perm.ce{1}{1}.adata.bins{:});

perm.ceMat = permute(RectFilter(permute(perm.ceMat,...
                                        [2,1,3,4,5]),...
                                3,1,'circular'),...
                     [2,1,3,4,5]);

[perm.phzRateMax,perm.phzRateMaxPosInd] = ...
    max(sq(mean(perm.ceMat,4,'omitnan')));
perm.phzRateMax = sq(perm.phzRateMax);
perm.phzRateMaxPosInd = sq(perm.phzRateMaxPosInd);

perm.ceMatSize = size(perm.ceMat);
perm.phzRateMeanPos = sq(sum(repmat(perm.drzBins,[1,perm.ceMatSize(2:5)]) ... 
                                 .*bsxfun(@rdivide,perm.ceMat,sum(perm.ceMat,'omitnan')),...
                             'omitnan'));

perm.phzRateMeanPosInd = discretize(sq(perm.phzRateMeanPos(:,:,1,:)),                      ...  Position
                               [perm.drzBins-diff(perm.drzBins(1:2));            ...
                                perm.drzBins(end)+diff(perm.drzBins(1:2))]);%Position Bins
perm.phzRateMeanPosInd = reshape(perm.phzRateMeanPosInd,[],1);
perm.phzRateMeanAll = nan(perm.ceMatSize([2,3,5,4]));
for f = 1:perm.ceMatSize(4),
    ratesm = sq(perm.ceMat(:,:,:,f,:));
    perm.phzRateMeanAll(:,:,:,f) = ...
        reshape(ratesm([0:numel(perm.phzRateMeanPosInd)-1]'.*size(ratesm,1)+perm.phzRateMeanPosInd),...
                perm.ceMatSize([2,3,5]));
end
perm.phzRateMean = mean(perm.phzRateMeanAll,4,'omitnan');
perm.phzRateStd  = std(perm.phzRateMeanAll,[],4,'omitnan');
%%%>>>

%%%<<< COMPUTE inter state permuted difference
% rear-high
% rear-low
% low-high
for s = 1:numel(perm.stateComparisons),
    sti = all(~cellfun(@isempty,regexp(perm.statePairs,...
                                      ['(',perm.stateComparisons{s}{1},')' ...
                                      '|(',perm.stateComparisons{s}{2},')'])),...
              2);
    stif = sti & ~cellfun(@isempty,regexp(perm.statePairs(:,1),...
                                          ['(',perm.stateComparisons{s}{1},')']));

    perm.phzRateMeanDiff(:, = reshape( perm.phzRateMeanAll(:,:,stif,:)...
                                   -permute(perm.phzRateMeanAll(:,:,sti&~stif,:),[1,2,4,3]),...
                                   numel(perm.phzBins),[]);
end
%dratesa = reshape(ratesa(:,9,1,:)-permute(ratesa(:,9,2,:),[1,2,4,3]),16,[]);
%%%>>>

%%%<<< COMPUTE complex value vector of phzRateMax
perm.phzRateMeanCpx = bsxfun(@times,                                                     ...
                perm.phzRateMean,                                                     ...
                exp(i.*(perm.phzBins+double(perm.phzBins<0).*2*pi)));
perm.phzRateMeanCpxMean = sq(sum(perm.phzRateMeanCpx)./sum(perm.phzRateMean,1,'omitnan'));
%%%>>>

%%%<<< GET sorted equivalents
% COMPUTE mean angle from weighted mean of the complex phase rates
perm.prmPhzAng = angle(perm.phzRateMeanCpxMean); 
[~,perm.prmPhzInd] = min(abs(bsxfun(@circ_dist,reshape(perm.prmPhzAng,1,[]),perm.phzBins)));
perm.prmPhzInd = reshape(perm.prmPhzInd,[],3);
perm.prmPhzRln = abs(perm.phzRateMeanCpxMean);
perm.prmMeanRate = sq(mean(perm.phzRateMean,'omitnan'));
[perm.prmMaxRate,perm.prmMaxRateInd] = max(phzRateMean,[],'omitnan');
[perm.prmMinRate,perm.prmMinRateInd] = min(phzRateMean,[],'omitnan');
[perm.prmMinRate,perm.prmMinRateInd,perm.prmMaxRate,perm.prmMaxRateInd] = ...
    deal(sq(perm.prmMinRate),sq(perm.prmMinRateInd),sq(perm.prmMaxRate),sq(perm.prmMaxRateInd));

perm.prmPhzAngSrt = nan(size(perm.prmMeanRate));
perm.prmPhzRlnSrt = nan(size(perm.prmMeanRate));
perm.prmMeanRateSrt   = nan(size(perm.prmMeanRate));
perm.prmMaxRateSrt    = nan(size(perm.prmMeanRate));
perm.prmMinRateSrt    = nan(size(perm.prmMeanRate));
perm.prmMaxRateIndSrt    = nan(size(perm.prmMeanRate));
perm.prmMinRateIndSrt    = nan(size(perm.prmMeanRate));
perm.prmPosSrt        = nan(size(perm.phzRateMeanPos));
perm.prmRateSrt       = nan(size(perm.phzRateMean));
for j = 1:size(perm.prmMeanRate,1);
    for s = 1:3,
        perm.prmPhzAngSrt    (j,s) = perm.prmPhzAng(j,ssi(j,s));
        perm.prmPhzRlnSrt    (j,s) = perm.prmPhzRln(j,ssi(j,s));        
        perm.prmMeanRateSrt  (j,s) = perm.prmMeanRate  (j,ssi(j,s));
        perm.prmMaxRateSrt   (j,s) = perm.prmMaxRate   (j,ssi(j,s));    
        perm.prmMinRateSrt   (j,s) = perm.prmMinRate   (j,ssi(j,s));            
        perm.prmMaxRateIndSrt(j,s) = perm.prmMaxRateInd(j,ssi(j,s));
        perm.prmMinRateIndSrt(j,s) = perm.prmMinRateInd(j,ssi(j,s));
        perm.prmPhzIndSrt    (j,s) = perm.prmPhzInd    (j,ssi(j,s));
        perm.prmPosSrt     (:,j,s) = perm.phzRateMeanPos(:,j,ssi(j,s));
        perm.prmRateSrt    (:,j,s) = perm.phzRateMean (:,j,ssi(j,s));
    end
end

perm.prmPhzAngSrtShift = perm.prmPhzAngSrt + double(perm.prmPhzAngSrt<0).*2.*pi;
%%%>>>

%%%<<< COLLECT rateMaxDcpp slope analysis 
%     THRESHOLD the rate to find the start, middle, and end of phase precession
%     REPORT the phase and position of each.
%
%     Compute start,middle,end of p,d ∀ u {u ∈ Units, p ∈ Phz, d ∈ Drz}
%     

smePhz = nan([size(prmMaxRateSrt,1),2,3]);
smePhzInd = nan([size(prmMaxRateSrt,1),2,3]);
smeDrz = nan([size(prmMaxRateSrt,1),2,3]);
errCnt = 0;
for u = stind
    for s = 1:2,
        if prmMaxRateSrt(u,s) > 4,
            rateThr = (prmMeanRateSrt(u,s)-prmMinRateSrt(u,s)).*0.5;
            pktrInds = [prmMinRateIndSrt(u,s),prmMaxRateIndSrt(u,s)];

            if pktrInds(1) > pktrInds(2),          
                indR = pktrInds(2)+1:pktrInds(1)-1;
                indL = [pktrInds(1)+1:size(prmRateSrt,1),1:pktrInds(2)-1];
            else
                indL = pktrInds(1)+1:pktrInds(2)-1;
                indR = [pktrInds(2)+1:size(prmRateSrt,1),1:pktrInds(1)-1];
            end
               

            try
            [~,pBR] = min(abs(prmRateSrt(indR,u,s)-prmMinRateSrt(u,s)-rateThr));
            [~,pBL] = min(abs(prmRateSrt(indL,u,s)-prmMinRateSrt(u,s)-rateThr));
            smePhzInd(u,s,:) = [indL(pBL),prmMaxRateIndSrt(u,s),indR(pBR)];
            smeDrz(u,s,:) = prmPosSrt(sq(smePhzInd(u,s,:)),u,s);
            smePhz(u,s,:) = phzBins(smePhzInd(u,s,:));
            catch
                errCnt = errCnt + 1;
            end

        end
    end
end

%%%>>>


