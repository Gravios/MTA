% MjgER2016_figure_thetaPhasePrecession_correctedPhasePreferenceDecomposition_statePermutations.m
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
% 
% drzphz ratemaps are computed in req20190527_perm.m

%%%<<< DEFINE default vars

perm.states = statesTPP;
perm.sigma = sigma;
perm.statePairs = {};
perm.numIter = 101;
nSts = numel(perm.states);
% GET all pair-wise state permutations
perm.statePairs = cat(1,...
                      perm.states(nonzeros(repmat(1:nSts,   nSts,1).*double(~eye(nSts)))),...
                      perm.states(nonzeros(repmat([1:nSts]',1,nSts).*double(~eye(nSts)))))';
perm.stateComparisons = {{'rear','high'},{'rear','low'},{'low','high'}};

%%%>>>

%%%<<< COLLECT all pair-wise permutation rate maps

% LOAD drzphz fields

perm.ce = cf(@(s1,s2) ...
             cf(@(T,u) ...
                MTAApfs(T,u,'tag',['ddtp2-','s',num2str(perm.sigma),'-',s1,'-',s2,'-1']), ...
                Trials, units), ...
             perm.statePairs(:,1),perm.statePairs(:,2));
for n = 2:perm.numIter,
    ceTemp = {};
    for p = 1:size(perm.statePairs,1)
        for tt = 1:numel(Trials)
            ceTemp{p}{tt} = MTAApfs(Trials{tt},units{tt},...
                    'tag',['ddtp2-s150-',perm.statePairs{p,1},'-',perm.statePairs{p,2},'-',num2str(n)]);
        end
    end
    

% $$$     ceTemp = cf(@(s1,s2) ...
% $$$                 cf(@(T,u) ...
% $$$                    MTAApfs(T,u,'tag',['ddtp2-','s',num2str(perm.sigma),'-',s1,'-',s2,'-',num2str(n)]), ...
% $$$                    Trials, units), ...
% $$$                 perm.statePairs(:,1),perm.statePairs(:,2));    
    for p = 1:size(perm.statePairs,1)
        for t = 1:numel(Trials)
            perm.ce{p}{t}.data.rateMap = cat(3,perm.ce{p}{t}.data.rateMap,ceTemp{p}{t}.data.rateMap);
        end
    end
end



%%%>>>

%%%<<< collection diagnostics
% $$$ for t = 1:23;
% $$$     ceTemp = cf(@(s1,s2) ...
% $$$                 cf(@(T,u) ...
% $$$                    MTAApfs(T,u,'tag',['ddtp-','s',num2str(perm.sigma),'-',s1,'-',s2,'-',num2str(n)]), ...
% $$$                    Trials(t), units(t)), ...
% $$$                 perm.statePairs(:,1),perm.statePairs(:,2));    
% $$$ end
% $$$ t = 19;
% $$$ n = 81;
% $$$ n = 65;
% $$$ n = 100;
% $$$ ppp = {};
% $$$ for s = 1:6,
% $$$ ppp{s} = MTAApfs(Trials{t},units{t},...
% $$$               'tag',['ddtp-','s',num2str(perm.sigma),'-',perm.statePairs{s,1},'-',perm.statePairs{s,2},'-',num2str(n)]);
% $$$ end
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

% $$$ perm.ceMat = permute(RectFilter(permute(perm.ceMat,...
% $$$                                         [2,1,3,4,5]),...
% $$$                                 3,1,'circular'),...
% $$$                      [2,1,3,4,5]);


for u = 1:size(perm.ceMat,3),
    for n = 1:size(perm.ceMat,4),
        for s = 1:size(perm.ceMat,5),
            for p = 1:size(perm.ceMat,2),
                nind = find(isnan(perm.ceMat(:,p,u,n,s)));
                gind = nind( nind<30 | nind>10 );
                if numel(gind)>4, continue; end
                bind = ~ismember(perm.drzBins,perm.drzBins(nind));
                if ~isempty(gind)
                    perm.ceMat(gind,p,u,n,s) =interp1(perm.drzBins(bind),          ...
                                                      perm.ceMat(bind,p,u,n,s),    ...
                                                      perm.drzBins(gind));
                end
            end
        end
    end
end



% $$$ phzCorrection = [6,6,              ... er01
% $$$                  4,4,4,            ... ER06
% $$$                  6,6,              ... Ed10
% $$$                  0,0,0,0,0,0,0,0,0,... jg04                 
% $$$                  2,2,2,2,2,2,2,2,2 ... jg05
% $$$                  6,6,6]; % new units - ER06, Ed10, er01
% $$$ 
% $$$ 
% $$$ for p = 1:numel(phzCorrection)
% $$$     perm.ceMat(:,:,cluSessionMapSubset(:,1)==p,:,:) = circshift(perm.ceMat(:,:,cluSessionMapSubset(:,1)==p,:,:),phzCorrection(p),2);
% $$$ end    



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

%%%<<< COMPUTE complex value vector of phzRateMax

perm.phzRateMeanCpx = bsxfun(@times, perm.phzRateMean, exp(i.*perm.phzBins));
perm.phzRateMeanCpxMean = sq(sum(perm.phzRateMeanCpx)./sum(perm.phzRateMean,1,'omitnan'));

%%%>>>

%%%<<< GET sorted equivalents

% COMPUTE mean angle from weighted mean of the complex phase rates
perm.prmPhzAng = angle(perm.phzRateMeanCpxMean); 
[~,perm.prmPhzInd] = min(abs(bsxfun(@circ_dist,reshape(perm.prmPhzAng,1,[]),perm.phzBins)));
perm.prmPhzInd = reshape(perm.prmPhzInd,[],3);
perm.prmPhzRln = abs(perm.phzRateMeanCpxMean);
perm.prmMeanRate = sq(mean(perm.phzRateMean,'omitnan'));
[perm.prmMaxRate,perm.prmMaxRateInd] = max(perm.phzRateMean,[],'omitnan');
[perm.prmMinRate,perm.prmMinRateInd] = min(perm.phzRateMean,[],'omitnan');
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

perm.smePhz = nan([size(perm.prmMaxRateSrt,1),2,3]);
perm.smePhzInd = nan([size(perm.prmMaxRateSrt,1),2,3]);
perm.smeDrz = nan([size(perm.prmMaxRateSrt,1),2,3]);
errCnt = 0;
for u = stind
    for s = 1:2,
        if perm.prmMaxRateSrt(u,s) > 4,
            rateThr = (perm.prmMeanRateSrt(u,s)-perm.prmMinRateSrt(u,s)).*0.5;
            pktrInds = [perm.prmMinRateIndSrt(u,s),perm.prmMaxRateIndSrt(u,s)];

            if pktrInds(1) > pktrInds(2),          
                indR = pktrInds(2)+1:pktrInds(1)-1;
                indL = [pktrInds(1)+1:size(perm.prmRateSrt,1),1:pktrInds(2)-1];
            else
                indL = pktrInds(1)+1:pktrInds(2)-1;
                indR = [pktrInds(2)+1:size(perm.prmRateSrt,1),1:pktrInds(1)-1];
            end
               

            try
            [~,pBR] = min(abs(perm.prmRateSrt(indR,u,s)-perm.prmMinRateSrt(u,s)-rateThr));
            [~,pBL] = min(abs(perm.prmRateSrt(indL,u,s)-perm.prmMinRateSrt(u,s)-rateThr));
            perm.smePhzInd(u,s,:) = [indL(pBL),perm.prmMaxRateIndSrt(u,s),indR(pBR)];
            perm.smeDrz(u,s,:) = perm.prmPosSrt(sq(smePhzInd(u,s,:)),u,s);
            perm.smePhz(u,s,:) = perm.phzBins(perm.smePhzInd(u,s,:));
            catch
                errCnt = errCnt + 1;
            end

        end
    end
end

%%%>>>

%%%<<< COMPUTE inter state permuted difference

% rear-high
% rear-low
% low-high              
perm.phzRateMeanDiff = zeros([size(perm.phzRateMeanAll,1),...
                              size(perm.phzRateMeanAll,2),...
                              (size(perm.phzRateMeanAll,4).^2-size(perm.phzRateMeanAll,4))./2,...
                              size(perm.statePairs,1)]);
tind = size(perm.phzRateMeanAll,4) >= rot90(bsxfun(@plus,...
                                                   1:size(perm.phzRateMeanAll,4),...
                                                   [1:size(perm.phzRateMeanAll,4)]')...
                                           );
for s = 1:size(perm.statePairs,1)
    sti = all(~cellfun(@isempty,regexp(perm.statePairs,...
                                      ['(',perm.statePairs{s,1},')' ...
                                      '|(',perm.statePairs{s,2},')'])),...
              2);
    stif = sti & ~cellfun(@isempty,regexp(perm.statePairs(:,1),...
                                          ['(',perm.statePairs{s,1},')']));

    tmat= reshape(perm.phzRateMeanAll(:,:,stif,:)...
                      -permute(perm.phzRateMeanAll(:,:,sti&~stif,:),[1,2,4,3]),...
                  numel(perm.phzBins),size(perm.phzRateMeanAll,2),[]);
    perm.phzRateMeanDiff(:,:,:,s) = tmat(:,:,tind(:));
end



%dratesa = reshape(ratesa(:,9,1,:)-permute(ratesa(:,9,2,:),[1,2,4,3]),16,[]);

%%%>>>

%%%<<< COMPUTE inter state difference

% rear-high
% rear-low
% low-high              
tpp.phzRateMeanDiff = zeros([size(tpp.phzRateMeanAll,1),...
                             size(tpp.phzRateMeanAll,2),...
                             (size(tpp.phzRateMeanAll,4).^2+size(tpp.phzRateMeanAll,4))./2,...
                             size(perm.statePairs,1)]);
tind = size(tpp.phzRateMeanAll,4) >= rot90(bsxfun(@plus,...
                                                  1:size(tpp.phzRateMeanAll,4),...
                                                  [1:size(tpp.phzRateMeanAll,4)]')-1 ...
                                           );
for s = 1:size(perm.statePairs,1)
    sti = ~cellfun(@isempty,regexp(statesTPP,['(',perm.statePairs{s,1},')']));
    stif = ~cellfun(@isempty,regexp(statesTPP,['(',perm.statePairs{s,2},')']));    
    tmat = reshape(tpp.phzRateMeanAll(:,:,sti,:)...
                       -permute(tpp.phzRateMeanAll(:,:,stif,:),[1,2,4,3]),...
                   numel(tpp.phzBins),size(tpp.phzRateMeanAll,2),[]);
    tpp.phzRateMeanDiff(:,:,:,s) = tmat(:,:,tind(:));
end

%%%>>>

%%%<<< COMPUTE complex value vector of zscore
tpp.zscore = (tpp.phzRateMeanDiff-mean(perm.phzRateMeanDiff,3,'omitnan'))  ...
                 ./std(perm.phzRateMeanDiff,[],3,'omitnan');

zScrPos = tpp.zscore;
zScrPos(zScrPos<0) = 0;
zScrPosCpx = bsxfun(@times,                                                      ...
                    zScrPos,                                                     ...
                    exp(i.*(perm.phzBins+double(perm.phzBins<0).*2*pi)));
perm.zScrPosCpxMean = sq(sum(zScrPosCpx)./sum(zScrPos,1,'omitnan'));

zScrNeg = -tpp.zscore;
zScrNeg(zScrNeg<0) = 0;
zScrNegCpx = bsxfun(@times,                                                      ...
                    zScrNeg,                                                     ...
                    exp(i.*(perm.phzBins+double(perm.phzBins<0).*2*pi)));
perm.zScrNegCpxMean = sq(sum(zScrNegCpx)./sum(zScrNeg,1,'omitnan'));

zScrShf = bsxfun(@plus,tpp.zscore,abs(min(tpp.zscore)));
zScrShf(zScrShf<0) = 0;
zScrShfCpx = bsxfun(@times,                                                      ...
                    zScrShf,                                                     ...
                    exp(i.*(perm.phzBins+double(perm.phzBins<0).*2*pi)));
perm.zScrShfCpxMean = sq(sum(zScrShfCpx)./sum(zScrShf,1,'omitnan'));





%%%>>>


% $$$ spb = perm.phzBins;
% $$$ 
% $$$ u = ismember(cluSessionMapSubset,[5,31],'rows');
% $$$ 
% $$$ u = ismember(cluSessionMapSubset,[21,22],'rows');
% $$$ figure();
% $$$ subplot(211);
% $$$     hold('on');
% $$$     plot(spb,sq( tpp.phzRateMeanDiff(phzOrder,u,1:150,4)),'r');
% $$$     plot(spb,sq(perm.phzRateMeanDiff(phzOrder,u,1:150,4)),'b');
% $$$     title('hloc-lloc (red)   hlperm (blue)')
% $$$     ylabel('rate difference');
% $$$     xlim([0,2*pi]);
% $$$ subplot(212);
% $$$     hold('on');
% $$$     plot(spb,sq(tpp.zscore(phzOrder,u,1:150,4)),'c');
% $$$     plot(spb,sq(tpp.zscore(phzOrder,u,1,4)),'k','LineWidth',5);
% $$$     ylabel('z-score')
% $$$     xlabel('Theta phase');
% $$$     xlim([0,2*pi]);
% $$$ 
% $$$     
% $$$     
% $$$     
% $$$ for p = 1:16;
% $$$ [hsig(p),pval(p)] =ttest2(sq(tpp.phzRateMeanDiff(p,u,:,4)),sq(perm.phzRateMeanDiff(p,u,:,4)));
% $$$ end
% $$$ 
% $$$ figure,plot(sq(tmap(:,u,tind(:))))
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,
% $$$ hold('on'),
% $$$ plot(sq(tpp.phzRateMeanDiff(:,u,:,4)),'r');
% $$$ plot(reshape(tpp.phzRateMeanAll(:,u,2,:)...
% $$$                             -permute(tpp.phzRateMeanAll(:,u,3,:),[1,2,4,3]),...
% $$$                                    numel(tpp.phzBins),[]))
% $$$ 
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ plot(tpp.phzRateMean(:,u,3,1),'b')
% $$$ plot(tpp.phzRateMean(:,u,2,1),'g');
% $$$ 
% $$$ 
% $$$ figure,
% $$$ hold('on');
% $$$ plot(reshape(tpp.phzRateMeanAll(:,u,2,:)-permute(tpp.phzRateMeanAll(:,u,3,:),[1,2,4,3]),...
% $$$              numel(tpp.phzBins),[]),'r')
% $$$ plot(reshape(perm.phzRateMeanAll(:,u,4,:)-permute(perm.phzRateMeanAll(:,u,6,:),[1,2,4,3]),...
% $$$              numel(tpp.phzBins),[]),'b')
% $$$ 
% $$$ 
% $$$ stind = find(validPosOccRlx & sesInds(unitSubset) & prmMaxRateSrt(:,2) > 4 & validUnits)';
% $$$ 
% $$$ set(0,'defaultAxesFontSize',10);
% $$$ 
% $$$ figure();
% $$$ for ind = 306;%stind(150:end),
% $$$     clf();
% $$$     t = cluSessionMapSubset(ind,1);
% $$$     u = cluSessionMapSubset(ind,2);
% $$$     mrateHZTPD = max(cell2mat(cf(@(p) p{t}.maxRate(u,false), pftHZTPD(statesIndsGPE(:)))));    
% $$$     subplot(241);
% $$$         plot(pfts{t},u,'mean','text',[],false,[],[],[],@jet,[],nanColor);
% $$$     for s = 1:3,
% $$$         sp2 = subplot(4,4,s+1);
% $$$         plot(pftHZTPD{statesIndsGPE(s)}{t},u,'mean','',[0,mrateHZTPD],false,[],[],[],@jet,[],nanColor); 
% $$$         sp1 = subplot(4,4,s+1+4);
% $$$         plot(pftHZTPD{statesIndsGPE(s)}{t},u,'mean','',[0,mrateHZTPD],false,[],[],[],@jet,[],nanColor); 
% $$$         sp1.Position(2) =         sp1.Position(2)+0.05;
% $$$     end
% $$$     subplot(245);
% $$$         hold('on');        
% $$$         for s = 1:3,
% $$$             plot(sq(tpp.phzRateMeanAll(phzOrder,ind,s,1:101)),sclr(s));
% $$$         end    
% $$$         title(num2str([cluSessionMapSubset(ind,:),find(stind==ind), ind]))
% $$$     subplot(2,4,6);
% $$$         hold('on')
% $$$         plot(spb,sq( tpp.zscore(phzOrder,ind,1:150,1)),'r');
% $$$         ylabel('z-score')        
% $$$         xlabel('theta phase')                
% $$$         title('rear vs high')
% $$$         xlim([0,2*pi])
% $$$         ylim([-10,10]);
% $$$     subplot(2,4,7);
% $$$         hold('on')
% $$$         plot(spb,sq( tpp.zscore(phzOrder,ind,1:150,4)),'g');
% $$$         ylabel('z-score')        
% $$$         xlabel('theta phase')                
% $$$         title('high vs low')        
% $$$         xlim([0,2*pi])        
% $$$         ylim([-10,10]);        
% $$$     subplot(2,4,8);
% $$$         hold('on')
% $$$         plot(sq( tpp.zscore(phzOrder,ind,1:150,5)),'b');
% $$$         ylabel('z-score')        
% $$$         xlabel('theta phase')                
% $$$         title('low vs rear')                
% $$$         xlim([0,2*pi])        
% $$$         ylim([-10,10]);        
% $$$     waitforbuttonpress();
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==2 & ssi(:,2)==3 )';
% $$$ 
% $$$ 
% $$$ % mean phz vs res
% $$$ figure();
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==3 & ssi(:,2)==2 )';
% $$$ plot(angle(perm.zScrPosCpxMean(stind,1,6)),abs(perm.zScrPosCpxMean(stind,1,6)),'.b','MarkerSize',15);
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==2 & ssi(:,2)==3 )';
% $$$ hold('on');
% $$$ plot(angle(perm.zScrPosCpxMean(stind,1,4)),abs(perm.zScrPosCpxMean(stind,1,4)),'.g','MarkerSize',15);
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==1 & ssi(:,2)==2 )';
% $$$ hold('on');
% $$$ plot(angle(perm.zScrPosCpxMean(stind,1,1)),abs(perm.zScrPosCpxMean(stind,1,1)),'.r','MarkerSize',15);
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure();
% $$$ subplot(325);
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==3 & ssi(:,2)==2 )';
% $$$ plot(angle(perm.zScrPosCpxMean(stind,1,6)),tpp.prmPhzAngSrt(stind,3),'.');
% $$$ title('LOW, high') 
% $$$ subplot(323);
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==2 & ssi(:,2)==3 )';
% $$$ hold('on');
% $$$ histogram(angle(perm.zScrPosCpxMean(stind,1,4)),linspace(-pi,pi,17));
% $$$ title('HIGH, low') 
% $$$ subplot(321);
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==1 & ssi(:,2)==2 )';
% $$$ hold('on');
% $$$ histogram(angle(perm.zScrPosCpxMean(stind,1,1)),linspace(-pi,pi,17));
% $$$ title('REAR, high') 
% $$$ 
% $$$ 
% $$$ figure();
% $$$     stind = find(    validPosOccRlx             ...
% $$$                    & sesInds(unitSubset)        ...
% $$$                    & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$                    & validUnits                 ...
% $$$                    & ssi(:,1)==3 & ssi(:,2)==2 )';                         
% $$$ %                   & (ssi(:,1)==3 | ssi(:,2)==2) )';                             
% $$$     %plot(max(tpp.zscore(:,stind,1,6)),std(tpp.zscore(:,stind,1,6)),'.')
% $$$     %plot(sum(tpp.zscore(:,stind,1,6)>3.1)+randn([1,numel(stind)])/5,mean(tpp.zscore(:,stind,1,6)),'.');
% $$$ subplot(121);
% $$$     plot(abs(perm.zScrShfCpxMean(stind,1,6)),sum(tpp.zscore(:,stind,1,6)>3.1)+randn([1,numel(stind)])/5,'.');
% $$$ subplot(122);    
% $$$     plot(abs(perm.zScrShfCpxMean(stind,1,6)),max(tpp.zscore(:,stind,1,6)),'.');
% $$$ 
% $$$ % phase of positive zscore of units which fire secondary state
% $$$ figure();
% $$$ subplot(313);
% $$$     hold('on')
% $$$     stind = find(  validPosOccRlx             ...
% $$$                    & sesInds(unitSubset)        ...
% $$$                    & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$                    & validUnits                 ...
% $$$                    & ssi(:,1)==3 & ssi(:,2)==2 )';
% $$$     histogram(angle(perm.zScrShfCpxMean(stind,1,6)),linspace(-pi,pi,17));
% $$$     stind = find(  validPosOccRlx             ...
% $$$                    & sesInds(unitSubset)        ...
% $$$                    & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$                    & validUnits                 ...
% $$$                    & ssi(:,1)==3 & ssi(:,2)==1 )';
% $$$     histogram(angle(perm.zScrShfCpxMean(stind,1,5)),linspace(-pi,pi,17));
% $$$     title('LOW, high,rear') 
% $$$ 
% $$$ subplot(312);
% $$$     stind = find(  validPosOccRlx             ...
% $$$                    & sesInds(unitSubset)        ...
% $$$                    & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$                    & validUnits                 ...
% $$$                    & ssi(:,1)==2 & ssi(:,2)==3 )';
% $$$     hold('on');
% $$$     histogram(angle(perm.zScrShfCpxMean(stind,1,4)),linspace(-pi,pi,17));
% $$$     stind = find(  validPosOccRlx             ...
% $$$                    & sesInds(unitSubset)        ...
% $$$                    & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$                    & validUnits                 ...
% $$$                    & ssi(:,1)==2 & ssi(:,2)==1 )';
% $$$     hold('on');
% $$$     histogram(angle(perm.zScrShfCpxMean(stind,1,3)),linspace(-pi,pi,17));
% $$$     title('HIGH, low','rear') 
% $$$ 
% $$$ subplot(311);
% $$$     stind = find(  validPosOccRlx             ...
% $$$                    & sesInds(unitSubset)        ...
% $$$                    & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$                    & validUnits                 ...
% $$$                    & ssi(:,1)==1 & ssi(:,2)==2 )';
% $$$     hold('on');
% $$$     histogram(angle(perm.zScrShfCpxMean(stind,1,1)),linspace(-pi,pi,17));
% $$$     stind = find(  validPosOccRlx             ...
% $$$                    & sesInds(unitSubset)        ...
% $$$                    & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$                    & validUnits                 ...
% $$$                    & ssi(:,1)==1 & ssi(:,2)==3 )';
% $$$     hold('on');
% $$$     histogram(angle(perm.zScrShfCpxMean(stind,1,2)),linspace(-pi,pi,17));
% $$$     title('REAR, high, low') 
% $$$ 
% $$$ ForAllSubplots('ylim([0,10]);')
% $$$ ForAllSubplots('xlim([-pi,pi]);')
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) <= 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==s)';
% $$$ histogram(tpp.prmPhzAngSrt(stind,s),linspace(-pi,pi,17));
% $$$ title('LOW') ;
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure();
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==2 & ssi(:,2)==3 )';
% $$$ plot(angle(perm.zScrPosCpxMean(stind,1,4)),...
% $$$      angle(perm.zScrNegCpxMean(stind,1,4)),...
% $$$      '.g','MarkerSize',15);
% $$$ hold('on');
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==3 & ssi(:,2)==2 )';
% $$$ plot(angle(perm.zScrPosCpxMean(stind,1,6)),...
% $$$      angle(perm.zScrNegCpxMean(stind,1,6)),...
% $$$      '.b','MarkerSize',15);
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==1 & ssi(:,2)==2 )';
% $$$ plot(angle(perm.zScrPosCpxMean(stind,1,1)),...
% $$$      angle(perm.zScrNegCpxMean(stind,1,1)),...
% $$$      '.r','MarkerSize',15);
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$  
% $$$ % Sort the z-scores 
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==2 & ssi(:,2)==3 )';
% $$$ figure();
% $$$ imagesc(tpp.zscore(phzOrder,stind,1,4)');
% $$$ caxis([-5.2,5.2]);
% $$$ 
% $$$ 
% $$$ % High Dominant -> Low
% $$$ figure();
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==2 & ssi(:,2)==3 )';
% $$$ [~,sind] = sortrows([angle(perm.zScrPosCpxMean(stind,1,4)),abs(perm.zScrPosCpxMean(stind,1,4))],[1,2]);
% $$$ zscrs = tpp.zscore(phzOrder,stind,1,4);;
% $$$ imagesc(zscrs(:,sind)');
% $$$ caxis([-5.2,5.2]);
% $$$ colormap('jet');
% $$$ 
% $$$ % High Dominant -> Rear
% $$$ figure();
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==2 & ssi(:,2)==1 )';
% $$$ [~,sind] = sortrows([angle(perm.zScrPosCpxMean(stind,1,3)),abs(perm.zScrPosCpxMean(stind,1,3))],[1,2]);
% $$$ zscrs = tpp.zscore(phzOrder,stind,1,3);;
% $$$ imagesc(zscrs(:,sind)');
% $$$ caxis([-5.2,5.2]);
% $$$ colormap('jet');
% $$$ 
% $$$ 
% $$$ 
% $$$ % Rear Dominant -> High
% $$$ figure();
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)         ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==1 & ssi(:,2)==2  ...
% $$$              & sum(tpp.zscore(phzOrder,:,1,1)>2.1)'>4)';
% $$$ [~,sind] = sort(angle(perm.zScrPosCpxMean(stind,1,1)));
% $$$ zscrs = tpp.zscore(phzOrder,stind,1,1);
% $$$ imagesc(repmat(zscrs(:,sind)',[1,2]));
% $$$ caxis([-6.2,6.2]);
% $$$ colormap('jet');
% $$$ 
% $$$ % Low Dominant -> High
% $$$ figure();
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==3 & ssi(:,2)==2  ...
% $$$              & sum(tpp.zscore(phzOrder,:,1,6)>2.1)'>4)';
% $$$ [~,sind] = sort(angle(perm.zScrPosCpxMean(stind,1,6)));
% $$$ zscrs = tpp.zscore(phzOrder,stind,1,6);
% $$$ imagesc(repmat(zscrs(:,sind)',[1,2]));
% $$$ caxis([-6.2,6.2]);
% $$$ colormap('jet');
% $$$ % Low Dominant -> Rear
% $$$ figure();
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==3 & ssi(:,2)==1  ...
% $$$              & sum(tpp.zscore(phzOrder,:,1,5)>2.1)'>4)';
% $$$ [~,sind] = sort(angle(perm.zScrPosCpxMean(stind,1,5)));
% $$$ zscrs = tpp.zscore(phzOrder,stind,1,5);
% $$$ imagesc(repmat(zscrs(:,sind)',[1,2]));
% $$$ caxis([-6.2,6.2]);
% $$$ colormap('jet');
% $$$ 
% $$$ 
% $$$ % High Dominant -> Low
% $$$ figure();
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==2 & ssi(:,2)==3  ...
% $$$              & sum(tpp.zscore(phzOrder,:,1,4)>2.1)'>4)';
% $$$ [~,sind] = sort(angle(perm.zScrPosCpxMean(stind,1,4)));
% $$$ zscrs = tpp.zscore(phzOrder,stind,1,4);
% $$$ imagesc(repmat(zscrs(:,sind)',[1,2]));
% $$$ caxis([-6.2,6.2]);
% $$$ colormap('jet');
% $$$ 
% $$$ % High Dominant -> Rear
% $$$ figure();
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==2 & ssi(:,2)==1  ...
% $$$              & sum(tpp.zscore(phzOrder,:,1,3)>2.1)'>4)';
% $$$ [~,sind] = sort(angle(perm.zScrPosCpxMean(stind,1,3)));
% $$$ zscrs = tpp.zscore(phzOrder,stind,1,3);
% $$$ imagesc(repmat(zscrs(:,sind)',[1,2]));
% $$$ caxis([-6.2,6.2]);
% $$$ colormap('jet');
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ % Sort the z-scores 
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==3 & ssi(:,2)==2 )';
% $$$ figure();
% $$$ imagesc(tpp.zscore(phzOrder,stind,1,6)');
% $$$ caxis([-5.2,5.2]);
% $$$ colormap('jet');
% $$$ 
% $$$ % Sort the z-scores 
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==2 & ssi(:,2)==3 )';
% $$$ tppz = tpp.zscore(phzOrder,stind,1,4)';
% $$$ figure();
% $$$ imagesc(tppz);
% $$$ caxis([-5.2,5.2]);
% $$$ colormap('jet');
% $$$ 
% $$$ % Sort the z-scores 
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==3 & ssi(:,2)==2 )';
% $$$ tppz = sortrows(tpp.zscore(phzOrder,stind,1,6)',5);
% $$$ figure();
% $$$ imagesc(tppz);
% $$$ caxis([-5.2,5.2]);
% $$$ colormap('jet');
% $$$ colorbar();
% $$$ 
% $$$ 
% $$$ 
% $$$ % Sort the z-scores 
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==2 & ssi(:,2)==3 )';
% $$$ 
% $$$ [~,rsi] = sort(mean(tpp.phzRateMean(:,stind,2),1,'omitnan')-mean(tpp.phzRateMean(:,stind,3),1,'omitnan'));
% $$$ figure();
% $$$ imagesc(tpp.zscore(phzOrder,stind(rsi),1,4)');
% $$$ caxis([-3.2,3.2])
% $$$ 
% $$$ figure();
% $$$ subplot(121)
% $$$ imagesc(tpp.phzRateMean(phzOrder,stind,3)');
% $$$ caxis([0,15.2])
% $$$ subplot(122)
% $$$ imagesc(tpp.phzRateMean(phzOrder,stind,2)');
% $$$ caxis([0,15.2])
% $$$ 
% $$$ figure();
% $$$ % z-score rear high
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==1 & ssi(:,2)==2 )';
% $$$ imagesc(tpp.zscore(phzOrder,stind,1,1)');
% $$$ caxis([-5.2,5.2]);
% $$$ % z-score high rear
% $$$ stind = find(  validPosOccRlx             ...
% $$$              & sesInds(unitSubset)        ...
% $$$              & tpp.prmMaxRateSrt(:,2) > 4 ...
% $$$              & validUnits                 ...
% $$$              & ssi(:,1)==1 & ssi(:,2)==2 )';
% $$$ imagesc(tpp.zscore(phzOrder,stind,1,1)');
% $$$ caxis([-5.2,5.2]);
% $$$ 
% $$$ % z-score high low
% $$$ % z-score low  rear
% $$$ 
% $$$ % rear only 
% $$$ % high only 
% $$$ % low  only