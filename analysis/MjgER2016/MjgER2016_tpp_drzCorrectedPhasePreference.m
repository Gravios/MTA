% MjgER2016_figure_thetaPhasePrecession_correctedPhasePreferenceDecomposition.m
%
%  for each drzphz rate map and each phase, find the rate weighted mean drz.
%
%
% DRZ := Directed Rate Zone (gaussian)
% TPP := Theta Phase Preference

% LOAD drzphz fields
statesIndsGPE = [2,3,4];
statesTPP = stateLabels(statesIndsGPE);
tpp.phzBins = pftHZTPD{1}{1}.adata.bins{2};
tpp.drzBins = pftHZTPD{1}{1}.adata.bins{1};
%%%<<< LOAD rate maps -> ce{states}{Trials}
tpp.ce = cf(@(s)                                                             ...
            cf(@(T,u) MTAApfs(T,u,'tag',['ddtp-','s',num2str(sigma),'-',s]), ...
               Trials,                                                       ...
               units),                                                       ...
            statesTPP);
%%%>>>

%%%<<< DECAPSULATE drzphz fields -> ceMat[drz,phz,unit,iter,state]
tpp.ceMat = cf(@(pfs)                                                         ...
               cf(@(p,u)                                                      ...
                  p.data.rateMap(:,ismember(p.data.clu,u),:),                 ...
                  pfs,units), ...
               tpp.ce);
% COLLECT cluster Ids
tpp.clu     = cf(@(pfs) cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),pfs,units), tpp.ce);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
tpp.ceMat = cf(@(pfs)   cat(2,pfs{:}), tpp.ceMat);
tpp.clu     = cf(@(c)   cat(2,c{:}),   tpp.clu);
assert(all(cell2mat(cf(@(c) all(c), cf(@(c,o) c==o, tpp.clu,circshift(tpp.clu,1,2))))),'Clus do not match');
tpp.ceMat = cat(4,tpp.ceMat{:});
tpp.clu = tpp.clu{1};
tlu=cat(2,tlu{:});
tpp.clu=[tlu',tpp.clu'];
% SORT drzphz field maps by the unitSubset
[tpp.clu,rind] = sortrows(tpp.clu);
tpp.clu     = tpp.clu(unitSubset,:);
tpp.ceMat = tpp.ceMat(:,rind,:,:);
tpp.ceMat = tpp.ceMat(:,unitSubset,:,:);
tpp.ceMat = reshape(permute(tpp.ceMat,[1,5,2,3,4]),...
                    [tpp.ce{1}{1}.adata.binSizes',...
                     numel(unitSubset),...
                     tpp.ce{1}{1}.parameters.numIter,...
                     numel(statesTPP)]);
%%%>>>

%%%<<< COLLECT drz corrected rates -> phzRateMean[phz,unit,state]
% FOR each unit, for each state, for each phase, find the maximum rate and position along the drz (gdz)
% GET phzRateMax
%     phzRateMaxPos
%     phzRateMaxPosMean
%     phzRateMaxPosStd
%     phzRateMean
%     phzRateStd
%     phzRateMeanPos

% $$$ tpp.ceMat = permute(RectFilter(permute(tpp.ceMat,...
% $$$                                         [2,1,3,4,5]),...
% $$$                                 3,1,'circular'),...
% $$$                      [2,1,3,4,5]);

for u = 1:size(tpp.ceMat,3),
    for n = 1:size(tpp.ceMat,4),
        for s = 1:size(tpp.ceMat,5),
            for p = 1:size(tpp.ceMat,2),
                nind = [];
                gind = [];
                nind = find(isnan(tpp.ceMat(:,p,u,n,s)));
                gind = nind( nind<30 & nind>10 );
                if numel(gind)>4, continue; end
                bind = ~ismember(tpp.drzBins,tpp.drzBins(nind));
                if ~isempty(gind)
                    tpp.ceMat(gind,p,u,n,s) = interp1(tpp.drzBins(bind),          ...
                                                      tpp.ceMat(bind,p,u,n,s),    ...
                                                      tpp.drzBins(gind));
                end
            end
        end
    end
end



[tpp.phzRateMax,tpp.phzRateMaxPosInd] = ...
    max(sq(mean(tpp.ceMat,4,'omitnan')));
tpp.phzRateMax = sq(tpp.phzRateMax);
tpp.phzRateMaxPosInd = sq(tpp.phzRateMaxPosInd);

tpp.ceMatSize = size(tpp.ceMat);
tpp.phzRateMeanPos = sq(sum(repmat(tpp.drzBins,[1,tpp.ceMatSize(2:5)]) ... 
                                 .*bsxfun(@rdivide,tpp.ceMat,sum(tpp.ceMat,'omitnan')),...
                             'omitnan'));

tpp.phzRateMeanPosInd = discretize(sq(tpp.phzRateMeanPos(:,:,1,:)),                      ...  Position
                               [tpp.drzBins-diff(tpp.drzBins(1:2));            ...
                                tpp.drzBins(end)+diff(tpp.drzBins(1:2))]);%Position Bins
tpp.phzRateMeanPosInd = reshape(tpp.phzRateMeanPosInd,[],1);
tpp.phzRateMeanAll = nan(tpp.ceMatSize([2,3,5,4]));
for f = 1:tpp.ceMatSize(4),
    ratesm = sq(tpp.ceMat(:,:,:,f,:));
    tpp.phzRateMeanAll(:,:,:,f) = ...
        reshape(ratesm([0:numel(tpp.phzRateMeanPosInd)-1]'.*size(ratesm,1)+tpp.phzRateMeanPosInd),...
                tpp.ceMatSize([2,3,5]));
end
tpp.phzRateMean = mean(tpp.phzRateMeanAll,4,'omitnan');
tpp.phzRateStd  = std(tpp.phzRateMeanAll,[],4,'omitnan');
tpp.phzRateMeanPos = sq(tpp.phzRateMeanPos(:,:,1,:));
%%%>>>

%%%<<< GET BS std phzRateMax (skip)
% $$$ [phzRateMaxStd]   = permute(sq(std(RectFilter(permute(tpp.ceMat,[2,1,3,4,5]), ...
% $$$                                  3,1,'circular'),[],4,'omitnan')),[2,1,3,4,5]);
% $$$ mprSize = size(phzRateMaxStd);
% $$$ phzRateMaxStd = reshape(phzRateMaxStd,mprSize(1),[]);
% $$$ mprPos = reshape(mprPos,[],1);
% $$$ phzRateMaxStda = [];
% $$$ for s = 1:size(phzRateMaxStd,2),
% $$$     phzRateMaxStda(s) = phzRateMaxStd(gpi(s),s);
% $$$ end
% $$$ phzRateMaxStd = phzRateMaxStda;
% $$$ phzRateMaxStd = reshape(sq(phzRateMaxStd),gsaSize(2:end));
% $$$ 
% $$$ [phzRateMaxBS,phzRateMaxPos] = max(permute(sq(RectFilter(permute(tpp.ceMat,[2,1,3,4,5]),3,1,'circular')),[2,1,3,4,5]));
% $$$ phzRateMaxPos = sq(phzRateMaxPos(1,:,:,:,statesIndsGPE));
% $$$ phzRateMaxPos = permute(phzRateMaxPos,[1,2,4,3]);
% $$$ phzRateMaxBS = sq(phzRateMaxBS(1,:,:,:,statesIndsGPE));
% $$$ phzRateMaxBS = permute(phzRateMaxBS,[1,2,4,3]);
% $$$ phzRateMaxPosMean = mean(pftHZTPD{1}{1}.adata.bins{1}(phzRateMaxPos),ndims(phzRateMaxPos),'omitnan');
% $$$ phzRateMaxPosStd = std(phzRateMaxPos,[],ndims(phzRateMaxPos),'omitnan');
% $$$ 
% $$$ phzRateMax = sq(phzRateMax(:,:,:,statesIndsGPE));
% $$$ phzRateMaxStd = sq(phzRateMaxStd(:,:,phzRateMaxPos));
% $$$ %     phzRateMaxPosMean
% $$$ %     phzRateMaxPosStd
% $$$ %     phzRateMean
% $$$ %     phzRateStd
% $$$ %     phzRateMeanPos
%%%>>>

%%%<<< COMPUTE complex value vector of phzRateMax
tpp.phzRateMeanCpx = bsxfun(@times,                                                     ...
                tpp.phzRateMean,                                                     ...
                exp(i.*(tpp.phzBins+double(tpp.phzBins<0).*2*pi)));
tpp.phzRateMeanCpxMean = sq(sum(tpp.phzRateMeanCpx)./sum(tpp.phzRateMean,1,'omitnan'));
%%%>>>

%%%<<< GET indices' order of dominant vs subdominant states 
msr = sq(max(sq(tpp.phzRateMean)));
[msr,msi] = sort(msr,2,'descend');
dsr = sq(mean(sq(tpp.phzRateMean)));
[dsr,dsi] = sort(dsr,2,'descend');
csr = sq(abs(sum(tpp.phzRateMeanCpx)));
[csr,csi] = sort(csr,2,'descend');

ssi = dsi;

validPosOccRlx =  all(sq(sum(any(isnan(tpp.ceMat(10:30,:,:,1,2:3)),2)))<5,2);
validPosOccRr =  all(sq(sum(any(isnan(tpp.ceMat(10:30,:,:,1,1)),2)))<5,2);
validUnits = true(size(validPosOccRr));

% drop index:dind:  bad occ in high or low bad
%            rind:  bad rate est in rear
%            bind:  bad
%            mind:  messy field
%            eind:  bad edge est
%            lind:  low rate
dind = [85,52,91,82,85,89,110,152,205,246,248,308,367];
rind = [44,46,49,59,62,106,112,157,158,186,187,188,191,200,201,202,205,217,220,230,243,253,256,...
        258,259,262,272,275,286,287,317,330,339,348,357,364,371,374,377,398,402];
bind = [25,34,40,57,64,65,69,77,87,117,161,199,204,208,209,213,268,285,308,350,361,376,377,381,...
        383,388];
mind = [18,37,48,121,153,156,196,210,301,305,324,361,375];
eind = [191,196,202,372,406];
lind = [27,72,76,305];
               
validPosOccRr(rind) = false;
validUnits(bind) = false;
validUnits(eind) = false;
validUnits(lind) = false;
validUnits(dind) = false;

% CORRECT bhv pref oder where rear cannot be compared 
tssi = ssi(validPosOccRlx & ~validPosOccRr,[1,2]);
tssi(tssi(:,1)==1&tssi(:,2)==3,1:2) = [3,2];
tssi(tssi(:,1)==1&tssi(:,2)==2,1:2) = [2,3];
tssi(tssi(:,2)==1&tssi(:,1)==3,1:2) = [3,2];
tssi(tssi(:,2)==1&tssi(:,1)==2,1:2) = [2,3];
ssi(validPosOccRlx & ~validPosOccRr,[1,2]) = tssi;
ssi(validPosOccRlx & ~validPosOccRr,3) = 1;

%%%>>>

%%%<<< GET sorted equivalents
tpp.prmPhzAng = angle(tpp.phzRateMeanCpxMean); 
[~,tpp.prmPhzInd] = min(abs(bsxfun(@circ_dist,reshape(tpp.prmPhzAng,1,[]),tpp.phzBins)));
tpp.prmPhzInd = reshape(tpp.prmPhzInd,[],3);
tpp.prmPhzRln = abs(tpp.phzRateMeanCpxMean);
tpp.prmMeanRate = sq(mean(tpp.phzRateMean,'omitnan'));
[tpp.prmMaxRate,tpp.prmMaxRateInd] = max(tpp.phzRateMean,[],'omitnan');
[tpp.prmMinRate,tpp.prmMinRateInd] = min(tpp.phzRateMean,[],'omitnan');
[tpp.prmMinRate,tpp.prmMinRateInd,tpp.prmMaxRate,tpp.prmMaxRateInd] = ...
    deal(sq(tpp.prmMinRate),sq(tpp.prmMinRateInd),sq(tpp.prmMaxRate),sq(tpp.prmMaxRateInd));

tpp.prmPhzAngSrt = nan(size(tpp.prmMeanRate));
tpp.prmPhzRlnSrt = nan(size(tpp.prmMeanRate));
tpp.prmMeanRateSrt   = nan(size(tpp.prmMeanRate));
tpp.prmMaxRateSrt    = nan(size(tpp.prmMeanRate));
tpp.prmMinRateSrt    = nan(size(tpp.prmMeanRate));
tpp.prmMaxRateIndSrt    = nan(size(tpp.prmMeanRate));
tpp.prmMinRateIndSrt    = nan(size(tpp.prmMeanRate));
tpp.prmPosSrt        = nan(size(tpp.phzRateMeanPos));
tpp.prmRateSrt       = nan(size(tpp.phzRateMean));
for j = 1:size(tpp.prmMeanRate,1);
    for s = 1:3,
        tpp.prmPhzAngSrt    (j,s) = tpp.prmPhzAng(j,ssi(j,s));
        tpp.prmPhzRlnSrt    (j,s) = tpp.prmPhzRln(j,ssi(j,s));        
        tpp.prmMeanRateSrt  (j,s) = tpp.prmMeanRate  (j,ssi(j,s));
        tpp.prmMaxRateSrt   (j,s) = tpp.prmMaxRate   (j,ssi(j,s));    
        tpp.prmMinRateSrt   (j,s) = tpp.prmMinRate   (j,ssi(j,s));            
        tpp.prmMaxRateIndSrt(j,s) = tpp.prmMaxRateInd(j,ssi(j,s));
        tpp.prmMinRateIndSrt(j,s) = tpp.prmMinRateInd(j,ssi(j,s));
        tpp.prmPhzIndSrt    (j,s) = tpp.prmPhzInd    (j,ssi(j,s));
        tpp.prmPosSrt     (:,j,s) = tpp.phzRateMeanPos(:,j,ssi(j,s));
        tpp.prmRateSrt    (:,j,s) = tpp.phzRateMean (:,j,ssi(j,s));
    end
end

tpp.prmPhzAngSrtShift = tpp.prmPhzAngSrt + double(tpp.prmPhzAngSrt<0).*2.*pi;
%%%>>>

%%%<<< COLLECT rateMaxDcpp slope analysis 
%     THRESHOLD the rate to find the start, middle, and end of phase precession
%     REPORT the phase and position of each.
%
%     Compute start,middle,end of p,d ∀ u {u ∈ Units, p ∈ Phz, d ∈ Drz}
%     
stind = find(validPosOccRlx & sesInds(unitSubset) & tpp.prmMaxRateSrt(:,2) > 4 & validUnits)';
tpp.smePhz = nan([size(tpp.prmMaxRateSrt,1),2,3]);
tpp.smePhzInd = nan([size(tpp.prmMaxRateSrt,1),2,3]);
tpp.smeDrz = nan([size(tpp.prmMaxRateSrt,1),2,3]);
errCnt = 0;
for u = stind
    for s = 1:2,
        if tpp.prmMaxRateSrt(u,s) > 4,
            rateThr = (tpp.prmMeanRateSrt(u,s)-tpp.prmMinRateSrt(u,s)).*0.5;
            pktrInds = [tpp.prmMinRateIndSrt(u,s),tpp.prmMaxRateIndSrt(u,s)];

            if pktrInds(1) > pktrInds(2),          
                indR = pktrInds(2)+1:pktrInds(1)-1;
                indL = [pktrInds(1)+1:size(tpp.prmRateSrt,1),1:pktrInds(2)-1];
            else
                indL = pktrInds(1)+1:pktrInds(2)-1;
                indR = [pktrInds(2)+1:size(tpp.prmRateSrt,1),1:pktrInds(1)-1];
            end
               

            try
            [~,pBR] = min(abs(tpp.prmRateSrt(indR,u,s)-tpp.prmMinRateSrt(u,s)-rateThr));
            [~,pBL] = min(abs(tpp.prmRateSrt(indL,u,s)-tpp.prmMinRateSrt(u,s)-rateThr));
            tpp.smePhzInd(u,s,:) = [indL(pBL),tpp.prmMaxRateIndSrt(u,s),indR(pBR)];
            tpp.smeDrz(u,s,:) = tpp.prmPosSrt(sq(tpp.smePhzInd(u,s,:)),u,s);
            tpp.smePhz(u,s,:) = tpp.phzBins(tpp.smePhzInd(u,s,:));
            catch
                errCnt = errCnt + 1;
            end

        end
    end
end
%%%>>>

%%%<<< DIAGNOSTIC figures
%% EXAMINE units -----------------------------------------------------------------------------------
% $$$ stind = find(validPosOccRlx & sesInds(unitSubset) & prmMaxRateSrt(:,2) > 4 & validUnits)';
% $$$ figure();
% $$$ for ind = stind(1:end),
% $$$     clf();
% $$$     t = cluSessionMapSubset(ind,1);
% $$$     u = cluSessionMapSubset(ind,2);
% $$$     mrateHZTPD = max(cell2mat(cf(@(p) p{t}.maxRate(u,false), pftHZTPD(statesIndsGPE(:)))));    
% $$$     subplot(1,5,1);
% $$$         plot(pfts{t},u,'mean','text',[],false,[],[],[],@jet,[],nanColor);
% $$$     for s = 1:3,
% $$$         subplot(1,5,s+1);
% $$$         plot(pftHZTPD{statesIndsGPE(s)}{t},u,'mean','',[0,mrateHZTPD],false,[],[],[],@jet,[],nanColor); 
% $$$     end
% $$$     subplot(155);
% $$$         hold('on');        
% $$$         for s = 1:3,
% $$$             plot(sq(tpp.phzRateMean(phzOrder,ind,s)),sclr(s));
% $$$         end    
% $$$         title(num2str([cluSessionMapSubset(ind,:),find(stind==ind), ind]))
% $$$     waitforbuttonpress();
% $$$ end
% $$$ 
% $$$ 

            
figure();
for s = 1:3;
    subplot(1,3,s);
    hold('on');
    for j = nonzeros(double(s~=[1:3]).*[1:3])',
        ind =  ssi(:,1)==s & ssi(:,2)==j        ...
               & validPosOccRlx                 ...
               & sesInds(unitSubset)            ...
               & validUnits;
        plot(tpp.prmPhzAngSrtShift(ind,1),...
             tpp.prmPhzAngSrtShift(ind,2),'.');
    end
    line([0,2*pi],[0,2*pi])
end

% $$$ 
% $$$ sesIds = [3:5,8:12,17:23];
% $$$ %sesIds = [8:12,17:23];
% $$$ sesInds = ismember(cluSessionMap(:,1),sesIds);
% $$$ 
% $$$ figure();
% $$$ for s = [3,2,1]
% $$$     ind =  ssi(:,1) == s                    ...
% $$$            & validPosOccRlx                 ...
% $$$            & sesInds(unitSubset)            ...
% $$$            & tpp.prmMaxRateSrt(:,2) > 4         ...
% $$$            & validUnits;
% $$$ 
% $$$ % $$$     hax=bar(linspace(0,2*pi,25),...
% $$$ % $$$             histc([prmPhzAngSrtShift(ind,1);prmPhzAngSrtShift(ind,1)+2*pi],...
% $$$ % $$$                   linspace(0,4*pi,25)),...
% $$$ % $$$             'histc');
% $$$ % $$$ % $$$     hax=polarhistogram(prmPhzAngSrtShift(ind,1),                 ... theta
% $$$ % $$$ % $$$                        linspace(0,2*pi,13),                          ... bins
% $$$ % $$$ % $$$                        'Normalization','probability'                 ... normalization
% $$$ % $$$ % $$$                       );
% $$$ % $$$     hold('on');    
% $$$ % $$$     hax.FaceColor = sclr(s);
% $$$ % $$$     hax.EdgeColor = sclr(s);    
% $$$ % $$$     hax.FaceAlpha = 0.2;
% $$$ % $$$     hax.EdgeAlpha = 0.4;
% $$$     subplot(3,1,s);
% $$$     hax=bar(linspace(0,4*pi,25),...
% $$$             histc([tpp.prmPhzAngSrtShift(ind,1);tpp.prmPhzAngSrtShift(ind,1)+2*pi],...
% $$$                   linspace(0,4*pi,25)),...
% $$$             'histc');
% $$$     xlim([0,4*pi]);
% $$$     xlabel('theta phase')
% $$$ 
% $$$         
% $$$ % $$$     polarplot([0,circ_mean(prmPhzAngSrtShift(ind,1))],               ... theta
% $$$ % $$$               [0,circ_r(prmPhzAngSrtShift(ind,1))],                  ... res len
% $$$ % $$$               sclr(s),                                               ... color
% $$$ % $$$               'LineWidth',1);
% $$$ end
% $$$ 
% $$$ 
% $$$ figure();hold('on');
% $$$ for s = [1,2,3]
% $$$     
% $$$     ind =    ssi(:,1) == s                  ...
% $$$            & validPosOccRlx                 ...
% $$$            & sesInds(unitSubset)            ...
% $$$            & tpp.prmMaxRateSrt(:,2) > 4         ...
% $$$            & validUnits;
% $$$     
% $$$     %subplot(3,1,s);
% $$$     hax=bar(linspace(-pi,pi,17),...
% $$$             histc(-circ_dist(prmPhzAngSrt(ind,1)+double(prmPhzAngSrt(ind,1)<0).*2.*pi,...
% $$$                             prmPhzAngSrt(ind,2)+double(prmPhzAngSrt(ind,2)<0).*2.*pi),...
% $$$                   linspace(-pi,pi,17)),...
% $$$             'histc');
% $$$     hax.FaceColor = sclr(3);
% $$$     hax.EdgeColor = sclr(3);    
% $$$     hax.FaceAlpha = 0.2;
% $$$     hax.EdgeAlpha = 0.4;
% $$$ end
% $$$ xlim([-pi,pi])
% $$$ 
% $$$ figure();
% $$$ ind = validPosOccRlx                 ...
% $$$      & sesInds(unitSubset)            ...
% $$$      & tpp.prmMaxRateSrt(:,2) > 4         ...
% $$$      & validUnits;
% $$$ h = boxplot(-circ_dist(tpp.prmPhzAngSrtShift(ind,1),...
% $$$                        tpp.prmPhzAngSrtShift(ind,2)),...
% $$$             ssi(ind,1),...
% $$$             'colors','rgb');
% $$$ 
% $$$ 
% $$$ 
% $$$ ind = validPosOccRlx                 ...
% $$$      & sesInds(unitSubset)            ...
% $$$      & tpp.prmMaxRateSrt(:,2) > 4         ...
% $$$      & validUnits;
% $$$ 
% $$$ paH = tpp.prmPhzAngSrtShift(ind&ssi(:,1)==2,1);
% $$$ paL = tpp.prmPhzAngSrtShift(ind&ssi(:,1)==3,1);
% $$$ paH(isnan(paH)) = [];
% $$$ paL(isnan(paL)) = [];
% $$$ 
% $$$ [pConc,f] = circ_ktest(paH,paL)
% $$$ 
% $$$ [pMedian,m,P] = circ_cmtest([paH;paL],...
% $$$                       [ones([numel(paH),1]);2*ones([numel(paL),1])])
% $$$ 
% $$$ [pMeanp,tab] = circ_wwtest([paH;paL],...
% $$$                       [ones([numel(paH),1]);2*ones([numel(paL),1])])
% $$$ 
% $$$ 
% $$$ 
% $$$ figure();
% $$$     plot(abs(diff(smeDrz(ind,1,[1,2]),1,3)),abs(diff(smeDrz(ind,2,[1,2]),1,3)),'.');
% $$$ 
% $$$ figure();
% $$$     hold('on');
% $$$     plot(smeDrz(ind,1,1),smeDrz(ind,2,1),'.');
% $$$     plot(smeDrz(ind,1,3),smeDrz(ind,2,3),'.');
% $$$ 
% $$$     
% $$$ ind = validPosOccRlx                 ...
% $$$      & sesInds(unitSubset)            ...
% $$$      & prmMaxRateSrt(:,2) > 4         ...
% $$$      & validUnits;
% $$$ 
% $$$ figure();
% $$$     hold('on');
% $$$     plot(smeDrz(ind,1,1),smePhz(ind,1,1)+randn(sum(ind),1)./25,'.');
% $$$     plot(smeDrz(ind,2,1),smePhz(ind,2,1)+randn(sum(ind),1)./25,'.');    
% $$$ 
% $$$     
% $$$ ind = validPosOccRlx                 ...
% $$$      & sesInds(unitSubset)            ...
% $$$      & prmMaxRateSrt(:,2) > 4         ...
% $$$      & validUnits                     ...
% $$$      & (ssi(:,1) == 2 | ssi(:,1) == 3);
% $$$ 
% $$$ figure();
% $$$     plot(prmPhzRlnSrt(ind,2)-prmPhzRlnSrt(ind,1),...
% $$$          circ_dist(smePhz(ind,2,2),smePhz(ind,1,2))+randn(sum(ind),1)./25,...
% $$$          '.');
% $$$ 
% $$$ 
% $$$ ind =     validPosOccRlx                   ...
% $$$           & sesInds(unitSubset)            ...
% $$$           & prmMaxRateSrt(:,2) > 4         ...
% $$$           & validUnits                     ...
% $$$           & ( ssi(:,1) == 2 | ssi(:,1) == 3);
% $$$ 
% $$$ figure();% DRZ vs TPP
% $$$     hold('on');
% $$$     plot(repmat(smeDrz(ind,1,2),[2,1]),...
% $$$          [smePhz(ind,1,2)+randn(sum(ind),1)./25;smePhz(ind,1,2)+randn(sum(ind),1)./25+2*pi],'.');
% $$$     plot(repmat(smeDrz(ind,2,2),[2,1]),...
% $$$          [smePhz(ind,2,2)+randn(sum(ind),1)./25;smePhz(ind,2,2)+randn(sum(ind),1)./25+2*pi],'.');
% $$$ 
% $$$ figure();% DRZ vs TPP
% $$$     for s = 1:3
% $$$         ind =     validPosOccRlx                 ...
% $$$                 & sesInds(unitSubset)            ...
% $$$                 & prmMaxRateSrt(:,2) > 4         ...
% $$$                 & validUnits                     ...
% $$$                 & ( ssi(:,1) == s);
% $$$         for p = 1:3,
% $$$             subplot2(3,3,s,p);
% $$$             plot(smeDrz(ind,2,p)-smeDrz(ind,1,p),...
% $$$                  circ_dist(smePhz(ind,2,p),smePhz(ind,1,p))+randn(sum(ind),1)./25,...
% $$$                  '.');
% $$$             xlim([-0.5,0.5]);
% $$$             ylim([-pi,pi]);
% $$$         end
% $$$     end
% $$$ 
% $$$ figure();% DRZ vs TPP
% $$$     for s = 1:3
% $$$         ind =     validPosOccRlx                 ...
% $$$                 & sesInds(unitSubset)            ...
% $$$                 & prmMaxRateSrt(:,2) > 4         ...
% $$$                 & validUnits                     ...
% $$$                 & ( ssi(:,1) == s);
% $$$             subplot2(3,1,s,1);
% $$$             plot(circ_dist(smePhz(ind,2,1),smePhz(ind,1,1))+randn(sum(ind),1)./25,...
% $$$                  circ_dist(smePhz(ind,2,3),smePhz(ind,1,3))+randn(sum(ind),1)./25,...
% $$$                  '.');
% $$$             xlim([-pi,pi]);
% $$$             ylim([-pi,pi]);
% $$$     end
% $$$ 
% $$$ 
% $$$     
% $$$ figure(); % rate vs TPP
% $$$     for s = 1:3,
% $$$         ind =     validPosOccRlx                 ...
% $$$                 & sesInds(unitSubset)            ...
% $$$                 & prmMaxRateSrt(:,2) > 4         ...
% $$$                 & validUnits                     ...
% $$$                 & ( ssi(:,1) == s);
% $$$         subplot2(1,3,1,s);
% $$$         plot([prmMaxRateSrt(ind,1);prmMaxRateSrt(ind,1)],...
% $$$              [prmPhzAngSrt(ind,1);prmPhzAngSrt(ind,1)+2*pi], ...
% $$$              '.')
% $$$         Lines([],pi,'k');
% $$$         title(statesTPP(s))
% $$$     end
% $$$ 
% $$$ dDrz = repmat(smeDrz(ind,1,2),[2,1]);
% $$$ dPhz = [smePhz(ind,1,2)+randn(sum(ind),1)./25;smePhz(ind,1,2)+randn(sum(ind),1)./25+2*pi];
% $$$ dsDrzVec = repmat(smeDrz(ind,2,2),[2,1])-repmat(smeDrz(ind,1,2),[2,1]);
% $$$ dsPhzVec = -circ_dist([smePhz(ind,2,2)+randn(sum(ind),1)./25;smePhz(ind,2,2)+randn(sum(ind),1)./25+2*pi],...
% $$$                      [smePhz(ind,1,2)+randn(sum(ind),1)./25;smePhz(ind,1,2)+randn(sum(ind),1)./25+2*pi]);
% $$$ dsVecMag = sqrt(sum([dsDrzVec,dsPhzVec].^2,2));
% $$$ 
% $$$ figure();
% $$$ smePhz(
% $$$ 
% $$$ figure();
% $$$     hold('on');
% $$$     h = quiver(dDrz,...
% $$$                dPhz,...
% $$$                dsDrzVec./dsVecMag,...
% $$$                dsPhzVec./dsVecMag,...
% $$$                0);
% $$$     
% $$$ 
% $$$     
% $$$ figure();
% $$$     subplot(131);
% $$$     plot(smePhz(ind,1,1)+randn(sum(ind),1)./10,smePhz(ind,2,1)+randn(sum(ind),1)./10,'.');    
% $$$     subplot(132);
% $$$     plot(smePhz(ind,1,2)+randn(sum(ind),1)./10,smePhz(ind,2,2)+randn(sum(ind),1)./10,'.');    
% $$$     subplot(133);
% $$$     plot(smePhz(ind,1,3)+randn(sum(ind),1)./10,smePhz(ind,2,3)+randn(sum(ind),1)./10,'.');    
% $$$     
% $$$     figure();
% $$$     subplot(131);
% $$$     plot(smePhz(ind,1,1)+randn(sum(ind),1)./10,smePhz(ind,2,1)+randn(sum(ind),1)./10,'.');    
% $$$ 
% $$$     
% $$$     figure();
% $$$         plot(abs(circ_dist(smePhz(ind,1,2),smePhz(ind,1,1)))+abs(circ_dist(smePhz(ind,1,2),smePhz(ind,1,3)))+randn(sum(ind),1)./10,...
% $$$              abs(circ_dist(smePhz(ind,2,2),smePhz(ind,2,1)))+abs(circ_dist(smePhz(ind,2,2),smePhz(ind,2,3)))+randn(sum(ind),1)./10,'.')
% $$$         
% $$$     figure();
% $$$     plot(circ_dist(smePhz(ind,1,3),smePhz(ind,1,2)).*(-double(smePhz(ind,1,3)<=0)+double(smePhz(ind,1,3)>0))...
% $$$              +circ_dist(smePhz(ind,1,2),smePhz(ind,1,1))...
% $$$              +randn(sum(ind),1)./10,...
% $$$          circ_dist(smePhz(ind,2,3),smePhz(ind,2,2)).*(-double(smePhz(ind,2,3)<=0)+double(smePhz(ind,2,3)>0))...
% $$$          +circ_dist(smePhz(ind,2,2),smePhz(ind,2,1))...
% $$$          +randn(sum(ind),1)./10,'.');
% $$$         
%%%>>>    
