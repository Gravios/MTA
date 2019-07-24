% MjgER2016_figure_thetaPhasePrecession_correctedPhasePreferenceDecomposition.m
% DRZ := Directed Rate Zone (gaussian)
% TPP := Theta Phase Preference

% LOAD drzphz fields
statesIndsGPE = [2,3,4];
statesTPP = stateLabels(statesIndsGPE);
tpp.phzBins = pftHZTPD{1}{1}.adata.bins{2};

%%%<<< LOAD rate maps
tpp.ce = cf(@(s)                                                             ...
            cf(@(T,u) MTAApfs(T,u,'tag',['ddtp-','s',num2str(sigma),'-',s]), ...
               Trials,                                                       ...
               units),                                                       ...
            stateLabels);
%%>>>

%%%<<< DECAPSULATE drzphz fields
tpp.ceMat = cf(@(pfs)                                                         ...
               cf(@(p,u)                                                      ...
                  p.data.rateMap(:,ismember(p.data.clu,u),:),                 ...
                  pfs,units), ...
               tpp.ce);
% COLLECT cluster Ids
clu     = cf(@(pfs) cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),pfs,units), pftHZTPD);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
pftTPZDa = cf(@(pfs) cat(2,pfs{:}), pftTPZDa);
clu     = cf(@(c)   cat(2,c{:}),clu);
assert(all(cell2mat(cf(@(c) all(c), cf(@(c,o) c==o, clu,circshift(clu,1,2))))),'Clus do not match');
pftTPZDa = cat(4,pftTPZDa{:});
clu = clu{1};
tlu=cat(2,tlu{:});
clu=[tlu',clu'];
% SORT drzphz field maps by the unitSubset
[clu,rind] = sortrows(clu);
clu     = clu(unitSubset,:);
pftTPZDa = pftTPZDa(:,rind,:,:);
pftTPZDa = pftTPZDa(:,unitSubset,:,:);
pftTPZDa = reshape(permute(pftTPZDa,[1,5,2,3,4]),[pftHZTPD{1}{1}.adata.binSizes',...
                    numel(unitSubset),pftHZTPD{1}{1}.parameters.numIter,numel(stateLabels)]);
%%%>>>

% FOR each unit, for each state, for each phase, find the maximum rate and position along the gdz

% GET phzRateMax
%     phzRateMaxPos
%     phzRateMaxPosMean
%     phzRateMaxPosStd
%     phzRateMean
%     phzRateStd
%     phzRateMeanPos


[phzRateMax,mprPos] = max(permute(sq(mean(RectFilter(permute(pftTPZDa,[2,1,3,4,5]),3,1,'circular'),4,'omitnan')),[2,1,3,4,5]));

rates = permute(RectFilter(permute(pftTPZDa,...
                                   [2,1,3,4,5]),...
                           3,1,'circular'),...
                [2,1,3,4,5]);
rateSize = size(rates);
dweights = bsxfun(@rdivide,rates,sum(rates,'omitnan'));
phzRateMeanPos = sq(sum(repmat(dbins,[1,rateSize(2:5)]).*dweights,'omitnan'));
phzRateMeanPosInd = discretize(sq(phzRateMeanPos(:,:,1,:)),                      ...  Position
                               [dbins-diff(dbins(1:2));dbins(end)+diff(dbins(1:2))]);%Position Bins
phzRateMeanPosInd = reshape(phzRateMeanPosInd,[],1);
ratesa = nan(rateSize([2,3,5,4]));
for f = 1:size(rates,4),
    ratesm = sq(rates(:,:,:,f,:));
    ratesa(:,:,:,f) = reshape(ratesm([0:numel(phzRateMeanPosInd)-1]'.*size(ratesm,1)+phzRateMeanPosInd),...
                              rateSize([2,3,5]));
end

phzRateMean = mean(ratesa,4,'omitnan');
phzRateMean = phzRateMean(:,:,statesIndsGPE);
phzRateStd  = std(ratesa,[],4,'omitnan');
phzRateStd  = phzRateStd(:,:,statesIndsGPE);
phzRateMeanPos = sq(phzRateMeanPos(:,:,1,statesIndsGPE));



% GET BS std phzRateMax
[phzRateMaxStd]   = permute(sq(std(RectFilter(permute(pftTPZDa,[2,1,3,4,5]), ...
                                 3,1,'circular'),[],4,'omitnan')),[2,1,3,4,5]);
mprSize = size(phzRateMaxStd);
phzRateMaxStd = reshape(phzRateMaxStd,mprSize(1),[]);
mprPos = reshape(mprPos,[],1);
phzRateMaxStda = [];
for s = 1:size(phzRateMaxStd,2),
    phzRateMaxStda(s) = phzRateMaxStd(gpi(s),s);
end
phzRateMaxStd = phzRateMaxStda;
phzRateMaxStd = reshape(sq(phzRateMaxStd),gsaSize(2:end));

[phzRateMaxBS,phzRateMaxPos] = max(permute(sq(RectFilter(permute(pftTPZDa,[2,1,3,4,5]),3,1,'circular')),[2,1,3,4,5]));
phzRateMaxPos = sq(phzRateMaxPos(1,:,:,:,statesIndsGPE));
phzRateMaxPos = permute(phzRateMaxPos,[1,2,4,3]);
phzRateMaxBS = sq(phzRateMaxBS(1,:,:,:,statesIndsGPE));
phzRateMaxBS = permute(phzRateMaxBS,[1,2,4,3]);
phzRateMaxPosMean = mean(pftHZTPD{1}{1}.adata.bins{1}(phzRateMaxPos),ndims(phzRateMaxPos),'omitnan');
phzRateMaxPosStd = std(phzRateMaxPos,[],ndims(phzRateMaxPos),'omitnan');



phzRateMax = sq(phzRateMax(:,:,:,statesIndsGPE));
phzRateMaxStd = sq(phzRateMaxStd(:,:,phzRateMaxPos));
%     phzRateMaxPosMean
%     phzRateMaxPosStd
%     phzRateMean
%     phzRateStd
%     phzRateMeanPos



% COMPUTE complex value vector of phzRateMax
prmCpx = bsxfun(@times,                                                     ...
                phzRateMean,                                                     ...
                exp(i.*(pftHZTPD{1}{1}.adata.bins{2}+double(pftHZTPD{1}{1}.adata.bins{2}<0).*2*pi)));
prmPhz = sq(sum(prmCpx)./sum(phzRateMean,1,'omitnan'));


msr = sq(max(sq(phzRateMean)));
[msr,msi] = sort(msr,2,'descend');
dsr = sq(mean(sq(phzRateMean)));
[dsr,dsi] = sort(dsr,2,'descend');
csr = sq(abs(sum(prmCpx)));
[csr,csi] = sort(csr,2,'descend');

ssi = dsi;


validPosOccRlx =  all(sq(sum(any(isnan(pftTPZDa(10:30,:,:,1,3:4)),2)))<5,2);
validPosOccRr =  all(sq(sum(any(isnan(pftTPZDa(10:30,:,:,1,2)),2)))<5,2);
validUnits = true(size(validPosOccRr));

% drop index: bad occ in high or low bad
%             bad rate est in rear
%             bad
%             messy field
%             bad edge est
%             low rate
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
tssi(tssi(:,1)==1&tssi(:,2)==3,1) = 2;
tssi(tssi(:,1)==1&tssi(:,2)==2,1) = 3;
tssi(tssi(:,2)==1&tssi(:,1)==3,2) = 2;
tssi(tssi(:,2)==1&tssi(:,1)==2,2) = 3;
ssi(validPosOccRlx & ~validPosOccRr,[1,2]) = tssi;
ssi(validPosOccRlx & ~validPosOccRr,3) = 1;

stsColor = 'rgb';
prmPhzAng = angle(prmPhz);
[~,prmPhzInd] = min(abs(bsxfun(@circ_dist,reshape(prmPhzAng,1,[]),pftHZTPD{1}{1}.adata.bins{2})));
prmPhzInd = reshape(prmPhzInd,[],3);
prmPhzRln = abs(prmPhz);
prmMeanRate = sq(mean(phzRateMean,'omitnan'));
[prmMaxRate,prmMaxRateInd] = max(phzRateMean,[],'omitnan');
[prmMinRate,prmMinRateInd] = min(phzRateMean,[],'omitnan');
[prmMinRate,prmMinRateInd,prmMaxRate,prmMaxRateInd] = ...
    deal(sq(prmMinRate),sq(prmMinRateInd),sq(prmMaxRate),sq(prmMaxRateInd));

prmPhzAngSrt = nan(size(prmMeanRate));
prmPhzRlnSrt = nan(size(prmMeanRate));
prmMeanRateSrt   = nan(size(prmMeanRate));
prmMaxRateSrt    = nan(size(prmMeanRate));
prmMinRateSrt    = nan(size(prmMeanRate));
prmMaxRateIndSrt    = nan(size(prmMeanRate));
prmMinRateIndSrt    = nan(size(prmMeanRate));
prmPosSrt        = nan(size(phzRateMeanPos));
prmRateSrt       = nan(size(phzRateMean));
for j = 1:size(prmMeanRate,1);
    for s = 1:3,
        prmPhzAngSrt(j,s) = prmPhzAng(j,ssi(j,s));
        prmPhzRlnSrt(j,s) = prmPhzRln(j,ssi(j,s));        
        prmMeanRateSrt  (j,s) = prmMeanRate  (j,ssi(j,s));
        prmMaxRateSrt   (j,s) = prmMaxRate   (j,ssi(j,s));    
        prmMinRateSrt   (j,s) = prmMinRate   (j,ssi(j,s));            
        prmMaxRateIndSrt(j,s) = prmMaxRateInd(j,ssi(j,s));
        prmMinRateIndSrt(j,s) = prmMinRateInd(j,ssi(j,s));
        prmPhzIndSrt    (j,s) = prmPhzInd    (j,ssi(j,s));
        prmPosSrt     (:,j,s) = phzRateMeanPos(:,j,ssi(j,s));
        prmRateSrt    (:,j,s) = phzRateMean (:,j,ssi(j,s));
    end
end

prmPhzAngSrtShift = prmPhzAngSrt + double(prmPhzAngSrt<0).*2.*pi;







%% EXAMINE units -----------------------------------------------------------------------------------
stind = find(validPosOccRlx & sesInds(unitSubset) & prmMaxRateSrt(:,2) > 4 & validUnits)';
figure();
for ind = stind(1:end),
    clf();
    t = cluSessionMapSubset(ind,1);
    u = cluSessionMapSubset(ind,2);
    mrateHZTPD = max(cell2mat(cf(@(p) p{t}.maxRate(u,false), pftHZTPD(statesIndsGPE(:)))));    
    subplot(1,5,1);
        plot(pfts{t},u,'mean','text',[],false,[],[],[],@jet,[],nanColor);
    for s = 1:3,
        subplot(1,5,s+1);
        plot(pftHZTPD{statesIndsGPE(s)}{t},u,'mean','',[0,mrateHZTPD],false,[],[],[],@jet,[],nanColor); 
    end
    subplot(155);
        hold('on');        
        for s = 1:3,
            plot(sq(phzRateMean(phzOrder,ind,s)),sclr(s));
        end    
        title(num2str([cluSessionMapSubset(ind,:),find(stind==ind), ind]))
    waitforbuttonpress();
end



%% rateMaxDcpp slope analysis ----------------------------------------------------------------------
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




            
            
figure();
for s = 1:3;
    subplot(1,3,s);
    hold('on');
    for j = nonzeros(double(s~=[1:3]).*[1:3])',
        ind =  ssi(:,1)==s & ssi(:,2)==j        ...
               & validPosOccRlx                 ...
               & sesInds(unitSubset)            ...
               & prmMaxRateSrt(:,2) > 4         ...
               & validUnits;
        plot(prmPhzAngSrt(ind,1)+double(prmPhzAngSrt(ind,1)<0).*2.*pi,...
             prmPhzAngSrt(ind,2)+double(prmPhzAngSrt(ind,2)<0).*2.*pi,smarkers(j));
    end
    line([0,2*pi],[0,2*pi])
end

sesIds = [3:5,8:12,17:23];
%sesIds = [8:12,17:23];
sesInds = ismember(cluSessionMap(:,1),sesIds);

figure();
for s = [3,2,1]
    ind =  ssi(:,1) == s                    ...
           & validPosOccRlx                 ...
           & sesInds(unitSubset)            ...
           & prmMaxRateSrt(:,2) > 4         ...
           & validUnits;

% $$$     hax=bar(linspace(0,2*pi,25),...
% $$$             histc([prmPhzAngSrtShift(ind,1);prmPhzAngSrtShift(ind,1)+2*pi],...
% $$$                   linspace(0,4*pi,25)),...
% $$$             'histc');
% $$$ % $$$     hax=polarhistogram(prmPhzAngSrtShift(ind,1),                 ... theta
% $$$ % $$$                        linspace(0,2*pi,13),                          ... bins
% $$$ % $$$                        'Normalization','probability'                 ... normalization
% $$$ % $$$                       );
% $$$     hold('on');    
% $$$     hax.FaceColor = sclr(s);
% $$$     hax.EdgeColor = sclr(s);    
% $$$     hax.FaceAlpha = 0.2;
% $$$     hax.EdgeAlpha = 0.4;
    subplot(3,1,s);
    hax=bar(linspace(0,4*pi,25),...
            histc([prmPhzAngSrtShift(ind,1);prmPhzAngSrtShift(ind,1)+2*pi],...
                  linspace(0,4*pi,25)),...
            'histc');
    
% $$$     polarplot([0,circ_mean(prmPhzAngSrtShift(ind,1))],               ... theta
% $$$               [0,circ_r(prmPhzAngSrtShift(ind,1))],                  ... res len
% $$$               sclr(s),                                               ... color
% $$$               'LineWidth',1);
end


figure();hold('on');
for s = [1,2,3]
    
    ind =    ssi(:,1) == s                  ...
           & validPosOccRlx                 ...
           & sesInds(unitSubset)            ...
           & prmMaxRateSrt(:,2) > 4         ...
           & validUnits;
    
    %subplot(3,1,s);
    hax=bar(linspace(-pi,pi,17),...
            histc(-circ_dist(prmPhzAngSrt(ind,1)+double(prmPhzAngSrt(ind,1)<0).*2.*pi,...
                            prmPhzAngSrt(ind,2)+double(prmPhzAngSrt(ind,2)<0).*2.*pi),...
                  linspace(-pi,pi,17)),...
            'histc');
    hax.FaceColor = sclr(3);
    hax.EdgeColor = sclr(3);    
    hax.FaceAlpha = 0.2;
    hax.EdgeAlpha = 0.4;
end
xlim([-pi,pi])

figure();
ind = validPosOccRlx                 ...
     & sesInds(unitSubset)            ...
     & prmMaxRateSrt(:,2) > 4         ...
     & validUnits;
h = boxplot(-circ_dist(prmPhzAngSrt(ind,1)+double(prmPhzAngSrt(ind,1)<0).*2.*pi,...
                       prmPhzAngSrt(ind,2)+double(prmPhzAngSrt(ind,2)<0).*2.*pi),...
            ssi(ind,1),...
            'colors','rgb');



circ_ktest(prmPhzAngSrtShift(ind&ssi(:,1)==2,1),...
           prmPhzAngSrtShift(ind&ssi(:,1)==3,2));

[p,m,P] = circ_cmtest([prmPhzAngSrtShift(ind&ssi(:,1)==2,1);...
                      prmPhzAngSrtShift(ind&ssi(:,1)==3,1)],...
                      [ones([sum(ind&ssi(:,1)==2),1]);ones([sum(ind&ssi(:,1)==3),1])])


paH = prmPhzAngSrtShift(ind&ssi(:,1)==2,1);
paL = prmPhzAngSrtShift(ind&ssi(:,1)==3,1);
paH(isnan(paH)) = [];
paL(isnan(paL)) = [];

[p,tab] = circ_wwtest([paH;paL],...
                      [ones([numel(paH),1]);2*ones([numel(paL),1])])

circ_ktest(paH,paL);


[p,m,P] = kstest2(prmPhzAngSrtShift(ind&ssi(:,1)==2,1),...
                      prmPhzAngSrtShift(ind&ssi(:,1)==3,1))


[h,p] = kstest2(-circ_dist(prmPhzAngSrtShift(ind&ssi(:,1)==2,1),...
                           prmPhzAngSrtShift(ind&ssi(:,1)==2,2)),...
                -circ_dist(prmPhzAngSrtShift(ind&ssi(:,1)==3,1),...
                           prmPhzAngSrtShift(ind&ssi(:,1)==3,2)));

[p,h] = ranksum(-circ_dist(prmPhzAngSrtShift(ind&ssi(:,1)==2,1),...
                           prmPhzAngSrtShift(ind&ssi(:,1)==2,2)),...
                -circ_dist(prmPhzAngSrtShift(ind&ssi(:,1)==3,1),...
                           prmPhzAngSrtShift(ind&ssi(:,1)==3,2)));

% test high->low phase shift distribution is different from low->high
[h,p] = kstest2(-circ_dist(prmPhzAngSrtShift(ind&ssi(:,1)==2&ssi(:,2)==3,1),...
                           prmPhzAngSrtShift(ind&ssi(:,1)==2&ssi(:,2)==3,2)),...
                -circ_dist(prmPhzAngSrtShift(ind&ssi(:,1)==3&ssi(:,2)==2,1),...
                           prmPhzAngSrtShift(ind&ssi(:,1)==3&ssi(:,2)==2,2)));


figure();
    plot(abs(diff(smeDrz(ind,1,[1,2]),1,3)),abs(diff(smeDrz(ind,2,[1,2]),1,3)),'.');

figure();
    hold('on');
    plot(smeDrz(ind,1,1),smeDrz(ind,2,1),'.');
    plot(smeDrz(ind,1,3),smeDrz(ind,2,3),'.');

    
ind = validPosOccRlx                 ...
     & sesInds(unitSubset)            ...
     & prmMaxRateSrt(:,2) > 4         ...
     & validUnits;

figure();
    hold('on');
    plot(smeDrz(ind,1,1),smePhz(ind,1,1)+randn(sum(ind),1)./25,'.');
    plot(smeDrz(ind,2,1),smePhz(ind,2,1)+randn(sum(ind),1)./25,'.');    

    
ind = validPosOccRlx                 ...
     & sesInds(unitSubset)            ...
     & prmMaxRateSrt(:,2) > 4         ...
     & validUnits                     ...
     & (ssi(:,1) == 2 | ssi(:,1) == 3);

figure();
    plot(prmPhzRlnSrt(ind,2)-prmPhzRlnSrt(ind,1),...
         circ_dist(smePhz(ind,2,2),smePhz(ind,1,2))+randn(sum(ind),1)./25,...
         '.');


ind =     validPosOccRlx                 ...
          & sesInds(unitSubset)            ...
          & prmMaxRateSrt(:,2) > 4         ...
          & validUnits                     ...
          & ( ssi(:,1) == 2 | ssi(:,1) == 3);

figure();% DRZ vs TPP
    hold('on');
    plot(repmat(smeDrz(ind,1,2),[2,1]),...
         [smePhz(ind,1,2)+randn(sum(ind),1)./25;smePhz(ind,1,2)+randn(sum(ind),1)./25+2*pi],'.');
    plot(repmat(smeDrz(ind,2,2),[2,1]),...
         [smePhz(ind,2,2)+randn(sum(ind),1)./25;smePhz(ind,2,2)+randn(sum(ind),1)./25+2*pi],'.');

figure();% DRZ vs TPP
    for s = 1:3
        ind =     validPosOccRlx                 ...
                & sesInds(unitSubset)            ...
                & prmMaxRateSrt(:,2) > 4         ...
                & validUnits                     ...
                & ( ssi(:,1) == s);
        for p = 1:3,
            subplot2(3,3,s,p);
            plot(smeDrz(ind,2,p)-smeDrz(ind,1,p),...
                 circ_dist(smePhz(ind,2,p),smePhz(ind,1,p))+randn(sum(ind),1)./25,...
                 '.');
            xlim([-0.5,0.5]);
            ylim([-pi,pi]);
        end
    end

figure();% DRZ vs TPP
    for s = 1:3
        ind =     validPosOccRlx                 ...
                & sesInds(unitSubset)            ...
                & prmMaxRateSrt(:,2) > 4         ...
                & validUnits                     ...
                & ( ssi(:,1) == s);
            subplot2(3,1,s,1);
            plot(circ_dist(smePhz(ind,2,1),smePhz(ind,1,1))+randn(sum(ind),1)./25,...
                 circ_dist(smePhz(ind,2,3),smePhz(ind,1,3))+randn(sum(ind),1)./25,...
                 '.');
            xlim([-pi,pi]);
            ylim([-pi,pi]);
    end


    
figure(); % rate vs TPP
    for s = 1:3,
        ind =     validPosOccRlx                 ...
                & sesInds(unitSubset)            ...
                & prmMaxRateSrt(:,2) > 4         ...
                & validUnits                     ...
                & ( ssi(:,1) == s);
        subplot2(1,3,1,s);
        plot([prmMaxRateSrt(ind,1);prmMaxRateSrt(ind,1)],...
             [prmPhzAngSrt(ind,1);prmPhzAngSrt(ind,1)+2*pi], ...
             '.')
        Lines([],pi,'k');
        title(statesTPP(s))
    end

dDrz = repmat(smeDrz(ind,1,2),[2,1]);
dPhz = [smePhz(ind,1,2)+randn(sum(ind),1)./25;smePhz(ind,1,2)+randn(sum(ind),1)./25+2*pi];
dsDrzVec = repmat(smeDrz(ind,2,2),[2,1])-repmat(smeDrz(ind,1,2),[2,1]);
dsPhzVec = -circ_dist([smePhz(ind,2,2)+randn(sum(ind),1)./25;smePhz(ind,2,2)+randn(sum(ind),1)./25+2*pi],...
                     [smePhz(ind,1,2)+randn(sum(ind),1)./25;smePhz(ind,1,2)+randn(sum(ind),1)./25+2*pi]);
dsVecMag = sqrt(sum([dsDrzVec,dsPhzVec].^2,2));

figure();
smePhz(

figure();
    hold('on');
    h = quiver(dDrz,...
               dPhz,...
               dsDrzVec./dsVecMag,...
               dsPhzVec./dsVecMag,...
               0);
    

    
figure();
    subplot(131);
    plot(smePhz(ind,1,1)+randn(sum(ind),1)./10,smePhz(ind,2,1)+randn(sum(ind),1)./10,'.');    
    subplot(132);
    plot(smePhz(ind,1,2)+randn(sum(ind),1)./10,smePhz(ind,2,2)+randn(sum(ind),1)./10,'.');    
    subplot(133);
    plot(smePhz(ind,1,3)+randn(sum(ind),1)./10,smePhz(ind,2,3)+randn(sum(ind),1)./10,'.');    
    
    figure();
    subplot(131);
    plot(smePhz(ind,1,1)+randn(sum(ind),1)./10,smePhz(ind,2,1)+randn(sum(ind),1)./10,'.');    

    
    figure();
        plot(abs(circ_dist(smePhz(ind,1,2),smePhz(ind,1,1)))+abs(circ_dist(smePhz(ind,1,2),smePhz(ind,1,3)))+randn(sum(ind),1)./10,...
             abs(circ_dist(smePhz(ind,2,2),smePhz(ind,2,1)))+abs(circ_dist(smePhz(ind,2,2),smePhz(ind,2,3)))+randn(sum(ind),1)./10,'.')
        
    figure();
    plot(circ_dist(smePhz(ind,1,3),smePhz(ind,1,2)).*(-double(smePhz(ind,1,3)<=0)+double(smePhz(ind,1,3)>0))...
             +circ_dist(smePhz(ind,1,2),smePhz(ind,1,1))...
             +randn(sum(ind),1)./10,...
         circ_dist(smePhz(ind,2,3),smePhz(ind,2,2)).*(-double(smePhz(ind,2,3)<=0)+double(smePhz(ind,2,3)>0))...
         +circ_dist(smePhz(ind,2,2),smePhz(ind,2,1))...
         +randn(sum(ind),1)./10,'.');
        