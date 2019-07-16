% MjgER2016_figure_thetaPhasePrecession_correctedPhasePreferenceDecomposition.m


sigma = 150;
pftHZTPD = cf(@(s) cf(@(T,u) MTAApfs(T,u,'tag',['ddtp-','s',num2str(sigma),'-',s]), Trials, units), stateLabels);
pftTPZDa = cf(@(pfs) ...
              cf(@(p,u) ...
                 p.data.rateMap(:,ismember(p.data.clu,u),:),  ...
                 pfs,units), ...
              pftHZTPD);

% $$$ u =3;figure;plot(sq(max(reshape(pftTPZDa{3}{20}(:,u,:),[40,16,101]))),'c');hold on, ...
% $$$    plot(sq(max(reshape(pftTPZDa{4}{20}(:,u,:),[40,16,101]))),'b');
% $$$ u =3;figure;plot(mean(RectFilter(sq(max(reshape(pftTPZDa{3}{20}(:,u,:),[40,16,101]))),3,1,'circular'),2),'c');
% $$$ hold('on'); plot(mean(RectFilter(sq(max(reshape(pftTPZDa{4}{20}(:,u,:),[40,16,101]))),3,1,'circular'),2),'b');

clu     = cf(@(pfs) cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),                        pfs,units), pftHZTPD);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
pftTPZDa = cf(@(pfs) cat(2,pfs{:}), pftTPZDa);
clu     = cf(@(c)   cat(2,c{:}),clu);
assert(all(cell2mat(cf(@(c) all(c), cf(@(c,o) c==o, clu,circshift(clu,1,2))))),'Clus do not match');
pftTPZDa = cat(4,pftTPZDa{:});
clu = clu{1};
tlu=cat(2,tlu{:});
clu=[tlu',clu'];
[clu,rind] = sortrows(clu);
clu     = clu(unitSubset,:);
pftTPZDa = pftTPZDa(:,rind,:,:);
pftTPZDa = pftTPZDa(:,unitSubset,:,:);
pftTPZDa = reshape(permute(pftTPZDa,[1,5,2,3,4]),[pftHZTPD{1}{1}.adata.binSizes',...
                    numel(unitSubset),pftHZTPD{1}{1}.parameters.numIter,numel(stateLabels)]);

% For each unit, for each state, for each phase, find the maximum rate and position along the gdz

% GET BS mean maxPhzRate
[maxPhzRate,mprPos] = max(permute(sq(mean(RectFilter(permute(pftTPZDa,[2,1,3,4,5]),3,1,'circular'),4,'omitnan')),[2,1,3,4,5]));

% GET BS std maxPhzRate
[maxPhzRateStd]   = permute(sq(std(RectFilter(permute(pftTPZDa,[2,1,3,4,5]), ...
                                 3,1,'circular'),[],4,'omitnan')),[2,1,3,4,5]);
mprSize = size(maxPhzRateStd);
maxPhzRateStd = reshape(maxPhzRateStd,mprSize(1),[]);
mprPos = reshape(mprPos,[],1);
maxPhzRateStda = [];
for s = 1:size(maxPhzRateStd,2),
    maxPhzRateStda(s) = maxPhzRateStd(gpi(s),s);
end
maxPhzRateStd = maxPhzRateStda;
maxPhzRateStd = reshape(sq(maxPhzRateStd),gsaSize(2:end));
maxPhzRateStd = sq(maxPhzRateStd(:,:,2:4));

[maxPhzRateBS,mprPosBS] = max(permute(sq(RectFilter(permute(pftTPZDa,[2,1,3,4,5]),3,1,'circular')),[2,1,3,4,5]));
mprPosBS = sq(mprPosBS(1,:,:,:,2:4));
mprPosBS = permute(mprPosBS,[1,2,4,3]);
maxPhzRateBS = sq(maxPhzRateBS(1,:,:,:,2:4));
maxPhzRateBS = permute(maxPhzRateBS,[1,2,4,3]);

mprPosBSMean = mean(pftHZTPD{1}{1}.adata.bins{1}(mprPosBS),ndims(mprPosBS),'omitnan');
mprPosBSStd = std(mprPosBS,[],ndims(mprPosBS),'omitnan');


% $$$ figure,plot(mprPosBSMean(:,314,3),pftHZTPD{1}{1}.adata.bins{2})
% $$$ figure,hold('on')
% $$$ u = 309;
% $$$ for s=1:3;
% $$$ plot(sq(maxPhzRateBS(phzOrder,u,s,:)),sclr(s))
% $$$ plot(mean(maxPhzRateBS(phzOrder,u,s),4,'omitnan'),'c','LineWidth',2)
% $$$ plot(maxPhzRate(phzOrder,u,s),'k','LineWidth',2)
% $$$ plot(maxPhzRate(phzOrder,u,s)+maxPhzRateStd(phzOrder,u,s),'m','LineWidth',2)
% $$$ plot(maxPhzRate(phzOrder,u,s)-maxPhzRateStd(phzOrder,u,s),'m','LineWidth',2)
% $$$ end

% SQUEEZE matrix 
% FILTER States
% GET max rate for each behavior
maxPhzRate = sq(maxPhzRate(:,:,:,2:4));
mprPos = sq(mprPos);
mprPos = mprPos(:,:,2:4);
mprPos = pftHZTPD{1}{1}.adata.bins{1}(mprPos);
[prate,pphz] = max(maxPhzRate);

% SELECT unit subset for phase preference decomposition
validPosOcc =  all(sq(sum(any(isnan(pftTPZDa(10:30,:,:,1,2:4)),2)))<5,2);
mprDecompUnits = validPosOcc & sesInds(unitSubset);

npftTPZa = reshape(maxPhzRate(:,mprDecompUnits,:),16,[]);
npftTPZa(~nniz(npftTPZa(:))) = 0;
nind = sum(npftTPZa==0)==16;
npftTPZa(:,nind) =[];

[W,H,D] = nnmf(npftTPZa,3);
% $$$ figure,plot(W(phzOrder,:))
fscrCPPD = zeros([numel(nind),3]);
fscrCPPD(~nind,:) = H';
fscrCPPD = reshape(fscrCPPD,[],3,3);

msr = sq(max(sq(maxPhzRate(phzOrder,:,:))));
[msr,msi] = sort(msr,2,'descend');
dsr = sq(mean(sq(maxPhzRate(phzOrder,:,:))));
[dsr,dsi] = sort(dsr,2,'descend');

nTPZmean = mean(npftTPZa,2);
nTPstd = std(npftTPZa,[],2);
DPhz = cov(npftTPZa');

pftTPZDs = sq(pftTPZDa(:,:,:,:,2:4));
W = bsxfun(@rdivide,W,sum(W));
FSCFrPhz = W * inv(W' * W);
rk = size(FSCFrPhz,2);%rank(D,1e-4);       % why not on D? would save on corrcoef(X) computation
FSCFrPhz = FSCFrPhz .* repmat(sqrt(diag(DPhz)),1,rk); % compute rotated factor scores from the normalized raw
FSrCPhz = multiprod(pftTPZDs,FSCFrPhz,[2],[1,2]);

fscrCPPD = reshape(fscrCPPD,[],3,3);


% PROJECT rateCPPD (maxPhzRate) onto nnmf vectors 
fscrCPPDAll = permute(multiprod(maxPhzRate,FSCFrPhz,[1],[1,2]),[2,3,1])./10000;

figure,plot(fscrCPPD(:,1,1),fscrCPPDAll(mprDecompUnits,1,1)/10000,'.');
line([0,0.2],[0,0.2])





sclr = 'rgb';
figure(1)
clf();
sp = tight_subplot(4,3,0,0.2);
u = 97;
while u~=-1,
    figure(1);
    mrate = prctile(nonzeros(pftTPZDs(:,:,u,:)),99);    
    for s = 1:3,
        axes(sp(s));cla();
            imagesc(mean(pftTPZDs(:,:,u,:,s),4,'omitnan')');
            caxis([0,mrate]);
            axis('xy')
            title(num2str(cluSessionMapSubsetCPPD(u,:)));
        axes(sp(s+3));
            imagesc(mean(pftTPZDs(:,:,u,:,s),4,'omitnan')');
            caxis([0,mrate]);        
            axis('xy')
        axes(sp(s+6));cla();hold('on')
            for g = 1:101
                for p = 1:3
                    plot(FSrC(:,p,u,g,s),sclr(p))
                end
            end

        axes(sp(s+9));        
            plot(maxPhzRate(:,u,s),sclr(s));
    end
    mscr = max(nonzeros(FSrC(:,:,u,:,:)));
    af(@(x) ylim(x,[0,mscr]),sp(s+4:s+6));
    axes('Position',[0.2,0.05,0.1,0.15]);
    plot(W);
    axes('Position',[0.4,0.05,0.1,0.15]);    
    plot(pfts{cluSessionMapSubsetCPPD(u,1)},cluSessionMapSubsetCPPD(u,2),'mean','text',[],true);
    axes('Position',[0.6,0.05,0.1,0.15]);        
    plot(pfbs{cluSessionMapSubsetCPPD(u,1)},cluSessionMapSubsetCPPD(u,2),'mean','text',[],false);
    
    figure(2);
    for s = 1:3,
        for p = 1:3,
            subplot2(3,3,s,p);cla();hold('on')
            for g = 1:3

            plot((mean(FSrC(:,g,u,:,s),4,'omitnan')-mean(FSrC(:,g,u,:,p),4,'omitnan'))...
                 ./sqrt(0.5*(var(FSrC(:,g,u,:,s),[],4,'omitnan')...
                             +var(FSrC(:,g,u,:,p),[],4,'omitnan'))),sclr(g));
            end
        end
    end
    
    u = figure_controls(figure(2),u,1:size(pftTPZDs,3));
end


FSrC(FSrC<0|isnan(FSrC)) = 0;
FSrCN = bsxfun(@rdivide,FSrC,sum(FSrC));
comCPPD = sq(sum(bsxfun(@times,pftHZTPD{1}{1}.adata.bins{1},FSrCN)));


figure,
for s = 1:3,
    subplot(2,3,s),
        hold('on');
        plot(comCPPD(:,:,s));
        plot(mean(comCPPD(:,:,s),2,'omitnan'),'k--','LineWidth',2); 
        ylim([-1,1]);
    
    subplot(2,3,s+3),
        hold('on');
        plot(diff(comCPPD(:,:,s)));
        plot(mean(diff(comCPPD(:,:,s)),2,'omitnan'),'k--','LineWidth',2); 
        ylim([-1,1]);
end


figure,
for s = 1:3,
    subplot(1,3,s);
        plot(comCPPD(2,:,s)-comCPPD(1,:,s),comCPPD(3,:,s)-comCPPD(2,:,s),'.');
        ylim([-1,1]);
        xlim([-1,1]);
end


figure,
for s = 1:3,
    for j = 1:3,
        subplot2(3,3,s,j);
            hold('on');
            plot(comCPPD(:,:,s)-comCPPD(:,:,j))
    end
end





% COMPUTE complex value vector of maxPhzRate
mprCpx = bsxfun(@times,                                                     ...
                maxPhzRate,                                                 ...
                exp(i.*(pftHZTPD{1}{1}.adata.bins{2}+double(pftHZTPD{1}{1}.adata.bins{2}<0).*2*pi)));
%exp(-i.*pftHZTPD{1}{1}.adata.bins{2}));
mprMeanPhz = sq(sum(mprCpx)./sum(maxPhzRate,1,'omitnan'));


msr = sq(max(sq(maxPhzRate)));
[msr,msi] = sort(msr,2,'descend');
dsr = sq(mean(sq(maxPhzRate)));
[dsr,dsi] = sort(dsr,2,'descend');
csr = sq(abs(sum(mprCpx)));
[csr,csi] = sort(csr,2,'descend');

ssi = dsi;


validPosOccRlx =  all(sq(sum(any(isnan(pftTPZDa(10:30,:,:,1,3:4)),2)))<5,2);
validPosOccRr =  all(sq(sum(any(isnan(pftTPZDa(10:30,:,:,1,2)),2)))<5,2);

% CORRECT bhv pref oder where rear cannot be compared 
tssi = ssi(validPosOccRlx & ~validPosOccRr,[1,2]);
tssi(tssi(:,1)==1&tssi(:,2)==3,1) = 2;
tssi(tssi(:,1)==1&tssi(:,2)==2,1) = 3;
tssi(tssi(:,2)==1&tssi(:,1)==3,2) = 2;
tssi(tssi(:,2)==1&tssi(:,1)==2,2) = 3;
ssi(validPosOccRlx & ~validPosOccRr,[1,2]) = tssi;
ssi(validPosOccRlx & ~validPosOccRr,3) = 1;

stsColor = 'rgb';
mprMeanPhzAng = angle(mprMeanPhz);
mprMeanPhzRln = abs(mprMeanPhz);
mprMeanRate = sq(mean(maxPhzRate,'omitnan'));
mprMaxRate = sq(max(maxPhzRate,[],'omitnan'));

mprMeanPhzAngSrt = nan(size(mprMeanRate));
mprMeanPhzRlnSrt = nan(size(mprMeanRate));
mprMeanRateSrt   = nan(size(mprMeanRate));
mprMaxRateSrt    = nan(size(mprMeanRate));
mprPosSrt        = nan(size(mprPos));
for j = 1:size(mprMeanRate,1);
    for s = 1:3,
        mprMeanPhzAngSrt(j,s) = mprMeanPhzAng(j,ssi(j,s));
        mprMeanPhzRlnSrt(j,s) = mprMeanPhzRln(j,ssi(j,s));        
        mprMeanRateSrt  (j,s) = mprMeanRate  (j,ssi(j,s));
        mprMaxRateSrt   (j,s) = mprMaxRate   (j,ssi(j,s));    
        mprPosSrt     (:,j,s) = mprPos     (:,j,ssi(j,s));
    end
end


mprMeanPhzAngSrtShift = mprMeanPhzAngSrt(:,1) + double(mprMeanPhzAngSrt(:,1)<0).*2.*pi;


