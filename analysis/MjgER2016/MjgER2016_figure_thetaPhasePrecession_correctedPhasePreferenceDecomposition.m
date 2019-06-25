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


pind =  all(sq(sum(any(isnan(pftTPZDa(10:30,:,:,1,2:4)),2)))<5,2);

grateUnits = pind&sesInds(unitSubset);
[grate,gpos] = max(permute(sq(mean(RectFilter(permute(pftTPZDa(:,:,grateUnits,:,:),[2,1,3,4,5]),3,1,'circular'),4,'omitnan')),[2,1,3,4,5]));

% $$$ [grate,gpos] = max(permute(sq(mean(permute(pftTPZDa(:,:,pind&sesInds(unitSubset),:,:),[2,1,3,4,5]),4,'omitnan')),[2,1,3,4,5]));

grate = sq(grate(:,:,:,2:4));
gpos = sq(gpos);
[prate,pphz] = max(grate);
clusub = clu(sesInds(unitSubset),:);

npftTPZa = reshape(grate,16,[]);
npftTPZa(~nniz(npftTPZa(:))) = 0;
nind = sum(npftTPZa==0)==16;
npftTPZa(:,nind) =[];

[W,H,D] = nnmf(npftTPZa,3);
% $$$ figure,plot(W(phzOrder,:))
fscrCPPD = zeros([numel(nind),3]);
fscrCPPD(~nind,:) = H';
fscrCPPD = reshape(fscrCPPD,[],3,3);

msr = sq(max(sq(grate(phzOrder,:,:))));
[msr,msi] = sort(msr,2,'descend');
dsr = sq(mean(sq(grate(phzOrder,:,:))));
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


[grate,gpos] = max(permute(sq(mean(RectFilter(permute(pftTPZDa,[2,1,3,4,5]),3,1,'circular'),4,'omitnan')),[2,1,3,4,5]));
grate = sq(grate(:,:,:,2:4));

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
            plot(grate(:,u,s),sclr(s));
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