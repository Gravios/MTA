%MjgER2016_supfig_bfrmCrossValidation

MjgER2016_load_data();
%  Variables:
%      Trials
%      units
%      cluSessionMap
%      pitchReferenceTrial
%      figBasePath
%      sessionListName
%      sessionList
%      states
%      numStates
%      interpParPfsp
%      interpParDfs
%
%  Functions:
%      reshape_eigen_vector

MjgER2016_figure5_args('section 1');
%  Defargs:
%      compute_bhv_ratemaps

% DEF place field rate maps
%pfts                  = cf(@(t,u)   pfs_2d_theta(t,u),                              Trials, units);
pfts   = cf(@(t,u)   compute_ratemaps(t,u),                          Trials, units);
pftCV1 = cf(@(t,u)   compute_ratemaps_crossval(t,u,'tag','cv1'),     Trials, units);
pftCV2 = cf(@(t,u)   compute_ratemaps_crossval(t,u,'tag','cv2'),     Trials, units);

brm    = cf(@(t,u)   compute_bhv_ratemaps(t,u),                      Trials, units);
brmCV1 = cf(@(t,u)   compute_bhv_ratemaps_crossval(t,u,'tag','cv1'), Trials, units);
brmCV2 = cf(@(t,u)   compute_bhv_ratemaps_crossval(t,u,'tag','cv2'), Trials, units);

overwrite = false;
[pfss,metaData]       = req20181106(Trials,[],[],units,'overwrite',overwrite);
[pfssCV1,metaDataCV1] = req20181106_crossval_1st(Trials,[],[],units,'overwrite',overwrite);
[pfssCV2,metaDataCV2] = req20181106_crossval_2nd(Trials,[],[],units,'overwrite',overwrite);

[eigVecs,eigScrs,eigVars,unitSubset,validDims] = compute_bhv_ratemaps_erpPCA(brm,units,'overwrite',overwrite);

pfta    = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'),  pfts  ,units);
clu     = cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),                        pfts  ,units);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])) ,units);
pfta =cat(2,pfta{:});clu=cat(2,clu{:});tlu=cat(2,tlu{:});clu=[tlu',clu'];
[clu,rind] = sortrows(clu);
clu     = clu(unitSubset,:);
pfta    = pfta(:,rind);
pfta    = pfta(:,unitSubset);
pfta(isnan(pfta)) = 0;

brma    = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'),   brm  ,units);
clu     = cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),                         brm  ,units);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])) ,units);
brma    =cat(2,brma{:});clu=cat(2,clu{:});tlu=cat(2,tlu{:});clu=[tlu',clu'];
[clu,rind] = sortrows(clu);
clu     = clu(unitSubset,:);
brma    = brma(:,rind);
brma    = brma(:,unitSubset);
brma(isnan(brma)) = 0;

pftCV1a = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'),  pftCV1,units);
clu     = cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),                        pftCV1,units);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
pftCV1a =cat(2,pftCV1a{:});clu=cat(2,clu{:});tlu=cat(2,tlu{:});clu=[tlu',clu'];
[clu,rind] = sortrows(clu);
clu     = clu(unitSubset,:);
pftCV1a = pftCV1a(:,rind);
pftCV1a = pftCV1a(:,unitSubset);
pftCV1a(isnan(pftCV1a)) = 0;

pftCV2a = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'),  pftCV2,units);
clu     = cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),                        pftCV2,units);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])) ,units);
pftCV2a =cat(2,pftCV2a{:});clu=cat(2,clu{:});tlu=cat(2,tlu{:});clu=[tlu',clu'];
[clu,rind] = sortrows(clu);
clu     = clu(unitSubset,:);
pftCV2a = pftCV2a(:,rind);
pftCV2a = pftCV2a(:,unitSubset);
pftCV2a(isnan(pftCV2a)) = 0;

brmCV1a = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'),  brmCV1,units);
clu     = cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),                        brmCV1,units);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])) ,units);
brmCV1a =cat(2,brmCV1a{:});clu=cat(2,clu{:});tlu=cat(2,tlu{:});clu=[tlu',clu'];
[clu,rind] = sortrows(clu);
clu     = clu(unitSubset,:);
brmCV1a = brmCV1a(:,rind);
brmCV1a = brmCV1a(:,unitSubset);
brmCV1a(isnan(brmCV1a)) = 0;

brmCV2a = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'),  brmCV2,units);
clu     = cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),                        brmCV2,units);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
brmCV2a =cat(2,brmCV2a{:});clu=cat(2,clu{:});tlu=cat(2,tlu{:});clu=[tlu',clu'];
[clu,rind] = sortrows(clu);
clu     = clu(unitSubset,:);
brmCV2a = brmCV2a(:,rind);
brmCV2a = brmCV2a(:,unitSubset);
brmCV2a(isnan(brmCV2a)) = 0;


maskS = create_tensor_mask(pftCV1{1}.adata.bins,struct('shape','circular','radius',440));
maskS = logical(maskS(:));

maskB = false([brm{1}.adata.binSizes']);
maskB(validDims) = true;
maskB = maskB(:);

% COMPUTE crossval of space and behavior
rhoB  = nan([repmat(size(brmCV1a,2),[1,2])]);
rhoS  = nan([repmat(size(pftCV1a,2),[1,2])]);
pvalB  = nan([repmat(size(brmCV1a,2),[1,2])]);
pvalS  = nan([repmat(size(pftCV1a,2),[1,2])]);
nindB = false([size(brmCV1a,1),1]);
nindS = false([size(pftCV1a,1),1]);
for j = 1:size(brmCV1a,2),
    for k = 1:size(brmCV1a,2),
        nindB =   ~isnan(brmCV1a(:,j)) & ~isnan(brmCV2a(:,k)) ...
                  & (brmCV1a(:,j)>0.25 | brmCV2a(:,k)>0.25)     ...
                  & maskB;
        if ~isempty(nindB)&&sum(nindB)>=20,    
            [rhoB(j,k),pvalB(j,k)] = corr(brmCV1a(nindB,j),brmCV2a(nindB,k));
        end
        nindS = ~isnan(pftCV1a(:,j)) & ~isnan(pftCV2a(:,k)) ...
                & (pftCV1a(:,j)>0.25 | pftCV2a(:,k)>0.25)   ... 
                & maskS;
        if ~isempty(nindS)&&sum(nindS)>=20,    
            [rhoS(j,k),pvalB(j,k)] = corr(pftCV1a(nindS,j),pftCV2a(nindS,k));
        end
    end
end

% PLOT spatial vs behavioral correlation between fist and second half of session
figure();
plot(diag(rhoS),diag(rhoB),'.r');
xlim([-1,1]);
ylim([-1,1]);
Lines([],0.5,'r');
Lines(0.5,[],'r');

% COLLECT neuron quality information from each Trial
cf(@(t) t.load('nq'), Trials);
nq = cf(@(t) StructArray(t.nq,1), Trials)
nq = cf(@(n,u) n(u), nq,units);
nq = cat(1,nq{:});

% PLOT cluster quality vs CV corr in {space,behavior}
figure();
subplot(121);
plot(([nq(unitSubset).eDist]),diag(rhoS),'.');
xlabel('eDist');
ylabel('spatial correlation');
title({'crossvalidation of unit stability','within session'});
xlim([10,70]);
subplot(122);
plot(([nq(unitSubset).eDist]),diag(rhoB),'.');
xlabel('eDist');
ylabel('behavioral correlation');
xlim([10,70]);



% PLOT distribution of correlation against random CV unit pairs
figure();
subplot(211)
hold('on');
hax = bar(-1:0.05:1,histc(diag(rhoS),-1:0.05:1)./size(rhoS,1),'histc');
hax.FaceColor = 'r';
hax.EdgeColor = 'r';
hax.FaceAlpha = 0.3
hax.EdgeAlpha = 0.3
hax = bar(-1:0.05:1,histc(nonzeros(triu(rhoS,1)),-1:0.05:1)./(numel(rhoS)/2),'histc');
hax.FaceColor = 'b';
hax.EdgeColor = 'b';
hax.FaceAlpha = 0.3
hax.EdgeAlpha = 0.3
Lines(0.5,[],'r');
subplot(212);
hold('on');
hax = bar(-1:0.05:1,histc(diag(rhoB),-1:0.05:1)./size(rhoB,1),'histc');
hax.FaceColor = 'r';
hax.EdgeColor = 'r';
hax.FaceAlpha = 0.3
hax.EdgeAlpha = 0.3
hax = bar(-1:0.05:1,histc(nonzeros(triu(rhoB,1)),-1:0.05:1)./(numel(rhoB)/2),'histc');
hax.FaceColor = 'b';
hax.EdgeColor = 'b';
hax.FaceAlpha = 0.3
hax.EdgeAlpha = 0.3
Lines(0.5,[],'r');



% SOMETHING ELSE is below


mask = nan(brm{1}.adata.binSizes');
mask(validDims) = 1;
smask = mask;
mask = permute(mask,[1,3,2]);


t = 20;
figPath = create_directory(figBasePath,'figure5','supfig','thetaPhaseDecomposition_CrossValidation');
% $$$ figure();
% $$$ for t = 1:23;
% $$$ 
% $$$ figure


hfig = figure();
for t = 1:23
clf();
Trial = Trials{t};
unitSubset = units{t};
create_directory(fullfile(figPath,Trial.filebase));
%unitSubset = ismember(
sp = reshape(tight_subplot(3,11,0,0.1,0.01),11,3);
for u = unitSubset,    
    if ~ismember([t,u],cluSessionMapSubset,'rows'), continue, end
    
    rmapAll = bsxfun(@times,plot(pfss{t}   ,u,1,11,[],false),mask);

    rmapCV1 = bsxfun(@times,plot(pfssCV1{t},u,1,'',[],false),mask);
    rmapCV2 = bsxfun(@times,plot(pfssCV2{t},u,1,'',[],false),mask);
   
    rmax = prctile([brmCV1{t}.data.rateMap(:,u==brmCV1{t}.data.clu);...
                    brmCV2{t}.data.rateMap(:,u==brmCV2{t}.data.clu)],95);
    %rmax = prctile([rmapCV1(:);rmapCV2(:)],99);
    rmaxCV1 = prctile([rmapCV1(:)],99);
    rmaxCV2 = prctile([rmapCV2(:)],99);    
    if ~nniz(rmaxCV1), rmaxCV1 = 1; end
    if ~nniz(rmaxCV2), rmaxCV2 = 1; end    
    
% THETA place field
    axes(sp(1,1));cla();        
    plot(pfts{t},u,'mean','text',[0,rmax]);
    title(num2str(u));
    axis('xy');  axis('tight');  set(gca,'YTick',[]);  set(gca,'XTick',[]);
    axes(sp(1,2));cla();        
    plot(pftCV1{t},u,'mean','text',[0,rmax]);
    title(num2str(u));
    axis('xy');  axis('tight');  set(gca,'YTick',[]);  set(gca,'XTick',[]);
    axes(sp(1,3));cla();        
    plot(pftCV2{t},u,'mean','text',[0,rmax]);
    title(num2str(u));
    axis('xy');  axis('tight');  set(gca,'YTick',[]);  set(gca,'XTick',[]);
% BRM both halves
    axes(sp(2,1));cla();        
    plot(brm{t},u,'mean','text',[0,rmax],false,'mazeMask',smask);
    axis('xy');  axis('tight');  set(gca,'YTick',[]);  set(gca,'XTick',[]);
% BRM first half
    axes(sp(2,2));cla();        
    plot(brmCV1{t},u,'mean','text',[0,rmax],false,'mazeMask',smask);
    axis('xy');  axis('tight');  set(gca,'YTick',[]);  set(gca,'XTick',[]);
% BRM second half
    axes(sp(2,3));cla();        
    plot(brmCV2{t},u,'mean','text',[0,rmax],false,'mazeMask',smask);
    axis('xy');  axis('tight');  set(gca,'YTick',[]);  set(gca,'XTick',[]);
    
    for p = 1:8,    
% BOTH halves
        axes(sp(p+2,1));  cla();  hold(gca,'on');
        imagescnan({pfd{1}.adata.bins{:},sq(sum(rmapAll(:,phzOrder([p*2-1,p*2]),:),2)./2)'},...
                   [0,rmax],[],false,nanColor,'colorMap',@parula);
        text(pfs.adata.bins{1}(end)-0.45*diff(pfs.adata.bins{1}([1,end])),...
             pfs.adata.bins{3}(end)-0.10*diff(pfs.adata.bins{3}([1,end])),...
             sprintf('%2.0f',rmax),...
             'Color','w','FontWeight','bold','FontSize',8);
        axis('xy');  axis('tight');  set(gca,'YTick',[]);  set(gca,'XTick',[]);                
        
% FIRST half
        axes(sp(p+2,2));  cla();  hold(sp(p+2,2),'on');
        imagescnan({pfd{1}.adata.bins{:},sq(sum(rmapCV1(:,phzOrder([p*2-1,p*2]),:),2)./2)'},...
                   [0,rmaxCV1],[],false,nanColor,'colorMap',@parula);
        text(pfs.adata.bins{1}(end)-0.45*diff(pfs.adata.bins{1}([1,end])),...
             pfs.adata.bins{3}(end)-0.10*diff(pfs.adata.bins{3}([1,end])),...
             sprintf('%2.0f',rmaxCV1),...
             'Color','w','FontWeight','bold','FontSize',8);
        axis('xy');  axis('tight');  set(gca,'YTick',[]);  set(gca,'XTick',[]);        
    
% SECOND half
        axes(sp(p+2,3));  cla();  hold(sp(p+2,3),'on');
        imagescnan({pfd{1}.adata.bins{:},sq(sum(rmapCV2(:,phzOrder([p*2-1,p*2]),:),2)./2)'},...
                   [0,rmaxCV2],[],false,nanColor,'colorMap',@parula);
        text(pfs.adata.bins{1}(end)-0.45*diff(pfs.adata.bins{1}([1,end])),...
             pfs.adata.bins{3}(end)-0.10*diff(pfs.adata.bins{3}([1,end])),...
             sprintf('%2.0f',rmaxCV2),...
             'Color','w','FontWeight','bold','FontSize',8);
        axis('xy');  axis('tight');  set(gca,'YTick',[]);  set(gca,'XTick',[]);
    end
    
    
    
    axes(sp(p+3,1));  cla();  hold(gca,'on');    
    rmapComplex = sum(reshape(bsxfun(@times,                                                     ...
                                     reshape(permute(rmapAll,[1,3,2]),[],size(rmapAll,2)),       ...
                                     exp(-i.*pfssCV2{t}.adata.bins{2}(phzOrder))'),              ...
                              [size(rmapAll,1),size(rmapAll,3),size(rmapAll,2)]),                ...
                      3)                                                                         ...
        ./sq(sum(rmapAll,2));

    imagescnan({pfd{1}.adata.bins{:},angle(rmapComplex)'.*sq(mask)'},...
               [-pi,pi],[],false,nanColor,'colorMap',@hsv);
    axis('xy');  axis('tight');  set(gca,'YTick',[]);  set(gca,'XTick',[]);

    axes(sp(p+3,2));  cla();  hold(gca,'on');    
    rmapComplex = sum(reshape(bsxfun(@times,                                                     ...
                                     reshape(permute(rmapCV1,[1,3,2]),[],size(rmapCV1,2)),       ...
                                     exp(-i.*pfssCV2{t}.adata.bins{2}(phzOrder))'),              ...
                              [size(rmapCV1,1),size(rmapCV1,3),size(rmapCV1,2)]),                ...
                      3)                                                                         ...
        ./sq(sum(rmapCV1,2));

    imagescnan({pfd{1}.adata.bins{:},angle(rmapComplex)'.*sq(mask)'},...
               [-pi,pi],[],false,nanColor,'colorMap',@hsv);
    axis('xy');  axis('tight');  set(gca,'YTick',[]);  set(gca,'XTick',[]);

    axes(sp(p+3,3));  cla();  hold(gca,'on');    
    rmapComplex = sum(reshape(bsxfun(@times,                                                     ...
                                     reshape(permute(rmapCV2,[1,3,2]),[],size(rmapCV2,2)),       ...
                                     exp(-i.*pfssCV2{t}.adata.bins{2}(phzOrder))'),              ...
                              [size(rmapCV2,1),size(rmapCV2,3),size(rmapCV2,2)]),                ...
                      3)                                                                         ...
        ./sq(sum(rmapCV2,2));
    imagescnan({pfd{1}.adata.bins{:},angle(rmapComplex)'.*sq(mask)'},...
               [-pi,pi],[],false,nanColor,'colorMap',@hsv);
    axis('xy');  axis('tight');  set(gca,'YTick',[]);  set(gca,'XTick',[]);        
    
    drawnow();
    pause(0.1);
    
    FigName = ['bfrm_',Trial.filebase,'_unit-',num2str(u)];
    print(hfig,'-dpng',  fullfile(figPath,Trial.filebase,[FigName,'.png']));
    %print(hfig,'-depsc2',fullfile(figPath,Trial.filebase,[FigName,'.eps']));    
end

end