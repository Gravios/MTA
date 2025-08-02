function [eigVecs, eigScrs, eigVars, eSpi, FSrBhvThp, metaData] = req20181119_nnmf(varargin);
%  Tags: behavior ratemap theta phase
%  Status: active
%  Type: Analysis
%  Author: Justin Graboski
%  Final_Forms: NA
%  Project: MjgER2016:placefields
%  Description: Theta phase resolved behavior rate maps
% 
%  Dependencies: req20181106 - generation of behavioral tuning maps 
%
%  Protocol: 
%    1. select interneurons
%    2  compute ratemaps
% 
%  Figures: 
%    1. 
% 
%  OUTPUT : eigVals, eigVecs, eigScrs, eigVars, eSpi
%

% GLOBALS ------------------------------------------------------------------------------------------
global MTA_PROJECT_PATH;
%---------------------------------------------------------------------------------------------------


% $$$ %varargin = {};
% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('Trials',                        {{}},                                          ...
                 'sessionListName',               'MjgER2016',                                   ...
                 'tag',                           'hbpptbpFS1v4',                                ...
                 'units',                         [],                                            ...
                 'sampleRate',                    250,                                           ...
                 'stcMode',                       'msnn_ppsvd_raux',                             ...
                 'states',                        {{'theta','rear','hloc','hpause',              ...
                                                   'lloc','lpause','groom','sit'}},              ...
                 'feature_fcn',                   @fet_HB_pitchB,                                ...
                 'overwrite',                     false                                          ...
);
[Trials,sessionListName,tag,units,sampleRate,stcMode,states,feature_fcn,overwrite] =             ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% SET Meta Data ------------------------------------------------------------------------------------
analysisHash = DataHash({sessionListName,tag,sampleRate,stcMode,states,func2str(feature_fcn)});
metaFilePath = fullfile(MTA_PROJECT_PATH,'analysis',[mfilename '-meta-',analysisHash,'.mat']);
dataFilePath = fullfile(MTA_PROJECT_PATH,'analysis',[mfilename '-data-',analysisHash,'.mat']);
metaVars = {'sessionListName','tag','sampleRate','stcMode','states','feature_fcn','analysisHash'};
dataVars = {'eigVecs','eigScrs','eigVars','eSpi','FSrBhvThp'};
%---------------------------------------------------------------------------------------------------

MjgER2016_load_bhv_erpPCA_scores();

if ~exist(dataFilePath,'file')  ||  overwrite,
    if isempty(Trials),
        Trials = af(@(s) MTATrial.validate(s), get_session_list(sessionListName));
    end
    if isempty(units),
        units = cf(@(T)  select_placefields(T),  Trials); 
        units = req20180123_remove_bad_units(units);
    end
    cluSessionMap = [];
    for u = 1:numel(units),
        cluSessionMap = cat(1,cluSessionMap,[u*ones([numel(units{u}),1]),units{u}(:)]);
    end
    

    bfrm = cf(@(t,u)  compute_bhv_ratemaps(t,u), Trials, units);
    [~,~,~,~,validDims] = compute_bhv_ratemaps_erpPCA(bfrm);
    
    % OLD
    %[~,~,~,~,validDims,~,~] = req20180123_pfd_erpPCA([],[],'HBPITCHxBPITCH_v13',[],[],[],false);

    
    ucnt = 1;

    nfac = 3;
    
    eigVecs = [];
    eigScrs = [];
    eigVars = [];
    eSpi = [];
    
    for t = 1:numel(Trials),
        Trial = Trials{t}; 
        unitSubset = units{t};
        %pfs = MTAApfs(Trial,'tag',[']);
        pfs = req20181106(Trial,sessionListName,tag,unitSubset);

        mask = zeros(pfs.adata.binSizes([1,3])');
        mask(validDims) = 1;
        mask = repmat(permute(mask,[1,3,2]),[1,pfs.adata.binSizes(2),1]);
        
        for u = unitSubset,
            rmap = plot(pfs,u,1,'colorbar',[],true,0.5,[],[],[],[],mask);
            rmap(~mask) = nan;
            if sum(rmap(:),'omitnan')~=0,
                rmapNind = sq(sum(reshape(~isnan(rmap(:)),size(rmap)),2))>=pfs.adata.binSizes(2);
                zmap = reshape(permute(rmap,[1,3,2]),[],pfs.adata.binSizes(2));
                zmap(repmat(~rmapNind(:),[1,pfs.adata.binSizes(2)])) = 0;            
                zmap(zmap(:)<0) = 0;            
                rmapNind = rmapNind&reshape(sum(reshape(zmap,[],pfs.adata.binSizes(2))','omitnan')~=0,size(rmap,1),size(rmap,3));
                % COMPUTE erpPCA 
                pmap = zmap;
                pmap(isnan(pmap(:,1))|sum(pmap','omitnan')'==0,:) = [];
                if size(pmap,1)>20    
                    mrate= prctile(rmap(:),99);
                    %[LU,LR,FSr,VT] = erpPCA(pmap',nfac);
                    [W,H,D] = nnmf(pmap,nfac);
                    for v = 1:nfac,
                        evec = nan([size(rmap,1),size(rmap,3)]);
                        evec(rmapNind) = W(:,v);
                        eigVecs(ucnt,v,:,:) = evec;
                        pmapHat = (W(:,v)*H(v,:))';
                        eigVars(ucnt,v) = sum(var(pmapHat))./sum(var(pmap'));
                    end
                    %eigVars(ucnt,:) = VT(1:nfac,4);                
                    eigScrs(ucnt,:,:) = H';
                end

                for p = 1:pfs.adata.binSizes(2),
                    mrate = mean(reshape(zmap(:,p),[],1),'omitnan');
                    pmap = zmap(:,p);
                    pmap = pmap(rmapNind);
                    eSpi(ucnt,1,p) = sum(pmap./mean(pmap,'omitnan').*log2(pmap./mean(pmap,'omitnan')),'omitnan')./numel(pmap);
                end

            end
            ucnt = ucnt+1;
        end%for u
    end%for t

    cEVecs = reshape(eigVecs,size(eigVecs,1),size(eigVecs,2),[]);
    cEVecs = cEVecs(:,:,validDims);
    cEVecs(isnan(cEVecs(:))) = 0;
    
    FSrBhvThp = [];
% $$$     for v = 1:size(cEVecs,2),
% $$$         FSrBhvThp(:,:,v) = sq(cEVecs(:,v,:)) * FSCFr;
% $$$     end

    save(dataFilePath,'-v7.3',dataVars{:});
    save(metaFilePath,'-v7.3',metaVars{:});
else
    load(dataFilePath);
end

metaData = load(metaFilePath);




% $$$ uind = {ismember(cluSessionMap(:,1),[1:2]);...
% $$$         ismember(cluSessionMap(:,1),[3:5]);...
% $$$         ismember(cluSessionMap(:,1),[6:7]);...
% $$$         ismember(cluSessionMap(:,1),[8:12]);...
% $$$         ismember(cluSessionMap(:,1),[13:16]);...
% $$$         ismember(cluSessionMap(:,1),[17:23])};
% $$$ sid = {'er01';'ER06';'Ed10';'jg04';'jg04';'jg05'};
% $$$ 
% $$$ csmUnitSubset = cluSessionMap(unitSubsets{1},:);
% $$$ csmEigScrs = eigScrs(unitSubsets{1},:,:);
% $$$ csmEigVecs = eigVecs(unitSubsets{1},:,:,:);
% $$$ csmEigVars = eigVars(unitSubsets{1},:,:);
% $$$ 
% $$$ 
% $$$ 
% $$$ % dim reduced csmEigVecs
% $$$ cEVecs = reshape(csmEigVecs,size(csmEigVecs,1),size(csmEigVecs,2),[]);
% $$$ cEVecs = cEVecs(:,:,validDims{1});
% $$$ cEVecs(isnan(cEVecs(:))) = 0;
% $$$ 
% $$$ FSrBhvThp = [];
% $$$ for v = 1:size(cEVecs,2),
% $$$     FSrBhvThp(:,:,v) = sq(cEVecs(:,v,:)) * FSCFr;
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ clear('i');
% $$$ % $$$ uind = ismember(csmUnitSubset(:,1),[3:5]);
% $$$ % $$$ uind = ismember(csmUnitSubset(:,1),[6:7]);
% $$$ 
% $$$ [~,csmPhzPrefMax] = max(sq(csmEigScrs(:,:,:)),[],3);
% $$$ 
% $$$ csmPhzPrefMean = angle(sum(reshape(bsxfun(@times,reshape(bsxfun(@plus,csmEigScrs,abs(min(csmEigScrs,[],3))),[],size(csmEigScrs,3)),exp(-i.*pfs.adata.bins{2})'),size(csmEigScrs)),3)./sum(bsxfun(@plus,csmEigScrs,abs(min(csmEigScrs,[],3))),3));
% $$$ 
% $$$ figure,plot(pfs.adata.bins{2}(abs(csmPhzPrefMax(:,1)))+randn(size(csmPhzPrefMax(:,1)))/10,...
% $$$             csmPhzPrefMean(:,1),'.');
% $$$ % $$$ figure,plot(sq(csmEigScrs(:,1,:))')
% $$$ 
% $$$ figure,plot(log10(csmEigVars(:,1)),log10(csmEigVars(:,2)),'.')
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ plot(circ_dist(csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)),csmEigVars(:,1),'.b')
% $$$ plot(circ_dist(csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)),csmEigVars(:,2),'.r')
% $$$ 
% $$$ 
% $$$ figure,
% $$$ plot(circ_mean([csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)]')+pi,csmEigVars(:,1),'.b')
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure();
% $$$ ind = csmEigVars(:,1)>55|csmEigVars(:,2)>30;
% $$$ %ind = csmEigVars(:,1)>60;
% $$$ %ind = ':';
% $$$ subplot(121);
% $$$ plot(csmPhzPrefMean(ind,1),circ_dist(csmPhzPrefMean(ind,1),csmPhzPrefMean(ind,2)),'.b')
% $$$ subplot(122);
% $$$ plot(csmPhzPrefMean(ind,2),circ_dist(csmPhzPrefMean(ind,1),csmPhzPrefMean(ind,2)),'.r')
% $$$ 
% $$$ 
% $$$ figure();
% $$$ scatter3(csmPhzPrefMean(ind,1),csmPhzPrefMean(ind,2),circ_dist(csmPhzPrefMean(ind,1),csmPhzPrefMean(ind,2)),10,log10(csmEigVars(ind,1)./csmEigVars(ind,2)),'filled');
% $$$ colormap(gca,'jet');
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ ind = csmEigVars(:,1)>50;
% $$$ plot(csmPhzPrefMean(ind,1),csmEigVars(ind,1),'.b')
% $$$ ind = csmEigVars(:,2)>20;
% $$$ plot(csmPhzPrefMean(ind,2),csmEigVars(ind,2),'.r')
% $$$ 
% $$$ figure();
% $$$ bar(linspace(-pi,pi,16),...
% $$$     histc([csmPhzPrefMean(csmEigVars(:,1)>55,1)],...
% $$$           linspace(-pi,pi,16)),...
% $$$     'histc');
% $$$ figure();
% $$$ bar(linspace(-pi,pi,16),...
% $$$     histc([csmPhzPrefMean(csmEigVars(:,1)>55&csmEigVars(:,2)>20,1)],...
% $$$           linspace(-pi,pi,16)),...
% $$$     'histc');
% $$$ 
% $$$ figure();
% $$$ bar(linspace(-pi,pi,16),...
% $$$     histc([csmPhzPrefMean(csmEigVars(:,1)>55&csmEigVars(:,2)<40,1)],...
% $$$           linspace(-pi,pi,16)),...
% $$$     'histc');
% $$$ 
% $$$ 
% $$$ figure();
% $$$ ind = csmEigVars(:,1)>55&csmEigVars(:,2)>20;
% $$$ hist(circ_dist(csmPhzPrefMean(ind,1),csmPhzPrefMean(ind,2)))
% $$$ 
% $$$ % from the preferred phases of the two main factors find a value alhpa, by which the phases are shifted,
% $$$ % which minimizes the variance of one of the two groups.
% $$$ 
% $$$ %csmPhzPrefMeanShifted
% $$$ 
% $$$ ind = csmEigVars(:,1)>40&csmEigVars(:,2)>30;
% $$$ ind = csmEigVars(:,1)>65;
% $$$ ind = nniz(csmPhzPrefMean(:,1:2));
% $$$ alpha = linspace(-pi,pi,100);
% $$$ cppmsVar = zeros([numel(alpha),2]);
% $$$ for s = 1:numel(alpha),
% $$$     cppms = [circ_dist(csmPhzPrefMean(ind,1),alpha(s)),circ_dist(csmPhzPrefMean(ind,2),alpha(s))]; 
% $$$     swpInd = bsxfun(@plus,[1,2],bsxfun(@times,double(gt(abs(cppms(:,1)),abs(cppms(:,2)))),[1,-1]));
% $$$     for a = 1:size(cppms,1),
% $$$         cppms(a,:) = cppms(a,swpInd(a,:));
% $$$     end
% $$$     cppmsVar(s,:) = circ_var(cppms);
% $$$ end
% $$$ 
% $$$ figure,plot(cppmsVar)
% $$$ 
% $$$ [~,minCppmsVarInd] = min(cppmsVar(:,1));
% $$$ 
% $$$ cppms = [circ_dist(csmPhzPrefMean(ind,1),alpha(minCppmsVarInd)),...
% $$$          circ_dist(csmPhzPrefMean(ind,2),alpha(minCppmsVarInd))]; 
% $$$ swpInd = bsxfun(@plus,[1,2],bsxfun(@times,double(gt(abs(cppms(:,1)),abs(cppms(:,2)))),[1,-1]));
% $$$ 
% $$$ cppm = csmPhzPrefMean(ind,1:2);
% $$$ for a = 1:size(cppm,1),
% $$$     cppm(a,:) = cppm(a,swpInd(a,:));
% $$$ end
% $$$ 
% $$$ figure();
% $$$ scatter(cppm(:,1),cppm(:,2),10,csmEigVars(ind,1),'filled');
% $$$ xlim([-pi,pi])
% $$$ ylim([-pi,pi])
% $$$ colormap(gca,'jet');
% $$$ 
% $$$ figure();
% $$$ scatter(csmPhzPrefMean(ind,1),csmPhzPrefMean(ind,2),10,csmEigVars(ind,1),'filled');xlim([-pi,pi]);ylim([-pi,pi]);
% $$$ colormap(gca,'jet');
% $$$ 
% $$$ figure();
% $$$ bar(linspace(-pi,pi,16),...
% $$$     histc([cppm(:,1)],...
% $$$           linspace(-pi,pi,16)),...
% $$$     'histc');
% $$$ 
% $$$ 
% $$$ figure,
% $$$ scatter(csmEigVars(ind,1),csmEigVars(ind,2),10,cppm(:,1),'filled');
% $$$ colormap(gca(),'hsv');
% $$$ caxis([-pi,pi]);
% $$$ 
% $$$ 
% $$$ % PLOT 1st vs 2nd EigenVectorVariance colored by 1st EigenVectorProjection onto bhvpfsEigenVectors
% $$$ figure,
% $$$ scatter(csmEigVars(ind,1),csmEigVars(ind,2),10,FSrBhvThp(ind,1,3),'filled');
% $$$ colormap(gca(),'jet');
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ plot(csmPhzPrefMean(:,1),log10(csmEigVars(:,1)),'.b')
% $$$ plot(csmPhzPrefMean(:,2),log10(csmEigVars(:,2)),'.r')
% $$$ 
% $$$ 
% $$$ figure();
% $$$ hist2([[csmPhzPrefMean(:,1);csmPhzPrefMean(:,2)],[csmEigVars(:,1);csmEigVars(:,2)]],10,10);
% $$$ 
% $$$ 
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ plot(csmPhzPrefMean(:,1),csmEigVars(:,1)./csmEigVars(:,2),'.b')
% $$$ %plot(csmPhzPrefMean(:,2),csmEigVars(:,2),'.b')
% $$$ 
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ for v = 1:size(csmPhzPrefMean),
% $$$     plot([csmPhzPrefMean(v,1),csmPhzPrefMean(v,2)],[csmEigVars(v,1),csmEigVars(v,2)])
% $$$ end
% $$$ 
% $$$ 
% $$$ figure();
% $$$ hist2([csmPhzPrefMean(:,1),csmEigVars(:,1);csmPhzPrefMean(:,2),csmEigVars(:,2)],10,10);
% $$$ 
% $$$ figure
% $$$ for u = 1:numel(uind)
% $$$     subplot(numel(uind),3,u*3-2);
% $$$     ind = uind{u}(unitSubsets{1})&(csmEigVars(:,1)>40&csmEigVars(:,2)>20)&FSrC(:,2)<0;
% $$$     %bar(pfs.adata.bins{2},histc(pfs.adata.bins{2}(abs(csmPhzPrefMax(ind,1))),pfs.adata.bins{2}),'histc');
% $$$     bar(pfs.adata.bins{2},histc(csmPhzPrefMean(ind,1),pfs.adata.bins{2}),'histc');    
% $$$     axis('tight');    
% $$$     title(sid{u});
% $$$     subplot(numel(uind),3,u*3-1);
% $$$     %bar(pfs.adata.bins{2},histc(pfs.adata.bins{2}(abs(csmPhzPrefMax(ind,2))),pfs.adata.bins{2}),'histc');    
% $$$     bar(pfs.adata.bins{2},histc(csmPhzPrefMean(ind,2),pfs.adata.bins{2}),'histc');
% $$$     axis('tight');
% $$$     title(sid{u});
% $$$     subplot(numel(uind),3,u*3);    
% $$$     hist2([csmPhzPrefMean(ind,1),...
% $$$            csmPhzPrefMean(ind,2)], ...
% $$$           pfs.adata.bins{2},pfs.adata.bins{2});
% $$$ % $$$     hist2([pfs.adata.bins{2}(abs(csmPhzPrefMax(ind,1))),...
% $$$ % $$$            pfs.adata.bins{2}(abs(csmPhzPrefMax(ind,2)))], ...
% $$$ % $$$           pfs.adata.bins{2},pfs.adata.bins{2});
% $$$     
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,
% $$$ ind = uind&csmEigVars(:,1)>50;
% $$$ ind = uind&all(csmEigVars(:,1:2)>20,2);
% $$$ subplot(121);
% $$$ scatter3(FSrBhvThp(ind,1,1),FSrBhvThp(ind,2,1),FSrBhvThp(ind,3,1),20,pfs.adata.bins{2}(csmPhzPrefMax(ind,1)),'filled');
% $$$ colormap(gca,'hsv');
% $$$ caxis([-pi,pi]);
% $$$ subplot(122);
% $$$ scatter3(FSrBhvThp(ind,1,2),FSrBhvThp(ind,2,2),FSrBhvThp(ind,3,2),20,pfs.adata.bins{2}(csmPhzPrefMax(ind,2)),'filled');
% $$$ colormap(gca,'hsv');
% $$$ caxis([-pi,pi]);
% $$$ 
% $$$ figure,
% $$$ labels = {'low','rear','high'};
% $$$ for v = 1:3,
% $$$ subplot(1,3,v);
% $$$ hist(pfs.adata.bins{2}(csmPhzPrefMax(FSrBhvThp(uind,v,1)>0.5,1)),linspace(-pi,pi,16))
% $$$ title(labels{v});
% $$$ end




% 1 find the index of minima for each phase
% $$$ [mcv,mcvi] = max(sq(-cvals(:,:,1,:)),[],2);
% $$$ 
% $$$ 
% $$$ 
% $$$ scvals = zeros(size(cvals));
% $$$ scmins = zeros(size(cmins));
% $$$ for u = 1:size(cvals,1),
% $$$ for p = 1:pfs.adata.binSizes(2),
% $$$ [~,scvmi] = sort(-cvals(:,:,1,p),2,'descend');    
% $$$ scvals(u,:,1,p) = -cvals(u,scvmi(u,:,1,1),1,p);
% $$$ scmins(u,:,:,p) = cmins(u,scvmi(u,:,1,1),:,p); 
% $$$ end
% $$$ end
% $$$ 
% $$$ % find bhv fields with only one main mode
% $$$ u1 = all(sum(scvals(uind,:,1,:)~=0,2)==1,4);
% $$$ % find bhv fields with two modes 
% $$$ u2 = all(sum(scvals(uind,:,1,:)~=0,2)==2,4);
% $$$ 
% $$$ figure,
% $$$ plot3(reshape(scmins(uind,1,1,:),[],1)+randn([sum(uind)*16,1])/5,...
% $$$       reshape(scmins(uind,1,2,:),[],1)+randn([sum(uind)*16,1])/5,...
% $$$       reshape(eSpi(uind,1,:),[],1)+randn([sum(uind)*16,1])/15,'.');
% $$$ 
% $$$ 
% $$$ p =  5
% $$$ figure,
% $$$ plot3(reshape(scmins(uind,1,1,p),[],1)+randn([sum(uind),1])/5,...
% $$$       reshape(scmins(uind,1,2,p),[],1)+randn([sum(uind),1])/5,...
% $$$       reshape(eSpi(uind,1,p),[],1)+randn([sum(uind),1])/15,'.');
% $$$ %      reshape(scvals(uind,1,1,:),[],1)+randn([sum(uind)*16,1])/5,'.');
% $$$ 
% $$$ 
% $$$ 
% $$$ uscmins = scmins(uind,:,:,:);
% $$$ uscvals = scvals(uind,:,:,:);
% $$$ 
% $$$ tind = 20;
% $$$ unit = 21;
% $$$ cind = ismember(cluSessionMap,[tind,unit],'rows');;             
% $$$ figure();
% $$$ u = find(cind(uind));
% $$$ hold('on');
% $$$ plot3(pfs.adata.bins{1}(nonzeros(uscmins(u,:,1,:))),...
% $$$       pfs.adata.bins{3}(nonzeros(uscmins(u,:,2,:))),...
% $$$       nonzeros(uscvals(u,:,1,:)),'*');
% $$$ 
% $$$ for i = 1:3,
% $$$     gvals = sq(uscvals(u,i,1,:)~=0);
% $$$     
% $$$     if sum(gvals)>0,
% $$$         plot3(pfs.adata.bins{1}(reshape(uscmins(u,i,1,:),[],1)),...
% $$$               pfs.adata.bins{3}(reshape(uscmins(u,i,2,:),[],1)),...
% $$$               reshape(uscvals(u,i,1,:),[],1),'*-');
% $$$     end
%end
% $$$ xlim([-2,0.5]);
% $$$ ylim([-0.5,1.5]);
% $$$ zlim([0,max(zlim)]);
% $$$ 
% $$$ peaks = [pfs.adata.bins{1}(nonzeros(uscmins(u,:,1,:))),...
% $$$           pfs.adata.bins{3}(nonzeros(uscmins(u,:,2,:)))];
% $$$ pkvls = nonzeros(reshape(uscvals(u,:,1,:),[],1));
% $$$ kls = kmeans(peaks,3);
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ plot3(peaks(kls==1,1),peaks(kls==1,2),pkvls(kls==1),'.')
% $$$ plot3(peaks(kls==2,1),peaks(kls==2,2),pkvls(kls==2),'.')
% $$$ plot3(peaks(kls==3,1),peaks(kls==3,2),pkvls(kls==3),'.')
% $$$ xlim([-2,0.5]);
% $$$ ylim([-0.5,1.5]);
% $$$ zlim([0,max(zlim)]);
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ usd1 = sq(sqrt(sum((scmins(uind,1,:,:)-circshift(scmins(uind,1,:,:),-1,3)).^2,3)));
% $$$ usd2 = sq(sqrt(sum((scmins(uind,2,:,:)-circshift(scmins(uind,2,:,:),-1,3)).^2,3)));
% $$$ usd12 = sq(sqrt(sum((scmins(uind,1,:,:)-circshift(scmins(uind,2,:,:),-1,3)).^2,3)));
% $$$ 
% $$$ figure,plot3(usd1(1,:),usd2(1,:),sq(uscvals(1,1,1,:)),'*-');
% $$$ figure,imagesc(u1d)
% $$$ 
% $$$ 
% $$$ [mcv,mcvipp] = max(sq(mcv),[],2);
% $$$ [bsi,bsipp] = max(sq(eSpi),[],2);
% $$$ 
% $$$ figure();
% $$$ hist2([bsipp(uind),mcvipp(uind)]+0.5,1:17,1:17);
% $$$ set(gca,'XTick',1:2:17)
% $$$ set(gca,'YTick',1:2:17)
% $$$ set(gca,'XTickLabels',round(linspace(-pi,pi,8),2));
% $$$ set(gca,'YTickLabels',round(linspace(-pi,pi,8),2));
% $$$ axis('xy');
% $$$ xlabel('BHV Info PP');
% $$$ ylabel('Max Rate PP');
% $$$ 
% $$$ [phzVals,phzPref] = max(sq(-cvals(:,1,1,:)),[],2);
% $$$ 
% $$$ muind = uind(ismember(cluSessionMap(:,1),17:23));
% $$$ myUDP = unitDepths(muind);
% $$$ 
% $$$ hpm = sq(cmins(:,1,1,:));
% $$$ hpm = hpm(sub2ind([size(cmins,1),10],(1:size(cmins,1))',phzPref));
% $$$ hpm = hpm(uind);
% $$$ 
% $$$ fhpm = sq(cmins(:,1,1,:));
% $$$ fhpm = fhpm(sub2ind([size(cmins,1),10],(1:size(cmins,1))',phzPref-4+10*double((phzPref-4)<=0)));
% $$$ %fhpm = fhpm(sub2ind([size(cmins,1),10],(1:size(cmins,1))',phzPref+3-10*double((phzPref+3)>=11)));
% $$$ fhpm = fhpm(uind);
% $$$ 
% $$$ bpm = sq(cmins(:,1,2,:));
% $$$ bpm = bpm(sub2ind([size(cmins,1),10],(1:size(cmins,1))',phzPref));
% $$$ bpm = bpm(uind);
% $$$ 
% $$$ 
% $$$ [bsi,bsipp] = max(sq(eSpi),[],2);
% $$$ tcvals = -cvals(:,1,1,:);
% $$$ bsippVals = tcvals(sub2ind([size(cmins,1),10],(1:size(cmins,1))',bsipp));
% $$$ bsippVals(uind);
% $$$ bsi = bsi(uind);
% $$$ bsipp = bsipp(uind);
% $$$ 


% ATTEMPT to regress pca components ------------ STATUS:MEH ------------
% $$$ t = 20;
% $$$ Trial = Trials{t}; 
% $$$ unitSubset = units{t};
% $$$ pfs = MTAApfs(Trial,'tag','hbpptbpFS1v3');
% $$$ pft = pfs_2d_theta(Trial);
% $$$ 
% $$$ 
% $$$ %imagesc(pfs.adata.bins{[1,3]},angle(ppmap)');
% $$$ 
% $$$ % VISUALIZE eigenvectors
% $$$ 
% $$$ clear('eigVecs');
% $$$ nfac = 5;
% $$$ hfig = figure();
% $$$ u = unitSubset(1);
% $$$ while u~=-1,
% $$$     hax = tight_subplot(2,pfs.adata.binSizes(2)+1,0,0.1,0.01);    
% $$$     axes(hax(1));
% $$$     plot(pft,u,'mean','text');
% $$$     title(['Unit: ',num2str(u)]);
% $$$ 
% $$$     rmap = plot(pfs,u,1,'colorbar',[],false);
% $$$     rmapNind = sq(sum(reshape(~isnan(rmap(:)),size(rmap)),2))>=pfs.adata.binSizes(2);
% $$$     zmap = reshape(permute(rmap,[1,3,2]),[],pfs.adata.binSizes(2));
% $$$     zmap(repmat(~rmapNind(:),[1,pfs.adata.binSizes(2)])) = 0;            
% $$$     zmap(zmap(:)<0) = 0;            
% $$$     rmapNind = rmapNind&reshape(sum(reshape(zmap,[],pfs.adata.binSizes(2))','omitnan')~=0,size(rmap,1),size(rmap,3));
% $$$ 
% $$$     % COMPUTE erpPCA 
% $$$     pmap = zmap;
% $$$     pmap(isnan(pmap(:,1))|sum(pmap','omitnan')'==0,:) = [];
% $$$     if size(pmap,1)>20    
% $$$         mrate= prctile(rmap(:),99);
% $$$         [LU,LR,FSr,VT] = erpPCA(pmap',nfac);
% $$$         for v = 1:nfac,
% $$$             evec = nan([size(rmap,1),size(rmap,3)]);
% $$$             evec(rmapNind) = LR(:,v);%*(-double(max(FSr(:,4))<abs(min(FSr(:,4)))));
% $$$             eigVecs(ucnt,v,:,:) = evec;
% $$$         end
% $$$ 
% $$$         for p = 1:pfs.adata.binSizes(2),
% $$$             axes(hax(p+1));
% $$$             hold('on');
% $$$             imagesc(pfd{1}.adata.bins{:},sq(rmap(:,p,:))');
% $$$             %imagesc(pfd{1}.adata.bins{:},sq(rmap(:,p,:))'-sq(mean(rmap(:,:,:),2,'omitnan'))');    
% $$$             axis('xy');
% $$$             axis('tight');
% $$$             caxis([0,mrate]);
% $$$             %caxis([-mrate/2,mrate/2]);
% $$$         end
% $$$         axes(hax(pfs.adata.binSizes(2)+2));
% $$$         plot(VT(:,4),'-+');
% $$$         for v = 1:nfac,  
% $$$             axes(hax(pfs.adata.binSizes(2)+2+v));
% $$$             hold('on');
% $$$             imagesc(pfs.adata.bins{[1,3]},sq(eigVecs(ucnt,v,:,:))');
% $$$             axis('xy');
% $$$             axis('tight');
% $$$         end
% $$$         axes(hax(pfs.adata.binSizes(2)+2+v+1));
% $$$ 
% $$$         rmap = reshape(permute(rmap,[1,3,2]),[],pfs.adata.binSizes(2));;
% $$$         ppmap = sum(bsxfun(@times,rmap,exp(i.*pfs.adata.bins{2}')),2)./sum(rmap,2);
% $$$         ppmap = reshape(ppmap,pfs.adata.binSizes([1,3])');
% $$$         imagescnan({pfs.adata.bins{[1,3]},angle(ppmap)'},[-pi,pi],'circular',true,'colorMap',@hsv);
% $$$         axis('xy');
% $$$         axes(hax(pfs.adata.binSizes(2)+2+v+3));
% $$$         imagescnan({pfs.adata.bins{[1,3]},abs(ppmap)'},[],'circular',true,'colorMap',@hsv);
% $$$         axis('xy');
% $$$         axes(hax(pfs.adata.binSizes(2)+2+v+5));
% $$$         imagescnan({1:size(FSr,1),1:size(FSr,2),FSr'},[],'circular',true,'colorMap',@jet);
% $$$     end
% $$$     u = figure_controls(hfig,u,unitSubset,false,[],[]);
% $$$     clf;
% $$$ end
% END Data exploration figure


% FIGURE : Behavioral theta phase decomposition Examples -------------------------------------------%

% bhvMulti  mod 52, 60, 61, 80,
% bhvSingle mod 81, 86, 103, 104 
% bhvProg   mod 63, 34, 83, 72, 145
% $$$ 
% $$$ req20190122();% -> bhvInfoTS
% $$$ 
% $$$ unitsExampleSubset = [81,86,104];
% $$$ unitsExampleSubset = [52,61,72];
% $$$ unitsExampleSubset = [63, 34, 83];
% $$$ 
% $$$ nx = pfs.adata.binSizes(2)/2+2;
% $$$ ny = numel(unitsExampleSubset);
% $$$ 
% $$$ hfig = figure();
% $$$ fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters','FontSize',8,'LineWidth',1);
% $$$ 
% $$$ hfig.Units = 'centimeters';
% $$$ set(hfig,'Position',[0,0,18,6]);
% $$$ hfig.Units = 'normalized';
% $$$ hax = tight_subplot(ny,nx,0,0.1,0.1);
% $$$ hfig.Units = 'centimeters';
% $$$ af(@(h) set(h,'Units','centimeters'), hax);
% $$$ af(@(h) set(h,'FontSize',8), hax);
% $$$ set(hfig,'Position',[0,0,18,10]);
% $$$ set(fax,'Position',[0,0,18,10]);
% $$$ phzOrder = [10,12,14,16,2,4,6,8];
% $$$ phzLables = {'','90','','180','','270','',''};
% $$$ 
% $$$ for y = 1:ny,
% $$$ % GET space and theta phase jointly restricted bhv rate maps
% $$$     rmap = plot(pfs,unitsExampleSubset(y),1,'colorbar',[],false);
% $$$     rmapNind = sq(sum(reshape(~isnan(rmap(:)),size(rmap)),2))>=pfs.adata.binSizes(2);
% $$$     zmap = reshape(permute(rmap,[1,3,2]),[],pfs.adata.binSizes(2));
% $$$     zmap(repmat(~rmapNind(:),[1,pfs.adata.binSizes(2)])) = 0;            
% $$$     zmap(zmap(:)<0) = 0;            
% $$$     rmapNind = rmapNind&reshape(sum(reshape(zmap,[],pfs.adata.binSizes(2))','omitnan')~=0,...
% $$$                                 size(rmap,1),size(rmap,3));      
% $$$     mrate= prctile(rmap(:),99.9);    
% $$$     
% $$$ % PLOT theta restricted spatial rate maps
% $$$     axes(hax((y-1)*nx+1));
% $$$     hx = gca();
% $$$     plot(pft,unitsExampleSubset(y),'mean','text',mrate,'colorMap',@jet);
% $$$     hx.Position(1) = hx.Position(1)-1;    
% $$$     hx.Position(2) = hx.Position(2)+(ny-y).*0.5;    
% $$$     hx.Position(4) = hx.Position(3);
% $$$     ylabel(['Unit: ',num2str(unitsExampleSubset(y))]);
% $$$     set(hx,'YTickLabels',{});
% $$$     set(hx,'XTickLabels',{});    
% $$$     if y == 1,
% $$$         title({'Spatial','Rate Map'});
% $$$     end
% $$$     
% $$$ % PLOT spatially restricted bhv rate maps
% $$$     axes(hax((y-1)*nx+2));
% $$$     hx = gca();
% $$$     plot(pfd{20,1},unitsExampleSubset(y),'mean','text',mrate,false,'colorMap',@jet);
% $$$     hx.Position(1) = hx.Position(1)-0.25;    
% $$$     hx.Position(2) = hx.Position(2)+(ny-y).*0.5;    
% $$$     hx.Position(4) = hx.Position(3);    
% $$$     set(hx,'YTickLabels',{});
% $$$     set(hx,'XTickLabels',{});    
% $$$     if y == 1,
% $$$         ylabel('Body Pitch');
% $$$         xlabel('Head Pitch');
% $$$         title({'Behavioral','Rate Map'});
% $$$     end
% $$$     
% $$$ % PLOT space and theta phase jointly restricted bhv rate maps
% $$$     for p = 1:(pfs.adata.binSizes(2)/2),
% $$$         axes(hax((y-1)*nx+p+2));
% $$$         hx = gca();
% $$$         hold(hx,'on');
% $$$         imagescnan({pfd{1}.adata.bins{:},sq(rmap(:,phzOrder(p),:))'},[0,mrate],[],false,[0,0,0],'colorMap',@jet);
% $$$         axis('xy');
% $$$         axis('tight');
% $$$         hx.Position(2) = hx.Position(2)+(ny-y).*0.5;
% $$$         hx.Position(4) = hx.Position(3);
% $$$         hx.Position(1) = hx.Position(1)+0.5;
% $$$     end 
% $$$    text(pfs.adata.bins{1}(end)-0.5*diff(pfs.adata.bins{1}([1,end])),...
% $$$          pfs.adata.bins{3}(end)-0.10*diff(pfs.adata.bins{3}([1,end])),...
% $$$          sprintf('%2.1f',mrate),...
% $$$         'Color','w','FontWeight','bold','FontSize',8)
% $$$ 
% $$$ end
% $$$ 
% $$$ % CREATE reference bars
% $$$ axes(fax)
% $$$ xlim([0,fax.Position(3)]);
% $$$ ylim([0,fax.Position(4)]);
% $$$ % SUPERAXIS Theta phase 
% $$$ line([hax(3).Position(1),sum(hax(nx).Position([1,3]))],...
% $$$      sum(hax(3).Position([2,4])).*[1,1]+0.5,'LineWidth',1)
% $$$ text(mean([hax(3).Position(1),sum(hax(nx).Position([1,3]))]),...
% $$$      sum(hax(3).Position([2,4]))+0.75,'Theta Phase','HorizontalAlignment','center','FontSize',8)
% $$$ % TICKS theta phase
% $$$ text(hax(3).Position(1),sum(hax(3).Position([2,4]))+0.25,'0^{o}','HorizontalAlignment','center','FontSize',8)
% $$$ text(mean([hax(3).Position(1),sum(hax(nx).Position([1,3]))]),...
% $$$      sum(hax(3).Position([2,4]))+0.25,'180^{o}','HorizontalAlignment','center','FontSize',8)
% $$$ text(sum(hax(nx).Position([1,3])),...
% $$$      sum(hax(3).Position([2,4]))+0.25,'360^{o}','HorizontalAlignment','center','FontSize',8)
% $$$ % REFERENCE bar space
% $$$ line([hax(2+nx*2).Position(1),sum(hax(2+nx*2).Position([1,3]).*[1,11/28])],...
% $$$      hax(2+nx*2).Position(2).*[1,1]-0.1,'LineWidth',1);
% $$$ text(hax(2+nx*2).Position(1),hax(2+nx*2).Position(2)-0.3,'1 rad','FontSize',8);
% $$$ % REFERENCE bar space
% $$$ line([hax(1+nx*2).Position(1),sum(hax(1+nx*2).Position([1,3]).*[1,0.5])],...
% $$$      hax(1+nx*2).Position(2).*[1,1]-0.1,'LineWidth',1);
% $$$ text(hax(1+nx*2).Position(1),hax(1+nx*2).Position(2)-0.3,'50 cm','FontSize',8);










% $$$ u = 61;    
% $$$ % RESTRICT map elements to bins with no nans across theta bins
% $$$ rmap = plot(pfs,u,'mean','colorbar',[],false,0.5);            
% $$$ rmapNind = sq(sum(reshape(~isnan(rmap(:)),size(rmap)),2))>=pfs.adata.binSizes(2);
% $$$ %zmap = permute(rmap,[1,3,2]);
% $$$ zmap = reshape(permute(rmap,[1,3,2]),[],pfs.adata.binSizes(2));
% $$$ % $$$             zmap(repmat(rmapNind,[1,1,pfs.adata.binSizes(2)]) & isnan(zmap)) = 0;    
% $$$ % $$$             zmap(repmat(~rmapNind,[1,1,pfs.adata.binSizes(2)])) = 0;
% $$$ zmap(repmat(~rmapNind(:),[1,pfs.adata.binSizes(2)])) = 0;            
% $$$ zmap(zmap(:)<0) = 0;            
% $$$ rmapNind = rmapNind&reshape(sum(reshape(zmap,[],pfs.adata.binSizes(2))','omitnan')~=0,size(rmap,1),size(rmap,3));
% $$$ 
% $$$ 
% $$$ % COMPUTE erpPCA 
% $$$ pmap = zmap;
% $$$ pmap(isnan(pmap(:,1))|sum(pmap','omitnan')'==0,:) = [];
% $$$ 
% $$$ [LU,LR,FSr,VT] = erpPCA(pmap',6);
% $$$ for v = 1:6,
% $$$     evec = nan([size(rmap,1),size(rmap,3)]);
% $$$     evec(rmapNind) = LR(:,v);
% $$$     eigVecs(ucnt,v,:,:) = evec;
% $$$ end
% $$$ 
% $$$ % VISUALIZE eigenvectors
% $$$ figure,
% $$$ hax = tight_subplot(2,pfs.adata.binSizes(2)+1,0,0.1,0.01);
% $$$ axes(hax(1));
% $$$ plot(pft,u);
% $$$ for p = 1:pfs.adata.binSizes(2),
% $$$     %zmapNormalizedPart = zmap(:,p)./sum(zmap(:,p));                
% $$$     %mBhvPos = [sum(hp.*zmapNormalizedPart),sum(bp.*zmapNormalizedPart)];
% $$$     axes(hax(p+1));
% $$$     hold('on');
% $$$     imagesc(pfd{1}.adata.bins{:},sq(rmap(:,p,:))');
% $$$     axis('xy');
% $$$     axis('tight');
% $$$     %plot(mBhvPos(1),mBhvPos(2),'*m');                
% $$$ end
% $$$ axes(hax(pfs.adata.binSizes(2)+2));
% $$$ plot(VT(:,4),'-+');
% $$$ for v = 1:6,  
% $$$     axes(hax(pfs.adata.binSizes(2)+2+v));
% $$$     hold('on');
% $$$     imagesc(pfs.adata.bins{[1,3]},sq(eigVecs(ucnt,v,:,:))');
% $$$     axis('xy');
% $$$     axis('tight');
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,
% $$$ subplot(131);
% $$$ evec = nan([size(rmap,1),size(rmap,3)]);
% $$$ evec(rmapNind) = efit(:,1);
% $$$ imagesc(pfs.adata.bins{[1,3]},evec');
% $$$ axis('xy');
% $$$ subplot(132);
% $$$ evec = nan([size(rmap,1),size(rmap,3)]);
% $$$ evec(rmapNind) = pmap(:,p);
% $$$ imagesc(pfs.adata.bins{[1,3]},evec');
% $$$ axis('xy');
% $$$ subplot(133);
% $$$ evec = nan([size(rmap,1),size(rmap,3)]);
% $$$ evec(rmapNind) = pmap(:,p)-efit;
% $$$ imagesc(pfs.adata.bins{[1,3]},evec');
% $$$ axis('xy');
% $$$ 
% $$$ 
% $$$ figure,
% $$$ subplot(121);
% $$$ evec = nan([size(rmap,1),size(rmap,3)]);
% $$$ evec(rmapNind) = PCAScores(:,1);
% $$$ imagesc(pfs.adata.bins{[1,3]},evec');
% $$$ subplot(122);
% $$$ evec = nan([size(rmap,1),size(rmap,3)]);
% $$$ evec(rmapNind) = PCAScores(:,2);
% $$$ imagesc(pfs.adata.bins{[1,3]},evec');
% $$$ 
% $$$ figure();
% $$$ hax = tight_subplot(1,pfs.adata.binSizes(2),0,0.1,0);
% $$$ for p = 1:pfs.adata.binSizes(2),
% $$$     axes(hax(p));
% $$$     evec = nan([size(rmap,1),size(rmap,3)]);
% $$$     evec(rmapNind) = pmap(:,p)-yfitPCR;
% $$$     imagesc(pfs.adata.bins{[1,3]},evec');
% $$$     axis('xy');
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ t = 20;
% $$$ Trial = Trials{t}; 
% $$$ unitSubset = units{t};
% $$$ pfs = MTAApfs(Trial,'tag','hbpptbpFS1');
% $$$ 
% $$$ % Plot bhv fields and their maxima
% $$$ figure();
% $$$ for u = 1:numel(unitSubset)
% $$$ rmap = plot(pfs,unitSubset(u),'','',[],false);
% $$$ rmax = prctile(rmap(nniz(rmap(:))),99);
% $$$ subplot(3,8,1);
% $$$ plot(pft,unitSubset(u),'mean','colorbar',[0,rmax],true);
% $$$ subplot(3,8,2);
% $$$ plot(pfd{1},unitSubset(u),'mean','colorbar',[0,rmax],false);
% $$$ for p = 1:pfs.adata.binSizes(2),
% $$$     mrate = mean(reshape(sq(rmap(:,p,:)),[],1),'omitnan');
% $$$     subplot(3,8,p+8);
% $$$     %imagescnan({pfs.adata.bins{[1,3]},sq(rmap(:,p,:))'-mrate},[0,mrate]);    
% $$$     imagescnan({pfs.adata.bins{[1,3]},sq(rmap(:,p,:))'},[0,rmax]);
% $$$     axis('xy');
% $$$     hold('on');
% $$$     for m = 1:sum(smins(u,:,1,p)~=0)
% $$$     %for m = 1:sum(sq(peakPos(u,p,1,:))~=0)
% $$$         %plot(peakPos(u,p,1,m),peakPos(u,p,2,m),'*m');
% $$$         plot(smins(u,m,1,p),smins(u,m,2,p),'*m');
% $$$         plot(pfs.adata.bins{1}(cmins(u,m,1,p)),...
% $$$              pfs.adata.bins{3}(cmins(u,m,2,p)),...
% $$$              '*g');        
% $$$     end
% $$$ end
% $$$ waitforbuttonpress();
% $$$ end




% TESTING GAUSSIAN FITTING -------------- STATUS:FAILED --------------------
% $$$ bins = linspace(-400,400,spcBinCnt);
% $$$ 
% $$$ % GET centers of the bins 
% $$$ binc = {(pfs.adata.bins{1}(1:end-1)+pfs.adata.bins{1}(2:end))./2, ...
% $$$         (pfs.adata.bins{3}(1:end-1)+pfs.adata.bins{3}(2:end))./2, };
% $$$ binc = pfs.adata.bins([1,3]);
% $$$ 
% $$$ % GET grid coordinates of all bins
% $$$ gridBins = cell([1,2]);
% $$$ [gridBins{:}] = ndgrid(binc{:});
% $$$ gridBins = reshape(cat(3,gridBins{:}),[],2);
% $$$ 
% $$$ 
% $$$ unit = unitSubset(4);
% $$$ 
% $$$ rmap = plot(pfs,unit,'','',[],false);
% $$$ 
% $$$ center = [0,0];
% $$$ weightsFunction = @(x,center) 1./sqrt(sum(bsxfun(@minus,x,center).^2,2)+1);
% $$$ 
% $$$ % TEST weight function 
% $$$ weights = weightsFunction(gridBins,[0,0]);
% $$$ weights = weights./max(weights(:));
% $$$ figure,imagesc(binc{:},reshape(weights',pfs.adata.binSizes([1,3])')');axis('xy');
% $$$ % TEST RMAP
% $$$ zmap = sq(rmap(:,7,:));
% $$$ figure,imagesc(binc{:},reshape(zmap,pfs.adata.binSizes([1,3])')');axis('xy');
% $$$ 
% $$$ 
% $$$ 
% $$$ % Gaussian fit 
% $$$ gaussianFitFunction = fittype( @(A,xa,ya,xya,xo,yo,x,y) ...
% $$$                                A.*exp(-(xa.*(x-xo).^2+xya.*(x-xo).*(y-yo)+ya.*(y-yo).^2)),...
% $$$                               'independent',{'x', 'y'},'dependent', 'z'); 
% $$$ gaussianFitStartPoint = [5,1,1,1,0,0];
% $$$ 
% $$$ 
% $$$ fo{p} = fit(gridBins,outIn(:),fitModel,'StartPoint',fitStartPoint);axis('xy');
% $$$ outModel = reshape(fo{p}(gridBins(:,1),gridBins(:,2)),size(outIn));
% $$$ 
% $$$ peakPos = zeros([numel(unitSubset),pfs.adata.binSizes(2),2,5]);
% $$$ 
% $$$ % for u = 1:numel(unitSubset)
% $$$ u = 3;
% $$$ rmap = plot(pfs,unitSubset(u),'','',[],false);
% $$$ % TEST Gaussian fit
% $$$ for p = 1:pfs.adata.binSizes(2),
% $$$     zmap = sq(rmap(:,p,:));
% $$$     mrate = mean(reshape(zmap(nind),[],1),'omitnan');
% $$$     fo = {};
% $$$     nind = nniz(zmap(:));
% $$$     for x = 1:pfs.adata.binSizes(1),
% $$$         for y = 1:pfs.adata.binSizes(3),            
% $$$             gaussianFitStartPoint = [mrate,1,1,1,pfs.adata.bins{1}(x),pfs.adata.bins{1}(y)];
% $$$             try,
% $$$                 fo{x,y} = fit(gridBins(nind,:),zmap(nind)-mrate,gaussianFitFunction,'StartPoint', ...
% $$$                               gaussianFitStartPoint);
% $$$             end
% $$$         end                                        
% $$$     end
% $$$     pos = cell2mat(cf(@(fo) cat(3,fo.xo,fo.yo), fo(~isempty(fo))));
% $$$     pos = round(reshape(pos,[],2),2);    
% $$$     pos(sqrt(sum(pos.^2,2))>4,:) = nan;                                
% $$$     tPeakPos = unique(pos,'rows');
% $$$     peakPos(u,p,:,1:sum(nniz(tPeakPos))) = tPeakPos(nniz(tPeakPos),:)';
% $$$ end
% $$$ amps = cell2mat(cf(@(fo) fo.A, fo));
% $$$ % CHECK distribution of amps
% $$$ figure,imagesc(binc{:},log10(amps)');axis('xy');
% $$$ pos = cell2mat(cf(@(fo) cat(3,fo.xo,fo.yo), fo));
% $$$ pos = reshape(pos,[],2);
% $$$ pos(sqrt(sum(pos.^2,2))>4) = nan;
% $$$ figure,plot(pos(:,1),pos(:,2),'.')
% $$$ figure,hist2(pos,100,100)



% FIGURES -----------------------------------------------------------------------------
% $$$ t = 20;
% $$$ Trial = Trials{t}; 
% $$$ unitSubset = units{t};
% $$$ pfs = MTAApfs(Trial,'tag','hbpptbpFS1');
% $$$ pfd = req20180123_ver5(Trial,[],'10');
% $$$ pft = pfs_2d_theta(Trial);
% $$$ unit = unitSubset(1);
% $$$ 
% $$$ figure,
% $$$ hold('on');
% $$$ plot(pfd{1},unit,'mean','colorbar',[],false);
% $$$ %uind = ismember(cluSessionMap,[t,unit],'rows');
% $$$ uind = ismember(cluSessionMap(:,1),17:23);
% $$$ uind = ismember(cluSessionMap(:,1),[3:5,17:23]);
% $$$ 
% $$$ % SVD  
% $$$ [U,S,V] = (sq(eSpi(uind,:,:,:)));
% $$$ figure();imagesc(V);
% $$$ pc2 = V(:,2)'*sq(eSpi(uind,:,:))';
% $$$ pc3 = V(:,3)'*sq(eSpi(uind,:,:))';
% $$$ figure();plot(pc2,pc3,'.');
% $$$ 
% $$$ 
% $$$ [Coeff,Scr,Latent] = pca(sq(eSpi(uind,:,:,:)));
% $$$ figure();
% $$$ subplot(131),imagesc(Coeff');
% $$$ subplot(132),imagesc(Scr');
% $$$ subplot(133),plot(Latent,'+-')
% $$$ 
% $$$ 
% $$$ [Coeff,Scr,Latent] = pca(sq(eSpi(uind,:,:,:)));
% $$$ figure();
% $$$ subplot(131),imagesc(Coeff');
% $$$ subplot(132),imagesc(Scr');
% $$$ subplot(133),plot(Latent,'+-')
% $$$ 
% $$$ 
% $$$ [Coeff,Scr,Latent] = princomp(bsxfun(@minus,sq(eSpi(uind,:,:,:)),mean(sq(eSpi(uind,:,:,:)),2)));
% $$$ figure();
% $$$ subplot(141),imagesc(Coeff');
% $$$ subplot(142),imagesc(Scr');
% $$$ subplot(143),plot(Latent,'+-')
% $$$ 
% $$$ figure,hist(Scr(:,1),20);
% $$$ 
% $$$ [Coeff,Scr,Latent] = princomp(bsxfun(@minus,sq(mcv(uind,:,:,:)),mean(sq(mcv(uind,:,:,:)),2)));
% $$$ figure();
% $$$ subplot(141),imagesc(Coeff');
% $$$ subplot(142),imagesc(Scr');
% $$$ subplot(143),plot(Latent,'+-')
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure(),imagesc(bsxfun(@rdivide,bsxfun(@minus,sq(eSpi(uind,:,:,:)),mean(sq(eSpi(uind,:,:,:)),2)),std(sq(eSpi(uind,:,:,:)),[],2)))
% $$$ 
% $$$ 
% $$$ 
% $$$ [~,mcvi] = max(sq(cvals(uind,:,1,:)),[],2);
% $$$ 
% $$$ figure();plot(Scr(:,1),Scr(:,2),'.');
% $$$ 
% $$$ 
% $$$ 
% $$$ pc2 = V(:,2)'*sq(eSpi(uind,:,:))';
% $$$ pc3 = V(:,3)'*sq(eSpi(uind,:,:))';
% $$$ figure();plot(pc2,pc3,'.');





% $$$ plot(pfd{1}.adata.bins{1}(nonzeros(reshape(cmins(uind,:,1,2),[],1))),...
% $$$      pfd{1}.adata.bins{2}(nonzeros(reshape(cmins(uind,:,2,2),[],1))),...
% $$$      '.')
% $$$ axis('xy');


% $$$ unit = unitSubset(1);
% $$$ hfig = figure();
% $$$ while unit ~=-1;
% $$$ 
% $$$ uind = ismember(cluSessionMap,[t,unit],'rows');
% $$$ 
% $$$ rmap = plot(pfs,unit,'mean','colorbar',[],false,0.5);
% $$$             rmapNind = sq(sum(reshape(~isnan(rmap(:)),size(rmap)),2))>=10;
% $$$             zmap = permute(rmap,[1,3,2]);
% $$$             %zmap = reshape(permute(rmap,[1,3,2]),[],10);
% $$$             zmap(repmat(rmapNind,[1,1,10]) & isnan(zmap)) = 0;    
% $$$             zmap(repmat(~rmapNind,[1,1,10])) = 0;
% $$$             %zmap(repmat(~rmapNind(:),[1,10])) = 0;            
% $$$             zmap(zmap(:)<0) = 0;            
% $$$ clim = [0,max(rmap(:))];
% $$$ subplot(2,14,[1:3]);
% $$$ plot(pft,unit,'mean','colorbar',clim);
% $$$ subplot(2,14,4);
% $$$ plot(pfd{1},unit,'mean','',clim,false);
% $$$ for p = 1:10
% $$$     subplot(2,14,p+4);hold('on');
% $$$     imagescnan({pfd{1}.adata.bins{:},sq(rmap(:,p,:))'},clim);axis('xy');axis('tight');
% $$$     for m = 1:sum(all(cmins(uind,:,:,p),3)),
% $$$         plot(pfd{1}.adata.bins{1}(cmins(uind,m,1,p)),...
% $$$              pfd{1}.adata.bins{2}(cmins(uind,m,2,p)),...
% $$$              '*m');
% $$$     end    
% $$$ % $$$     subplot(2,14,p+4+14);hold('on');    
% $$$ % $$$     
% $$$ % $$$     imagescnan({pfd{1}.adata.bins{:},convn(zmap(:,:,p),Smoother,'same')'},clim);axis('xy');axis('tight');
% $$$ end
% $$$             
% $$$ 
% $$$             
% $$$ 
% $$$             subplot(2,14,19);plot(sq(iva(uind,1,:)),'-+');
% $$$             subplot(2,14,20);plot(sq(iva(uind,1,:)),sq(eSpi(uind,:,:,:))','-+');
% $$$             subplot(2,14,21);plot(sq(iva(uind,1,:)),-sq(cvals(uind,1,1,:))','-+');
% $$$ 
% $$$ subplot(2,14,18);
% $$$ plot(sq(eSpi(uind,:,:,:)),'+-');
% $$$ unit = figure_controls(hfig,unit,unitSubset,false,[],[]);
% $$$ end
% $$$ 
% $$$ figure,plot(bsxfun(@rdivide,sq(eSpi(uind,:,:))',mean(sq(eSpi(uind,:,:)),[],2)'))
% $$$ figure,plot(bsxfun(@rdivide,sq(eSpi(uind,:,:))',mean(sq(eSpi(uind,:,:)),2)'))
% $$$ figure,plot(mean(bsxfun(@rdivide,sq(eSpi(uind,:,:))',mean(sq(eSpi(uind,:,:)),2)'),2))
% $$$ 
% $$$ figure,imagesc(sq(cvals(:,1,1,:))');axis('xy');
% $$$ 
% $$$ ncval = sq(bsxfun(@rdivide,-cvals(:,1,1,:),max(-cvals(:,1,1,:),[],4)))';
% $$$ 
% $$$ uind = ismember(cluSessionMap(:,1),17:23);
% $$$ seslins = find(diff(cluSessionMap(ismember(cluSessionMap(:,1),17:23),1))==1);
% $$$ figure,imagesc(ncval([6:10,1:5],uind));axis('xy');
% $$$ Lines(seslins,[],'r',[],4)
% $$$ caxis([0,1])
% $$$ 
% $$$ 
% $$$ cf(@(t) t.load('nq'), Trials(17:23));
% $$$ nq = cf(@(t) t.nq, Trials(17:23));
% $$$ 
% $$$ allUnitDepths = cf(@(t,n,a)  n.maxAmpChan-a(n.ElNum)',  Trials(17:23),nq,anatGrpSWPCenter);
% $$$ unitDepths   = cf(@(a,u)  a(u),                       allUnitDepths,units(17:23));
% $$$ unitDepths = cat(1,unitDepths{:});
% $$$ %elnum        = cf(@(t,u)  t.spk.map(u,2),             Trial(17:23),units(17:23));
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ uind = ismember(cluSessionMap(:,1),1:23);
% $$$ 
% $$$ uind = ismember(cluSessionMap(:,1),[17:23]);
% $$$ 
% $$$ % 1 find the index of minima for each phase
% $$$ [mcv,mcvi] = max(sq(-cvals(:,:,1,:)),[],2);
% $$$ [mcv,mcvipp] = max(sq(mcv),[],2);
% $$$ [bsi,bsipp] = max(sq(eSpi),[],2);
% $$$ 
% $$$ figure();
% $$$ hist2([bsipp(uind),mcvipp(uind)]+0.5,1:17,1:17);
% $$$ set(gca,'XTick',1:2:17)
% $$$ set(gca,'YTick',1:2:17)
% $$$ set(gca,'XTickLabels',round(linspace(-pi,pi,8),2));
% $$$ set(gca,'YTickLabels',round(linspace(-pi,pi,8),2));
% $$$ axis('xy');
% $$$ xlabel('BHV Info PP');
% $$$ ylabel('Max Rate PP');
% $$$ 
% $$$ [phzVals,phzPref] = max(sq(-cvals(:,1,1,:)),[],2);
% $$$ 
% $$$ muind = uind(ismember(cluSessionMap(:,1),17:23));
% $$$ myUDP = unitDepths(muind);
% $$$ 
% $$$ hpm = sq(cmins(:,1,1,:));
% $$$ hpm = hpm(sub2ind([size(cmins,1),10],(1:size(cmins,1))',phzPref));
% $$$ hpm = hpm(uind);
% $$$ 
% $$$ fhpm = sq(cmins(:,1,1,:));
% $$$ fhpm = fhpm(sub2ind([size(cmins,1),10],(1:size(cmins,1))',phzPref-4+10*double((phzPref-4)<=0)));
% $$$ %fhpm = fhpm(sub2ind([size(cmins,1),10],(1:size(cmins,1))',phzPref+3-10*double((phzPref+3)>=11)));
% $$$ fhpm = fhpm(uind);
% $$$ 
% $$$ bpm = sq(cmins(:,1,2,:));
% $$$ bpm = bpm(sub2ind([size(cmins,1),10],(1:size(cmins,1))',phzPref));
% $$$ bpm = bpm(uind);
% $$$ 
% $$$ 
% $$$ [bsi,bsipp] = max(sq(eSpi),[],2);
% $$$ tcvals = -cvals(:,1,1,:);
% $$$ bsippVals = tcvals(sub2ind([size(cmins,1),10],(1:size(cmins,1))',bsipp));
% $$$ bsippVals(uind);
% $$$ bsi = bsi(uind);
% $$$ bsipp = bsipp(uind);
% $$$ 
% $$$ 
% $$$ 
% $$$ phzPref = phzPref(uind);
% $$$ nind = nniz([phzPref,hpm,bpm])&phzVals(uind)>1&bpm<=6;%(bpm==3|bpm==4);%
% $$$ %nind = nniz([phzPref,hpm,bpm])&phzVals(uind)>1.5&bpm>6;%(bpm==3|bpm==4);
% $$$ out = accumarray([phzPref(nind),hpm(nind)],ones([sum(nind),1]),[10,13],@sum);
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,
% $$$ for ses = 1:23
% $$$ muind = ismember(cluSessionMap(uind,1),ses)&nind;
% $$$ subplot(5,5,ses);
% $$$ hist(phzPref(muind)-5+10*double(phzPref(muind)<=5),1:10)
% $$$ title(Trials{ses}.name)
% $$$ end
% $$$ 
% $$$ figure,
% $$$ for ses = 1:23
% $$$ muind = ismember(cluSessionMap(uind,1),ses)&nind;
% $$$ subplot(5,5,ses);
% $$$ hist(bsipp(muind)-5+10*double(bsipp(muind)<=5),1:10)
% $$$ title(Trials{ses}.name)
% $$$ end
% $$$ 
% $$$ figure,
% $$$ subplot(221);barh((pfs.adata.bins{2}+pi)/pi*180,histc(phzPref(nind)-5+10*double((phzPref(nind)<=5)),1:10),'histc');
% $$$ ylabel('theta phase');
% $$$ subplot(222);imagesc(pfs.adata.bins{1},(phzBins+pi)/pi*180,out([6:10,1:5],:));axis('xy');
% $$$ subplot(224);bar(pfs.adata.bins{1},histc(hpm(nind),1:13),'histc');
% $$$ xlabel('head pitch');
% $$$ xlim(gca,pfs.adata.bins{1}([1,end]));
% $$$ suptitle('Preferred theta and head pitch for units of subject jg05');
% $$$ %caxis([0,10])
% $$$ %u,m,1,p
% $$$ 
% $$$ figure,hist2([cmins(uind,1,1,3),cmins(uind,1,1,8)]+0.5,1:14,1:14);line([1,14],[1,14]);  
% $$$ 
% $$$ figure,
% $$$ hist2([bsipp(nind)-5+10*double((bsipp(nind)<=5)),phzPref(nind)-5+10*double((phzPref(nind)<=5))],10,10);
% $$$ 
% $$$ figure();
% $$$ subplot(151);ind = nind&ismember(phzPref,[8,9,10]);hist2([hpm(ind),fhpm(ind)]+0.5,1:14,1:14);line([1,14],[1,14]);  
% $$$ subplot(152);ind = nind&ismember(phzPref,[6,7,8]);hist2([hpm(ind),fhpm(ind)]+0.5,1:14,1:14);line([1,14],[1,14]);
% $$$ subplot(153);ind = nind&ismember(phzPref,[4,5,6]);hist2([hpm(ind),fhpm(ind)]+0.5,1:14,1:14);line([1,14],[1,14]);
% $$$ subplot(154);ind = nind&ismember(phzPref,[2,3,4]);hist2([hpm(ind),fhpm(ind)]+0.5,1:14,1:14);line([1,14],[1,14]);
% $$$ subplot(155);ind = nind&ismember(phzPref,[10,1,2]);hist2([hpm(ind),fhpm(ind)]+0.5,1:14,1:14);line([1,14],[1,14]);
% $$$ 
% $$$ figure,
% $$$ scatter(bsippVals(nind),phzVals(nind),10,bsipp(nind)-5+10*double((bsipp(nind)<=5)),'filled')
% $$$ colormap(gca,'hsv')
% $$$ 
% $$$ figure,
% $$$ scatter(phzPref(nind)-5+10*double((phzPref(nind)<=5)),phzVals(nind),10,myUDP(nind),'filled')
% $$$ colormap('jet');
% $$$ 
% $$$ figure();
% $$$ plot(phzPref(nind)-5+10*double((phzPref(nind)<=5))+randn([sum(nind),1])/8,myUDP(nind)+randn([sum(nind),1])/8,'.');
% $$$ figure();
% $$$ hist2([phzPref(nind)-5+10*double((phzPref(nind)<=5)),myUDP(nind)],10,15);
% $$$ 
% $$$ pcc = [phzPref(nind)-5+10*double((phzPref(nind)<=5)), ...
% $$$        myUDP(nind)];
% $$$ [R,P] = corrcoef(pcc(nniz(pcc),:))
% $$$ 
% $$$ figure();
% $$$ subplot(121);
% $$$ hist2([phzPref(nind)-5+10*double((phzPref(nind)<=5)),myUDP(nind)],10,10);
% $$$ ylabel('depth');
% $$$ xlabel('preferred theta phase');
% $$$ subplot(122);
% $$$ hist2([hpm(nind)+0.5,myUDP(nind)],1:14,10);
% $$$ ylabel('depth');
% $$$ xlabel('head pitch');
% $$$ 
% $$$ pcc = [hpm(nind&bsi>1),myUDP(nind&bsi>1)];
% $$$ [R,P] = corrcoef(pcc(nniz(pcc),:))


%subplot(142);ind = nind&phzPref<=5;hist2([hpm(ind),fhpm(ind)]+0.5,1:14,1:14);line([1,14],[1,14]);






% $$$ 
% $$$             figure,
% $$$             subplot(1,11,1);                            
% $$$             plot(pft,u);
% $$$             for p = 1:10,
% $$$                 subplot(2,11,p+1);hold('on');imagescnan({pfd{1}.adata.bins{:},szmap(:,:,p)'./mean(nonzeros(reshape(szmap(:,:,p),[],1)))},[0,max(szmap(:))],[],true);axis('xy');axis('tight');
% $$$                 subplot(2,11,p+12);hold('on');imagescnan({pfd{1}.adata.bins{:},zmap(:,:,p)'},[0,max(szmap(:))],[],true);axis('xy');axis('tight');
% $$$                 for m = 1:sum(all(cmins(:,:,p),2)),
% $$$                     plot(pfd{1}.adata.bins{1}(cmins(m,1,p)),pfd{1}.adata.bins{2}(cmins(m,2,p)),'*m');
% $$$                 end
% $$$             end
% $$$  
% $$$             
% $$$             
% $$$             for i= 1:10,
% $$$                 for j = 1:10,
% $$$                 iva(i,j) = sum(reshape(zmap(:,:,i),[],1).*reshape(zmap(:,:,mod(j,10)+1),[],1))./ ...
% $$$                             (sqrt(sum(reshape(zmap(:,:,i),[],1).^2)).*sqrt(sum(reshape(zmap(:,:,mod(j,10)+1),[],1).^2)));
% $$$                 end
% $$$             end
% $$$ 
% $$$             
% $$$             figure,
% $$$             subplot(131);imagesc(iva)
% $$$             subplot(132);imagesc(diff(iva,1,1))
% $$$             subplot(133);imagesc(diff(iva,2,1))            
% $$$             
% $$$ 
% $$$             
% $$$             figure,
% $$$             subplot(141);imagesc(srmap'),axis('xy')
% $$$             subplot(142);imagesc(Smoother'),axis('xy')
% $$$             subplot(143);imagesc(zmap(:,:,5)'),axis('xy')
% $$$             subplot(144);imagesc(zmap(:,:,9)'),axis('xy')
% $$$ 
% $$$             
% $$$ % $$$             for v = 1:5
% $$$ % $$$                 eScra(ucnt,:,v) = sum(bsxfun(@times,repmat(eigVec(:,v),[1,10]),zmap),'omitnan');
% $$$ % $$$                 eScrn(ucnt,:,v) = sum(bsxfun(@times,...
% $$$ % $$$                                             bsxfun(@rdivide,zmap,mean(zmap,'omitnan')),... normalize map
% $$$ % $$$                                             repmat(eigVec(:,v),[1,10])),'omitnan');
% $$$ % $$$                 eReg(ucnt,:,v) = sum(bsxfun(@minus,...
% $$$ % $$$                                             bsxfun(@times,max(zmap),...
% $$$ % $$$                                                    repmat(eigVec(:,v)./max(eigVec(:,v)),[1,10])),...
% $$$ % $$$                                             zmap).^2,'omitnan');
% $$$ % $$$             end            
% $$$ 
% $$$             eVar(ucnt,:) = bsxfun(@rdivide,...
% $$$                                     sum(bsxfun(@rdivide,zmap,mean(zmap,'omitnan')),'omitnan'),...
% $$$                                     sum(~isnan(zmap)));
% $$$             eSpi(ucnt,:) = bsxfun(@rdivide,sum(bsxfun(@rdivide,zmap,mean(zmap,'omitnan')).*log2(bsxfun(@rdivide,zmap,mean(zmap,'omitnan'))),'omitnan'),sum(~isnan(zmap)));
% $$$ 
% $$$             
% $$$             zmap(isnan(zmap(:,1))|sum(zmap','omitnan')'==0,:) = [];
% $$$             
% $$$             [LU,LR,FSr,VT] = erpPCA(zmap',6);
% $$$             for v = 1:6,
% $$$                 evec = nan([size(rmap,1),size(rmap,3)]);
% $$$                 evec(rmapNind) = LR(:,v);
% $$$                 eigVecs(ucnt,v,:,:) = evec;
% $$$             end
% $$$ 
% $$$ 
% $$$             figure,
% $$$             subplot(2,11,1);                            
% $$$             plot(pft,u);
% $$$             for p = 1:10,
% $$$                 %zmapNormalizedPart = zmap(:,p)./sum(zmap(:,p));                
% $$$                 %mBhvPos = [sum(hp.*zmapNormalizedPart),sum(bp.*zmapNormalizedPart)];
% $$$                 subplot(2,11,p+1);hold('on');
% $$$                 imagesc(pfd{1}.adata.bins{:},sq(rmap(:,p,:))');
% $$$                 axis('xy');
% $$$                 axis('tight');
% $$$                 %plot(mBhvPos(1),mBhvPos(2),'*m');                
% $$$             end
% $$$             subplot(2,11,12),
% $$$             plot(VT(:,4),'-+');
% $$$             for v = 1:6,                     
% $$$                 subplot(2,11,11+1+v);hold('on');                
% $$$                 imagesc(pfs.adata.bins{[1,3]},sq(eigVecs(ucnt,v,:,:))');
% $$$             end
% $$$             
% $$$             eigScrs(ucnt,:,:) = FSr(:,1:2);
% $$$             eigVals(ucnt,:) = VT(1:2,4);
% $$$         end
% $$$         sessionCluMap(ucnt,:) = [t,Trial.spk.map(u,:)];
% $$$         ucnt = ucnt+1;
% $$$     end
% $$$ end

%eReg = sum(eReg(:,:,[1:3]),3);
%eReg = sum(eReg(:,:,[1]),3)
% $$$ 
% $$$ 
% $$$ hfig = figure();
% $$$ unitInds = 1:numel(unitSubset)
% $$$ u = 1;
% $$$ 
% $$$ %while u ~=-1
% $$$ for u = unitInds,
% $$$     clf();
% $$$     for i= 1:4
% $$$         subplot2(10,10,1:2,i*2-1:i*2);    
% $$$         imagescnan({pfs.adata.bins{[1,3]},reshape(eigVec(:,i),[pfd{1}.adata.binSizes'])'}),axis xy;
% $$$     end
% $$$     subplot2(10,10,1:2,5*2-1:5*2);        
% $$$     hold('on');
% $$$     plot3(reshape(eScrn(u,:,2),[],1)./reshape(eScrn(u,:,3),[],1),...
% $$$           reshape(eScrn(u,:,1),[],1)./reshape(eScrn(u,:,3),[],1),...
% $$$           eScrn(u,:,3));
% $$$     scatter3(reshape(eScrn(u,:,2),[],1)./reshape(eScrn(u,:,3),[],1),...
% $$$              reshape(eScrn(u,:,1),[],1)./reshape(eScrn(u,:,3),[],1),...
% $$$              eScrn(u,:,3),...
% $$$              20,...
% $$$              reshape(repmat(permute(hsv(10),[1,3,2]),[1,1,1]),[],3),...
% $$$              'filled');
% $$$     daspect([0.01,0.01,1]);    
% $$$     subplot2(10,10,3:8,1:5);
% $$$     plot(pfd{1},unitSubset(u),'mean','colorbar',[],false);
% $$$     subplot2(10,10,3:8,6:10);    
% $$$     scatter3(reshape(eScrn(u,:,2),[],1),...
% $$$              reshape(eScrn(u,:,3),[],1),...
% $$$              reshape(eScrn(u,:,1),[],1),...
% $$$              30,...
% $$$              reshape(repmat(permute(hsv(10),[1,3,2]),[1,1,1]),[],3),...
% $$$              'filled'...
% $$$              );
% $$$     hold('on');
% $$$     plot3(reshape(eScrn(u,:,2),[],1),...
% $$$           reshape(eScrn(u,:,3),[],1),...
% $$$           reshape(eScrn(u,:,1),[],1));    
% $$$     daspect([1,1,1]);
% $$$     colorbar();
% $$$     caxis([-pi,pi]);
% $$$     colormap(gca(),'hsv');
% $$$     rmap = plot(pfs,unitSubset(u),1,'colorbar',[],false);    
% $$$     rmax = prctile(rmap(nniz(rmap(:))),98)+2;
% $$$     for i = 1:10,
% $$$         subplot2(10,10,[9,10],i);
% $$$         imagescnan({pfs.adata.bins{[1,3]},sq(rmap(:,mod(i+4,10)+1,:))'},[0,rmax],'colorMap',@jet);
% $$$         axis('xy');
% $$$         title(num2str(pfs.adata.bins{2}(mod(i+4,10)+1)));
% $$$     end
% $$$     
% $$$     waitforbuttonpress();
% $$$     %u = figure_controls(hfig,u,unitInds,false,[],[]);
% $$$ end
% $$$ 
% $$$   
% $$$ figure
% $$$     plot3(log10(eReg(:,7,1)),...
% $$$           log10(eReg(:,7,3)),...
% $$$           log10(eReg(:,7,2)),'.');
% $$$ 
% $$$ figure,imagesc(eigScrs(:,:,1))
% $$$ 
% $$$ [~,pc1m] = max(abs(eigScrs(:,:,1)),[],2);
% $$$ [~,pc2m] = max(abs(eigScrs(:,:,2)),[],2);
% $$$ 
% $$$ figure,hist2(eigVals,30,30)
% $$$ figure,hist(diff(eigVals,1,2),100)
% $$$ phzBins = linspace(-pi,pi,10)';
% $$$ 
% $$$ ind = diff(eigVals,1,2)>-20;
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,plot(phzBins(pc1m(ind))+randn([sum(ind),1])/3,phzBins(pc2m(ind))+randn([sum(ind),1])/3,'.')
% $$$ 
% $$$ figure,rose(phzBins(pc1m(ind))+randn([sum(ind),1])/3-phzBins(pc2m(ind))+randn([sum(ind),1])/3)
% $$$ 
% $$$ figure,
% $$$ hold('on');plot(phzBins(pc1m(ind))+randn([sum(ind),1])/3,phzBins(pc1m(ind))+randn([sum(ind),1])/3- ...
% $$$                 phzBins(pc2m(ind))+randn([sum(ind),1])/3,'.b')
% $$$ 
% $$$ 
% $$$ 
% $$$ figure();hold('on');
% $$$ plot3(reshape(eScrn(:,:,2),[],1),...
% $$$       reshape(eScrn(:,:,3),[],1),...
% $$$       reshape(eScrn(:,:,1),[],1),'.b');
% $$$ 
% $$$ 
% $$$ figure();hold('on');
% $$$ scatter3(reshape(eScrn(:,:,2),[],1),...
% $$$       reshape(eScrn(:,:,3),[],1),...
% $$$       reshape(eScrn(:,:,1),[],1),30,reshape(eSpi,[],1),'filled');
% $$$ 
% $$$ ind = 20:40
% $$$ figure();hold('on');
% $$$ plot3(reshape(eScra(ind,[2:2:8,2],2),[],5)',...
% $$$       reshape(eScra(ind,[2:2:8,2],3),[],5)',...
% $$$       reshape(eScra(ind,[2:2:8,2],1),[],5)');
% $$$ 
% $$$ 
% $$$ ind = 20:40
% $$$ figure();hold('on');
% $$$ plot3(reshape(eScrn(ind,[2:2:8,2],2),[],5)',...
% $$$       reshape(eScrn(ind,[2:2:8,2],3),[],5)',...
% $$$       reshape(eScrn(ind,[2:2:8,2],1),[],5)');
% $$$ 
% $$$ plot3(reshape(eScr(:,7,2),[],1),...
% $$$       reshape(eScr(:,7,3),[],1),...
% $$$       reshape(eScr(:,7,1),[],1),'.r');
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ for i = 1:10
% $$$ plot3(eScr(:,i,2)',eScr(:,i,3)',eScr(:,i,1)','.');
% $$$ waitforbuttonpress();
% $$$ end
% $$$ 
% $$$ ind = mean(eScr(:,7,1),2)<18;
% $$$ figure,plot(mean(eVar(ind,:,1))),
% $$$ hold('on');plot(mean(eVar(ind,:,1))+std(eVar(ind,:,1)),'r'),
% $$$ plot(mean(eVar(ind,:,1))-std(eVar(ind,:,1)),'r')
% $$$ 
% $$$ 
% $$$ %ind = mean(eScr(:,9,1),2)<18;
% $$$ ind = ':';
% $$$ %ind = any(eSpi>0.5,2);
% $$$ mes = mean(eSpi(ind,:));
% $$$ ses = std(eSpi(ind,:));
% $$$ 
% $$$ figure,plot(phzBins,mes);
% $$$ hold('on');
% $$$ plot(phzBins,mes+ses,'r'),
% $$$ plot(phzBins,mes-ses,'r')
% $$$ 
% $$$ ind = ':';
% $$$ mdes = mean(diff(eSpi(ind,:,1),1,2));
% $$$ sdes = std(diff(eSpi(ind,:,1),1,2));
% $$$ phzBinsd = phzBins(1:end-1)+0.5*diff(phzBins(1:2));
% $$$ figure();
% $$$ plot(phzBinsd,mdes)
% $$$ hold('on');
% $$$ plot(phzBinsd,mdes+sdes,'r');
% $$$ plot(phzBinsd,mdes-sdes,'r');
% $$$ Lines(phzBinsd([1,end]),0,'k');
% $$$ 
% $$$ 
% $$$ co = [];
% $$$ sc = [];
% $$$ la = [];
% $$$ for u= 1:size(eScrn,1),
% $$$     [co(u,:,:),sc(u,:,:),la(u,:)] = pca(sq(eScrn(u,:,1:3)));
% $$$ end
% $$$ 
% $$$ figure,plot(log10(la(:,1)./la(:,2)),eSpi(:,7),'.')
% $$$ 
% $$$ 
% $$$ sng = (co(:,2,1)-polyval([0,0.5],co(:,1,1)))<0;
% $$$ mrr = [-1,0,0;0,-1,0,;0,0,-1];
% $$$ sco = co(:,:,1);
% $$$ for u = find(sng),
% $$$     sco(u,:) = sco(u,:)*mrr;
% $$$ end
% $$$ 
% $$$ figure();
% $$$ scatter3(sco(:,1,1),sco(:,2,1),sco(:,3,1),20,log10(la(:,1)),'filled');


% $$$                 %[mins,vals] = LocalMinimaN(convn(-zmap(:,:,p),Smoother,'same'),-2,7);
% $$$                 [mins,vals] = LocalMinimaN(-zmap(:,:,p)+mrate,-0.5,5);
% $$$                 vals =  vals-mrate;
% $$$                 cmins(ucnt,1:size(mins,1),:,p) = mins;
% $$$                 cvals(ucnt,1:size(mins,1),1,p) = vals;
% $$$                 for c = 1:size(mins,1)
% $$$                     weights = weightsFunction(gridBins,...
% $$$                                               gridBins(bsxfun(@plus,sub2ind(cellfun(@numel,binc),...
% $$$                                                                     cmins(ucnt,c,1,p),...
% $$$                                                                     cmins(ucnt,c,2,p)),[0,size(gridBins,1)])));
% $$$                     weights = weights./sum(weights(:),'omitnan');                                        
% $$$                     weigths = weights.*reshape(zmap(:,:,p),[],1)./max(reshape(zmap(:,:,p),[],1));
% $$$                     weights(~rmapNind) = 0;
% $$$                     weights = weights./sum(weights(:),'omitnan');                    
% $$$                     smins(ucnt,c,:,p) = sum(bsxfun(@times,weights,gridBins));
% $$$                 end
