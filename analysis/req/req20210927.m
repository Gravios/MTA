;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.


% center1 center2
% area around center1 as state


cf(@(T) T.load('nq'), Trials);

% bhv field parameters
bfs      = cf(@(t,u)  compute_bhv_ratemaps(t,u),          Trials, units);
[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfs, units, [], [], false);
numComp = size(eigVecs,2);



overwrite = false;
threshDist = 150;
sampleRate = 16;
pfsArgs = struct('states',           'theta-groom-sit',  ...
                 'binDims',          [0.1,0.1],          ...
                 'SmoothingWeights', [1.8,1.8],          ...
                 'numIter',          200,                ...                
                 'boundaryLimits',   [-2,0.8;-0.8,2],    ...
                 'mask',             validDims,          ...
                 'halfsample',       false);

unit = 34;
%% compute_permuted_patch_ratemap
%uids = find(ismember(cluSessionMap(:,1),tid));

rmapA = zeros([784,numUnits,6,6]);
rmapB = zeros([784,numUnits,6,6]);
rmapZscr = zeros([784,numUnits,6,6]);
rmapCorr = zeros([200,numUnits,6,6]);
for uid = 1:numUnits
    if tid ~= cluSessionMap(uid,1) | uid==1
        tid = cluSessionMap(uid,1);
        Trial = Trials{tid};
        xyz = preproc_xyz(Trial,'trb',sampleRate);
        fet = fet_HB_pitchB(Trial,sampleRate);
        spk  = create(copy(Trial.spk),Trial,sampleRate,pfsArgs.states,units{tid},'');
    end
        
    unit = cluSessionMap(uid,2);
    numPatches = sum(~isnan(patchCntrF(uid,:,1)));
    tic
    if numPatches > 1
        for p1 = 1:numPatches-1
            for p2 = p1+1:numPatches
                tag = ['pfsPerm_',num2str([unit,p1,p2],'%d_%d-%d')];
% $$$                 [rmapDiff,ratemapA,ratemapB] = ...
                [rmapZscr(:,uid,p1,p2),                                         ...
                 rmapA(:,uid,p1,p2),rmapB(:,uid,p1,p2),                         ...
                 rmapCorr(:,uid,p1,p2)] =                                       ...
                    compute_permuted_patch_ratemap(Trial,                       ...
                                                   unit,                        ...
                                                   [],                          ... fetset
                                                   sampleRate,                  ...
                                                   patchCntrF(uid,[p1,p2],:),   ...
                                                   'hcom',                      ... marker
                                                   pfsArgs,                     ... 
                                                   [],                          ... threshRate
                                                   threshDist,                  ...
                                                   xyz,                         ...
                                                   fet,                         ...
                                                   spk,                         ...
                                                   tag,                         ...
                                                   overwrite);

            end
        end
    end
    toc
end


%%%<<< COMPUTE the patch bhv rate maps -------------------------------------------------------------
rmapP = zeros([784,numUnits,6]);
for uid = 1:numUnits
    if tid ~= cluSessionMap(uid,1)
        tid = cluSessionMap(uid,1);
        Trial = Trials{tid};
        xyz = preproc_xyz(Trial,'trb',sampleRate);
        fet = fet_HB_pitchB(Trial,sampleRate);
        spk  = create(copy(Trial.spk),Trial,sampleRate,pfsArgs.states,units{tid},'');
    end
    unit = cluSessionMap(uid,2);
    numPatches = sum(~isnan(patchCntrF(uid,:,1)));
    tic
    for p = 1:numPatches
        [rmapP(:,uid,p)] =                                         ...
                compute_patch_ratemap(Trial,                       ...
                                      unit,                        ...
                                      [],                          ... fetset
                                      sampleRate,                  ...
                                      patchCntrF(uid,p,:),         ...
                                      'hcom',                      ... marker
                                      pfsArgs,                     ... 
                                      [],                          ... threshRate
                                      threshDist,                  ...
                                      xyz,                         ...
                                      fet,                         ...
                                      spk,                         ...
                                      tag,                         ...
                                      overwrite);
    end
    toc
end
%%%>>>


figure,
for p = 1:6,
subplot(1,6,p);
imagesc(reshape(rmapP(:,ismember(cluSessionMap,[20,31],'rows'),p).*validDims,[28,28])');axis('xy');colormap('jet');colorbar();
end



stateLabels = {'Rear','High','Low'};
tid = 20;
u = 34;
p1 = 1; p2 = 2;
%for u = units{tid};
    uid = find(ismember(cluSessionMap,[tid,u],'rows'));
    mrate = max(apfstats.peakFR(uid,:));
    clf();
    for s = 1:3
        subplot2(2,3,1,s);
        plot(pfsr{tid}{s},u,1,'colorbar',[0,mrate],'colorMap',@parula);
        hold('on');
        for p = 1:2
            if apfstats.patchPFR(uid,s,p)>rthresh
                rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
                rinds = rinds(nniz(rinds),:);
                if size(rinds,1)>athresh/400
                    plot(pfsr{tid}{s}.adata.bins{1}(rinds(:,1)),...
                         pfsr{tid}{s}.adata.bins{2}(rinds(:,2)),...
                         'm.');                    
                end        
            end
        end
        for p = 1:sum(nniz(patchCntrF(uid,:,1)'))
            circle(patchCntrF(uid,p,1),patchCntrF(uid,p,2),150);
            text(patchCntrF(uid,p,1),patchCntrF(uid,p,2),num2str(p));
        end
        title({['Unit: ',num2str(u)],[stateLabels{s}, ' State Ratemap']});
    end
    subplot2(2,3,2,1);
        imagesc(bins{:},reshape(rmapA(:,uid,p1,p2).*validDims,cellfun(@numel,bins))');
        title(['Patch ',num2str(p1),' bhv ratemap']);
        axis('xy');
        cax = colorbar();
        ylabel(cax,'Hz');                
        colormap(gca(),'jet');
        xlabel('Head-Body pitch (rad)');
        ylabel('Body pitch (rad)');    
    subplot2(2,3,2,2);
        imagesc(bins{:},reshape(rmapB(:,uid,p1,p2).*validDims,cellfun(@numel,bins))');
        title(['Patch ',num2str(p2),' bhv ratemap']);
        xlabel('Head-Body pitch (rad)');
        ylabel('Body pitch (rad)');
        axis('xy');
        cax = colorbar();
        ylabel(cax,'Hz');                
        colormap(gca(),'jet');
    subplot2(2,3,2,3);
        imagesc(bins{:},reshape(rmapZscr(:,uid,p1,p2).*validDims,cellfun(@numel,bins))');
        axis('xy');
        cax = colorbar();
        ylabel(cax,'z-score');
        colormap(gca(),'jet');    
        title(num2str([p1,p2],'%d <-> %d'));
        xlabel('Head-Body pitch (rad)');
        ylabel('Body pitch (rad)');
        title('Permuted Patch Ratemap Z-Score');
        
    
% $$$     waitforbuttonpress();
% $$$ end



%%%<<< COMPUTE inter placecell bhv ratemap correlations
% only compare units on separate electrodes
rmapIPCorr = [];
rmapIPDist = [];
rmapIPBcnt = [];
for tid = numel(Trials),
    for u1 = units{tid}
        uid1 = find(ismember(cluSessionMap,[tid,u1],'rows'));
        for u2 = units{tid}
            if u1 == u2 || (Trials{tid}.nq.ElNum(u1)==Trials{tid}.nq.ElNum(u2) && abs(Trials{tid}.nq.maxAmpChan(u1)-Trials{tid}.nq.maxAmpChan(u2))<=2)
                continue
            end
            uid2 = find(ismember(cluSessionMap,[tid,u2],'rows'));
            for gptch1 = find(~isnan(patchCntrF(uid1,:,2)))
                for gptch2 = find(~isnan(patchCntrF(uid2,:,2)))
                    rmapPC = [rmapP(validDims,uid1,gptch1),rmapP(validDims,uid2,gptch2)];
                    if any(sum(~isnan(rmapPC)) > 200)
                        continue
                    end
                    c = corr(rmapPC(~isnan(rmapPC(:,1))&~isnan(rmapPC(:,2)),:),'type','Spearman');
                    rmapIPCorr(end+1) = c(2);
                    rmapIPDist(end+1) = sqrt(sum((patchCntrF(uid1,gptch1,:)-patchCntrF(uid2,gptch2,:)).^2,3));
                    rmapIPBcnt(end+1) = sum(double(~isnan(rmapPC(:,1))&~isnan(rmapPC(:,2))));
                end
            end
        end
    end
end

rmapPCA = nan([784,1]);
rmapPCB = nan([784,1]);
rmapPCA(validDims) = rmapPC(:,1);
rmapPCB(validDims) = rmapPC(:,2);
figure,
subplot(121);
pcolor(reshape(rmapPCA,[28,28])');axis('xy');
subplot(122);
pcolor(reshape(rmapPCB,[28,28])');axis('xy');

figure,
uid = find(ismember(cluSessionMap,[20,20],'rows'));
imagesc(reshape(rmapP(:,uid,1),[28,28])');
axis('xy');

figure,
plot(rmapIPDist/10,rmapIPCorr,'.')
xlabel('Inter Patch Distance (cm)')
ylabel('Bhv Ratemap Correlation');

patchCntF = sum(~isnan(patchCntrF(:,:,1)),2);
out = histcounts(patchCntF,0.5:5.5);
figure,pie(out)
%%%>>>





        
stateLabels = {'Rear','High','Low'};
tid = 20;
u = 25;
p1 = 1; p2 = 2;
%for u = units{tid};
    uid = find(ismember(cluSessionMap,[tid,u],'rows'));
    mrate = max(apfstats.peakFR(uid,:));
    clf();
    for s = 1:3
        subplot2(2,3,1,s);
        plot(pfsr{tid}{s},u,1,'colorbar',[0,mrate],'colorMap',@parula);
        hold('on');
        for p = 1:2
            if apfstats.patchPFR(uid,s,p)>rthresh
                rinds = sq(apfstats.patchRateInd(uid,s,p,:,:))';
                rinds = rinds(nniz(rinds),:);
                if size(rinds,1)>athresh/400
                    plot(pfsr{tid}{s}.adata.bins{1}(rinds(:,1)),...
                         pfsr{tid}{s}.adata.bins{2}(rinds(:,2)),...
                         'm.');                    
                end        
            end
        end
        for p = 1:sum(nniz(patchCntrF(uid,:,1)'))
            circle(patchCntrF(uid,p,1),patchCntrF(uid,p,2),150);
            text(patchCntrF(uid,p,1),patchCntrF(uid,p,2),num2str(p));
        end
        title({['Unit: ',num2str(u)],[stateLabels{s}, ' State Ratemap']});
    end
        
set(gcf(),'PaperOrientation','landscape');
set(gcf(),'PaperType','A4');
        
nrA = ratemapA(validDims,uid,p1,p2);
nrB = ratemapB(validDims,uid,p1,p2);
nind = nniz(nrA) & nniz(nrB);
nrA = nrA(nind);
nrB = nrB(nind);

corr(nrA,nrB,'type','Spearman')

uid = find(ismember(cluSessionMap,[20,44],'rows'));
patchCenter = patchCntrF(uid,[1,2],:);

% REPORT intra cell patch distance vs correlation 
figure
hold('on');
ph = patch([10.1,10.1,30.0,30.0],[-1,1,1,-1],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.3);
pd = reshape(patchDistF(unitSubset,:,:),[],1)/10;
rcorr = reshape(rmapCorr(1,unitSubset,:,:),[],1);
nind = nniz(pd) & nniz(rcorr) & pd>10;
plot(pd(nind),rcorr(nind),'.');
box('on');
xlabel('Patch Distance (cm)')
ylabel('Bhv Correlation');





polyfit(pd(nind),rcorr(nind),'.');

ind = unitSubset;
ind = ':';
figure,
subplot(131);
plot(reshape(patchPFRF(ind,2,:),[],1),reshape(patchPFRF(ind,1,:),[],1),'.');xlim([0,25]);ylim([0,25]);
line([0,25],[0,25],'Color','k');
line([0,12.5],[0,25],'Color','r');
line([0,25],[0,12.5],'Color','r');
xlabel('Peak Patch Rate Hz (high)');
ylabel('Peak Patch Rate Hz (rear)');
subplot(132);
plot(reshape(patchMFRF(ind,3,:),[],1),reshape(patchMFRF(ind,1,:),[],1),'.');xlim([0,25]);ylim([0,25]);
line([0,25],[0,25],'Color','k');
line([0,12.5],[0,25],'Color','r');
line([0,25],[0,12.5],'Color','r');
xlabel('Peak Patch Rate Hz (low)');
ylabel('Peak Patch Rate Hz (rear)');
subplot(133);
plot(reshape(patchPFRF(ind,3,:),[],1),reshape(patchPFRF(ind,2,:),[],1),'.');xlim([0,25]);ylim([0,25]);
line([0,25],[0,25],'Color','k');
line([0,12.5],[0,25],'Color','r');
line([0,25],[0,12.5],'Color','r');
xlabel('Peak Patch Rate Hz (low)');
ylabel('Peak Patch Rate Hz (high)');


figure();
subplot(131);
violin((reshape(patchPFRF(ind,1,:),[],1)-reshape(patchPFRF(ind,2,:),[],1))./(reshape(patchPFRF(ind,1,:),[],1)+reshape(patchPFRF(ind,2,:),[],1)));
subplot(132);
violin((reshape(patchPFRF(ind,1,:),[],1)-reshape(patchPFRF(ind,3,:),[],1))./(reshape(patchPFRF(ind,1,:),[],1)+reshape(patchPFRF(ind,3,:),[],1)));
subplot(133);
violin((reshape(patchPFRF(ind,3,:),[],1)-reshape(patchPFRF(ind,2,:),[],1))./(reshape(patchPFRF(ind,3,:),[],1)+reshape(patchPFRF(ind,2,:),[],1)));

figure
histogram(patchCNTF(ind),0.5:6.5);