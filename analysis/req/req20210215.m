

MjgER2016_load_data();
%  Variables:
%      Trials
%      units
%      cluSessionMap
%      pitchReferenceTrial
%      figBasePath
%      sessionListName
%      sessionList
%      stcMode
%      states
%      numStates
%      interpParPfsp
%      interpParDfs
%
%  Functions:
%      reshape_eigen_vector

MjgER2016_load_interneurons();
%  Variables:
%      unitsInts

configure_default_args();
%  GENERAL 
%      sessionListName
%      activeState 
%      fet_HB_pitchB.referenceTrial = 'Ed05-20140529.ont.all';
%  FUNCTIONS 
%      compute_bhv_ratemaps
%      compute_bhv_ratemaps_shuffled
%      compute_ratemaps
%      compute_bhv_ratemap_erpPCA
%      compute_bhv_ratemap_nnmf
%      compute_xyhb_ratemaps


% LOAD place restricted behavior fields
bfs   = cf(@(t,u)   compute_bhv_ratemaps(t,u),  Trials, units);
%pfsa = cf(@(s,t) cat(2,{t},s), pfss,pfts);
% COMPUTE bfs erpPCA
[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfs, units);

bhvMask = false(size(validDims));
bhvMask(validDims) = true;
bhvMask = reshape_eigen_vector(bhvMask,bfs)';

% fet
t = 20;
clearvars('-GLOBAL','AP');

pfsArgs = struct('states',           'theta-groom-sit',      ...
                 'binDims',          [0.1,0.1,0.8],          ...
                 'SmoothingWeights', [1.8,1.8,0.5],                ...
                 'numIter',          1,                      ...
                 'boundaryLimits',   [-2,0.8;-0.8,2;-9,-4],   ...
                 'halfsample',       false);

ifr = cf(@(T,U)  compute_bhv_ratemaps(T,U,'fet_HPBPHR',[],'none',pfsArgs,1,inf,true), Trials(t),unitsInts(t));

rhm = fet_rhm(Trial,[],'mtchglong','newSR',16);

pfsArgs = struct('states',           'theta-groom-sit',      ...
                 'binDims',          [0.1,0.1,0.5],          ...
                 'SmoothingWeights', [1.8,1.8,0.5],                ...
                 'numIter',          1,                      ...
                 'boundaryLimits',   [-2,0.8;-0.8,2;-2,2],   ...
                 'halfsample',       false);

ifs = cf(@(T,U)  compute_bhv_ratemaps(T,U,@fet_HPBPHS,[],'none',pfsArgs,1,inf,false), Trials,unitsInts);

[irmaps,cluSessionMap] = decapsulate_and_concatenate_mtaapfs(ifs,unitsInts);

ibins = ifs{1}.adata.bins;
ibinDims = ifs{1}.adata.binSizes';

% RESHAPE rmaps into [ xPosition, yPosition, headPitch, bodyPitch, unit ]
irmapa = nan([ibinDims,size(irmaps,2)]);
for u = 1:size(irmaps,2),
    irmapa(:,:,:,u) = reshape(irmaps(:,u),ibinDims);
end
irmapa = irmapa.*repmat(maskBhv,[1,1,size(irmapa,3),size(irmapa,4)]);
irmapa(irmapa==0) = nan;
sbins =     -2:0.5:2;




cu = [27,25];
cu = [20,45];
cu = [20,92];
cu = [20,106];
cu = [4,15];
cu = [27,31];
figure,
u = find(ismember(cluSessionMap,cu,'rows'));%unitsInts{20} == 106;%find(ismember(cluSessionMap,[20,8],'rows'));
cmax = prctile(nonzeros(irmapa(:,:,:,u)),99);

subplot2(4,ibinDims(3),[2:4],[1,ibinDims(3)]);
hold('on');
for y = 1:28
    for x = 1:28
        plot(sq(irmapa(x,y,:,u)),'+-');
    end
end
xlim([0.5,8.5]);

for s = 1:ibinDims(3)
    subplot2(4,ibinDims(3),1,s);
    imagesc(ibins{1},ibins{2},irmapa(:,:,s,u)');%.*bhvMask')
    axis('xy');
    colormap('jet');
    caxis([0,cmax]);
    title(num2str(10.^sbins([s,s+1]),'%3.2f - %3.2f'));
    Lines([-1.75,0.5],0.65,'m');
    Lines(-0.8,[-0.5,0.65],'m');
    
end
cax = colorbar();
cax.Position(1) = cax.Position(1) +0.05;



figure,        
hold('on');
plot(mean(reshape(irmapa(:,ibins{2}>0.65,:,u),[],8),'omitnan'),'+-r');
plot(mean(reshape(irmapa(ibins{1}<-0.8,ibins{2}<0.65,:,u),[],8),'omitnan'),'+-c');
plot(mean(reshape(irmapa(ibins{1}>-0.8,ibins{2}<0.65,:,u),[],8),'omitnan'),'+-g');

for u = 1:size(irmapa,4)
P(u,:,1) = polyfit(ibins{3}(2:end-1),mean(reshape(irmapa(ibins{1}<-0.8,ibins{2}<0.65,2:end-1,u),[],6),'omitnan')',1);
P(u,:,2) = polyfit(ibins{3}(2:end-1),mean(reshape(irmapa(ibins{1}>-0.8,ibins{2}<0.65,2:end-1,u),[],6),'omitnan')',1);
P(u,:,3) = polyfit(ibins{3}(2:end-1),mean(reshape(irmapa(:,ibins{2}>0.65,2:end-1,u),[],6),'omitnan')',1);
end


figure,
subplot(131);plot(P(:,1,1),P(:,2,1),'.');title('low')
subplot(132);plot(P(:,1,2),P(:,2,2),'.');title('high')
subplot(133);plot(P(:,1,3),P(:,2,3),'.');title('rear')


figure();
plot(P(:,2,2)-P(:,2,1),P(:,1,2),'.');

figure();
plot(P(:,1,1),P(:,1,3),'.');
grid('on');
line([-4,10],[-4,10]);
xlabel('slope (low)')
ylabel('slope (high)')


figure();
ind = all(P(:,1,[1,2])<0,3);
plot(P(ind,1,3),P(ind,2,3),'.')



figure();
hist2([P(:,1,1),P(:,2,1)],linspace(-5,10,20),linspace(0,50,20));

ind = P(:,2,1)>10;
figure,
subplot(121)
plot(P(ind,1,1),P(ind,1,2),'.')
subplot(122)
plot(P(ind,1,1),P(ind,1,3),'.')