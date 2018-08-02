


% LOAD session data
MjgER2016_load_data();

% LOAD behavioral data
[pfd ,tags ,eigVec, eigVar, eigScore, validDims,...
 unitSubsets, unitIntersection, zrmMean, zrmStd] = req20180123_ver5(Trials);
numComp = size(eigVec{1},2);
pfindex = 1;

% LOAD behavioral scores
MjgER2016_load_bhv_erpPCA_scores();

% GENERATE gaussian smoothing kernel with std of # bins <- smoothingWeight
smoothingWeights = [0.05,0.05];
sind = cell(1,ndims(bins));
[sind{:}] = ndgrid(bins{:});
for i = 1:ndims(bins),  sind{i} = sind{i}.^2/smoothingWeights(i)^2/2;  end
smootherNarrow = exp(sum(-cat(ndims(bins)+1,sind{:}),ndims(bins)+1));
smootherNarrow = smootherNarrow./sum(smootherNarrow(:));

smoothingWeights = [0.15,0.15];
sind = cell(1,ndims(bins));
[sind{:}] = ndgrid(bins{:});
for i = 1:ndims(bins),  sind{i} = sind{i}.^2/smoothingWeights(i)^2/2;  end
smootherWide = exp(sum(-cat(ndims(bins)+1,sind{:}),ndims(bins)+1));
smootherWide = smootherWide./sum(smootherWide(:));

shuffle_bmap = @(m) circshift(circshift(m,randi(size(m,1)),1),randi(size(m,2)),2);


% GET Pfd information
binSizes = pfd{1}.adata.binSizes';
bins =  pfd{1}.adata.bins;

% SELECT behavior field
VTR = [];
LRR = [];
numIter = 100;
for n = 1:numIter,
    bmapRandShift = zeros([size(rmaps)]);
    for u = 1:size(rmaps,2)
        bmap = zeros([binSizes]); 
        bmap(validDims{1}) = rmaps(:,u);
        bmapSmooth  = convn(bmap, smootherWide,'same');
        bmapSmooth  = convn(bmapSmooth, smootherWide,'same');
        bmapMerge = bmapSmooth;
        bmapMerge(validDims{1}) = rmaps(:,u);
        bmapMergeSubset = bmapMerge;
        bmapMergeSubset(validDims{1}&bmap(:)==0) = bmapSmooth(validDims{1}&bmap(:)==0);

        % FIND the rows and columns of matrx 
        ginds = cell([1,ndims(bins)]);
        ginds{1} = ~all(~bmap,2);
        ginds{2} = ~all(~bmap,1);

        bmapMergeSubset = bmapMerge(ginds{:});
        binSubset = cf(@(b,i) b(i), bins, ginds);

        bmapMergeSubsetShuffle = shuffle_bmap(bmapMergeSubset);               
        bmapMergeSubsetShuffleSmooth  = convn(bmapMergeSubsetShuffle, smootherNarrow,'same');
        bmapMergeShuffleSmoothFull = zeros([binSizes]); 
        bmapMergeShuffleSmoothFull(ginds{:}) = bmapMergeSubsetShuffleSmooth;
        bmapMergeShuffleSmoothFull(~validDims{1}) = 0;

% $$$         while sum(bmapMergeShuffleSmoothFull(:))<sum(bmap(:)).*0.9
% $$$             bmapMergeSubsetShuffle = shuffle_bmap(bmapMergeSubset);               
% $$$             bmapMergeSubsetShuffleSmooth  = convn(bmapMergeSubsetShuffle, smootherNarrow,'same');
% $$$             bmapMergeShuffleSmoothFull = zeros([binSizes]); 
% $$$             bmapMergeShuffleSmoothFull(ginds{:}) = bmapMergeSubsetShuffleSmooth;
% $$$             bmapMergeShuffleSmoothFull(~validDims{1}) = 0;
% $$$ % $$$             
% $$$ % $$$             bmapMergeShuffleSmoothFull(:,:,u) = zeros([binSizes]); 
% $$$ % $$$             bmapMergeShuffleSmoothFull(ginds{:},u) = bmapMergeSubsetShuffleSmooth(:,:,u);
% $$$ % $$$             bmapMergeShuffleSmoothFull(~validDims{1},:) = 0;
% $$$         end

        bmapRandShift(:,u) = bmapMergeShuffleSmoothFull(validDims{1});
    end
    
    [~,LR,FSr,VT] = erpPCA(bmapRandShift',numComp,200);            % COMPUTE covariance-based PCA with Varimax rotation
    VTR = cat(2,VTR,VT(1:50,4));
    LRR = cat(3,LRR,LR);
    
end



% SELECT behavior field
bmapOriginal = zeros(size(rmaps));    
for u = 1:size(rmaps,2)
    bmap = zeros([binSizes]); 
    bmap(validDims{1}) = rmaps(:,u);
    bmapSmooth  = convn(bmap, smootherWide,'same');
    bmapSmooth  = convn(bmapSmooth, smootherWide,'same');
    bmapMerge = bmapSmooth;
    bmapMerge(validDims{1}) = rmaps(:,u);
    bmapMergeSubset = bmapMerge;
    bmapMergeSubset(validDims{1}&bmap(:)==0) = bmapSmooth(validDims{1}&bmap(:)==0);
    ginds = cell([1,ndims(bins)]);
    ginds{1} = ~all(~bmap,2);
    ginds{2} = ~all(~bmap,1);
    bmapMergeSubset = bmapMerge(ginds{:});
    binSubset = cf(@(b,i) b(i), bins, ginds);
    bmapMergeSubsetShuffle = bmapMergeSubset;                       
    bmapMergeSubsetShuffleSmooth  = convn(bmapMergeSubsetShuffle, smootherNarrow,'same');
    bmapMergeShuffleSmoothFull = zeros([binSizes]); 
    bmapMergeShuffleSmoothFull(ginds{:}) = bmapMergeSubsetShuffleSmooth;
    bmapMergeShuffleSmoothFull(~validDims{1}) = 0;
    bmapOriginal(:,u) = bmapMergeShuffleSmoothFull(validDims{1});
end
VTO = [];
LRO = [];
for n = 1:numIter,
    randPermInd = randperm(size(bmapOriginal,2));
    [~,LR,FSr,VT] = erpPCA(bmapOriginal(:,randPermInd(1:400))',numComp);            % COMPUTE covariance-based PCA with Varimax rotation
    VTO = cat(2,VTO,VT(1:50,4));
    LRO = cat(3,LRO,LR);
end

% $$$ randPermInd = randperm(size(bmapOriginal,2));
% $$$ [U,S,V] = svd(rmaps(:,randPermInd(1:500))',0);            % COMPUTE covariance-based PCA with Varimax rotation

fpcs  = cell([1,10]);
for i = 1:10,
    fpcs{i} = nan([numel(validDims{pfindex}),1]);
    fpcs{i}(validDims{pfindex}) = V(:,i);
end


figure;
for i = 1:10,  
    subplot(1,11,i);  imagesc(bins{:},reshape_eigen_vector(fpcs{i},pfd));axis('xy');  title(['Original PC',num2str(i)]);
end
subplot(1,11,11);hold('on');plot(diag(S(1:10,1:10)),'r-+');





LRop = [];
VTop = [];
for n = 1:100,
randPermInd = randperm(size(bmapOriginal,2));
[~,LR,FSr,VT] = erpPCA(bmapOriginal(:,randPermInd(1:400))',numComp);            % COMPUTE covariance-based PCA with Varimax rotation
LRop = cat(3,LRop,LR);
VTop = cat(2,VTop,VT(1:50,4));
end

figure,
for n =  1:5,
    subplot(5,5,n);hist(VTop(n,:),25)
end


k = 5;
fpct  = cell([1,numComp]);
for n = 1:4
    for i = 1:numComp,
        fpct{i} = nan([numel(validDims{pfindex}),1]);
        fpct{i}(validDims{pfindex}) = LRop(:,i,n+10);
        subplot(5,5,k+i)
        imagesc(bins{:},reshape_eigen_vector(fpct{i},pfd));axis('xy');  title(['Original PC',num2str(i)]);    
    end
    k = k+5;
end





fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan([numel(validDims{pfindex}),1]);
    fpc{i}(validDims{pfindex}) = eigVec{1}(:,i);
end

fpcOriBs  = cell([1,numComp]);
for i = 1:numComp,
    fpcOriBs{i} = nan([numel(validDims{pfindex}),numIter]);
    for n = 1:numIter,
        fpcOriBs{i}(validDims{pfindex},n) = LRO(:,i,n);        
    end
end


fpcShuff  = cell([1,numComp]);
for i = 1:numComp,
    fpcShuff{i} = nan([numel(validDims{pfindex}),numIter]);
    for n = 1:numIter,
        fpcShuff{i}(validDims{pfindex},n) = LRR(:,i,n);        
    end
end


figure();
subplot(4,5,1);  imagesc(bins{:},bmap');axis('xy');                               title('Original Map')
subplot(4,5,2);  imagesc(smootherWide); axis('xy');                               title('Smoothing Kernel')
subplot(4,5,3);  imagesc(bins{:},bmapSmooth');axis('xy');                         title('Smoothed Background Map');
subplot(4,5,4);  imagesc(bins{:},bmapMerge');axis('xy');                          title('Merged Original and Background Map');
subplot(4,5,5);  imagesc(binSubset{:},bmapMergeSubset');  axis('xy');             title('Merged Map Subset');
subplot(4,5,10);  imagesc(binSubset{:},bmapMergeSubsetShuffle');  axis('xy');     title('Random Circular Shift');
subplot(4,5,9);  imagesc(smootherNarrow);                                         title('Smoothing Kernel');
subplot(4,5,8);  imagesc(binSubset{:},bmapMergeSubsetShuffleSmooth');  axis('xy');title('Smoothed Shifted Map')
subplot(4,5,7);  imagesc(bins{:},bmapMergeShuffleSmoothFull(:,:,1)');  axis('xy');title('Full Shifted Map')


for i = 1:4,  subplot(4,5,10+i);  imagesc(bins{:},reshape_eigen_vector(fpc{i},pfd));axis('xy');  title(['Original PC',num2str(i)]);end
subplot(4,5,15);hold('on');plot(VT(1:10,end),'r-+');plot(mean(VTR(1:10,:),2),'b-+'); title('Eigenvalues');
legend({'Original','Shifted'});

for i = 1:4,  subplot(4,5,15+i);  imagesc(bins{:},reshape_eigen_vector(fpce{i},pfd));axis('xy');  title(['Shifted PC',num2str(i)]);end
subplot(4,5,20); plot((VT(1:20,4)-mean(VTR(i,1:32)))./std(VTR(i,1:32)),'b-+');                    title('z-score oriXshft');


print(gcf,'-dpng',fullfile('/storage/share/Projects/BehaviorPlaceCode/bhv_decomp',...
                             'eigenspectrum_rmap_random_shifts.png'));
print(gcf,'-depsc2',fullfile('/storage/share/Projects/BehaviorPlaceCode/bhv_decomp',...
                             'eigenspectrum_rmap_random_shifts.eps'));



reshape_eigen_vector_bs = @(vec) permute(reshape(vec(:,:),pfd{1}.adata.binSizes(1),pfd{1}.adata.binSizes(2),[]),[2,1,3,4]);
 
ovec = [];
for i = 1:5,
    ovec = cat(3,ovec,reshape_eigen_vector(fpc{i},pfd));
end
ovec(isnan(ovec)) = 0;
 
evec = [];
for i = 1:5,
    evec = cat(4,evec,reshape_eigen_vector_bs(fpcOriBs{i}));
end
evec(isnan(evec)) = 0;

svec = [];
for i = 1:5,
    svec = cat(4,svec,reshape_eigen_vector_bs(fpcShuff{i}));
end
svec(isnan(svec)) = 0;





vscale = abs(diff(bins{1}(1:2)));
cOEvec = [];
mOEvec = [];
for n = 1:numIter,
    for j = 1:5,
        for k = 1:5,
            [lm,m] = LocalMinimaN(-xcorr2(ovec(:,:,j),evec(:,:,n,k)),0,100);
            mOEvec(j,k,n) = -m;
            cOEvec(j,k,n) = sqrt(sum([vscale*(lm(1)-40),vscale*(lm(2)-40)].^2));
        end
    end
end




figure,
subplot(221);imagesc(mvec(:,:,end)');
subplot(222);imagesc(cvec(:,:,end)');
subplot(223);imagesc(mOEvec(:,:,end)');
subplot(224);imagesc(cOEvec(:,:,end)');


% minimize 
bmatch = [];
cOEvecSorted = [];
evecSorted = [];
for n = 1:numIter,
    for j = 1:5,    
        if j == 1,
            [~,bmatch(j,n)] = min(cOEvec(j,:,n));
        elseif j == 5,
            bmatch(j,n) = find(~ismember(1:5,bmatch(1:j-1,n)));
        else,
            [~,m] = min(cOEvec(j,~ismember(1:5,bmatch(1:j-1,n)),n));
            % m is a subset
            msubset = find(~ismember(1:5,bmatch(1:j-1,n)));
            bmatch(j,n) = msubset(m);
        end
        cOEvecSorted(j,n) = cOEvec(j,bmatch(j,n),n);
        evecSorted(:,:,n,j) = evec(:,:,n,bmatch(j,n));
    end
end



vscale = abs(diff(bins{1}(1:2)));
cvec = [];
mvec = [];
for n = 1:numIter,
    for j = 1:5,
        for k = 1:5,
            [lm,m] = LocalMinimaN(-xcorr2(evecSorted(:,:,n,j),svec(:,:,n,k)),0,100);
            mvec(j,k,n) = -m;
            cvec(j,k,n) = sqrt(sum([vscale*(lm(1)-40),vscale*(lm(2)-40)].^2));
        end
    end
end




% minimize 
bmatch = [];
sOEvecSorted = [];
for n = 1:numIter,
    for j = 1:5,    
        if j == 1,
            [~,bmatch(j,n)] = min(cvec(j,:,n));
        elseif j == 5,
            bmatch(j,n) = find(~ismember(1:5,bmatch(1:j-1,n)));
        else,
            [~,m] = min(cvec(j,~ismember(1:5,bmatch(1:j-1,n)),n));
            % m is a subset
            msubset = find(~ismember(1:5,bmatch(1:j-1,n)));
            bmatch(j,n) = msubset(m);
        end
        cvecSorted(j,n) = cvec(j,bmatch(j,n),n);
    end
end

eds = linspace(0,1.2,100);
figure,
bar(eds,histc(cvecSorted(3,:),eds),'histc');




figure,
for n = 1:5, subplot(3,5,n);  imagesc(evec(:,:,1,n));end
for n = 1:5, subplot(3,5,n+5);imagesc(svec(:,:,1,n));end
for n = 1:5, subplot(3,5,n+10);imagesc(ovec(:,:,n));end
    

testbinSizes = [20,20];
testbins =  {linspace(-1,1,20),linspace(-1,1,20)};

ssi = 0.1:0.025:0.5;

shuffle_testmap = @(m) m(randperm(numel(m)));

for s = 1:numel(ssi)
smoothingWeights = [ssi(s),ssi(s)];
sind = cell(1,ndims(testbins));
[sind{:}] = ndgrid(testbins{:});
for i = 1:ndims(bins),  sind{i} = sind{i}.^2/smoothingWeights(i)^2/2;  end
smootherTest = exp(sum(-cat(ndims(bins)+1,sind{:}),ndims(bins)+1));
smootherTest = smootherTest./sum(smootherTest(:));

testshift = zeros([400,10000]);
for n = 1:10000,
    testShift(:,n) = reshape(shuffle_bmap(smootherTest),[],1);
end

[~,LR,FSr,VT] = erpPCA(testShift',24);            % COMPUTE covariance-based PCA with Varimax rotation

pcCount(s) = sum(VT(:,4)>1);
end

figure
for i = 1:24,
    subplot(5,5,i);
    imagesc(reshape(LR(:,i),[20,20]));    
end
subplot(5,5,25);
plot(VT(1:24,4));

    
    
testshift = zeros([400,10000]);
for n = 1:10000,
    testShift(:,n) = shuffle_testmap(reshape(shuffle_bmap(smootherTest),[],1));
end
[~,LR,FSr,VT] = erpPCA(testShift',24,500);            % COMPUTE covariance-based PCA with Varimax rotation


figure
for i = 1:24,
    subplot(5,5,i);
    imagesc(reshape(LR(:,i),[20,20]));    
end
subplot(5,5,25);
plot(VT(1:24,4));


[U,S,V] = svd(testShift','econ');            % COMPUTE covariance-based PCA with Varimax rotation

figure
for i = 1:24,
    subplot(5,5,i);
    imagesc(reshape(V(:,i),[20,20]));    
end
subplot(5,5,25);
plot(diag(S(1:24,1:24)));

