% req20180630 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Author: Justin Graboski
%  Final_Forms: NA
%  Description: validate erpPCA eigenvectors of behavioral placefields
%  Bugs: NA

% Load and preprocess data
% Generate resampled and surrogate data sets
% Concatenation of eigenvectors from different data treatments
% Matching and sorting eigenvectors 



%% Load and preprocess data

% LOAD session data
MjgER2016_load_data();

% LOAD behavioral data
[pfd ,tags ,eigVec, eigVar, eigScore, validDims,...
 unitSubsets, unitIntersection, zrmMean, zrmStd] = req20180123_ver5(Trials);
numComp = size(eigVec{1},2);
pfindex = 1;

% LOAD behavioral scores
MjgER2016_load_bhv_erpPCA_scores();


%% Generate resampled and surrogate data sets

% GET behavior field information
binSizes = pfd{1}.adata.binSizes';
bins =  pfd{1}.adata.bins;

% HELPER function for sample randomization
shuffle_bmap = @(m) circshift(circshift(m,randi(size(m,1)),1),randi(size(m,2)),2);


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



% GENERATE bootstrapped eigenvectors based on two-thirds sample
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


% GENERATE bootstrapped eigenvectors based on two-thirds sample
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


% $$$ % PLOT eigenvectors
% $$$ k = 5;
% $$$ fpct  = cell([1,numComp]);
% $$$ for n = 1:4
% $$$     for i = 1:numComp,
% $$$         fpct{i} = nan([numel(validDims{pfindex}),1]);
% $$$         fpct{i}(validDims{pfindex}) = LRop(:,i,n+10);
% $$$         subplot(5,5,k+i)
% $$$         imagesc(bins{:},reshape_eigen_vector(fpct{i},pfd));axis('xy');  title(['Original PC',num2str(i)]);    
% $$$     end
% $$$     k = k+5;
% $$$ end



% RESHAPE eigenvectors into behavioral space
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



% DIAGNOSTIC figure
% $$$ figure();
% $$$ subplot(4,5,1);  imagesc(bins{:},bmap');axis('xy');                               title('Original Map')
% $$$ subplot(4,5,2);  imagesc(smootherWide); axis('xy');                               title('Smoothing Kernel')
% $$$ subplot(4,5,3);  imagesc(bins{:},bmapSmooth');axis('xy');                         title('Smoothed Background Map');
% $$$ subplot(4,5,4);  imagesc(bins{:},bmapMerge');axis('xy');                          title('Merged Original and Background Map');
% $$$ subplot(4,5,5);  imagesc(binSubset{:},bmapMergeSubset');  axis('xy');             title('Merged Map Subset');
% $$$ subplot(4,5,10);  imagesc(binSubset{:},bmapMergeSubsetShuffle');  axis('xy');     title('Random Circular Shift');
% $$$ subplot(4,5,9);  imagesc(smootherNarrow);                                         title('Smoothing Kernel');
% $$$ subplot(4,5,8);  imagesc(binSubset{:},bmapMergeSubsetShuffleSmooth');  axis('xy');title('Smoothed Shifted Map')
% $$$ subplot(4,5,7);  imagesc(bins{:},bmapMergeShuffleSmoothFull(:,:,1)');  axis('xy');title('Full Shifted Map')
% $$$ 
% $$$ for i = 1:4,  subplot(4,5,10+i);  imagesc(bins{:},reshape_eigen_vector(fpc{i},pfd));axis('xy');  title(['Original PC',num2str(i)]);end
% $$$ subplot(4,5,15);hold('on');plot(VT(1:10,end),'r-+');plot(mean(VTR(1:10,:),2),'b-+'); title('Eigenvalues');
% $$$ legend({'Original','Shifted'});
% $$$ 
% $$$ for i = 1:4,  subplot(4,5,15+i);  imagesc(bins{:},reshape_eigen_vector(fpce{i},pfd));axis('xy');  title(['Shifted PC',num2str(i)]);end
% $$$ subplot(4,5,20); plot((VT(1:20,4)-mean(VTR(i,1:32)))./std(VTR(i,1:32)),'b-+');                    title('z-score oriXshft');
% $$$ 
% $$$ print(gcf,'-dpng',fullfile('/storage/share/Projects/BehaviorPlaceCode/bhv_decomp',...
% $$$                              'eigenspectrum_rmap_random_shifts.png'));
% $$$ print(gcf,'-depsc2',fullfile('/storage/share/Projects/BehaviorPlaceCode/bhv_decomp',...
% $$$                              'eigenspectrum_rmap_random_shifts.eps'));





%% Concatenation of eigenvectors from different data treatments

reshape_eigen_vector_bs = @(vec) permute(reshape(vec(:,:),pfd{1}.adata.binSizes(1),pfd{1}.adata.binSizes(2),[]),[2,1,3,4]);

% CONCATENATE eigenvectors of original data
eigenVecOri = [];
for i = 1:5,
    eigenVecOri = cat(3,eigenVecOri,reshape_eigen_vector(fpc{i},pfd));
end
eigenVecOri(isnan(eigenVecOri)) = 0;
% CONCATENATE eigenvectors of bootstrapped original data
eigenVecBs = [];
for i = 1:5,
    eigenVecBs = cat(4,eigenVecBs,reshape_eigen_vector_bs(fpcOriBs{i}));
end
eigenVecBs(isnan(eigenVecBs)) = 0;
% CONCATENATE eigenvectors of shuffled data
eigenVecShuff = [];
for i = 1:5,
    eigenVecShuff = cat(4,eigenVecShuff,reshape_eigen_vector_bs(fpcShuff{i}));
end
eigenVecShuff(isnan(eigenVecShuff)) = 0;




%% Matching and sorting eigenvectors 

% FAILED METHOD
% TITLE       : Match eigenvectors XCORR
% DESCRIPTION : cross correlation to match similar eigenvectors
% ASSUMPTIONS : All major eigenvectors have most of their weights loaded locally in the form of a gaussian distribution
% $$$ vscale = abs(diff(bins{1}(1:2)));
% $$$ distanceOriBs = [];
% $$$ xcorrOriBs = [];
% $$$ for n = 1:numIter,
% $$$     for j = 1:5,
% $$$         for k = 1:5,
% $$$ % FIND shift where xcorr is maximized 
% $$$             [lm,m] = LocalMinimaN(-xcorr2(eigenVecOri(:,:,j),eigenVecBs(:,:,n,k)),0,100);
% $$$             xcorrOriBs(j,k,n) = -m;
% $$$ % COMPUTE distance of shift
% $$$             distanceOriBs(j,k,n) = sqrt(sum([vscale*(lm(1)-40),vscale*(lm(2)-40)].^2));
% $$$         end
% $$$     end
% $$$ end
% $$$ figure,
% $$$ subplot(221);imagesc(mvec(:,:,end)');
% $$$ subplot(222);imagesc(cvec(:,:,end)');
% $$$ subplot(223);imagesc(xcorrOriBs(:,:,end)');
% $$$ subplot(224);imagesc(distanceOriBs(:,:,end)');
% POF        : non peak values were large enough to influence the cross correlation, such that the
%            | resulting shift did not reflect the relationship between eigenvector peaks.


% ACTIVE METHOD
% TITLE       : Match eigenvectors LM
% DESCRIPTION : eigenvector Local Maxima to match similar eigenvectors
% ASSUMPTIONS : All major eigenvectors have most of their weights loaded locally in the form of a gaussian distribution
distanceOriBs = [];
vscale = abs(diff(bins{1}(1:2)));
for n = 1:numIter,
    for j = 1:5,
        for k = 1:5,
% FIND shift where xcorr is maximized 
            [lmOri] = LocalMinima2(-eigenVecOri(:,:,j),0,100);
            [lmBs]  = LocalMinima2(-eigenVecBs(:,:,n,k),0,100);
% COMPUTE distance of shift
            distanceOriBs(j,k,n) = sqrt(sum([lmOri-lmBs].^2)).*vscale;
        end
    end
end





% SORT bootstrapped eigenvectors to match original eigenvectors
bmatch = [];
distanceOriBsSorted = [];
eigenVecBsSorted = [];
for n = 1:numIter,
    for j = 1:5,    
        if j == 1,
            [~,bmatch(j,n)] = min(distanceOriBs(j,:,n));
        elseif j == 5,
            bmatch(j,n) = find(~ismember(1:5,bmatch(1:j-1,n)));
        else,
            [~,m] = min(distanceOriBs(j,~ismember(1:5,bmatch(1:j-1,n)),n));
            % m is a subset
            msubset = find(~ismember(1:5,bmatch(1:j-1,n)));
            bmatch(j,n) = msubset(m);
        end
        distanceOriBsSorted(j,n) = distanceOriBs(j,bmatch(j,n),n);
        eigenVecBsSorted(:,:,n,j) = eigenVecBs(:,:,n,bmatch(j,n));
    end
end




vscale = abs(diff(bins{1}(1:2)));
cvec = [];
mvec = [];
for n = 1:numIter,
    for j = 1:5,
        for k = 1:5,
            [lm,m] = LocalMinimaN(-xcorr2(eigenVecOri(:,:,j),eigenVecShuff(:,:,n,k)),0,100);
            mvec(j,k,n) = -m;
            cvec(j,k,n) = sqrt(sum([vscale*(lm(1)-40),vscale*(lm(2)-40)].^2));
        end
    end
end




% SORT shuffled eigenvectors according to best fit to original eigenvector
bmatch = [];
eigenVecShuffSorted = [];
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
        sOEvecSorted(j,n) = cvec(j,bmatch(j,n),n);        
        eigenVecShuffSorted(:,:,j,n) = eigenVecShuff(:,:,n,bmatch(j,n));
    end
end

eds = linspace(0,1.2,100);
figure,
bar(eds,histc(cvecSorted(3,:),eds),'histc');

figure();
hist(sOEvecSorted(3,:),100)

figure();
subplot(131);imagesc(eigenVecOri(:,:,1)');
subplot(132);imagesc(eigenVecShuff(:,:,1)');




j = 97;
figure,
for n = 1:5, subplot(3,5,n);  imagesc(eigenVecBs(:,:,j,n));end
for n = 1:5, subplot(3,5,n+5);imagesc(eigenVecShuff(:,:,j,n));end
for n = 1:5, subplot(3,5,n+10);imagesc(eigenVecOri(:,:,n));end
    

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

