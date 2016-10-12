
% MTAAknnpfs_bs Parameter search to 




% KNNPFS DEFARGS ---------------------------------------------------------------
defargs = struct('units',              [],                                   ...
                 'states',             {{'walk'}},                           ...
                 'overwrite',          false,                                ...
                 'tag',                [],                                   ...
                 'ufr',                [],                                   ...
                 'binDims',            [20,20],                              ...
                 'nNearestNeighbors' , 30,                                   ...
                 'distThreshold',      125,                                  ...
                 'type',               'xy',                                 ...
                 'ufrShufBlockSize',   1,                                    ... 
                 'numIter',            1,                                    ...
                 'pos',                [],                                   ...
                 'sampleRate',         4,                                    ...
                 'absTimeSubSample',   300                                   ...
);%-----------------------------------------------------------------------------



Trial = MTATrial('jg05-20120310');

xyz = Trial.load('xyz');
xyz.resample(defargs.sampleRate);

pos = xyz(Trial.stc{'w'},5,:);


%% Check change in sample variance with samples---------------------------------




apf = struct('overwrite',          false, ...
             'binDims',            [20,20],...
             'nNearestNeighbors',  60,...
             'ufrShufBlockSize',   1,...
             'distThreshold',      125,...
             'numIter',            1001,...
             'sampleRate',         10);
apf = struct2varargin(apf);
pf = MTAAknnpfs_bs(Trial, apf{:});

units = pf.data.clu;

pfstats = {};
pfshuff = {};
for unit = units;
    [pfstats{unit},pfshuff{unit}] = PlaceFieldStats(Trial,pf,unit);
end


figure,hist(pfshuff{10}.patchPFR(1,:,1),100)



figure,hist(pfshuff{10}.patchCOM(1,:,1,2),100)


accVar = nan([1001,2000]);
for i = 1:2000,
for n = 10:1001,
    accVar(n,i) = var(pfshuff{10}.patchCOM(1,randi(1001,[n,1]),1,2));
end
end


accVarStd = std(accVar(10:100:end,:),[],2);
figure,errorbar(1:100:1001-9,mean(accVar(10:100:end,:),2), ...
                accVarStd*2,-accVarStd*2);




%% Start parameter search

Trial = MTATrial('jg05-20120310');

% KNNPFS DEFARGS ---------------------------------------------------------------
defargs = struct('units',              10,                                   ...
                 'states',             {{'walk'}},                           ...
                 'overwrite',          false,                                ...
                 'tag',                [],                                   ...
                 'ufr',                [],                                   ...
                 'binDims',            [20,20],                              ...
                 'nNearestNeighbors' , 30,                                   ...
                 'distThreshold',      125,                                  ...
                 'type',               'xy',                                 ...
                 'ufrShufBlockSize',   1,                                    ... 
                 'numIter',            10001,                                ...
                 'pos',                [],                                   ...
                 'sampleRate',         4,                                    ...
                 'absTimeSubSample',   300                                   ...
);%-----------------------------------------------------------------------------


defargs = struct2varargin(defargs);


pf = MTAAknnpfs_bs(Trial,defargs{:});



[pfst,pfsf] = PlaceFieldStats(Trial,pf,10);


width = pf.adata.binSizes(1);
height = pf.adata.binSizes(2);
radius = round(pf.adata.binSizes(1)/2)-find(pf.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;



ratemap = reshape(pf.data.rateMap,[pf.adata.binSizes',10001]);
ratemap = bsxfun(@times,ratemap,mask);
ratemap(isnan(ratemap)) = -1;



figure;

subplot(231);
imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap(:,:,1)');
axis xy
colormap([0,0,0;parula]);
caxis([-1,10]);
hold on
plot(mean(pfsf.patchCOM(1,:,1,1)),mean(pfsf.patchCOM(1,:,1,2)),'*')

subplot(232);
imagesc(pf.adata.bins{1},pf.adata.bins{2},median(ratemap,3)');
axis xy
colormap([0,0,0;parula]);
caxis([-1,10]);
hold on
plot(mean(pfsf.patchCOM(1,:,1,1)),mean(pfsf.patchCOM(1,:,1,2)),'*')

subplot(233);
ratemap(ratemap<prctile(ratemap(nniz(ratemap(:))),95)) = -1;
imagesc(pf.adata.bins{1},pf.adata.bins{2},median(ratemap,3)');
axis xy
colormap([0,0,0;parula]);
caxis([-1,10]);
hold on
plot(mean(pfsf.patchCOM(1,:,1,1)),mean(pfsf.patchCOM(1,:,1,2)),'*')

subplot(234);
hist(pfsf.patchCOM(1,:,1,1),100)

subplot(235)
hist(pfsf.patchCOM(1,:,1,2),100)

subplot(236);
plot(pfsf.patchCOM(1,:,1,1),pfsf.patchCOM(1,:,1,2),'.')
xlim([-500,500])
ylim([-500,500])