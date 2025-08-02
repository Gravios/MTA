

% Pitch related firing

Trial = MTATrial('jg05-20120310');
Trial.load('stc','NN0317R');
xyz = Trial.load('xyz');
fxyz = xyz.copy();
fxyz.filter('ButFilter',3,0.1,'low');
fang = create(MTADang,Trial,fxyz);



figure,hold on,
ind = 1:1000;
plot(xyz(ind,1,1),xyz(ind,1,2))
plot(fxyz(ind,1,1),fxyz(ind,1,2),'.')

figure,
plot(fang(:,1,5,1))

pos = xyz.copy();
pos.data = cat(3,xyz(:,7,[1,2]),fang(:,1,7,1));

Trial.maze.boundaries(3,:) = [-pi,pi];

Pfs = MTAApfs(Trial,         ... MTATrial
                 [],         ... units
                 'loc&theta',... states
                 true,       ... overwrite
                 [],         ... tag
                 [20,20,0.2],... BinDims
                 [2.2,2.2,2.2], ... SmoothingWeights
                 'xya',      ... type
                 [],         ... spkShuffle
                 [],         ... posShuffle
                 1,          ... numIter
                 pos,         ... pos
                 [-500,500;-500,500;-pi,pi]... bound_lims
);


bin1 = Pfs.adata.bins{1};
bin2 = Pfs.adata.bins{2};
bin3 = Pfs.adata.bins{3};

width = Pfs.adata.binSizes(1);
height = Pfs.adata.binSizes(2);
radius = round(Pfs.adata.binSizes(1)/2)-find(Pfs.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;

unit = 107;
rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),3);

rateMap = reshape(rateMap,numel(bin1),numel(bin2),numel(bin3)).*repmat(mask,[1,1,numel(bin3)]);

pfSize = [numel(bin1),numel(bin2)];

figure, hold on
for i = 1:4:numel(bin3);
    surf(bin3(i)*ones(pfSize),rateMap(:,:,i)','EdgeColor','none')
end
view([50,10])
whitebg('k')
grid on
colorbar



% Height related

Trial = MTATrial('jg05-20120310');
Trial.load('stc','NN0317R');
xyz = Trial.load('xyz');
fxyz = xyz.copy();
fxyz.filter('ButFilter',3,1,'low');

pos = xyz.copy();
pos.data = cat(3,xyz(:,7,[1,2]),fxyz(:,7,3));

Pfs = MTAApfs(Trial,         ... MTATrial
                 [],         ... units
                 'loc&theta',... states
                 false,       ... overwrite
                 [],         ... tag
                 [20,20,5],... BinDims
                 [2.2,2.2,3.2], ... SmoothingWeights
                 'xya',      ... type
                 [],         ... spkShuffle
                 [],         ... posShuffle
                 1,          ... numIter
                 pos,         ... pos
                 [-500,500;-500,500;30,160]... bound_lims
);


bin1 = Pfs.adata.bins{1};
bin2 = Pfs.adata.bins{2};
bin3 = Pfs.adata.bins{3};

width = Pfs.adata.binSizes(1);
height = Pfs.adata.binSizes(2);
radius = round(Pfs.adata.binSizes(1)/2)-find(Pfs.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;

figure,
for unit = 1:20,
    clf,hold on,
rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),3);
rateMap = reshape(rateMap,numel(bin1),numel(bin2),numel(bin3)).*repmat(mask,[1,1,numel(bin3)]);
pfSize = [numel(bin1),numel(bin2)];


for i = 1:4:numel(bin3);
    surf(bin3(i)*ones(pfSize),rateMap(:,:,i)','EdgeColor','none')
end
view([50,10])
whitebg('k')
grid on
colorbar
view([-70,10])
pause(0.5)
end



Pfs = MTAApfs(Trial,            ... MTATrial
                 [7:14],     ... units
                 'theta-groom-sit',   ... states
                 false,         ... overwrite
                 [],            ... tag
                 [20,20,20],     ... BinDims
                 [2.2,2.2,2.2], ... SmoothingWeights
                 'xya',         ... type
                 [],            ... spkShuffle
                 true,          ... posShuffle
                 1000,          ... numIter
                 pos,           ... pos
                 [-500,500;-500,500;0,250]... bound_lims
);




bin1 = Pfs.adata.bins{1};
bin2 = Pfs.adata.bins{2};
bin3 = Pfs.adata.bins{3};

width = Pfs.adata.binSizes(1);
height = Pfs.adata.binSizes(2);
radius = round(Pfs.adata.binSizes(1)/2)-find(Pfs.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;

figure,
clf,hold on,
rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),3);
rateMap = reshape(rateMap,numel(bin1),numel(bin2),numel(bin3)).*repmat(mask,[1,1,numel(bin3)]);
pfSize = [numel(bin1),numel(bin2)];

for i = 1:2:numel(bin3);
    surf(bin3(i)*ones(pfSize),rateMap(:,:,i)','EdgeColor','none')
end
view([50,10])
whitebg('k')
grid on
colorbar
view([-70,10])
pause(0.5)

LocalMinimaN(-rateMap,-4,3)
