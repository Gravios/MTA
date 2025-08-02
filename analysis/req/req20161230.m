
Trial = MTATrial('Ed10-20140905');

Pfs = MTAApfs(Trial,[],'theta',...
              false,...
              'numIter',1,...
              'binDims',[20,20,20],...
              'type','xyz',...
              'SmoothingWeights',[2.2,2.2,2.2]);




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

unit = 36
rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),3);

rateMap = reshape(rateMap,numel(bin1),numel(bin2),numel(bin3)).*repmat(mask,[1,1,numel(bin3)]);

pfSize = [numel(bin1),numel(bin2)];

figure, hold on
for i = 1:3:18;
    surf(bin3(i)*ones(pfSize),rateMap(:,:,i)','EdgeColor','none')
end
view([50,10])
whitebg('k')
grid on
colorbar