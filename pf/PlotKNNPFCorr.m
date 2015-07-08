function [RateMap,Bins,distdw,distIndw]= PlotKNNPFCorr(Session,ufr,pos,varargin)
[binDims,nnn,dthresh,downSampleRate,type,distdw,distIndw] = DefaultArgs(varargin,{50,70,60,20,'xy',[],[]});

ndims = numel(binDims);
bound_lims = Session.maze.boundaries(1:ndims,:);
Nbin = round(diff(bound_lims,1,2)./binDims');
nbins = prod(Nbin);

Bins = cell(1,ndims);
k = Nbin./diff(bound_lims,1,2);
msize = round(sum(abs(bound_lims),2).*k);
for i = 1:ndims
    Bins{i} = ([1:msize(i)]'-1)./repmat(k(i)',msize(i),1)+repmat(bound_lims(i,1),msize(i),1)+round(repmat(k(i)',msize(i),1).^-1/2);
end

spind = 1:round(Session.xyz.sampleRate/downSampleRate):size(pos,1);
mywu = repmat(ufr(spind,:),[1,1,nbins]);

if isempty(distdw)&&isempty(distIndw)
    spind = 1:round(Session.xyz.sampleRate/downSampleRate):size(pos,1);
    mywx = repmat(pos(spind,:),[1,1,nbins]);
    xy = cell(2,1);
    [xy{:}]= meshgrid(Bins{:});
    xy = repmat(permute(reshape(cat(3,xy{:}),nbins,2),fliplr(1:3)),size(mywx,1),1);
    distw = sqrt(sum((mywx-xy).^2,2));
    [distdw,distIndw] = sort(distw,1,'ascend');  
end

pfknnmdw = permute(median(distdw(1:nnn,:,:,:),1),[2,3,1]);

s = size(mywu);
cm = reshape(repmat([0:s(1):prod(s)-1],s(1),1),s);     
smywut = mywu(repmat(distIndw,1,size(ufr,2))+cm);

RateMap = permute(1/(nnn-1).*sum((distdw(1:nnn,:,:,:)-repmat(mean(distdw(1:nnn,:,:,:),1),[nnn,1,1,1])).*(smywut(1:nnn,:,:,:)-repmat(mean(smywut(1:nnn,:,:,:),1),[nnn,1,1,1])),1)./sqrt(var(distdw(1:nnn,:,:,:),[],1).*var(smywut(1:nnn,:,:,:),[],1)),[ 3,2,1]);



RateMap(repmat(pfknnmdw,[1,1,size(ufr,2)])>dthresh) = nan;
