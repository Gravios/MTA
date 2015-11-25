function [RateMap,Bins,distdw,distIndw]= PlotKNNPF(Session,ufr,pos,varargin)
%function [RateMap,Bins,distdw,distIndw]= PlotKNNPF(Session,ufr,pos,varargin)
%[binDims,nnn,dthresh,downSampleRate,type,distdw,distIndw] = DefaultArgs(varargin,{50,70,.5*Session.xyz.sampleRate,20,'xy',[],[]});
[binDims,nnn,dthresh,type,distdw,distIndw,bound_lims,stat_fun] = DefaultArgs(varargin,{50,70,.5*Session.xyz.sampleRate,'xy',[],[],[],@nanmean},1);


ndims = numel(binDims);
if isempty(bound_lims),bound_lims = Session.maze.boundaries(1:ndims,:);end
Nbin = round(diff(bound_lims,1,2)./binDims');
nbins = prod(Nbin);

Bins = cell(1,ndims);
k = Nbin./diff(bound_lims,1,2);
msize = round(sum(abs(bound_lims),2).*k);
for i = 1:ndims
    Bins{i} = ([1:msize(i)]'-1)./repmat(k(i)',msize(i),1)+repmat(bound_lims(i,1),msize(i),1)+round(repmat(k(i)',msize(i),1).^-1/2);
end



mywu = ufr.copy;
%mywu.data = mywu(:,1);
mywu.data = repmat(mywu.data,[1,1,nbins]);

if isempty(distdw)&&isempty(distIndw)
    mywx = repmat(pos.data,[1,1,nbins]);
    xy = cell(ndims,1);
    [xy{:}]= meshgrid(Bins{:});
    xy = repmat(permute(reshape(cat(3,xy{:}),nbins,ndims),fliplr(1:3)),size(mywx,1),1);
    distw = sqrt(nansum((mywx-xy).^2,2));
    [distdw,distIndw] = sort(distw,1,'ascend');  
end


if size(distdw,1)>nnn,
    pfknnmdw = permute(nanmedian(distdw(1:nnn,:,:,:),1),[2,3,1]);
    s = size(mywu);
    cm = reshape(repmat([0:s(1):prod(s)-1],s(1),1),s);     
    smywut = mywu.data(repmat(distIndw,[1,size(mywu,2)])+cm);
    RateMap = permute(feval(stat_fun,smywut(1:nnn,:,:,:)),[3,2,1]);
    RateMap(repmat(pfknnmdw,[1,1,size(ufr,2)])>dthresh) = nan;
else
    RateMap = nan([size(distw,3),size(ufr,2)]);
end