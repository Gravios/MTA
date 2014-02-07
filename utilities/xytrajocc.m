function [xyto,Bins] = xytrajocc(Session,pos,varargin)

[binDims,dthresh,type] = DefaultArgs(varargin,{[20,20],80,'xy'});

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

spind = 1:round(Session.xyz.sampleRate/60):size(pos,1);
mywx = repmat(pos(spind,:),[1,1,nbins]);


xy = cell(2,1);
[xy{:}]= meshgrid(Bins{:});
xy = repmat(permute(reshape(cat(3,xy{:}),nbins,2),fliplr(1:3)),size(mywx,1),1);

distw = sqrt(sum((mywx-xy).^2,2));

xyto = zeros(nbins,1);
for i = 1:size(distw,3),
    tsi = find(diff(distw(:,1,i)<dthresh)>0);
    tei = find(diff(distw(:,1,i)<dthresh)<0);
    xyto(i) = sum(abs(diff([tsi(2:numel(tsi)),tei(1:numel(tsi)-1)],1,2))>120);
end

xyto = reshape(xyto,Nbin');