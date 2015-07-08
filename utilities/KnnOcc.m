function [xyto,Bins] = KnnOcc(Session,state,varargin)
[binDims,nnn,dthresh,downSampleRate,type,distdw,distIndw] = DefaultArgs(varargin,{50,70,60,20,'xy',[],[]});

if Session.xyz.isempty, 
    Session.load('xyz');
end

pos = sq(Session.xyz(Session.stc{state},Session.trackingMarker,1:numel(binDims))); 


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
mywx = repmat(pos(spind,:),[1,1,nbins]);


xy = cell(2,1);
[xy{:}]= meshgrid(Bins{:});
xy = repmat(permute(reshape(cat(3,xy{:}),nbins,2),fliplr(1:3)),size(mywx,1),1);
distw = sqrt(sum((mywx-xy).^2,2));
[distdw,distIndw] = sort(distw,1,'ascend');  
pfknnmdw = permute(median(distdw(1:nnn,:,:,:),1),[2,3,1]);

% $$$ xyto = zeros(nbins,1);
% $$$ for i = 1:size(distw,3),
% $$$     tsi = find(diff(distw(:,1,i)<dthresh)>0);
% $$$     tei = find(diff(distw(:,1,i)<dthresh)<0);
% $$$     xyto(i) = sum(abs(diff([tsi(2:numel(tsi)),tei(1:numel(tsi)-1)],1,2))>120);
% $$$ end
% $$$ 
% $$$ xyto = reshape(xyto,Nbin');