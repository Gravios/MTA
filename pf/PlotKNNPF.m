function [RateMap,Bins,MRate,SI,Spar]= PlotKNNPF(Session,ufr,pos,varargin)
[binDims,SmoothingWeights,type] = DefaultArgs(varargin,{50,[],'xy'});

% Session = MTATrial('jg05-20120317');
% Session.xyz.load(Session);
% Session = Session.filter('xyz');
% ufr = Session.ufr.create(Session,Session.xyz,[],[],0.2);
% 
% mywx = Session.xyz(Session.stc{'w'},7,[1,2]);
% mywu = ufr(Session.stc{'w'},:);
% spind = 1:10:size(mywx,1);
% mywx = mywx(spind,:);
% mywu = mywu(spind,:);
ndims = numel(binDims);
bound_lims = Session.maze.boundaries(1:ndims,:);
Nbin = round(diff(bound_lims,1,2)./binDims');


Bins = cell(1,ndims);
k = Nbin./diff(bound_lims,1,2);
msize = round(sum(abs(bound_lims),2).*k);
for i = 1:ndims
    Bins{i} = ([1:msize(i)]'-1)./repmat(k(i)',msize(i),1)+repmat(bound_lims(i,1),msize(i),1)+round(repmat(k(i)',msize(i),1).^-1/2);
end

spind = 1:10:size(mywx,1);
mywx = pos(spind,:);
mywu = ufr(spind,:);

nxbins = 50;
nybins = 50;
xbins = linspace(Session.maze.boundaries(1,1),Session.maze.boundaries(1,2),nxbins)';
ybins = linspace(Session.maze.boundaries(2,2),Session.maze.boundaries(2,1),nxbins)';
nnn = 70;
pfknnmrw = nan(,ufs.size(2));
pfknnmdw = nan(nxbins,nybins);
dthresh = 40;

xy = cell(2,1);
[xy{:}]= meshgrid(xbins,ybins);

xy = repmat(permute(cat(3,xy{:}),[4,3,1,2]),size(mywx,1),1);
mywxt = repmat(mywx,[1,1,nxbins,nybins]);
mywut = repmat(mywu,[1,1,nxbins,nybins]);
distw = sqrt(sum((mywxt-xy).^2,2));
[distdw,distIndw ] = sort(distw,1,'ascend');  

pfknnmdw = sq(median(distdw(1:nnn,:,:,:),1));

s = size(mywut);
cm = reshape(repmat([0:s(1):prod(s)-1],s(1),1),s);     
smywut = mywut(repmat(distIndw,1,ufr.size(2))+cm);
pfknnmrw = permute(mean(smywut(1:nnn,:,:,:),1),[3,4,2,1]);
pfknnmrw(repmat(pfknnmdw,[1,1,ufr.size(2)])>dthresh) = nan;
