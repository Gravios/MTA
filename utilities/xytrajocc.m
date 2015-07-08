function [xyto,Bins] = xytrajocc(Session,state,varargin)
%function [xyto,Bins] = xytrajocc(Session,pos,varargin)
%[binDims,dthresh,type] = DefaultArgs(varargin,{[20,20],80,'xy'});
[overwrite,binDims,dthresh,type] = DefaultArgs(varargin,{false,[20,20],120,'xy'});



if ischar(state), 
    state = Session.stc{state,Session.xyz.sampleRate}; 
end

if ~exist(fullfile(Session.spath,[Session.filebase,'xytrajocc','.',state.label,'.mat']),'file')||overwrite

    xyz = Session.load('xyz');
    xyz.filter(gtwin(.4,xyz.sampleRate));
    pos = sq(xyz(state,Session.trackingMarker,[1,2]));

    stateBoundaries = cumsum(diff(state.data,1,2));

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

    mywx = repmat(pos,[1,1,nbins]);

    xy = cell(2,1);
    [xy{:}]= meshgrid(Bins{:});
    xy = repmat(permute(reshape(cat(3,xy{:}),nbins,2),fliplr(1:3)),size(mywx,1),1);


    distw = sqrt(sum((mywx-xy).^2,2));

    clear('xy')
    clear('mywx')

    distw(distw<dthresh) = 1;
    distw(distw>dthresh) = 0;
    distw(stateBoundaries,:,:) = 0;

    xyto = zeros([size(distw,3),1]);
    for i = 1:size(distw,3),
        xyto(i) = size(ThreshCross(distw(:,1,i),.5,round(.25*xyz.sampleRate)),1);
    end


    xyto = reshape(xyto,Nbin');

    save(fullfile(Session.spath,[Session.filebase,'xytrajocc','.',state.label,'.mat']),...
         'xyto','Bins');
else
    load(fullfile(Session.spath,[Session.filebase,'xytrajocc','.',state.label,'.mat']));
end