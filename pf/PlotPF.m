function [RateMap, Bins, MRate, SI, Spar] = PlotPF(Session,spkpos,pos,varargin)
[binDims,SmoothingWeights,type,bound_lims] = DefaultArgs(varargin,{50,[],'xy',[]});

ndims = numel(binDims);

if isempty(bound_lims),
    bound_lims = bound_lims(1:ndims,:);
end

Nbin = round(abs(diff(bound_lims,1,2))./binDims');

if isempty(SmoothingWeights)
  SmoothingWeights = Nbin./30;
end

%% Constraint to maze is forced
switch type
    case 'pfcrz'
        switch Session.maze.shape
            case 'circle'
                Xmin = 0;
                Xmax = Session.maze.boundaries(1,2)-Session.maze.boundaries(1,1);
                Ymin = Session.maze.boundaries(3,1);
                Ymax = Session.maze.boundaries(3,2);
                %% scaling factor for rounding position
                dx = Xmax - Xmin;
                dy = Ymax - Ymin;
                k = Nbin./[dx;dy];
                %% matrix size
                msize = round([sum(abs([Xmin,Xmax]))*k(1) sum(abs([Ymin,Ymax]))*k(2)])';
                Bin1 = ([1:msize(1)]-1)/k(1);
                Bin2 = ([1:msize(2)]-1)/k(2);
                Bins = {Bin1(:),Bin2(:)};
                bound_lims = [Xmin,Xmax;Ymin,Ymax];
            case 'square'
        end
        
    otherwise
        Bins = cell(1,ndims);
        k = Nbin./abs(diff(bound_lims,1,2));
        msize = round(abs(diff(bound_lims,1,2)).*k);
        for i = 1:ndims
            Bins{i} = ([1:msize(i)]'-1)./repmat(k(i)',msize(i),1)+repmat(bound_lims(i,1),msize(i),1)+round(repmat(k(i)',msize(i),1).^-1/2);
        end
end



%% rounded position and removal of bins outside the computational volume 
Pos = round((pos-repmat(bound_lims(:,1)',size(pos,1),1)).*repmat(k',size(pos,1),1))+1;
for i = 1:ndims   
    Pos(Pos(:,i)<1|Pos(:,i)>Nbin(i),:) = [];
end
Occupancy = accumarray(Pos,1,msize')./Session.xyz.sampleRate;

spkpos = round((spkpos-repmat(bound_lims(:,1)',size(spkpos,1),1)).*repmat(k',size(spkpos,1),1))+1;
for i = 1:ndims
    spkpos(spkpos(:,i)<1|spkpos(:,i)>Nbin(i),:) = [];
end
SpikeCount = accumarray(spkpos,1,msize');



sind = cell(1,ndims);
for i = 1:ndims,
    sind{i} = linspace(-round(msize(i)/2),round(msize(i)/2),msize(i));
end
[sind{:}] = ndgrid(sind{:});
for i = 1:ndims,
    sind{i} = sind{i}.^2/SmoothingWeights(i)^2/2;
end
Smoother = exp(sum(-cat(ndims+1,sind{:}),ndims+1));
Smoother = Smoother./sum(Smoother(:));


SOcc = convn(Occupancy,Smoother,'same');
SCount = convn(SpikeCount,Smoother,'same');


OccThresh = 0.06;%0.12;%
%% Find the total occupancy and each pixels 
%% probability of occupancy
gtind = SOcc>OccThresh;
TotOcc = sum(SOcc(gtind));
POcc = SOcc./TotOcc;
%% Rate Map
RateMap = NaN(prod(Nbin),1);
RateMap(gtind) = SCount(gtind)./SOcc(gtind);
RateMap = RateMap(:);
RateMap(~gtind) = NaN;
%% Find the units overall mean rate given the 
%% current state
MRate = sum(SCount(gtind))/TotOcc;
SI = nansum(POcc(gtind).*(RateMap(gtind)./MRate).*log2(RateMap(gtind)./MRate));
Spar = 1/nansum(POcc(gtind).*RateMap(gtind).^2./MRate.^2);

