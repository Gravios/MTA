function [RateMap, Bins, SI, Spar] = PlotPF(Session,spkpos,pos,varargin)

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('binDims',                                50,                                   ...
                 'SmoothingWeights',                       [],                                   ...
                 'type',                                   'xy',                                 ...
                 'bound_lims',                             [],                                   ...
                 'posSampleRate',                          Session.xyz.sampleRate                ...
);
[binDims,SmoothingWeights,type,bound_lims,posSampleRate] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

SI = [];
Spar = [];

% SET 
ndims = numel(binDims);
if isempty(bound_lims),  bound_lims = Trial.maze.boundaries(ismember('xyz',type),:);  end
Nbin = round(abs(diff(bound_lims,1,2))./binDims');
if isempty(SmoothingWeights),  SmoothingWeights = Nbin./30;  end


% INITIALIZE bins
Bins = cell(1,ndims);
k = Nbin./abs(diff(bound_lims,1,2));
msize = round(abs(diff(bound_lims,1,2)).*k);
for i = 1:ndims
    Bins{i} = ([1:msize(i)]'-1)./repmat(k(i)',msize(i),1) ... 
              +repmat(bound_lims(i,1),msize(i),1)+round(repmat(k(i)',msize(i),1).^-1/2);
end


% ROUND position
% REMOVE bins outside the computational region
Pos = round((pos-repmat(bound_lims(:,1)',size(pos,1),1)).*repmat(k',size(pos,1),1))+1;
for i = 1:ndims   
    Pos(Pos(:,i)<1|Pos(:,i)>Nbin(i)|~nniz(Pos),:) = [];
end

if ndims==1,
    amsize = [msize',1];
else
    amsize = msize';
end

% ACCUMULATE feature space occupancy in seconds
Occupancy = accumarray(Pos,1,amsize)./posSampleRate;

% ACCUMULATE spikes binned by feature space (e.g. xy coordinates)
if ~isempty(spkpos),
    spkpos = round((spkpos-repmat(bound_lims(:,1)',size(spkpos,1),1)).*repmat(k',size(spkpos,1),1))+1;
    for i = 1:ndims
        spkpos(spkpos(:,i)<1|spkpos(:,i)>Nbin(i)|~nniz(spkpos),:) = [];
    end
    SpikeCount = accumarray(spkpos,1,amsize);
else
    SpikeCount = zeros([amsize,1]);
end


% GENERATE gaussian smoothing kernel with std of # bins <- SmoothingWeight
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

% SMOOTH Occupancy
SOcc   = convn(Occupancy, Smoother,'same');
% SMOOTH Spike count
SCount = convn(SpikeCount,Smoother,'same');

% $$$ SOcc   = convn(Occupancy, Smoother,'same');
% $$$ SCount = convn(SpikeCount,Smoother,'same');


OccThresh = 4*(10/200).^numel(binDims);

%% Find the total occupancy and each pixels 
%% probability of occupancy
gtind = SOcc>OccThresh;
% TotOcc = sum(SOcc(gtind));
% POcc = SOcc./TotOcc;
TotOcc = sum(Occupancy(gtind));
POcc = Occupancy./TotOcc;


% COMPUTE RateMap
RateMap(gtind) = SCount(gtind)./SOcc(gtind);
RateMap = RateMap(:);
RateMap(~gtind) = nan;

% COMPUTE unit mean rate given state
MRate = sum(SpikeCount(gtind))/TotOcc;
% COMPUTE Spatial Information
if nargout >= 3,  SI = nansum(POcc(gtind).*(RateMap(gtind)./MRate).*log2(RateMap(gtind)./MRate));end
% COMPUTE Spatial Sparsity
if nargout >= 4,  Spar = 1/nansum(POcc(gtind).*RateMap(gtind).^2./MRate.^2);  end

% END MAIN -----------------------------------------------------------------------------------------

