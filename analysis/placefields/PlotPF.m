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


% SET bins
Bins = cell(1,ndims);
k = Nbin./abs(diff(bound_lims,1,2));
msize = round(abs(diff(bound_lims,1,2)).*k);
for i = 1:ndims
    Bins{i} = ([1:msize(i)]'-1)./repmat(k(i)',msize(i),1)+repmat(bound_lims(i,1),msize(i),1)+round(repmat(k(i)',msize(i),1).^-1/2);
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

Occupancy = accumarray(Pos,1,amsize)./posSampleRate;

if ~isempty(spkpos),
    spkpos = round((spkpos-repmat(bound_lims(:,1)',size(spkpos,1),1)).*repmat(k',size(spkpos,1),1))+1;
    for i = 1:ndims
        spkpos(spkpos(:,i)<1|spkpos(:,i)>Nbin(i)|~nniz(spkpos),:) = [];
    end
    SpikeCount = accumarray(spkpos,1,amsize);
else
    SpikeCount = zeros([amsize,1]);
end


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
OccThresh = 0.1.^numel(binDims);%0.06;0.12;%
                 %OccThresh = .03;%0.06;0.12;%
%% Find the total occupancy and each pixels 
%% probability of occupancy
gtind = SOcc>OccThresh;
% TotOcc = sum(SOcc(gtind));
% POcc = SOcc./TotOcc;
TotOcc = sum(Occupancy(gtind));
POcc = Occupancy./TotOcc;


%% Rate Map
%RateMap = NaN(prod(Nbin),1);
%RateMap(gtind) = SCount(gtind)./SOcc(gtind);
%RateMap = RateMap(:);
%RateMap(~gtind) = nan;
RateMap = SpikeCount./Occupancy;
RateMap(isnan(RateMap)) = 0;
RateMap = convn(RateMap,Smoother,'same');
%RateMap = nanconv(RateMap,Smoother);
RateMap = RateMap(:);
RateMap(~gtind) = nan;
%% Find the units overall mean rate given the 
%% current state
MRate = sum(SpikeCount(gtind))/TotOcc;
if nargout >= 3,  SI = nansum(POcc(gtind).*(RateMap(gtind)./MRate).*log2(RateMap(gtind)./MRate));  end
if nargout >= 4,  Spar = 1/nansum(POcc(gtind).*RateMap(gtind).^2./MRate.^2);  end

% END MAIN -----------------------------------------------------------------------------------------

