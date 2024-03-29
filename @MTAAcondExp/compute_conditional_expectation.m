function [RateMap, Bins, SI, Spar] = compute_conditional_expectation(FetData,SpkData,varargin)

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('binDims',                                50,                                   ...
                 'SmoothingWeights',                       [],                                   ...
                 'type',                                   'xy',                                 ...
                 'bound_lims',                             [],                                   ...
);
[binDims,SmoothingWeights,type,bound_lims,posSampleRate] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

SI = [];
Spar = [];

% SET
boundaries = Data.model.domain;

ndims = numel(binDims);
Nbin = round(abs(diff(boundaries,1,2))./binDims');
if isempty(SmoothingWeights),  SmoothingWeights = Nbin./30;  end


% SET bins
Bins = cell(1,ndims);
k = Nbin./abs(diff(boundaries,1,2));
msize = round(abs(diff(boundaries,1,2)).*k);
for i = 1:ndims
    Bins{i} = ([1:msize(i)]'-1)./repmat(k(i)',msize(i),1)+repmat(boundaries(i,1),msize(i),1)+round(repmat(k(i)',msize(i),1).^-1/2);
end


% ROUND position
% REMOVE bins outside the computational region
Pos = round((pos-repmat(boundaries(:,1)',size(pos,1),1)).*repmat(k',size(pos,1),1))+1;
for i = 1:ndims   
    Pos(Pos(:,i)<1|Pos(:,i)>Nbin(i)|~nniz(Pos),:) = [];
end

if ndims==1,
    amsize = [msize',1];
else
    amsize = msize';
end

Occupancy = accumarray(Pos,1,amsize)./Data.sampleRate

if ~isempty(spkpos),
    spkpos = round((spkpos-repmat(boundaries(:,1)',size(spkpos,1),1)).*repmat(k',size(spkpos,1),1))+1;
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
if nargout >= 3,  SI = nansum(POcc(gtind).*(RateMap(gtind)./MRate).*log2(RateMap(gtind)./MRate));  end
if nargout >= 4,  Spar = 1/nansum(POcc(gtind).*RateMap(gtind).^2./MRate.^2);  end

% END MAIN -----------------------------------------------------------------------------------------

