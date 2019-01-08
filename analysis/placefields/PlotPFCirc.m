function [ratemap, Bins, spatialInformation, sparsity] = ...
    PlotPFCirc(Session,spkpos,pos,binDims,SmoothingWeights,type,bound_lims,posSampleRate)

% $$$ % DEFARGS ------------------------------------------------------------------------------------------
% $$$ defargs = struct('binDims',                                50,                                   ...
% $$$                  'SmoothingWeights',                       [],                                   ...
% $$$                  'type',                                   'xy',                                 ...
% $$$                  'bound_lims',                             [],                                   ...
% $$$                  'posSampleRate',                          Session.xyz.sampleRate                ...
% $$$ );
% $$$ [binDims,SmoothingWeights,type,bound_lims,posSampleRate] = DefaultArgs(varargin,defargs,'--struct');
% $$$ %---------------------------------------------------------------------------------------------------





% MAIN ---------------------------------------------------------------------------------------------

spatialInformation = [];
sparsity = [];

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
    Bins{i} = ([1:msize(i)]'-1)./repmat(k(i)',msize(i),1)+...
               repmat(bound_lims(i,1),msize(i),1)+...
               repmat(k(i)',msize(i),1).^-1/2;
    %Bins{i} = linspace(bound_lims(i,1),bound_lims(i,2),Nbin(i));
end


% ROUND position
% REMOVE bins outside the computational region
pos = round((pos-repmat(bound_lims(:,1)',size(pos,1),1)).*repmat(k',size(pos,1),1))+1;
for i = 1:ndims   
    pos(pos(:,i)<1|pos(:,i)>Nbin(i)|~nniz(pos),:) = [];
end

if ndims==1,
    amsize = [msize',1];
else
    amsize = msize';
end

% ACCUMULATE feature space occupancy in seconds
occupancy = accumarray(pos,1,amsize)./posSampleRate;

% ACCUMULATE spikes binned by feature space (e.g. xy coordinates)
if ~isempty(spkpos),
    spkpos = round((spkpos-repmat(bound_lims(:,1)',size(spkpos,1),1)).*repmat(k',size(spkpos,1),1))+1;
    for i = 1:ndims
        spkpos(spkpos(:,i)<1|spkpos(:,i)>Nbin(i)|~nniz(spkpos),:) = [];
    end
    spikeCount = accumarray(spkpos,1,amsize);
else
    spikeCount = zeros([amsize,1]);
end


 % GENERATE gaussian smoothing kernel with std of # bins <- SmoothingWeight
sind = cell(1,ndims);
for i = 1:ndims,
    sind{i} = linspace(-round(msize(i)/2),round(msize(i)/2),msize(i)+double(mod(msize(i),2)==0));
end
[sind{:}] = ndgrid(sind{:});
for i = 1:ndims,
    sind{i} = sind{i}.^2/SmoothingWeights(i)^2/2;
end
Smoother = exp(sum(-cat(ndims+1,sind{:}),ndims+1));
Smoother = Smoother./sum(Smoother(:));

% SMOOTH occupancy
% SMOOTH Spike count
ROcc   = repmat(occupancy,[1,3,1]);
SOcc   = convn(ROcc, Smoother,'same');
SOcc   = SOcc(:,Nbin(2)+1:Nbin(2)*2,:,:,:,:);
RCount = repmat(spikeCount,[1,3,1]);
SCount = convn(RCount,Smoother,'same');
SCount = SCount(:,Nbin(2)+1:Nbin(2)*2,:,:,:,:);

% CREATE minimally smoothed occupancy map
soc = occupancy;
soc(isnan(soc)) = 0;

if ndims==1,
    soc = RectFilter(soc,3,1);
elseif ndims>=2,    
    dshift = circshift(1:ndims,-1);
    for t = 1:2,
    for s =1:ndims    
        soc = RectFilter(permute(soc,dshift),3,1);
    end
    end
end

%OccThresh = 4*(mean(binDims)/200).^numel(binDims);

%gtind = SOcc > 2*(10/100).^numel(binDims)*numel(Bins{2});
gtind = soc > 2*(1/15).^numel(binDims)*numel(Bins{2});



%% Find the total occupancy and each pixels 
%% probability of occupancy

% totalOcc = sum(SOcc(gtind));
% pOcc = SOcc./totalOcc;
totalOcc = sum(occupancy(gtind));
pOcc = occupancy./totalOcc;


% COMPUTE ratemap from smoothed spike and occupancy hists
ratemap = SCount./SOcc;
ratemap(~gtind) = nan;
ratemap = ratemap(:);

% COMPUTE unit mean rate given state
MRate = sum(spikeCount(gtind))/totalOcc;

% COMPUTE smoothed rate map 
% $$$ srmap = RCount./ROcc;           
% $$$ srmap(~gtind) = nan;
% $$$ %srmap(~gtind) = nan;
% $$$ %srmap(isnan(srmap)) = 10.^mean(log10(ratemap(:)),'omitnan');
% $$$ srmap(isnan(srmap)) = 0
% $$$ srmap = convn(srmap,Smoother,'same');
% $$$ srmap = srmap(:,Nbin(2)+1:Nbin(2)*2,:,:,:,:);
% $$$ 
% $$$ % $$$ % AVERAGE 
% $$$ ratemap = reshape(mean(cat(ndims+1,srmap,ratemap),ndims+1),[],1);

if nargout >= 3,  
% COMPUTE Spatial Information  
    spatialInformation = nansum(1./sum(double(gtind(:))).*(ratemap(gtind)./MRate)...
                                .*log2(ratemap(gtind)./MRate));
    if nargout == 4,
% COMPUTE Sparsity
        sparsity = nansum(1./sum(double(gtind(:))).*ratemap(gtind)).^2 ...
            ./nansum(1./sum(double(gtind(:))).*ratemap(gtind).^2);
    end    
end



% END MAIN -----------------------------------------------------------------------------------------

