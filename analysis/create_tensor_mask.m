function [mask] = create_tensor_mask(bins,varargin)


if ~isempty(varargin),
    boundary = varargin{1};
else
    boundary = struct('shape',  'circular',   ...
                      'radius', 440);
end

% $$$ defargs = struct('boundary',               struct('shape',  'circular',   ...
% $$$                                                   'radius', 440)        ...
% $$$ );
% $$$ [boundary] = DefaultArgs(varargin,defargs,'--struct');


switch boundary.shape
  case 'circular',
    width = numel(bins{1});    
    height = numel(bins{2});    
    radius = round(numel(bins{1})/2)-find(bins{1}<-boundary.radius,1,'last');
    centerW = width/2;
    centerH = height/2;
    [W,H] = meshgrid(1:width,1:height);
    mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);

  case 'circular_bhv'
    
% $$$     width = numel(bins{1});    
% $$$     height = numel(bins{2});    
% $$$     radius = round(numel(bins{1})/2)-find(bins{1}<-boundary.radius,1,'last');
% $$$     centerW = width/2;
% $$$     centerH = height/2;
% $$$     [W,H] = meshgrid(1:width,1:height);
% $$$     circMask = repmat(double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius),[1,1,pfsBinsDims(3:4)]);

% $$$         interpParPfb = struct('bins',{ {pitchSpaceDomain,...
% $$$                             pitchSpaceDomain}},...
% $$$                               'nanMaskThreshold', 0.01,...
% $$$                               'methodNanMap',     'linear',...
% $$$                               'methodRateMap',    'linear');
% $$$         bhvMask = zeros([cellfun(@numel,pfb.adata.bins)]);
% $$$         bhvMask(ds.vDims) = 1;
% $$$         interpGrids = cell([1,numel(interpParPfb.bins)]);
% $$$         [interpGrids{:}] = ndgrid(interpParPfb.bins{:});
% $$$         bhvMask = interpn(pfb.adata.bins{:},bhvMask,interpGrids{:},interpParPfb.methodRateMap);
% $$$         bhvMask(isnan(bhvMask)) = 0;
% $$$         sinterpGrids = interpGrids;
% $$$         SmoothingWeights = [0.2,0.2];
% $$$         for i = 1:ndims(sinterpGrids),
% $$$             sinterpGrids{i} = sinterpGrids{i}.^2/SmoothingWeights(i)^2/2;
% $$$         end
% $$$         Smoother = exp(sum(-cat(ndims(sinterpGrids)+1,sinterpGrids{:}),ndims(sinterpGrids)+1));
% $$$         Smoother = Smoother./sum(Smoother(:));
% $$$         bhvMask   = convn(bhvMask, Smoother,'same');
% $$$         bhvMask = double(bhvMask>0.01);
% $$$         if strcmp(mode,'xyhb')
% $$$             bhvMask = repmat(permute(bhvMask,[3,4,1,2]),[pfsBinsDims(1:2),1,1]);
% $$$         elseif strcmp(mode,'xyh')
% $$$             bhvMask = repmat(permute(~all(bhvMask==0,2),[3,4,1,2]),[pfsBinsDims(1:2),1,1]);
% $$$         elseif strcmp(mode,'xyb')            
% $$$             bhvMask = repmat(permute(~all(bhvMask==0,2),[3,4,1,2]),[pfsBinsDims(1:2),1,1]);
% $$$             %bhvMask = repmat(permute(~all(bhvMask==0),[3,4,2,1]),[pfsBinsDims(1:2),1,1]);
% $$$         end
    
% $$$     bhvMask = repmat(permute(bhvMask,[3,4,1,2]),[pfsBinsDims(1:2),1,1]);    
% $$$     mask = logical(circMask.*bhvMask);
    
  case 'HB'
    
    

end
