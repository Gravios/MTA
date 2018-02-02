function [mxr,mxp] = maxRate(Pfs,varargin)
% function [mxr,mxp] = maxRate(Pfs,varargin)
% return the maximum firing rate and location of a unit(s)


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',                  [],                                                   ...
                 'isCircular',             true,                                                 ...
                 'mode',                   'mean'                                                ...
);
[units,isCircular,mode] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------
if isempty(units),
    units = Pfs.data.clu;
end
if isCircular,
    width = Pfs.adata.binSizes(1);
    height = Pfs.adata.binSizes(2);
    radius = round(Pfs.adata.binSizes(1)/2)-find(Pfs.adata.bins{1}<-420,1,'last');
    centerW = width/2;
    centerH = height/2;
    [W,H] = meshgrid(1:width,1:height);           
    mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
    mask(mask==0)=nan;
    if numel(Pfs.parameters.type)>2,
        mask = repmat(mask,[1,1,Pfs.adata.binSizes(3)]);
    end
    mask = reshape(mask,[],1);
else
    mask = 1;
end

mxr = nan(numel(units),1);
mxp = nan(numel(units),1);

for u = units(:)',
    switch mode
      case 'mean'        
        rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==u,:),3,'omitnan');
      case 'std'        
        rateMap = std(Pfs.data.rateMap(:,Pfs.data.clu==u,:),[],3,'omitnan');
      otherwise
        rateMap = Pfs.data.rateMap(:,Pfs.data.clu==u,1);
    end
    if size(rateMap,2)==0,
        rateMap = nan([size(Pfs.data.rateMap,1),1]);
    end

    [mxr(u==units),mxp(u==units)] = max(rateMap.*mask);
end
mxp = Ind2Sub(Pfs.adata.binSizes',mxp);
if numel(Pfs.parameters.type)>2,
    mxp = [Pfs.adata.bins{1}(mxp(:,1)), ...
           Pfs.adata.bins{2}(mxp(:,2)), ...
           Pfs.adata.bins{3}(mxp(:,3))];                                
else
    mxp = [Pfs.adata.bins{1}(mxp(:,1)), ...
           Pfs.adata.bins{2}(mxp(:,2))];
end

% END MAIN -----------------------------------------------------------------------------------------