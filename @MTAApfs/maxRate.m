function [mxr,mxp] = maxRate(Pfs,varargin)
% function [mxr,mxp] = maxRate(Pfs,varargin)
% return the maximum firing rate and location of a unit(s)
%
%  VARARGIN:
%    units                   []
%    isCircular              true
%    mode                    'mean'
%


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',                  [],                                                   ...
                 'mazeMaskFlag',           true,                                                 ...
                 'mode',                   'mean',                                               ...
                 'sigThresh',              [],                                                   ...
                 'interpPar',              struct('bins',{{linspace(-500,500,200)',              ...
                                                           linspace(-500,500,200)'}},            ...
                                                  'nanMaskThreshold', 0,                         ...
                                                  'methodNanMap',     'linear',                  ...
                                                  'methodRateMap',    'linear')                  ...
);
[units,mazeMaskFlag,mode,sigThresh,interPar] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------
if isempty(units),
    units = sort(Pfs.data.clu');
end
if mazeMaskFlag,
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

        sigMask = [];
        if ~isempty(sigThresh) && Pfs.parameters.numIter > 1,
            sigMask = sum(~isnan(Pfs.data.rateMap(:,Pfs.data.clu==u,:)),3)>Pfs.parameters.numIter*sigThresh;
        end 
        
        rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==u,:),3,'omitnan');
        
        rateMap(~sigMask(:)) = nan;
            
        if size(rateMap,2)==0,  rateMap = nan([size(Pfs.data.rateMap,1),1]);  end
        [mxr(u==units),mxp(u==units)] = max(rateMap.*mask);
        
      case 'std'        
        if ~isempty(sigThresh) && Pfs.parameters.numIter > 1,
            sigMask = sum(~isnan(Pfs.data.rateMap(:,Pfs.data.clu==u,:)),3)>Pfs.parameters.numIter*sigThresh;
        end         
        
        rateMap = std(Pfs.data.rateMap(:,Pfs.data.clu==u,:),[],3,'omitnan');
        
        rateMap(~sigMask(:)) = nan;
        
        if size(rateMap,2)==0,  rateMap = nan([size(Pfs.data.rateMap,1),1]);  end
        [mxr(u==units),mxp(u==units)] = max(rateMap.*mask);
        
      case 'prctile99'
        rateMap = Pfs.data.rateMap(:,Pfs.data.clu==u,:);
        if ~isempty(sigThresh) && Pfs.parameters.numIter > 1,
            sigMask = sum(~isnan(Pfs.data.rateMap(:,Pfs.data.clu==u,:)),3)>Pfs.parameters.numIter*sigThresh;
            rateMap(~sigMask(:)) = nan;
        end 
        
        mxr(u==units) = prctile(rateMap(nniz(rateMap)),99.9);
        return
        
      otherwise
        rateMap = Pfs.data.rateMap(:,Pfs.data.clu==u,1);
        if size(rateMap,2)==0,  rateMap = nan([size(Pfs.data.rateMap,1),1]);  end
        [mxr(u==units),mxp(u==units)] = max(rateMap.*mask);
        
    end
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