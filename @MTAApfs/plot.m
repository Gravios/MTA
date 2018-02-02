
function rateMap = plot(Pfs,varargin)
% function rateMap = plot(Pfs,varargin)
% plot placefield of the specified unit
% varargin:
%    unit                   [],                                                   
%    nMode                  'mean',                                               
%    ifColorbar             false,                                                
%    maxRate                [],                                                   
%    isCircular             true
%    reqSigVals             false
%


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('unit',                   [],                                                   ...
                 'nMode',                  'mean',                                               ...
                 'ifColorbar',             false,                                                ...
                 'maxRate',                [],                                                   ...
                 'isCircular',             true,                                                 ...
                 'sigThresh',              []                                                    ...
);
[unit,nMode,ifColorbar,maxRate,isCircular,sigThresh] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

if isempty(unit),unit=Pfs.data.clu(1);end
switch numel(Pfs.parameters.type)
  case 2
    bin1 = Pfs.adata.bins{1};
    bin2 = Pfs.adata.bins{2};

    
    mask = 1;
    if isCircular,
% RESTRICT rateMap area to within the maze radius
        width = Pfs.adata.binSizes(1);
        height = Pfs.adata.binSizes(2);
        radius = round(Pfs.adata.binSizes(1)/2)-find(Pfs.adata.bins{1}<-420,1,'last');
        centerW = width/2;
        centerH = height/2;
        [W,H] = meshgrid(1:width,1:height);           
        mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
        mask(mask==0)=nan;
    end

% COMPUTE rateMap output type
    if isnumeric(nMode)
        rateMap = Pfs.data.rateMap(:,Pfs.data.clu==unit,nMode);    
    else
        switch nMode
          case 'mean',    rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),3,'omitnan');
          case 'std' ,    rateMap = std(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),[],3,'omitnan');
          case 'snr',     rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),3,'omitnan')...
                                    ./std(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),[],3,'omitnan');
                          rateMap(isinf(rateMap)) = nan;                          
          case 'sig'
            rateMap = 1./sum((repmat(max(Pfs.data.rateMap(:,Pfs.data.clu==unit,:)),[size(Pfs.data.rateMap,1),1,1])...
                              -repmat(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),[1,1,Pfs.parameters.numIter]))<0,3)';
          otherwise
            rateMap = Pfs.data.rateMap(:,Pfs.data.clu==unit,1);
        end                    
    end
    
    if ~isempty(sigThresh) && Pfs.parameters.numIter > 1,
        sigMask = sum(~isnan(Pfs.data.rateMap(:,Pfs.data.clu==unit,:)),3)>Pfs.parameters.numIter*sigThresh;
        rateMap(~sigMask(:)) = nan;
    end 

    
% RESHAPE rateMap from 1D to 2D
    rateMap = reshape(rateMap,numel(bin1),numel(bin2)).*mask;
    
    if nargout>0,return,end

% ASSIGN maximum rate of colormap
    if isempty(maxRate), maxRate = max(rateMap(:)); end
    if isnan(maxRate), maxRate = 0; end    

% ASSIGN elements with nans a black color value
    rateMap(isnan(rateMap)) = -1;
    imagesc(bin1,bin2,rateMap');

% APPEND text which contains the maximum rate of the map
% $$$     text(Pfs.adata.bins{1}(end)-250,Pfs.adata.bins{2}(end)-50,...
% $$$          sprintf('%2.1f',max(rateMap(:))),'Color','w','FontWeight','bold','FontSize',8)
% SET colormap with lowest value as black
    colormap([0,0,0;parula]);
% ADJUST color scale
    caxis([-1,maxRate]);
    axis('xy');
  
  case 3
    c = eye(3);
    r = [1.2,3,6];
    rateMap = permute(reshape(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),Pfs.adata.binSizes'),[1,2,3]);
    if nargout>0,return,end
    hrate = max(Pfs.data.rateMap(:,Pfs.data.clu==unit,1))/2;
    [mind] = LocalMinimaN(-rateMap,-hrate,9);

    m = find(ismember('xyz',nMode));
    o = find(~ismember([1,2,3],m));
    ss = {};
    if ~isempty(mind),
        mind = mind(1,:);
        ss{m(1)} = ':';
        ss{m(2)} = ':';
        ss{o}    = mind(o);
        imagescnan({Pfs.adata.bins{m(1)},Pfs.adata.bins{m(2)},sq(nanmean(subsref(rateMap,substruct('()',ss)),o))'},colorLimits,[],ifColorbar,[0,0,0]);
        axis xy
    end

end

% END MAIN -----------------------------------------------------------------------------------------