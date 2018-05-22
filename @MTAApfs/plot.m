
function rateMap = plot(Pfs,varargin)
% function rateMap = plot(Pfs,varargin)
% plot placefield of the specified unit
% varargin:
%    unit - NumericArray: {[]}
%
%    mode - String: {x,'mean'} ploting method
%
%    rateReportMethod - String or Cellstr: {'','colorbar','text'}
%
%    maxRate - NumericArray: {[]}                                                   
%
%    mazeMaskFlag - Logical: {true} set bins to NAN based on mask
%
%    sigThresh - Numeric: {[]} proportion of non nan values accepted in the resampled estimates
%
%    flipAxesFlag - Logical: {false}
%
%    interpPar - Struct: {[]}
%
%        bins:             cellArray[1,2]( Numeric [N,1] ) 
%                Example: {linspace(-500,500,200)',linspace(-500,500,200)'}
%
%        nanMaskThreshold: Numeric,
%                Example: 0
%
%        methodNanMap: String
%                Example: 'linear'
%                Options: see interp2
%
%        methodRateMap: String,
%                Example: 'linear'
%                Options: see interp2
%
%    colorMap:     FuncHandle, @parula function handle for the colormap 




% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('unit',                   [],                                                   ...
                 'mode',                  'mean',                                                ...
                 'rateReportMethod',       '',                                                   ...
                 'maxRate',                [],                                                   ...
                 'mazeMaskFlag',           true,                                                 ...
                 'sigThresh',              [],                                                   ...
                 'flipAxesFlag',           false,                                                ...
                 'interpPar',              [],                                                   ...
                 'colorMap',               @parula,                                              ...
                 'mazeMask',               1,                                                    ...
                 'nanColor',               [0,0,0]                                               ...
);
[unit,mode,rateReportMethod,maxRate,mazeMaskFlag,sigThresh,flipAxesFlag,interpPar,colorMap,      ...
     mazeMask,nanColor] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

if isempty(unit),unit=Pfs.data.clu(1);end
switch numel(Pfs.parameters.type)
  case 2


% COMPUTE rateMap output type
    if isnumeric(mode)
        rateMap = Pfs.data.rateMap(:,Pfs.data.clu==unit,mode);    
    else
        if Pfs.parameters.numIter>1, 
            ind = 2:Pfs.parameters.numIter;
        else,
            ind = 1;
        end
        switch mode            
          case 'mean',    rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind),3,'omitnan');
          case 'std' ,    rateMap = std(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind),[],3,'omitnan');
          case 'snr',     rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind),3,'omitnan')...
                                    ./std(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind),[],3,'omitnan');
                          rateMap(isinf(rateMap)) = nan;                          
          case 'snrs',    
            tStd = std(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind),[],3,'omitnan');
            tStd(tStd<1) = 1;
            rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind),3,'omitnan')./tStd;
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
    rateMap = reshape(rateMap,Pfs.adata.binSizes(1),Pfs.adata.binSizes(2));

    if ~isempty(interpPar),
% INTERPOLATE ratemap 
        interpGrids= cell([1,numel(interpPar.bins)]);
        [interpGrids{:}] = ndgrid(interpPar.bins{:});

        nanMask = double(isnan(rateMap));
% Shold I use mean rate again for zeros interpolation?        
        rateMap(isnan(rateMap)) = 0;


        rateMap        = interp2(Pfs.adata.bins{:},rateMap',interpGrids{:},interpPar.methodRateMap);
        interpdNanMask = interp2(Pfs.adata.bins{:},nanMask',interpGrids{:},interpPar.methodNanMap);

% SMOOTH edges with interpolated nan mask
        rateMap(interpdNanMask>interpPar.nanMaskThreshold) = nan;
% CORRECT for cubic undershoot 
        rateMap(rateMap<0) = 0;

        bins = interpPar.bins;
        binSizes = cellfun(@numel,interpPar.bins);
    else
        bins = Pfs.adata.bins;
        binSizes = Pfs.adata.binSizes;
    end
    


    if mazeMaskFlag,
% RESTRICT rateMap area to within the maze radius
        if mazeMask==1,
            width = binSizes(1);
            height =binSizes(2);
            radius = round(binSizes(1)/2)-find(bins{1}<-450,1,'last');
            centerW = width/2;
            centerH = height/2;
            [W,H] = meshgrid(1:width,1:height);           
            mazeMask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
            mazeMask(mazeMask==0)=nan;
        end
    end
    rateMap = rateMap.*mazeMask;
    
    
    
    if nargout>0,return,end

% ASSIGN maximum rate of colormap
    if isempty(maxRate), maxRate = max(rateMap(:)); end
    if isnan(maxRate), maxRate = 0; end    
    if numel(maxRate)<2, maxRate = [0,maxRate]; end

% ASSIGN elements with nans a black color value
    %rateMap(isnan(rateMap)) = -1;
    if flipAxesFlag,
        bins = fliplr(bins(:)');
        imagescnan({bins{:},fliplr(rot90(rateMap',-1))},maxRate,'linear',false,nanColor,[],[],colorMap);
    else
        imagescnan({bins{:},rateMap'},maxRate,'linear',false,nanColor,[],[],colorMap);
    end
    
    if islogical(rateReportMethod) && rateReportMethod,
        rateReportMethod = 'colorbar';
    elseif islogical(rateReportMethod) && ~rateReportMethod,
        rateReportMethod = '';
    end
    
    if ~iscell(rateReportMethod)
        rateReportMethod = {rateReportMethod};
    end
    
    for method = rateReportMethod,
        switch method{1}
          case 'colorbar'
            colorbar();
            colormap(func2str(colorMap));
            caxis([maxRate]);
          case 'text'
            text(Pfs.adata.bins{1}(end)-0.25*diff(Pfs.adata.bins{1}([1,end])),...
                 Pfs.adata.bins{2}(end)-0.10*diff(Pfs.adata.bins{2}([1,end])),...
                 sprintf('%2.1f',max(rateMap(:))),...
                 'Color','w','FontWeight','bold','FontSize',8)
        end
    end

    axis('xy');
  
% $$$   case 3
% $$$     c = eye(3);
% $$$     r = [1.2,3,6];
% $$$     rateMap = permute(reshape(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),Pfs.adata.binSizes'),[1,2,3]);
% $$$     if nargout>0,return,end
% $$$     hrate = max(Pfs.data.rateMap(:,Pfs.data.clu==unit,1))/2;
% $$$     [mind] = LocalMinimaN(-rateMap,-hrate,9);
% $$$ 
% $$$     m = find(ismember('xyz',mode));
% $$$     o = find(~ismember([1,2,3],m));
% $$$     ss = {};
% $$$     if ~isempty(mind),
% $$$         mind = mind(1,:);
% $$$         ss{m(1)} = ':';
% $$$         ss{m(2)} = ':';
% $$$         ss{o}    = mind(o);
% $$$         imagescnan({Pfs.adata.bins{m(1)},                                ... xaxis labels
% $$$                     Pfs.adata.bins{m(2)},                                ... yaxis labels
% $$$                     sq(nanmean(subsref(rateMap,substruct('()',ss)),o))'},... imagedata
% $$$                    colorLimits,                                          ... color scale
% $$$                    [],                                                   ... color scale mode
% $$$                    colorbarFlag,                                         ... flag to add colorbar
% $$$                    [0,0,0]                                               ... Nan color value
% $$$         );
% $$$         axis xy
% $$$     end

end

% END MAIN -----------------------------------------------------------------------------------------