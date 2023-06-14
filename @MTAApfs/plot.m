function [ratemap] = plot(Pfs,varargin)
% function ratemap = plot(Pfs,varargin)
% plot placefield of the specified unit
%
% IN:
%    Pfs - MTAApfs: ratemap analysis object
% 
% VARARGIN:
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
%    colorMap - FuncHandle: {@parula} function handle for the colormap 
%
%    mazeMask - NumericArray: {1} mask to nullifiy elements out of bounds
%
%    nanColor - NumericArray: {[0,0,0]} normalized rgb values for nan elements
%
% OUT:
%    ratemap - NumericMatrix: matrix containing the ratemap of the requested unit




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

if isnumeric(mode)
    ratemap = Pfs.data.rateMap(:,Pfs.data.clu==unit,mode);    
else
    if Pfs.parameters.numIter>1, 
        ind = 2:Pfs.parameters.numIter;
    else,
        ind = 1;
    end
    
    % yep.... i'm going there
    if ischar(mode),
        switch mode            
          case 'mean',    ratemap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind),3,'omitnan');
          case 'mean_circ_smooth', 
            ratemap = mean(RectFilter(sq(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind)),3,1,'circular'),2,'omitnan');
          case 'all_circ_smooth', 
            ratemap = RectFilter(sq(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind)),3,1,'circular');
          case 'std' ,    ratemap = std(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind),[],3,'omitnan');
          case 'snr',     ratemap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind),3,'omitnan')...
                ./std(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind),[],3,'omitnan');
            ratemap(isinf(ratemap)) = nan;                          
          case 'snrs',    
            tStd = std(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind),[],3,'omitnan');
            tStd(tStd<1) = 1;
            ratemap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,ind),3,'omitnan')./tStd;
            ratemap(isinf(ratemap)) = nan;                          
          case 'sig'
            ratemap = 1./sum((repmat(max(Pfs.data.rateMap(:,Pfs.data.clu==unit,:)),[size(Pfs.data.rateMap,1),1,1])...
                              -repmat(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),[1,1,Pfs.parameters.numIter]))<0,3)';
          otherwise
            ratemap = Pfs.data.rateMap(:,Pfs.data.clu==unit,1);
        end                    
    elseif iscell(mode),
        ratemap = Pfs.data.rateMap(:,Pfs.data.clu==unit,mode{end});
    else
        ratemap = Pfs.data.rateMap(:,Pfs.data.clu==unit,1);
    end
    
end

% COMPUTE ratemap output type    
if ~isempty(sigThresh) && Pfs.parameters.numIter > 1,
    sigMask = sum(~isnan(Pfs.data.rateMap(:,Pfs.data.clu==unit,:)),3)>Pfs.parameters.numIter*sigThresh;
    ratemap(~sigMask(:)) = nan;
end 

% RESHAPE ratemap from 1D to ND
if numel(Pfs.adata.binSizes) > 1,
    ratemap = reshape(ratemap,Pfs.adata.binSizes');    
    switch numel(Pfs.adata.binSizes)
      case 2

        

        if ~isempty(interpPar),
            % INTERPOLATE ratemap 
            interpGrids= cell([1,numel(interpPar.bins)]);
            [interpGrids{:}] = ndgrid(interpPar.bins{:});

            nanMask = double(isnan(ratemap));
            % Shold I use mean rate again for zeros interpolation?        
            ratemap(isnan(ratemap)) = 0;


            ratemap        = interp2(Pfs.adata.bins{:},ratemap',interpGrids{:},interpPar.methodRateMap);
            interpdNanMask = interp2(Pfs.adata.bins{:},nanMask',interpGrids{:},interpPar.methodNanMap);

            if numel(mazeMask)>1,
                mazeMask = double(interp2(Pfs.adata.bins{:},...
                                          double(mazeMask)',...
                                          interpGrids{:},...
                                          interpPar.methodRateMap)>0.5);
            end
            % SMOOTH edges with interpolated nan mask
            ratemap(interpdNanMask>interpPar.nanMaskThreshold) = nan;
            % CORRECT for cubic undershoot 
            ratemap(ratemap<0) = 0;

            bins = interpPar.bins;
            binSizes = cellfun(@numel,interpPar.bins);
        else
            bins = Pfs.adata.bins;
            binSizes = Pfs.adata.binSizes;
        end
        
        if mazeMaskFlag,
            % RESTRICT ratemap area to within the maze radius
            if numel(mazeMask)==1,
                width = binSizes(1);
                height =binSizes(2);
                radius = round(binSizes(1)/2)-find(bins{1}<=-450,1,'last');
                centerW = width/2;
                centerH = height/2;
                [W,H] = meshgrid(1:width,1:height);           
                mazeMask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
            end
        end
        if islogical(mazeMask),
            mazeMask = double(mazeMask);
        end
        mazeMask(mazeMask==0)=nan;

        if ~isempty(mazeMask)
            ratemap = ratemap.*mazeMask;
        end
        
        
        
        
        if nargout>0,return,end

        % ASSIGN maximum rate of colormap
        if isempty(maxRate), maxRate = max(ratemap(:)); end
        if isnan(maxRate), maxRate = 0; end    
        if numel(maxRate)<2, maxRate = [0,maxRate]; end

        % ASSIGN elements with nans a black color value
        %ratemap(isnan(ratemap)) = -1;
        if flipAxesFlag,
            bins = fliplr(bins(:)');
            imagescnan({bins{:},fliplr(rot90(ratemap',-1))},maxRate,'linear',false,nanColor,[],[],colorMap);
        else
            imagescnan({bins{:},ratemap'},maxRate,'linear',false,nanColor,[],[],colorMap);
        end
        
        if islogical(rateReportMethod) && rateReportMethod,
            rateReportMethod = 'colorbar';
        elseif islogical(rateReportMethod) && ~rateReportMethod,
            rateReportMethod = '';
        end
        

        switch rateReportMethod,
          case 'colorbar'
            colorbar();
            colormap(gca(),func2str(colorMap));
            caxis([maxRate]);
          case 'text'
            mrate = max(ratemap(:));
            hax = text;
            
            if mrate < 10
                hax.String   = sprintf('%2.1f ',max(ratemap(:)));                    
            else
                hax.String   = sprintf('%2.0f ',max(ratemap(:)));
            end
            %hax.Position = [Pfs.adata.bins{1}(end)-0.45*diff(Pfs.adata.bins{1}([1,end])),...
            %                Pfs.adata.bins{2}(end)-0.10*diff(Pfs.adata.bins{2}([1,end]))];
            hax.Position(1:2) = [Pfs.adata.bins{1}(end) - hax.Extent(3),...
                                 Pfs.adata.bins{2}(end) - hax.Extent(4)];
            hax.Color    = 'w';
            hax.FontWeight = 'bold';
            hax.FontSize = 8;
        end


        axis('xy');
        
      case 3
        if iscell(mode),
            mode = mode(1:end-1);
            ratemap = sq(ratemap(mode{:}));
            primaryDims = find(cellfun(@ischar,mode));
            
            if ~isempty(interpPar),
                % INTERPOLATE ratemap 
                interpGrids= cell([1,numel(interpPar.bins)]);
                [interpGrids{:}] = ndgrid(interpPar.bins{:});

                nanMask = double(isnan(ratemap));
                % Shold I use mean rate again for zeros interpolation?        
                ratemap(isnan(ratemap)) = 0;


                ratemap        = interp2(Pfs.adata.bins{:},ratemap',interpGrids{:},interpPar.methodRateMap);
                interpdNanMask = interp2(Pfs.adata.bins{:},nanMask',interpGrids{:},interpPar.methodNanMap);

                if numel(mazeMask)>1,
                    mazeMask = double(interp2(Pfs.adata.bins{:},...
                                              double(mazeMask)',...
                                              interpGrids{:},...
                                              interpPar.methodRateMap)>0.5);
                end
                % SMOOTH edges with interpolated nan mask
                ratemap(interpdNanMask>interpPar.nanMaskThreshold) = nan;
                % CORRECT for cubic undershoot 
                ratemap(ratemap<0) = 0;

                bins = interpPar.bins;
                binSizes = cellfun(@numel,interpPar.bins);
            else
                bins = Pfs.adata.bins;
                binSizes = Pfs.adata.binSizes;
            end
            
            if mazeMaskFlag,
                % RESTRICT ratemap area to within the maze radius
                if numel(mazeMask)==1,
                    width = binSizes(1);
                    height =binSizes(2);
                    radius = round(binSizes(1)/2)-find(bins{1}<-450,1,'last');
                    centerW = width/2;
                    centerH = height/2;
                    [W,H] = meshgrid(1:width,1:height);           
                    mazeMask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
                end
            end
            if islogical(mazeMask),
                mazeMask = double(mazeMask);
            end
            mazeMask(mazeMask==0)=nan;

            if ~isempty(mazeMask)
                ratemap = ratemap.*mazeMask;
            end
            
            
            
            
            if nargout>0,return,end

            % ASSIGN maximum rate of colormap
            if isempty(maxRate), maxRate = max(ratemap(:)); end
            if isnan(maxRate), maxRate = 0; end    
            if numel(maxRate)<2, maxRate = [0,maxRate]; end

            % ASSIGN elements with nans a black color value
            %ratemap(isnan(ratemap)) = -1;
            if flipAxesFlag,
                bins = fliplr(bins(:)');
                imagescnan({bins{primaryDims},...
                            fliplr(rot90(ratemap',-1))},...
                           maxRate,'linear',false,nanColor,[],[],colorMap);
            else
                imagescnan({bins{primaryDims},ratemap'},maxRate,'linear',false,nanColor,[],[],colorMap);
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
                    mrate = max(ratemap(:));
                    if mrate < 10
                        tax = text(Pfs.adata.bins{1}(end)-0.45*diff(Pfs.adata.bins{1}([1,end])),...
                             Pfs.adata.bins{2}(end)-0.10*diff(Pfs.adata.bins{2}([1,end])),...
                             sprintf('%2.1f',max(ratemap(:))),...
                             'Color','w','FontWeight','bold','FontSize',8)
                    else
                        tax = text(Pfs.adata.bins{1}(end)-0.45*diff(Pfs.adata.bins{1}([1,end])),...     
                             Pfs.adata.bins{2}(end)-0.10*diff(Pfs.adata.bins{2}([1,end])),...
                             sprintf('%2.0f',max(ratemap(:))),...
                             'Color','w','FontWeight','bold','FontSize',8)
                    end
                end
            end

            axis('xy');
        else,
            if ~isempty(interpPar),
                % INTERPOLATE ratemap 
                interpGrids= cell([1,numel(interpPar.bins)]);
                % WARNING!!! Double check for non symetric placefields due to the incongruence between
                %            meshgrid and ndgrid.
                [interpGrids{:}] = ndgrid(interpPar.bins{:});

                nanMask = double(isnan(ratemap));
                % Shold I use mean rate again for zeros interpolation?        
                ratemap(isnan(ratemap)) = 0;

                ratemap        = interpn(Pfs.adata.bins{:},ratemap,interpGrids{:},interpPar.methodRateMap);
                interpdNanMask = interpn(Pfs.adata.bins{:},nanMask,interpGrids{:},interpPar.methodNanMap);

                % SMOOTH edges with interpolated nan mask
                ratemap(interpdNanMask>interpPar.nanMaskThreshold) = nan;
                % CORRECT for cubic undershoot 
                ratemap(ratemap<0) = 0;
            end
            
            ratemap = ratemap.*mazeMask;
        end
        
% $$$     c = eye(3);
% $$$     r = [1.2,3,6];
% $$$     ratemap = permute(reshape(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),Pfs.adata.binSizes'),[1,2,3]);
% $$$     if nargout>0,return,end
% $$$     hrate = max(Pfs.data.rateMap(:,Pfs.data.clu==unit,1))/2;
% $$$     [mind] = LocalMinimaN(-ratemap,-hrate,9);
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
% $$$                     sq(nanmean(subsref(ratemap,substruct('()',ss)),o))'},... imagedata
% $$$                    colorLimits,                                          ... color scale
% $$$                    [],                                                   ... color scale mode
% $$$                    colorbarFlag,                                         ... flag to add colorbar
% $$$                    [0,0,0]                                               ... Nan color value
% $$$         );
% $$$         axis xy
% $$$     end
      case 4

        if ~isempty(interpPar),
            % INTERPOLATE ratemap 
            interpGrids= cell([1,numel(interpPar.bins)]);
            % WARNING!!! Double check for non symetric placefields due to the incongruence between
            %            meshgrid and ndgrid.
            [interpGrids{:}] = ndgrid(interpPar.bins{:});

            nanMask = double(isnan(ratemap));
            % Shold I use mean rate again for zeros interpolation?        
            ratemap(isnan(ratemap)) = 0;

            ratemap        = interpn(Pfs.adata.bins{:},ratemap,interpGrids{:},interpPar.methodRateMap);
            interpdNanMask = interpn(Pfs.adata.bins{:},nanMask,interpGrids{:},interpPar.methodNanMap);

            % SMOOTH edges with interpolated nan mask
            ratemap(interpdNanMask>interpPar.nanMaskThreshold) = nan;
            % CORRECT for cubic undershoot 
            ratemap(ratemap<0) = 0;
        end

        ratemap = ratemap.*mazeMask;
        
    end% switch
end% if

% END MAIN -----------------------------------------------------------------------------------------



