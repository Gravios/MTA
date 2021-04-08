function [mxr,mxp,fo] = get_max_rate(Pfs,varargin)
% function [mxr,mxp] = maxRate(Pfs,varargin)
% return the maximum firing rate and location of a unit(s)
%
%  IN:
%
%    units - NumericArray: { [] } list of unit clu numbers associated with the rate maps in Pfs
%    mazeMaskFlag - Logical: { true } Apply 
%    mode - String: { 'mean' }                                               ...
%                 'sigThresh',              [],                                                   ...
%                 'interpPar',              [],                                                   ...
%                 'fitStartPoint',          [],                                                   ...
%                 'mask',                   1                                                     ...
%
% interpolaiton only makes sense with extra filtering
%
% interpPar - struct: ('bins',{{linspace(-500,500,200)',              ...
%                                                           linspace(-500,500,200)'}},            ...
%                                                  'nanMaskThreshold', 0,                         ...
%                                                  'methodNanMap',     'linear',                  ...
%                                                  'methodRateMap',    'linear')
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',                  [],                                                   ...
                 'mazeMaskFlag',           true,                                                 ...
                 'mode',                   'mean',                                               ...
                 'sigThresh',              [],                                                   ...
                 'interpPar',              [],                                                   ...
                 'fitStartPoint',          [],                                                   ...
                 'mask',                   1                                                     ...
);
[units,mazeMaskFlag,mode,sigThresh,interpPar,fitStartPoint,mask] =                               ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

if isempty(interpPar),
    bins = Pfs.adata.bins;
else
    bins = interpPar.bins;
end

if isempty(units),
    units = sort(Pfs.data.clu');
end
if mazeMaskFlag,
    width = numel(bins{1});
    height = numel(bins{2});
    radius = round(numel(bins{1})/2)-find(bins{1}<-420,1,'last');
    centerW = width/2;
    centerH = height/2;
    [W,H] = meshgrid(1:width,1:height);           
    mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
    mask(mask==0)=nan;
    if numel(Pfs.parameters.type)>2,
        mask = repmat(mask,[1,1,numel(bins{3})]);
    end
    mask = reshape(mask,[],1);
end

mxr = nan(numel(units),1);
mxp = nan(numel(units),1);
fo  = cell(numel(units),1);
for u = units(:)',
    switch mode
% $$$       case 'mean-gaussfit'
% $$$         switch numel(Pfs.adata.bins),
% $$$           case 1
% $$$             g = fittype( @(A,xa,xo,x) A.*exp(-(xa.*(x-xo).^2)),'independent',{'x', 'y'},'dependent', 'z' ); 
% $$$             fitStartPoint = [8,0.00001,0];
% $$$           case 2
% $$$             g = fittype( @(A,xa,ya,xya,xo,yo,x,y) A.*exp(-(xa.*(x-xo).^2+xya.*(x-xo).*(y-yo)+ya.*(y-yo).^2)),...
% $$$                        'independent',{'x', 'y'},'dependent', 'z' ); 
% $$$             fitStartPoint = [100,0.00001,0.0000001,0.00001,0,0];
% $$$           case 3                        
% $$$             g = fittype( @(A,xa,ya,za,xya,xza,zya,xyza,xo,yo,zo,x,y,z) ...
% $$$                          A.*exp(-(xa.*(x-xo).^2+xya.*(x-xo).*(y-yo)+xza.*(x-xo).*(z-zo)+zya.*(z-zo).*(y-yo) ...
% $$$                                   +xyza.*(x-xo).*(z-zo).*(y-yo)+ya.*(y-yo).^2)),...
% $$$                        'independent',{'x', 'y','z'},'dependent', 'f'); 
% $$$         end
% $$$         
% $$$         [mxp,mxr] = LocalMinimaN();
% $$$ 
% $$$         sigMask = [];
% $$$         if ~isempty(sigThresh) && Pfs.parameters.numIter > 1,
% $$$             sigMask = sum(~isnan(Pfs.data.rateMap(:,Pfs.data.clu==u,:)),3)>Pfs.parameters.numIter*sigThresh;
% $$$         end 
% $$$ 
% $$$         %rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==u,:),3,'omitnan');
% $$$         rateMap = Pfs.plot(u,'mean',false,[],false,0.99,false,interpPar);
% $$$ 
% $$$         rateMap = RectFilter(rateMap');
% $$$         rateMap = RectFilter(rateMap');        
% $$$         rateMap = rateMap(:);
% $$$         
% $$$         rateMap(~sigMask(:)) = nan;
% $$$             
% $$$         if size(rateMap,2)==0,  rateMap = nan([size(Pfs.data.rateMap,1),1]);  end
% $$$         [mxr(u==units),mxp(u==units)] = max(rateMap.*mask);
        
      case 'mean'        

        sigMask = [];
        if ~isempty(sigThresh) && Pfs.parameters.numIter > 1,
            sigMask = sum(~isnan(Pfs.data.rateMap(:,Pfs.data.clu==u,:)),3)>Pfs.parameters.numIter*sigThresh;
        end 
        

        %rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==u,:),3,'omitnan');
        rateMap = Pfs.plot(u,'mean',false,[],false,0.99,false,interpPar);

        rateMap = RectFilter(rateMap');
        rateMap = RectFilter(rateMap');        
        rateMap = rateMap(:);
        
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

mxp = Ind2Sub(cellfun(@numel,bins),mxp);
if numel(Pfs.parameters.type)>2,
    mxp = [bins{1}(mxp(:,1)), ...
           bins{2}(mxp(:,2)), ...
           bins{3}(mxp(:,3))];                                
else
    mxp = [bins{1}(mxp(:,1)), ...
           bins{2}(mxp(:,2))];
end

% END MAIN -----------------------------------------------------------------------------------------