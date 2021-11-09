function [pfstats,pfbstats,pfmstats] = PlaceFieldStats(Trial,pf,unit,varargin)


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('maxNumPatches',                    2,                                          ...
                 'verbose',                          true,                                        ...
                 'maskFlag',                         false                                        ...
);
[maxNumPatches,verbose,maskFlag] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

VERSION = '0.1';
%thresholdMethod = 'shuffle';
thresholdMethod = 'percentile';

if verbose,
    fprintf('PlaceFieldStats v%s \n\n',version)
    fprintf('Trial: %s\n',Trial.filebase)
    fprintf('Unit: %i\n',unit)
end
    
map = reshape(pf.plot(unit,'mazeMaskFlag',maskFlag),[],1);
%map = sq(pf.data.rateMap(:,unit==pf.data.clu,:));
pkfr = max(map(:,1));

rateThreshold = select_rate_threshold(Trial,pf,unit,thresholdMethod); %helper function (see end)

%maxNumPatches =2;
maxBinCount = round(prod(pf.adata.binSizes)*0.22);


%%!!!! Discrepency between MTAApfs and MTAAknnpfs_bs
%ratemap = reshape(map(:,1)',fliplr(pf.adata.binSizes'))';


ratemap = pf.plot(unit,'mazeMaskFlag',maskFlag);


% Maximum firing rate found within the normal place field
pfstats.peakFR = pkfr;
pfstats.rateThreshold = rateThreshold;
pfstats.spatialCoherence = pf.spatialCoherence(unit,false);
%pfstats.spatialInformation = nansum((map(~isnan(map))./nanmean(map(:))).*log2(map(~isnan(map))./nanmean(map(:))));


pfstats.patchArea = nan([1,1,maxNumPatches]);
pfstats.patchCOM  = nan([1,1,maxNumPatches,2]);
pfstats.patchPFR = nan([1,1,maxNumPatches]);
pfstats.patchMFR = nan([1,1,maxNumPatches]);
pfstats.patchCnt = nan([1,1,maxNumPatches]);
pfstats.patchRateInd = nan([1,1,maxNumPatches,2,maxBinCount]);
pfstats.patchRateMap = nan([1,1,maxNumPatches,maxBinCount]);  
 

if verbose,
    fprintf('patch detection:\n')
end


[B,L] = bwboundaries(ratemap>rateThreshold,'noholes');

if ~isempty(B),
 
    
    patchArea = nan([1,1,maxNumPatches]);                                        
    patchCOM  = nan([1,1,maxNumPatches,2]);                                      
    patchPFR = nan([1,1,maxNumPatches]);                                         
    patchMFR = nan([1,1,maxNumPatches]);                                         
    patchCnt = nan([1,1,maxNumPatches]);                                         
    patchRateInd = nan([1,1,maxNumPatches,2,maxBinCount]);
    patchRateMap = nan([1,1,maxNumPatches,maxBinCount]);  
    
    clear nB    
    if verbose,
        fprintf('%i patches detected\n')
    end

    for i = 1:numel(B)
        [nB{i}(:,1),nB{i}(:,2)] = ind2sub(size(L),find(L==i));
        npix = size(nB{i},1);

        patchArea(1,1,i) = prod(pf.parameters.binDims)*npix;
        patchRateMap(1,1,i,1:npix) = ratemap(sub2ind(pf.adata.binSizes',nB{i}(:,1),nB{i}(:,2)));

        patchRateInd(1,1,i,:,1:npix) = nB{i}';

        if isa(pf,'MTAApfs'),
            pxy = [pf.adata.bins{1}(nB{i}(:,1)),pf.adata.bins{2}(nB{i}(:,2))];
        elseif isa(pf,'MTAAknnpfs_bs')
            pxy = [pf.adata.bins{2}(nB{i}(:,2)),pf.adata.bins{1}(nB{i}(:,1))];
        end
        
        patchCOM(1,1,i,:)  = [sq(patchRateMap(1,1,i,1:npix))'*pxy(:,1),sq(patchRateMap(1,1,i,1:npix))'*pxy(:,2)]/nansum(patchRateMap(1,1,i,1:npix));
        patchPFR(1,1,i)    = max    (patchRateMap(1,1,i,1:npix),[],4);
        patchMFR(1,1,i)    = nanmean(patchRateMap(1,1,i,1:npix),4);

    end
    


    [~,pind] = sort(sq(patchArea),'descend');
    if numel(nB)>=maxNumPatches,
        for i = 1:maxNumPatches,
            nNotNan = sum(~isnan(sq(patchRateMap(1,1,pind(i),:))));
            pfstats.patchRateInd(1,1,i,:,1:nNotNan) = patchRateInd(1,1,pind(i),:,1:nNotNan);
            pfstats.patchRateMap(1,1,i,1:nNotNan) = patchRateMap(1,1,pind(i),1:nNotNan);
        end
        pfstats.patchCOM(1,1,:,:)     = patchCOM(1,1,pind(1:maxNumPatches),:);
        pfstats.patchArea(1,1,:)    = patchArea(1,1,pind(1:maxNumPatches));
        pfstats.patchPFR(1,1,:)     = patchPFR(1,1,pind(1:maxNumPatches));
        pfstats.patchMFR(1,1,:)     = patchMFR(1,1,pind(1:maxNumPatches));
    else
        for i = 1:maxNumPatches,
            nNotNan = sum(~isnan(sq(patchRateMap(1,1,1,:))));
            pfstats.patchRateInd(1,1,1,:,1:nNotNan) = patchRateInd(1,1,1,:,1:nNotNan);
            pfstats.patchRateMap(1,1,1,1:nNotNan) = patchRateMap(1,1,1,1:nNotNan);
        end
        pfstats.patchCOM(1,1,1,:)   = patchCOM (1,1,1,:);
        pfstats.patchArea(1,1,1)      = patchArea(1,1,1);
        pfstats.patchPFR(1,1,1)       = patchPFR (1,1,1);
        pfstats.patchMFR(1,1,1)       = patchMFR (1,1,1);
    end

    
end






if verbose,
    fprintf('Boundary detection on bootstraped placefields:\n')
end


if nargout>1

    
    
    % Maximum firing rate found within the normal place field
    pfmstats.peakFR = max(nanmean(map,2));
    pfmstats.rateThreshold = rateThreshold;
    pfmstats.spatialCoherence = pf.spatialCoherence(unit,false);
    %pfmstats.spatialInformation = nansum((map(~isnan(map))./nanmean(map(:))).*log2(map(~isnan(map))./nanmean(map(:))));


    pfmstats.patchArea = nan([1,1,maxNumPatches]);
    pfmstats.patchCOM  = nan([1,1,maxNumPatches,2]);
    pfmstats.patchPFR = nan([1,1,maxNumPatches]);
    pfmstats.patchMFR = nan([1,1,maxNumPatches]);
    pfmstats.patchCnt = nan([1,1,maxNumPatches]);
    pfmstats.patchRateInd = nan([1,1,maxNumPatches,2,maxBinCount]);
    pfmstats.patchRateMap = nan([1,1,maxNumPatches,maxBinCount]);


    if verbose,
        fprintf('patch detection:\n')
    end

    ratemap = pf.plot(unit,'mean','mazeMaskFlag',maskFlag);
        
    [B,L] = bwboundaries(ratemap>rateThreshold,'noholes');

    if ~isempty(B),        
        clear nB    
        if verbose,
            fprintf('%i patches detected\n')
        end

        
        patchArea = nan([1,1,maxNumPatches]);                                        
        patchCOM  = nan([1,1,maxNumPatches,2]);                                      
        patchPFR = nan([1,1,maxNumPatches]);                                         
        patchMFR = nan([1,1,maxNumPatches]);                                         
        patchCnt = nan([1,1,maxNumPatches]);                                         
        patchRateInd = nan([1,1,maxNumPatches,2,maxBinCount]);
        patchRateMap = nan([1,1,maxNumPatches,maxBinCount]);          
        
        for i = 1:numel(B)
            [nB{i}(:,1),nB{i}(:,2)] = ind2sub(size(L),find(L==i));
            npix = size(nB{i},1);

            patchArea(1,1,i) = prod(pf.parameters.binDims)*npix;
            patchRateMap(1,1,i,1:npix) = ratemap(sub2ind(pf.adata.binSizes',nB{i}(:,1),nB{i}(:,2)));

            patchRateInd(1,1,i,:,1:npix) = nB{i}';

            if isa(pf,'MTAApfs'),
                pxy = [pf.adata.bins{1}(nB{i}(:,1)),pf.adata.bins{2}(nB{i}(:,2))];
            elseif isa(pf,'MTAAknnpfs_bs')
                pxy = [pf.adata.bins{2}(nB{i}(:,2)),pf.adata.bins{1}(nB{i}(:,1))];
            end
            
            patchCOM(1,1,i,:)  = [sq(patchRateMap(1,1,i,1:npix))'*pxy(:,1),sq(patchRateMap(1,1,i,1:npix))'*pxy(:,2)]/nansum(patchRateMap(1,1,i,1:npix));
            patchPFR(1,1,i)    = max    (patchRateMap(1,1,i,1:npix),[],4);
            patchMFR(1,1,i)    = nanmean(patchRateMap(1,1,i,1:npix),4);

        end
        


        [~,pind] = sort(sq(patchArea),'descend');
        if numel(nB)>=maxNumPatches,
            for i = 1:maxNumPatches,
                nNotNan = sum(~isnan(sq(patchRateMap(1,1,pind(i),:))));
                pfmstats.patchRateInd(1,1,i,:,1:nNotNan) = patchRateInd(1,1,pind(i),:,1:nNotNan);
                pfmstats.patchRateMap(1,1,i,1:nNotNan) = patchRateMap(1,1,pind(i),1:nNotNan);
            end
            pfmstats.patchCOM(1,1,:,:)     = patchCOM(1,1,pind(1:maxNumPatches),:);
            pfmstats.patchArea(1,1,:)    = patchArea(1,1,pind(1:maxNumPatches));
            pfmstats.patchPFR(1,1,:)     = patchPFR(1,1,pind(1:maxNumPatches));
            pfmstats.patchMFR(1,1,:)     = patchMFR(1,1,pind(1:maxNumPatches));
        else
            for i = 1:maxNumPatches,
                nNotNan = sum(~isnan(sq(patchRateMap(1,1,1,:))));
                pfmstats.patchRateInd(1,1,1,:,1:nNotNan) = patchRateInd(1,1,1,:,1:nNotNan);
                pfmstats.patchRateMap(1,1,1,1:nNotNan) = patchRateMap(1,1,1,1:nNotNan);
            end
            pfmstats.patchCOM(1,1,1,:)   = patchCOM (1,1,1,:);
            pfmstats.patchArea(1,1,1)      = patchArea(1,1,1);
            pfmstats.patchPFR(1,1,1)       = patchPFR (1,1,1);
            pfmstats.patchMFR(1,1,1)       = patchMFR (1,1,1);
        end

        
    end

    
    
    
    
    pfbstats.patchArea =    nan([1,pf.parameters.numIter,maxNumPatches]);
    pfbstats.patchCOM  =    nan([1,pf.parameters.numIter,maxNumPatches,2]);
    pfbstats.patchPFR  =    nan([1,pf.parameters.numIter,maxNumPatches]);
    pfbstats.patchMFR  =    nan([1,pf.parameters.numIter,maxNumPatches]);
    pfbstats.patchCnt  =    nan([1,pf.parameters.numIter,maxNumPatches]);
    pfbstats.patchRateInd = nan([1,pf.parameters.numIter,maxNumPatches,2,maxBinCount]);
    pfbstats.patchRateMap = nan([1,pf.parameters.numIter,maxNumPatches,maxBinCount]);
% $$$ 
% $$$     for k = 1:pf.parameters.numIter,    
% $$$         if verbose,
% $$$             if ~mod(k,100),
% $$$                 fprintf('.\n')
% $$$             elseif mod(k,100)==1,
% $$$                 fprintf('bs iter:%i.',k)
% $$$             else
% $$$                 fprintf('.',k)
% $$$             end
% $$$         end        
% $$$         
% $$$         
% $$$         bsPatchArea =    nan([1,1,maxNumPatches]);                                     
% $$$         bsPatchCOM  =    nan([1,1,maxNumPatches,2]);                                   
% $$$         bsPatchPFR  =    nan([1,1,maxNumPatches]);                                     
% $$$         bsPatchMFR  =    nan([1,1,maxNumPatches]);                                     
% $$$         bsPatchCnt  =    nan([1,1,maxNumPatches]);                                     
% $$$         bsPatchRateInd = nan([1,1,maxNumPatches,2,maxBinCount]);
% $$$         bsPatchRateMap = nan([1,1,maxNumPatches,maxBinCount]);  
% $$$         
% $$$         
% $$$         % if k == 86,keyboard,end
% $$$         %MTAApfs
% $$$         %ratemap = reshape(map(:,k),numel(bin1),numel(bin2)).*mask;    
% $$$ 
% $$$         ratemap = pf.plot(unit,k,'mazeMaskFlag',maskFlag);
% $$$         
% $$$         %MTAAknnpfs_bs
% $$$         %ratemap = reshape(map(:,k)',fliplr(pf.adata.binSizes'))';
% $$$         clear nB B L npix bsPatchArea bsPatchRateMap bsPatchRateInd ...
% $$$               pxy bsPatchCOM bsPatchPFR bsPatchMFR pind 
% $$$         [B,L] = bwboundaries(ratemap>rateThreshold,'noholes');
% $$$ 
% $$$         if ~isempty(B),
% $$$                 
% $$$             for i = 1:numel(B),
% $$$                 [nB{i}(:,1),nB{i}(:,2)] = ind2sub(size(L),find(L==i));
% $$$                 npix = size(nB{i},1);
% $$$                 
% $$$                 bsPatchArea(1,1,i) = prod(pf.parameters.binDims)*npix;
% $$$                 bsPatchRateMap(1,1,i,1:npix) = ratemap(sub2ind(pf.adata.binSizes',nB{i}(:,1),nB{i}(:,2)));
% $$$ 
% $$$                 bsPatchRateInd(1,1,i,:,1:npix) = nB{i}';
% $$$ 
% $$$                 %pxy = [pf.adata.bins{1}(nB{i}(:,1)),pf.adata.bins{2}(nB{i}(:,2))];
% $$$                 pxy = [pf.adata.bins{2}(nB{i}(:,2)),pf.adata.bins{1}(nB{i}(:,1))];
% $$$ 
% $$$                 bsPatchCOM(1,1,i,:)  = [sq(bsPatchRateMap(1,1,i,1:npix))'*pxy(:,1),sq(bsPatchRateMap(1,1,i,1:npix))'*pxy(:,2)]/nansum(bsPatchRateMap(1,1,i,1:npix));
% $$$                 bsPatchPFR(1,1,i)    = max    (bsPatchRateMap(1,1,i,1:npix),[],4);
% $$$                 bsPatchMFR(1,1,i)    = nanmean(bsPatchRateMap(1,1,i,1:npix),4);
% $$$                 %%%%%% start here
% $$$             end
% $$$             
% $$$             [~,pind] = sort(sq(bsPatchArea),'descend');
% $$$             if numel(nB)>=maxNumPatches,
% $$$                 for i = 1:maxNumPatches,
% $$$                     nNotNan = sum(~isnan(sq(bsPatchRateMap(1,1,pind(i),:))));
% $$$                     pfbstats.patchRateInd(1,k,i,:,1:nNotNan) = bsPatchRateInd(1,1,pind(i),:,1:nNotNan);
% $$$                     pfbstats.patchRateMap(1,k,i,1:nNotNan) = bsPatchRateMap(1,1,pind(i),1:nNotNan);
% $$$                 end
% $$$                 pfbstats.patchCOM(1,k,:,:)     = bsPatchCOM(1,1,pind(1:maxNumPatches),:);
% $$$                 pfbstats.patchArea(1,k,:)    = bsPatchArea(1,1,pind(1:maxNumPatches));
% $$$                 pfbstats.patchPFR(1,k,:)     = bsPatchPFR(1,1,pind(1:maxNumPatches));
% $$$                 pfbstats.patchMFR(1,k,:)     = bsPatchMFR(1,1,pind(1:maxNumPatches));
% $$$             else
% $$$                 for i = 1:maxNumPatches,
% $$$                     nNotNan = sum(~isnan(sq(bsPatchRateMap(1,1,1,:))));
% $$$                     pfbstats.patchRateInd(1,k,1,:,1:nNotNan) = bsPatchRateInd(1,1,1,:,:);
% $$$                     pfbstats.patchRateMap(1,k,1,1:nNotNan) = bsPatchRateMap(1,1,1,:);
% $$$                 end
% $$$                 pfbstats.patchCOM(1,k,1,:)     = bsPatchCOM (1,1,1,:);
% $$$                 pfbstats.patchArea(1,k,1)      = bsPatchArea(1,1,:);
% $$$                 pfbstats.patchPFR(1,k,1)       = bsPatchPFR (1,1,:);
% $$$                 pfbstats.patchMFR(1,k,1)       = bsPatchMFR (1,1,:);
% $$$             end
% $$$         end
% $$$     end
end

if verbose,
    fprintf('Status: Complete\n')
end



function threshold = select_rate_threshold(Trial,pf,unit,method)
switch method
  case 'shuffle'
    defargs = get_default_args('MjgER2016','MTAAknnpfs','struct');
    defargs.units = unit;
    defargs.states = 'theta';
    defargs.overwrite = false;
    defargs = struct2varargin(defargs);        
    pf = MTAAknnpfs(Trial,defargs{:});      
    
    mmap = pf.plot(unit,'mean','mazeMaskFlag',maskFlag);
    smap = pf.plot(unit,'std','mazeMaskFlag',maskFlag);
    threshold = nanmean(mmap(:))+nanmean(smap(:))*3;
  case 'percentile'    
    map = sq(pf.data.rateMap(:,unit==pf.data.clu,:));    
    threshold = prctile(sq(map(:)),90);
end
