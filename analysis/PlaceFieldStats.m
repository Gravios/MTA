function [pfstats,pfbstats] = PlaceFieldStats(Trial,pf,unit)

version = '0.1';
verbose = true;

if verbose,
    fprintf('PlaceFieldStats v%s \n\n',version)
    fprintf('Trial: %s\n',Trial.filebase)
    fprintf('Unit: %i\n',unit)
end
    
map = sq(pf.data.rateMap(:,unit==pf.data.clu,:));
pkfr = max(map(:,1));
rateThreshold = prctile(sq(map(:)),95);
maxNumPatches =2;


%%!!!! Discrepency between MTAApfs and MTAAknnpfs_bs
%%!!!! Redo plot function for MTAAknpfs_bs to match MTAApfs
ratemap = reshape(map(:,1)',fliplr(pf.adata.binSizes'))';
%ratemap = pf.plot(unit,'isCircular',false);
%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


% Maximum firing rate found within the normal place field
pfstats.peakFR = pkfr;
pfstats.rateThreshold = rateThreshold;
pfstats.spatialCoherence = pf.spatialCoherence(unit);
%pfstats.spatialInformation = nansum((map(~isnan(map))./nanmean(map(:))).*log2(map(~isnan(map))./nanmean(map(:))));


pfstats.patchArea = nan([1,1,maxNumPatches]);
pfstats.patchCOM  = nan([1,1,maxNumPatches,2]);
pfstats.patchPFR = nan([1,1,maxNumPatches]);
pfstats.patchMFR = nan([1,1,maxNumPatches]);
pfstats.patchCnt = nan([1,1,maxNumPatches]);
pfstats.patchRateInd = nan([1,1,maxNumPatches,2,round(prod(pf.adata.binSizes)*0.1)]);
pfstats.patchRateMap = nan([1,1,maxNumPatches,round(prod(pf.adata.binSizes)*0.1)]);


if verbose,
    fprintf('patch detection:\n')
end


[B,L] = bwboundaries(ratemap>rateThreshold,'noholes');

if ~isempty(B),
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

        %MTAApfs        
        %pxy = [pf.adata.bins{1}(nB{i}(:,1)),pf.adata.bins{2}(nB{i}(:,2))];

        %MTAAknnpfs_bs
        pxy = [pf.adata.bins{2}(nB{i}(:,2)),pf.adata.bins{1}(nB{i}(:,1))];

        patchCOM(1,1,i,:)  = [sq(patchRateMap(1,1,i,1:npix))'*pxy(:,1),sq(patchRateMap(1,1,i,1:npix))'*pxy(:,2)]/nansum(patchRateMap(1,1,i,1:npix));
        patchPFR(1,1,i)    = max    (patchRateMap(1,1,i,1:npix),[],4);
        patchMFR(1,1,i)    = nanmean(patchRateMap(1,1,i,1:npix),4);

    end
    
    %keyboard


    [~,pind] = sort(sq(patchArea),'descend');
    if numel(nB)>=maxNumPatches,
        for i = 1:maxNumPatches,
            nNotNan = sum(~isnan(sq(patchRateMap(1,1,pind(i),:))));
            pfstats.patchRateInd(1,1,i,:,1:nNotNan) = patchRateInd(1,1,pind(i),:,:);
            pfstats.patchRateMap(1,1,i,1:nNotNan) = patchRateMap(1,1,pind(i),:);
        end
        pfstats.patchCOM(1,1,:,:)     = patchCOM(1,1,pind(1:maxNumPatches),:);
        pfstats.patchArea(1,1,:)    = patchArea(1,1,pind(1:maxNumPatches));
        pfstats.patchPFR(1,1,:)     = patchPFR(1,1,pind(1:maxNumPatches));
        pfstats.patchMFR(1,1,:)     = patchMFR(1,1,pind(1:maxNumPatches));
    else
        for i = 1:maxNumPatches,
            nNotNan = sum(~isnan(sq(patchRateMap(1,1,1,:))));
            pfstats.patchRateInd(1,1,1,:,1:nNotNan) = patchRateInd(1,1,1,:,:);
            pfstats.patchRateMap(1,1,1,1:nNotNan) = patchRateMap(1,1,1,:);
        end
        pfstats.patchCOM(1,1,1,:)   = patchCOM (1,1,1,:);
        pfstats.patchArea(1,1,1)      = patchArea(1,1,:);
        pfstats.patchPFR(1,1,1)       = patchPFR (1,1,:);
        pfstats.patchMFR(1,1,1)       = patchMFR (1,1,:);
    end

    
end

if verbose,
    fprintf('Boundary detection on bootstraped placefields:\n')
end


if nargout>1
    
    pfbstats.patchArea =    nan([1,pf.parameters.numIter,maxNumPatches]);
    pfbstats.patchCOM  =    nan([1,pf.parameters.numIter,maxNumPatches,2]);
    pfbstats.patchPFR  =    nan([1,pf.parameters.numIter,maxNumPatches]);
    pfbstats.patchMFR  =    nan([1,pf.parameters.numIter,maxNumPatches]);
    pfbstats.patchCnt  =    nan([1,pf.parameters.numIter,maxNumPatches]);
    pfbstats.patchRateInd = nan([1,pf.parameters.numIter,maxNumPatches,2,round(prod(pf.adata.binSizes)*0.1)]);
    pfbstats.patchRateMap = nan([1,pf.parameters.numIter,maxNumPatches,round(prod(pf.adata.binSizes)*0.1)]);

    for k = 1:pf.parameters.numIter,    
        if verbose,
            if ~mod(k,100),
                fprintf('.\n')
            elseif mod(k,100)==1,
                fprintf('bs iter:%i.',k)
            else
                fprintf('.',k)
            end
        end        
        
        % if k == 86,keyboard,end
        %MTAApfs
        %ratemap = reshape(map(:,k),numel(bin1),numel(bin2)).*mask;    

        %MTAAknnpfs_bs
        ratemap = reshape(map(:,k)',fliplr(pf.adata.binSizes'))';
        clear nB B L npix bsPatchArea bsPatchRateMap bsPatchRateInd ...
              pxy bsPatchCOM bsPatchPFR bsPatchMFR pind 
        [B,L] = bwboundaries(ratemap>rateThreshold,'noholes');

        if ~isempty(B),
                
            for i = 1:numel(B),
                [nB{i}(:,1),nB{i}(:,2)] = ind2sub(size(L),find(L==i));
                npix = size(nB{i},1);
                
                bsPatchArea(1,1,i) = prod(pf.parameters.binDims)*npix;
                bsPatchRateMap(1,1,i,1:npix) = ratemap(sub2ind(pf.adata.binSizes',nB{i}(:,1),nB{i}(:,2)));

                bsPatchRateInd(1,1,i,:,1:npix) = nB{i}';

                %pxy = [pf.adata.bins{1}(nB{i}(:,1)),pf.adata.bins{2}(nB{i}(:,2))];
                pxy = [pf.adata.bins{2}(nB{i}(:,2)),pf.adata.bins{1}(nB{i}(:,1))];

                bsPatchCOM(1,1,i,:)  = [sq(bsPatchRateMap(1,1,i,1:npix))'*pxy(:,1),sq(bsPatchRateMap(1,1,i,1:npix))'*pxy(:,2)]/nansum(bsPatchRateMap(1,1,i,1:npix));
                bsPatchPFR(1,1,i)    = max    (bsPatchRateMap(1,1,i,1:npix),[],4);
                bsPatchMFR(1,1,i)    = nanmean(bsPatchRateMap(1,1,i,1:npix),4);
                %%%%%% start here
            end
            
            [~,pind] = sort(sq(bsPatchArea),'descend');
            if numel(nB)>=maxNumPatches,
                for i = 1:maxNumPatches,
                    nNotNan = sum(~isnan(sq(bsPatchRateMap(1,1,pind(i),:))));
                    pfbstats.patchRateInd(1,k,i,:,1:nNotNan) = bsPatchRateInd(1,1,pind(i),:,:);
                    pfbstats.patchRateMap(1,k,i,1:nNotNan) = bsPatchRateMap(1,1,pind(i),:);
                end
                pfbstats.patchCOM(1,k,:,:)     = bsPatchCOM(1,1,pind(1:maxNumPatches),:);
                pfbstats.patchArea(1,k,:)    = bsPatchArea(1,1,pind(1:maxNumPatches));
                pfbstats.patchPFR(1,k,:)     = bsPatchPFR(1,1,pind(1:maxNumPatches));
                pfbstats.patchMFR(1,k,:)     = bsPatchMFR(1,1,pind(1:maxNumPatches));
            else
                for i = 1:maxNumPatches,
                    nNotNan = sum(~isnan(sq(bsPatchRateMap(1,1,1,:))));
                    pfbstats.patchRateInd(1,k,1,:,1:nNotNan) = bsPatchRateInd(1,1,1,:,:);
                    pfbstats.patchRateMap(1,k,1,1:nNotNan) = bsPatchRateMap(1,1,1,:);
                end
                pfbstats.patchCOM(1,k,1,:)     = bsPatchCOM (1,1,1,:);
                pfbstats.patchArea(1,k,1)      = bsPatchArea(1,1,:);
                pfbstats.patchPFR(1,k,1)       = bsPatchPFR (1,1,:);
                pfbstats.patchMFR(1,k,1)       = bsPatchMFR (1,1,:);
            end
        end
    end
end

if verbose,
    fprintf('Status: Complete\n')
end


% $$$ 
% $$$ if nargout>1,        
% $$$     pfshuff.peakFR = zeros(pf.parameters.numIter,1);
% $$$     pfshuff.patchArea = zeros(pf.parameters.numIter,1);
% $$$     pfshuff.patchPFR = zeros(pf.parameters.numIter,1);
% $$$     pfshuff.patchMFR = zeros(pf.parameters.numIter,1);
% $$$     pfshuff.patchCnt = zeros(pf.parameters.numIter,1);;
% $$$     %% Patch stats
% $$$     bsB = cell(pf.parameters.numIter,1);
% $$$     bsPFR = zeros(pf.parameters.numIter,1);
% $$$     bsMFR = zeros(pf.parameters.numIter,1);
% $$$ 
% $$$     for i = 1:pf.parameters.numIter,
% $$$         bsB{i} = bwboundaries(reshape(map(:,i),pf.adata.binSizes')'>pfstats.rateThreshold);
% $$$         if ~isempty(bsB{i}),
% $$$             [~,tmpi] = max(cellfun(@numel,bsB{i}));
% $$$             bsB{i} = bsB{i}{tmpi};
% $$$ 
% $$$             bsMFR(i) = mean(map(sub2ind(pf.adata.binSizes',bsB{i}(:,2),bsB{i}(:,1)),i));
% $$$             bsPFR(i) = max(map(sub2ind(pf.adata.binSizes',bsB{i}(:,2),bsB{i}(:,1)),i));
% $$$         end
% $$$     end
% $$$ 
% $$$     pfshuff.peakFR = sq(max(map))';
% $$$     pfshuff.patchArea = cellfun(@size,bsB,repmat({1},pf.parameters.numIter,1));
% $$$     pfshuff.patchPFR = bsPFR;
% $$$     pfshuff.patchMFR = bsMFR;
% $$$ end
% $$$ 

% $$$ if nargout>1,        
