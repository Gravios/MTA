function [dc] = decode_ufr_boxcar(Trial,varargin)
% function [dc] = decode_ufr(Trial,varargin)
% 
% 
%
% Input:
%    Trial
%    sampleRate
%    units
%    mode
%    overwrite
%
% Output:
%    decEstCom
%    decEstMax
%    decteriorMax
%

global MTA_PROJECT_PATH

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',                         [],                                            ...
                 'sampleRate',                    [20],                                         ...
                 'ufr',                           [],                                            ...
                 'pfs',                           [],                                            ...
                 'mask',                          [],                                            ...
                 'halfSpkWindow',                 0.15,                                          ...
                 'smoothingWeights',              [250.^2,250.^2],                               ...
                 'tag',                           '',                                            ...
                 'overwrite',                     false                                          ...
);
[units, sampleRate, ufr, pfs, mask, halfSpkWindow, smoothingWeights, tag, overwrite] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


if isempty(tag)
% GENERATE unique hash for save file
    tag = DataHash({Trial.filebase,       ...
                    units,                ...                    
                    sampleRate,           ...
                    pfs.filename,         ...
                    mask,                 ...
                    halfSpkWindow,        ...                    
                    smoothingWeights      ...
                   });
end

% GENERATE path for save file 
filepath = fullfile(Trial.spath,[Trial.filebase,'.',mfilename,'.',tag,'.mat']);

if exist(filepath,'file') && ~overwrite,    
    load(filepath);
else

% SET parameters
    smoothingWeights = diag(smoothingWeights);
    pfsBins = pfs.adata.bins;
    ndims = numel(pfs.adata.bins);
    
    spkWindow = round(halfSpkWindow.*sampleRate).*[-1,1]; %s*i/s
    
% ACCUMULATE masked rate maps
    ratemap = zeros([sum(mask(:)),0]);
    for u = 1:numel(units),
        if isa(pfs,'MTAApfs');
            trm = pfs.plot(units(u),1,false,[],false,0.25,false);
            ratemap = cat(2,ratemap,trm(logical(mask(:))));
        else
            ratemap = pfs.data.rateMap(logical(mask(:)),:);
        end
    end
    
% COMPUTE Posterior distribution base on ratemaps and unit firing rates
    ratemap(isnan(ratemap)) = 0;
    ratemap = ratemap+1e-3;    
    logRatemap = log(ratemap);
    priorRatemap = -sum(ratemap,2)*ufr.spikeWindow;

% ESTIMATE number of valid samples
    esamp = size(nonzeros(sum(ufr.data>=1/(ufr.spikeWindow*sampleRate),2)),1).*3;
    
% INSTANTIATE output 
    dc.ind  = nan([esamp*2,1]);    
    dc.max  = nan([esamp*2,ndims]);    
    dc.com  = nan([esamp*2,ndims]);
    dc.sax  = nan([esamp*2,ndims]);
    dc.lom  = nan([esamp*2,ndims]);
    dc.lax  = nan([esamp*2,ndims]);
    dc.post = nan([esamp*2,1]);
    dc.ucnt = nan([esamp*2,1]);    
    dc.uinc = nan([esamp*2,size(ufr,2)]);    

% SELECT ratemap bins within mask
    binGrid = cell([1,numel(pfsBins)]); 
    [binGrid{:}] = ndgrid(pfsBins{:});
    gbinm = nan([size(ratemap,1),numel(pfsBins)]);
    for d = 1:ndims,
        gbinm(:,d) = cat(2,binGrid{d}(logical(mask(:))));
    end

    cindex = 1;    
    tic        
    for tind = (1+abs(spkWindow(1))):(size(ufr,1)-abs(spkWindow(2))-10)

        if mod(cindex,10000)==0,toc,disp(cindex);end
        twind = tind+spkWindow;

        cufr = sum(ufr(twind(1):twind(2),:));

        gind = cufr>0.5;
        
        if any(gind)
            dc.ind (cindex,1) = tind;
            dc.ucnt(cindex,1) = sum(gind);            
            dc.uinc(cindex,:) = gind;

% COMPUTE posterior
            E = exp(priorRatemap+logRatemap*(cufr+eps)');
            E = E./sum(E,'omitnan');

% LOCATE peak of posterior
            [dc.post(cindex,1),tbin] = max(E);
            dc.max(cindex,:) = gbinm(tbin,:);  

% COMPUTE COM and weighted max of posterior
            wbinm = bsxfun(@minus,gbinm,gbinm(tbin,:));                
            weights = exp(multiprod(-wbinm,multiprod(inv(smoothingWeights),wbinm,[1,2],[2]),[2],[2]));
            weights = weights./sum(weights,'omitnan');
            weights = bsxfun(@times,weights,repmat(E(:),[1,ndims]));
            weights = weights./sum(weights,'omitnan');
            dc.sax(cindex,:) = sum(gbinm.*weights,'omitnan');
            dc.com(cindex,:) = sum(gbinm.*repmat(E(:),[1,ndims]),'omitnan');     
            
            lpost = log10(E(:));
            lpost = lpost+8;
            lpost(lpost<=0) = eps;
            lpost = lpost./sum(lpost,'omitnan');
            dc.lom(cindex,:) = sum(gbinm.*repmat(lpost(:),[1,ndims]),'omitnan');
    
            wbinm = bsxfun(@minus,gbinm,gbinm(tbin,:));
            weights = exp(multiprod(-wbinm,multiprod(inv(smoothingWeights),wbinm,[1,2],[2]),[2],[2]));
            weights = weights./sum(weights,'omitnan');
            weights = bsxfun(@times,weights,repmat(lpost(:),[1,ndims]));
            weights = weights./sum(weights,'omitnan');
            dc.lax(cindex,:) = sum(gbinm.*weights,'omitnan');
            
            cindex = cindex + 1;
            
        end%if
        
    end%for;tind
    

% REMOVE nan data points
    dc.max (isnan(dc.ind),:) = [];
    dc.com (isnan(dc.ind),:) = [];
    dc.sax (isnan(dc.ind),:) = [];
    dc.lom (isnan(dc.ind),:) = [];
    dc.lax (isnan(dc.ind),:) = [];
    dc.post(isnan(dc.ind),:) = [];
    dc.ucnt(isnan(dc.ind),:) = [];
    dc.uinc(isnan(dc.ind),:) = [];
    dc.ind (isnan(dc.ind),:) = [];    

% SAVE reconstruction statistics
    dc.smoothingWeights = diag(smoothingWeights);    
    dc.window = 2*halfSpkWindow;
    dc.sampleRate = sampleRate;

    save(filepath,'dc');
end


