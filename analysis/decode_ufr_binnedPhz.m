function [dc] = decode_ufr_binnedPhz(Trial,varargin)
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
                 'sampleRate',                    [250],                                         ...
                 'ufr',                           [],                                            ...
                 'phz',                           [],                                            ...
                 'pfs',                           [],                                            ...
                 'mask',                          [],                                            ...
                 'halfSpkWindow',                 0.15,                                          ...
                 'halfPhzWindow',                 0.25,                                          ...
                 'phzBins',                       linspace(0,2*pi,25),                           ...
                 'smoothingWeights',              [250.^2,250.^2],                               ...
                 'tag',                           '',                                            ...
                 'overwrite',                     false                                          ...
);
[units, sampleRate, ufr, phz, pfs, mask, halfSpkWindow, halfPhzWindow, phzBins, smoothingWeights,...
 tag, overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


if isempty(tag)
% GENERATE unique hash for save file
    tag = DataHash({Trial.filebase,       ...
                    units,                ...                    
                    sampleRate,           ...
                    phz.hash,             ...
                    phzBins,              ...
                    pfs.filename,         ...
                    mask,                 ...
                    halfSpkWindow,        ...                    
                    halfPhzWindow,        ...                                        
                    phzBins,              ...
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
    
    phzWindow = round(halfPhzWindow./diff(phzBins(1:2)));
    phzWindow = -phzWindow:1:phzWindow;

    spkWindow = round(halfSpkWindow.*sampleRate).*[-1,1]; %s*i/s
    
    iphz = discretize(phz.data,phzBins);    
    
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

% ESTIMATE number of valid samples
    esamp = size(nonzeros(sum(ufr.data>=1/(ufr.spikeWindow*sampleRate),2)),1).*3;
    
% INSTANTIATE output 
    dc.ind  = nan([esamp*2,1]);    
    dc.max  = nan([esamp*2,ndims]);    
    dc.com  = nan([esamp*2,ndims]);
    dc.sax  = nan([esamp*2,ndims]);
    dc.post = nan([esamp*2,1]);
    dc.phz  = nan([esamp*2,1]);
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
        cphz = phz(tind);
        pwind = mod(iphz(tind)-1+phzWindow,numel(phzBins)-1)+1;

        tufr = ufr(twind(1):twind(2),:);
        tphz = iphz(twind(1):twind(2));
        
        cufr = sum(tufr(ismember(tphz,pwind),:));
        gind = cufr>0.5;
        
        if any(gind)
            dc.ind (cindex,1) = tind;
            dc.phz (cindex,1) = cphz;
            dc.ucnt(cindex,1) = sum(gind);            
            dc.uinc(cindex,:) = gind;

% COMPUTE posterior
            E = exp(-sum(ratemap,2)*ufr.spikeWindow+log(ratemap)*(cufr+eps)');
            E = E./sum(E,'omitnan');

% LOCATE peak of posterior
            [dc.post(cindex,1),tbin] = max(E);
            dc.max(cindex,:) = gbinm(tbin,:);  

% COMPUTE COM and weighted max of posterior
            wbinm = bsxfun(@minus,gbinm,gbinm(tbin,:));                
            weights = exp(multiprod(-wbinm,multiprod(inv(smoothingWeights),wbinm,[1,2],[2]),[2],[2]));
            weights = weights./sum(weights);
            weights = bsxfun(@times,weights,repmat(E(:),[1,ndims]));
            weights = weights./sum(weights);
            dc.sax(cindex,:) = sum(gbinm.*weights,'omitnan');
            dc.com(cindex,:) = sum(gbinm.*repmat(E(:),[1,ndims]),'omitnan');     
            
            cindex = cindex + 1;
            
        end%if
        
    end%for;tind
    

% REMOVE nan data points
    dc.max (isnan(dc.ind),:) = [];
    dc.com (isnan(dc.ind),:) = [];
    dc.sax (isnan(dc.ind),:) = [];
    dc.post(isnan(dc.ind),:) = [];
    dc.phz (isnan(dc.ind),:) = [];
    dc.ucnt(isnan(dc.ind),:) = [];
    dc.uinc(isnan(dc.ind),:) = [];
    dc.ind (isnan(dc.ind),:) = [];    

% SAVE reconstruction statistics
    dc.smoothingWeights = diag(smoothingWeights);    
    dc.window = 2*halfSpkWindow;
    dc.sampleRate = sampleRate;
    dc.phzBins = phzBins;

    save(filepath,'dc');
end


