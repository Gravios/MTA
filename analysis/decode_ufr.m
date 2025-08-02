function [decEstCom,decEstMax,decEstSax,posteriorMax] = decode_ufr(Trial,varargin)
% function [decEstCom,decEstMax,decEstSax,decteriorMax] = decode_ufr(Trial,varargin)
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
                 'sampleRate',                    [30],                                          ...
                 'ufr',                           [],                                            ...
                 'pfs',                           [],                                            ...
                 'interpParPfs',                  [],                                            ...
                 'mask',                          [],                                            ...
                 'smoothingWeights',              [250.^2,250.^2],                               ...
                 'posteriorThresh',               [1e-5],                                        ...
                 'bufferSize',                    2^11,                                          ...
                 'tag',                           '',                                            ...
                 'overwrite',                     false                                          ...
);
[units,sampleRate,ufr,pfs,interpParPfs,mask,smoothingWeights,posteriorThresh,bufferSize,tag,overwrite] =...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


if isempty(tag)
% GENERATE unique hash for save file
    tag = DataHash({Trial.filebase,pfs.filename,sampleRate,units,spikeWindow});
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
    ufrLength = size(ufr,1);
    spikeWindow = ufr.spikeWindow;

% SET interpolation parameters
% $$$     pfsBins = cell([1,4]);
% $$$     if ~isempty(interpParPfs),
% $$$         pfsBins(1:numel(interpParPfs.bins)) = interpParPfs.bins;
% $$$     else,
% $$$         pfsBins = pfs.adata.bins;
% $$$     end
% $$$     pfsBinsDims = cellfun(@numel,pfsBins);
% $$$     pfsBinsDims(pfsBinsDims==0) = 1;
% $$$     pfsBins(cellfun(@isempty,pfsBins)) = [];
    

    
% ACCUMULATE masked rate maps
    ratemap = zeros([sum(mask(:)),0]);
    for u = 1:numel(units),
        if isa(pfs,'MTAApfs');
            trm = pfs.plot(units(u),1,false,[],false,0.25,false,interpParPfs);
            ratemap = cat(2,ratemap,trm(logical(mask(:))));
        else
            ratemap = pfs.data.rateMap(logical(mask(:)),:);
        end
    end
    
    
% COMPUTE Posterior distribution base on ratemaps and unit firing rates
    ratemap(isnan(ratemap)) = 0;
    ratemap = ratemap+1e-3;    

    decEstCom = nan([ufrLength,ndims]);
    decEstMax = decEstCom;
    decEstSax = decEstCom;
    posteriorMax = nan([ufrLength,1]);

    numDimPos = numel(pfsBins);

    binGrid = cell([1,numel(pfsBins)]); 
    [binGrid{:}] = ndgrid(pfsBins{:});

    gbinm = nan([size(ratemap,1),numel(pfsBins)]);
    for d = 1:ndims,
        gbinm(:,d) = cat(2,binGrid{d}(logical(mask(:))));
    end
    %clear('binGrid','mask');

    
    mpos = nan([bufferSize,numDimPos]);
    apos = nan([bufferSize,1]);            
    
    for i = 1:ceil(size(ufr,1)/bufferSize)
        tic
        disp([num2str(i),' of ' num2str(ceil(size(ufr,1)/bufferSize))]);

        
% SELECT data subset (bufferSize)
        if i == ceil(size(ufr,1)/bufferSize),
            ind = (1+bufferSize*(i-1)) : (bufferSize*(i-1)+mod(size(ufr,1),bufferSize));
            bufferSize = numel(ind);                        
            mpos = nan([bufferSize,numDimPos]);
            tpos = nan([bufferSize,numDimPos]);
            apos = nan([bufferSize,1]);            
        else        
            ind = (1+bufferSize*(i-1)) : (bufferSize*i);
        end
        
% COMPUTE posterior
        E = exp(-sum(ratemap,2)*ones([1,bufferSize])*spikeWindow+log(ratemap)*(ufr.data(ind,:)+eps)');
        E = bsxfun(@rdivide,E,sum(E));

% LOCATE peak of posterior
        [apos,tbin] = max(E);
        mpos = gbinm(tbin,:);  

% COMPUTE COM and weighted max of posterior
        tpos = nan([bufferSize,numDimPos]);
        cpos = nan([bufferSize,numDimPos]);        
        for tind = 1:bufferSize
            if apos(tind)<posteriorThresh,continue;end;
            wbinm = bsxfun(@minus,gbinm,gbinm(tbin(tind),:));                
            weights = exp(multiprod(-wbinm,multiprod(inv(smoothingWeights),wbinm,[1,2],[2]),[2],[2]));
            weights = weights./sum(weights);
            weights = bsxfun(@times,weights,repmat(E(:,tind),[1,numDimPos]));
            weights = weights./sum(weights);
            tpos(tind,:) = sum(gbinm.*weights,'omitnan');
            cpos(tind,:) = sum(gbinm.*repmat(E(:,tind),[1,numDimPos]),'omitnan');            
        end

% ACCUMULATE position estimates                
        decEstCom(ind,:)      = cpos;
        decEstSax(ind,:)      = tpos;
        decEstMax(ind,:)      = mpos;
        posteriorMax(ind)     = apos;

        toc
        

% $$$         if i == 4, 
% $$$ 
% $$$         
% $$$         % DIAGNOSTIC figure -----------------
% $$$         xyz = preproc_xyz(Trial,'trb');
% $$$         xyz.resample(sampleRate);
% $$$         hxyz = xyz(:,'hcom',1:2);
% $$$         nxyz = xyz(:,'nose',1:2);        
% $$$ 
% $$$         keyboard();
% $$$ 
% $$$         figure();
% $$$         hold('on');
% $$$         pax = pcolor(pfsBins{:},mask);
% $$$         pax.EdgeColor = 'none';
% $$$         lax = plot(0,0,'*m','MarkerSize',10);
% $$$         lax2 = plot(0,0,'*g','MarkerSize',10);
% $$$         lax3 = plot(0,0,'*c','MarkerSize',10);        
% $$$         for j = 1:1:2000,
% $$$             rmap = nan(size(mask));
% $$$             rmap(logical(mask(:))) = E(:,j);
% $$$             pax.CData = rmap';
% $$$             lax.XData = hxyz(ind(j),1);
% $$$             lax.YData = hxyz(ind(j),2);            
% $$$             lax2.XData = nxyz(ind(j),1);
% $$$             lax2.YData = nxyz(ind(j),2);            
% $$$             lax3.XData = decEstSax(ind(j),1);
% $$$             lax3.YData = decEstSax(ind(j),2);            
% $$$             drawnow();
% $$$             pause(.05);
% $$$         end
% $$$ 
% $$$         figure();
% $$$         histogram(sqrt(sum((decEstSax(ind,:)-nxyz(ind,:)).^2,2))./10,linspace(0,100,50));
% $$$ 
% $$$         % -----------------------------------
% $$$         end
        
        
        
    end% i

    
    % SAVE reconstruction statistics
    smoothingWeights = diag(smoothingWeights);    
    save(filepath,'decEstCom',   'decEstSax', 'decEstMax',...
                  'posteriorMax','sampleRate','spikeWindow','smoothingWeights');
end


% $$$ ufr = Trial.load('ufr',xyz,spk,unitSubset,spikeWindow.*2,true,'gauss');
% $$$ xyz = preproc_xyz(Trial,'trb');
% $$$ xyz.resample(sampleRate);
% $$$ fet = fet_HB_pitchB(Trial,sampleRate);
% $$$ 
% $$$ figure,
% $$$ unitInclusion = sum(ufr.data>0.2,2);    
% $$$ ind = nan([size(xyz,1),1]);
% $$$ ind( posteriorMax>0.0001&unitInclusion>=3) = 1;
% $$$ subplot(4,1,1);plot(bsxfun(@times,[xyz(:,5,1),decEstCom(:,1)],ind));    
% $$$ subplot(4,1,2);plot(bsxfun(@times,[xyz(:,5,2),decEstCom(:,2)],ind));        
% $$$ subplot(4,1,3);plot(bsxfun(@times,[fet(:,1),  decEstCom(:,3)],ind));
% $$$ subplot(4,1,4);plot(bsxfun(@times,[fet(:,2),  decEstCom(:,4)],ind));
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
