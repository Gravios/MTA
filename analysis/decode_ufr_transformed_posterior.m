function [decEstCom,decEstMax,decEstSax,posteriorMax,posteriorPatch] = decode_ufr_transformed_posterior(Trial,varargin)
% function [decEstCom,decEstMax,decEstSax,decteriorMax] = decode_ufr_transformed_posterior(Trial,varargin)
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
    tag = DataHash({Trial.filebase,pfs.filename,sampleRate,units,ufr.spikeWindow});
end

filepath = fullfile(Trial.spath,[Trial.filebase,'.',mfilename(),'.',tag,'.mat']);

if exist(filepath,'file') && ~overwrite,    
    load(filepath);
else
% SET helper vars

    smoothingWeights = diag(smoothingWeights);
    pfsBins = pfs.adata.bins;
    ndims = numel(pfs.adata.bins);
    ufrLength = size(ufr,1);
    spikeWindow = ufr.spikeWindow;

% NEW VARS NEEDED
    xyz = {resample(preproc_xyz(Trial,'trb'),30)};
    hvec = cf(@(x)    x(:,'head_front',[1,2])-x(:,'head_back',[1,2]),                    xyz          );
    hvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                          hvec         );
    hvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                         hvec         );
    hvec = hvec{1};
    xyz = xyz{1};

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
        trm = pfs.plot(units(u),1,false,[],false,0.25,false,interpParPfs);
        ratemap = cat(2,ratemap,trm(logical(mask(:))));
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

    gbina = nan([size(ratemap,1),numel(pfsBins)]);
    for d = 1:ndims,
        gbina(:,d) = cat(2,binGrid{d}(:));
    end
    gbinPatch = gbina(sqrt(sum(gbina.^2,2))<300,:);
    posteriorPatch = nan([ufrLength,size(gbinPatch,1)]);
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
        E = exp(-sum(ratemap,2)*ones([1,bufferSize])*spikeWindow+log(ratemap)*(ufr.data(ind,:)-eps)');
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
        
        patchTF = bsxfun(@plus,...
                         sq(multiprod(permute(gbinPatch,[4,2,1,3]),hvec(ind,:,:),2,[2,3])),...
                         sq(xyz(ind,'hcom',[1,2])));
        
        pPatch = [];
        for tind = 1:bufferSize        
            Et = zeros([size(gbina,1),1]);
            Et(logical(mask(:))) = E(:,tind);
            
            F = scatteredInterpolant(gbinm(:,1),...
                                     gbinm(:,2),...
                                     E(:,tind),...
                                     'linear');
            pPatch(tind,:) = F(sq(patchTF(tind,1,:)),...
                               sq(patchTF(tind,2,:)));
        end

% $$$ hpos =     xyz(:,{'head_back','head_front'},[1,2]);
% $$$ hcom =     sq(xyz(:,{'hcom'},[1,2]));
% $$$ figure, 
% $$$ sp = gobjects([0,1]);
% $$$ sp(end+1) = subplot(121);
% $$$ sp(end+1) = subplot(122);
% $$$ for tind = 1:bufferSize,
% $$$     axes(sp(1));cla();hold('on');
% $$$     scatter(gbinm(:,1),gbinm(:,2),100,E(:,tind),'filled');
% $$$     circle(hcom(ind(tind),1),hcom(ind(tind),2),200,'r');
% $$$     line(hpos(ind(tind),:,1),hpos(ind(tind),:,2));
% $$$     axes(sp(2));cla();
% $$$     scatter(gbinPatch(:,1),gbinPatch(:,2),100,pPatch(tind,:),'filled');
% $$$     drawnow();
% $$$ end
        
% ACCUMULATE position estimates                
        decEstCom(ind,:)      = cpos;
        decEstSax(ind,:)      = tpos;
        decEstMax(ind,:)      = mpos;
        posteriorMax(ind)     = apos;
        posteriorPatch(ind,:) = pPatch;

        toc
    end% i

    
    % SAVE reconstruction statistics
    smoothingWeights = diag(smoothingWeights);    
    save(filepath,'decEstCom',   'decEstSax', 'decEstMax',...
                  'posteriorMax','sampleRate','spikeWindow','smoothingWeights',...
                  'gbinPatch','posteriorPatch');
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
