function [posEstCom,posEstMax,posteriorMax] = bhv_decode(Trial,varargin);
% function [posEstCom,posEstMax,posteriorMax] = bhv_decode(Trial,varargin);
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
%    posEstCom
%    posEstMax
%    posteriorMax
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('sampleRate',                    40,                                            ...
                 'units',                         [],                                            ...
                 'mode',                          'xy',                                          ...
                 'pfsArgs',                       '',                                            ...
                 'interpParPfs',                  [],                                            ...
                 'spikeWindow',                   0.02,                                          ...
                 'overwrite',                     false                                          ...
);
[sampleRate,units,mode,pfsArgs,interpParPfs,spikeWindow,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


global MTA_PROJECT_PATH
addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/fwdrebayesiandecodingtools/
addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/tempScripts/
addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/aux/

hashTag = DataHash({Trial.filebase,mode,sampleRate,units,spikeWindow});


fileName = fullfile(Trial.spath,[Trial.filebase,'.decoded_',mode,'_',hashTag,'.mat']);

if ~exist(fileName,'file') || overwrite,

    % COMPUTE / LOAD ND place fields
    if strcmp(mode,'xyhb'),
        xyz = preproc_xyz(Trial,'SPINE_SPLINE_HEAD_EQI');
        xyz.resample(16);
        fet = fet_HB_pitchB(Trial,16);
        xyzp = copy(fet);
        xyzp.data = cat(2,sq(xyz(:,'nose',[1,2])),xyzp.data);
        if ~isempty(pfsArgs),
            pfsArgs = struct('units',              units,                                ...
                             'states',             'theta-groom-sit',                    ...
                             'overwrite',          false,                                ...
                             'tag',                '',                                   ...
                             'binDims',            [ 100, 100,0.4,0.4],                  ...
                             'SmoothingWeights',   [0.8,0.8,0.8,0.8],                    ...
                             'type',               'xyhb',                               ...
                             'spkShuffle',         false,                                ...
                             'posShuffle',         false,                                ...
                             'numIter',            1000,                                 ...
                             'xyzp',               xyzp,                                 ...
                             'boundaryLimits',     [-500,500;-500,500;-2,2;-2,2],        ...
                             'bootstrap',          false,                                ...
                             'halfsample',         true                                  ...
                             );
        end
        
        pfsArgs = struct2varargin(pfsArgs);
        pfs = MTAApfs(Trial,pfsArgs{:});    
        % SET interp parameters                
        interpParPfs = struct('bins',{{linspace(-500,500,50),...
                                       linspace(-500,500,50),...
                                       linspace(  -2,  2,50),...
                                       linspace(  -2,  2,50)}},...
                              'nanMaskThreshold', 0.1,...
                              'methodNanMap',     'linear',...
                              'methodRateMap',    'linear');
        
    elseif strcmp(mode,'xy')
        % COMPUTE / LOAD 2d place fields
        pfs = pfs_2d_theta(Trial,units);
        % SET interp parameters
        interpParPfs = struct('bins',{{linspace(-500,500,100),...
                                       linspace(-500,500,100)}},...
                              'nanMaskThreshold', 0.1,...
                              'methodNanMap',     'linear',...
                              'methodRateMap',    'linear');
    end


    % SET interpolation parameters
    pfsBins = cell([1,4]);
    if ~isempty(interpParPfs),
        pfsBins(1:numel(interpParPfs.bins)) = interpParPfs.bins;
    else,
        pfsBins = pfs.adata.bins;
    end
    pfsBinsDims = cellfun(@numel,pfsBins);
    pfsBinsDims(pfsBinsDims==0) = 1;
    pfsBins(cellfun(@isempty,pfsBins)) = [];


    % LOAD xyz data
    xyz = filter(resample(load(Trial,'xyz'),sampleRate),'ButFilter',3,3);
    ufr = load(Trial,'ufr',xyz,[],units,spikeWindow,true,'gauss');


    % CREATE spatial mask
    width = numel(pfsBins{1});
    height = numel(pfsBins{2});
    radius = round(numel(pfsBins{1})/2)-find(pfsBins{1}<-440,1,'last');
    centerW = width/2;
    centerH = height/2;
    [W,H] = meshgrid(1:width,1:height);           
    circMask = repmat(double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius),[1,1,pfsBinsDims(3:4)]);
    if strcmp(mode,'xyhb'),
        % CREATE behavioral mask
        % LOAD bhv analysis vars
        ds = load(fullfile(MTA_PROJECT_PATH,'analysis','req20180123_pfd_erpPCA-HBPITCHxBPITCH_v7.mat'));
        if ~exist('pfb','var'),  pfb = MTAApfs(Trial,'tag',['HBPITCHxBPITCH_v7']);  end
        interpParPfb = struct('bins',{{pfsBins{3},...
                            pfsBins{4}}},...
                              'nanMaskThreshold', 0.01,...
                              'methodNanMap',     'linear',...
                              'methodRateMap',    'linear');
        bhvMask = zeros([cellfun(@numel,pfb.adata.bins)]);
        bhvMask(ds.vDims) = 1;
        interpGrids = cell([1,numel(interpParPfb.bins)]);
        [interpGrids{:}] = ndgrid(interpParPfb.bins{:});
        bhvMask = interpn(pfb.adata.bins{:},bhvMask,interpGrids{:},interpParPfb.methodRateMap);
        bhvMask(isnan(bhvMask)) = 0;
        sinterpGrids = interpGrids;
        SmoothingWeights = [0.1,0.1];
        for i = 1:ndims(sinterpGrids),
            sinterpGrids{i} = sinterpGrids{i}.^2/SmoothingWeights(i)^2/2;
        end
        Smoother = exp(sum(-cat(ndims(sinterpGrids)+1,sinterpGrids{:}),ndims(sinterpGrids)+1));
        Smoother = Smoother./sum(Smoother(:));
        bhvMask   = convn(bhvMask, Smoother,'same');
        bhvMask = double(bhvMask>0.05);
        bhvMask = repmat(permute(bhvMask,[3,4,1,2]),...
                         [pfsBinsDims(1:2),1,1]);
    else
        bhvMask = 1;    
    end

    % CREATE rate map mask
    mazeMask = logical(circMask.*bhvMask);
    %mazeMask(mazeMask==0)=nan;

    % ACCUMULATE rate maps for bayesian decoding
% $$$     rateMap = zeros([cell2mat(cf(@numel,pfsBins)),0]);
% $$$     for u = 1:numel(units),
% $$$         trm = pfs.plot(units(u),'mean',false,[],false,0.25,false,interpParPfs,[],mazeMask);
% $$$         rateMap = cat(numel(pfsBins)+1,rateMap,trm);
% $$$     end
% $$$     clear('mazeMask','bhvMask','circMask','trm');
    rateMap = zeros([sum(mazeMask(:)),0]);
    for u = 1:numel(units),
        trm = pfs.plot(units(u),'mean',false,[],false,0.25,false,interpParPfs);
        rateMap = cat(2,rateMap,trm(logical(mazeMask(:))));
    end
    clear('bhvMask','circMask','trm');
    

    % REDUCE prior dimensionality
% $$$     mapSize = size(rateMap);
% $$$     erows = {};
% $$$     grows = {};
% $$$     eeCount = [];
% $$$     for i = 1:ndims(rateMap)-1;
% $$$         eeCount(:,i) = sum(reshape(permute(isnan(rateMap),[i,find(~ismember(1:ndims(rateMap),i))]),mapSize(i),[]),2);
% $$$         erows{i}(:) = eeCount(:,i)==prod(mapSize(~ismember(1:ndims(rateMap),i)));
% $$$         grows{i}(:) = eeCount(:,i)~=prod(mapSize(~ismember(1:ndims(rateMap),i)));
% $$$     end
% $$$     grows{end+1} = ones([1,numel(units)]);
% $$$     erows{end+1} = ones([1,numel(units)]);
% $$$ 
% $$$     gpfsBins = {};
% $$$     for i = 1:ndims(rateMap)-1;
% $$$         gpfsBins{i} = pfsBins{i}(grows{i});
% $$$     end
% $$$ 
% $$$     % REDUCE the dimensions of rateMap to save memory 
% $$$     smap = rateMap(grows{:});
% $$$     smap = rateMap;
% $$$     for i = 1:ndims(rateMap)-1,
% $$$         smap(erows{i},:,:,:,:) = [];
% $$$         smap = permute(smap,[2:ndims(rateMap),1]);
% $$$     end
% $$$     smap = permute(smap,[2:ndims(rateMap),1]);
    smap = rateMap;

    %clear('rateMap');

    % COMPUTE Posterior distribution base on ratemaps and unit firing rates
    smap(isnan(smap)) = 0;

% $$$     posEstCom = nan([size(xyz,1),ndims(smap)-1]);
% $$$     posEstMax = nan([size(xyz,1),ndims(smap)-1]);
% $$$     posteriorMax = nan([size(xyz,1),1]);
    posEstCom = nan([size(xyz,1),numel(pfsBins)]);
    posEstMax = nan([size(xyz,1),numel(pfsBins)]);
    posteriorMax = nan([size(xyz,1),1]);

    numDimPos = numel(pfsBins);

% $$$     binsE = {};
% $$$     for bin = 1:numel(pfs.adata.bins),  
% $$$         binsE{bin} = pfsBins{bin}(grows{bin});  
% $$$     end
% $$$     gbins = cell([1,numel(pfsBins)]); 
% $$$     [gbins{:}] = ndgrid(binsE{:});
% $$$     gbinm = cat(numDimPos+1,gbins{:});

    binGrid = cell([1,numel(pfsBins)]); 
    [binGrid{:}] = ndgrid(pfsBins{:});

    gbinm = nan([size(smap,1),2]);
    for d = 1:numDimPos,
        gbinm(:,d) = cat(2,binGrid{d}(mazeMask(:)));
    end
    
    bufferSize = 10000;
    for i = 1:ceil(size(ufr,1)/bufferSize)
        tic
        disp([num2str(i),' of ' num2str(ceil(size(ufr,1)/bufferSize))]);
        if i == ceil(size(ufr,1)/bufferSize),
            ind = (1+bufferSize*(i-1)) : (bufferSize*(i-1)+mod(size(ufr,1),bufferSize));
            bufferSize = numel(ind);                        
        else        
            ind = (1+bufferSize*(i-1)) : (bufferSize*i);
        end
        
% $$$         subind = nniz(xyz.data(ind,1,1));
% $$$         nsubsamp = sum(subind);
        
% $$$         if nsubsamp>0,
% $$$             tE = decode_bayesian_poisson(smap,ufr.data(ind(subind),:)'-eps,'bins',ones([1,bufferSize])/250);
            tE = decode_bayesian_poisson(smap,ufr.data(ind,:)'-eps,'bins',ones([1,bufferSize])/sampleRate);
% $$$             tESize = size(tE);

            % ESTIMATE positions from posterior center of mass
            % ESTIMATE positions from posterior max        
% $$$             ss = substruct('()',[repmat({':'},[1,numel(gbins)]),{1}]);
            mpos = nan([bufferSize,numDimPos]);
            tpos = nan([bufferSize,numDimPos]);
            apos = nan([bufferSize,1]);            
% $$$             mpos = nan([nsubsamp,numDimPos]);
% $$$             apos = nan([nsubsamp,1]);
% $$$             tpos = nan([nsubsamp,numDimPos]);
            for tind = 1:size(tE,ndims(tE)),
% $$$                 ss.subs{end} = tind;
% $$$                 [tval,tbin] = max(reshape(subsref(tE,ss),[],1));
                [tval,tbin] = max(tE(:,tind));
                if ~isempty(tbin),        
% $$$                     maxInd       = cell([1,numDimPos]);
% $$$                     [maxInd{:}]  = ind2sub(tESize(1:numDimPos),tbin);
% $$$                     mpos(tind,:) = cell2mat(cf(@(gbins,mbin) gbins(mbin), gpfsBins,maxInd));
% $$$                     apos(tind)   = tE(maxInd{:},tind);
                    mpos(tind,:) = gbinm(tbin,:);
                    apos(tind) = tval;
                    tpos(tind,:) = sum(gbinm.*repmat(tE(:,tind),[1,numDimPos]));
                end
                
% $$$                 tpos(tind,:) = sum(reshape(gbinm.*repmat(subsref(tE,ss),[ones([1,numDimPos]),numDimPos]),[],numDimPos));

            end
            
            % ACCUMULATE position estimates
            posEstCom(ind,:)      = tpos;
            posEstMax(ind,:)      = mpos;
            posteriorMax(ind)     = apos;
% $$$             posEstCom(ind(subind),:)      = tpos;
% $$$             posEstMax(ind(subind),:)      = mpos;
% $$$             posteriorMax(ind(subind))     = apos;
% $$$         end

        toc
    end

    % SAVE reconstruction statistics
    save(fileName,'posEstCom','posEstMax','posteriorMax','sampleRate','spikeWindow');

else
    load(fileName);
end



