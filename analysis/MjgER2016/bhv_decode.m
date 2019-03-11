function [posEstCom,posEstMax,posEstSax,posteriorMax] = bhv_decode(Trial,varargin);
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
                 'ufr',                           [],                                            ...
                 'units',                         [],                                            ...
                 'mode',                          'xy',                                          ...
                 'pfsArgs',                       '',                                            ...
                 'interpParPfs',                  [],                                            ...
                 'spikeWindow',                   0.02,                                          ...
                 'overwrite',                     false                                          ...
);
[sampleRate,ufr,units,mode,pfsArgs,interpParPfs,spikeWindow,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


global MTA_PROJECT_PATH

hashTag = DataHash({Trial.filebase,mode,sampleRate,units,spikeWindow});


fileName = fullfile(Trial.spath,[Trial.filebase,'.decoded_',mode,'_',hashTag,'.mat']);

pitchSpaceDomain = linspace(  -2,  2,50);


if ~exist(fileName,'file') || overwrite,

    % COMPUTE / LOAD ND place fields
    if strcmp(mode,'xyhb'),
        xyz = preproc_xyz(Trial,'SPINE_SPLINE_HEAD_EQI');
        xyz.resample(16);
        fet = fet_HB_pitchB(Trial,16);
        xyzp = copy(fet);
        xyzp.data = cat(2,sq(xyz(:,'nose',[1,2])),xyzp.data);
        if isempty(pfsArgs),
            pfsArgs = struct('units',              units,                                ...
                             'states',             'theta-groom-sit',                    ...
                             'overwrite',          false,                                ...
                             'tag',                '',                                   ...
                             'binDims',            [ 100, 100,0.1,0.1],                  ...
                             'SmoothingWeights',   [0.8,0.8,2.5,2.5],                    ...
                             'type',               'xyhb',                               ...
                             'spkShuffle',         false,                                ...
                             'posShuffle',         false,                                ...
                             'numIter',            1,                                    ...
                             'xyzp',               xyzp,                                 ...
                             'boundaryLimits',     [-500,500;-500,500;-2,2;-2,2],        ...
                             'bootstrap',          false,                                ...
                             'halfsample',         false                                 ...
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
    elseif strcmp(mode,'xyh'),
        xyz = preproc_xyz(Trial,'SPINE_SPLINE_HEAD_EQI');
        xyz.resample(16);
        fet = fet_HB_pitchB(Trial,16);
        xyzp = copy(fet);
        xyzp.data = cat(2,sq(xyz(:,'nose',[1,2])),xyzp.data(:,1));
        if isempty(pfsArgs),
            pfsArgs = struct('units',              units,                                ...
                             'states',             'theta-groom-sit-rear',               ...
                             'overwrite',          false,                                ...
                             'tag',                '',                                   ...
                             'binDims',            [ 50, 50,0.1],                        ...
                             'SmoothingWeights',   [1.5,1.5,2],                          ...
                             'type',               'xyh',                                ...
                             'spkShuffle',         false,                                ...
                             'posShuffle',         false,                                ...
                             'numIter',            1,                                    ...
                             'xyzp',               xyzp,                                 ...
                             'boundaryLimits',     [-500,500;-500,500;-2,0.7],           ...
                             'bootstrap',          false,                                ...
                             'halfsample',         false                                 ...
                             );
        end
        
        pfsArgs = struct2varargin(pfsArgs);
        pfs = MTAApfs(Trial,pfsArgs{:});    
        % SET interp parameters                
        interpParPfs = struct('bins',{{linspace(-500,500,50),...
                                       linspace(-500,500,50),...
                                       linspace(  -2,  2,50)}},...
                              'nanMaskThreshold', 0.1,...
                              'methodNanMap',     'cubic',...
                              'methodRateMap',    'cubic');
    elseif strcmp(mode,'xyhi'),
        xyz = preproc_xyz(Trial,'SPINE_SPLINE_HEAD_EQI');
        xyz.resample(32);
        fet = fet_HB_pitch(Trial);
        fet.map_to_reference_session(Trial,'Ed05-20140529.ont.all');
        fet.resample(32);

        xyzp = copy(fet);
        xyzp.data = cat(2,sq(xyz(:,'nose',[1,2])),xyzp.data(:,3));
        if isempty(pfsArgs),
            pfsArgs = struct('units',              units,                                ...
                             'states',             'theta-groom-sit-rear',               ...
                             'overwrite',          false,                                ...
                             'tag',                'xyhi',                               ...
                             'binDims',            [ 100, 100,0.1],                      ...
                             'SmoothingWeights',   [1.5,1.5,1.1],                        ...
                             'type',               'xyh',                                ...
                             'spkShuffle',         false,                                ...
                             'posShuffle',         false,                                ...
                             'numIter',            1,                                    ...
                             'xyzp',               xyzp,                                 ...
                             'boundaryLimits',     [-500,500;-500,500;-2,2],             ...
                             'bootstrap',          false,                                ...
                             'halfsample',         false                                 ...
                             );
        end
        
        pfsArgs = struct2varargin(pfsArgs);
        %MTAApfs.purge_tagged_savefile(Trial,'xyhi');
        pfs = MTAApfs(Trial,pfsArgs{:});    
        % SET interp parameters                
        interpParPfs = struct('bins',{{linspace(-500,500,50),...
                                       linspace(-500,500,50),...
                                       linspace(  -2,  2,50)}},...
                              'nanMaskThreshold', 0.1,...
                              'methodNanMap',     'cubic',...
                              'methodRateMap',    'cubic');
        
    elseif strcmp(mode,'xyb'),
        xyz = preproc_xyz(Trial,'SPINE_SPLINE_HEAD_EQI');
        xyz.resample(16);
        fet = fet_HB_pitchB(Trial,16);
        xyzp = copy(fet);
        xyzp.data = cat(2,sq(xyz(:,'nose',[1,2])),xyzp.data(:,2));
        if isempty(pfsArgs),
            pfsArgs = struct('units',              units,                                ...
                             'states',             'theta-groom-sit',                    ...
                             'overwrite',          false,                                ...
                             'tag',                '',                                   ...
                             'binDims',            [ 100, 100,0.1],                    ...
                             'SmoothingWeights',   [0.8,0.8,2.5],                    ...
                             'type',               'xyb',                               ...
                             'spkShuffle',         false,                                ...
                             'posShuffle',         false,                                ...
                             'numIter',            1,                                 ...
                             'xyzp',               xyzp,                                 ...
                             'boundaryLimits',     [-500,500;-500,500;-2,2],        ...
                             'bootstrap',          false,                                ...
                             'halfsample',         false                                  ...
                             );
        end
        
        pfsArgs = struct2varargin(pfsArgs);
        pfs = MTAApfs(Trial,pfsArgs{:});    
        % SET interp parameters                
        interpParPfs = struct('bins',{{linspace(-500,500,50),...
                                       linspace(-500,500,50),...
                                       linspace(  -2,  2,50)}},...
                              'nanMaskThreshold', 0.1,...
                              'methodNanMap',     'linear',...
                              'methodRateMap',    'linear');
    elseif strcmp(mode,'xyz'),
        if isempty(pfsArgs),
            pfsArgs = struct('units',              units,                                ...
                             'states',             'theta-groom-sit',                    ...
                             'overwrite',          false,                                ...
                             'tag',                '',                                   ...
                             'binDims',            [ 100, 100, 30],                      ...
                             'SmoothingWeights',   [0.8,0.8,0.8],                        ...
                             'type',               'xyz',                                ...
                             'spkShuffle',         false,                                ...
                             'posShuffle',         false,                                ...
                             'numIter',            1,                                    ...
                             'xyzp',               [],                                   ...
                             'boundaryLimits',     [-500,500;-500,500;-60,360],          ...
                             'bootstrap',          false,                                ...
                             'halfsample',         false                                 ...
                             );
        end
        
        pfsArgs = struct2varargin(pfsArgs);
        pfs = MTAApfs(Trial,pfsArgs{:});    
        % SET interp parameters                
        interpParPfs = struct('bins',{{linspace(-500,500,50),...
                                       linspace(-500,500,50),...
                                       linspace( -60,360,60)}},...
                              'nanMaskThreshold', 0.1,...
                              'methodNanMap',     'cubic',...
                              'methodRateMap',    'cubic');
        
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
    if isempty(ufr),
        ufr = load(Trial,'ufr',xyz,spk,units,spikeWindow,true,'gauss');
    end
    
    %fet = fet_HB_pitchB(Trial,sampleRate);

% CREATE spatial mask
    width = numel(pfsBins{1});
    height = numel(pfsBins{2});
    radius = round(numel(pfsBins{1})/2)-find(pfsBins{1}<-440,1,'last');
    centerW = width/2;
    centerH = height/2;
    [W,H] = meshgrid(1:width,1:height);           
    circMask = repmat(double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius),[1,1,pfsBinsDims(3:4)]);
    if strcmp(mode,'xyhb')||strcmp(mode,'xyh')||strcmp(mode,'xyb'),
% CREATE behavioral mask
% LOAD bhv analysis vars
        ds = load(fullfile(MTA_PROJECT_PATH,'analysis','req20180123_pfd_erpPCA-HBPITCHxBPITCH_v7.mat'));
        if ~exist('pfb','var'),  pfb = MTAApfs(Trial,'tag',['HBPITCHxBPITCH_v7']);  end
        interpParPfb = struct('bins',{ {pitchSpaceDomain,...
                            pitchSpaceDomain}},...
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
        SmoothingWeights = [0.2,0.2];
        for i = 1:ndims(sinterpGrids),
            sinterpGrids{i} = sinterpGrids{i}.^2/SmoothingWeights(i)^2/2;
        end
        Smoother = exp(sum(-cat(ndims(sinterpGrids)+1,sinterpGrids{:}),ndims(sinterpGrids)+1));
        Smoother = Smoother./sum(Smoother(:));
        bhvMask   = convn(bhvMask, Smoother,'same');
        bhvMask = double(bhvMask>0.01);
        if strcmp(mode,'xyhb')
            bhvMask = repmat(permute(bhvMask,[3,4,1,2]),[pfsBinsDims(1:2),1,1]);
        elseif strcmp(mode,'xyh')
            bhvMask = repmat(permute(~all(bhvMask==0,2),[3,4,1,2]),[pfsBinsDims(1:2),1,1]);
        elseif strcmp(mode,'xyb')            
            bhvMask = repmat(permute(~all(bhvMask==0,2),[3,4,1,2]),[pfsBinsDims(1:2),1,1]);
            %bhvMask = repmat(permute(~all(bhvMask==0),[3,4,2,1]),[pfsBinsDims(1:2),1,1]);
        end
    elseif strcmp(mode,'xyz'),
        bhvMask = pfsBins{3}>0&pfsBins{3}<300;
        bhvMask = repmat(permute(bhvMask,[3,4,2,1]),[pfsBinsDims(1:2),1,1]);        
            
        
    else
        bhvMask = ones(pfsBinsDims);
    end

% CREATE rate map mask
    mazeMask = logical(circMask.*bhvMask);

    rateMap = zeros([sum(mazeMask(:)),0]);
    for u = 1:numel(units),
        trm = pfs.plot(units(u),1,false,[],false,0.25,false,interpParPfs);
        rateMap = cat(2,rateMap,trm(mazeMask(:)));
    end
    clear('bhvMask','circMask','trm');
    

% COMPUTE Posterior distribution base on ratemaps and unit firing rates
    rateMap(isnan(rateMap)) = 0;
    rateMap = rateMap+1e-3;    

    posEstCom = nan([size(xyz,1),numel(pfsBins)]);
    posEstMax = nan([size(xyz,1),numel(pfsBins)]);
    posEstSax = nan([size(xyz,1),numel(pfsBins)]);    
    posteriorMax = nan([size(xyz,1),1]);

    numDimPos = numel(pfsBins);

    binGrid = cell([1,numel(pfsBins)]); 
    [binGrid{:}] = ndgrid(pfsBins{:});

    gbinm = nan([size(rateMap,1),2]);
    for d = 1:numDimPos,
        gbinm(:,d) = cat(2,binGrid{d}(mazeMask(:)));
    end
    clear('binGrid');
    
    switch mode
      case 'xyhb'
        bufferSize = 2^11;
        posteriorThresh = 0.00001;        
        A = [250.^2,       0,      0,      0;...
                  0,  250.^2,      0,      0;...
                  0,       0, 0.4.^2,      0;...
                  0,       0,      0, 0.4.^2];
      case 'xyh'
        bufferSize = 2^16;
        posteriorThresh = 0.0001;
        A = [250.^2,     0,      0;...
                0,  250.^2,      0;...
                0,       0, 0.4.^2];

      case 'xyhi'
        bufferSize = 2^16;
        posteriorThresh = 0.0001;
        A = [250.^2,     0,      0;...
                0,  250.^2,      0;...
                0,       0, 0.4.^2];
        
      case 'xyb'
        bufferSize = 2^16;
        posteriorThresh = 0.0001;        
        A = [250.^2,       0,      0;...
                  0,  250.^2,      0;...
                  0,       0, 0.4.^2];
        
      case 'xyz'
        bufferSize = 2^14;
        posteriorThresh = 0.00001;
        A = [200.^2,      0,     0;...
                  0, 200.^2,     0;...
                  0,      0, 20.^2];
        
      case 'xy'
        bufferSize = 2^16;
        posteriorThresh = 0.001;        
        A = [250.^2,      0;...
                  0, 250.^2];
        g = fittype( @(A,xa,ya,xya,xo,yo,x,y)                                   ...
                     A.*exp(-(xa.*(x-xo).^2+xya.*(x-xo).*(y-yo)+ya.*(y-yo).^2)),...
                     'independent',{'x', 'y'},'dependent', 'z' ); 
        
      otherwise
        bufferSize = 2^10;
        posteriorThresh = 0.001;        
    end
    
    mpos = nan([bufferSize,numDimPos]);
    apos = nan([bufferSize,1]);            



    
    for i = 1:ceil(size(ufr,1)/bufferSize)
        tic
        disp([num2str(i),' of ' num2str(ceil(size(ufr,1)/bufferSize))]);
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
        E = exp(-sum(rateMap,2)*ones([1,bufferSize])*spikeWindow+log(rateMap)*(ufr.data(ind,:)-eps)');
        E = bsxfun(@rdivide,E,sum(E));

        [apos,tbin] = max(E);
        mpos = gbinm(tbin,:);  

% $$$         wbinm = bsxfun(@minus,repmat(permute(gbinm,[3,2,1]),[size(tbin,2),1,1]),gbinm(tbin',:));
% $$$         weights = exp(multiprod(-wbinm,multiprod(inv(A),wbinm,[1,2],[2]),[2],[2]));        
% $$$                 
        tpos = nan([bufferSize,numDimPos]);
        cpos = nan([bufferSize,numDimPos]);        
        for tind = 1:bufferSize
            if apos(tind)<posteriorThresh,continue;end;
            wbinm = bsxfun(@minus,gbinm,gbinm(tbin(tind),:));                
            weights = exp(multiprod(-wbinm,multiprod(inv(A),wbinm,[1,2],[2]),[2],[2]));
            weights = weights./sum(weights);
            weights = bsxfun(@times,weights,repmat(E(:,tind),[1,numDimPos]));
            weights = weights./sum(weights);

            tpos(tind,:) = sum(gbinm.*weights,'omitnan');
            cpos(tind,:) = sum(gbinm.*repmat(E(:,tind),[1,numDimPos]),'omitnan');            
        end


        
        
% ACCUMULATE position estimates                
        posEstCom(ind,:)      = cpos;
        posEstSax(ind,:)      = tpos;
        posEstMax(ind,:)      = mpos;
        posteriorMax(ind)     = apos;

        toc
    end

    % SAVE reconstruction statistics
    save(fileName,'posEstCom','posEstSax','posEstMax','posteriorMax','sampleRate','spikeWindow');

else
    load(fileName);
end




