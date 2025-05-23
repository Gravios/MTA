classdef MTAApfs < hgsetget %< MTAAnalysis
% function Pfs = MTAApfs(Obj,varargin)            
%
% Computes the rate map conditioned on the provided covariates
%
% 
%
% Obj (MTATrial/MTASession/MTAApfs):
%
%    MTATrial & MTASession -> Create rate maps for neural units
%
%    MTAApfs -> update object
%
%    otherwise -> return empty MTAApfs object
%
%-----------------------------------------------------------------------------
%  varargin:
% 
%    units      -  Array,        list of units
%    states     -  String,       computation periods
%    overwrite  -  Logical,      overwrite rate maps
%    tag        -  String,       optional string for filename 
%    binDims    -  Array,        dimensions of ratemap bins (mm)
%    SmoothingWeights - Array,   std of gaussian kernel (bins)
%    type       -  String,       dimension discriptive characters (e.g. xy, xyz)
%    spkShuffle -  Logical,      FLAG, shuffle spiketimes  randomly
%    posShuffle -  Logical,      FLAG, shuflle position (blocks) randomly
%    numIter    -  Integer,      number of bootstrap iterations
%    xyzp       -  MTAData,      Alternative data object for ratemap space
%    boundaryLimits - Matrix,    Computational boundaries (default->Trial.maze.boundaries)
%    bootstrap  -  Logical,      FLAG, Compute bootstrap ratemaps
%    halfsample -  Logical,      FLAG, Compute each bootstrap with half the data
%    shuffleBlockSize - Double,  Data block size in seconds
%    trackingMarker   - String,  Marker name within MTADxyz object 
%                                (ignored if xyzp not empty)
%    autoSaveFlag - Logical,     FLAG save between unit computations
%    spkMode      - Logical,     FLAG remove burst (ignored if spk is provided)
%    spk          - MTASpk,      contains unit res and clu
%    compute_pfs  - function,    function handle to computational core
%    loadMethod   - RESERVED FOR FUTURE USE
%    posShuffleDims - Array(Numeric), specify shuffling dimentions
%
%-----------------------------------------------------------------------------
%  Description
%    Create an object which contains and manages rate maps,which
%    are spatially binned expected rates for individual neurons
%
%-----------------------------------------------------------------------------
%  Output:
%    Pfs (MTAApfs object) 
%  
%-----------------------------------------------------------------------------
%  Examples:    
%    pfs = MTAApfs(Trial,units,states);
%
%-----------------------------------------------------------------------------
%  NOTES:
%    hash - future version will include hash mod tracking
%    states - do not use TimeSeries yet
%
%-----------------------------------------------------------------------------
%  UPDATES:
%   X 20180515 (Justin Graboski) change preprocessing to always include first
%              steps for partitioning the data into disjoint blocks
%
        
    properties 
        %path - string: location of file containing MTAApfs object
        path = '';
        
        %filename - string: location of file containing MTAApfs object
        filename = '';
        
        %session - struct(sessionName; mazeName; trialName): requried to load parent session
        session = struct('sessionName',           '',                                            ...
                         'trialName',             '',                                            ...
                         'mazeName',              ''                                             ...
        );
        
        %tag - string: unique identifier
        tag = '';
        
        %ext - string: file extention
        ext = 'pfs';
        
        %parameters - struct:
        parameters
        
        %mdata - struct: metadata
        mdata
        
        %adata - struct: auxdata 
        adata = struct('trackingMarker',          'hcom',                                        ...
                       'bins',                    {{}},                                          ...
                       'binSizes',                [],                                            ...
                       'mask',                    []                                             ...
        );

        %data - struct: placefield data
        data = struct( 'clu',                     [],                                            ...
                       'elClu',                   [],                                            ...
                       'el',                      [],                                            ...
                       'maxRate',                 [],                                            ...
                       'maxRateInd',              [],                                            ...
                       'maxRatePos',              [],                                            ...
                       'rateMap',                 [],                                            ...
                       'meanRate',                [],                                            ...
                       'si',                      [],                                            ...
                       'spar',                    []                                             ...
        );
        
        % hash - future version will include hash mod tracking
    end

    methods

        function Pfs = MTAApfs(Obj, varargin)     
            if ischar(Obj)||isstruct(Obj),
               Obj = MTATrial.validate(Obj);
            end
% DEFARGS ------------------------------------------------------------------------------------------ 
            defargs = struct(                                                                    ...
                'units',                          [],                                            ...
                'states',                         'walk',                                        ...
                'overwrite',                      0,                                             ...
                'tag',                            [],                                            ...
                'binDims',                        [30,30],                                       ...
                'SmoothingWeights',               [],                                            ...
                'type',                           'xy',                                          ...
                'spkShuffle',                     0,                                             ...
                'posShuffle',                     0,                                             ...
                'numIter',                        1,                                             ...
                'xyzp',                           MTADxyz([]),                                   ...
                'boundaryLimits',                 [],                                            ...
                'bootstrap',                      0,                                             ...
                'halfsample',                     0,                                             ...
                'shuffleBlockSize',               1,                                             ...
                'trackingMarker',                 'hcom',                                        ...
                'autoSaveFlag',                   true,                                          ...
                'spkMode',                        'deburst',                                     ...
                'spk',                            [],                                            ...
                'compute_pfs',                    @PlotPF,                                       ...
                'loadMethod',                     '',                                            ...
                'posShuffleDims',                 [],                                            ...
                'purge',                          false                                           ...
            );            
            [units,states,overwrite,tag,binDims,SmoothingWeights,type,spkShuffle,posShuffle,     ...
             numIter,xyzp,boundaryLimits,bootstrap,halfsample,shuffleBlockSize,trackingMarker,   ...
             autoSaveFlag,spkMode,spk,compute_pfs,loadMethod,posShuffleDims, purge] = ...
                DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------
                


            if isa(Obj,'MTASession')
                Session = Obj;
                    
% $$$                     if ~isempty(tag),
% $$$                         filepath = fullfile(Session.spath,[Session.filebase,'.Pfs.' tag '.mat']);
% $$$                         if exist(filepath),
% $$$                             load(filepath);
% $$$                             if all(ismember(units,Pfs.data.clu)),
% $$$                                 return
% $$$                             else
% $$$                                 continue
% $$$                             end
% $$$                         end                        
% $$$                     end
                    

% SET path
                Pfs.path                  = Session.spath;
                Pfs.tag                   = tag;
% SET session
                Pfs.session.sessionName   = Session.name;
                Pfs.session.trialName     = Session.trialName;
                Pfs.session.mazeName      = Session.maze.name;
% SET parameters
                Pfs.parameters.binDims    = binDims;                    
                Pfs.parameters.type       = type;
                Pfs.parameters.spkMode    = spkMode;                    
                Pfs.parameters.spkShuffle = spkShuffle;
                Pfs.parameters.posShuffle = posShuffle;
                Pfs.parameters.numIter    = numIter;
                Pfs.parameters.bootstrap  = bootstrap;
                Pfs.parameters.smoothingWeights = SmoothingWeights;
                Pfs.parameters.states     = '';
                
                Pfs.adata.trackingMarker  = trackingMarker;
                    
% SET the epochs
                if ischar(states),
                    Pfs.parameters.states = states;
                    pfsState = Session.stc{states};
                elseif isa(states,'MTAData'),
                    pfsState = states.copy;
                    Pfs.parameters.states = pfsState.label;                       
                end
% ASSIGN computational boundaries
                if isempty(boundaryLimits),
                    boundaryLimits = Session.maze.boundaries(1:numel(type),:);
                end
                assert(all(diff(boundaryLimits,1,2)>0),'MTA:MTAApfs:BoundaryLimitRangeNegative');
% COMPUTE bin count based on bin dimension and computational boundaries
                Pfs.adata.binSizes = round(abs(diff(boundaryLimits(1:numel(binDims),:),1,2))./binDims');
                
                Pfs.update_filename(Session,tag);
                if exist(fpath(Pfs),'file') && ~overwrite
% LOAD place fields
                    disp(Pfs.fpath);
                    if purge,
                        delete(Pfs.fpath);
                    end
                    load(Pfs.fpath);
                    if all(ismember(units,Pfs.data.clu)),
                        return;
                    end
                end
                    
% SET the independent variabels
                if isempty(xyzp),
                    try,
                        xyz = preproc_xyz(Session,'trb');
                    catch err,
                        disp(err)                        
                        try,
                            xyz = preproc_xyz(Session,'crb');
                            %xyz = preproc_xyz(Session);
                            trackingMarker = 'hcom';
                        catch
                            xyz = Session.load('xyz');
                            trackingMarker = 'hcom';
                        end
                    end
                    xyz.resample(16);
                    xyz.data = sq(xyz(:,trackingMarker,1:numel(binDims)));
                else
                    xyz = xyzp;
                end
% RESAMPLE state to new sampleRate
                resample(pfsState,xyz);                    
                
                
            elseif isa(Obj,'MTAApfs'),
                Pfs = Obj;
                Session = MTATrial(Obj.session.sessionName,...
                                   Obj.session.mazeName,...
                                   Obj.session.trialName);
% RESET the independent variabels
                if xyzp.isempty,
                    try,
                        xyz = preproc_xyz(Session,'trb');
                    catch err,
                        disp(err)
                        xyz = preproc_xyz(Session);
                        trackingMarker = 'hcom';
                    end
                    xyz.resample(16);
                    xyz.data = sq(xyz(:,trackingMarker,1:numel(binDims)));
                else
                    xyz = xyzp;
                end
                    
% RESET the epochs
                if ischar(states),
                    Pfs.parameters.states = states;
                    pfsState = Session.stc{states};
                    resample(pfsState,xyz);
                elseif isa(states,'MTAData'),
                    pfsState = states.copy;
                    resample(pfsState,xyz);
                    Pfs.parameters.states = pfsState.label;                       
                end
% ASSIGN computational boundaries
                if isempty(boundaryLimits),
                    boundaryLimits = Session.maze.boundaries(1:numel(type),:);
                end
% COMPUTE bin count based on bin dimension and computational boundaries
                Pfs.adata.binSizes = round(abs(diff(boundaryLimits(1:numel(binDims),:),1,2))./binDims');
            else,
                return
            end%if isa
                        
            numUnits = numel(units);
            selected_units = units;

            if purge,
                delete(Pfs.fpath);
            end
            
            pf_tmpfile = Pfs.fpath;

% LOAD existing data
            epftmp = exist(pf_tmpfile,'file');

            if epftmp && ~overwrite,
                if isempty(Pfs.data.clu),
                    load(pf_tmpfile);
                end
                % Load specific units
                if isempty(units)
                    return
                else
                    unprocessed_units = ~ismember(units,Pfs.data.clu);
                    if sum(unprocessed_units)==0
% LOAD and RETURN 
                        selected_units_ind = ismember(Pfs.data.clu,units);
                        field = fieldnames(Pfs.data);
                        for f = 1:numel(field);
                            Pfs.data.(field{f}) = Pfs.data.(field{f})(:,selected_units_ind,:,:,:);
                        end
                        return
                    else
% ADD entrys in Pfs.data for uncomputed units
                        oldUnitInds = find(ismember(Pfs.data.clu,units));
                        numOldUnits = numel(oldUnitInds);
                        numNewUnits = numUnits - numOldUnits;
                        tnumUnits = numNewUnits + numel(Pfs.data.clu);
                        if numNewUnits>0,
                            newdata =struct( 'clu',        zeros([1,numNewUnits]),...
                                'elClu',      zeros([1,numNewUnits]),...
                                'el',         zeros([1,numNewUnits]),...
                                'maxRate',    zeros([5,numNewUnits]),...
                                'maxRateInd', zeros([5,numNewUnits,numel(binDims)]),...
                                'maxRatePos', zeros([5,numNewUnits,numel(binDims)]),...
                                'rateMap',    zeros([prod(Pfs.adata.binSizes),numNewUnits,numIter]),...
                                'meanRate',   zeros([1,numNewUnits]),...
                                'si',         zeros([1,numNewUnits]),...
                                'spar',       zeros([1,numNewUnits]));
                            field = fieldnames(newdata);
                            for f = 1:numel(field);
                                Pfs.data.(field{f}) = cat(2,Pfs.data.(field{f}),newdata.(field{f}));
                            end
                        end
                        selected_units = units(~ismember(units,Pfs.data.clu));
                        dind = [tnumUnits-numNewUnits+1:tnumUnits];
                    end
                end

            elseif epftmp && overwrite,
            %% Extend Pfs data for additional units
                numUnits = numel(units);
                if numUnits ==0
                    numUnits = size(Session.spk.map,1);
                    units    = 1:size(Session.spk.map,1);
                end
                if isempty(Pfs.data.clu),
                    load(pf_tmpfile);
                end
                oldUnitInds = find(ismember(Pfs.data.clu,units));
                numOldUnits = numel(oldUnitInds);
                numNewUnits = numUnits - numOldUnits;
                tnumUnits = numNewUnits + numel(Pfs.data.clu);
                if numNewUnits>0,
                    newdata =struct( 'clu',        zeros([1,numNewUnits]),...
                                     'elClu',      zeros([1,numNewUnits]),...
                                     'el',         zeros([1,numNewUnits]),...
                                     'maxRate',    zeros([5,numNewUnits]),...
                                     'maxRateInd', zeros([5,numNewUnits,numel(binDims)]),...
                                     'maxRatePos', zeros([5,numNewUnits,numel(binDims)]),...
                                     'rateMap',    zeros([prod(Pfs.adata.binSizes),numNewUnits,numIter]),...
                                     'meanRate',   zeros([1,numNewUnits]),...
                                     'si',         zeros([1,numNewUnits]),...
                                     'spar',       zeros([1,numNewUnits]));
                    field = fieldnames(newdata);
                    for f = 1:numel(field);
                        Pfs.data.(field{f}) = cat(2,Pfs.data.(field{f}),newdata.(field{f}));
                    end
                end
                selected_units = [units(ismember(units,Pfs.data.clu)),units(~ismember(units,Pfs.data.clu))];
                dind = [oldUnitInds(:)',[tnumUnits-numNewUnits+1:tnumUnits]];                
            elseif ~epftmp
            %% Instantiate Pfs Data Variables if DeNovoCalc
                numUnits = numel(units);
                if numUnits ==0
                    numUnits = size(Session.spk.map,1);
                    units    = 1:size(Session.spk.map,1);
                    selected_units = units;                    
                end
                dind = 1:numUnits;
                Pfs.data =struct('clu',        zeros([1,numUnits]),...
                                 'elClu',      zeros([1,numUnits]),...
                                 'el',         zeros([1,numUnits]),...
                                 'maxRate',    zeros([5,numUnits]),...
                                 'maxRateInd', zeros([5,numUnits,numel(binDims)]),...
                                 'maxRatePos', zeros([5,numUnits,numel(binDims)]),...
                                 'rateMap',    zeros([prod(Pfs.adata.binSizes),numUnits,numIter]),...
                                 'meanRate',   zeros([1,numUnits]),...
                                 'si',         zeros([1,numUnits]),...
                                 'spar',       zeros([1,numUnits]));

            end
            
            
% GET State Positions

            returnFlag = false;
            if ~isempty(pfsState),
                switch pfsState.type
                  case 'TimePeriods'
                    asstpos = sq(xyz(pfsState,:));
                  case 'TimeSeries'
                    asstpos = sq(xyz(logical(pfsState(:)),:));
                end
                
% LOAD Units into spk object;
                if isempty(spk)
                    spk = Session.spk.copy;
                    spk.create(Session, xyz.sampleRate, pfsState, ...
                               units, spkMode);
                end
            else
                spk = Session.spk.copy;
                returnFlag = true;
            end

            

            i = 1;
            for unit=selected_units(:)',

                Pfs.data.clu(dind(i))            = spk.map(unit,1);
                Pfs.data.el(dind(i))             = spk.map(unit,2);
                Pfs.data.elClu(dind(i))          = spk.map(unit,3);
                Pfs.data.maxRate(:,dind(i))      = zeros([5,1]);
                Pfs.data.maxRateInd(:,dind(i),:) = zeros([5,1,numel(binDims)]);
                Pfs.data.maxRatePos(:,dind(i),:) = zeros([5,1,numel(binDims)]);
                Pfs.data.rateMap(:,dind(i),:)    = zeros([prod(Pfs.adata.binSizes),1,numIter]);
                Pfs.data.meanRate(:,dind(i))     = zeros([1,1]);
                Pfs.data.si(:,dind(i))           = zeros([1,1]);
                Pfs.data.spar(:,dind(i))         = zeros([1,1]);
                
                if returnFlag,
                    continue
                end
                
                
                res = spk(unit);
                %% Skip unit if too few spikes
                if numel(res)>21
                    
                    % SUPER annoying fix for probe shifts
                    if false%~isempty(spk.per) && any(spk.perInd(unit,:))
                        % Get the periods
                        spkPer = spk.per.copy();
                        spkPer.data(~spk.perInd(unit,:),:) = [];
                        resync  (spkPer, Session); 
                        resample(spkPer, xyz    ); % match spkPer to xyz timeframe.

                        psPer = pfsState&spkPer;
                        sstres = SelectPeriods(res,psPer.data,'d',1,1);
                        cast(spkPer, 'TimeSeries');
                        spkPer = logical(SelectPeriods(spkPer.data,pfsState,'c',1,1));
                        sstpos = asstpos(spkPer,:);
                    else
                        sstpos = asstpos;
                        switch pfsState.type
                          case 'TimePeriods'
                            sstend = round((mean(res(end-10:end)) - Session.sync.data(end)) .* xyz.sampleRate);
                            sstbeg = round((mean(res(1:10)) - Session.sync.data(1)).*xyz.sampleRate);
                            
                            spkPer = copy(pfsState);
                            spkPer.label = 'spikes';
                            spkPer.sampleRate = xyz.sampleRate;
                            spkPer.data = [sstbeg, sstend];
                            psPer = pfsState&spkPer;
                            
                            sstres = SelectPeriods(res,psPer.data,'d',1,1);
                            cast(spkPer, 'TimeSeries');

                            spkPer = logical(SelectPeriods(spkPer.data,pfsState,'c',1,1));
                            sstpos = asstpos(spkPer,:);
                            
                          case 'TimeSeries'
                            warning('Add previous cases modification');
                            sstres = res(pfsState.data(res));
                        end
                        
                    end
                    nSpk = size(sstres,1);
                    sresind = repmat(sstres,1,numIter);

                    if spkShuffle
                        % implement spkshuffle
                    end

                    if bootstrap,
                        sresind = reshape(sresind(randi(nSpk,[round(nSpk.*bootstrap),numIter]),1),[],numIter);
                    end
                    

% CHUNK data into uniform blocks of specified temporal length
                    samplesPerBlock = round(xyz.sampleRate*shuffleBlockSize);
                    blockCount      = (size(sstpos,1)-mod(size(sstpos,1),samplesPerBlock))/samplesPerBlock;
% REMOVE positions and spikes which fall outside of the final block
                    sstpos = sstpos(1:size(sstpos,1)-mod(size(sstpos,1),samplesPerBlock),:);
                    sresind(any(sresind>size(sstpos,1),2),:) = [];
                    sstres(sstres>size(sstpos,1)) = [];

                    
                    if posShuffle,                        
% $$$                         % shifts the position of the res along the sstpos
% $$$                         sresind = bsxfun(@plus,...
% $$$                                          sstres,...
% $$$                                          repmat(randi(round([-size(sstpos,1)/2,size(sstpos,1)/2]),1,numIter),...
% $$$                                                 [numel(sstres),1]));
% $$$                         % Wraps negative res to end of the sstpos vector
% $$$                         sresind(sresind<=0) = sresind(sresind<=0)+size(sstpos,1);                           
% $$$                         % Wraps res greater than the size of the sstpos vector
% $$$                         sresind(sresind>size(sstpos,1)) = sresind(sresind>size(sstpos,1))-size(sstpos,1);
                        shuffle_positions = @(p) permute(reshape(subsref(reshape(permute(p,[3,1,2]),...
                                                                                 samplesPerBlock,[],size(p,2)),...
                                                                         substruct('()',{':',randperm(size(p,1)./samplesPerBlock)',':'})),...
                                                         [],1,size(p,2)),...
                                                         [1,3,2]);
                    else,
                        shuffle_positions = @(p) p;
                    end                    

                    if halfsample,
% INDEX blocks for halfsampling
                        halfBlockCount  = (blockCount-mod(blockCount,2))/2;
                        shuffleBlockInd = reshape(1:size(sstpos,1),[],blockCount);
                        halfSampleBlockInd = zeros([numIter,halfBlockCount]);
                        halfSampleInd = zeros([halfBlockCount*samplesPerBlock,numIter]);
                        for j = 1:2:numIter,
% SPLIT blocks randomly into two groups as random half samples
                            rp = randperm(blockCount);
                            halfSampleBlockInd(j  ,:) = rp(1:halfBlockCount);
                            halfSampleBlockInd(j+1,:) = rp(halfBlockCount+1:2*halfBlockCount);
                        end
% REBUILD half samples from blocks
                        for j = 1:numIter
                            halfSampleInd(:,j) = reshape(shuffleBlockInd(:,halfSampleBlockInd(j,:)),[],1);
                        end
                    else
                        halfSampleInd = repmat([1:size(sstpos,1)]',[1,numIter]);
                    end

                    
                    
                    if isempty(posShuffleDims),
                        posShuffleDims = 1:size(sstpos,2);
                    end
                    
% COMPUTE Place Field, the first always uses all available data
                    tic; %disp(unit);
                    sstposf = sstpos;
                    sstposf(:,posShuffleDims) = shuffle_positions(sstpos(:,posShuffleDims));
                    [Pfs.data.rateMap(:,dind(i),1), ... Rate Map
                     Pfs.adata.bins,                ... Bins
                     Pfs.data.si(:,dind(i)),        ... Spatial Information
                     Pfs.data.spar(:,dind(i))]   =  ... Sparsity
                        compute_pfs(Session,                         ... MTASession Object
                                    sstposf(sresind(:,1),:),         ... Spike Postion
                                    sstposf,                         ... Marker Postion
                                    binDims,                         ... Bin Dimensions
                                    SmoothingWeights,                ... Weights
                                    type,                            ... Type {'xy','xyz' ...}
                                    boundaryLimits,                  ... Computational Boundaries
                                    xyz.sampleRate                   ... Sample Rate
                    );
                    toc

% COMPUTE Bootstrap
                    if numIter>1,
                        tic
                        for bsi = 2:numIter
                            sstposf = sstpos;
                            sstposf(:,posShuffleDims) = shuffle_positions(sstpos(:,posShuffleDims));
                            Pfs.data.rateMap(:,dind(i),bsi) = ...
                                compute_pfs(Session,...
                                            sstposf(sresind(ismember(sresind(:,bsi),halfSampleInd(:,bsi)),bsi),:),...
                                            sstposf(halfSampleInd(:,bsi),:),...
                                            binDims,...
                                            SmoothingWeights,...
                                            type,...
                                            boundaryLimits,...
                                            xyz.sampleRate...
                            );
                        end%for bsi
                        toc
                    end%if
                    
                    
                end%if too few res
                i = i+1;
            end%for unit


            if autoSaveFlag,  
% SAVE Data after Calculations                
                Pfs.save();  
% SELECT units to be returned
                field = fieldnames(Pfs.data);
                Clu = Pfs.data.clu;
                for f = 1:numel(field);
                    Pfs.data.(field{f}) = Pfs.data.(field{f})(:,ismember(Clu,units),:,:,:);
                end
                
            end
            
        end%function MTAApfs

end%methods

methods (Static)
    function status = exist(Trial,varargin)
        if isa(Trial,'MTASession'),
            assert(numel(varargin)==1,'MTAApfs:exist:MissingPfsTag');
            status = exist(fullfile(Trial.spath,[Trial.filebase,'.Pfs.' varargin{1} '.mat']),'file');
        elseif isa(Trial,'MTAApfs'),
            status = exist(fullfile(pfs.path,pfs.filename),'file');
        end
    end
end % methods (Static)

%methods
%purge_tagged_savefile
%end
    
end%classdef MTAApfs