classdef MTAApfs < hgsetget %< MTAAnalysis
% function Pfs = MTAApfs(Obj,varargin)            
%
%  Obj (MTATrial/MTASession/MTAApfs):
%
%    MTATrial & MTASession -> Create rate maps for neural units
%
%    MTAApfs -> update object???
%
%    otherwise -> return empty MTAApfs object
%
%-----------------------------------------------------------------------------
%
%  varargin:
%    units
%    states
%    overwrite
%    tag
%    binDims
%    SmoothingWeights
%    type
%    spkShuffle
%    posShuffle
%    numIter
%    xyzp
%    bound_lims
%    bootstrap
%
%-----------------------------------------------------------------------------
%
%  Description
%    Create an object which contains and manages rate maps,which
%    are spatially binned expected rates for individual neurons
%
%-----------------------------------------------------------------------------
%  Output:
%    Pfs (MTAApfs object) - 
%  
%
%  Examples:
    
%  MTAApfs(Obj,units,states,overwrite,tag,binDims,SmoothingWeights,type,spkShuffle,posShuffle,numIter)

    properties 
        %path - string: location of file containing MTAApfs object
        path
        
        %filename - string: location of file containing MTAApfs object
        filename = '';
        
        %session - struct(sessionName; mazeName; trialName): requried to load parent session
        session
        
        %tag - string: unique identifier
        tag
        
        %ext - string: file extention
        ext
        
        %parameters - struct:
        parameters
        
        %mdata - struct: metadata
        mdata
        
        %adata - struct: auxdata 
        adata

        %data - struct: placefield data
        data
    end

    methods

        function Pfs = MTAApfs(Obj, varargin)     
            [units,states,overwrite,tag,binDims,SmoothingWeights,type,...
             spkShuffle,posShuffle,numIter,xyzp,boundaryLimits,bootstrap]=...
            DefaultArgs(varargin,{[],'walk',0,[],[30,30],[1.2,1.2],'xy',0,0,1,MTADxyz([]),[],0});

            units = units(:)';            

            
            switch class(Obj)
                case 'MTATrial'
                    Session = Obj;
                    SessionName = Session.name;
                    MazeName    = Session.maze.name;
                    TrialName   = Session.trialName;

                    % Update map - need better method
                    Session.spk.create(Session);
                    
                    
                    if xyzp.isempty,
                        xyz = Session.xyz.copy;
                        xyz.load(Session);
                        xyz.data = sq(xyz(:,Session.trackingMarker,1:numel(binDims)));
                    else
                        xyz = xyzp;
                    end

                    
                    Pfs.path = Session.spath;
                    Pfs.tag  = tag;
                    
                    Pfs.session.sessionName = SessionName;
                    Pfs.session.trialName   = TrialName;
                    Pfs.session.mazeName    = MazeName;

                    if ischar(states),
                        Pfs.parameters.states = states;
                        pfsState = Session.stc{states,xyz.sampleRate}.copy;
                    elseif isa(states,'MTAData'),
                        pfsState = states.copy;
                        pfsState.resample(xyz);
                        Pfs.parameters.states = pfsState.label;                       
                    end
                    
                    Pfs.parameters.type   = type;
                    Pfs.parameters.spkShuffle = spkShuffle;
                    Pfs.parameters.posShuffle = posShuffle;
                    Pfs.parameters.numIter  = numIter;
                    if isempty(SmoothingWeights)&&numel(binDims)==numel(SmoothingWeights)
                        SmoothingWeights = Nbin./30;
                    end
                    Pfs.parameters.smoothingWeights   = SmoothingWeights;
                    Pfs.parameters.binDims = binDims;
                    Pfs.parameters.bootstrap = bootstrap;
                    
                    Pfs.adata.trackingMarker = Session.trackingMarker;
                    Pfs.adata.bins = [];
                    if isempty(boundaryLimits),boundaryLimits = Session.maze.boundaries(1:numel(type),:);end
                    Pfs.adata.binSizes = round(abs(diff(boundaryLimits(1:numel(binDims),:),1,2))./binDims');
                    
                    Pfs.data =struct( 'clu',        [],...
                        'elClu',      [],...
                        'el',         [],...
                        'maxRate',    [],...
                        'maxRateInd', [],...
                        'maxRatePos', [],...
                        'rateMap',    [],...
                        'meanRate',   [],...
                        'si',         [],...
                        'spar',       []);
                    Pfs.update_filename(Session,tag);
                    Pfs.ext = 'pfs';
                    
                case 'MTAApfs'
                    Pfs = Obj;
                    Session = MTASession(Obj.session.sessionName,Obj.session.mazeName);
                    Session = MTATrial(Session,Obj.session.mazeName,Obj.session.trialName);
                    pfsState = Session.stc{Pfs.states,Session.xyz.sampleRate}.copy;
                    
                otherwise
                    prop = properties('MTAAPfs');
                    for i = 1:length(prop)
                        Pfs.(prop{i})=[];
                    end
                    return
            end
                        

            numUnits = numel(units);
            selected_units = units;
            
            pf_tmpfile = Pfs.fpath;
            %% load existing data
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
                        selected_units_ind = ismember(Pfs.data.clu,units);
                        field = fieldnames(Pfs.data);
                        for f = 1:numel(field);
                            Pfs.data.(field{f}) = Pfs.data.(field{f})(:,selected_units_ind,:,:,:);
                        end
                        return
                    else
                        oldUnitInds = find(ismember(Pfs.data.clu,units));
                        numOldUnits = numel(oldUnitInds);
                        numNewUnits = numUnits - numOldUnits;
                        tnumUnits = numNewUnits + numel(Pfs.data.clu);
                        if numNewUnits>0,
                            newdata =struct( 'clu',        zeros(1,numNewUnits),...
                                'elClu',      zeros(1,numNewUnits),...
                                'el',         zeros(1,numNewUnits),...
                                'maxRate',    zeros(5,numNewUnits),...
                                'maxRateInd', zeros(5,numNewUnits,numel(binDims)),...
                                'maxRatePos', zeros(5,numNewUnits,numel(binDims)),...
                                'rateMap',    zeros(prod(Pfs.adata.binSizes),numNewUnits,numIter),...
                                'meanRate',   zeros(1,numNewUnits),...
                                'si',         zeros(1,numNewUnits),...
                                'spar',       zeros(1,numNewUnits));
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
                    newdata =struct( 'clu',        zeros(1,numNewUnits),...
                                     'elClu',      zeros(1,numNewUnits),...
                                     'el',         zeros(1,numNewUnits),...
                                     'maxRate',    zeros(5,numNewUnits),...
                                     'maxRateInd', zeros(5,numNewUnits,numel(binDims)),...
                                     'maxRatePos', zeros(5,numNewUnits,numel(binDims)),...
                                     'rateMap',    zeros(prod(Pfs.adata.binSizes),numNewUnits,numIter),...
                                     'meanRate',   zeros(1,numNewUnits),...
                                     'si',         zeros(1,numNewUnits),...
                                     'spar',       zeros(1,numNewUnits));
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
                Pfs.data =struct('clu',        zeros(1,numUnits),...
                                 'elClu',      zeros(1,numUnits),...
                                 'el',         zeros(1,numUnits),...
                                 'maxRate',    zeros(5,numUnits),...
                                 'maxRateInd', zeros(5,numUnits,numel(binDims)),...
                                 'maxRatePos', zeros(5,numUnits,numel(binDims)),...
                                 'rateMap',    zeros(prod(Pfs.adata.binSizes),numUnits,numIter),...
                                 'meanRate',   zeros(1,numUnits),...
                                 'si',         zeros(1,numUnits),...
                                 'spar',       zeros(1,numUnits));

            end
            
            
            %% Get State Positions
            if pfsState.isempty,return,end
            sstpos = sq(xyz(pfsState,:));

            
            %% load Units into spk object;
            Session.spk.create(Session,xyz.sampleRate,pfsState,units);


            i = 1;
            for unit=selected_units,

                Pfs.data.clu(dind(i)) = Session.spk.map(unit,1);
                Pfs.data.el(dind(i)) = Session.spk.map(unit,2);
                Pfs.data.elClu(dind(i)) = Session.spk.map(unit,3);
                res = Session.spk(unit);
                %% Skip unit if too few spikes
                if numel(res)>10,

                    nSpk = numel(res);
                    sstres = SelectPeriods(res,pfsState.data,'d',1,1);

                    if ~spkShuffle
                        res = repmat(res,1,numIter);
                    else
% $$$                         res = repmat(sstres,1,numIter);
% $$$                         spkswind = round(spkShuffle*xyz.sampleRate);
% $$$                         spkWindInd =  [1:spkswind:res(end)-spkswind;(spkswind+1):spkswind:res(end)];
% $$$                         for s = 2:numIter,
% $$$                             tres = [];
% $$$                             for t = spkWindInd
% $$$                             tres = vertcat(tres,resSelectPeriods(sstres,t,'d',1,0));
% $$$                         end

% $$$                         shufSpkInd = zeros([nSpk,numIter]);
% $$$                         spkswind = round(spkShuffle*xyz.sampleRate);
% $$$                         startres = res(1);
% $$$                         spkTSeries = false([res(end)-startres,1]);
% $$$                         spad = spkswind-mod(numel(spkTSeries),spkswind);
% $$$                         spkTSeries = [spkTSeries;false([spad,1])];
% $$$                         spkTSeries(res(:,1)-res(1)+1) = true;
% $$$                         spkTSeriesInd =  1:numel(spkTSeries);
% $$$                         spkTSeriesInd = reshape(spkTSeriesInd,[],spkswind);
% $$$                         nsbins = size(spkTSeriesInd,1);
% $$$                         
% $$$                         for s = 2:numIter,
% $$$                             shufInd = reshape(spkTSeriesInd(randperm(nsbins),:),[],1);
% $$$                             
% $$$                             res(:,s) = find(spkTSeries(shufInd))+startres;
% $$$                             %res(:,s) = find(reshape(spkTSeries(:,spkTSeries(randperm(nsbins),:),s),[],1))+startres;
% $$$                         end
                    end
                    
                    % shifts the position of the res along the sstpos
                    sresind = bsxfun(@plus,sstres,repmat(randi([-posShuffle,posShuffle],1,numIter),numel(sstres),1));
                    % Wraps negative res to end of the sstpos vector
                    sresind(sresind<=0) = sresind(sresind<=0)+size(sstpos,1);                           
                    % Wraps res greater than the size of the sstpos vector
                    sresind(sresind>size(sstpos,1)) = sresind(sresind>size(sstpos,1))-size(sstpos,1);   
                    
                    if bootstrap,
                        sresind = reshape(sresind(randi(nSpk,[round(nSpk.*bootstrap),numIter]),1),[],numIter);
                    end

                    %% Caluculate Place Fields
                    tic
                    [Pfs.data.rateMap(:,dind(i),1), Pfs.adata.bins] =  ...
                        PlotPF(Session,sstpos(sresind(:,1),:),sstpos,binDims,SmoothingWeights,type,boundaryLimits,xyz.sampleRate);
                    toc
                    if numIter>1,
                         for bsi = 2:numIter
                             Pfs.data.rateMap(:,dind(i),bsi) = PlotPF(Session,sstpos(sresind(:,bsi),:),sstpos,binDims,SmoothingWeights,type,boundaryLimits);
                         end
                     end
%                         Pfs.rateMap{unit} = sq(bsMap);
%                         PlaceField.stdMap{unit} = sq(std(bsMap,0,3));
%                     else
%                         PlaceField.rateMap{unit} = sq(bsMap);
%                         PlaceField.stdMap{unit} = [];
%                     end
                    
%                     try
%                         if isempty(PlaceField.rateMap{unit}), continue, end,
%                         PlaceField.rateMap{unit}(isnan(PlaceField.rateMap{unit})) = 0;
%                         PlaceField.maxRateInd{unit} = LocalMinima2(-PlaceField.rateMap{unit},-0.2,12);
%                         PlaceField.rateMap{unit}(PlaceField.rateMap{unit}(:)==0) = nan;
%                         if isempty(PlaceField.maxRateInd{unit}), continue, end,
%                         PlaceField.maxRatePos{unit} = [PlaceField.ybin(PlaceField.maxRateInd{unit}(:,2));PlaceField.xbin(PlaceField.maxRateInd{unit}(:,1))]'*[0 1; 1 0];
%                         PlaceField.maxRate{unit} = PlaceField.rateMap{unit}(round(size(PlaceField.rateMap{unit},1)*[PlaceField.maxRateInd{unit}(:,2)-1]+PlaceField.maxRateInd{unit}(:,1)));
%                         [~,PlaceField.maxRateMax{unit}] = max(PlaceField.maxRate{unit});
%                     end
                    
                end  
                i = i+1;
            end
            % Save Data after Calculations
            save(pf_tmpfile,'Pfs','-v7.3')
            field = fieldnames(Pfs.data);
            Clu = Pfs.data.clu;
            for f = 1:numel(field);
                Pfs.data.(field{f}) = Pfs.data.(field{f})(:,ismember(Clu,units),:,:,:);
            end
        end


        
        

    end
    
end