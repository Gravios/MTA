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
             spkShuffle,posShuffle,numIter,xyzp,bound_lims,bootstrap]=...
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
                    if isempty(bound_lims),bound_lims = Session.maze.boundaries(1:numel(type),:);end
                    Pfs.adata.binSizes = round(abs(diff(bound_lims(1:numel(binDims),:),1,2))./binDims');
                    
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
                    Pfs.updateFilename(Session,tag);
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

            if epftmp&&overwrite~=1,
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

            elseif epftmp&&overwrite,
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
                        PlotPF(Session,sstpos(sresind(:,1),:),sstpos,binDims,SmoothingWeights,type,bound_lims,xyz.sampleRate);
                    toc
                    if numIter>1,
                         for bsi = 2:numIter
                             Pfs.data.rateMap(:,dind(i),bsi) = PlotPF(Session,sstpos(sresind(:,bsi),:),sstpos,binDims,SmoothingWeights,type,bound_lims);
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
        
        function rateMap = plot(Pfs,varargin)
            [unit,nMode,ifColorbar,maxRate,isCircular] = DefaultArgs(varargin,{[],'mean',0,[],true});

            if isempty(unit),unit=Pfs.data.clu(1);end
            switch numel(Pfs.parameters.type)
                case 2
                   bin1 = Pfs.adata.bins{1};
                   bin2 = Pfs.adata.bins{2};
                   %MTATrial.validate(Pfs.session).maze.shape
% $$$                    switch 'bob'
% $$$                      case 'circle'
% $$$                        width = Pfs.adata.binSizes(1);
% $$$                        height = Pfs.adata.binSizes(2);
% $$$                        radius = round(Pfs.adata.binSizes(1)/2)-...
% $$$                                      find(Pfs.adata.bins{1}<-420,1,'last');
% $$$                        centerW = width/2;
% $$$                        centerH = height/2;
% $$$                        [W,H] = meshgrid(1:width,1:height);           
% $$$                        mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
% $$$                        mask(mask==0)=nan;
% $$$                      otherwise
% $$$                        mask = 1;
% $$$                    end
                   %% This is only valid for the circular maze
                   if isCircular,
                       width = Pfs.adata.binSizes(1);
                       height = Pfs.adata.binSizes(2);
                       radius = round(Pfs.adata.binSizes(1)/2)-find(Pfs.adata.bins{1}<-420,1,'last');
                       centerW = width/2;
                       centerH = height/2;
                       [W,H] = meshgrid(1:width,1:height);           
                       mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
                       mask(mask==0)=nan;
                   else
                       mask = 1;
                   end

                  

                    switch nMode
                        case 'mean'
                            rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),3);
                        case 'std'
                            rateMap = std(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),[],3);
                        case 'sig'
                            rateMap = 1./sum((repmat(max(Pfs.data.rateMap(:,Pfs.data.clu==unit,:)),[size(Pfs.data.rateMap,1),1,1])...
                                                     -repmat(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),[1,1,Pfs.parameters.numIter]))<0,3)';
                        otherwise
                            rateMap = Pfs.data.rateMap(:,Pfs.data.clu==unit,1);
                    end
                    
                    
                    rateMap = reshape(rateMap,numel(bin1),numel(bin2)).*mask;
                    
                    if nargout>0,return,end
                    
                    if isempty(maxRate),
                        maxRate = max(rateMap(:));
                    end

                    rateMap(isnan(rateMap)) = -1;
                    imagesc(bin1,bin2,rateMap');

                    text(Pfs.adata.bins{1}(end)-250,Pfs.adata.bins{2}(end)-50,...
                        sprintf('%2.1f',max(rateMap(:))),'Color','w','FontWeight','bold','FontSize',8)
                    colormap([0,0,0;parula]);
                    caxis([-1,maxRate]);        
                    
        % $$$                     imagescnan({bin1,bin2,rateMap'},colorLimits,[],ifColorbar,[0,0,0]);
        % $$$                     
        % $$$                     if ~isempty(rateMap)&&~isempty(bin1)&&~isempty(bin2),
        % $$$                         text(bin1(1)+30,bin2(end)-50,sprintf('%2.1f',max(rateMap(:))),'Color','w','FontWeight','bold','FontSize',10)
        % $$$                     end
                    ;                    axis xy
              case 3
                c = eye(3);
                r = [1.2,3,6];
                rateMap = permute(reshape(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),Pfs.adata.binSizes'),[1,2,3]);
                if nargout>0,return,end
                %var = cat(2,Pfs.adata.bins,{permute(reshape(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),Pfs.adata.binSizes'),[2,1,3])},{[]});
                hrate = max(Pfs.data.rateMap(:,Pfs.data.clu==unit,1))/2;
                [mind] = LocalMinimaN(-rateMap,-hrate,9);

                    
                    m = find(ismember('xyz',nMode));
                    o = find(~ismember([1,2,3],m));
                    ss = {};
                    if ~isempty(mind),
                        mind = mind(1,:);
                        ss{m(1)} = ':';
                        ss{m(2)} = ':';
                        ss{o}    = mind(o);
                        imagescnan({Pfs.adata.bins{m(1)},Pfs.adata.bins{m(2)},sq(nanmean(subsref(rateMap,substruct('()',ss)),o))'},colorLimits,[],ifColorbar,[0,0,0]);
% $$$                         if ~isempty(rateMap)&&~isempty(Pfs.adata.bins{m(1)})&&~isempty(Pfs.adata.bins{m(2)}),
% $$$                             text(Pfs.adata.bins{m(1)}+30,Pfs.adata.bins{m(1)}-50,sprintf('%2.1f',max(rateMap(:))),'Color','w','FontWeight','bold','FontSize',10)
% $$$                         end
                        axis xy
                    end

        % $$$                     [X,Y,Z] = meshgrid(Pfs.adata.bins{:});
        % $$$                     var = permute(reshape(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),Pfs.adata.binSizes'),[2,1,3]);
        % $$$                     if nargout==0,                        
        % $$$                         XS = Pfs.adata.bins{1}(round(Pfs.adata.binSizes(1)/4):round(Pfs.adata.binSizes(1)/4):Pfs.adata.binSizes(1));
        % $$$                         YS = Pfs.adata.bins{2}(round(Pfs.adata.binSizes(2)/4):round(Pfs.adata.binSizes(2)/4):Pfs.adata.binSizes(2));
        % $$$                         slice(X,Y,Z,var,XS,YS,Pfs.adata.bins{3}(round(Pfs.adata.binSizes(3)/2)))
        % $$$                         xlim([min(Pfs.adata.bins{1}),max(Pfs.adata.bins{1})]);
        % $$$                         ylim([min(Pfs.adata.bins{2}),max(Pfs.adata.bins{2})]);
        % $$$                         zlim([min(Pfs.adata.bins{3}),max(Pfs.adata.bins{3})]);
        % $$$                         %view(3)
        % $$$                         if ~isempty(colorLimits), caxis(colorLimits); end
        % $$$                         if ifColorbar, colorbar; end
        % $$$                     else
        % $$$                         rateMap = var;
        % $$$                     end


% $$$                     for i = 1:3,
% $$$                         var(end) = {max(Pfs.data.rateMap(:,Pfs.data.clu==unit,1))/r(i)};
% $$$                         fv = isosurface(var{:});
% $$$                         patch(fv,'facecolor',c(i,:),'edgecolor','none');
% $$$                         alpha(1/r(i)*r(1));
% $$$                     end
% $$$                     xlim([min(Pfs.adata.bins{1}),max(Pfs.adata.bins{1})]);
% $$$                     ylim([min(Pfs.adata.bins{2}),max(Pfs.adata.bins{2})]);
% $$$                     zlim([min(Pfs.adata.bins{3}),max(Pfs.adata.bins{3})]);
% $$$                     view(3)
            end
        end

        function [mxr,mxp] = maxRate(Pfs,varargin)
            [units,isCircular] = DefaultArgs(varargin,{[],true},false);

            if isempty(units),
                units = Pfs.data.clu;
            end
            if isCircular,
                width = Pfs.adata.binSizes(1);
                height = Pfs.adata.binSizes(2);
                radius = round(Pfs.adata.binSizes(1)/2)-find(Pfs.adata.bins{1}<-420,1,'last');
                centerW = width/2;
                centerH = height/2;
                [W,H] = meshgrid(1:width,1:height);           
                mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
                mask(mask==0)=nan;
                if numel(Pfs.parameters.type)>2,
                    mask = repmat(mask,[1,1,Pfs.adata.binSizes(3)]);
                end
                mask = reshape(mask,[],1);
            else
                mask = 1;
            end

            mxr = nan(numel(units),1);
            mxp = nan(numel(units),1);
            
            for u = units(:)',
                rateMap = Pfs.data.rateMap(:,Pfs.data.clu==u,1);
                if size(rateMap,2)==0,
                    rateMap = nan([size(Pfs.data.rateMap,1),1]);
                end

                [mxr(u==units),mxp(u==units)] = max(rateMap.*mask);
            end
            mxp = Ind2Sub(Pfs.adata.binSizes',mxp);
            if numel(Pfs.parameters.type)>2,
                mxp = [Pfs.adata.bins{1}(mxp(:,1)), ...
                       Pfs.adata.bins{2}(mxp(:,2)), ...
                       Pfs.adata.bins{3}(mxp(:,3))];                                
            else
                mxp = [Pfs.adata.bins{1}(mxp(:,1)), ...
                       Pfs.adata.bins{2}(mxp(:,2))];
            end
            
        end

        function rho = spatialCoherence(Pfs,units)
        % rho = spatialCoherence(Pfs,units);
        % The correlation between a list of firing rates in each
        % pixel and a corresponding list of firing rates averaged
        % over the 8 nearestneighbors of each pixel. 
        % muller_kubie_1989
            if isempty(units),
                units = Pfs.data.clu;
            end
            rho = nan(numel(units),1);
            for u = units',
                rm = Pfs.plot(u);
                cmat = [1,1,1;1,0,1;1,1,1];
                c = conv2(rm,cmat,'valid')./8;
                rmss = rm(2:end-1,2:end-1);
                ind = ~isnan(c)&~isnan(rmss);
                rho(u==units) = corr(c(ind),rmss(ind));
            end
        end
        
        
        function Pfs = updateFilename(Pfs,Session,varargin)
            Pfs.tag= varargin{1};

            pfsState = Session.stc{Pfs.parameters.states,Session.xyz.sampleRate};
            if isempty(Pfs.tag)
                binDimTag = num2str(Pfs.parameters.binDims);
                binDimTag(isspace(binDimTag)&isspace(circshift(binDimTag',1)'))=[];
                binDimTag(isspace(binDimTag)) = '_';
                smwTag = num2str(Pfs.parameters.smoothingWeights*100);
                smwTag(isspace(smwTag)&isspace(circshift(smwTag',1)'))=[];
                smwTag(isspace(smwTag)) = '_';
                Pfs.filename = [Session.filebase ...
                    '.pfs.' Pfs.parameters.type '.' Session.trackingMarker '.' Session.stc.mode '.' pfsState.label '.' ...
                    'ss' num2str(Pfs.parameters.spkShuffle) 'ps' num2str(Pfs.parameters.posShuffle) ...
                    'ni' num2str(Pfs.parameters.numIter) 'bs' num2str(Pfs.parameters.bootstrap) 'sm' smwTag...
                    'bd' binDimTag '.mat'];
            else
                Pfs.filename = [Session.filebase '.Pfs.' Pfs.tag '.mat'];
            end
        end
        


        function path = fpath(Pfs)
            path = fullfile(Pfs.path,Pfs.filename);
        end

    end
    
end