classdef MTAApfs < hgsetget %< MTAAnalysis

    properties 
        path
        filename = '';
        session
        tag
        ext
        parameters
        mdata
        adata
        data
    end

    methods

        function Pfs = MTAApfs(Obj, varargin)     
            [units,states,overwrite,tag,binDims,SmoothingWeights,type,spkShuffle,posShuffle,numIter]=...
            DefaultArgs(varargin,{[],'walk',0,[],[20,20],[],'xy','n',0,1,});
            
            switch class(Obj)
                case 'MTATrial'
                    Session = Obj;
                    SessionName = Session.name;
                    MazeName    = Session.maze.name;
                    TrialName   = Session.trialName;
                    
                    Pfs.path = Session.spath;
                    Pfs.tag  = tag;
                    
                    Pfs.session.name = SessionName;
                    Pfs.session.trialName   = TrialName;
                    Pfs.session.mazeName    = MazeName;
                    
                    pfsState = Session.stc{states,Session.xyz.sampleRate}.copy;

                    Pfs.parameters.states = states;
                    Pfs.parameters.type   = type;
                    Pfs.parameters.spkShuffle = spkShuffle;
                    Pfs.parameters.posShuffle = posShuffle;
                    Pfs.parameters.numIter  = numIter;
                    if isempty(SmoothingWeights)&&numel(binDims)==numel(SmoothingWeights)
                        SmoothingWeights = Nbin./30;
                    end
                    Pfs.parameters.smoothingWeights   = SmoothingWeights;
                    Pfs.parameters.binDims = binDims;
                    
                    Pfs.adata.trackingMarker = Session.trackingMarker;
                    Pfs.adata.bins = [];
                    Pfs.adata.binSizes = round(diff(Session.maze.boundaries(1:numel(binDims),:),1,2)./binDims');
                    
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
                    Session = MTASession(Obj.session.name,Obj.session.mazeName);
                    Session = MTATrial(Session,Obj.session.trialName);
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

            if epftmp&&~overwrite,
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
                if numUnits == 0
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
                if numUnits == 0
                    numUnits = size(Session.spk.map,1);
                    units    = 1:numUnits;
                    selected_units    = units;
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
            

            
            %% load Units into spk object;
            Session.spk.create(Session,Session.xyz.sampleRate,pfsState.label,units);
            
            %% Get State Positions
            switch type
                case 'pfcrz'
                    pfc = MTAApfs(Session,units,'theta',0,[],[30,30],[1.5,1.5],'xy');
                    tsstpos = Session.xyz(pfsState,Session.trackingMarker,:);
                otherwise
                    sstpos = Session.xyz(pfsState,Session.trackingMarker,1:numel(binDims));
            end
            

            i = 1;
            for unit=selected_units,
                Pfs.data.clu(dind(i)) = Session.spk.map(unit,1);
                Pfs.data.el(dind(i)) = Session.spk.map(unit,2);
                Pfs.data.elClu(dind(i)) = Session.spk.map(unit,3);

                spkres = Session.spk(unit);
                nSpk = numel(spkres);
                
                %% Skip unit if too few spikes
                if nSpk>10,
                    
                    
                    switch type
                        case 'pfcrz'
                            mySpkPos = sq(Session.xyz(spkres,Session.trackingMarker,:));
                            
                            [~,mind] = max(pfc.data.rateMap(:,pfc.data.clu==unit));
                            tind = cell(1,numel(pfc.adata.binSizes));
                            [tind{:}] = ind2sub(pfc.adata.binSizes',mind);
                            posind = cell2mat(tind);
                            for j = 1:numel(posind)
                                posind(j) = pfc.adata.bins{j}(posind(j));
                            end
                            mySpkPos = [sqrt(sum((mySpkPos(:,1:2) - repmat(posind,size(mySpkPos,1),1)).^2,2)), mySpkPos(:,3)];
                            sstpos = [sqrt(sum((tsstpos(:,1:2) - repmat(posind,size(tsstpos,1),1)).^2,2)), tsstpos(:,3)];
                        otherwise
                            mySpkPos = sq(Session.xyz(spkres,Session.trackingMarker,1:numel(binDims)));
                    end
                    %% CircShift position data
                    shuffled_Pos = @(posShuffle,stspos) circshift(stspos,randi([-posShuffle,posShuffle]));
                    switch spkShuffle
                        case 'n'
                            shuffled_Spk = @(nSpk) mySpkPos(1:nSpk,:);
                        case 'r'
                            shuffled_Spk = @(nSpk) mySpkPos(randi(nSpk,1,nSpk),:);
                    end
                    
                    [Pfs.data.rateMap(:,dind(i),1), Pfs.adata.bins, Pfs.data.meanRate(dind(i)), Pfs.data.si(dind(i)), Pfs.data.spar(dind(i))] =  ...
                        PlotPF(Session,shuffled_Spk(nSpk),shuffled_Pos(posShuffle,sstpos),binDims,SmoothingWeights,type);

%                     if numBSiterations>1,
%                         for bsi = 2:numBSiterations,
%                             Pfs.data.rateMap(:,dind(i),bsi) = PlotPF(Session,shuffled_Spk(nSpk),shuffled_Pos(posShuffle,sstpos),binDims,SmoothingWeights,type);
%                         end
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
                save(pf_tmpfile,'Pfs','-v7.3')
            end

            field = fieldnames(Pfs.data);
            Clu = Pfs.data.clu;
            for f = 1:numel(field);
                Pfs.data.(field{f}) = Pfs.data.(field{f})(:,ismember(Clu,units),:,:,:);
            end
        end
        
        function plot(Pfs,unit,varargin)
            [ifColorbar,colorLimits] = DefaultArgs(varargin,{0,[]});
            switch Pfs.parameters.type
                case 'xy'
                    bin1 = Pfs.adata.bins{1};
                    bin2 = Pfs.adata.bins{2};
                    rateMap = reshape(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),numel(bin1),numel(bin2));
                    imagescnan({bin1,bin2,rateMap'},colorLimits,[],ifColorbar,[0,0,0]);
                    if ~isempty(rateMap)&&~isempty(bin1)&&~isempty(bin2),
                        text(bin1(1)+30,bin2(end)-50,sprintf('%2.1f',max(rateMap(:))),'Color','w','FontWeight','bold','FontSize',10)
                    end
                    axis xy
                case 'pfcrz'
                    bin1 = Pfs.adata.bins{1};
                    bin2 = Pfs.adata.bins{2};
                    rateMap = reshape(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),numel(bin1),numel(bin2));
                    imagescnan({bin1,bin2,rateMap'},colorLimits,[],ifColorbar,[0,0,0]);
                    if ~isempty(rateMap)&&~isempty(bin1)&&~isempty(bin2),
                        text(bin1(1)+30,bin2(end)-50,sprintf('%2.1f',max(rateMap(:))),'Color','w','FontWeight','bold','FontSize',10)
                    end
                    axis xy
                case 'xyz'
                    c = eye(3);
                    r = [1.2,3,6];
                    var = cat(2,Pfs.adata.bins,{permute(reshape(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),Pfs.adata.binSizes'),[2,1,3])},{[]});
                    for i = 1:3,
                        var(end) = {max(Pfs.data.rateMap(:,Pfs.data.clu==unit,1))/r(i)};
                        fv = isosurface(var{:});
                        patch(fv,'facecolor',c(i,:),'edgecolor','none');
                        alpha(1/r(i)*r(1));
                    end
                    xlim([min(Pfs.adata.bins{1}),max(Pfs.adata.bins{1})]);
                    ylim([min(Pfs.adata.bins{2}),max(Pfs.adata.bins{2})]);
                    zlim([min(Pfs.adata.bins{3}),max(Pfs.adata.bins{3})]);
                    view(3)
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
                    '.pfs.' Pfs.parameters.type '.' Session.trackingMarker '.' pfsState.label '.' ...
                    Pfs.parameters.spkShuffle 'ps' num2str(Pfs.parameters.posShuffle) ...
                    'bs' num2str(Pfs.parameters.numIter) 'sm' smwTag...
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