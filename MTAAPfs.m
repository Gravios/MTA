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

        function Pfs = MTAApfs(Session_args, parameters, varargin)     
            %Pfs = MTAAnalysis();
            prop = properties('MTAAPfs');
            for i = 1:length(prop)
                Pfs.(prop{i})=[];
            end
            
            %if isempty(Session_args),return,end
            
            %if isa(Session_args,'MTATrial'),
                Session = Session_args;
                SessionName = Session.name;
                MazeName    = Session.maze.name;
                TrialName   = Session.trialName;
%             elseif iscell(Session_args)
%                 [SessionName,MazeName,TrialName] = DefaultArgs(Session_args,{'jg05-20120317','cof','all'});
%                 Session = MTASession(SessionName,[],MazeName);
%                 Session = MTATrial(Session,{{'Spk',S}},TrialName);
%             end
            
            Pfs.ext = 'pfs';
            
            [units,states,overwrite,tag,type,spkShuffle,posShuffle,numIter,smooth,binDims]=DefaultArgs(varargin,{[],'walk',0,[],'xy','n',0,1,0.03,[20,20]});
            % Use parameters to populate database            
            
            Pfs.path = Session.spath;
            Pfs.tag  = tag;
            
            Pfs.session.name = SessionName;
            Pfs.session.trialName   = TrialName;
            Pfs.session.mazeName    = MazeName;

            pfsState = Session.stc{states,Session.xyz.sampleRate}.copy;

            Pfs.parameters.states = pfsState.label;
            Pfs.parameters.type   = type;
            Pfs.parameters.spkShuffle = spkShuffle;
            Pfs.parameters.posShuffle = posShuffle;
            Pfs.parameters.numIter  = numIter;
            Pfs.parameters.smooth   = smooth;
            Pfs.parameters.binDims = binDims;

            Pfs.adata.trackingMarker = Session.trackingMarker;
            Pfs.adata.bins = [];
            Pfs.adata.binSizes = diff(Session.maze.boundaries(1:numel(binDims),:),1,2)./binDims';
                        
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
            
            numUnits = numel(units);
            
            
            pf_tmpfile = Pfs.fpath;
            %% load existing data
            epftmp = exist(pf_tmpfile,'file');
            if epftmp&&overwrite~=1,
                load(pf_tmpfile);
                % Load specific units
                if isempty(units)
                    return
                else
                    selected_units = ismember(Pfs.data.clu,units);
                    field = fieldnames(Pfs.data);
                    for f = 1:numel(field);
                        Pfs.data.(field{f}) = Pfs.data.(field{f})(:,selected_units,:,:,:);
                    end
                    return
                end
                
            elseif ~epftmp
            %% Instantiate Pfs Data Variables if DeNovoCalc
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
                
            elseif epftmp&&overwrite==1,
            %% Extend Pfs data for additional units
                load(pf_tmpfile);
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
                end
                field = fieldnames(newdata);
                for f = 1:numel(field);
                    Pfs.data.(field{f}) = cat(2,Pfs.data.(field{f}),newdata.(field{f}));
                end
                dind = [oldUnitInds(:)',[tnumUnits-numNewUnits+1:tnumUnits]];
            end
            
            %% load Units into spk object;
            Session.spk.create(Session,Session.xyz.sampleRate,pfsState.label,units);
            
            %% Get State Positions
            sstpos = Session.xyz(pfsState,Session.trackingMarker,:);

            i = 1;
            for unit=units,
                Pfs.data.clu(dind(i)) = Session.spk.map(unit,1);
                Pfs.data.el(dind(i)) = Session.spk.map(unit,2);
                Pfs.data.elClu(dind(i)) = Session.spk.map(unit,3);
                mySpkPos = Session.xyz(Session.spk(unit),Session.trackingMarker,:);
                nSpk = size(mySpkPos,1);
                
                %% Skip unit if too few spikes
                if nSpk>10,
                    %% CircShift position data
                    shuffled_Pos = @(posShuffle,stspos) circshift(stspos,randi([-posShuffle,posShuffle]));
                    switch spkShuffle
                        case 'n'
                            shuffled_Spk = @(nSpk) mySpkPos(1:nSpk,:);
                        case 'r'
                            shuffled_Spk = @(nSpk) mySpkPos(randi(nSpk,1,nSpk),:);
                    end
                    
                    [Pfs.data.rateMap(:,dind(i),1), Pfs.adata.bins, Pfs.data.meanRate, Pfs.data.si, Pfs.data.spar] =  PlotPF(Session,shuffled_Spk(nSpk),shuffled_Pos(posShuffle,sstpos),binDims,smooth,type);
%                     if numBSiterations>1,
%                         for bsi = 2:numBSiterations,
%                             bsMap(:,:,bsi) = PlotPF(Session,shuffled_Spk(nSpk),shuffled_Pos(pos_shuffle,pos),binDims,smooth,type);
%                         end
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
                save(pf_tmpfile,'Pfs','-v7.3')
            end
            
            
        end
        
        function plot(Pfs,unit,varargin)
            [ifColorbar,colorLimits] = DefaultArgs(varargin,{0,[]});
            switch Pfs.parameters.type
                case 'xy'
                    bin1 = Pfs.adata.bins(:,1);
                    bin2 = Pfs.adata.bins(:,2);
                    rateMap = reshape(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),numel(bin1),numel(bin2));
                    imagescnan({bin1,bin2,rateMap'},colorLimits,[],ifColorbar,[0,0,0]);
                    if ~isempty(rateMap)&&~isempty(bin1)&&~isempty(bin2),
                        text(bin1(1)+30,bin2(end)-50,sprintf('%2.1f',max(rateMap(:))),'Color','w','FontWeight','bold','FontSize',10)
                    end
                    axis xy
            end
        end

        function Pfs = updateFilename(Pfs,Session,varargin)
            Pfs.tag= varargin{1};

            pfsState = Session.stc{Pfs.parameters.states,Session.xyz.sampleRate};
            if isempty(Pfs.tag)
                binDimTag = num2str(Pfs.parameters.binDims);
                binDimTag(isspace(binDimTag)&isspace(circshift(binDimTag',1)'))=[];
                binDimTag(isspace(binDimTag)) = '_';
                Pfs.filename = [Session.filebase ...
                    '.pfs.' Pfs.parameters.type '.' Session.trackingMarker '.' pfsState.label '.' ...
                    Pfs.parameters.spkShuffle 'ps' num2str(Pfs.parameters.posShuffle) ...
                    'bs' num2str(Pfs.parameters.numIter) 'sm' num2str(Pfs.parameters.smooth*100)...
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