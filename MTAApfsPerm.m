classdef MTAApfsPerm < hgsetget %< MTAAnalysis

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

        function Pfs = MTAApfsPerm(Obj, varargin)     
        % MTAApfs(Obj,{units,states,overwrite,tag,binDims,SmoothingWeights,type,spkShuffle,posShuffle,numIter})
            [units,states,overwrite,tag,binDims,SmoothingWeights,type,spkShuffle,posShuffle,numIter,xyz,bound_lims]=...
            DefaultArgs(varargin,{[],{'rear','walk'},0,[],[30,30],[1.2,1.2],'xy',0,0,1,MTADxyz([]),[]});

            units = units(:)';            

                
            switch class(Obj)
                case 'MTATrial'
                    Session = Obj;
                    SessionName = Session.name;
                    MazeName    = Session.maze.name;
                    TrialName   = Session.trialName;
                    
                    if xyz.isempty,
                        xyz = Session.xyz.copy;
                        xyz.load(Session);
                        xyz.data = sq(xyz(:,Session.trackingMarker,1:numel(binDims)));
                    end

                    
                    Pfs.path = Session.spath;
                    Pfs.tag  = tag;
                    
                    Pfs.session.name = SessionName;
                    Pfs.session.trialName   = TrialName;
                    Pfs.session.mazeName    = MazeName;


                    if prod(cellfun(@ischar,states)),
                        Pfs.parameters.states = states;
                        pfsState{1} = Session.stc{states{1},xyz.sampleRate}.copy;
                        pfsState{2} = Session.stc{states{2},xyz.sampleRate}.copy;
                    elseif isa(states,'MTAData'),
                        pfsState{1} = states{1}.copy;
                        pfsState{2} = states{2}.copy;                        
                        pfsState{1}.resample(xyz);
                        pfsState{2}.resample(xyz);                        
                        Pfs.parameters.states = {pfsState{1}.label,pfsState{2}.label};
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
            %sstpos = sq(xyz(pfsState,:));

            
            %% load Units into spk object;
            Session.spk.create(Session,xyz.sampleRate,'theta',units);




            
            pfsState{1}.cast('TimeSeries');
            pfsState{2}.cast('TimeSeries');

            durStateA = sum(pfsState{1}.data);
            durStateB = sum(pfsState{2}.data);
            
            pooledState = pfsState{1}.data+pfsState{2}.data;
            pooledState = find(pooledState);
            
            assert(durStateA+durStateB==length(pooledState),['States ', ...
                               'labeling must be non degenerate'])
            
            i = 1;
            for unit=selected_units,
                disp(['Calculating pfsPerm for unit: ' num2str(unit)]);

                Pfs.data.clu(dind(i)) = Session.spk.map(unit,1);
                Pfs.data.el(dind(i)) = Session.spk.map(unit,2);
                Pfs.data.elClu(dind(i)) = Session.spk.map(unit,3);
                res = Session.spk(unit);

                if  numel(res)< 10,continue,end
                sstposA = xyz(logical(pfsState{1}.data),:);
                sstposB = xyz(logical(pfsState{2}.data),:);                    

                nSpk = numel(res);

                spkindA = ismember(res,find(pfsState{1}.data==1));
                if sum(spkindA)~=0, 
                    spkposA = xyz(res(spkindA),:);
                else
                    spkposA = [];
                end
                
                spkindB = ismember(res,find(pfsState{2}.data==1));
                if sum(spkindB)~=0,
                    spkposB = xyz(res(spkindB),:);
                else
                    spkposB = [];
                end                
                
                tic
                %% Caluculate Place Fields
                [rateMapA, Pfs.adata.bins] = PlotPF(Session,spkposA,sstposA,binDims,SmoothingWeights,type,bound_lims,xyz.sampleRate);
                [rateMapB, Pfs.adata.bins] = PlotPF(Session,spkposB,sstposB,binDims,SmoothingWeights,type,bound_lims,xyz.sampleRate);
                Pfs.data.rateMap(:,dind(i),1) = rateMapA-rateMapB;

                if numIter>1,
                    for bsi = 2:numIter
                        rind = randperm(length(pooledState));
                        permStateA = pooledState(rind(1:durStateA));
                        permStateB = pooledState(rind(durStateA+1:end));

                        sstposA = xyz(permStateA,:);
                        sstposB = xyz(permStateB,:);                    

                        spkindA = ismember(res,permStateA);
                        if sum(spkindA)~=0,
                            spkposA = xyz(res(spkindA),:);
                        else
                            spkposA = [];
                        end
                        
                        spkindB = ismember(res,permStateB);
                        if sum(spkindB)~=0,
                            spkposB = xyz(res(spkindB),:);
                        else
                            spkposB = [];
                        end
                        

                        [rateMapA, Pfs.adata.bins] = PlotPF(Session,spkposA,sstposA,binDims,SmoothingWeights,type,bound_lims,xyz.sampleRate);
                        [rateMapB, Pfs.adata.bins] = PlotPF(Session,spkposB,sstposB,binDims,SmoothingWeights,type,bound_lims,xyz.sampleRate);
                        Pfs.data.rateMap(:,dind(i),bsi) = rateMapA-rateMapB;
                    end
                end
                toc

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
            [unit,nMode,ifColorbar,colorLimits] = DefaultArgs(varargin,{[],'mean',0,[]});

            if isempty(unit),unit=Pfs.data.clu(1);end
            switch numel(Pfs.parameters.type)
              case 2
                bin1 = Pfs.adata.bins{1};
                bin2 = Pfs.adata.bins{2};
                
                switch nMode
                  case 'mean'
                    rateMap = mean(Pfs.data.rateMap(:,Pfs.data.clu==unit,2:end),3);
                  case 'std'
                    rateMap = std (Pfs.data.rateMap(:,Pfs.data.clu==unit,2:end),[],3);
                  case 'sig'
        % $$$                     rateMapT = 1./sum(repmat(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),[1,1,Pfs.parameters.numIter])>Pfs.data.rateMap(:,Pfs.data.clu==unit,:),3);
        % $$$                     rateMap = nan([Pfs.adata.binSizes']);
        % $$$                     try,rateMap(rateMapT<0.05&nniz(rateMapT)) = rateMapT(rateMapT<0.05&nniz(rateMapT));end
        % $$$                     try,rateMap((1-rateMapT)<0.05&nniz(rateMapT)) = -(1-rateMapT((1-rateMapT)<0.05&nniz(rateMapT)));end
                    rateMap = zscore(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),[],3);
                    rateMap = rateMap(:,:,1);
                    %rateMap = 1-normcdf(abs(rateMap(:,1,1)),0,1);
                  case 'sig'
        % $$$                     rateMapT = 1./sum(repmat(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),[1,1,Pfs.parameters.numIter])>Pfs.data.rateMap(:,Pfs.data.clu==unit,:),3);
        % $$$                     rateMap = nan([Pfs.adata.binSizes']);
        % $$$                     try,rateMap(rateMapT<0.05&nniz(rateMapT)) = rateMapT(rateMapT<0.05&nniz(rateMapT));end
        % $$$                     try,rateMap((1-rateMapT)<0.05&nniz(rateMapT)) = -(1-rateMapT((1-rateMapT)<0.05&nniz(rateMapT)));end
                    rateMap = zscore(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),[],3);
                    rateMap = rateMap(:,:,1);
                  case 'sigks'
        % $$$                     rateMapT = 1./sum(repmat(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),[1,1,Pfs.parameters.numIter])>Pfs.data.rateMap(:,Pfs.data.clu==unit,:),3);
        % $$$                     rateMap = nan([Pfs.adata.binSizes']);
        % $$$                     try,rateMap(rateMapT<0.05&nniz(rateMapT)) = rateMapT(rateMapT<0.05&nniz(rateMapT));end
        % $$$                     try,rateMap((1-rateMapT)<0.05&nniz(rateMapT)) = -(1-rateMapT((1-rateMapT)<0.05&nniz(rateMapT)));end
                    rateMap = zscore(Pfs.data.rateMap(:,Pfs.data.clu==unit,:),[],3);
                    rateMap = rateMap(:,:,1);
                    %rateMap = 1-normcdf(abs(rateMap(:,1,1)),0,1);
                    for u = 1:numel(unit),
                        mask_kstest = zeros([size(rateMap,1),numel(unit)]);
                        for i = 1:size(rateMap,1),
                            try,mask_kstest(i,u) = kstest(rateMap(i,u,:),'alpha',.00001);end
                        end
                    end
                    rateMap(~repmat(mask_kstest,[1,1,size(rateMap,3)]))=nan;
                    %rateMap = 1-normcdf(abs(rateMap(:,1,1)),0,1);
                  otherwise
                    rateMap = Pfs.data.rateMap(:,Pfs.data.clu==unit,1);
                end
                
                
                rateMap = reshape(rateMap,numel(bin1),numel(bin2));
                if nargout>0,return,end
                
                imagescnan({bin1,bin2,rateMap'},colorLimits,[],ifColorbar,[0,0,0]);
                
                if ~isempty(rateMap)&&~isempty(bin1)&&~isempty(bin2),
                    text(bin1(1)+30,bin2(end)-50,sprintf('%2.1f',max(rateMap(:))),'Color','w','FontWeight','bold','FontSize',10)
                end
                axis xy
              case 3
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

            pfsState{1} = Session.stc{Pfs.parameters.states{1},Session.xyz.sampleRate};
            pfsState{2} = Session.stc{Pfs.parameters.states{2},Session.xyz.sampleRate};
            pfsStateLabel = [pfsState{1}.key,pfsState{2}.key];
            if isempty(Pfs.tag)
                binDimTag = num2str(Pfs.parameters.binDims);
                binDimTag(isspace(binDimTag)&isspace(circshift(binDimTag',1)'))=[];
                binDimTag(isspace(binDimTag)) = '_';
                smwTag = num2str(Pfs.parameters.smoothingWeights*100);
                smwTag(isspace(smwTag)&isspace(circshift(smwTag',1)'))=[];
                smwTag(isspace(smwTag)) = '_';
                Pfs.filename = [Session.filebase ...
                    '.pfsPerm.' Pfs.parameters.type '.' Session.trackingMarker '.' Session.stc.mode '.' pfsStateLabel '.' ...
                    'ss' num2str(Pfs.parameters.spkShuffle) 'ps' num2str(Pfs.parameters.posShuffle) ...
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