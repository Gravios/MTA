classdef MTAAknnpfs < hgsetget %< MTAAnalysis
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

        function Pfs = MTAAknnpfs(Obj, varargin)     
        %            [units,states,overwrite,tag,binDims,nNearestNeighbors,distThreshold,...
        %             type,ufrShufBlockSize,numIter, downSampleRate]=...

            % DEFARGS ----------------------------------------------------------------------
            defargs = struct('units',              [],                                   ...
                             'states',             {{'walk'}},                           ...
                             'overwrite',          false,                                ...
                             'tag',                [],                                   ...
                             'ufr',                [],                                   ...
                             'binDims',            [20,20],                              ...
                             'nNearestNeighbors' , 60,                                   ...
                             'distThreshold',      125,                                  ...
                             'type',               'xy',                                 ...
                             'ufrShufBlockSize',   1,                                    ... 
                             'numIter',            1,                                    ...
                             'pos',                [],                                   ...
                             'sampleRate',         10,                                   ...
                             'absTimeSubSample',   300,                                  ...
                             'nNodes',             4                                     ...
            );            
            [units,states,overwrite,tag,ufr,binDims,nNearestNeighbors,distThreshold,...
                type,ufrShufBlockSize,numIter,pos,sampleRate,absTimeSubSample,nNodes]=...
            DefaultArgs(varargin,defargs,'--struct');
            %-----------------------------------------------------------------------------             
            
            units = units(:)';
        
            switch class(Obj)
                case 'MTATrial'
                    %% Base values for every analysis -> perhaps create superclass later
                    Session = Obj;

                    Pfs.path = Session.spath;
                    Pfs.tag  = tag;
                    
                    % Session struct holding the basic information about
                    % the session used to create the analysis object.
                    Pfs.session.name = Session.name;
                    Pfs.session.trialName   = Session.trialName;
                    Pfs.session.mazeName    = Session.maze.name;
                    
                    pfsState = Session.stc{states};
                    
                    Pfs.parameters.states = states;
                    Pfs.parameters.type   = type;
                    Pfs.parameters.ufrShufBlockSize= ufrShufBlockSize;
                    Pfs.parameters.numIter  = numIter;
                    Pfs.parameters.nNearestNeighbors   = nNearestNeighbors;
                    Pfs.parameters.distThreshold = distThreshold;
                    Pfs.parameters.binDims = binDims;
                    Pfs.parameters.sampleRate = sampleRate;
                    
                    Pfs.adata.trackingMarker = Session.trackingMarker;
                    Pfs.adata.bins = [];
                    Pfs.adata.binSizes = round(diff(Session.maze.boundaries(1:numel(binDims),:),1,2)./binDims');
                    
                    Pfs.data =struct( 'clu',        [],...
                                      'elClu',      [],...
                                      'el',         [],...
                                      'rateMap',    []);
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
                    Pfs.rmClu(0);
                end
                % Load specific units
                if isempty(units)
                    return
                else
                    unprocessed_units = ~ismember(units,Pfs.data.clu);
                    if sum(unprocessed_units)==0
                        selected_units_ind = ismember(Pfs.data.clu,units);
                        field = fieldnames(Pfs.data);
                        for f = 1:numel(field)
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
                                'rateMap',    zeros(prod(Pfs.adata.binSizes),numNewUnits,numIter));
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
                    Pfs.rmClu(0);
                end
                oldUnitInds = find(ismember(Pfs.data.clu,units));
                numOldUnits = numel(oldUnitInds);
                numNewUnits = numUnits - numOldUnits;
                tnumUnits = numNewUnits + numel(Pfs.data.clu);
                if numNewUnits>0,
                    newdata =struct( 'clu',        zeros(1,numNewUnits),...
                                     'elClu',      zeros(1,numNewUnits),...
                                     'el',         zeros(1,numNewUnits),...
                                     'rateMap',    zeros(prod(Pfs.adata.binSizes),numNewUnits,numIter));
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
                                 'rateMap',    zeros(prod(Pfs.adata.binSizes),numUnits,numIter));
            end
            

            %% Get State Positions
            if isempty(pos), 
                pos = Session.load('xyz');
                pos.resample(sampleRate);                
                sstpos = sq(pos(pfsState,Session.trackingMarker,1:numel(binDims)));                 
            elseif size(pos,2)==1,
                pos.resample(sampleRate);
                sstpos = sq(pos(pfsState,1,1:numel(binDims)));
            end


            %% load unit firing rate


            if isempty(ufr),
                ufr = Session.ufr.copy;
                ufr = ufr.create(Session,pos,pfsState.label,selected_units,0.8);
            end
            sstufr = ufr(pfsState,:);

            
           
            %% Trim ufr and xyz
            if numIter>1
                samplesPerBlock = round(Session.xyz.sampleRate*ufrShufBlockSize);
                vlen = size(sstpos,1);
                trim = mod(vlen,samplesPerBlock);
                sstpos = sstpos(1:vlen-trim,:);
                sstufr = sstufr(1:vlen-trim,:);
                ufrShufBlockCount = (vlen-trim)/samplesPerBlock;
                ufrBlockInd = reshape(1:size(sstufr,1),[],ufrShufBlockCount);

                ufrShufPermIndices = zeros(numIter,ufrShufBlockCount);
                ufrShufPermIndices(1,:) = 1:ufrShufBlockCount;
                for i = 2:numIter,
                    ufrShufPermIndices(i,:) = randperm(ufrShufBlockCount);
                end
            else 
                ufrBlockInd = ':';
                ufrShufPermIndices = ':';
            end
            

            
            parp = [];
            if numIter>1,
                try,delete(gcp('nocreate')),end
                parp = parpool(nNodes);
            end
            

            i = 1;
            for unit=selected_units(:)',
                Pfs.data.clu(dind(i)) = Session.spk.map(unit,1);
                Pfs.data.el(dind(i)) = Session.spk.map(unit,2);
                Pfs.data.elClu(dind(i)) = Session.spk.map(unit,3);
                
                %% Calculate Place Fields
                tic
                [Pfs.data.rateMap(:,dind(i),1), Pfs.adata.bins,dw,diw] =  ...
                PlotKNNPF(Session,...
                          sstufr(reshape(ufrBlockInd(:,ufrShufPermIndices(1,:)),[],1),unit==selected_units),...
                          sstpos,...
                          binDims,...
                          nNearestNeighbors,...
                          distThreshold,...
                          type);
                toc

                tic
                if numIter>1,
                    RandRateMap = zeros([size(Pfs.data.rateMap,1),numIter-1]);

                    
                    parfor bsi = 1:numIter-1,
                        RandRateMap(:,bsi) = ...
                        PlotKNNPF...
                        (Session,sstufr(reshape(ufrBlockInd(:,ufrShufPermIndices(bsi+1,:)),[],1),unit==selected_units),...
                         sstpos,binDims,nNearestNeighbors,distThreshold,type);%,dw,diw);
                    end
                    Pfs.data.rateMap(:,dind(i),2:end) = RandRateMap;
                end
                toc
                
                i = i+1;
                %save(pf_tmpfile,'Pfs','-v7.3')
            end

            save(pf_tmpfile,'Pfs','-v7.3')

            % PAR POOL close 
            if ~isempty(parp), 
                delete(parp);
            end

            
            field = fieldnames(Pfs.data);
            Clu = Pfs.data.clu;

            for f = 1:numel(field);
                Pfs.data.(field{f}) = Pfs.data.(field{f})(:,ismember(Clu,units),:,:,:);
            end
        end
        
        function rateMap = plot(Pfs,varargin)
        % rateMap = plot(Pfs,unit,varargin)
        % [nMode,ifColorbar,colorLimits,sigDistr] = DefaultArgs(varargin,{'',0,[],[]});
            [unit,mode,ifColorbar,colorLimits,sigDistr,isCircular] = DefaultArgs(varargin,{[],'slice',false,[],[],true});
            if isempty(unit),unit=Pfs.data.clu(1);end
            switch Pfs.parameters.type
                case 'xy'

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
                  
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
                  bin1 = Pfs.adata.bins{1};
                  bin2 = Pfs.adata.bins{2};
                  
                  switch mode
                    case 'field'
                      rateMap = Pfs.data.rateMap(:,Pfs.data.clu==unit,1);
                    case 'mean'
                      rateMap = nanmean(Pfs.data.rateMap(:,Pfs.data.clu==unit,2:end),3);
                    case 'std'
                      rateMap = nanstd(Pfs.data.rateMap(:,Pfs.data.clu==unit,2:end),[],3);
                    case 'sig'
                      rateMap = 1./nansum((repmat(max(Pfs.data.rateMap(:,Pfs.data.clu==unit,:)),...
                                                  [size(Pfs.data.rateMap,1),1,1])...
                                          -repmat(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),...
                                                  [1,1,Pfs.parameters.numIter]))<0,3)';
                    otherwise
                      if isnumeric(mode)
                          rateMap = Pfs.data.rateMap(:,Pfs.data.clu==unit,mode);
                      else
                          rateMap = Pfs.data.rateMap(:,Pfs.data.clu==unit,1);
                      end
                      
                  end
                  
                  
                  %% Normal rateMap = reshape(rateMap',numel(bin1),numel(bin2))';
                  rateMap = reshape(rateMap,numel(bin2),numel(bin1))'.*mask;
                  
                  if nargout==0,                    

                      if isempty(colorLimits),
                          colorLimits = [-1,max(rateMap(:))];
                      end

                      rateMap(isnan(rateMap)) = -1;
                      imagesc(bin1,bin2,rateMap');
                      colormap([0,0,0;parula]);
                      caxis(colorLimits);        

                      axis xy
                      text(Pfs.adata.bins{1}(end)-250,Pfs.adata.bins{2}(end)-50,...
                          sprintf('%2.1f',max(rateMap(:))),'Color','w','FontWeight','bold','FontSize',10)
                      
                  end
                  return
                  
              case 'xyz'
                switch mode
                  case 'isosurface'
                    c = eye(3);
                    r = [1.2,3,6];
                    var = cat(2,Pfs.adata.bins, ...
                              {permute(reshape(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),Pfs.adata.binSizes'),[2,1,3])},{[]});
                    if nargout==0,                    
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
                    else
                        rateMap = var;
                    end
                  case 'slice'
                    
                    [X,Y,Z] = meshgrid(Pfs.adata.bins{:});
                    var = permute(reshape(Pfs.data.rateMap(:,Pfs.data.clu==unit,1),Pfs.adata.binSizes'),[2,1,3]);
                    if nargout==0,                        
                        XS = Pfs.adata.bins{1}(round(Pfs.adata.binSizes(1)/4):round(Pfs.adata.binSizes(1)/4):Pfs.adata.binSizes(1));
                        YS = Pfs.adata.bins{2}(round(Pfs.adata.binSizes(2)/4):round(Pfs.adata.binSizes(2)/4):Pfs.adata.binSizes(2));
                        slice(X,Y,Z,var,XS,YS,Pfs.adata.bins{3}(round(Pfs.adata.binSizes(3)/2)))
                        xlim([min(Pfs.adata.bins{1}),max(Pfs.adata.bins{1})]);
                        ylim([min(Pfs.adata.bins{2}),max(Pfs.adata.bins{2})]);
                        zlim([min(Pfs.adata.bins{3}),max(Pfs.adata.bins{3})]);
                        %view(3)
                        if ~isempty(colorLimits), caxis(colorLimits); end
                        if ifColorbar, colorbar; end
                    else
                        rateMap = var;
                    end
                    

                end
            end
        end

        function Pfs = updateFilename(Pfs,Session,varargin)
            Pfs.tag= varargin{1};

            pfsState = Session.stc{Pfs.parameters.states};
            if isempty(Pfs.tag)
                binDimTag = num2str(Pfs.parameters.binDims);
                binDimTag(isspace(binDimTag)&isspace(circshift(binDimTag',1)'))=[];
                binDimTag(isspace(binDimTag)) = '_';
                nnnTag = num2str(Pfs.parameters.nNearestNeighbors);
                Pfs.filename = [Session.filebase ...
                    '.pfk.' Pfs.parameters.type '.' Session.trackingMarker '.' pfsState.label '.' ...
                    'sr' num2str(Pfs.parameters.sampleRate) ...
                    'us' num2str(Pfs.parameters.ufrShufBlockSize) ...
                    'bs' num2str(Pfs.parameters.numIter) ...
                    'nnn' nnnTag ...
                    'dt' num2str(Pfs.parameters.distThreshold) ...
                    'bd' binDimTag ...
                    '.mat'];
            else
                Pfs.filename = [Session.filebase '.pfk.' Pfs.tag '.mat'];
            end
        end
        


        function path = fpath(Pfs)
            path = fullfile(Pfs.path,Pfs.filename);
        end

        function Pfs = rmClu(Pfs,clu)
            cind = Pfs.data.clu==clu;
            switch class(Pfs.data)
              case 'struct'
                fnames = fieldnames(Pfs.data);
                for f = 1:numel(fnames),
                    Pfs.data.(fnames{f})(:,cind,:) = [];
                end
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
            for u = units,
                rm = Pfs.plot(u);
                cmat = [1,1,1;1,0,1;1,1,1];
                c = conv2(rm,cmat,'valid')./8;
                rmss = rm(2:end-1,2:end-1);
                ind = ~isnan(c)&~isnan(rmss);
                rho(u==units) = corr(c(ind),rmss(ind));
            end
        end
        
        function [mxr,mxp] = maxRate(Pfs,units,varargin)
            [mode,isCircular] = DefaultArgs(varargin,{'first',true},1);

            % create mask to remove regions outside accessible area
            if isCircular,
                width = Pfs.adata.binSizes(1);
                height = Pfs.adata.binSizes(2);
                radius = round(Pfs.adata.binSizes(1)/2)-find(Pfs.adata.bins{1}<-420,1,'last');
                centerW = width/2;
                centerH = height/2;
                [W,H] = meshgrid(1:width,1:height);           
                mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
                mask(mask==0)=nan;
                if numel(Pfs.adata.binSizes)>2,
                    mask = repmat(mask,[1,1,Pfs.adata.binSizes(3)]);
                end
                mask = repmat(mask(:),[1,1,Pfs.parameters.numIter]);                
            else
                mask = 1;
            end
            
            
            if isempty(units),
                units = Pfs.data.clu;
            end
            mxr = nan(numel(units),1);
            mxp = nan(numel(units),1);
            switch mode
              case 'first';
                ind = 1;
                mfunc = @(x,y) x;
              case 'mean';
                ind = ':';
                mfunc = @nanmean;
            end
            
            for u = units,
                [mxr(u==units),mxp(u==units)] = max(mfunc(Pfs.data.rateMap(:,Pfs.data.clu==u,ind).*mask(:,1,ind),3));
            end
            mxp = Ind2Sub(Pfs.adata.binSizes',mxp);
            mxp = [Pfs.adata.bins{2}(mxp(:,2)),Pfs.adata.bins{1}(mxp(:,1))];
        end
    end
    
end