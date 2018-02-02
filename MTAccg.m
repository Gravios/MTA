classdef  MTAccg

    properties 
        name
        method
        filebase
        mazeName
        trialName

        % analysis
        binSize
        halfBins
        normalization
        numIterations
        surrogate_sample_size
        partitions
        Description

        %sampleRate
        tbin 

        %data
        ccg

        %auxData
        CluTags
        posind
        cluMap
        partition_feature
        partition_boundaries
        calculation_completion_map   
    end

    methods
        
        function Bccg = MTAccg(Trial,name,Description,ResTrain,CluTags,varargin)

            Trial = MTATrial.validate(Trial);
            
            prop = properties('MTAccg');
            %comment
            for i = 1:length(prop)
                Bccg.(prop{i})=[];
            end

            if isempty(Trial),return,end

% DEFARGS ------------------------------------------------------------------------------------------
            defargs = struct('units',                    [],                                     ...
                             'overwrite',                false,                                  ...
                             'rand_tag',                 1,                                      ...
                             'method',                   'normal',                               ...
                             'partitions',               1,                                      ...
                             'partition_feature',        [],                                     ...
                             'surrogate_sample_size',    [],                                     ...
                             'numIterations',            1,                                      ...
                             'binSize',                  100,                                    ... 
                             'halfBins',                 20,                                     ...
                             'normalization',            'hz'                                    ...
            );
             [units,overwrite,rand_tag, method, partitions,  partition_feature, ...
             surrogate_sample_size, numIterations, binSize, halfBins, normalization] = ...
                 DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

            %% Load Trial data
            if partitions == 1, method = 'normal';end
            Trial.spk.create(Trial,Trial.lfp.sampleRate);
            numClu = size(Trial.spk.map,1);

            Bccg.mazeName = Trial.maze.name;
            Bccg.trialName = Trial.trialName; 
            Bccg.name = name;
            Bccg.filebase = Trial.filebase;
            Bccg.method = method;
            Bccg.cluMap = Trial.spk.map;     
            Bccg.partitions = partitions;
            Bccg.CluTags = CluTags;
            Bccg.Description = Description;
            Bccg.binSize = binSize;
            Bccg.halfBins = halfBins;
            Bccg.normalization = normalization;
            Bccg.surrogate_sample_size =surrogate_sample_size;
            if isempty(ResTrain),
                Bccg.ccg = [];
                Bccg.partition_feature = [];
                Bccg.partition_boundaries = [];
                Bccg.posind = [];       
                Bccg.numIterations = [];
                return
            end        



            %% Check if Trial already contains pf with given parameters
            %% If terminal error occurred during a previous run pick up from the last unit.    
            pf_tmpfile = fullfile(Trial.spath ,[Trial.filebase '.ccg.' name '.mat']);

            if ~exist(pf_tmpfile,'file'),
                Bccg.calculation_completion_map = false(numClu,1);
            end

            startAtUnit = 1;
            if overwrite&isempty(units),
                Bccg.calculation_completion_map = false(size(Bccg.calculation_completion_map));
            elseif overwrite&~isempty(units),
                Bccg.calculation_completion_map(units) = false;
            elseif ~isempty(units),
                if length(units)==sum(Bccg.calculation_completion_map(units)),
                    return
                else
                    startAtUnit = find(Bccg.calculation_completion_map,1,'first');
                    if isempty(startAtUnit),
                        startAtUnit=1;
                    end
                end
            end


            %% If no specific unit(s), calculate all pfs
            if isempty(units),
                units = startAtUnit:numClu;
            end


            numClu = size(Trial.spk.map,1);
            numResTrains = length(ResTrain);
            Bccg.ccg = zeros(halfBins*2+1,numClu,length(ResTrain),partitions,numIterations);
            posind = {};
            partBoundary = {};
            partition_distMap = {};

            uRes = Trial.spk.res;
            xyz = preproc_xyz(Trial,'trb');


            %% Change res sampling rate to lfp.sampleRate and create
            %% surrogates where ResTrain is a list of periods.
            surrogateIndex = {};
            for nr = 1:numResTrains,               
                surrogateIndex{nr} = [];
                if size(ResTrain{nr},2)==2,
                    for i = 1:size(ResTrain{nr},1)
                        surrogateIndex{nr} = cat(1,surrogateIndex{nr},[ResTrain{nr}(i,1):ResTrain{nr}(i,2)]');
                    end
                end

                if ~isempty(surrogateIndex{nr}),
                    if ~isempty(surrogate_sample_size),
                        rpi = randi(size(surrogateIndex{nr},1),surrogate_sample_size,numIterations);
                    else
                        rpi = randi(size(surrogateIndex{nr},1),length(ResTrain{1}),numIterations);
                    end
                    

                    for i = 1:numIterations,
                        posind{nr,i} = round((surrogateIndex{nr}(rpi(:,i))-1)./Trial.lfp.sampleRate*xyz.sampleRate)+1;
                        tRes{nr,i} = surrogateIndex{nr}(rpi(:,i));
                        posind{nr,i}(posind{nr,i}>size(xyz,1)) = [];
                        tRes{nr,i}(posind{nr,i}>size(xyz,1)) = [];
                    end
                else
                    for i = 1:numIterations,
                        tRes{nr,i} = ResTrain{nr};
                        posind{nr,i} = round((ResTrain{nr}-1)./Trial.lfp.sampleRate.*xyz.sampleRate)+1;
                        posind{nr,i}(posind{nr,i}>size(xyz,1)) = [];
                        tRes{nr,i}(posind{nr,i}>size(xyz,1)) = [];
                        
                    end
                end%if ~isempty(surrogateIndex{nr}),
                
                


                %% Partition Res into groups based on a distance metric
                if partitions>1,
                    switch method,                        
                      case 'prctile_dist'
                        for unit = units,
                            for i = 1:numIterations,
                                try
                                    partition_distMap{nr,i}(unit,:) = ...
                                        sqrt(sum((repmat(partition_feature(unit,:),size(posind{nr,i},1),1)...
                                                  -sq(xyz(posind{nr,i},Trial.trackingMarker,[1,2]))).^2,2));
                                catch err,
                                    disp(err);
                                end
                            end
                        end
                        
                      case 'abs_dist'
                        for unit = units,
                            for iter = 1:numIterations,                     
                                try
                                    partition_distMap{nr,iter}(unit,:) = ...
                                        sqrt(sum((repmat(partition_feature(unit,:),size(posind{nr,iter},1),1)...
                                                  -sq(xyz(posind{nr,iter},Trial.trackingMarker,[1,2]))).^2, ...
                                                 2));
                                catch err,
                                    disp(err);
                                end
                            end
                        end%for unit
                    end%switch method
                end%if partitions>1
            end%for nr

            
            
            %% CCG of non-rearing periods and units
            for i=1:numIterations,
                Clu = num2cell([1:length(ResTrain)]);

                for unit = units,
                    for nr = 1:numResTrains,
                        if partitions>1,
                            switch method,
                              case 'prctile_dist'                    
                                partBoundary{nr,i}(unit,:) = prctile(partition_distMap{nr,i}(unit,:),[0:100/partitions:100]);
                              case 'abs_dist'
                                partBoundary{nr,i}(unit,:) = [0:max(abs(diff(Trial.maze.boundaries,1,2)))/partitions:max(abs(diff(Trial.maze.boundaries,1,2)))];
                            end
                        end
                    end
                    for np = 1:partitions,
                        if partitions>1,
                            for nr = 1:numResTrains,
                                Res{nr} = tRes{nr,i}((partition_distMap{nr,i}(unit,:)>=partBoundary{nr,i}(unit,np))'&(partition_distMap{nr,i}(unit,:)<=partBoundary{nr,i}(unit,np+1))');
                            end
                        else 
                            Res = tRes;
                        end

                        Res{numResTrains+1} = uRes(Trial.spk.clu==unit);
                        Clu{numResTrains+1} = Trial.spk.clu(Trial.spk.clu==unit);

                        [rccg Bccg.tbin ] = Trains2CCG(Res(:),Clu(:),binSize,halfBins,Trial.lfp.sampleRate,normalization);
                        emptyResInd = cellfun(@isempty,Res);

                        for nr = 1:numResTrains,
                            if ~emptyResInd(nr)&~emptyResInd(end),
                                Bccg.ccg(:,unit,nr,np,i) = sq(rccg(:,nr-sum(emptyResInd(1:nr)),end));
                            else 
                                Bccg.ccg(:,unit,nr,np,i) = false(halfBins*2+1,1);
                            end
                        end


                    end
                end
            end


            Bccg.partition_feature = partition_distMap;
            Bccg.partition_boundaries = partBoundary;
            Bccg.posind = posind;       
            Bccg.numIterations = numIterations;


            if rand_tag,
                save(fullfile(Trial.spath ,[Trial.filebase '.ccg.' Bccg.name '.' num2str(randi([100000,999999],1)) '.mat']),'Bccg','-v7.3');
            else
                save(fullfile(Trial.spath ,[Trial.filebase '.ccg.' Bccg.name '.mat']),'Bccg','-v7.3');
            end

            Bccg.save(Trial,rand_tag);
        end

        function save(Bccg,Trial,rand_tag)
            if rand_tag,
                %% Seed RandStream with system clock
                RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
                save(fullfile(Trial.spath, [Trial.filebase '.ccg.' Bccg.name '.' num2str(randi([100000,999999],1)) '.mat']),'Bccg','-v7.3');
            else
                save(fullfile(Trial.spath, [Trial.filebase '.ccg.' Bccg.name '.mat']),'Bccg','-v7.3');
            end
        end

        function fccg = filter(Bccg,kernel)
            padding_length = length(kernel);
            fccg = zeros(size(Bccg.ccg));
            dimpad = zeros(1,ndims(Bccg.ccg));
            dimpad(1) = padding_length*2;
            tccg = zeros(size(Bccg.ccg)+dimpad);

            tccg(1:padding_length,:,:,:,:) = flipdim(Bccg.ccg(1:padding_length,:,:,:,:),1);
            tccg(padding_length+1:size(fccg,1)+padding_length,:,:,:,:) = Bccg.ccg;
            tccg(end-padding_length:end,:,:,:,:) = flipdim(Bccg.ccg(end-padding_length:end,:,:,:,:),1);
            tccg = reshape(Filter0(kernel./sum(kernel),tccg),size(tccg));
            fccg = tccg(padding_length+1:size(fccg,1)+padding_length,:,:,:,:);
        end


        function axis_handel = plot(Bccg,unit,varargin)
            [partitions,filterKernel,ifColorBar] = DefaultArgs(varargin,{1,gausswin(5)./sum(gausswin(5)),1});
            if length(filterKernel) > 1,
                fccg = Bccg.filter(filterKernel);
            elseif filterKernel~=0,
                fccg = Bccg.filter(gausswin(filterKernel));
            else
                fccg = Bccg.ccg;
            end

            if length(partitions) > 1 
                axis_handel = imagescnan({Bccg.tbin,Bccg.partition_boundaries});
            else
                bar(Bccg.tbin/1000,Bccg.ccg(:,unit,partitions));
            end
        end
    end

end


