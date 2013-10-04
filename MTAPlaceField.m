classdef MTAPlaceField
%PlaceField = MTAPlaceField(Session,varargin) 
%
%  Description: Data structure used to organize unit rate maps. 
%
%  varargin:
%  
%    [units,states,overwrite,trialName,type,spk_shuffle,pos_shuffle,numBSiterations,
%     smooth,nbins,mazeName,constrained_to_maze]
%                
%-------------------------------------------------------------------------------------------
%
%  Example
% 
%    Placefield = MTAPlaceField(Session,[],{'head','theta'},1)


    properties (SetAccess = public)

        %name - string: describes the placefield file 
        name

        %type - string: describes place field type (eg 'xy','rz')
        type

        %filebase - string: full file head in the format SessionName.(Maze.name).trialName
        filebase

        %mazeName - string: 3-4 letter name of the 
        %  testing arena (e.g. 'rof' := rectangular open-field)
        mazeName

        %trialName - string: designation of trial, 
        %  the full Session has the default name 'all'
        trialName

        %trackingMarker - string: Marker name used for 
        %  place field calculations
        trackingMarker

        %stateLabel - string: state(s) used to define periods 
        %  used to calculate place fields
        stateLabel

        %rateMap - cellArray -> numericArray: {clu}(xbin,ybin) 
        %  unit firing rate map
        rateMap

        %stdMap -  cellArray -> numericArray: {clu}(xbin,ybin) 
        %  standard deviation of unit firing rate map
        stdMap

        %xbin - cellArray -> numericArray: {clu}(x-position) 
        %  spatial bins associated with rateMap and stdMap
        xbin

        %ybin - cellArray -> numericArray: {clu}(y-position) 
        %  spatial bins associated with rateMap and stdMap
        ybin

        %nbins - int: number of bins with which to calculate 
        %  the place field (i.e. the place field resolution)
        nbins

        %cluMap - numericArray: (clu,electrode_clu,eletrode) mapping of unit clusters 
        %  to electrode cluster and electrode
        cluMap

        %smooth - double: smoothing factor for gaussian 
        %  smoothing of occupancy and rate maps
        smooth

        %maxRate - cellArray -> numericArray: {clu}(peak_rates_of_placefield)
        maxRate

        %maxRateMax -cellArray -> int: {clu}(index_maxRate) index of maxRate
        %  corresponding to the overall maximum rate of the place field
        maxRateMax

        %rateMap - cellArray -> numericArray: {clu}(maxRate,[xbin,ybin]) 
        maxRateInd

        %rateMap - cellArray -> numericArray: {clu}(maxRate,[x-position,y-position]) 
        maxRatePos

        %numBSiterations - int: number of random interations
        numBSiterations

        %calculation_completion_map - booleanArray: indicates if
        %  place field for a unit has been calculated
        calculation_completion_map    

        %spk_shuffle - string: shuffling method
        %                      'n': normal indexing
        %                      'r': randomly shuffled spike time index
        spk_shuffle

        %pos_shuffle - int: range of random position shuffle 
        %  (e.g. pos_shuffle = 40000 -> [-40000,40000]
        pos_shuffle
    end

    methods
        function PlaceField = MTAPlaceField(Session,varargin)
            prop = properties('MTAPlaceField');
                for i = 1:length(prop)
                    PlaceField.(prop{i})=[];
                end

                if isempty(Session),return,end
                
               [units,states,overwrite,name,trialName,type,spk_shuffle,pos_shuffle,numBSiterations,smooth,nbins,mazeName,constrained_to_maze]=DefaultArgs(varargin,{[],'walk',0,[],'all','xy','n',0,1,0.03,50,'cof',1});

                %% Load MTASession object if Session is type char
                if ischar(Session)|~isa(Session,'MTATrial')&isa(Session,'MTASession'),
                    Session = MTATrial(Session,{},trialName,[],[],mazeName);
                end

                %% Load State specific periods
                [stsp stateLabel] = Session.statePeriods(states);
                if isempty(name)
                    pf_tmpfile = [Session.spath.analysis Session.filebase ...
                                  '.Pfs.' type '.' Session.trackingMarker '.' stateLabel '.' ...
                                  spk_shuffle num2str(pos_shuffle) 'bs' num2str(numBSiterations) 'sm' num2str(smooth*100) 'bn' num2str(nbins) '.mat'];
                else
                    pf_tmpfile = [Session.spath.analysis Session.filebase '.Pfs.' name '.mat'];
                end

                %% Load PlaceField if file exists
                if exist(pf_tmpfile,'file')&overwrite~=1,
                    load(pf_tmpfile);                    
                    if isempty(units)
                        %% if no units were spcified and overwrite
                        %% was not selected, return PlaceField
                        return
                    elseif min(ismember(find(PlaceField.calculation_completion_map),units)),
                        %% if all specified units have already been
                        %% computed then return the loaded
                        %% PlaceField object
                        return
                    end
                    %% Load Units and rescale sampling freq
                    Session = Session.load_CluRes(Session.xyzSampleRate,units);
                    numClu = size(Session.map,1);
                else
                    %% Load Units and rescale sampling freq
                    Session = Session.load_CluRes(Session.xyzSampleRate,units);
                    numClu = size(Session.map,1);
                    PlaceField.calculation_completion_map = false(numClu,1);
                end
                
                PlaceField.type = type;
                PlaceField.filebase = Session.filebase;
                PlaceField.mazeName = Session.Maze.name;
                PlaceField.trialName = Session.trialName;

                %% Constrain the computational volume
                if ~constrained_to_maze,
                    maxz = max(sq(Session.xyz(:,Session.Model.gmi(Session.trackingMarker),3)));
                    minz = min(sq(Session.xyz(:,Session.Model.gmi(Session.trackingMarker),3)));
                else
                    minz = Session.Maze.boundaries(3,1);
                    maxz = Session.Maze.boundaries(3,2);
                end

                %% Instantiate PlaceField Variables
                PlaceField.trackingMarker = Session.trackingMarker;
                PlaceField.stateLabel = stateLabel;

                PlaceField.maxRate    = cell(numClu,1);
                PlaceField.maxRateMax = cell(numClu,1);
                PlaceField.maxRateInd = cell(numClu,1);
                PlaceField.maxRatePos = cell(numClu,1);
                PlaceField.xbin       = zeros(nbins,1);
                PlaceField.ybin       = zeros(nbins,1);
                PlaceField.stdMap     = cell(numClu,1);  

                if numBSiterations>1,
                    PlaceField.rateMap    = cell(numClu,1);
                else
                    PlaceField.rateMap    = repmat({zeros(nbins,nbins)},numClu,1);
                end

                PlaceField.nbins = nbins;
                PlaceField.smooth = smooth;
                PlaceField.numBSiterations = numBSiterations;
                PlaceField.spk_shuffle = spk_shuffle;
                PlaceField.pos_shuffle = pos_shuffle;

                PlaceField.cluMap = Session.map;

                %% Return dummy PlaceField for instances where a particular state(s) was never occupied
                if isempty(stsp),
                    PlaceField.calculation_completion_map = true(numClu,1);
                    return
                end




                startAtUnit = 1;
                if overwrite&isempty(units),
                    PlaceField.calculation_completion_map = false(size(PlaceField.calculation_completion_map));
                elseif ~isempty(units),
                    PlaceField.calculation_completion_map(units) = false;
                else
                    startAtUnit = find(PlaceField.calculation_completion_map,1,'first');
                    if isempty(startAtUnit),startAtUnit=1;,end
                end

                %% If no specific unit(s), calculate all pfs
                if isempty(units),
                    units = startAtUnit:numClu;
                end

                %% get State Positions and Units
                stsp = round(stsp./Session.lfpSampleRate.*Session.xyzSampleRate)+1;
                [myRes,ind] = SelectPeriods(Session.res,stsp,'d',1,1);
                myClu = Session.clu(ind);
                stspos = SelectPeriods(sq(Session.xyz(:,Session.Model.gmi(Session.trackingMarker),:)),stsp, 'c', 1);

                switch type
                  case 'xy'
                    pos = stspos;
                  case 'pfcrz'
                    pf_search = MTAPlaceField([]);
                    pf_search.type = 'xy';
                    pf_search.mazeName = 'cof';
                    pf_search.trialName = Session.trialName;
                    pf_search.trackingMarker = 'head_front';
                    pf_search.stateLabel = 'theta';
                    pf_search.spk_shuffle = 'n';
                    pf_search.pos_shuffle = 0;
                    pf_search.numBSiterations = 1;
                    pf_search.nbins = 50;
                    pf_search.smooth = 0.03;

                    Session = Session.load_Pfs(pf_search);
                    newOrigin = zeros(numClu,2);


                    for unit=units,
                        if ~isempty(Session.Pfs.maxRatePos{unit}),
                            newOrigin(unit,:) = Session.Pfs.maxRatePos{unit}(Session.Pfs.maxRateMax{unit},:);
                        end
                        PlaceField.rateMap{unit}    = zeros(nbins,nbins,numBSiterations);
                        PlaceField.stdMap{unit}     = zeros(nbins,nbins,numBSiterations);
                    end
                end

                bsMap = zeros(nbins,nbins,numBSiterations);
                %% Calculate place fields for specified units
                for unit=units,
                    %% Skip unit if already calculated
                    if PlaceField.calculation_completion_map(unit),continue,end
                    mySpk = myRes(myClu==unit);
                    nSpk = length(mySpk);

                    %% Skip unit if too few spikes
                    if nSpk>10,
                        %% CircShift position data                        

                        switch type
                          case 'pfcrz'
                            pos = cat(2,sqrt(sum((stspos(:,[1,2])-repmat(newOrigin(unit,:),size(stspos,1),1)).^2,2)),stspos(:,3));
                        end

                        shuffled_Pos = @(pos_shuffle,stspos) circshift(stspos,randi([-pos_shuffle,pos_shuffle]));
                        
                        switch spk_shuffle
                          case 'n'
                            shuffled_Spk = @(nSpk) mySpk(1:nSpk);
                          case 'r'
                            shuffled_Spk = @(nSpk) mySpk(randi(nSpk,1,nSpk));
                        end
                        
                        [bsMap(:,:,1) PlaceField.xbin PlaceField.ybin] =  PlotPF(Session,shuffled_Spk(nSpk),shuffled_Pos(pos_shuffle,pos),nbins,smooth,type);
                        if numBSiterations>1,
                            for bsi = 2:numBSiterations,
                                bsMap(:,:,bsi) = PlotPF(Session,shuffled_Spk(nSpk),shuffled_Pos(pos_shuffle,pos),nbins,smooth,type);
                            end
                            PlaceField.rateMap{unit} = sq(bsMap);
                            PlaceField.stdMap{unit} = sq(std(bsMap,0,3));
                        else
                            PlaceField.rateMap{unit} = sq(bsMap);
                            PlaceField.stdMap{unit} = [];
                        end

                        try
                            if isempty(PlaceField.rateMap{unit}), continue, end,
                            PlaceField.rateMap{unit}(isnan(PlaceField.rateMap{unit})) = 0; 
                            PlaceField.maxRateInd{unit} = LocalMinima2(-PlaceField.rateMap{unit},-0.2,12);
                            PlaceField.rateMap{unit}(PlaceField.rateMap{unit}(:)==0) = nan; 
                            if isempty(PlaceField.maxRateInd{unit}), continue, end,
                            PlaceField.maxRatePos{unit} = [PlaceField.ybin(PlaceField.maxRateInd{unit}(:,2));PlaceField.xbin(PlaceField.maxRateInd{unit}(:,1))]'*[0 1; 1 0];
                            PlaceField.maxRate{unit} = PlaceField.rateMap{unit}(round(size(PlaceField.rateMap{unit},1)*[PlaceField.maxRateInd{unit}(:,2)-1]+PlaceField.maxRateInd{unit}(:,1)));
                            [~,PlaceField.maxRateMax{unit}] = max(PlaceField.maxRate{unit});
                        end
                        
                    end

                PlaceField.calculation_completion_map(unit) = true;
                save(pf_tmpfile,'PlaceField','-v7.3')
            end
        end 

        
        function plot(PlaceField,unit,varargin)
            [ifColorbar,colorLimits] = DefaultArgs(varargin,{0,[]});
            bin1 = PlaceField.xbin;
            bin2 = PlaceField.ybin;
            rateMap = PlaceField.rateMap{unit};
            imagescnan({bin1,bin2,rateMap'},colorLimits,[],ifColorbar,[0,0,0]);
            if ~isempty(rateMap)&~isempty(bin1)&~isempty(bin2),
                text(bin1(1)+30,bin2(end)-50,sprintf('%2.1f',max(rateMap(:))),'Color','w','FontWeight','bold','FontSize',10)
            end
            axis xy
        end

        function report(PlaceField)
            filename = [PlaceField.filebase '.Pfs.' PlaceField.type '.' PlaceField.trackingMarker '.' PlaceField.stateLabel '.' ...
                        PlaceField.spk_shuffle num2str(PlaceField.pos_shuffle) 'bs' num2str(PlaceField.numBSiterations) 'sm' ...
                        num2str(PlaceField.smooth*100) 'bn' num2str(PlaceField.nbins)];
            figure(1001)
            for unit = 1:length(PlaceField.rateMap),
                if isempty(PlaceField.maxRateMax{unit}),continue,end
                if max(PlaceField.maxRate{unit})<2,continue,end
                clf
                PlaceField.plot(unit)
                reportfig(1001,filename,[],['unit: ' num2str(unit)],[],0);                
            end
        end


    end
end
