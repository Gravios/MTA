classdef MTASession < hgsetget
% MTASession(name,varargin) 
% Data structure to organize the analysis of neural and spatial data.
%
%   name - string: Same name as the directory of the session
%
%   varargin:
%     [mazeName, overwrite, TTLValue, xyzSystem, ephySystem, xyzSampleRate]
%
%     mazeName:        string, 3-4 letter name of the testing arena 
%                      (e.g. 'rof' := rectangular open-field)
%
%     overwrite:       boolean, flag to overwrite saved Sessions
%
%     TTLValue:        string, used to synchronize position and electrophysiological data
%
%     xyzSystem:       string, name/id of the system to record position
%    
%     ephySystem:      string, name/id of the system to record neural activity
%
%     xyzSampleRate:   numeric, samples per second for xyz tracking
%
%---------------------------------------------------------------------------------------------------------
%   General Loading:
%     
%     Load from saved Session,
%     Session = MTASession(name,mazeName);
%     
%     Create new session,
%     Session = MTASession(name,mazeName,overwrite,TTLValue,xyzSystem,ephySystem);
%
%---------------------------------------------------------------------------------------------------------
%     examples:
%       load saved session,
%         Session = MTASession('jg05-20120309','rof');
%
%       Create New Session
%         Session = MTASession('jg05-20120309','rof',1,'0x0040','vicon','nlx',119.881035);
%   
%---------------------------------------------------------------------------------------------------------

    properties (SetAccess = public)
        
        %path - struct: holds all paths of the constructed data tree created by MTAConfiguration.m
        path          

        %spath - struct: same as path but with Session name appended to the end
        spath

        %filebase - string: full file head in the format name.(Maze.name).trialName
        filebase       

        %name - string: name of the directory holding session information
        name           

        %trialName - string: designation of trial the full Session has the default name 'all'
        trialName      

        %Maze - MTAMaze: Object containing all maze information
        maze           

        %Model - MTAModel: Object contianing all marker information
        model          

        %sync - MTADepoch: Loading periods relative to the primary recording system
        sync

        %sampleRate - double: Sample Rate of electrophysiological recording system
        sampleRate     

        %trackingMarker - string: Marker name used for place field calculations
        trackingMarker = 'head_front';

        %xyz - MTADxyz: (Time,Marker,Dimension) XYZ position of each marker 
        xyz 

        %ang - MTADang: (Time,Marker,Marker,Dimension) marker to marker angles 
        ang 

        %stc - MTAStateCollection: Stores the periods of events
        stc
        
        %spk - MTASpk: Stores neurons' action potential timing and features
        spk

        %Fet - MTAFet: Object containing behavioral features
        fet

        %lfp - MTADlfp: (Time,channel) local field potential
        lfp

        %ufr - MTADufr: (Time,ClusterId) unit firing rates with lfpSampleRate
        ufr = [];      

        %nq - struct: clustering quality and spike waveform characteristics
        nq = {};

    end

    methods

        %% Session Constructor - Creation & Loading-----------------------------------------------------%

        function Session = MTASession(varargin)
            [name,mazeName,overwrite,TTLValue,xyzSystem,ephySystem,xyzSampleRate] = ...
             DefaultArgs(varargin,{[],'cof',0,'0x0040','vicon','nlx',[]});
            Session.path = load('MTAPaths.mat');

            if isempty(name),
                Session.filebase = '';
                Session.name = '';
            elseif isa(name,'MTASession')
                prop = properties('MTASession');
                    for i = 1:length(prop)
                        if ismethod(name.(prop{i}),'copy')
                            Session.(prop{i})=name.(prop{i}).copy;
                        else
                        Session.(prop{i})=name.(prop{i});
                        end
                    end
            else
                Session.name = name;
                Session.spath = fullfile(Session.path.data,Session.name);
                Session.trialName = 'all';
                Session.maze = MTAMaze(mazeName);
                Session.filebase = [Session.name '.' Session.maze.name '.' Session.trialName];
                if exist(fullfile(Session.spath, [Session.filebase, '.ses.mat']),'file')&&~overwrite
                    Session = Session.load();
                    Session = Session.updatePaths();
                    Session.xyz.load(Session.sync);
                elseif overwrite
                    warning(['Overwriting Session: ' fullfile(Session.spath, [Session.filebase, '.ses.mat'])])
                    Session.create(TTLValue,xyzSystem,ephySystem,xyzSampleRate);
                else
                    warning(['Subsession with maze, ' Session.maze.name ', does not exist: creating session']);
                    Session.create(TTLValue,xyzSystem,ephySystem,xyzSampleRate);
                end  
            end
        end

        function Session = create(Session,varargin)
        % Session = create(Session,varargin)
        % Wrapper function for functions used to synchronize experimental
        % data. The choice of function is dependent upon the combination of
        % recording systems used in the session.
        %  
            [TTLValue,xyzSystem,ephySystem,xyzSampleRate] = DefaultArgs(varargin,{'0x8000','vicon','nlx',[]});
            switch ephySystem,

              case 'nlx',
                switch xyzSystem,
                  case 'vicon',
                    Session = syncViconNlx(Session,TTLValue,xyzSampleRate);
                end

              case 'blackrock'
                warning(['Session creation routine does not exist ' ...
                         'for Blackrock, thank you and have a nice day'])

              otherwise
                switch xyzSystem,
                  case 'vicon'
                    Session = loadVicon(Session,xyzSampleRate);
                end
            end
        end

        function Session = updatePaths(Session,varargin)
        % Session = updatePaths(Session)
        % Change Session.path & Session.spath to the current
        % MTAPath.mat configuration found in the matlab path.
            Session.path = load('MTAPaths.mat');
            Session.spath = fullfile(Session.path.data, Session.name);
            propList = properties(Session);
            for i = 1:numel(propList),
                if isa(Session.(propList{i}),'MTAData')||isa(Session.(propList{i}),'MTAStateCollection')
                    Session.(propList{i}).updatePath(Session.spath);
                end
            end
        end

        %%---------------------------------------------------------------------------------%


        %% Save and Load-------------------------------------------------------------------%

        function save(Session)
        %save(Session)
            save(fullfile(Session.spath, [Session.filebase '.ses.mat']),'Session','-v7.3');
        end

        function Session = load(Session,varargin)
        %Session = load(Session,varargin)
        %load the session file
        %
        %  Session:     MTASession, The data object which synchronizes and
        %                           holds all experimental data.
        %
        %  field:       string,     The name of a field which belongs to
        %                           session, which will be loaded from the
        %                           data file and synchronized with the
        %                           current sync.
        %
            fvarargin = {};
            if numel(varargin)>1,
                fvarargin = varargin(2:end);
            end
            [field] = DefaultArgs(varargin,{[]});
            if ~isempty(field),
                pattern = {'MTAData';'MTASpk';'MTAStateCollection'};
                F_classes = cat(2,superclasses(Session.(field))',{class(Session.(field))});
                if any(subsref(~cellfun(@isempty,regexp(repmat(F_classes,[numel(pattern),1]),...
                                                repmat(pattern,[1,numel(F_classes)]))),...
                               substruct('()',{':'}))),
                    %if isa(Session.(field),'MTAData')||isa(Session.(field),'MTASpk'),      
                    if nargout==1,
                        Data = Session.(field).copy;
                        Data = Data.load(Session,fvarargin{:}); 
                        Session = Data;
                    else
                        Session.(field).load(Session,fvarargin{:}); 
                    end
                else
                    switch varargin{1}
                      case 'nq'
                        ds = load(fullfile(Session.spath, [Session.name '.NeuronQuality.mat']));
                        Session.nq = ds.nq;
                    end
                end
            else
                load(fullfile(Session.spath, [Session.filebase '.ses.mat']));    
            end
        end

        %%---------------------------------------------------------------------------------%

        
        
        function Data = resync(Session,Data,varargin)
        % Data = resync(Session,Data,varargin)
        % Resync uses the sync object of a session or one given through the
        % varargin option. 
        %
        %  Session: MTASession, Data object holding all session information
        %
        %  Data:    MTAData,    Data object targeted for resynchronization.
        %
        %  sync:    MTADepoch,  A set of periods defining what data needs to
        %                       be loaded. 
        %           double,     Set of periods which have to be specifed in
        %                       indecies in the sampling rate of the Data
        %                       object
        %
        [sync] = DefaultArgs(varargin,{[]});
        
        %global diagnostic;
        
        %diagnostic,
        %Data = Data.copy;
        
        switch class(Data)
          case 'MTADepoch'
            if ~strcmp(Data.path,Session.spath),
                Data.updatePath(Session.spath);
            end
            
            if exist(Data.fpath,'file')&&~isempty(Data.filename),
                Data.load;
            end
        end
        
        if ~isempty(sync)
            switch class(sync)
                case 'double',
                    msg.message = 'The provided synchronization periods are empty.';
                    msg.identifier = 'MTASession:resync:EmptySync';
                    msg.stack = dbstack;
                    if isempty(sync),error(msg),end
                    Data.sync.data = sync./Data.sampleRate+Session.sync.origin;
                    
                case 'MTADepoch'
                    Data.sync.sync = sync.copy;
            end
        else
            Data.sync.sync = Session.sync.copy;
        end
        


        
        if isa(Data,'MTADepoch'),
            %%%%%%%%%%%%%%%%%%% REDO %%%%%%%%%%%%%%%%%%%
            %Data.resample(1);
        % $$$        Data.data = IntersectRanges(Data.data+Data.origin,Data.sync.sync.data+Data.sync.sync.origin-1)-Data.sync.sync(1);
        % $$$             Data.origin = Data.sync.sync(1)+1;
            if Data.origin ~= Data.sync.sync.data(1),
                indShift = round((Data.origin - ...
                                  Data.sync.sync.data(1))*Data.sampleRate);            
            else
                indShift = 0;
            end
            syncp = Data.sync.sync.copy;
            syncp.resample(Data.sampleRate);
            
            Data.data = IntersectRanges(Data.data+indShift,syncp.data-syncp.data(1)+1);
            Data.origin = Data.sync.sync.data(1);

            
            return
            %%%%%%%%%%%%%%%%%%% REDO %%%%%%%%%%%%%%%%%%%

        elseif isa(Data,'MTAData'),
            % TODO: Need to deal with lfp resync eventually.
            if isa(Data,'MTADlfp')||Data.isempty,
                Data.sync.sync = Session.sync.copy;
                return
            end
            
            % The periods when the data was recorded
            dataEpoch = Data.sync.copy;
            dataEpoch.cast('TimeSeries',Data.sampleRate,'absolute');
            dataOrigin = round(Data.sync(1)*Data.sampleRate);
            
            % The periods of data which are already loaded
            loadedData = ones(Data.size(1),1);
            try
            loadedData(Data.data(:,1,1,1,1)==0) = 0;
            end
            loadedData = cat(1,zeros(dataOrigin,1),loadedData);
            tailbuff = dataEpoch.size(1)-size(loadedData,1);
            if tailbuff<0,
                loadedData = cat(1,loadedData,zeros(dataEpoch.size(1)-size(loadedData,1),1));
            else
                loadedData = loadedData(1:end+tailbuff);
            end
            
            loadedDataEnd = find(loadedData==1,1,'last');
            
            % The desired synchronization periods
            syncEpoch = Data.sync.sync.copy;
            syncEpoch.cast('TimeSeries',Data.sampleRate,'absolute');
            try
                syncEpoch.data = syncEpoch.data(1:dataEpoch.size(1));
            catch err
                warning('Resync may be shifted 1 index');
                syncEpoch.data = [0;syncEpoch.data];
            end

            newOrigin = find(syncEpoch.data==1,1,'first');


            %Diagnostic
            %diagnostic,
% $$$                   figure,
% $$$                   plot(dataEpoch.data)
% $$$                   ylim([-3,3])
% $$$                   hold on
% $$$                   plot(loadedData-dataEpoch.data,'r')
% $$$                   plot(loadedData-syncEpoch.data,'g')
% $$$                   plot(dataEpoch.data-syncEpoch.data,'c')
% $$$                   plot(loadedData-syncEpoch.data+dataEpoch.data,'m')
% $$$ 
% $$$                   figure,
% $$$                   plot(dataEpoch.data)
% $$$                   ylim([-3,3])
% $$$                   hold on
% $$$                   plot(loadedData+.1,'r')
% $$$                   plot(syncEpoch.data+.2,'g')
% $$$                   plot(dataEpoch.data,'c')
% $$$                   plot(loadedData-syncEpoch.data+dataEpoch.data,'m')
% $$$ 
% $$$                   ds = load(Data.fpath)
% $$$                   figure,
% $$$                   plot(ds.data(:,1,1)==eps)
% $$$                   hold on
% $$$                   ,plot(ds.data(:,1,1)==0,'r')+.5
% $$$                   ylim([-3,3])
% $$$ 
% $$$                   figure
% $$$                   plot(Data.data(:,1,1)==eps)
% $$$                   hold on
% $$$                   plot((Data.data(:,1,1)==0)+.5,'r')
% $$$                   ylim([-3,3])
% $$$                   
                  %%Trim ends
            endSync = syncEpoch.size(1);
            endShiftIndex = endSync - loadedDataEnd;
            endShiftIndex(endShiftIndex==0)=1;           
            startShiftIndex = newOrigin-dataOrigin;
            startShiftIndex(startShiftIndex==0)=1;           
            if endShiftIndex < 0,
                if startShiftIndex < 0,
                    Data.data = cat(1,zeros([abs(startShiftIndex),Data.size(2:end)]),Data.data(1:abs(endSync-dataOrigin+1),:,:,:,:));
                    dataEpoch.data = dataEpoch.data(newOrigin:endSync);
                    loadedData = loadedData(newOrigin:endSync);
                    syncEpoch.data = syncEpoch.data(newOrigin:endSync);
                else
                    Data.data = Data.data(startShiftIndex:endSync-newOrigin,:,:,:,:);
                    dataEpoch.data = dataEpoch.data(newOrigin:endSync);
                    loadedData = loadedData(newOrigin:endSync);
                    syncEpoch.data = syncEpoch.data(newOrigin:endSync);
                end
            else
                if endShiftIndex == 1, endShiftIndex = 0; end
                if startShiftIndex < 0,
                    Data.data = cat(1,zeros([abs(startShiftIndex),Data.size(2:end)]),Data.data,zeros([abs(endShiftIndex),Data.size(2:end)]));
                    dataEpoch.data = dataEpoch.data(newOrigin:endSync);
                    loadedData = loadedData(newOrigin:endSync);
                    syncEpoch.data = syncEpoch.data(newOrigin:endSync);
                else
                    Data.data = cat(1,Data.data(startShiftIndex:end,:,:,:,:),zeros([abs(endShiftIndex),Data.size(2:end)]));
                    dataEpoch.data = dataEpoch.data(newOrigin:endSync);
                    loadedData = loadedData(newOrigin:endSync);
                    syncEpoch.data = syncEpoch.data(newOrigin:endSync);
                end
            end
            
            %Diagnostic
            %diagnostic,
% $$$                  figure,
% $$$                  plot(loadedData)
% $$$                  ylim([-2,2])
% $$$                  hold on
% $$$                  plot(loadedData-dataEpoch.data,'r')
% $$$                  plot(loadedData-syncEpoch.data,'g')
% $$$                  plot(dataEpoch.data-syncEpoch.data,'c')
            %%keyboard
            
            syncFeature = (loadedData-dataEpoch.data).*syncEpoch.data-syncEpoch.data;
            syncDataIndex = syncFeature==-2;
            syncDataPeriods = ThreshCross(syncDataIndex,0.5,3);
            syncZeroIndex = syncFeature==0;
            
            if ~isempty(syncDataPeriods),
                keyboard
                %syncshift = 0;
                %syncshift = round(Data.sync(1).*Data.sampleRate)-newOrigin-1;
                %syncshift = Data.sync(1)-newOrigin-1;
                if syncshift ==-2,syncshift=0;end
                Data.load(syncDataPeriods,[],syncshift);
% $$$                 Data.load@MTAData(syncDataPeriods,[],syncshift);
% $$$                 fh = @(Data,syncDataPeriods,syncshift)load(Data,syncDataPeriods,[],syncshift);
% $$$                 feval(fh,Data,syncDataPeriods,syncshift)
                if ~isempty(find(syncZeroIndex,1)),
                    Data.data(syncZeroIndex,:,:,:,:) = 0;
                end
            else
                if ~isempty(find(syncZeroIndex,1)),
                    Data.data(syncZeroIndex,:,:,:,:) = 0;
                end
                nper = syncEpoch.copy;
                 
                if (nper.size(1)-Data.size(1))==1, 
                    nper.data = nper.data(1:Data.size(1));
                    warning(['Synchronization offset 1 index larger than data']);
                end
                
                
                Data.data = Data.data(find(nper==1,1,'first'):find(nper==1,1,'last'),:,:,:,:);
            end
            
        end
        
        % Is this right? I think it is
        if round(Data.origin*Data.sampleRate) ~= newOrigin;
            Data.origin = newOrigin/Data.sampleRate;
        end

        end
        
        
        
        

        %% Variables from XYZ -------------------------------------------------------------%    

        function angles = transformOrigin(Session,varargin)               
            %angles = transformOrigin(Session, xyz, origin, orientationVector, vectorTranSet)   
            [xyz,origin,orientationVector,vectorTranSet] = DefaultArgs(varargin,{Session.xyz.copy,'head_back','head_front',{'head_left','head_right'}});
            
            if xyz.isempty,
                xyz.load(Session);
                xyz.filter('ButFilter',3,50,'low');
            end
            
            diffMat = markerDiffMatrix(xyz);
            mdvlen = size(diffMat,1);
            origin = xyz.model.gmi(origin);
            orientationVector = xyz.model.gmi(orientationVector);
            
            % Get transformation Matricies
            [rz,     rzMat, direction] = rotZAxis(squeeze(diffMat(:,origin,orientationVector,:)));
            [oriVector, ryMat, pitch ] = rotYAxis(rz);
            
            tCoordinates = [];
            roll = [];
            tCMarkers = [];
            
            % Transform other marker difference vectors
            if ~isempty(vectorTranSet),
                vecTSet = zeros(mdvlen,numel(vectorTranSet),3);
                for i = 1:length(vectorTranSet),
                    rztrans = sum(rzMat.* repmat(permute(shiftdim(squeeze(diffMat(:,origin,xyz.model.gmi(vectorTranSet(i)),:)),-1),[2 1 3]),[1 3 1]),3);
                    rytrans = sum(ryMat.* repmat(permute(shiftdim(rztrans,-1),[2 1 3]),[1 3 1]),3);
                    vecTSet(:,i,:) = rytrans;
                    tCMarkers(end+1) = xyz.model.gmi(vectorTranSet(i)); %#ok<*AGROW>
                end
                % detect head roll and remove by transformation
                [tCoordinates, roll] = detectRoll(vecTSet);
            end
            angles = struct('oriVect',oriVector,'transVec',tCoordinates,'transMarkers',tCMarkers,'direction',direction,'pitch',pitch,'roll',roll);
        end
        
        
        %%---------------------------------------------------------------------------------%

        


        %% Load NLX Related data ---------------------------------------------------------%


        function units = selectUnits(Session,query) %#ok<INUSL>
        %units = selectUnits(Session,query)
        %Return unit ids based on a query   
        %
        %  query: cellArray(formated msg) 
        %
        %  query example:
        % 
        %    query = {{PlaceField_walk,'maxRate',@gt,3},@or,{PlaceField_rear,'maxRate',@gt,3}}
       
            unit_status = [];
            fh_list = {};
            % If the query is a single criteria wrap in a cell
            if ~iscell(query{1})
                query = {query};
            end
            % parse query
            while numel(query)~=0,
                if iscell(query{1}),
                    obj = query{1}{1};                
                    objmeta = metaclass(obj);

                    if numel(query{1})>1,
                        prop = query{1}{2};
                        operator = query{1}{3};
                        selection_value = query{1}{4};
                    end
                else
                    obj = query{1};                
                    objmeta = metaclass(obj);
                end

                switch objmeta.Name
                  case 'MTAApfs'
                    % only really meant for the property max rate
                    % at the moment
                    target_prop = obj.data.(prop);
                    value = cellfun(@max,target_prop,'UniformOutput',false);
                    value(cellfun(@isempty,value))={[0]};
                    value = cell2mat(value);                    
                    unit_status = cat(2,unit_status,operator(value,selection_value));
                  case 'struct'
                    value = obj.(prop)
                    unit_status = cat(2,unit_status,operator(value(:),selection_value));
                  case 'function_handle'
                    % stores the logical operators in a cell array
                    fh_list{end+1} = obj;
                end
                query(1) = [];
            end
            
            while size(unit_status,2)>1,
                unit_status(:,2) = fh_list{1}(unit_status(:,1),unit_status(:,2));
                unit_status(:,1) = [];
                fh_list(1) = [];
            end

            units = find(unit_status);
        end


        %%---------------------------------------------------------------------------------%



        %% Model Statistics and Corrections based on Session Data--------------------------%

%% REDO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         function Session = updateModel(Session,varargin)
         %Session = updateModel(Session,depth)
         %update or create model, used for error correction
         %
         %  depth - int: default is 3, leave it alone
         %
             [depth] = DefaultArgs(varargin,{3});
             if Session.ang.isempty, Session.ang.create(Session);end 
             dist = Session.ang(:,:,:,3);
             dist_std = zeros(Session.model.N,Session.model.N,depth);
             dist_mean = zeros(Session.model.N,Session.model.N,depth);
             dist_error = zeros(size(Session.xyz,1),Session.model.N,Session.model.N,depth);
             fdist = sum(sum(dist,2),3);
             zdist = (fdist-mean(fdist(~isnan(fdist))))/std(fdist(~isnan(fdist)));
             for i = 1:Session.model.N,
                 for j = 1:Session.model.N,
                     tdist = sq(dist(:,i,j));
                     dist_std(i,j,1)  =  std(tdist(~dist_error(:,i,j,1)&~isnan(tdist)&zdist<1));
                     dist_mean(i,j,1) = mean(tdist(~dist_error(:,i,j,1)&~isnan(tdist)&zdist<1));
                     dist_error((tdist(~isnan(tdist))>(dist_std(i,j,1)+dist_mean(i,j,1)))|(tdist(~isnan(tdist))<(-dist_std(i,j,1)+dist_mean(i,j,1))),i,j,1)=1;
                     for k = 2:depth,
                         dist_std(i,j,k)  =  std(tdist(~dist_error(:,i,j,k-1)&~isnan(tdist)&zdist<1));
                         dist_mean(i,j,k) = mean(tdist(~dist_error(:,i,j,k-1)&~isnan(tdist)&zdist<1));
                         dist_error((tdist(~isnan(tdist))>(dist_std(i,j,k)+dist_mean(i,j,k)))|(tdist(~isnan(tdist))<(-dist_std(i,j,k)+dist_mean(i,j,k))),i,j,k)=1;
                     end
                 end
             end
             Session.model.imdMean = dist_mean;
             Session.model.imdStd = dist_std;
 
             % Find Index within the Session with Best Fit - Though
             % best to select it by eye
             model_index = 0;
             for i = 0:depth-1,
                 model_index = find(sum(sum(dist_error(:,:,:,depth-i),2),3)==0);
                 model_index = model_index(randi(length(model_index)));
                 if model_index, break, end
             end    
             Session.model.index = model_index;                    
         end
 

%% REDO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
         function Session = correctRigidBody(Session,varargin)
         %Session = correctRigidBody(Session,markerSubset,depth,modelIndex,display)
         %find rigid body errors and correct them, if possible  
         %uses model frame as a basis for the geometrical
         %organization of the rigid body
         %
         %  markerSubset - cellArray: group of markers included in
         %                            the rigid body
         %                            
         %    (e.g. {'head_back','head_left','head_front','head_right','head_top'})
         %
         %  depth - int: just leave this alone, since the default works fine
         %  modelIndex - int: xyz index which contains an ideal constellation
         %                    of markers in the rigid body 
         %  display - boolean: display rb errors (only for head)
         %
 
             [markerSubset,depth,modelIndex,display] = DefaultArgs(varargin,{{'head_back','head_left','head_front','head_right','head_top'},3,[],0});
             if size(Session.model.imdMean,3)~=3|Session.model.index==0,
                 Session = Session.updateModel(depth);
             end
             rb = Session.model.rb(markerSubset);
             if ~isempty(modelIndex), Session.model.index = modelIndex; end
             pose = sq(Session.xyz(Session.model.index,Session.model.gmi(markerSubset),:));
             for i = 1:size(Session.xyzPeriods,1),
                 Temp = MTATrial(Session,{},[],Session.xyzPeriods(i,:));
                 %%REDO this part  Temp = CorrectRigidBody(Temp,rb,pose);
                 Session.xyz(Temp.xyzPeriods(1,1):Temp.xyzPeriods(1,2),:,:) = Temp.xyz;
             end
             Session.ang = Session.load_ang(1);
             if display, PlotSessionErrors(Session),end
         end
 
         function Session = correctPointErrors(Session,varargin)
         %Session = correctPointErrors(Session,varargin)
         %

             [markerSubset] = DefaultArgs(varargin,{{'head_back','head_left','head_front','head_right','head_top'}});
             rb = Session.model.rb(markerSubset);
             Session.xyz.sync.resample(Session.xyz.sampleRate);
             for i = 1:size(Session.xyz.sync,1),
                 %%REDO this part  Temp = MTATrial(Session,[],Session.xyz.sync(i,:));
                 Temp = CorrectPointErrors(Temp,rb);
                 Session.xyz(Temp.xyz.sync(1,1):Temp.xyz.sync(1,2),:,:) = Temp.xyz;
             end
             Session.ang = Session.load_ang(1);
         end

         function Session = addMarker(Session,DataObj,name,color,sticks,data)
         %Session = addMarker(Session,name,color,sticks,xyz)
         %modify Model and xyz by adding an additional marker
         %
         %  name - string: marker descriptor (e.g.'bodypart_position')
         %  color - numericArray: rgb values range - [0,1] (e.g.[0.4,1,0])
         %  sticks - cellArray: information used to visualize connections/spatial
         %                      organizations of markers
         %                      {'marker_name1','marker_name2',[stick_color]}
         %                      (e.g. {'head_front','head_back',[0,0,1]})
         %
         %             marker_name1 - string: marker descriptor (e.g.'bodypart_position')
         %             marker_name2 - string: marker descriptor (e.g.'bodypart_position')
         %             stick_color - numericArray: rgb values range - [0,1] (e.g.[0.4,1,0])
         %
             
             Marker = MTAMarker(name,color);
             DataObj.model.Markers{end+1} = Marker;
             DataObj.model.N = DataObj.model.N + 1;
             for i = 1:length(sticks),
                 Session.model.Connections{end+1} = MTAStick(sticks{i}{1},sticks{i}{2},sticks{i}{3});
             end
             DataObj.data = cat(2,DataObj.data,data);
         end

        %%---------------------------------------------------------------------------------%


        %% Place Fields -------------------------------------------------------------------%
% CHECK IF STILL FUNCTIONAL
        function Session = load_Pfs(Session,varargin)
        %Session = load_Pfs(Session,pf_search)
        %
        %loads and populates Pfs field of a Session
        %loads all saved place fields by default
        %
        %pf_search - MTAPlaceField: contains ratemap, bins and
        %                           all calculation parameters
        %
        %pf_search - cellArray: contains PlaceField states as found
        %                       in the third varargin option of the
        %                       MTAPlaceField Class
        %            e.g. 'rear', {'rear','theta'}
        %
        %-----------------------------------------------------------
        %
        %  Example:
        %    load all placefields calculated to date
        %
        %    Session = Session.load_Pfs();
        %
        %-----------------------------------------------------------
        %
        %  load placefields with specified parameters
        %
        %    pf_search = MTAPlaceField([]);
        %    pf_search.mazeName = 'cof';
        %    pf_search.trialName = 'all';
        %    pf_search.trackingMarker = 'head_front';
        %    pf_search.stateLabel = 'head.theta';
        %    pf_search.spk_shuffle = 'n';
        %    pf_search.pos_shuffle = 0;
        %    pf_search.numBSiterations = 1;
        %    pf_search.numZslices = 1;
        %    pf_search.nbins = 50;
        %    pf_search.smooth = 0.03;
        %    
        %    Session = Session.load_Pfs(pf_search);
        %
        %-----------------------------------------------------------
        %
        %  load placefields specified by state and with default
        %  parameters 
        %
        %    Session = Session.load_pfs({'rear',{'rear','theta'}}); 
        %
        %-----------------------------------------------------------

            [Pfs_search]=DefaultArgs(varargin,{{}});
            if ischar(Pfs_search),
                for i = 1:numel(Pfs_search)                    
                    Session.Pfs{1} = MTAPlaceField(Session,[],Pfs_search);
                end
                return

            elseif iscell(Pfs_search),
                Pfs = {};
                for i = 1:numel(Pfs_search)
                    Session.Pfs{i} = MTAPlaceField(Session,[],Pfs_search{i});
                end
                return
            elseif ~isempty(Pfs_search),
                Session = Session.load_MTAObject('Pfs');
                Session.Pfs = Session.getPfs(Pfs_search);
                return
            end
        end
% CHECK IF STILL FUNCTIONAL
        function Pfs = getPfs(Session,Pfs_search)            
        %Pfs = getPfs(Session,pf_search)       
        %returns existing place fields loaded in Session
        %
        %pf_search - cellArray: contains PlaceField states as found
        %                       in the third varargin option of the
        %                       MTAPlaceField Class
        %            e.g. {'rear',{'rear','theta'}};
        %                                  
        %
        %pf_search - MTAPlaceField: contains ratemap, bins and
        %                                  all calculation parameters
        %
        %See MTASession.load_Pfs() for further help
        %
            if ischar(Pfs_search),
                for i = 1:numel(Pfs_search)                    
                    Pfs = MTAPlaceField(Session,[],Pfs_search);
                end
                return

            elseif iscell(Pfs_search),
                Pfs = {};
                for i = 1:numel(Pfs_search)                    
                    Pfs{i} = MTAPlaceField(Session,[],Pfs_search{i});
                end
                return

            elseif isa(Pfs_search,'MTAPlaceField')
                Pfs = Session.findObj('Pfs',Pfs_search);
                return
            end
        end

        %%---------------------------------------------------------------------------------%

% CHECK IF STILL FUNCTIONAL
        function Session = load_MTAObject(Session,Object_File_Tag,varargin)
            [MTAObject_search,safe_load]=DefaultArgs(varargin,{{},true});
            Session.(Object_File_Tag) = {};
            files = dir(Session.spath);
            re = ['\.' Object_File_Tag '\.'];
            objFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};
            objFileSize = {files(~cellfun(@isempty,regexp({files.name},re))).bytes};
            for i = 1:length(objFileList),
                if objFileList{i},
                    if objFileSize{i}>100000000,continue,end
                    obj = load([Session.spath objFileList{i}]);
                    type = fieldnames(obj);
                    if strcmp(obj.(type{1}).mazeName,Session.maze.name),
                        Session.(Object_File_Tag){end+1} = obj.(type{1});
                    end
                end
            end
            if ~isempty(MTAObject_search),
                Session.(Object_File_Tag) = Session.findObj(MTAObject_search);
            end
        end
% CHECK IF STILL FUNCTIONAL
        function Objs = findObj(Session,Object_File_Tag,MTAObject_search)
            Objs = {};
            prop = properties(MTAObject_search);
            for j = 1:length(Session.(Object_File_Tag)),
                matchCnt = 0;
                numProp = length(prop);
                for g = 1:length(prop),

                    if isequal(Session.(Object_File_Tag){j}.(prop{g}),getfield(MTAObject_search,prop{g}))&~isempty(getfield(MTAObject_search,prop{g}))
                        matchCnt = matchCnt + 1;
                    elseif isempty(getfield(MTAObject_search,prop{g}))
                        numProp = numProp-1;
                    end
                end
                if matchCnt == numProp & matchCnt>0,
                    Objs{end+1} = Session.(Object_File_Tag){j};
                end
            end
            if length(Objs)==1,
                Objs = Objs{1};
            end
        end

        function trialNames = listTrialNames(Session)
        %trialNames = listTrialNames(Session)
        %returns a cell array of trial names associated with the Session
        
            trialNames = {};
            files = dir(Session.spath);
            re = ['\.trl\.'];
            trlFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};
            for file = 1:length(trlFileList),
                points = regexp(trlFileList{file},'[.]');
                fileparts = [1 points+1; points-1 length(trlFileList{file})]';
                if size(fileparts,1)==5,
                    trialNames{end+1} = trlFileList{file}(fileparts(3,1):fileparts(3,2));
                end
            end
        end

        function Session = filter(Session,varargin)
        %Session = filter(Session,field_name,kernel)       
        %filters a field with specified kernel
        %
        %  field_name - string/cellArray(string): string of the property to be filtered
        %
        %  kernel - numericArray: the weighted kernal used in the convolution
        %
        %  example:
        %    
        %    Session = Session.filter('xyz',gausswin(9)./sum(gausswin(9))
        %
            [field_name,kernel] = DefaultArgs(varargin,{'xyz',gausswin(9)./sum(gausswin(9))});
            padding_length = length(kernel);
            fFieldData = zeros(Session.(field_name).size);
            tFieldData = fFieldData;
            tFieldData = cat(1,tFieldData,tFieldData(1:2*padding_length,:,:,:,:));

            tFieldData(1:padding_length,:,:,:,:) = flipdim(Session.(field_name).data(1:padding_length,:,:,:,:),1);
            tFieldData(padding_length+1:size(fFieldData,1)+padding_length,:,:,:,:) = Session.(field_name).data;
            
            tFieldData(end-padding_length:end,:,:,:,:) = flipdim(Session.(field_name).data(end-padding_length:end,:,:,:,:),1);
            
            tFieldData = reshape(Filter0(kernel,tFieldData),size(tFieldData));

            Session.(field_name).data = tFieldData(padding_length+1:size(fFieldData,1)+padding_length,:,:,:,:);

        end

        function printFig(Session,varargin)
        %printFig(Session,varargin)
        %saves an image of a target figure. 
        %
        %  handle     - figureHandel: default current figure handel
        %
        %  imFileType - string: image file type, supports 'eps' and 'png'
        %
        %  id         - string: identifier for the figure, (e.g. unit_number )
        %
        %  name       - string: name, default - random number between 
        %
        %  path       - string: path where the images will be saved -
        %                      default ~/figures/f_"datestr"/
        %
        %  example:
        %
        %    Session.printFig(gcf,'png',num2str(1),[],[]);
        %
        %    Session.printFig(gcf,'png',num2str(1),'PlaceFields',Session.spath);
        %
            [handle,imFileType,id,name,path] = DefaultArgs(varargin,{gcf,'png',[],num2str(randi(100000,1)),['~/figures/f' datestr(now,29)]});
            if ~exist(path,'dir'),
                mkdir(path),
            end
            if ~isempty(id), id = [ '.' num2str(id)];end
            switch imFileType
              case 'eps'
                print(handle,'-dpsc2',[path '/' Session.filebase '.' name id '.' imFileType]);
              case 'png'
                print(handle,'-dpng',[path '/' Session.filebase '.' name id '.' imFileType]);
            end
        end

%% REDO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function Session = updateBhvSource(Session,varargin)
%             [triallist] = DefaultArgs(varargin,{{Session.trialName}});
%             tempSession = Session;
%             for i = 1:length(triallist),
%                 if ~strcmp(Session.trialName,triallist{i})
%                     Session = MTATrial(Session.name,{},triallist{i},[],0,Session.maze.name);
%                 end
%                 bhvModeNames = {};
%                 files = dir(Session.spath);
%                 re = [Session.name '\.' Session.maze.name '\.' Session.trialName '\.bhv\.'];
%                 bhvFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};
%                 bhvDnumList = {files(~cellfun(@isempty,regexp({files.name},re))).datenum};
%                 re = ['\.mat']; 
%                 bhvFileListBadInd = find(cellfun(@isempty,regexp(bhvFileList,re)));
%                 bhvFileList(bhvFileListBadInd) = [];
%                 bhvDnumList(bhvFileListBadInd) = [];
%                 [~,newestBhvFileInd] = max(cell2mat(bhvDnumList));
%                 if isempty(newestBhvFileInd),continue,end
%                 load([Session.spath bhvFileList{newestBhvFileInd}])
%                 bhvmode = Bhv.mode;
%                 Session.Bhv = {};
%                 clear('Bhv');
%                 Session = Session.load_Bhv(bhvmode);
%                 Session.save;
%             end
%             Session = tempSession;;
%         end
%
%% REDO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function Trial = consolidateTrials(Session,varargin)
%             if ~isa(Session,'MTASession')&isa(Session,'MTATrial'),
%                 Session = MTASession(Session.name,{},Session.maze.name);
%             end
%             [newTrialName,trialnames] = DefaultArgs(varargin,{'cnsldtd',Session.list_trialNames});
%             re = 'all';
%             allInd = find(~cellfun(@isempty,regexp(trialnames,re)));
%             if ~isempty(allInd),
%                 trialnames(allInd) = [];                        
%             end
%             re = newTrialName;
%             allInd = find(~cellfun(@isempty,regexp(trialnames,re)));
%             if ~isempty(allInd),
%                 trialnames(allInd) = [];                        
%             end
%             Bhv = MTABhv([],'cnsldtd');
%             newXYZPeriods = [];
%             for i = 1:length(trialnames),
%                 Trial = MTATrial(Session,{},trialnames{i},[],0,Session.maze.name,'minimal');
%                 newXYZPeriods = cat(1,newXYZPeriods,Trial.xyzPeriods);
%                 if i==1,
%                     Bhv.States = Trial.Bhv.States;
%                     for j = 1:length(Bhv.States),
%                         Bhv.States{j}.state = Bhv.States{j}.state+Trial.xyzPeriods(1);
%                     end
%                 else
%                     for j = 1:length(Bhv.States),
%                         Bhv.States{j}.state = cat(1,Bhv.States{j}.state,Trial.Bhv.getState(Bhv.States{j}.label).state+Trial.xyzPeriods(1));
%                         if i==length(trialnames),
%                             Bhv.States{j}.state = sort(Bhv.States{j}.state)-min(newXYZPeriods(:));
%                         end
%                     end
%                 end                                
%             end
%             newXYZPeriods = sort(newXYZPeriods);
%             Trial = MTATrial(Session,{},newTrialName,newXYZPeriods,1,Session.maze.name,'minimal');            
%             Trial.Bhv = Bhv;
%             Bhv.save(Trial,1);
%             Trial.save
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
        
    end %methods

end %MTASession

