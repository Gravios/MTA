classdef MTASession < hgsetget
% MTASession(name,varargin) 
% Data structure to organize the analysis of neural and spatial data.
%
%   name - string: Same name as the directory of the session
%
%   varargin:
%     [mazeName,overwrite,TTLValue,xyzSampleRate,xyzSystem,ephySystem]
%
%     mazeName:        string, 3-4 letter name of the testing arena 
%                      (e.g. 'rof' := rectangular open-field)
%
%     overwrite:       boolean, flag to overwrite saved Sessions
%
%     TTLValue:        string, used to synchronize position and electrophysiological data
%
%     xyzSampleRate:   double, sample rate of position tracking system
%
%     xyzSystem:       string, name/id of the system to record position
%    
%     ephySystem:      string, name/id of the system to record neural activity
%
%---------------------------------------------------------------------------------------------------------
%   General Loading:
%     
%     Load from saved Session,
%     Session = MTASession(name,mazeName);
%     
%     Create new session,
%     Session = MTASession(name,mazeName,overwrite,TTLValue,xyzSampleRate,xyzSystem,ephySystem);
%
%---------------------------------------------------------------------------------------------------------
%     examples:
%       load saved session,
%         Session = MTASession('jg05-20120309','rof');
%
%       Create New Session
%         Session = MTASession('jg05-20120309','rof',1,'0x0040',119.880135,'vicon','nlx');
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

        function Session = MTASession(name,varargin)
            [mazeName,overwrite,TTLValue,xyzSampleRate,xyzSystem,ephySystem] = ...
             DefaultArgs(varargin,{'cof',0,'0x0040',119.881035,'vicon','nlx'});
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
                    Session.create(xyzSampleRate,TTLValue,xyzSystem,ephySystem);
                else
                    warning(['Subsession with maze, ' Session.maze.name ', does not exist: creating session']);
                    Session.create(xyzSampleRate,TTLValue,xyzSystem,ephySystem);
                end  
            end
        end

        function Session = create(Session,varargin)
        % Session = create(Session,varargin)
        % Wrapper function for functions used to synchronize experimental
        % data. The choice of function is dependent upon the combination of
        % recording systems used in the session.
        %  
            [xyzSampleRate,TTLValue,xyzSystem,ephySystem] = DefaultArgs(varargin,{[],'0x8000','vicon','nlx'});
            switch ephySystem,

              case 'nlx',
                switch xyzSystem,
                  case 'vicon',
                    Session = syncViconNlx(Session,xyzSampleRate,TTLValue);
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
            [field] = DefaultArgs(varargin,{[]});
            if ~isempty(field),
                if isa(Session.(field),'MTAData')
                    Session.(field).load(Session); 
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
                if exist(Data.fpath,'file'),
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
                    Data.sync.sync.resample(Data.sampleRate);
            end
        else
            Data.sync.sync = Session.sync.copy;
            Data.sync.sync.resample(Data.sampleRate);
        end
        

        
        if isa(Data,'MTADepoch'),
            Data.data = IntersectRanges(Data.data+Data.origin,Data.sync.sync.data+Data.sync.sync.origin-1)-Data.sync.sync(1);
            Data.origin = Data.sync.sync(1)+1;
            return

        elseif isa(Data,'MTAData'),
            
            % The periods when the data was recorded
            dataEpoch = Data.sync.copy;
            dataEpoch.cast('TimeSeries',Data.sampleRate,'absolute');
            dataOrigin = Data.origin;
            
            % The periods of data which are already loaded
            loadedData = ones(Data.size(1),1);
            try
            loadedData(Data.data(:,1,1,1,1)==0|isnan(Data.data(:,1,1,1,1))|isnan(Data.data(:,1,1,1,1))) = 0;
            end
            loadedData = cat(1,zeros(dataOrigin,1),loadedData);
            tailbuff = dataEpoch.size(1)-size(loadedData,1);
            if tailbuff>0,
                loadedData = cat(1,loadedData,zeros(dataEpoch.size(1)-size(loadedData,1),1));
            else
                loadedData = loadedData(1:end+tailbuff);
            end
            
            loadedDataEnd = find(loadedData==1,1,'last');
            
            % The desired synchronization periods
            syncEpoch = Data.sync.sync.copy;
            syncEpoch.cast('TimeSeries',Data.sampleRate,'absolute');
            syncEpoch.data = syncEpoch.data(1:dataEpoch.size(1));
            newOrigin = find(syncEpoch.data==1,1,'first');
            
            %Diagnostic
            %diagnostic,
% $$$                  figure,
% $$$                  plot(dataEpoch.data)
% $$$                  ylim([-3,3])
% $$$                  hold on
% $$$                  plot(loadedData-dataEpoch.data,'r')
% $$$                  plot(loadedData-syncEpoch.data,'g')
% $$$                  plot(dataEpoch.data-syncEpoch.data,'c')
% $$$                  plot(loadedData-syncEpoch.data+dataEpoch.data,'m')
             
            
            %%Trim ends
            endSync = Data.sync.sync(end);
            endShiftIndex = endSync - loadedDataEnd+1;
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
                syncshift = Data.sync(1)-newOrigin-1;
                Data.load(syncDataPeriods,syncshift);
            end
            if ~isempty(find(syncZeroIndex,1)),
                Data.data(syncZeroIndex,:,:,:,:) = 0;
            end
            
        end
        
        Data.origin = newOrigin-1;

        end
        
        
        
        

        %% Variables from XYZ -------------------------------------------------------------%        
        function diffMat = markerDiffMatrix(Session,varargin)
        %diffMat = markerDiffMatrix(Session)
        %Create a time series where the position of every marker
        %has been substracted from one another
        %
        %  Output: 
        %
        %  diffMat - numericArray: (index,marker1,marker2,dim)
        %
            [xyz] = DefaultArgs(varargin,{Session.xyz});
            diffMat = zeros(xyz.size(1),xyz.model.N,xyz.model.N,3);
            for i=1:xyz.model.N,
                for j=1:xyz.model.N,
                    diffMat(:,i,j,:) = xyz(:,j,:)-xyz(:,i,:);
                end
            end
        end


        function angles = transformOrigin(Session, origin, orientationVector, vectorTranSet)   
            %angles = transformOrigin(Session, origin, orientationVector, vectorTranSet)   
            diffMat = Session.markerDiffMatrix();
            mdvlen = size(diffMat,1);
            origin = Session.model.gmi(origin);
            orientationVector = Session.model.gmi(orientationVector);
            
            % Get transformation Matricies
            [rz,     rzMat, direction] = rotZAxis(squeeze(diffMat(:,origin,orientationVector,:)));
            [oriVector, ryMat, pitch ] = rotYAxis(rz);
            
            tCoordinates = [];
            roll = [];
            tCMarkers = [];
            
            % Transform other marker difference vectors
            if ~isempty(vectorTranSet),
                vecTSet = zeros(mdvlen,length(vectorTranSet),3);
                for i = 1:length(vectorTranSet),
                    rztrans = sum(rzMat.* repmat(permute(shiftdim(squeeze(diffMat(:,origin,Session.model.gmi(vectorTranSet(i)),:)),-1),[2 1 3]),[1 3 1]),3);
                    rytrans = sum(ryMat.* repmat(permute(shiftdim(rztrans,-1),[2 1 3]),[1 3 1]),3);
                    vecTSet(:,i,:) = rytrans;
                    tCMarkers(end+1) = Session.model.gmi(vectorTranSet(i)); %#ok<*AGROW>
                end
                % detect head roll and remove by transformation
                [tCoordinates, roll] = detectRoll(vecTSet);
            end
            angles = struct('oriVect',oriVector,'transVec',tCoordinates,'transMarkers',tCMarkers,'direction',direction,'pitch',pitch,'roll',roll);
        end
        
        
        function v = vel(Session,varargin)
        %v = vel(Session,varargin)
        %calculate the speed of marker(s)
        %[marker,dim] = DefaultArgs(varargin,{[1:Session.model.N],[1:size(Session.xyz,3)]});
        if Session.xyz.isempty, Session.xyz.load(Session); end    
        [marker,dim] = DefaultArgs(varargin,{1:Session.model.N, 1:Session.xyz.size(3)});
            v = sqrt(sum(diff(Session.xyz(:,marker,dim),1).^2,3));
        end            

        function a = acc(Session,varargin)
        %a = acc(Session,varargin)
        %calculate the acceleration of marker(s)
        %[marker,dim] = DefaultArgs(varargin,{[1:Session.model.N],[1:size(Session.xyz,3)]});
            [marker,dim,padded] = DefaultArgs(varargin,{1:Session.model.N, 1:Session.xyz.size(3),1});
            a = diff(Session.vel(marker,dim),1);
            if padded==1
                a = cat(1,a(1,:),a,a(end,:));
            end
        end

        function center_of_mass = com(Session,Model)
        %center_of_mass = com(Session,Model)
        %
        %Model - MTAModel: MTA object holding marker information
        %
        %Examples:
        %  Find the center of mass of the session model
        %    center_of_mass = Session.com(Session.xyz.model);
        %
        %  Select a model based on a subset of markers from a larger model
        %    Model = Session.model.rb({'head_back','head_left','head_front','head_right'});
        %    center_of_mass = Session.com(Model);
        %
            center_of_mass = mean(Session.xyz(:,Session.model.gmi(Model.ml()),:),2);
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
             if Session.ang.isempty, Session.ang.load(Session);end             
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
            if ~isempty(id), id = [ '.' num2str(id)];,end
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

