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
                    %Session.xyz.load(Session.sync);
                elseif overwrite
                    warning(['Overwriting Session: ' fullfile(Session.spath, [Session.filebase, '.ses.mat'])])
                    Session.create(TTLValue,xyzSystem,ephySystem,xyzSampleRate);
                else
                    warning(['Subsession with maze, ' Session.maze.name ', does not exist: creating session']);
                    Session.create(TTLValue,xyzSystem,ephySystem,xyzSampleRate);
                end  
            end
        end



        
        
        
        
        
        
        
        
        
        
        
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

    methods (Static)
        function Session = validate(Session)
            
            if isa(Session,'MTASession'),
                return;
            elseif ischar(Session),
                pat =['(?<sessionName>[a-z_A-Z]+\d{2,2}[-]\d{8,8})\.'...
                      '(?<mazeName>\w+)\.'...
                      '(?<trialName>\w+)'];
                tok = regexp(Session,pat,'names');
                if ~isempty(tok),
                    Session = tok;
                    Session = MTASession.validate(Session);
                else
                    Session = MTASession(Session);                                    
                end

                
            elseif iscell(Session),
                Session = MTASession(Session{:});
                
            elseif isstruct(Session),
                stcMode = '';
                if isfield(Session,'stcMode'),
                    stcMode = Session.stcMode;
                end

                Session = MTASession(Session.sessionName,...
                                     Session.mazeName,...
                                     Session.trialName);

                if ~isempty(stcMode),
                    Session.load('stc',stcMode);
                end

            else
                error('MTA:validate_trial: unrecognized format');
            end

        end
    end

    
end %MTASession

