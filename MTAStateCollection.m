classdef MTAStateCollection < hgsetget
%MTAStateCollection(Session,varargin) 
%
%  Session - MTASession/string: Either an MTASession containing the paths to the Bhv file
%                               or a string with the complete path filebase, where the
%                               file base follows the convention
%                               [Session.name Session.Maze.name Session.trialName]
%  varargin:
%    [path,filename,mode,sync,origin,overwrite,ext]
%
%    label_mode:  string,  3-4 letter name of the labeling object
%    overwrite:   boolean, flag to overwrite an existing bhv file
%

    properties (SetAccess = public)
        mode = [];
        filename = [];
        path = [];
        sync = [];
        origin = [];
        ext = [];        
        %States - cellArray(MTADepoch): object containing State information         
        states = {};
    end
    
    methods
        function Stc = MTAStateCollection(varargin)
            [path,filename,mode,sync,origin,overwrite,ext] = DefaultArgs(varargin,{[],[],'manual',[],[],0,'stc'}); %#ok<*PROP>
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename = [filename '.' ext '.' mode '.mat'];
                end
            else
                return
            end
            Stc.filename = filename;
            Stc.path = path;
            
            if ~exist(Stc.fpath,'file')||overwrite,
                Stc.mode = mode;
                Stc.sync = sync;
                Stc.origin = origin;
                Stc.ext = ext;
                Stc.save(overwrite);
            else
                Stc.load;
            end
        end

        function Stc = create(Stc,Session,varargin)
            [method,training_set] = DefaultArgs(varargin,{'auto',[]});
            
            Stc.updateFilename([Session.filebase '.stc.' method '.mat']);
            Stc.mode = method;
            Stc.updatePath(Session.spath);
            Stc.ext = 'stc';
            Stc.sync = Session.sync.copy;
            Stc.origin = Session.sync.origin;
            
            switch method
                case 'auto'
                    bhv_auto(Session,Stc);
                case 'hmm'
                    %bhv_hmm(Session,Stc,training_set);
                    disp('bhv_hmm is not ready for use at this time');
                case 'qda'
                    %bhv_qda(Session,Stc,training_set);
                    disp('bhv_qda is not ready for use at this time');
                otherwise
                    error('MTAStateCollection: Unknown method for creating a state collection');
            end
            
        end

        function Stc = load(Stc,varargin)
            [nMode] = DefaultArgs(varargin,{[]});
            if isempty(nMode)&&exist(Stc.fpath,'file')
                ds = load(Stc.fpath);
                %Stc = ds.Stc.copy;
                sprop = properties(Stc);
                for s = 1:numel(sprop)
                    Stc.(sprop{s}) = ds.Stc.(sprop{s});
                end
            else
                Stc.updateMode(nMode);
                Stc.states = {};

                if exist(Stc.fpath,'file')
                sprop = properties(Stc);
                for s = 1:numel(sprop)
                    Stc.(sprop{s}) = ds.Stc.(sprop{s});
                end
                end
            end
        end

        function out = save(Stc,varargin)
            [overwrite] = DefaultArgs(varargin,{0});
            out = false;
            states = Stc.states;  %#ok<NASGU>
            if ~exist(Stc.fpath,'file')
                save( Stc.fpath,'Stc','-v7.3');
                out = true;
            elseif exist(Stc.fpath,'file')&&overwrite
                warning(['Overwriting: ' Stc.fpath]);
                out = true;
                save( Stc.fpath,'Stc','-v7.3');
            else
                warning(['File exists: ' Stc.fpath, ' - flag the overwrite option  to save']);
            end
        end
        
        function fpath = fpath(Stc)
            fpath = fullfile(Stc.path,Stc.filename);
        end
        
        function Stc = rmState(Stc,name)
            %Stc = rmState(Stc,name)
            %Removes a state from the collection 
            %
            % name: string, the label or key of the state to be removed
            %
            if isempty(Stc.states),
                return
            else
                for i = 1:numel(Stc.states)
                    if strcmp(name, Stc.states{i}.label),Stc.states(i) = []; return, end
                    if strcmp(name, Stc.states{i}.key),  Stc.states(i) = []; return, end
                end
            end
        end
        
        function Stc = updateSync(Stc,sync)
            Stc.sync = sync;
        end
        function Stc = updateOrigin(Stc,origin)
            Stc.origin = origin;
        end
        
        function Stc = addState(Stc,path,filename,data,sampleRate,sync,origin,varargin)            
        %Stc = addState(Stsc,key,label,state)
        %
        %   Stc - MTAStateCollection: object containing state epochs
        %
        %   key - char: single character used for keyboard control
        %               of labeling in MTABrowser 
        % 
        %   label - string: one word description of the state
        %
        %   data - numericArray: (event,[start,stop]) Periods
        %           in which the state occurs @ xyzSampleRate
        %
            if isa(path,'MTADepoch'),
                 if ~isempty(Stc.gsi(path.label))&&~isempty(Stc.gsi(path.key)),
                     return
                 end                                      
                Stc.states(end+1) = {path};                
                return
            end
            [label,key,type,ext] = DefaultArgs(varargin,{[],[],[],[]});
            assert(isempty(Stc.gsi(label)),...
                   'MTAStateCollection:addState:ExistingLabel',...
                   ['State: ' label ', already exists in this collection']);
            assert(isempty(Stc.gsi(key)),...
                   'MTAStateCollection:addState:ExistingKey',...
                   ['State: ' label ', already exists in this collection']);
            Stc.states{end+1} = MTADepoch(path,filename,data,sampleRate,sync,origin,label,key,type,ext);
        end
        
        function [Stc,varargout] = subsref(Stc,S)
            varargout = cell(nargout-1,1);
            ni = numel(S);
            for n = 1:ni,
                if isa(Stc,'MTAStateCollection'),
                    if strcmp(S(n).type,'{}'),
                        if ischar(S(n).subs{1}),

                            stsSampleRate = [];
                            if numel(S(n).subs)>=2
                                stsSampleRate = S(n).subs{2};
                                S(n).subs(2) = [];
                            end
                            
                            stsFuncs = regexp(S(n).subs{1},'\&*\^*\|*\+*','match');
                            stsNames = regexp(S(n).subs{1},'\&*\^*\|*\+*','split');
                            assert(numel(stsFuncs)+1==numel(stsNames),...
                                'MTAStateCollection:subsref:UnequalStatesAndFunctionCount',...
                                ['Each state name must be separated by a \n' ...
                                'join(''^'',''&'') or intersect(''|'',''+'') operator']);
                            sts = {};
                            for i = 1:numel(stsNames),
                                stci = Stc.gsi(stsNames{i});
                                if isempty(stci)
                                    if numel(stsNames{i})==1,
                                        sts{i} = MTADepoch(Stc.path,[],[],[],[],[],[],stsNames{i},[],[]);
                                    else
                                        sts{i} = MTADepoch(Stc.path,[],[],[],[],[],stsNames{i},[],[],[]);
                                    end
                                    try,
                                        sts{i} = sts{i}.load(Stc.sync);
                                        Stc.addState(sts{i});
                                    catch
                                        warning(['MTAStateCollection:subrefs:StateNotFound, state: ''' stsNames{i} ''' does not exist']);
                                        sts{i} = {};
                                    end

                                else
                                    sts{i} =  Stc.states{stci}.copy;
                                end
                            end
                            
                            while ~isempty(stsFuncs),
                                switch stsFuncs{1}
                                  case '&'
                                    sts{2} = MTADepoch.intersect(sts(1:2));
                                    sts(1) = [];
                                  case '^'
                                    sts{2} = MTADepoch.intersect(sts(1:2));
                                    sts(1) = [];
                                  case '+'
                                    sts{2} = MTADepoch.join(sts(1:2));
                                    sts(1) = [];
                                  case '|'
                                    sts{2} = MTADepoch.join(sts(1:2));
                                    sts(1) = [];
                                end
                                stsFuncs(1) = [];
                            end
                            
                            if ~isempty(sts{1}),
                                sts = sts{1}.copy;
                                if ~isempty(stsSampleRate)
                                    sts.resample(stsSampleRate);
                                end
                                if isempty(Stc.gsi(sts.label)),
                                    Stc.addState(sts);
                                end
                                Stc = sts;
                            else
                                Stc = {};
                            end

                        else
                            Stc = Stc.states{S(n).subs{1}}.copy;
                            Stc.resample(S(n).subs{2});
                        end
                    else
                        Stc = builtin('subsref',Stc,S(n:end));                       
                        return
                    end
                else
                    Stc = subsref(Stc,S(n:end));
                    return
                end
            end
        end

        function Stc = updatePath(Stc,path)
            Stc.path = path;
        end

        function Stc = updateMode(Stc,mode)
            Stc.mode = mode;
            [~,fname,fext] = fileparts(Stc.filename);
            pind = strfind(fliplr(fname),'.');
            fname = [fname(1:end-pind+1) mode fext];
            Stc.updateFilename(fname);
        end
        
        function Stc = updateFilename(Stc,filename)
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename = [filename '.' Stc.ext '.' Stc.mode '.mat'];
                end
            end
            Stc.filename = filename;
        end
        
        function DataCopy = copy(Data)
        % Make a copy of a handle object.
        % Instantiate new object of the same class.
            DataCopy = feval(class(Data),[]);
            % Copy all non-hidden properties.
            p = properties(Data);
            for i = 1:length(p)                
                if isa(Data.(p{i}),'MTAData'),
                    DataCopy.(p{i}) = Data.(p{i}).copy;
                else
                    DataCopy.(p{i}) = Data.(p{i});
                end
            end
            if ~isempty(Data.states)
                for s = 1:numel(Data.states)
                    DataCopy.states{s} = Data.states{s}.copy;
                end
            else
                DataCopy.states = Data.states;
            end    
        end
        
        function state_index = gsi(Stc,state)
            state_index = [];
            for i = 1:numel(Stc.states),
                if strcmp(Stc.states{i}.label,state)||strcmp(Stc.states{i}.key,state),
                    state_index = i;
                    return;
                end
            end
        end

        function state_attrib_list = list_state_attrib(Stc,varargin)
        % state_attrib_list = list_state_attrib(Bhv,attrib)
        %
        % available attributes to be listed
        % key - string: keyboard character associated with state
        % label - string: name describing state (e.g. 'rear' or 'walk')
        % state - numericArray<-double: (event,[time_onset,time_offset]) @xyzSampleRate
 
        [attrib] = DefaultArgs(varargin,{'label'});

            state_attrib_list = {};
            for i = 1:length(Stc.states),
                state_attrib_list{i} = Stc.states{i}.(attrib);
            end
        end

        function [sts_res,filterName] = filter(Stc,sampleRate,command_list)
        %  use cell arrays for message passing,
        % The first index of the command list i
        %  order of commands will be important: make rules
        % {'walk',{'exclude',{'rear','immobility'}},{'trim',2}}
            ENFT.type = 'MTABhv:filter:NoFilterTarget';
            ENFT.msg  = ['Behavior: %s, doesn''''t exist. No target ' ...
                         'periods for filtration.'];
            sts_periods = [];
            tname = '';

            if numel(command_list(1))>1
                switch command_list{1}{1}
                  case 'join'
                    tname = 'bhvj';
                    assert(~isempty(Stc.subsref(substruct('{}',{command_list{1}{2}{1},sampleRate}))),ENFT.type,ENFT.msg,command_list{1}{2}{1});
                    bhvp{1} = Stc.subsref(substruct('{}',{command_list{1}{2}{1},sampleRate},'.','copy'));
                    tname = strcat(tname,'_',command_list{1}{2}{1});
                    for s = 2:numel(command_list{1}{2}),                                          
                        assert(~isempty(Stc.subsref(substruct('{}',{command_list{1}{2}{s},sampleRate}))),ENFT.type,ENFT.msg,command_list{1}{2}{s});
                        bhvp{s} = Stc.subsref(substruct('{}',{command_list{1}{2}{s},sampleRate},'.','data'));                                                                        
                        tname = strcat(tname,'_',command_list{1}{2}{1});
                    end
                    sts_periods = JoinRanges(bhvp);

                  case 'intersect' 
                    tname = 'bhvi';         
                    assert(~isempty(Stc.subsref(substruct('{}',{command_list{1}{2}{1},sampleRate},'.','data'))),ENFT.type,ENFT.msg,command_list{1}{2}{1});
                    sts_periods = Stc.subsref(substruct('{}',{command_list{1}{2}{1},sampleRate},'.','data'));
                    tname = strcat(tname,'_',command_list{1}{2}{1});
                    for s = 2:numel(command_list{1}{2}),                                          
                        assert(~isempty(Stc.subsref(substruct('{}',{command_list{1}{2}{s},sampleRate},'.','data'))),ENFT.type,ENFT.msg,command_list{1}{2}{s});
                        sts_periods = IntersectRanges(sts_periods,Stc.subsref(substruct('{}',{command_list{1}{2}{s},sampleRate},'.','data')));
                        tname = strcat(tname,'_',command_list{1}{2}{1});
                    end
                    %% Subsref fixed upto this point

                  otherwise
                    error(['Commmand: ' command_list{1}{1} ' not found.'])

                end
            else 
                tname = strcat('Stcs_',command_list{1});
                assert(~isempty(Stc.subsref(substruct('{}',{command_list{1},sampleRate},'.','data'))),ENFT.type,ENFT.msg,command_list{1});
                sts_periods = Stc.subsref(substruct('{}',{command_list{1},sampleRate},'.','data'));
            end

            Sts_onset_ind  = true(size(sts_periods,1),1);
            Sts_offset_ind = true(size(sts_periods,1),1);

            %% Remove Target from the command list
            command_list(1)=[];

            filterName = ['F_' tname '_'];
            tshift = 0;

            while numel(command_list)~=0,
                switch command_list{1}{1},
                  
                  case 'exclusion',  filter_tag = 'e';                  
                  
                    if iscell(command_list{1}{2}),
                        ExState = Stc.JoinStates(command_list{1}{2});                       
                    else
                        ExState = Stc.subsref(substruct('{}',{command_list{1}{2}},'.','copy'));
                    end
                    ExState.resample(sampleRate);
                    filterName = strcat(filterName,filter_tag,ExState.key,num2str(command_list{1}{3}));                    
                    ExThresh = command_list{1}{3}*sampleRate;

                    temp_sts_onset_ind  = false(size(sts_periods,1),1);
                    temp_sts_offset_ind = false(size(sts_periods,1),1);
                    event_prox_ind = [];
                    for i = 1:length(Sts_onset_ind)
                        event_prox = ExState.data(:,2)-sts_periods(i,1);
                        event_prox_ind = find(event_prox>-ExThresh & event_prox < 0,1,'first');
                        if isempty(event_prox_ind),
                            temp_sts_onset_ind(i)  = true;
                        end
                    end
                    for i = 1:length(Sts_offset_ind)
                        event_prox = ExState.data(:,1)-sts_periods(i,2);
                        event_prox_ind = find(event_prox< ExThresh & event_prox > 0,1,'first');
                        if isempty(event_prox_ind),
                            temp_sts_offset_ind(i)  = true;
                        end
                    end
                    Sts_onset_ind  = Sts_onset_ind & temp_sts_onset_ind;
                    Sts_offset_ind = Sts_offset_ind & temp_sts_offset_ind;

                  case 'select_boarder_states', filter_tag = 's';
                    %% Spelect event onset/offset based on the
                    %% boardering states                    
                    if iscell(command_list{1}{2}),
                        ExState = Stc.JoinStates(command_list{1}{2});
                    else
                        ExState = Stc.subsref(substruct('{}',{command_list{1}{2}},'.','copy'));
                    end
                    ExState.resample(sampleRate);
                    filterName = strcat(filterName,filter_tag,ExState.key,num2str(command_list{1}{3}));
                    ExThresh = command_list{1}{3}*sampleRate;

                    temp_sts_onset_ind  = false(size(sts_periods,1),1);
                    temp_sts_offset_ind = false(size(sts_periods,1),1);                    
                    event_prox_ind = [];
                    for i = 1:length(Sts_onset_ind)
                        event_prox = ExState.data(:,2)-sts_periods(i,1);
                        event_prox_ind = find(event_prox>-ExThresh & event_prox < 0,1,'first');
                        if ~isempty(event_prox_ind),
                            temp_sts_onset_ind(i)  = true;
                        end
                    end
                    for i = 1:length(Sts_offset_ind)
                        event_prox = ExState.data(:,1)-sts_periods(i,2);
                        event_prox_ind = find(event_prox< ExThresh & event_prox > 0,1,'first');
                        if ~isempty(event_prox_ind),
                            temp_sts_offset_ind(i)  = true;
                        end
                    end
                    Sts_onset_ind  = Sts_onset_ind & temp_sts_onset_ind;
                    Sts_offset_ind = Sts_offset_ind & temp_sts_offset_ind;


                  case 'duration', filter_tag = 'd';

                    filterName = strcat(filterName,filter_tag,num2str(command_list{1}{2}));
                    bdur = diff(sts_periods,1,2)/sampleRate;

                    temp_sts_onset_ind  = false(size(sts_periods,1),1);
                    temp_sts_offset_ind = false(size(sts_periods,1),1);
                    temp_sts_onset_ind (bdur>command_list{1}{2}) = true; 
                    temp_sts_offset_ind(bdur>command_list{1}{2}) = true;
                    Sts_onset_ind  = Sts_onset_ind & temp_sts_onset_ind;
                    Sts_offset_ind = Sts_offset_ind & temp_sts_offset_ind;
                    

                  case 'trim',  filter_tag = 't';
              
                    filterName = strcat(filterName,filter_tag,num2str(command_list{1}{2}));
                    tshift = round(command_list{1}{2}*sampleRate);
                    
                  case 'complete', filter_tag = 'c';

                    filterName = strcat(filterName,filter_tag);
                    temp_sts_onset_ind  = Sts_offset_ind;
                    temp_sts_offset_ind = Sts_onset_ind;
                    Sts_onset_ind  = Sts_onset_ind & temp_sts_onset_ind;
                    Sts_offset_ind = Sts_offset_ind & temp_sts_offset_ind;
                                             
                  otherwise
                    error(['Commmand: ' command_list{1}{1} ' not found.'])
                end
    
                command_list(1)=[];
            end    
                        

            filterName(ismember(filterName,'.'))='_';
            sts_res{1}  = sts_periods(Sts_onset_ind, 1)+tshift;
            sts_res{2}  = sts_periods(Sts_offset_ind,2)-tshift;
            

        end

        function composite_state = JoinStates(Stc,varargin)
        % JoinStates(Bhv,varargin)
        % varargin can be a list of state labels
        % or a single cell array of state labels
        % The final sample rate is set by the first listed state.
            if numel(varargin)==1,
                states = varargin{1};
            else
                states = varargin;
            end
            
            if numel(states) == 1,
                composite_state = Stc.subsref(substruct('{}',states(1)));
            else
                sampleRate = Stc.subsref(substruct('{}',states(1),'.','sampleRate'));
                oper = cell(numel(states),1);
                keys = '';
                for i=1:numel(states)
                    oper{i} = Stc.subsref(substruct('{}',{states{i},sampleRate},'.','data'));
                    keys(end+1) = Stc.subsref(substruct('{}',states(i),'.','key'));
                end
                newStateName = ['COMP_' keys];
                uper = JoinRanges(oper);
                composite_state = MTADepoch([],uper,sampleRate, ...
                                            states{1}.sync.copy,states{1}.origin,newStateName,keys,[],[]);
            end
        end
        
        function out = isempty(Stc)
            out = isempty(Stc.states);
        end

    end
end
