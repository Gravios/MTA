classdef MTAStateCollection < hgsetget
%MTAStateCollection(Session,varargin) 
%
%  Session - MTASession/string: Either an MTASession containing the paths to the Bhv file
%                               or a string with the complete path filebase, where the
%                               file base follows the convention
%                               [Session.name Session.Maze.name Session.trialName]
%  varargin:
%    [label_mode,overwrite]
%
%    label_mode:  string,  3-4 letter name of the labeling object
%    overwrite:   boolean, flag to overwrite an existing bhv file
%

    properties (SetAccess = public)
        mode = [];
        filename = [];
        path = [];
        ext = [];        
    end
    
    properties (Transient=true)
        %States - cellArray(MTADepoch): object containing State information         
        states = {};
    end

    methods
        function Stc = MTAStateCollection(varargin)
            [path,filename,mode,overwrite,ext] = DefaultArgs(varargin,{[],[],'manual',0,'stc'}); %#ok<*PROP>
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
                Stc.states  = {};
                Stc.ext = ext;
                Stc.save(overwrite);
            else
                Stc = Stc.load;
            end
        end
      
        function Stc = load(Stc)
            ds = load(Stc.fpath);
            Stc = ds.Stc;
            Stc.states = {};
        end
        
        function out = save(Stc,overwrite)
            out = false;
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
        
        function Stc = addState(Stc,path,filename,data,sampleRate,varargin)            
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
                Stc.states(end+1) = {path};                
                path.save;
                return
            end
            [key,label,type,ext] = DefaultArgs(varargin,{[],[],[],[]});
            assert(isempty(Stc.gsi(label)),...
                   'MTAStateCollection:addState:ExistingLabel',...
                   ['State: ' label ', already exists in this collection']);
            assert(isempty(Stc.gsi(key)),...
                   'MTAStateCollection:addState:ExistingKey',...
                   ['State: ' label ', already exists in this collection']);
            Stc.states{end+1} = MTADepoch(path,filename,data,sampleRate,type,ext,label,key);
            Stc.states{end}.save;
        end
        
        function [Stc,varargout] = subsref(Stc,S)
            varargout = cell(nargout-1,1);
            ni = numel(S);
            for n = 1:ni,
                if isa(Stc,'MTAStateCollection'),
                    if strcmp(S(n).type,'{}'),
                        if ischar(S(n).subs{1}),
                            stci = Stc.gsi(S(n).subs{1});
                            if isempty(stci)
                                if numel(S(n).subs{1})==1,
                                    sts = MTADepoch(Stc.path,[],[],[],[],[],S(n).subs{1},[]);
                                else
                                    sts = MTADepoch(Stc.path,[],[],[],[],[],[],S(n).subs{1});
                                end
                                sts = sts.load;
                                Stc.addState(sts);
                                Stc = sts;
                                %return
                            else
                                Stc =  Stc.states{stci};
                                %return
                            end
                        else
                            Stc = Stc.states{S(n).subs{1}};
                        end
%                     elseif strcmp(S(n).type,'()'),
%                         if ischar(S(n).subs{1}),
%                             if strcmp(S(n).subs{1},':')
%                                 Stc = Stc.states;
%                                 return
%                             end
%                             stci = Stc.gsi(S(n).subs{1});
%                             if isempty(stci)
%                                 sts = MTADepoch(Stc.path,[],[],[],[],S(n).subs{1},[]);
%                                 sts = sts.load;
%                                 Stc.addState(sts);
%                                 Stc = sts;
%                             else
%                                 Stc = Stc.states{stci};
%                             end
%                         else
%                             Stc = Stc.states{S(n).subs{1}};
%                         end
                        
                    else
%                         dbreak = false;
%                         if ismethod(Stc,S(n).subs),dbreak = true;end
                        Stc = builtin('subsref',Stc,S(n:end));
                        return
%                         if dbreak,break,end
                    end
                else
                    Stc = subsref(Stc,S(n:end));
                end
            end
        end
        
        function Stc = updatePath(Stc,path)
            Stc.path = path;
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
                DataCopy.(p{i}) = Data.(p{i});
            end
            for s = 1:numel(Data.states)
                DataCopy.states{s} = Data.states{s}.copy;
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

        function [sts_res,filterName] = filter(Stc,command_list)
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
                    assert(~isempty(Stc(command_list{1}{2}{1})),ENFT.type,ENFT.msg,command_list{1}{2}{1});
                    bhvp{1} = Stc(command_list{1}{2}{1});
                    tname = strcat(tname,'_',command_list{1}{2}{1});
                    for s = 2:numel(command_list{1}{2}),                                          
                        assert(~isempty(Stc(command_list{1}{2}{s})),ENFT.type,ENFT.msg,command_list{1}{2}{s});
                        bhvp{s} = Stc(command_list{1}{2}{s}).data;                                                                        
                        tname = strcat(tname,'_',command_list{1}{2}{1});
                    end
                    sts_periods = JoinRanges(bhvp);

                  case 'intersect' 
                    tname = 'bhvi';         
                    assert(~isempty(Stc(command_list{1}{2}{1}).data),ENFT.type,ENFT.msg,command_list{1}{2}{1});
                    sts_periods = Stc(command_list{1}{2}{1}).data;
                    tname = strcat(tname,'_',command_list{1}{2}{1});
                    for s = 2:numel(command_list{1}{2}),                                          
                        assert(~isempty(Stc(command_list{1}{2}{s}).data),ENFT.type,ENFT.msg,command_list{1}{2}{s});
                        sts_periods = IntersectRanges(sts_periods,Stc(command_list{1}{2}{s}).data);
                        tname = strcat(tname,'_',command_list{1}{2}{1});
                    end

                  otherwise
                    error(['Commmand: ' command_list{1}{1} ' not found.'])

                end
            else 
                tname = strcat('Stcs_',command_list{1});
                assert(~isempty(Stc(command_list{1}).data),ENFT.type,ENFT.msg,command_list{1});
                sts_periods = Stc(command_list{1}).data;
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
                        ExState = Stc(command_list{1}{2});
                    end
                    filterName = strcat(filterName,filter_tag,ExState.key,num2str(command_list{1}{3}));                    
                    ExThresh = command_list{1}{3}*Stc.sampleRate;

                    temp_sts_onset_ind  = false(size(sts_periods,1),1);
                    temp_sts_offset_ind = false(size(sts_periods,1),1);
                    event_prox_ind = [];
                    for i = 1:length(Sts_onset_ind)
                        event_prox = ExState.state(:,2)-sts_periods(i,1);
                        event_prox_ind = find(event_prox>-ExThresh & event_prox < 0,1,'first');
                        if isempty(event_prox_ind),
                            temp_sts_onset_ind(i)  = true;
                        end
                    end
                    for i = 1:length(Sts_offset_ind)
                        event_prox = ExState.state(:,1)-sts_periods(i,2);
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
                        ExState = Stc.getState(command_list{1}{2});
                    end
                    filterName = strcat(filterName,filter_tag,ExState.key,num2str(command_list{1}{3}));
                    ExThresh = command_list{1}{3}*Stc.sampleRate;

                    temp_sts_onset_ind  = false(size(sts_periods,1),1);
                    temp_sts_offset_ind = false(size(sts_periods,1),1);                    
                    event_prox_ind = [];
                    for i = 1:length(Sts_onset_ind)
                        event_prox = ExState.state(:,2)-sts_periods(i,1);
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
                    bdur = diff(sts_periods,1,2)/Stc.sampleRate;

                    temp_sts_onset_ind  = false(size(sts_periods,1),1);
                    temp_sts_offset_ind = false(size(sts_periods,1),1);
                    temp_sts_onset_ind (bdur>command_list{1}{2}) = true; 
                    temp_sts_offset_ind(bdur>command_list{1}{2}) = true;
                    Sts_onset_ind  = Sts_onset_ind & temp_sts_onset_ind;
                    Sts_offset_ind = Sts_offset_ind & temp_sts_offset_ind;
                    

                  case 'trim',  filter_tag = 't';
              
                    filterName = strcat(filterName,filter_tag,num2str(command_list{1}{2}));
                    tshift = round(command_list{1}{2}*Stc.sampleRate);
                    
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
                        
            pind = find(filterName=='.');
            filterName(pind) = '_';
            sts_res{1}  = sts_periods(Sts_onset_ind, 1)+tshift;
            sts_res{2}  = sts_periods(Sts_offset_ind,2)-tshift;
            

        end

        function composite_state = JoinStates(Stc,varargin)
        % JoinStates(Bhv,varargin)
        % varargin can be a list of state labels
        % or a single cell array of state labels
            if numel(varargin)==1,
                states = varargin{1};
            else
                states = varargin;
            end

            oper = cell(numel(states),1);
            for i=1:numel(states)
                oper{i} = cat(1,sper,Stc(states{i}).data);
                keys{end+1} = Stc(states{i}).key;
            end
            newStateName = ['COMP_' [keys{:}]];
            uper = JoinRanges(oper);
            composite_state = MTADepoch([],uper,Stc.sampleRate,[],[],[keys{:}],newStateName);
        end
        
        function out = isempty(Stc)
            out = isempty(Stc.states);
        end

    end
end
