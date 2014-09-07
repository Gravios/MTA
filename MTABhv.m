classdef MTABhv < hgsetget
%MTABhv(Session,varargin) 
% DO NOT USE: OLD
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

        %mode - string: descriptor often associated with xml file, specifying heuristic functions and parameters, if no xml file exists then MTABhv object was manually labeled.
        mode

        %sampleRate - double: sampling rate of the state periods
        sampleRate

        %States - cellArray(MTAState): object containing State information 
        States        
    end

    methods
        function Bhv = MTABhv(Session,varargin)
            [label_mode,overwrite] = DefaultArgs(varargin,{'manual',0});
            if isa(Session,'MTASession'),
                bhvfile = [Session.filebase '.bhv.' label_mode '.mat'];
                bhvpath = Session.spath;
            end
            
            if isempty(bhvpath),                
                Bhv.mode = label_mode;
                Bhv.States  = {};            
            elseif ~exist(fullfile(bhvpath,bhvfile),'file')||overwrite,
                Bhv.mode = label_mode;
                Bhv.States  = {};
                Bhv.save(bhvpath,overwrite);
            else
                load(fullfile(bhvpath,bhvfile));
            end
        end

        function save(Bhv,Session,overwrite)
            if isa(Session,'MTASession'),                                
                bhvpath = fullfile(Session.spath, [Session.filebase '.bhv.' Bhv.mode '.mat']);
            end

            if ~exist(bhvpath,'file')
                save( bhvpath,'Bhv','-v7.3');
            elseif exist(bhvpath,'file')&&overwrite
                warning(['Overwriting: ' bhvpath]);
                save( bhvpath,'Bhv','-v7.3');
            else
                warning(['File exists: ' bhvpath, ' - flag the overwrite option  to save']);
            end
        end      

        function Bhv = addState(Bhv,key,label,state)
        %Bhv = addState(Bhv,key,label,state)
        %
        %   Bhv - MTABhv: object containing behavioral states
        %
        %   key - char: single character used for keyboard control
        %               of labeling in MTABrowser 
        % 
        %   label - string: one word description of the state
        %
        %   state - numericArray: (event,[start,stop]) Periods
        %           in which the state occurs @ xyzSampleRate
        %
            Bhv.States(end+1) = {MTAState(key,label,state)};
        end

        function Bhv = subsref(Bhv,S)
        %State = getState(Bhv,label)
        %
        %   Bhv - MTABhv: object containing behavioral states
        %
        %   index - string: one word description of the state which
        %                   will be retrieved 
        %
            ni = numel(S);
            for n = 1:ni,
                if isa(Bhv,'MTABhv'),
                    if strcmp(S(n).type,'()'),
                        if ischar(S(n).subs{1}),
                            Bhv = Bhv.States{Bhv.gsi(S(n).subs{1})};
                        else
                            Bhv = Bhv.States{S(n).subs{1}};
                        end
                    else
                        Bhv = builtin('subsref',Bhv,S(n));
                    end
                else
                    Bhv = subsref(Bhv,S(n));
                end
            end
        end
        
        function state_index = gsi(Bhv,state)
            state_index = [];
            for i = 1:length(Bhv.States),
                if strcmp(Bhv.States{i}.label,state),
                    state_index = i;
                    return;
                end
            end
        end

        function state_attrib_list = list_state_attrib(Bhv,varargin)
        % state_attrib_list = list_state_attrib(Bhv,attrib)
        %
        % available attributes to be listed
        % key - string: keyboard character associated with state
        % label - string: name describing state (e.g. 'rear' or 'walk')
        % state - numericArray<-double: (event,[time_onset,time_offset]) @xyzSampleRate
 
        [attrib] = DefaultArgs(varargin,{'label'});

            state_attrib_list = {};
            for i = 1:length(Bhv.States),
                state_attrib_list{i} = Bhv.States{i}.(attrib);
            end
        end

        function [bhv_res,filterName] = filter(Bhv,command_list)
        %  use cell arrays for message passing,
        % The first index of the command list i
        %  order of commands will be important: make rules
        % {'walk',{'exclude',{'rear','immobility'}},{'trim',2}}
            ENFT.type = 'MTABhv:filter:NoFilterTarget';
            ENFT.msg  = ['Behavior: %s, doesn''''t exist. No target ' ...
                         'periods for filtration.'];
            bhv_periods = [];
            tname = '';

            if numel(command_list(1))>1
                switch command_list{1}{1}
                  case 'join'
                    tname = 'bhvj';
                    assert(~isempty(Bhv.getState(command_list{1}{2}{1}).state),ENFT.type,ENFT.msg,command_list{1}{2}{1});
                    bhvp{1} = Bhv.getState(command_list{1}{2}{1}).state;
                    tname = strcat(tname,'_',command_list{1}{2}{1});
                    for s = 2:numel(command_list{1}{2}),                                          
                        assert(~isempty(Bhv.getState(command_list{1}{2}{s}).state),ENFT.type,ENFT.msg,command_list{1}{2}{s});
                        bhvp{s} = Bhv.getState(command_list{1}{2}{s}).state;                                                                        
                        tname = strcat(tname,'_',command_list{1}{2}{1});
                    end
                    bhv_periods = JoinRanges(bhvp);

                  case 'intersect' 
                    tname = 'bhvi';         
                    assert(~isempty(Bhv.getState(command_list{1}{2}{1}).state),ENFT.type,ENFT.msg,command_list{1}{2}{1});
                    bhv_periods = Bhv.getState(command_list{1}{2}{1}).state;
                    tname = strcat(tname,'_',command_list{1}{2}{1});
                    for s = 2:numel(command_list{1}{2}),                                          
                        assert(~isempty(Bhv.getState(command_list{1}{2}{s}).state),ENFT.type,ENFT.msg,command_list{1}{2}{s});
                        bhv_periods = IntersectRanges(bhv_periods,Bhv.getState(command_list{1}{2}{s}).state);
                        tname = strcat(tname,'_',command_list{1}{2}{1});
                    end

                  otherwise
                    error(['Commmand: ' command_list{1}{1} ' not found.'])

                end
            else 
                tname = strcat('bhvs_',command_list{1});
                assert(~isempty(Bhv.getState(command_list{1}).state),ENFT.type,ENFT.msg,command_list{1});
                bhv_periods = Bhv.getState(command_list{1}).state;
            end

            bhv_onset_ind  = true(size(bhv_periods,1),1);
            bhv_offset_ind = true(size(bhv_periods,1),1);

            %% Remove Target from the command list
            command_list(1)=[];

            filterName = ['f_' tname '_'];
            tshift = 0;


            while numel(command_list)~=0,
                switch command_list{1}{1},
                  
                  case 'exclusion',  filter_tag = 'e';                  
                  
                    if iscell(command_list{1}{2}),
                        ExState = Bhv.JoinStates(command_list{1}{2});
                    else
                        ExState = Bhv.getState(command_list{1}{2});
                    end
                    filterName = strcat(filterName,filter_tag,ExState.key,num2str(command_list{1}{3}));                    
                    ExThresh = command_list{1}{3}*Bhv.sampleRate;

                    temp_bhv_onset_ind  = false(size(bhv_periods,1),1);
                    temp_bhv_offset_ind = false(size(bhv_periods,1),1);
                    event_prox_ind = [];
                    for i = 1:length(bhv_onset_ind)
                        event_prox = ExState.state(:,2)-bhv_periods(i,1);
                        event_prox_ind = find(event_prox>-ExThresh & event_prox < 0,1,'first');
                        if isempty(event_prox_ind),
                            temp_bhv_onset_ind(i)  = true;
                        end
                    end
                    for i = 1:length(bhv_offset_ind)
                        event_prox = ExState.state(:,1)-bhv_periods(i,2);
                        event_prox_ind = find(event_prox< ExThresh & event_prox > 0,1,'first');
                        if isempty(event_prox_ind),
                            temp_bhv_offset_ind(i)  = true;
                        end
                    end
                    bhv_onset_ind  = bhv_onset_ind & temp_bhv_onset_ind;
                    bhv_offset_ind = bhv_offset_ind & temp_bhv_offset_ind;

                  case 'select_boarder_states', filter_tag = 's';
                    %% Spelect event onset/offset based on the
                    %% boardering states                    
                    if iscell(command_list{1}{2}),
                        ExState = Bhv.JoinStates(command_list{1}{2});
                    else
                        ExState = Bhv.getState(command_list{1}{2});
                    end
                    filterName = strcat(filterName,filter_tag,ExState.key,num2str(command_list{1}{3}));
                    ExThresh = command_list{1}{3}*Bhv.sampleRate;

                    temp_bhv_onset_ind  = false(size(bhv_periods,1),1);
                    temp_bhv_offset_ind = false(size(bhv_periods,1),1);                    
                    event_prox_ind = [];
                    for i = 1:length(bhv_onset_ind)
                        event_prox = ExState.state(:,2)-bhv_periods(i,1);
                        event_prox_ind = find(event_prox>-ExThresh & event_prox < 0,1,'first');
                        if ~isempty(event_prox_ind),
                            temp_bhv_onset_ind(i)  = true;
                        end
                    end
                    for i = 1:length(bhv_offset_ind)
                        event_prox = ExState.state(:,1)-bhv_periods(i,2);
                        event_prox_ind = find(event_prox< ExThresh & event_prox > 0,1,'first');
                        if ~isempty(event_prox_ind),
                            temp_bhv_offset_ind(i)  = true;
                        end
                    end
                    bhv_onset_ind  = bhv_onset_ind & temp_bhv_onset_ind;
                    bhv_offset_ind = bhv_offset_ind & temp_bhv_offset_ind;


                  case 'duration', filter_tag = 'd';

                    filterName = strcat(filterName,filter_tag,num2str(command_list{1}{2}));
                    bdur = diff(bhv_periods,1,2)/Bhv.sampleRate;

                    temp_bhv_onset_ind  = false(size(bhv_periods,1),1);
                    temp_bhv_offset_ind = false(size(bhv_periods,1),1);
                    temp_bhv_onset_ind (bdur>command_list{1}{2}) = true; 
                    temp_bhv_offset_ind(bdur>command_list{1}{2}) = true;
                    bhv_onset_ind  = bhv_onset_ind & temp_bhv_onset_ind;
                    bhv_offset_ind = bhv_offset_ind & temp_bhv_offset_ind;
                    

                  case 'trim',  filter_tag = 't';
              
                    filterName = strcat(filterName,filter_tag,num2str(command_list{1}{2}));
                    tshift = round(command_list{1}{2}*Bhv.sampleRate);
                    
                  case 'complete', filter_tag = 'c';

                    filterName = strcat(filterName,filter_tag);
                    temp_bhv_onset_ind  = bhv_offset_ind;
                    temp_bhv_offset_ind = bhv_onset_ind;
                    bhv_onset_ind  = bhv_onset_ind & temp_bhv_onset_ind;
                    bhv_offset_ind = bhv_offset_ind & temp_bhv_offset_ind;
                                             
                  otherwise
                    error(['Commmand: ' command_list{1}{1} ' not found.'])
                end
    
                command_list(1)=[];
            end    
                        
            pind = find(filterName=='.');
            filterName(pind) = '_';
            bhv_res{1}  = bhv_periods(bhv_onset_ind, 1)+tshift;
            bhv_res{2}  = bhv_periods(bhv_offset_ind,2)-tshift;
            

        end

        function composite_state = JoinStates(Bhv,varargin)
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
                oper{i} = cat(1,sper,Bhv.getState(states{i}).state);
                keys{end+1} = Bhv.getState(states{i}).key;
            end
            newStateName = ['u_' [keys{:}]];
            uper = JoinRanges(oper);
            composite_state = MTAState([keys{:}],newStateName,uper);
        end

    end
end
