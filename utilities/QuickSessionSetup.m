function QuickSessionSetup(SessionParm,varargin)
% function QuickSessionSetup(SessionParm,varargin)
% HostConf = load('MTAConf');
% [host,local] = DefaultArgs(varargin,{HostConf.host_server,false},1);

HostConf = load('MTAConf');
[host,local,link] = DefaultArgs(varargin,{HostConf.host_server,false,true},1);
    
        
    % Expecting Name of Session list, cell, or struct
    if ischar(SessionParm), 
        Sessions = get_session_list(SessionParm);
    elseif iscell(SessionParm),
        Sessions = get_session_list(SessionParm{:});
    else
        Sessions = SessionParm;
    end

    
    % SessionList Struct Type
    if isstruct(Sessions),
        for s = 1:numel(Sessions)
            if local,
                MTAstartup(host,HostConf.data_server);
            else
                MTAstartup(host,Sessions(s).host);
            end
                
            if all(isfield(Sessions(s),{'xyz_host','nlx_host'}))&&~ispc&&link,
                linkSession(Sessions(s).sessionName,...
                            Sessions(s).xyz_host,...
                            Sessions(s).nlx_host);
                dataLoggers = {'vicon','nlx'};
            elseif isfield(Sessions(s),'xyz_host')&&~ispc&&link,
                linkSession(Sessions(s).sessionName,...
                            Sessions(s).xyz_host);
                dataLoggers = {'vicon'};
            end
            
            % Check for MoCap SampleRate
            if isfield(Sessions(s),'xyzSampleRate'),
                xyzSampleRate = Sessions(s).xyzSampleRate;
            else
                xyzSampleRate = [];
            end
            
            % Create/Overwrite Session
            Session = MTASession(Sessions(s).sessionName,...
                                 Sessions(s).mazeName,...
                                 true,...
                                 Sessions(s).TTLValue,...
                                 dataLoggers,...
                                 'xyzSampleRate',xyzSampleRate);
        end
    
    % SessionList Cell Type
    elseif iscell(Sessions),    
        for s = 1:numel(Sessions)
            % Setup MTA
            MTAstartup(host,Sessions{s}{4});
            % Create/Overwrite Session
            Session = MTASession(Sessions{s}{1},Sessions{s}{2},true,Sessions{s}{5});
        end

    end


end
