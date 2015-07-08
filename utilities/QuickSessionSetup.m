function QuickSessionSetup(SessionParm,varargin)
HostConf = load('MTAConf');
[host,local] = DefaultArgs(varargin,{HostConf.host_server,false},1);
    
    
        
    if ischar(SessionParm), 
        Sessions = SessionList(SessionParm);
    else
        Sessions = {SessionParm};
    end

    if isstruct(Sessions),
        for s = 1:numel(Sessions)
            if local,
                MTAstartup(host,HostConf.data_server);
            else
                MTAstartup(host,Sessions(s).host);
            end
                
            if all(isfield(Sessions(s),{'xyz_host','nlx_host'}))&&~ispc,
                linkSession(Sessions(s).name,...
                            Sessions(s).xyz_host,...
                            Sessions(s).nlx_host);
            end
            
            if isfield(Sessions(s),'xyzSampleRate'),
                xyzSampleRate = Sessions(s).xyzSampleRate;
            else
                xyzSampleRate = [];
            end
            
            Session = MTASession(Sessions(s).name,...
                                 Sessions(s).mazeName,...
                                 true,...
                                 Sessions(s).TTLValue,...
                                 'xyzSampleRate',xyzSampleRate);
        end
    elseif iscell(Sessions),    
        for s = 1:numel(Sessions)
            MTAstartup(host,Sessions{s}{4});
            Session = MTASession(Sessions{s}{1},Sessions{s}{2},true,Sessions{s}{5});
        end

    end


end
