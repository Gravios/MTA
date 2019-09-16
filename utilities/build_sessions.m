function build_sessions(sessionListName,varargin)
% function build_sessions(sessionListName,overwrite)
%
%
% sessionList is a struct array which requires the * following
% fields to function
%
%-----------------------------------------------------------------------
%  req   Field            Type      Units          Example             |
%-----------------------------------------------------------------------      
%   *    sessionName:     String,   NA,            'jg05-20120311'     |
%   *    mazeName:        String,   NA,            'cof'               |
%   *    trialName:       String,   NA,            'all'               |
%   *    xyz_host:        String,   NA,            '/path/to/xyz/data' |
%   *    nlx_host:        String,   NA,            '/path/to/nlx/data' |
%   *    xyzSampleRate:   Numeric   frames/second   119.881035         |
%   *    host:            String,   NA,            'lmu'               |
%   *    project',        String,   NA,            'general'           |
%   *    TTLValue',       String,   NA,            'Vicon start'       |
%        includeSyncInd', Integers, NA,            []                  |
%        offsets',        Numeric,  Seconds,       [15,0]              |
%   *    xOffSet',        Numeric,  mm,             0                  |                    
%   *    yOffSet',        Numeric,  mm,             0                  |
%        stcMode',        String,   NA,            'default'           |
%-----------------------------------------------------------------------
%   
    
    
% DEFARGS ------------------------------------------------------------------------------------------
Config = load('MTAConf.mat');    
[overwrite] = DefaultArgs(varargin,{false},1);
%--------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% RETRIVE session list
sessionList = get_session_list(sessionListName);

% FOREACH Session in Sessions
%     CREATE new session if ses file is absent or overwrite is flagged
for session = sessionList
        
%   CONFIGURE MTA toolbox
    if ~all([strcmp(session.project,    Config.project_name ),...
             strcmp(session.hostServer, Config.host_server  ),...
             strcmp(session.dataServer, Config.data_server  )]),
        MTAstartup(session.project,session.hostServer,session.dataServer);
    end

%   LINK data to project folder
    link_session( session.sessionName, session.dPaths)
        
%   CHECK for existing session and overwrite flag
    if isempty(list_files(session.sessionName,['.',session.mazeName,'.ses.'])) || overwrite,        
%       CREATE Session
        MTASession(session.sessionName,                         ...
                   session.mazeName,                            ...
                   true,                                        ...
                   session.TTLValue,                            ...
                   session.dLoggers,                            ...
                   session.xyzSampleRate);
    end
end

%---------------------------------------------------------------------------------------------------