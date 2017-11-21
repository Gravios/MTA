function QuickSessionSetup(sessionListName,varargin)
% function QuickSessionSetup(SessionParm,varargin)
% HostConf = load('MTAConf');
% [host,local,link,overwrite] =
% DefaultArgs(varargin,{HostConf.host_server,false},1);
%
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
for sessionArgs = sessionList

        
%   CONFIGURE MTA toolbox
    if ~all([strcmp(sessionArgs.project,    Config.project_name ),...
             strcmp(sessionArgs.hostServer, Config.host_server  ),...
             strcmp(sessionArgs.dataServer, Config.data_server  )]),
        MTAstartup(sessionArgs.project,sessionArgs.hostServer,sessionArgs.dataServer);
    end
        

%   LINK data to project folder
    if all(isfield(sessionArgs,{'xyz_host','nlx_host'}))&&~ispc
        linkSession(sessionArgs.sessionName,...
                    sessionArgs.xyz_host,...
                    sessionArgs.nlx_host);
        dataLoggers = {'nlx','vicon'};
    elseif isfield(sessionArgs,'xyz_host')&&~ispc
        linkSession(sessionArgs.sessionName,...
                    sessionArgs.xyz_host);
        dataLoggers = {'vicon'};
    end

        
%   CHECK for existing session and overwrite flag        
    if isempty(list_files(sessionArgs.sessionName,['.',sessionArgs.mazeName,'.ses.'])) || overwrite,        

%       CHECK for MoCap SampleRate
        if isfield(sessionArgs,'xyzSampleRate'),
            xyzSampleRate = sessionArgs.xyzSampleRate;
        else
            xyzSampleRate = [];
        end
        
%       CREATE Session
        Session = MTASession(sessionArgs.sessionName,...
                             sessionArgs.mazeName,...
                             true,...
                             sessionArgs.TTLValue,...
                             dataLoggers,...
                             'xyzSampleRate',xyzSampleRate);
    end
end

%---------------------------------------------------------------------------------------------------