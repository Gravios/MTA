function link_sessions_to_project(sessionList,varargin)
%function link_sessions_to_project(sessionList,varargin)
%
% link data of multiple sessions to project folder
% 
% Assumes the original data you are attempting to link resides
% within your mta data folders
% 
% INPUTS: 
%     sessionList: string,      name of entry within get_session_list
%                  structArray, Structs

% DEFARGS ------------------------------------------------------------------------------------------
[projectName] = DefaultArgs(varargin,{'general'},true);
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

datPath = getenv('MTA_DATA')
projectPath = struct('xyz','','nlx','');
projectPath.xyz = fullfile(datPath,'processed','xyz');
projectPath.nlx = fullfile(datPath,'processed','nlx');


if ischar(slist),
    slist = get_session_list(slist,projectPath.xyz,projectPath.nlx);
end

assert(isstruct(slist),'MTA:utilities:link_sessions_to_project:SessionListNotFound');

for s = slist,
    liskSession(s.sessionName,s.xyz_host,s.nlx_host);
end

%---------------------------------------------------------------------------------------------------