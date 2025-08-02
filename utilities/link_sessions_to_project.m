function link_sessions_to_project(sessionlist,varargin)
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

if ischar(sessionlist),
    sessionlist = query_session_list(session);
end


for ses = sessionlist,
    link_session_Dpath(ses.sessionName, ses.dPaths);
end

%---------------------------------------------------------------------------------------------------
