function linkSessions(slist,varargin)
[prjName] = DefaultArgs(varargin,{'general'},true);
datPath = getenv('PROJECT')
prjPath_xyz = fullfile(datPath,'data','processed','xyz');
prjPath_nlx = fullfile(datPath,'data','processed','nlx');

slist = SessionList(slist,prjPath_xyz,prjPath_nlx);
for s = slist,
    liskSession(s.sessionName,s.xyz_host,s.nlx_host);
end
