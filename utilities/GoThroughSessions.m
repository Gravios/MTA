function GoThroughSessions(SessionListName,funcHandle,varargin)

Sessions = SessionList(SessionListName);


for s = 1:numel(Sessions)
MTAstartup([],Sessions{s}{4});
Session = MTASession(Sessions{s}{1},Sessions{s}{2});

feval(funcHandle,Session,varargin{:});

end
