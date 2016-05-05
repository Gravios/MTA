function GoThroughSessions(SessionListName,funcHandle,varargin)
% function GoThroughSessions(SessionListName,funcHandle,varargin)
% 
% run funcHandle over a list of trials
% 
% Note: assumes first argument of funcHandle to accept an MTASession
Sessions = SessionList(SessionListName);
for s = Sessions
    try,
        feval(funcHandle,MTASession.validate(s),varargin{:});
    catch err
        for e = err.stack'
            disp(e)
        end
    end

end
