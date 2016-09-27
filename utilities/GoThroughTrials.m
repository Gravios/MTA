function GoThroughTrials(TrialListName,funcHandle,varargin)
% function GoThroughTrials(TrialListName,funcHandle,varargin)
% 
% run funcHandle over a list of trials
% 
% Note: assumes first argument of funcHandle to accept an MTATrial/MTASession
Trials = get_session_list(TrialListName);
for t = Trials
    try,
        feval(funcHandle,MTATrial.validate(t),varargin{:});
    catch err
        for e = err.stack'
            disp(e)
        end
    end
end
