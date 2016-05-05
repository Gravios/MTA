function GoThroughTrials(TrialListName,funcHandle,varargin)
% function GoThroughTrials(TrialListName,funcHandle,varargin)
% 
% run funcHandle over a list of trials
% 
% Note: assumes first argument of funcHandle to accept an MTATrial/MTASession
Trials = SessionList(TrialListName);
for t = Trials
    try,
        feval(funcHandle,MTATrial.validate(Trials(t)),varargin{:});
    catch err
        for e = err.stack'
            disp(e)
        end
    end
end
