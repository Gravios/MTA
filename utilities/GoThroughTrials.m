function GoThroughTrials(TrialListName,funcHandle,varargin)
% function GoThroughTrials(TrialListName,funcHandle,varargin)
% 
% run funcHandle over a list of trials
% 
% Note: assumes first argument of funcHandle to accept an MTATrial/MTASession
Trials = SessionList(TrialListName);

for s = 1:numel(Trials)
    Trial = MTATrial.validate(Trials(s));
    try,
        feval(funcHandle,Trial,varargin{:});
    end

end
