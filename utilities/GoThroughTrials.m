function GoThroughTrials(TrialListName,funcHandle,varargin)

Trials = SessionList(TrialListName);


for s = 1:numel(Trials)
    %MTAstartup('lmu',Trials{s}{4});
    Trial = MTATrial(Trials{s}{1},Trials{s}{3},Trials{s}{2});

    try,
        feval(funcHandle,Trial,varargin{:});
    end

end
