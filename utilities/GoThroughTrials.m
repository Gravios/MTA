function GoThroughTrials(TrialListName,funcHandle)

Trials = SessionList(TrialListName);


for s = numel(Trials)
MTAstartup('cin',Trials{s}{4});
Trial = MTATrial(Trials{s}{1},Trials{s}{3},Trials{s}{2});

feval(funcHandle,Trial);

end
