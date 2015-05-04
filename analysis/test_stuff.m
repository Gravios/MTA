
SessionName = 'Ed05-20140529'; 
MazeName = 'ont';
TrialName = 't2'; 


s = MTASession(SessionName,MazeName);

%startStopShift = [20,-1]; 
%QuickTrialSetup(s,'all',startStopShift);

ignoredViconTrials = [1,3:10];
startStopShift = [20,-1]; 
QuickTrialSetup(s,TrialName,startStopShift,ignoredViconTrials,[],false);
Trial = MTATrial(SessionName,TrialName,MazeName);



