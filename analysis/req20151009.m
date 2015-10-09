function req20151009(Trial)


%% A1 #DESC
% JPDF of the XY distance between the lower spine marker and the
% upper spine marker with contours representing the occupancy of
% other sessions whos labels are derived from the 
%
sltag = 'req20151009';
figTitle = 'BLBU_{pitch} vs d(BMBU_{pitch})/dt';
figPath  = '/storage/gravio/figures/';
hfig = figure(2838293);clf
Trials = SessionList(sltag);
numOfTrials = numel(Trials);
Scolors = jet(numOfTrials);
tlist = {};
state = 'rear';
for s = 1:numOfTrials
    tlist{s} = Trials{s}{1};
    req20151009_A1(MTATrial(Trials{s}{1},Trials{s}{3},Trials{s}{2}),hfig,state,Scolors(s,:));
end
axis xy
caxis([0,200])
legend(tlist,'location','SouthEast')
hfig.Position  = [100   100   782   629];
reportfig(figPath, hfig, 'CJPDF', 'req', false,[sltag '_A1'],figTitle,[],false,'png',[],[],[],false);


%% A2 #DESC
%
%
sltag = 'req20151009';
figTitle = 'BLBU_{pitch} vs d(BMBU_{pitch})/dt';
figPath  = '/storage/gravio/figures/';
hfig = figure(2838293);clf
Trials = SessionList(sltag);
numOfTrials = numel(Trials);
Scolors = jet(numOfTrials);
tlist = {};
state = 'rear';
for s = 1:numOfTrials
    if s == 1, AddBreak = true; else  AddBreak = false; end
    clf(hfig);
    tlist{s} = Trials{s}{1};
    req20151009_A2(MTATrial(Trials{s}{1},Trials{s}{3},Trials{s}{2}),hfig,state);
    axis xy
    caxis([0,200])
    legend({state},'location','SouthEast')
    hfig.Position  = [100   100   782   629];
    reportfig(figPath, hfig, 'CJPDF', 'req', false,[sltag '_A2'],figTitle,[],false,'png',[],[],[],AddBreak);
end



%% B1 #DESC
% JPDF of the XY distance between the lower spine marker and the
% upper spine marker with contours representing the occupancy of
% other sessions whos labels are derived from the 
%
sltag = 'req20151009';
figTitle = 'BLBU_{pitch} vs d(BMBU_{pitch})/dt';
figPath  = '/storage/gravio/figures/';
hfig = figure(2838293);clf
Trials = SessionList(sltag);
numOfTrials = numel(Trials);
Scolors = jet(numOfTrials);
tlist = {};
state = 'rear';
for s = 1:numOfTrials
    tlist{s} = Trials{s}{1};
    req20151009_A1(MTATrial(Trials{s}{1},Trials{s}{3},Trials{s}{2}),hfig,state,Scolors(s,:));
end
caxis([0,200])
legend(tlist,'location','SouthEast')
hfig.Position  = [100   100   782   629];
reportfig(figPath, hfig, 'CJPDF', 'req', false,[sltag '_A1'],figTitle,[],false,'png',[],[],[],false);


%% B2 #DESC
%
%
sltag = 'req20151009';
figTitle = 'BLBU_{pitch} vs d(BMBU_{pitch})/dt';
figPath  = '/storage/gravio/figures/';
hfig = figure(2838293);clf
Trials = SessionList(sltag);
numOfTrials = numel(Trials);
Scolors = jet(numOfTrials);
tlist = {};
state = 'rear';
for s = 1:numOfTrials
    if s == 1, AddBreak = true; else  AddBreak = false; end
    clf(hfig);
    tlist{s} = Trials{s}{1};
    req20151009_A2(MTATrial(Trials{s}{1},Trials{s}{3},Trials{s}{2}),hfig,state);
    caxis([0,200])
    legend({state},'location','SouthEast')
    hfig.Position  = [100   100   782   629];
    reportfig(figPath, hfig, 'CJPDF', 'req', false,[sltag '_A2'],figTitle,[],false,'png',[],[],[],AddBreak);
end



