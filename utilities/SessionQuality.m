function s_quality = SessionQuality(Session)

%% check if there are 8 or more markers
if Session.Model.N<8,return,end

s_quality.sessionName = Session.filebase;

%% Number of Overall units
% number of spikes must be greater than 10
[Res Clu Map] = LoadCluRes([Session.spath.nlx Session.name]);
Res = round(Res*Session.lfpSampleRate/Session.sampleRate);
[~,ind] = SelectPeriods(Res,Session.syncPeriods,'d',1,1);
myClu = Clu(ind);

spknum = 0;
for i = 1:length(unique(myClu))
    if length(myClu(myClu==i))>10,
        spknum = spknum+1;
    end
end
s_quality.totalGCluCnt = spknum;


%% Total Distance Traveled
dxyz = diff(sq(Session.xyz(:,Session.Model.gmi(Session.trackingMarker),1:2)),2);
xySpeed = sqrt(sum(dxyz.^2,2));
s_quality.totalDistance = sum(xySpeed)/1000;    


%% State dependent unit count
% number of spikes must be greater than 10
num_bhv_states = length(Session.Bhv.States);

% Backup bhv xml if it exists
if exist([Session.spath.analysis Session.filebase '.bhv.auto.xml'],'file'),
    system(['mv ' Session.spath.analysis Session.filebase '.bhv.auto.xml ' Session.spath.analysis Session.filebase '.bhv.auto.bkp'])
end
system(['cp ' Session.path.root 'config/xml/bhv.auto.xml ' Session.spath.analysis Session.filebase '.bhv.auto.xml']);
Session = Session.autoLabelBhv('auto');
for i = 1:length(Session.Bhv.States)
    [myRes,ind] = SelectPeriods(Res,stsp,'d',1,1);
    myClu = Clu(ind);
    spknum = 0;
    for j = 1:length(unique(myClu))
        if length(myClu(myClu==j))>10,
            spknum = spknum+1;
        end
    end
    s_quality.stateLabels{i} = Session.Bhv.States{i}.label;
    s_quality.stateGCluCnt{i} = spknum;
end



%% Vissualization

%% xy occupancy


%% rearing
rear = Session.Bhv.getState('rear')
scatter(Session.xyz(rear.state(:,1),Session.Model.gmi('head_front'),1),Session.xyz(rear.state(:,1),Session.Model.gmi('head_front'),2))
xlim([-500,500])
ylim([-500,500])



Sessions =  {};
Session_Summary.name = Session.filebase


