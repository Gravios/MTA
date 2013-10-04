function  Bccg = gen_state_ccg(Session,varargin)
[state,dur_thresh,partitions,trialName,mazeName] = DefaultArgs(varargin,{'rear',180,10,'all','cof'});

if ~isa(Session,'MTATrial')
    Session = MTATrial(Session,{},trialName,[],[],mazeName);
end

sper = Session.Bhv.getState(state).state;
sname = Session.Bhv.getState(state).label;


sdur = diff(sper,1,2);
sper = sper(sdur>dur_thresh,:);


if partitions==1,

Bccg = MTAccg(Session,[sname 'par_1'] ,['CCG around' sname 'and offset'], ...
             {sper(:,1),sper(:,2)},{1,2},{[sname ' onset'],[sname ' offset']},0);

elseif partitions>1,
%% For partitioning ccgs
Session = Session.load_CluRes();
numClu = size(Session.map,1);

pf_search = MTAPlaceField([]);
pf_search.mazeName = Session.Maze.name;
pf_search.trialName = Session.trialName;
pf_search.trackingMarker = 'head_front';
pf_search.stateLabel = 'head.theta';
pf_search.spk_shuffle = 'n';
pf_search.pos_shuffle = 0;
pf_search.numBSiterations = 1;
pf_search.nbins = 50;
pf_search.smooth = 0.03;
Session = Session.load_Pfs(pf_search);
pfMaxPos = zeros(numClu,2);
for unit = 1:numClu,
    try
        pfMaxPos(unit,:) = Session.Pfs.maxRatePos{unit}(Session.Pfs.maxRateMax{unit},:);
    end
end

Bccg = MTAccg(Session,[sname '_par' num2str(partitions)] ,['CCG around' sname 'and offset'], ...
              {sper(:,1),sper(:,2)},{1,2},{[sname ' onset'],[sname ' offset']},...
              0,'abs_dist',partitions,pfMaxPos,1,[],[]);
end



