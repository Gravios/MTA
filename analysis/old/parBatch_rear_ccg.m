function parBatch_rear_ccg(Session,varargin)
[trialName,batch_mode] = DefaultArgs(varargin,{'all','run'});

%% Load Session 
Session = MTASession(Session,{'CluRes'});
if ~isempty(trialName)
    Session = MTATrial(Session,{},trialName);
end


%% Batch Information
time_stamp = '20121231T1730';

%name - string: unique descriptive name of the ccg
name = 'test_par1';

%Description - string: short statement describing the ccg
Description = 'testing';

% $$$ ['CCG of rearing onset/offset and non-rearing head periods'];

% $$$  partitioned into ' ...
% $$$                'groups based on their distance from the center of ' ...
% $$$                'each placefield'];

method = 'abs_dist';
surrogate_sample_size = [];
partitions = 1;
numIterations = 1;%1000;



%% Check if Backups Exist
switch batch_mode

  case 'setup'
    %% Save a script backup 
    oriScript = mfilename('fullpath');
    system(['cp ' oriScript '.m ' Session.path.root 'batch/scripts/' mfilename '-' time_stamp '.m']);

  case 'run'
    %% MTAccg Parameters
    numClu = size(Session.map,1);

    rper =   Session.Bhv.getState('rear').state;
    nrhper = Session.Bhv.getState('nrhp').state;
    rdur = diff(rper,1,2)/Session.xyzSampleRate;
    rpos = sq(Session.xyz(rper(:,1),Session.Model.gmi(Session.trackingMarker),[1,2]));
    rper = rper(rdur>1.5,:);
    rpos = rpos(rdur>1.5,:);

    Res = {rper(:,1),rper(:,2),nrhper};
    Clu = {1,2,3};
    CluTags = {'rear onset','rear offset','non-rearing head periods'};

    %% Create partition feature vectors if there is more that one partition
    if partitions >1;
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
    else 
        pfMaxPos = [];
    end

    %% Run MTACCG
    Bccg = MTAccg(Session,name,Description,Res,Clu,CluTags,0,method,partitions,pfMaxPos,1,surrogate_sample_size,numIterations);

end



