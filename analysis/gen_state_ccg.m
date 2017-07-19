function  Bccg = gen_state_ccg(Trial,varargin)
Trial = MTATrial.validate(Trial);

[state,durationThreshold,partitions] = DefaultArgs(varargin,{'rear',1800,1});


sper = [ Trial.Stc{state} ];

sper = [ mStc{state,Trial.lfp.sampleRate} ];


sdur = diff(sper.data,1,2);
sper.data = sper(sdur>durationThreshold,:);


if partitions==1,
    Bccg = MTAccg(Trial,...                                   % Trial
                  [sper.label 'par_1'] ,...                   % Name
                  ['CCG around' sper.label 'and offset'], ... % Description
                  {sper(:,1),sper(:,2)},...                   % ResTrains
                  {[sper.label ' onset'],[sper.label ' offset']},... % CluTags
                  [],true);

elseif partitions>1,
    %% For partitioning ccgs
    Trial = Trial.load_CluRes();
    numClu = size(Trial.map,1);

    pf_search = MTAPlaceField([]);
    pf_search.mazeName = Trial.maze.name;
    pf_search.trialName = Trial.name;
    pf_search.trackingMarker = 'head_front';
    pf_search.stateLabel = 'head.theta';
    pf_search.spk_shuffle = 'n';
    pf_search.pos_shuffle = 0;
    pf_search.numBSiterations = 1;
    pf_search.nbins = 50;
    pf_search.smooth = 0.03;
    Trial = Trial.load_Pfs(pf_search);
    pfMaxPos = zeros(numClu,2);
    for unit = 1:numClu,
        try
            pfMaxPos(unit,:) = Trial.Pfs.maxRatePos{unit}(Trial.Pfs.maxRateMax{unit},:);
        end
    end

    Bccg = MTAccg(Trial,[sper.label '_par' num2str(partitions)] ,...
                  ['CCG around' sper.label 'and offset'], ...
                  {sper(:,1),sper(:,2)},{1,2},...
                  {[sper.label ' onset'],[sper.label ' offset']},...
                  0,'abs_dist',partitions,pfMaxPos,1,[],[]);
end



