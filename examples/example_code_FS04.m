%% How to MTA your Data

% The Following topics will be discussed
% 
%    1. Configuration:      Setting up the MTA environment on linux
%    2. Data Arangment:     Where and how to store your data
%    3. Meta Data:          Meta data for generating a session object
%    4. Preprocess CSV:     parse_rbo_from_csv.m
%    5. Build Session:      generate the MTASession objects
%    6. Load Session:       Various methods for loading an MTASession object
%    7. Load subjects:      Load arena and subject rbo(s) 
%    8. Build Trial:        generate the MTATrial objects
%    9. Label Behaviors:    Label basic behaviors (theta periods, good periods, speed thresholded periods)
%   10. Compute ratemaps:   Do it.
%   11. Simple figure:      booo.
%   12. Fabulous figure:    wooo.
%   13. Sort Units:         by pyramidal and interneurons, and by center and border types
%   14. Decode position:    from placefields and unit firing rate
%   15. Decode position:    by theta phase

%%%<<< 1. Configuration ----------------------------------------------------------------------------
% NOTICE - if your system has never been configured for MTA follow the next step
% CONFIGURE your OS environmental variables ->  saved in ~/.bashrc 
%
% shell $> cd /storage/share/matlab/MTA
% shell $> ./configure -u fabian -s storage2 -d Data -t /storage2/fabian/Code 

%%%>>>


%%%<<< 2. Data Arangement --------------------------------------------------------------------------
% example for FS04/FS04-20210323a
%
% CREATE data directories 
% shell $> mkdir -p /storage2/fabian/Data/raw/xyz/FS04/FS04-20210323/vrr
% shell $> mkdir -p /storage2/fabian/Data/processed/xyz/FS04/FS04-20210323/vrr
% shell $> mkdir -p /storage2/fabian/Data/processed/ephys/FS04/FS04-20210323
%
% MOVE data
% shell $> mv /storage2/fabian/Data/raw/FS04/HSW_2021_03_23__20_07_03__30min_44sec__hsamp_64ch_25000sps /storage2/fabian/Data/processed/ephys/FS04/
% OR 
% SYMLINK data
% shell $> cp -Rs /storage2/fabian/Data/raw/FS04/HSW_2021_03_23__20_07_03__30min_44sec__hsamp_64ch_25000sps /storage2/fabian/Data/processed/ephys/FS04/
%
% RENAME from HSW_ to nlx format
% shell $> cd /storage2/fabian/Data/processed/FS04/
% shell $> hsw_rename_files 
%
% RENAME data files
% MOVE motion tracking csv files to Data/raw/xyz
% shell $> cp /storage2/fabian/Data/processed/ephys/FS04/FS04-20210322a/Take*.csv /storage2/fabian/Data/raw/xyz/FS04/FS04-20210323a/vrr
%
%%%>>>


%%%<<< 3. Meta Data --------------------------------------------------------------------------------
% SETUP struct which holds session metadata
% struct Session  
%
% WARNING : sessionName - Must match data directory names 
% WARNING : mazeName    - Explicitly register in MTA/MTAConfiguration.m
% WARNING : trialName   - 'all' is default where entire session is considered a trial
% WARNING : dLoggers    - see: MTA/utilities/synchronization/format.txt
% WARNING : thetaRef    - Only one channel per electrode shank
% WARNING : thetaRefGeneral - Only one channel
%
%--------------------------------------------------------------
%|  sessionName    |  String        |   'jg05-20120309'       | 
%|  subjects       |  CellStr       |   {{'jg05'}}            | 
%|  mazeName       |  String        |   'cof'                 | 
%|  trialName      |  String        |   'all'                 | 
%|  dLoggers       |  CellStr       |   {{'nlx','vicon'}}     |
%|  dPaths         |  struct        |   SEE dPaths bellow     |
%|  xyzSampleRate  |  Double        |   119.881035            |
%|  hostServer     |  String        |   'lmu'                 |
%|  dataServer     |  String        |   'lmu',                |
%|  project        |  String        |   'general'             |
%|  TTLValue       |  CONTEXT       |   SEE TTLValue bellow   |
%|  includeSyncInd |  Array[double] |   []                    |
%|  offsets        |  Array[double] |   [15,-15]              |
%|  xOffSet        |  Double        |   0                     |
%|  yOffSet        |  Double        |   0                     |
%|  stcMode        |  String        |   'default'             |
%|  rotation       |  Double        |   0                     |
%|  thetaRef       |  Array[double] |   [8:8:64]              | 
%| thetaRefGeneral |  Double        |   8                     | 
%--------------------------------------------------------------
%
% dPaths structure:
%----------------------------------------------------------------------------------------------
%|     FIELD       |      TYPE      |          example                                        |
%----------------------------------------------------------------------------------------------
%|  xyz            |  String        |   '/storage/<user>/data/processed/xyz/<subjectId>/'     |
%|  ephys          |  String        |   '/storage/<user>/data/processed/nlx/<subjectId>/'     |
%|  video          |  String        |   '/storage/<user>/data/processed/video/<subjectId>/'   |
%----------------------------------------------------------------------------------------------

% METADATA  ----------------------------------------------------------------------------------
% SET Session info
meta.sessionBase = 'FS04';
meta.sessionName = 'FS04-20210322a';
meta.mazeName = 'vrr';
meta.trialName = 'all';
meta.dLoggers = {'WHT','CSV'};
meta.dPaths.xyz   = fullfile('/storage/gravio/data/processed/xyz/',  meta.sessionBase);
meta.dPaths.ephys = fullfile('/storage/gravio/data/processed/ephys/',meta.sessionBase);
meta.xyzSampleRate = 120.00;% ???
meta.hostServer = 'lmu';
meta.dataServer = 'lmu';
meta.project    = 'general';
% SET subject info
meta.primarySubject = 'FS04';
meta.subjects(1).name = 'FS04';
meta.subjects(1).type = 'rbo';
meta.subjects(1).rb(1).name = 'Head';
meta.subjects(1).rb(1).alias = 'RatW';
meta.subjects(2).name = 'Arena';
meta.subjects(2).type = 'rbo';
meta.subjects(2).rb(1).name = 'Arena';
meta.subjects(2).rb(1).alias = 'Arena';
% SET Trials
meta.includeSyncInd = [];
meta.offsets  = [0,0];
% SET Maze corrections if the maze isn't already centered
meta.xOffSet  = 0;
meta.yOffSet  = 0;
meta.zOffSet  = 270;
meta.rotation = 0;
% SET State label collection
meta.stcMode  = 'default';
% SET LFP Theta reference channels
meta.thetaRef = [1:11:64];
meta.thetaRefGeneral = 1;


% CONSTRUCT paths
MTA_DATA_PATH = getenv('MTA_DATA_PATH');
meta.path.raw.ephys =                                                        ...
    create_directory(                                                        ...
        fullfile(MTA_DATA_PATH,'raw/ephys/',meta.sessionBase,meta.sessionName));

meta.path.raw.xyz =                                                                     ...
    create_directory(                                                                   ...
        fullfile(MTA_DATA_PATH,'raw/xyz/',meta.sessionBase,meta.sessionName,meta.mazeName));

meta.path.processed.xyz =                                                                       ...
    create_directory(                                                                           ...
        fullfile(MTA_DATA_PATH,'processed/xyz/',  meta.sessionBase,meta.sessionName,meta.mazeName));

meta.path.processed.ephys =                                                       ...
    create_directory(                                                             ...
        fullfile(MTA_DATA_PATH,'processed/ephys/',meta.sessionBase,meta.sessionName));


% Constructed as cell array with empty strings where no mocap occured.
% WARNING 1 take per record
% EXAMPLE : given 3 records exist in a session
%           and record 1 and 3 have associated csv files
%           the meta.csv field should look like this 
%
% meta.csv = {'Take 2021-03-22 08.24.46 PM.csv',     '', 'Take 2021-03-22 09.24.46 PM.csv'};
%                      record1,                 record2,              record3
%
meta.csv = {'Take 2021-03-22 08.24.46 PM.csv'};

% MOVE csv file from the raw ephys folder to raw xyz folder
cf(@(csv)                                                               ...
   movefile(fullfile(meta.path.processed.ephys,csv),meta.path.raw.xyz), ...
   meta.csv(~isempty(meta.csv)));


% NOTE : not really a TTLValue
%
meta.TTLValue = convert_srs_to_TTLValue(meta);    
% EXPECTED output
% meta.TTLValue = { ...                                record1
%                   {[46029977]},...                     record sample count
%                   {'FS04-20210322a.take_0001.mat'},... file with extracted rigidbodies
%                   {'FS04'}...                          session's primary subject
% }

%%%>>>


%%%<<< 4. Preprocess CSV ---------------------------------------------------------------------------
%
% EXTRACTS rigid bodies from the csv file generated by motive
%
parse_rbo_from_csv(meta);

%%%>>>


%%%<<< 5. Build Session ----------------------------------------------------------------------------
link_session( meta.sessionName, meta.dPaths)
build_sessions(meta);
% -> generates symbolic link file structure in project folder
% -> synchronizes lfp and rbo/xyz objecs
% -> creates default objects

%%%>>>


%%%<<< 6. Load Session -----------------------------------------------------------------------------
% NOTE the next 3 commands do the same thing
% Meta data struct
Session = MTASession.validate(meta);
% Piecing together the session name 
Session = MTASession.validate([meta.sessionName,'.',meta.mazeName,'.',meta.trialName]);
% writing out the session name 
Session = MTASession.validate('FS04-20210322a.vrr.all');

% $$$ Session.load('spk');
% $$$ Session.save();

%%%>>>


%%%<<< 7. Load Subjects ----------------------------------------------------------------------------
% problem: getting the modified subject objects passed to the pfs computation
% problem: loading Objects

rat   = Session.load('subject','FS04');
Arena = Session.load('subject','Arena');





% FILL gaps within recording periods
oriRat = copy(rat);

gapPer = ThreshCross(isnan(rat(:,1,1)),0.5,0);
gapLen = diff(gapPer,1,2);
gapPer = bsxfun(@plus,gapPer,[1,0]);
gapDstBA = [gapPer(:,1);size(rat,1)]-[1;gapPer(:,2)];

rat.fill_gaps(4,30);

rat.fill_gaps(20,50);

rat.fill_gaps(50,100);

rat.fill_gaps(300,inf);

figure();
hold('on');
plot(rat(:,1,3),'.r')
plot(oriRat(:,1,3),'.c')
plot(rat(:,1,2),'.r')
plot(oriRat(:,1,2),'.c')
plot(rat(:,1,1),'.r')
plot(oriRat(:,1,1),'.c')
Lines(gapPer(:),[],'k');


figure,
hold('on');
plot(rat(:,1,1),rat(:,1,2),'.r')
plot(oriRat(:,1,1),oriRat(:,1,2),'.c')


% $$$ figure();
% $$$ plot(Arena(:,'Arena',1),Arena(:,'Arena',2),'.')
% $$$ hold('on');
% $$$ plot(mean(Arena(:,'Arena',1),'omitnan'),mean(Arena(:,'Arena',2),'omitnan'),'or')
% $$$ plot(mean([min(Arena(nniz(Arena(:,1,2)),'Arena',1)),max(Arena(nniz(Arena(:,1,2)),'Arena',1))]),...
% $$$      mean([min(Arena(nniz(Arena(:,1,2)),'Arena',2)),max(Arena(nniz(Arena(:,1,2)),'Arena',2))]),...
% $$$      'dg');
% $$$ plot([min(Arena(nniz(Arena(:,1,2)),'Arena',1)),max(Arena(nniz(Arena(:,1,2)),'Arena',1))],...
% $$$      [min(Arena(nniz(Arena(:,1,2)),'Arena',2)),max(Arena(nniz(Arena(:,1,2)),'Arena',2))],...
% $$$      'dm');


mazeVec = [max(Arena(nniz(Arena(:,1,2)),'Arena',1))-min(Arena(nniz(Arena(:,1,2)),'Arena',1)),...
           max(Arena(nniz(Arena(:,1,2)),'Arena',2))-min(Arena(nniz(Arena(:,1,2)),'Arena',2))];
mazeVec = mazeVec./norm(mazeVec);


% $$$ o = atan2(mazeVec(1),mazeVec(2));
% $$$ nArena = multiprod([cos(o),-sin(o);sin(o),cos(o)],sq(Arena(:,1,1:2)),[1,2],[2]);
% $$$ figure();
% $$$ hold('on');
% $$$ plot(Arena(:,'Arena',1),Arena(:,'Arena',2),'.')
% $$$ plot(nArena(:,1),nArena(:,2),'.r')


% $$$ figure,
% $$$ plot(nratAC(:,1,1),nratAC(:,1,2),'.');
% $$$ figure,
% $$$ plot(rat(:,1,1),rat(:,1,2),'.');
% $$$ figure();
% $$$ plot(Arena(:,'Arena',1))
% $$$ figure();
% $$$ plot(Arena(:,'Arena',2))
% $$$ figure();
% $$$ plot(Arena(:,'Arena',3))
% $$$ figure();
% $$$ plot(ratAC(:,'Head',1),ratAC(:,'Head',2),'.');
% $$$ 
% $$$ xyzShift = [min(Arena(nniz(Arena(:,1,1)),'Arena',1)),min(Arena(nniz(Arena(:,1,2)),'Arena',2)),median(Arena(nniz(Arena(:,1,2)),'Arena',3))];
% $$$ 
% $$$ xyzShift = [min(Arena(nniz(Arena(:,1,1)),'Arena',1)),min(Arena(nniz(Arena(:,1,2)),'Arena',2)),median(Arena(nniz(Arena(:,1,2)),'Arena',3))];
% $$$ xyzShift = [mean([min(Arena(nniz(Arena(:,1,1)),'Arena',1)),max(Arena(nniz(Arena(:,1,1)),'Arena',1))]),...
% $$$             mean([min(Arena(nniz(Arena(:,1,2)),'Arena',2)),max(Arena(nniz(Arena(:,1,2)),'Arena',2))]),...
% $$$             median(Arena(nniz(Arena(:,1,2)),'Arena',3))];

xyzShift = [mean([min(Arena(nniz(Arena(:,1,1)),'Arena',1)),max(Arena(nniz(Arena(:,1,1)),'Arena',1))]),...
            mean([min(Arena(nniz(Arena(:,1,2)),'Arena',2)),max(Arena(nniz(Arena(:,1,2)),'Arena',2))]),...
            median(Arena(nniz(Arena(:,1,2)),'Arena',3))];


% CONVERT the xyz coordinates to Room frame of refrences centered on mean Arena position
ratRC = copy(rat);
% TRANSLATE 
ratRC.data(:,1,1:3) = bsxfun(@minus,ratRC.data(:,1,1:3),permute(xyzShift,[1,3,2]));
% ROTATE 
o = 0.0349065850398866;
ratRC.data(:,1,1:2) = multiprod([cos(o),-sin(o);sin(o),cos(o)],sq(ratRC.data(:,1,1:2)),[1,2],[2]);
ratRC.label = [ratRC.name,'_RC'];
ratRC.update_filename([Session.filebase,'.rbo.',rat.name,'_RC.s.mat']);
ratRC.save();


% CONVERT the xyz coordinates to Arena frame of refrences centered on the Arena
ratAC = copy(rat);
% TRANSLATE 
ratAC.data(:,1,1:3) = ratAC(:,'Head',1:3) - [Arena(:,'Arena',1:3);zeros([1,1,3])];
% ROTATE 
o = 0.0349065850398866;
ratAC.data(:,1,1:2) = multiprod([cos(o),-sin(o);sin(o),cos(o)],sq(ratAC.data(:,1,1:2)),[1,2],[2]);
ratAC.data(:,1,3) = ratAC.data(:,1,3)+meta.zOffSet;
ratAC.label = [ratAC.name,'_AC'];
ratAC.update_filename([Session.filebase,'.rbo.',rat.name,'_AC.s.mat']);
ratAC.save();

% NOTE arena rotation does not correct the rat head quaternion :(

%%%>>>


%%%<<< 8. Build Trial ------------------------------------------------------------------------------
QuickTrialSetup(meta);

%%%>>>


%%%<<< 9. Load Trials ------------------------------------------------------------------------------
% NOTE the next 3 commands do the same thing
% Meta data struct
Trial = MTATrial.validate(meta);
% Piecing together the session name 
Trial = MTATrial.validate([meta.sessionName,'.',meta.mazeName,'.',meta.trialName]);
% writing out the session name 
Trial = MTATrial.validate('FS04-20210322a.vrr.all');


%%%>>>


%%%<<< 9. Label Behaviors --------------------------------------------------------------------------
% LABEL 



rat   = Trial.load('subject','FS04_AC');
ratAC   = Trial.load('subject','FS04_AC');
Arena = Trial.load('subject','Arena');


% LABEL theta periods from lfp
label_theta(Session);

% LABEL Non nan periods
gper = ThreshCross(double(nniz(rat(:,1,1))),0.5,10);
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   gper,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'gper','a');
Trial.stc.states{end}.save(1);


% LABEL speed thresholded periods
ratFilt = filter(copy(ratAC),'ButFilter',4,1,'low');
%figure();plot(sqrt(sum(diff(ratFilt(:,1,1:2)).^2,3)).*ratFilt.sampleRate./10)
vper = ThreshCross(double(sqrt(sum(diff(ratFilt(:,1,1:2)).^2,3)).*ratFilt.sampleRate./10>2),0.5,10);
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   vper,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'vel','v');
Trial.stc.states{end}.save(1);



% STATE ArenaConfigA
xyzShift = [mean([min(Arena(nniz(Arena(:,1,1)),'Arena',1)),max(Arena(nniz(Arena(:,1,1)),'Arena',1))]),...
            mean([min(Arena(nniz(Arena(:,1,2)),'Arena',2)),max(Arena(nniz(Arena(:,1,2)),'Arena',2))]),...
            median(Arena(nniz(Arena(:,1,2)),'Arena',3))];
arenaConfA = ThreshCross(double(Arena(:,'Arena',2)>974),0.5,10);
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   arenaConfA,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'ArenaConfA','A');
Trial.stc.states{end}.save(1);
% STATE ArenaConfigB
arenaConfB = ThreshCross(double(Arena(:,'Arena',2)<700),0.5,10);
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   arenaConfB,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'ArenaConfB','B');
Trial.stc.states{end}.save(1);
% STATE ArenaConfigC
arenaConfC = ThreshCross(double(Arena(:,'Arena',2)<974&Arena(:,'Arena',2)>700),0.5,10);
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   arenaConfC,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'ArenaConfC','C');
Trial.stc.states{end}.save(1);

% STATE rear
rfet = filter(copy(ratAC),'ButFilter',4,0.5,'low');
rearPer = ThreshCross(double(rfet(:,1,3)>180),0.5,10);
rearPer(end-3:end,:) = [];
figure();
hold('on');
plot(ratAC(:,1,3),'.');
plot(rfet(:,1,3),'.r');
Lines(rearPer(:),[],'k');
Lines(rearPer(:,1)-round(ratAC.sampleRate*0.5),[],'g');
Lines(rearPer(:,2)+round(ratAC.sampleRate*0.5),[],'m');

rearPer(:,1) = rearPer(:,1) - round(ratAC.sampleRate*0.5);
rearPer(:,2) = rearPer(:,2) + round(ratAC.sampleRate*0.5);

Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   rearPer,...
                   ratAC.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'rear','r');
Trial.stc.states{end}.save(1);

figure,
rper = Trial.stc{'r'};
plot(ratAC(rper,1,1),ratAC(rper,1,2),'.');

%%%>>>


%%%<<< 10. Compute ratemaps ------------------------------------------------------------------------
% SUBSET of xyz within active stae
% p(s)
% 
%
activeState = 'theta&gper&vel'; 


% ROOM frame of reference  (X,Y)
pfsArgs = struct('states',           activeState,                  ... Computational Periods 
                 'binDims',          [30,30],                      ... Physical size of bins in milimeters
                 'SmoothingWeights', [2.5,2.5],                    ... Gaussian smoother prameters, std deviation in bins 
                 'numIter',          1,                            ... number of bootstraps
                 'boundaryLimits',   [-500,500;-1000,1000],        ... Computational domain
                 'halfsample',       false);                      %... throw out half of the data each iteration ... or don't

pfs = compute_ratemaps(Trial,                                      ... MTATrial
                       [],                                         ... Unit list (e.g. [1, 2, ... , N])
                       @fet_rbo_RC,                                ... Function handle, ratemap space
                       [],                                         ... sampleRate, default 16Hz to speed it up
                       pfsArgs,                                    ... Arguments to the MTAApfs 
                       'overwrite',true);

                       
% ARENA frame of reference (X,Y)
pfsArgs = struct('states',           activeState,                  ...
                 'binDims',          [30,30],                      ...
                 'SmoothingWeights', [2.5,2.5],                    ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-500,500;-1000,1000],        ...
                 'halfsample',       false);
pfa = compute_ratemaps(Trial,[],@fet_rbo_AC,[],pfsArgs,'overwrite',true);



% ARENA frame of reference (X,Y)
pfsArgs = struct('states',           [activeState, '&ArenaConfA'], ...
                 'binDims',          [30,30],                      ...
                 'SmoothingWeights', [2.5,2.5],                    ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-500,500;-1000,1000],        ...
                 'halfsample',       false);
pfs_A = compute_ratemaps(Trial,[],@fet_rbo_AC,[],pfsArgs,'overwrite',true);

pfsArgs = struct('states',           [activeState, '&ArenaConfB'], ...
                 'binDims',          [30,30],                      ...
                 'SmoothingWeights', [2.5,2.5],                    ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-500,500;-1000,1000],        ...
                 'halfsample',       false);
pfs_B = compute_ratemaps(Trial,[],@fet_rbo_AC,[],pfsArgs,'overwrite',true);

pfsArgs = struct('states',           [activeState, '&ArenaConfC'], ...
                 'binDims',          [30,30],                      ...
                 'SmoothingWeights', [2.5,2.5],                    ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-500,500;-1000,1000],        ...
                 'halfsample',       false);
pfs_C = compute_ratemaps(Trial,[],@fet_rbo_AC,[],pfsArgs,'overwrite',true);


% ARENA frame of reference (Z,Y)
pfsArgs = struct('states',           activeState,                  ...
                 'binDims',          [30,30],                      ...
                 'SmoothingWeights', [2.5,2.5],                    ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-1000,1000;-300,0],        ...
                 'halfsample',       false);
pfz = compute_ratemaps(Trial,[],@fet_rbo_YZ,[],pfsArgs,'overwrite',true);

pfsArgs = struct('states',           [activeState, '&ArenaConfA'], ...
                 'binDims',          [30,30],                      ...
                 'SmoothingWeights', [2.5,2.5],                    ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-1000,1000;-300,0],        ...
                 'halfsample',       false);
pfz_A = compute_ratemaps(Trial,[],@fet_rbo_YZ,[],pfsArgs,'overwrite',true);

pfsArgs = struct('states',           [activeState, '&ArenaConfB'], ...
                 'binDims',          [30,30],                      ...
                 'SmoothingWeights', [2.5,2.5],                    ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-1000,1000;-300,0],        ...
                 'halfsample',       false);
pfz_B = compute_ratemaps(Trial,[],@fet_rbo_YZ,[],pfsArgs,'overwrite',true);

pfsArgs = struct('states',           [activeState, '&ArenaConfC'], ...
                 'binDims',          [30,30],                      ...
                 'SmoothingWeights', [2.5,2.5],                    ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-1000,1000;-300,0],        ...
                 'halfsample',       false);
pfz_C = compute_ratemaps(Trial,[],@fet_rbo_YZ,[],pfsArgs,'overwrite',true);


%%%>>>


%%%<<< 11. Simple figure ---------------------------------------------------------------------------
% START SIMPLE FIGURE ------------------------------------------------------------------------------
% PLOT Placefields for each Unit
hfig = figure();
u = 1;
while u~=-1
    subplot(231);
        plot(pfs,u,1,'text',[],false);
        daspect([1,1,1]);
        title('room coord');
    subplot(232);
        plot(pfa,u,1,'text',[],false);    
        daspect([1,1,1]);
        title('arena coord');        
    subplot(233);
        plot(pfz,u,1,'text',[],false,'flipAxesFlag',true);    
        daspect([1,1,1]);
        xlabel('Low       High');
        title(num2str(u));
        
    subplot(234);
        plot(pfs_A,u,1,'text',[],false);
        daspect([1,1,1]);
        title(pfs_A.parameters.states);
    subplot(235);
        plot(pfs_B,u,1,'text',[],false);
        daspect([1,1,1]);
        title(pfs_B.parameters.states);
    subplot(236);
        hold('on');
        plot(pfs_C,u,1,'text',[],false);
        daspect([1,1,1]);
        title(pfs_C.parameters.states);
        plot(ratAC([Trial.stc{[activeState, '&ArenaConfC']}],1,1),...
             ratAC([Trial.stc{[activeState, '&ArenaConfC']}],1,2),...
             '.','MarkerSize',1);
    u = figure_controls(hfig,u,Trial.spk.map(:,1));
end
% END SIMPLE FIGURE --------------------------------------------------------------------------------

%%%>>>


%%%<<< 12. Fabulous figure -------------------------------------------------------------------------
% START fabulous figure  ---------------------------------------------

figpath = fullfile('/storage/gravio/figures/analysis/fabian/',Trial.filebase);
create_directory(figpath);
configAColor = 'm';
configBColor = 'g';

for unit = 1:74;

[hfig,fig,fax,sax] = set_figure_layout(figure(666001),'A4','landscape',[],10,3,0.1,0.2);
delete(sax);
sax = gobjects([0,1]);

text(fax,                ... axesHandle
     fig.page.marginLeft,... xpos
     fig.page.height-3,  ... ypos
     {'Vestibular',...
      '30cm shift',...
      Trial.filebase,...
      [' Unit: ',num2str(unit)]},...
     'FontSize',24);

ratemapXY   = plot(pfs,  unit,[],[],[],false);
ratemapXYA  = plot(pfs_A,unit,[],[],[],false);
ratemapXYB  = plot(pfs_B,unit,[],[],[],false);
clim = [0,max([ratemapXY(:);ratemapXYA(:);ratemapXYB(:)])];
% ROOM COORDINATES ---------------------------------------------------------------------------------
secXind = 1;
secYind = 2;

[yind, yOffSet, xind, xOffSet] = deal(secYind, 0, secXind, 1.5);
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);


binsXY      = cell([1,2]);
[binsXY{:}] = meshgrid(pfs.adata.bins{:});
binsXY      = cat(2,binsXY,{zeros(size(binsXY{1}))});

hold('on');
set( surf(binsXY{:},ratemapXY'), 'EdgeColor', 'none');
view([90,-90]);
daspect([1,1,1]);
caxis(clim);
cax = colorbar(sax(end));
cax.Units = 'centimeters';
drawnow();
pause(0.1);
cax.Position(1) = sum(sax(end).Position([1,3]))+0.5;
colormap('parula');

sax(end).XTick = [];
sax(end).YTick = [];
sax(end).Color = [0.9,0.9,0.9];
axis('off')

rax = rectangle('Position',[-510,-1000,1000,1700],'LineWidth',2,'EdgeColor',configAColor);
rax = rectangle('Position',[-490,-700,1000,1700],'LineWidth',2,'EdgeColor',configBColor);


% MAZE COORDINATES CONF A & B XY -----------------------------------------------------------------------
yOffSet = 2;
[yind, yOffSet, xind, xOffSet] = deal(secYind+2, yOffSet, secXind, 1);
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);


binsXY      = cell([1,2]);
[binsXY{:}] = meshgrid(pfs.adata.bins{:});
binsXY      = cat(2,binsXY,{zeros(size(binsXY{1}))});

hold('on');
set( surf(binsXY{:},ratemapXYA'), 'EdgeColor', 'none');
view([90,-90]);
daspect([1,1,1]);
caxis(clim);

sax(end).XTick = [];
sax(end).YTick = [];
sax(end).Color = [0.9,0.9,0.9];
axis('off')

rax = rectangle('Position',[-500,-850,1000,1700],'LineWidth',2,'EdgeColor',configAColor);


[yind, yOffSet, xind, xOffSet] = deal(secYind+3, yOffSet, secXind, 2);
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);

binsXY      = cell([1,2]);
[binsXY{:}] = meshgrid(pfs.adata.bins{:});
binsXY      = cat(2,binsXY,{zeros(size(binsXY{1}))});

hold('on');
set( surf(binsXY{:},ratemapXYB'), 'EdgeColor', 'none');
view([90,-90]);
daspect([1,1,1]);
caxis(clim);

sax(end).XTick = [];
sax(end).YTick = [];
sax(end).Color = [0.9,0.9,0.9];
axis('off')

rax = rectangle('Position',[-500,-850,1000,1700],'LineWidth',2,'EdgeColor',configBColor);

% YZ STUFF
ratemapYZ = plot(pfz,unit,[],[],[],false);
ratemapYZA = plot(pfz_A,unit,[],[],[],false);
ratemapYZB = plot(pfz_B,unit,[],[],[],false);
clim = [0,max([ratemapYZ(:);ratemapYZA(:);ratemapYZB(:)])];
% ROOM COORDINATES YZ -------------------------------------------------------------------------------
secXind = 2;
[yind, yOffSet, xind, xOffSet] = deal(secYind, 0, secXind, 1.5);
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);


binsYZ      = cell([1,2]);
[binsYZ{:}] = meshgrid(pfz.adata.bins{:});
binsYZ{2} = binsYZ{2}+300;
binsYZ      = cat(2,{pfs.adata.bins{1}(end).*ones(size(binsYZ{1}))},binsYZ);
set( surf(binsYZ{:},ratemapYZ'), 'EdgeColor', 'none');
view([90,0]);
hold('on');

zlim([-500,500]);

daspect([1,1,1]);
caxis(clim);

sax(end).XTick = [];
sax(end).YTick = [];
sax(end).Color = [0.9,0.9,0.9];
axis('off')

cax = colorbar(sax(end));
cax.Units = 'centimeters';
drawnow();
pause(0.1);
cax.Position(1) = sum(sax(end).Position([1,3]))+0.5;
colormap('parula');


%rax = rectangle('Position',[-510,-1000,1000,1700],'LineWidth',2,'EdgeColor','m');
%rax = rectangle('Position',[-490,-700,1000,1700], 'LineWidth',2,'EdgeColor','g');

line([0,0],[-1000,700],[-20,-20],'LineWidth',2,'Color',configAColor);
line([0,0],[-700,1000],[ 0, 0],'LineWidth',2,'Color',configBColor);
%line([sax(end).Position(1),sum(sax(end).Position([1,3]))]+[1,0],sax(end).Position(2).*[1,1],'EdgeColor','m')

%line([sax(end).Position(1),sum(sax(end).Position([1,3]))]+[0,-1],sax(end).Position(2).*[1,1],'Color','m')
%line([sax(end).Position(1),sum(sax(end).Position([1,3]))]+[1,0],sax(end).Position(2).*[1,1],'EdgeColor','m')

%rax = rectangle('Position',[-510,-1000,1000,1700],'LineWidth',2,'EdgeColor','m');
%rax = rectangle('Position',[-490,-700,1000,1700], 'LineWidth',2,'EdgeColor','g');


% MAZE COORDINATES CONF A & B YZ -----------------------------------------------------------------------
yOffSet = 2;
[yind, yOffSet, xind, xOffSet] = deal(secYind+2, yOffSet, secXind, 1);
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);



binsYZ      = cell([1,2]);
[binsYZ{:}] = meshgrid(pfz.adata.bins{:});
binsYZ{2} = binsYZ{2}+300;
binsYZ      = cat(2,{pfs.adata.bins{1}(end).*ones(size(binsYZ{1}))},binsYZ);
set( surf(binsYZ{:},ratemapYZA'), 'EdgeColor', 'none');
view([90,0]);
hold('on');

zlim([-500,500]);

daspect([1,1,1]);
caxis(clim);

sax(end).XTick = [];
sax(end).YTick = [];
sax(end).Color = [0.9,0.9,0.9];
axis('off')

line([0,0],[-900,900],[-40,-40],'LineWidth',2,'Color',configAColor);


[yind, yOffSet, xind, xOffSet] = deal(secYind+3, yOffSet, secXind, 2);
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);


binsYZ      = cell([1,2]);
[binsYZ{:}] = meshgrid(pfz.adata.bins{:});
binsYZ{2} = binsYZ{2}+300;
binsYZ      = cat(2,{pfs.adata.bins{1}(end).*ones(size(binsYZ{1}))},binsYZ);
set( surf(binsYZ{:},ratemapYZB'), 'EdgeColor', 'none');
view([90,0]);
hold('on');

zlim([-500,500]);

daspect([1,1,1]);
caxis(clim);

sax(end).XTick = [];
sax(end).YTick = [];
sax(end).Color = [0.9,0.9,0.9];
axis('off')

line([0,0],[-900,900],[-40,-40],'LineWidth',2,'Color',configBColor);

% PRINT it
figname = [Trial.filebase,'-vst_unit-',num2str(unit,'%04.f')];
print(gcf, '-depsc2', fullfile(figpath,[figname,'.eps']));
print(gcf,'-dpng',    fullfile(figpath,[figname,'.png']));

end

% END fabulous figure ------------------------------------------------------------------------------

%%%>>>


%%%<<< 13. Sort Units ------------------------------------------------------------------------------

Trial = MTATrial.validate(meta);

activeState = 'theta&gper&vel'; 
pfsArgs = struct('states',           activeState,          ... Computational Periods 
                 'binDims',          [30,30],              ... Physical size of bins in milimeters
                 'SmoothingWeights', [2.5,2.5],            ... Gaussian smoother prameters, std deviation in bins 
                 'numIter',          1,                    ... number of bootstraps
                 'boundaryLimits',   [-500,500;-1000,1000],... Computational domain
                 'halfsample',       false);              %... throw out half of the data for each iteration
pft = compute_ratemaps(Trial,                              ... MTATrial
                       [],                                 ... Unit list (e.g. [1, 2, ... , N])
                       @fet_rbo_AC,                        ... Function handle, ratemap space
                       [],                                 ... sampleRate, default 16Hz to speed it up
                       pfsArgs,                            ... Arguments to the MTAApfs 
                       'overwrite',true);

activeState = 'theta-rear&gper&vel'; 
pfsArgs = struct('states',           activeState,          ... Computational Periods 
                 'binDims',          [30,30],              ... Physical size of bins in milimeters
                 'SmoothingWeights', [2.5,2.5],            ... Gaussian smoother prameters, std deviation in bins 
                 'numIter',          1,                    ... number of bootstraps
                 'boundaryLimits',   [-500,500;-1000,1000],... Computational domain
                 'halfsample',       false);              %... throw out half of the data for each iteration
pfs = compute_ratemaps(Trial,                              ... MTATrial
                       [],                                 ... Unit list (e.g. [1, 2, ... , N])
                       @fet_rbo_AC,                        ... Function handle, ratemap space
                       [],                                 ... sampleRate, default 16Hz to speed it up
                       pfsArgs,                            ... Arguments to the MTAApfs 
                       'overwrite',true);

activeState = 'rear&theta&gper'; 
pfsArgs = struct('states',           activeState,          ... Computational Periods 
                 'binDims',          [30,30],              ... Physical size of bins in milimeters
                 'SmoothingWeights', [2.5,2.5],            ... Gaussian smoother prameters, std deviation in bins 
                 'numIter',          1,                    ... number of bootstraps
                 'boundaryLimits',   [-500,500;-1000,1000],... Computational domain
                 'halfsample',       false);              %... throw out half of the data for each iteration
pfr = compute_ratemaps(Trial,                              ... MTATrial
                       [],                                 ... Unit list (e.g. [1, 2, ... , N])
                       @fet_rbo_AC,                        ... Function handle, ratemap space
                       [],                                 ... sampleRate, default 16Hz to speed it up
                       pfsArgs,                            ... Arguments to the MTAApfs 
                       'overwrite',true);

% PLOT Placefields for each Unit
hfig = figure();
u = 1;
while u~=-1
    clf();
    subplot(131);
        plot(pft,u,1,'text',[],false);
        daspect([1,1,1]);
        title(num2str(u));
    subplot(132);
        plot(pfs,u,1,'text',[],false);
        daspect([1,1,1]);
        title(num2str(u));
    subplot(133);
        plot(pfr,u,1,'text',[],false);
        daspect([1,1,1]);
        title(num2str(u));
    u = figure_controls(hfig,u,Trial.spk.map(:,1));
end

% ASSIGN units manually  
unitsPyr = [1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,19,20,22,23,24,26,28,29,30,31,32,...
            33,34,36,37,39,40,44,45,46,48,49,50,51,53,54,55,56,57,58,61,62,63,64,66,68,...
            69,70,71];
set_unit_set(Trial.spk,Trial,'placefields',unitsPyr);

unitsCtr = [1,3,8,10,13,15,16,17,19,26,37,32,33,37,39,44,55,64,66,68,70];
unitsBdr = [2,4,6,7,9,11,12,14,20,22,23,24,28,29,30,31,34,36,40,45,46,47,48,49,50,51,53,54,...
            56,57,58,61,62,63,69,71];
unitsInt = [5,18,21,41,42,52,59,65,67,72,73,74];
set_unit_set(Trial.spk,Trial,'interneurons',unitsInt);

unitsSlt = [25,27,35,38,43,60];

units = unitsPyr;
%%%>>>


%%%<<< 14. Decode position -------------------------------------------------------------------------

units = get_unit_set(Trial.spk,Trial,'placefields');

% SET decoding parameters
tunits = units;
thetaRefChan = 22;
sampleRate = 30;
% $$$ phzCorr = 0;
halfSpkWindow = 0.15;
ufrWindow = 0.02;
phzBins = linspace(0,2*pi,25);
smoothingWeights = [250.^2,250.^2];

% LOAD subject
rat = Trial.load('subject','FS04_AC');
rat.resample(sampleRate);
% $$$ xyz = preproc_xyz(Trial,'trb');
% $$$ xyz.resample(sampleRate);

% $$$ phz = load_theta_phase(Trial,rat,thetaRefChan,phzCorr);


%pfs = pfs_2d_theta(Trial,tunits);
activeState = 'theta&gper&vel'; 
pfsArgs = struct('states',           activeState,          ... Computational Periods 
                 'binDims',          [30,30],              ... Physical size of bins in milimeters
                 'SmoothingWeights', [2.5,2.5],            ... Gaussian smoother prameters, std deviation in bins 
                 'numIter',          1,                    ... number of bootstraps
                 'boundaryLimits',   [-500,500;-1000,1000],... Computational domain
                 'halfsample',       false);              %... throw out half of the data for each iteration
pfs = compute_ratemaps(Trial, units, @fet_rbo_AC, [], pfsArgs, 'overwrite',true);


% GENERATE mask to remove poorly occupied spatial bins (i.e. the edges are noisy)
%mask = create_tensor_mask(pfs.adata.bins(1:2));
mask = true([pfs.adata.binSizes']);
mask([1:3,end-2:end],:) = false;
mask(:,[1:3,end-2:end]) = false;

% LOAD spikes and GENERATE unit firing rate matrix
spk = Trial.load('spk', sampleRate, 'gper', tunits, '');
ufr = Trial.load('ufr', rat,spk,tunits,ufrWindow,'boxcar',true);

% DECODE position
dc = decode_ufr_boxcar(Trial, ...
                       tunits, ...
                       sampleRate, ...
                       ufr, ...
                       pfs, ...
                       mask,...
                       halfSpkWindow,...
                       smoothingWeights,...
                       'overwrite',true);

% DECODE position with respect to theta phase
% $$$ dc = decode_ufr_binnedPhz(Trial, ...
% $$$                           tunits, ...
% $$$                           sampleRate, ...
% $$$                           ufr, ...
% $$$                           phz, ...
% $$$                           pfs, ...
% $$$                           mask,...
% $$$                           [],[],[],...
% $$$                           smoothingWeights,...
% $$$                           'overwrite',true);




% PLOT xy vs decoded xy
figure();
subplot(311); hold('on');
plot((1:size(rat))./dc.sampleRate,rat(:,1,1),'.b');
scatter(dc.ind./dc.sampleRate,dc.sax(:,1),10,dc.ucnt,'Filled');
ylim([-425,425]);
subplot(312);hold('on');
plot((1:size(rat))./dc.sampleRate,rat(:,1,2),'.b')
scatter(dc.ind./dc.sampleRate,dc.sax(:,2),10,dc.ucnt,'Filled');
ylim([-900,900]);
colormap('summer');
subplot(313);
plotSTC(Trial.stc,1,'text',{'t','A','B','C'},'rgbm');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


% PLOT error distributions 
figure();
norm = 'xprob';
subplot2(4,2,1,1);
    hist2([rat(dc.ind,1,1),dc.sax(:,1)],...
          linspace(-425,425,30),...
          linspace(-425,425,30),...
          norm);
    daspect([1,1,1]);
    xlabel('X mm');
    ylabel('X'' mm');
subplot2(4,2,2,1);
    hist2([rat(dc.ind,1,2),dc.sax(:,2)],...
          linspace(-900,900,30),...
          linspace(-900,900,30),...
          norm);
    daspect([1,1,1]);
    xlabel('Y mm');
    ylabel('Y'' mm');
subplot2(4,2,3,1);
    errorXY = sqrt(sum([sq(rat(dc.ind,1,1:2))-dc.sax].^2,2));
    set(histogram2(dc.ucnt,...
               log10(errorXY),...
               1:20,...
               linspace(0,3.2,50),...
               'DisplayStyle','tile'),...
        'EdgeColor','none');
    daspect([11,1,1]);
    ylabel('XY log10(rms error mm)');
    xlabel('decoding unit count');
subplot2(4,2,1,2);
    errorXY = sqrt(sum([sq(rat(dc.ind,1,1:2))-dc.sax].^2,2));
    xBinEds = linspace(-450,450,30);
    yBinEds = linspace(-900,900,60);
    xBinCtr = (xBinEds(2:end)+xBinEds(1:end-1))/2;
    yBinCtr = (yBinEds(2:end)+yBinEds(1:end-1))/2;
    xInd = discretize(rat(dc.ind,1,1),xBinEds);
    yInd = discretize(rat(dc.ind,1,2),yBinEds);
    nind = nniz(xInd) & nniz(yInd) & nniz(errorXY);
    errorXYgXY = accumarray([xInd(nind),yInd(nind)],errorXY(nind),[numel(xBinCtr),numel(yBinCtr)],@median);
    errorXYgXY(errorXYgXY==0) = nan;
    set(pcolor(xBinCtr,yBinCtr,log10(errorXYgXY)'),'EdgeColor','none');
    daspect([1,1,1]);
    colormap(gca(),'jet');
    colorbar(gca());
    caxis([0,3.2]);
    title('Error all Configs');
subplot2(4,2,2,2);
    confAper = [Trial.stc{'A'}];
    confAper.cast('TimeSeries',rat);
    confAper.data = logical(confAper.data(dc.ind));
    nind = nniz(xInd) & nniz(yInd) & nniz(errorXY) & confAper.data;
    errorXYgXY = accumarray([xInd(nind),yInd(nind)],errorXY(nind),[numel(xBinCtr),numel(yBinCtr)],@median);
    errorXYgXY(errorXYgXY==0) = nan;
    set(pcolor(xBinCtr,yBinCtr,log10(errorXYgXY)'),'EdgeColor','none');
    daspect([1,1,1]);
    colormap(gca(),'jet');
    colorbar(gca());
    caxis([0,3.2]);
    title('Error for Config A');
subplot2(4,2,3,2);
    confBper = [Trial.stc{'B'}];
    confBper.cast('TimeSeries',rat);
    confBper.data = logical(confBper.data(dc.ind));
    nind = nniz(xInd) & nniz(yInd) & nniz(errorXY) & confBper.data;
    errorXYgXY = accumarray([xInd(nind),yInd(nind)],errorXY(nind),[numel(xBinCtr),numel(yBinCtr)],@median);
    errorXYgXY(errorXYgXY==0) = nan;
    set(pcolor(xBinCtr,yBinCtr,log10(errorXYgXY)'),'EdgeColor','none');
    daspect([1,1,1]);
    colormap(gca(),'jet');
    colorbar(gca());
    caxis([0,3.2]);
    title('Error for Config B');    
subplot2(4,2,4,2);
    confCper = [Trial.stc{'C'}];
    confCper.cast('TimeSeries',rat);
    confCper.data = logical(confCper.data(dc.ind));
    nind = nniz(xInd) & nniz(yInd) & nniz(errorXY) & confCper.data;
    errorXYgXY = accumarray([xInd(nind),yInd(nind)],errorXY(nind),[numel(xBinCtr),numel(yBinCtr)],@median);
    errorXYgXY(errorXYgXY==0) = nan;
    set(pcolor(xBinCtr,yBinCtr,log10(errorXYgXY)'),'EdgeColor','none');
    daspect([1,1,1]);
    colormap(gca(),'jet');
    colorbar(gca());
    caxis([0,3.2]);
    title('Error for Config C');    


%%%>>>
    
%%%<<< 15. Decode position by theta phase ----------------------------------------------------------

tunits = get_unit_set(Trial.spk,Trial,'placefields');

sampleRate = 250;
thetaRefChan = 22;
phzCorr = 0
ufrWindow = 0.02;
halfSpkWindow = 0.15;
halfPhzWindow = 0.25;
phzBins = linspace(0,2*pi,25);
smoothingWeights = [250.^2,250.^2];

xyz = Trial.load('subject','FS04_AC');
xyz.resample(sampleRate);

phz = load_theta_phase(Trial,xyz,thetaRefChan,phzCorr);

%pfs = pfs_2d_theta(Trial,tunits);
activeState = 'theta&gper&vel'; 
pfsArgs = struct('states',           activeState,          ... Computational Periods 
                 'binDims',          [30,30],              ... Physical size of bins in milimeters
                 'SmoothingWeights', [2.5,2.5],            ... Gaussian smoother prameters, std deviation in bins 
                 'numIter',          1,                    ... number of bootstraps
                 'boundaryLimits',   [-500,500;-1000,1000],... Computational domain
                 'halfsample',       false);              %... throw out half of the data for each iteration
pfs = compute_ratemaps(Trial, units, @fet_rbo_AC, [], pfsArgs, 'overwrite',true);

% GENERATE mask to remove poorly occupied spatial bins (i.e. the edges are noisy)
%mask = create_tensor_mask(pfs.adata.bins(1:2));
mask = true([pfs.adata.binSizes']);
mask([1:3,end-2:end],:) = false;
mask(:,[1:3,end-2:end]) = false;

spk = Trial.load('spk', sampleRate, '', tunits, 'deburst');
ufr = Trial.load('ufr', xyz,spk,tunits,ufrWindow,'boxcar',true);

dc = decode_ufr_binnedPhz(Trial, ...
                          tunits, ...
                          sampleRate, ...
                          ufr, ...
                          phz, ...
                          pfs, ...
                          mask,...
                          halfSpkWindow,...
                          halfPhzWindow,...
                          phzBins,...
                          smoothingWeights,...
                          'overwrite',true);

dc.hRot = 0;
% DISCRETIZE theta phaze
dc.iphz = discretize(dc.phz,phzBins);
% GENERATE state collection matrix    
dc.stcm = stc2mat(Trial.stc,xyz,{'theta','vel','gper'});
dc.stcm = dc.stcm(dc.ind,:);
% COMPUTE head and body relative velocity vectors
% $$$ hvfl = fet_href_HXY(Trial,sampleRate,[],'trb');
% $$$ dc.hvfl = hvfl(dc.ind,:);



figure,
hold('on');
plot(dc.sax(:,1),dc.sax(:,2)
plot(dc.sax(:,1),dc.sax(:,2)

eang = quaternion2rad(sq(xyz(:,1,5:8)));

xind =  220300:220700;
figure,plot(real(eang(xind,2)))


rotm = quat2rotm(sq(xyz(:,1,5:8)));


nxyz = permute(multiprod(rotm,repmat([50, 0, 0],[size(rotm,1),1]),[2,3],[2]),[1,3,2])+xyz(:,1,1:3);
rxyz = permute(multiprod(rotm,repmat([ 0,-50, 0],[size(rotm,1),1]),[2,3],[2]),[1,3,2])+xyz(:,1,1:3);
zxyz = permute(multiprod(rotm,repmat([ 0, 0,50],[size(rotm,1),1]),[2,3],[2]),[1,3,2])+xyz(:,1,1:3);

figure,
hold('on')
for xind =  220300:220700;
plot3([xyz(xind,1,1)],...
     [xyz(xind,1,2)],...
     [xyz(xind,1,3)],...
      'm*');
% $$$ plot3([xyz(xind,1,1);zxyz(xind,1,1)],...
% $$$      [xyz(xind,1,2);zxyz(xind,1,2)],...
% $$$      [xyz(xind,1,3);zxyz(xind,1,3)],...
% $$$       '-b');
% $$$ plot3([xyz(xind,1,1);nxyz(xind,1,1)],...
% $$$      [xyz(xind,1,2);nxyz(xind,1,2)],...
% $$$      [xyz(xind,1,3);nxyz(xind,1,3)],...
% $$$       '-g');
plot3([xyz(xind,1,1);rxyz(xind,1,1)],...
     [xyz(xind,1,2);rxyz(xind,1,2)],...
     [xyz(xind,1,3);rxyz(xind,1,3)],...
      '-r');
plot3([xyz(xind,1,1);xyz(xind+100,1,1)],...
      [xyz(xind,1,2);xyz(xind+100,1,2)],...
      [xyz(xind,1,3);xyz(xind+100,1,3)],...      
      '-m');
end



figure,
plot(sqrt(sum((xyz(:,1,[1,2])-circshift(xyz(:,1,[1,2]),30)).^2,3)))


% CHECK velocity vector matches head direction vector
xvec = sq((xyz(:,1,[1,2])-circshift(xyz(:,1,[1,2]),30)));
xvec = bsxfun(@rdivide,xvec,sqrt(sum(xvec.^2,2)));
xang = atan2(xvec(:,1),xvec(:,2));

rvec = sq(rxyz(:,1,[1,2])-xyz(:,1,[1,2]));
rvec = bsxfun(@rdivide,rvec,sqrt(sum(rvec.^2,2)));
rang = atan2(rvec(:,1),rvec(:,2));


rvec = sq(rxyz(:,1,[1,2])-xyz(:,1,[1,2]));
rvec = bsxfun(@rdivide,rvec,sqrt(sum(rvec.^2,2)));
rang = atan2(rvec(:,1),rvec(:,2));

figure,
ind = nniz(mang) & nniz(xang);
subplot(121);
rose(circ_dist(mang(ind),xang(ind)));
subplot(122);
ind = nniz(rang) & nniz(xang);
rose(circ_dist(rang(ind),xang(ind)));

figure();
hold('on');
xind = 220300:220700;
plot3(xyz(xind,1,1),xyz(xind,1,2),xyz(xind,1,3),'.b')
plot3(nxyz(xind,1,1),nxyz(xind,1,2),nxyz(xind,1,3),'.r')
plot3(rxyz(xind,1,1),rxyz(xind,1,2),rxyz(xind,1,3),'.g')
plot3(xyz(xind(1),1,1),xyz(xind(1),1,2),xyz(xind(1),1,3),'*b')
plot3(nxyz(xind(1),1,1),nxyz(xind(1),1,2),nxyz(xind(1),1,3),'*r')
plot3(rxyz(xind(1),1,1),rxyz(xind(1),1,2),rxyz(xind(1),1,3),'*g')

eind = dc.ind>=xind(1)&dc.ind<=xind(end);
scatter3(dc.sax(eind,1),...
         dc.sax(eind,2),...
         -250.*ones([sum(eind),1]),...
         10,...
         dc.iphz(eind),...
         'Filled');
colormap('hsv');

figure,
plot(dc.esax(eind,2),dc.esax(eind,1),'.')



figure,
subplot(211);hold('on');
plot(xyz(:,1,1),'.b');
plot(dc.ind,dc.sax(:,1),'.r');
subplot(212);hold('on');
plot(xyz(:,1,2),'.b');
plot(dc.ind,dc.sax(:,2),'.r');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


% PLOT obsevered and decoded side by side
figure,
subplot(121);
plot(xyz(:,1,1),xyz(:,1,2),'.')
subplot(122);
plot(dc.sax(:,1),dc.sax(:,2),'.')
linkaxes(findobj(gcf(),'Type','Axes'),'xy');


dc.xyz =  xyz(dc.ind,1,1:3);
dc.hvec = rxyz(dc.ind,1,[1,2])-xyz(dc.ind,1,[1,2]);
dc.hvec = sq(bsxfun(@rdivide,dc.hvec,sqrt(sum(dc.hvec.^2,3))));
dc.hvec = cat(3,dc.hvec,sq(dc.hvec)*[0,-1;1,0]);
dc.hvec = multiprod(dc.hvec,                           ...
                          [cos(dc.hRot),-sin(dc.hRot); ...
                           sin(dc.hRot), cos(dc.hRot)],...
                          [2,3],                                   ...
                          [1,2]);


dc.ecom = cf(@(e,x,h,f)  ...
             sq(multiprod(bsxfun(@minus,...
                                 e(:,[1,2]),...
                                 sq(x(:,1,[1,2]))),...
                          h,...
                          2,...
                          [2,3])...
                ), ...
             {dc.com},{dc.xyz},{dc.hvec});
dc.ecom = dc.ecom{1};
% SAX 
dc.esax = cf(@(e,x,h,f)  ...
             sq(multiprod(bsxfun(@minus,...
                                 e(:,[1,2]),...
                                 sq(x(:,1,[1,2]))),...
                          h,...
                          2,...
                          [2,3])...
                ), ...
             {dc.sax},{dc.xyz},{dc.hvec});
dc.esax = dc.esax{1};
% MAX 
dc.emax = cf(@(e,x,h,f)  ...
             sq(multiprod(bsxfun(@minus,...
                                 e(:,[1,2]),...
                                 sq(x(:,1,[1,2]))),...
                          h,...
                          2,...
                          [2,3])...
                ), ...
             {dc.max},{dc.xyz},{dc.hvec});
dc.emax = dc.emax{1};

phzBinc = (phzBins(2:end)+phzBins(1:end-1))./2;
msaxf = [];
for p = 1:24,
    ind =   dc.iphz==p         ...
          & dc.ucnt>=3         ...
          & dc.stcm(:,3)==3    ...            
          & dc.stcm(:,2)==2    ...
          & dc.stcm(:,1)==1;
    msaxf(p) = median(dc.esax(ind,1));
end

figure();
plot(msaxf);



figure();
ind = dc.ucnt>=2 & ismember(dc.iphz,10:15)& dc.xyz(:,1,3)<-100;
%ind = dc.ucnt>=2 & ismember(dc.iphz,[1:8,18:24]) & dc.xyz(:,1,3)<-100;
hist2([dc.esax(ind,2),...
       dc.esax(ind,1)],...
      linspace(-300,300,30),...
      linspace(-300,300,30));



% COMPUTE the egocentric shift to match head position
ind = dc.stcm(:,1)==1 & dc.stcm(:,2)==2 & dc.stcm(:,3)==3  & dc.ucnt>=2;
xshifts = -80:5:80;
yshifts = -80:5:80;
medianPosError = nan([numel(xshifts),numel(yshifts)]);
for xshift = 1:numel(xshifts),
    for yshift = 1:numel(yshifts),
        medianPosError(xshift,yshift) = median(sqrt(sum([dc.esax(ind,1)+xshifts(xshift),...
                                                         dc.esax(ind,2)+yshifts(yshift)].^2,2)));
    end
end
mmpe = LocalMinimaN(medianPosError,100,10);
figure,
set(pcolor(yshifts,xshifts,medianPosError),'EdgeColor','none');
hold('on');
scatter(yshifts(mmpe(2)),xshifts(mmpe(1)),30,'r','Filled');
colorbar();

%%%>>>
