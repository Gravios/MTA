%% How to MTA your Data
% /storage/share/matlab/MTA/examples/example_code_FS03.m

% The Following topics will be discussed
% 
%    1. Configuration:  Setting up the MTA environment on linux
%    2. Data Arangment: Where and how to store your data
%    3. Meta Data:      Meta data for generating a session object
%    4. Preprocess CSV: parse_rbo_from_csv.m
%    5. Build Sessions: generate the MTASession and MTATrial objects

%% 1. Configuration ----------------------------------------------------------------------------------------------------
% NOTICE - if your system has never been configured for MTA follow the next step
% CONFIGURE your OS environmental variables ->  saved in ~/.bashrc 

% shell $> cd /storage/share/matlab/MTA
% shell $> ./configure -u fabian -s storage2 -d Data -t /storage2/fabian/Code 
% shell $> source ~/.bashrc

% Start Matlab
% matlab >> MTASession([]); 

%% 2. Data Arangement --------------------------------------------------------------------------------------------------
% example for FS03/FS03-20201223
%
% CREATE data directories 
% shell $> mkdir -p /storage2/fabian/Data/raw/xyz/FS03/FS03-20201223/vrr
% shell $> mkdir -p /storage2/fabian/Data/raw/xyz/FS03/FS03-20201223/vrr

%% START HERE
% $$$ MTA_DATA_PATH = getenv('MTA_DATA_PATH');
% $$$ meta.path.raw.ephys       = fullfile(MTA_DATA_PATH,'raw/ephys/',      meta.sessionBase,meta.sessionName);
% $$$ meta.path.raw.xyz         = fullfile(MTA_DATA_PATH,'raw/xyz/',        meta.sessionBase,meta.sessionName,meta.mazeName);
% $$$ meta.path.processed.xyz   = fullfile(MTA_DATA_PATH,'processed/xyz/',  meta.sessionBase,meta.sessionName);
% $$$ meta.path.processed.ephys = fullfile(MTA_DATA_PATH,'processed/ephys/',meta.sessionBase,meta.sessionName);

mkdir('/storage2/fabian/Data/raw/xyz/',
% shell $> mkdir -p /storage2/fabian/Data/processed/xyz/FS03/FS03-20201223/vrr
% shell $> mkdir -p /storage2/fabian/Data/processed/ephys/FS03/FS03-20201223
%
% shell $> mkdir -p /storage2/fabian/Data/processed/ephys/
% shell $> mv /storage2/fabian/Data/processed/FS03 /storage2/fabian/Data/processed/ephys/
% 
% shell $> cd /storage2/fabian/Data/processed/ephys/FS03
% shell $> mv ./FS03_2020_12_23 ./FS03-20201223
%
% MOVE data to data directories 
% shell $> mv /storage2/fabian/Data/processed/ephys/FS03/FS03_2020_12_23 /storage2/fabian/Data/processed/ephys/FS03/FS03-20201223
% shell $> cd /storage2/fabian/Data/processed/ephys/FS03/FS03-20201223
%
% RENAME data files
% NOTICE - the flag '-n' lists the filename changes without changing them
% shell $> rename -n 's/(\w{4})_(\d{4})_(\d{2})_(.*)/$1-$2$3$4/' ./*
%
% WARNING - The following line will rename all files in the current directory
% shell $> rename 's/(\w{4})_(\d{4})_(\d{2})_(.*)/$1-$2$3$4/' ./*
%
% MOVE motion tracking csv files to Data/raw/xyz
% NOTICE - assuming the maze vrr is apporiate for the virtual arena setup
% shell $> cp /storage2/fabian/Data/processed/ephys/FS03/FS03-20201223/Take*.csv /storage2/fabian/Data/raw/xyz/FS03/FS03-20201223/vrr



%% 3. Meta Data --------------------------------------------------------------------------------------------------------
% SETUP struct which holds session metadata
% struct Session  
%--------------------------------------------------------------
%|  sessionName    |  String        |   'jg05-20120309'       | Warning: sessionName must match data directory names 
%|  subjects       |  CellStr       |   {{'jg05'}}            | 
%|  mazeName       |  String        |   'cof'                 | Warning: Maze name must be explicitly registered in MTA/MTAConfiguration.m
%|  trialName      |  String        |   'all'                 | Warning: 'all' is default where the entire session is considered a trial
%|  dLoggers       |  CellStr       |   {{'nlx','vicon'}}     | Warning: MTA/utilities/synchronization must contain a sync_DataLogger1_DataLogger2.m
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
%|  thetaRef       |  Array[double] |   [8:8:64]              | WARNING: Only one channel per electrode shank
%| thetaRefGeneral |  Double        |   8                     | WARNING: Only one channel 
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

% DEFAULT stuff
meta.sessionBase = 'FS03';
meta.sessionName = 'FS03-20201223';

meta.primarySubject = 'FS03';

meta.subjects(1).name = 'FS03';
meta.subjects(1).type = 'rbo'
%meta.subjects(1).type = 'xyz'
meta.subjects(1).rb(1).name = 'Head';
meta.subjects(1).rb(1).alias = 'RatW';
meta.subjects(2).name = 'Arena';
meta.subjects(2).type = 'rbo'
%meta.subjects(1).type = 'xyz'
meta.subjects(2).rb(1).name = 'Arena';
meta.subjects(2).rb(1).alias = 'Arena';

%meta.subjects(1).rb(1).nMarkers = 4;
% $$$ meta.subjects(1).rb(2).name = 'body';
% $$$ meta.subjects(1).rb(2).alias = 'RatB';
% $$$ meta.subjects(1).rb(2).nMarkers = 4;
% $$$ 
% $$$ meta.subjects(2).name = 'FS04';
% $$$ meta.subjects(2).rb(1).name = 'head';
% $$$ meta.subjects(2).rb(1).alias = 'FS04_head';
% $$$ meta.subjects(2).rb(1).nMarkers = 4;

meta.mazeName = 'vrr';
meta.trialName = 'all';
meta.dLoggers = {'WHT','CSV'};
meta.dPaths.xyz   = fullfile('/storage/gravio/data/processed/xyz/',  meta.subjects(1).name);
meta.dPaths.ephys = fullfile('/storage/gravio/data/processed/ephys/',meta.subjects(1).name);
meta.xyzSampleRate = 120.00;% ???
meta.hostServer = 'lmu';
meta.dataServer = 'lmu';
meta.project    = 'general';
% $$$ meta.TTLValue   = [];
% $$$ meta.TTLValue = {                                 ...
% $$$     {48218111,''},                                ...
% $$$     {47050073,'Take 2020-12-23 12.51.38 AM.csv'}, ...
% $$$     {22460889,''},                                ...
% $$$     {58916569,''},                                ...
% $$$     {46969305,'Take 2020-12-23 08.55.44 PM.csv'}  ...
% $$$ };

meta.TTLValue = {                                   ...
    {48218111,'',''},                               ...
    {47050073,'FS03-20201223.take_0001.mat','FS03'},...
    {22460889,'',''},                               ...
    {58916569,'',''},                               ...
    {46969305,'FS03-20201223.take_0002.mat','FS03'} ...
};

% Trials
meta.includeSyncInd = [];
meta.offsets  = [0,0];
% Maze corrections if the maze isn't already centered
meta.xOffSet  = 0;
meta.yOffSet  = 0;
meta.rotation = 0;
% State label collection
meta.stcMode  = 'default';
% LFP Theta reference channels
meta.thetaRef = [1:11:64];
meta.thetaRefGeneral = 1;

meta.csv = {'','Take 2020-12-23 12.51.38 AM.csv','','','Take 2020-12-23 08.55.44 PM.csv'};

% WHITEHAT stuff
% Link WhiteHat -> motive
% Get automatically framo FS03-20201223.srs
MTA_DATA_PATH = getenv('MTA_DATA_PATH');
meta.path.raw.ephys       = fullfile(MTA_DATA_PATH,'raw/ephys/',      meta.sessionBase,meta.sessionName);
meta.path.raw.xyz         = fullfile(MTA_DATA_PATH,'raw/xyz/',        meta.sessionBase,meta.sessionName,meta.mazeName);
meta.path.processed.xyz   = fullfile(MTA_DATA_PATH,'processed/xyz/',  meta.sessionBase,meta.sessionName);
meta.path.processed.ephys = fullfile(MTA_DATA_PATH,'processed/ephys/',meta.sessionBase,meta.sessionName);

meta.TTLValue = convert_srs_to_TTLValue(meta);


%% 4. Preprocess CSV ---------------------------------------------------------------------------------------------------
parse_rbo_from_csv(meta);




% DEV NOTE ---------------------------------------------------------------------
% INTEGRITY CHECK - data may contain nans where model was temporarily lost
%    nanGaps = isnan(data(:,3));
%    nanTot = sum(nanGaps);

% ONCE extracted create MTADrbo.m data type for rigidBodyObjects
% TODO How to represent and setup rbo is the meta data
%      subject has one or more rbo.
%        subjecs may have rbo(s) with same FINAL name (e.g. head or body)
%        MAP motive rbo names to final common naming scheme 
%        OR decide upon naming scheme beforehand 
% FORNOW 
%      name RatW -> head
%
% TODO load/reference scheme
% rbo = [Trial.load('skeleton','FS03').rbo{'head'}]
% rat = Trial.load('skeleton','FS03');
%
% PROBLEM 
%    DESCRIPTION unless corrections are applied at data acquisition to correctly 
%        orient the rbo relative to the marker
%    CONDITONS 
%        if rbo construction originates from marker set?
%        if rbo data is provided but no correction at time of acquisition?
%
% END DEV NOTE -----------------------------------------------------------------





%% 5. Build Sessions ---------------------------------------------------------------------------------------------------

% CREATE session objects
build_sessions(meta);
% -> generates symbolic link file structure in project folder
% -> synchronizes lfp and rbo/xyz objecs
% -> creates default objects

%% 6. Loading Sessions -------------------------------------------------------------------------------------------------
% Meta data struct
Session = MTASession.validate(meta);

% Piecing together the session name 
Session = MTASession.validate([meta.sessionName,'.',meta.mazeName,'.',meta.trialName]);

% writing out the session name 
Session = MTASession.validate('FS03-20201223.vrr.all');

Session.model = Session.subject.model;
Session.save();

QuickTrialSetup(meta);

Trial = MTATrial.validate([meta.sessionName,'.',meta.mazeName,'.',meta.trialName]);

rat   = Session.load('subject','FS03');
Arena = Session.load('subject','Arena');

% CORRECT for arena position
rat.data(:,1,1:3) = bsxfun(@minus,rat.data(:,1,1:3),mean(Arena.data(nniz(Arena.data),1,1:3)));
save(rat);

% PLOT the height of the rat's head
figure,plot(rat(:,'Head',3));

% PLOT xy position of the rat's head
figure();
plot(rat.data(:,1,1),rat.data(:,1,2),'.');

% CHECK arena movement for coresponding rat behavior
figure();
plot(Arena(:,'Arena',3));

%label_theta(Session);
rat = Trial.load('subject','FS03');

gper = ThreshCross(double(nniz(rat(:,1,1))),0.5,10);
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   gper,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'gper','a');



ratFilt = filter(copy(rat),'ButFilter',4,1,'low');

figure();plot(sqrt(sum(diff(ratFilt(:,1,1:2)).^2,3)).*ratFilt.sampleRate./10)

vper = ThreshCross(double(sqrt(sum(diff(ratFilt(:,1,1:2)).^2,3)).*ratFilt.sampleRate./10>2),0.5,10);

Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   vper,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'vel','v');



% SUBSET of xyz within active stae
% p(s)
activeState = 'theta&gper&vel';
pfsArgs = struct('states',           activeState,                  ...
                 'binDims',          [20,20],                  ...
                 'SmoothingWeights', [3.5,3.5],                           ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-500,500;-1000,1000],          ...
                 'halfsample',       false);
pfs = compute_ratemaps(Trial,[],@fet_rbo,[],pfsArgs);

figure,
for u = 1:144
plot(pfs,u,1,'text',[],false);
title(num2str(u));
waitforbuttonpress();
end


xyz = Session.load('xyz');


    
%plots the height of the front marker on the rats head against time
plot(Session.xyz(:,Session.trackingMarker,3));
    %plots x versus y 
    plot(Session.xyz(:,Session.trackingMarker,1),Session.xyz(:,Session.trackingMarker,2),'.');

    
    % Returns a Trial object from the data of the Session object 
    % and do automatic behavioral segmentation
    Trial = QuickTrialSetup(Session,         ...
                            TrialName,       ...
                            startStopShift,  ...
                            ignoredViconTrials);

end

%% Loading stuff

% This loads a trial
%Trial = MTATrial.validate('jg05-20120312.cof.all');
Trial = MTATrial (SessionName,    TrialName,MazeName);


% Does the same thing as above in a single line.
xyz = Trial.load('xyz');


% You can filter data with an arbitrary window
%xyz.filter(ones([1,7])./7);
xyz.filter('ButFilter',3,50,'low'); %see help gtwin for gaussian kernals

% this is crapy though since the angles are not easily smoothed in
% the time domain, so we can create a new set based on a smoothed xyz
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,50,'low'); 
ang = create(MTADang,Trial,xyz);

% Marker segment angles are relative to the room coordinate system.
% The following code plots the distribution of head pitch during rearing
figure,hist(ang(Trial.stc{'r'},'head_back','head_front',2),100)

% Doing the same thing with numerical indexing
figure,hist(ang(Trial.stc{'r'},5,7,2),100)


%% Loading LFP

% for jg05 how to load the raw lfp of the hippocampus H64BUZ
lfp = Trial.load('lfp',1:64);

% for jg05 how to load the raw lfp of the hippocampus H32LIN
lfp = Trial.load('lfp' ,65:96);

lfp = Trial.load('lfp' ,70);

lfp.resample(xyz);

hang = create(MTADang,Trial,xyz);


pos = sq(bsxfun(@minus,sq(xyz(:,'head_back',[1:2])),[-100,10]));
ang = cell([1,2]);
[ang{:}] = cart2pol(pos(:,1),pos(:,2));
ang = cat(2,ang{:});

figure,plot(circ_dist(hang(:,'head_back','head_front',1),ang(:,1)),ang(:,2),'.')





%% Behavior Segmentation stuff - It's automatic if done by Quick trial setup
% most behaviors are stored as periods in the stc field of a Trial
% Stc should be a MTAStateCollection object.
% you can see which states are available by the following commands
disp(Trial.stc.list_state_attrib('label'))
disp(Trial.stc.list_state_attrib('key'))

% You can reference timeperiods of MTADepoch objects. These are
% stored in the Trial.stc property: stc = "State Collection".
% The periods will automatically be resampled if they don't match
% the sampling rate of the calling object.
figure,hist(ang(Trial.stc{'r'},'head_back','head_front',2),100)





%% Place Fields

% I'll comment on this another time
%units = select_units(Trial,18,'pyr');

pfs = MTAAknnpfs(Trial,...
                 'units',            [],...     % [] = all units (numbered by their clu
                 'states',           'walk',... % the state you're interested in
                 'overwrite',        true,...   % this will redo the calculations for
                                          ...     the specified units
                 'numIter',          1,...
                 'ufrShufBlockSize', 0,...
                 'binDims',          [30,30],...
                 'distThreshold',    125,...
                 'nNearestNeighbors',110);








