%% How to MTA your Data


% The Following topics will be discussed
% 
%    1. Configuration:      Setting up the MTA environment on linux
%    2. Data Arangment:     Where and how to store your data
%    3. Meta Data:          Meta data for generating a session object
%    4. Preprocess CSV:     parse_rbo_from_csv.m
%    5. Build Sessions:     generate the MTASession objects
%    6. Loading Sessions:   Various methods for loading an MTASession object
%    7. Build Sessions:     generate the MTATrial objects
%    8. Load Subjects:      Load arena and subject rbo(s) 
%    9. Label Behaviors:    Label basic behaviors (theta periods, good periods, speed thresholded periods)
%   10. Compute ratemaps:   Do it.

%% 1. Configuration ----------------------------------------------------------------------------------------------------
% NOTICE - if your system has never been configured for MTA follow the next step
% CONFIGURE your OS environmental variables ->  saved in ~/.bashrc 
%
% shell $> cd /storage/share/matlab/MTA
% shell $> ./configure -u fabian -s storage2 -d data -t /storage2/fabian/Code 



%% 2. Data Arangement --------------------------------------------------------------------------------------------------
% example for FS04/FS04-20210319a
%
% CREATE data directories 
!mkdir -p /storage2/fabian/data/raw/xyz/FS04/FS04-20210319a/vrr
!mkdir -p /storage2/fabian/data/processed/xyz/FS04/FS04-20210319a/vrr
!mkdir -p /storage2/fabian/data/processed/ephys/FS04/FS04-20210319a
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

% METADATA 
% SET Session info
meta.sessionBase = 'FS04';
meta.sessionName = 'FS04-20210319a';
meta.mazeName = 'viz';
meta.trialName = 'all';
meta.dLoggers = {'WHT','CSV'};
meta.dPaths.xyz   = fullfile('/storage2/fabian/data/processed/xyz/',  meta.sessionBase);
meta.dPaths.ephys = fullfile('/storage2/fabian/data/processed/ephys/',meta.sessionBase);
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
meta.rotation = 0;
% SET State label collection
meta.stcMode  = 'default';
% SET LFP Theta reference channels
meta.thetaRef = [1:11:64];
meta.thetaRefGeneral = 1;

% Constructed as cell array with empty strings where no mocap occured.
% WARNING 1 take per record
% EXAMPLE : given 3 records exist in a session
%           and record 1 and 3 have associated csv files
%           the meta.csv field should look like this 
%
% meta.csv = {'Take 2021-03-22 08.24.46 PM.csv',     '', 'Take 2021-03-22 09.24.46 PM.csv'};
%                      record1,                 record2,              record3
%
meta.csv = {'Take 2021-03-19 07.48.56 PM.csv',''};

MTA_DATA_PATH = getenv('MTA_DATA_PATH');
meta.path.raw.ephys       = create_directory(fullfile(MTA_DATA_PATH,'raw/ephys/',      meta.sessionBase,meta.sessionName));
meta.path.raw.xyz         = create_directory(fullfile(MTA_DATA_PATH,'raw/xyz/',        meta.sessionBase,meta.sessionName,meta.mazeName));
meta.path.processed.xyz   = create_directory(fullfile(MTA_DATA_PATH,'processed/xyz/',  meta.sessionBase,meta.sessionName,meta.mazeName));
meta.path.processed.ephys = create_directory(fullfile(MTA_DATA_PATH,'processed/ephys/',meta.sessionBase,meta.sessionName));



% MOVE csv file from the raw ephys folder to raw xyz folder
cf(@(csv) movefile(fullfile(meta.path.processed.ephys,csv),meta.path.raw.xyz), meta.csv(~isempty(meta.csv)));


% MOVE csv file from the raw ephys folder to raw xyz folder

% NOTE : not really a TTLValue
meta.TTLValue = convert_srs_to_TTLValue(meta);    
% EXPECTED output
% meta.TTLValue = { ...                                record1
%                   {[46029977]},...                     record sample count
%                   {'FS04-20210323a.take_0001.mat'},... file with extracted rigidbodies
%                   {'FS04'}...                          session's primary subject
% }



%% 4. Preprocess CSV ---------------------------------------------------------------------------------------------------
%
% EXTRACTS rigid bodies from the csv file generated by motive
%
parse_rbo_from_csv(meta);



%% 5. Build Sessions ---------------------------------------------------------------------------------------------------
build_sessions(meta);
% -> generates symbolic link file structure in project folder
% -> synchronizes lfp and rbo/xyz objecs
% -> creates default objects



%% 6. Loading Sessions -------------------------------------------------------------------------------------------------
% NOTE the next 3 commands do the same thing
% Meta data struct
Session = MTASession.validate(meta);
% Piecing together the session name 
Session = MTASession.validate([meta.sessionName,'.',meta.mazeName,'.',meta.trialName]);
% Writing out the session name 
Session = MTASession.validate('FS04-20210319a.vrr.all');




%% 7. Build Trials -----------------------------------------------------------------------------------------------------
QuickTrialSetup(meta);



%% 8. Load Trials ------------------------------------------------------------------------------------------------------
% NOTE the next 3 commands do the same thing
% Meta data struct
Trial = MTATrial.validate(meta);
% Piecing together the session name 
Trial = MTATrial.validate([meta.sessionName,'.',meta.mazeName,'.',meta.trialName]);
% writing out the session name 
Trial = MTATrial.validate('FS04-20210319a.vrr.all');



%% 8. Load Subjects ----------------------------------------------------------------------------------------------------
% problem: getting the modified subject objects passed to the pfs computation
% problem: loading Objects

rat   = Trial.load('subject','FS04');
Arena = Trial.load('subject','Arena');


% CONVERT the xyz coordinates to Arena frame of refrences centered on the Arena
ratAC = copy(rat);
% TRANSLATE 
ratAC.data(:,1,1:3) = ratAC(:,'Head',1:3) - Arena(:,'Arena',1:3);
% ROTATE 
o = 0.0349065850398866;
ratAC.data(:,1,1:2) = multiprod([cos(o),-sin(o);sin(o),cos(o)],sq(ratAC.data(:,1,1:2)),[1,2],[2]);
ratAC.label = [ratAC.label,'_AC'];
ratAC.update_filename([Trial.filebase,'.rbo.',rat.name,'_AC.s.mat']);
ratAC.save();


%% 9. Label Behaviors --------------------------------------------------------------------------------------------------
% LABEL 

% LABEL theta periods from lfp
label_theta(Session,[],32);

% LABEL Non nan periods
gper = ThreshCross(double(nniz(rat(:,1,1))),0.5,10);
Trial.stc.addState(Trial.spath,                     ... path to project folder
                   Trial.filebase,                  ... filebase
                   gper,                            ... state periods [n,2] 
                   rat.sampleRate,                  ... period sample rate
                   Trial.sync.copy,                 ... Trial synchronization object
                   Trial.sync.data(1),              ... Trial synchronization origin
                   'gper',                          ... state label
                   'a');


% LABEL speed thresholded periods
ratFilt = filter(copy(rat),'ButFilter',4,1,'low');
%figure();plot(sqrt(sum(diff(ratFilt(:,1,1:2)).^2,3)).*ratFilt.sampleRate./10)
vper = ThreshCross(double(sqrt(sum(diff(ratFilt(:,1,1:2)).^2,3)).*ratFilt.sampleRate./10>2),0.5,10);
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   vper,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'vel','v');



%% 10. Compute ratemaps ------------------------------------------------------------------------------------------------
% SUBSET of xyz within active stae
% p(s)
% 
%
activeState = 'theta&gper&vel'; 

% ROOM frame of reference  (X,Y)
pfsArgs = struct('states',           activeState,                  ... Computational Periods 
n                 'binDims',          [20,20],                      ... Physical size of bins in milimeters
                 'SmoothingWeights', [3.5,3.5],                    ... Gaussian smoother prameters, std deviation in bins 
                 'numIter',          1,                            ... number of bootstraps
                 'boundaryLimits',   [-500,500;-1000,1000],        ... Computational domain
                 'halfsample',       false);                      %... throw out half of the data each iteration ... or don't

pfs = compute_ratemaps(Trial,                                      ... MTATrial
                       [],                                         ... Unit list (e.g. [1, 2, ... , N])
                       @fet_rbo_AC,                                ... Function handle, ratemap space
                       [],                                         ... sampleRate, default 16Hz to speed it up
                       pfsArgs,                                    ... Arguments to the MTAApfs 
                       'overwrite',true);

                       
% ARENA frame of reference (Z,Y)
pfsArgs = struct('states',           activeState,                  ...
                 'binDims',          [20,20],                      ...
                 'SmoothingWeights', [3.5,3.5],                    ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-1000,1000;-300,0],        ...
                 'halfsample',       false);
pfz = compute_ratemaps(Trial,[],@fet_rbo_YZ,[],pfsArgs,'overwrite',true);



% LOAD spk object - holds the spike times and cluster ids ( see MTASpk )
% help MTASpk
spkw = Trial.spk.copy();
spkw.load_spk(Trial);
                       



% PLOT Placefields for each Unit
spkScale = 0.2;
hfig = figure(101);
subplotHandles = gobjects([1,3]);             % INITIALIZE empty GraphicsPlaceholder array (matlab::graphics::objectsystem)
subplotHandles = tight_subplot(1,3,0.1);      % GENERATE 1x2 subplot grid                  (labbox::Graphics)


unit = 1;
while unit~=-1,
    ax = subplotHandles(1);
    hfig.CurrentAxes = ax;                    % SET focus on first set of axes
    cla(ax);                                  % CLEAR the contents from axes (ax)          (matlab::graphics)
    ax.Title.String = '';
    plot   (pfs,unit,1,'text',[],false);      % PLOT the place field of unit               (MTA::MTAApfs)
    daspect(ax,[1,1,1]);                      % SET the aspect ratio                       (matlab::graph3d)
    title  (ax,num2str(unit));                % SET the title of the axes (ax)             (matlab::graph2d)

    ax = subplotHandles(2);
    hfig.CurrentAxes = ax;                    % SET focus on first set of axes
    cla(ax);                                  % CLEAR the contents from axes (ax)          (matlab::graphics)
    ax.Title.String = '';
    plot   (pfz,unit,1,'text',[],false,'flipAxesFlag',true);      % PLOT the place field of unit               (MTA::MTAApfs)
    daspect(ax,[1,1,1]);                      % SET the aspect ratio                       (matlab::graph3d)
    title  (ax,num2str(unit));                % SET the title of the axes (ax)             (matlab::graph2d)

    ax = subplotHandles(3);
    cla(ax);
    uResw = spkw.res(spkw.clu==unit);         % SELECT the spike times where the clusterId equals unit
    uSpkw = spkw.spk(spkw.clu==unit,:,:);     % SELECT the spike waveforms where the clusterId equals unit
    [~,sInd] = SelectPeriods(uResw,[Trial.stc{activeState,1}],'d',1,0); % SELECT spike times within the activeState

    if numel(sInd)<=1, mspkt = zeros([size(uSpkw,2),size(uSpkw,3)]);    % COMPUTE the mean spike waveform
    else,              mspkt = bsxfun(@plus,sq(mean(uSpkw(sInd,:,:)))'.*spkScale,fliplr(linspace(100,800,size( uSpkw,2))));
    end
    
    if numel(sInd)>5,
        plot(ax,mspkt,'b');                   % PLOT the mean spike waveform               (matlab::graph2d) 
    end
    xlim(ax,[0,size(mspk,1)]);                % SET the Y axis visible range               (matlab::graph3d)
    ylim(ax,[0,900]);                         % SET the Y axis visible range               (matlab::graph3d)
    set(ax,'YTickLabel',{});                  % REMOVE labels along the Y axis             (matlab::graphics::objectsystem)
    set(ax,'XTickLabel',{});                  % REMOVE labels along the X axis             (matlab::graphics::objectsystem)
    
    unit = figure_controls(hfig,unit,Trial.spk.map(:,1));
end


compute_neuron_quality(Session);

% NOTES ----------------------------------------------------------------------------------------------------------------
% graphical handles inherit the hgsetget superclass, 
% which allows you to use the get and set functions.
%
% Try using the get function on a handle from a figure/axes/plot/line 
%  EX: get(hfig)







