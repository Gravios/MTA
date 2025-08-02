function process_session_FS(meta,varargin)

% TESTVARS -----------------------------------------------------------------------------------------
meta.sessionBase = 'FS04';
meta.sessionName = 'FS04-20210318a';
meta.subjectName = 'FS04';
meta.mazeName    = 'vrr';
meta.mazeList    = {'vrr'};
meta.xml         = 'FS04.xml';
%---------------------------------------------------------------------------------------------------

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('test', 'echo '...
);
[test] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

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

% rawData : rename CSV files one per ephys
% rawData -> raw/ephys: link dat files
% rawData -> raw/xyz  rename CSV files one per ephys
%

meta.path = generate_mta_paths(meta);

hswfiles = importdata(fullfile(meta.path.raw.data,[meta.sessionName,'.txt']));


for file = reshape(hswfiles(3:end),1,[])
    file = file{1};
    % MOVE and RENAME CSV files to xyz.raw
    for maze = meta.mazeList
        maze = maze{1};
        if exist(fullfile(meta.path.raw.data,file,maze))
            csvfile = dir(fullfile(meta.path.raw.data,file,'*.csv'));
            movefile(fullfile(csvfile.folder,csvfile.name),fullfile(csvfile.folder,[file,'.csv']));
            system([test,'ln -s ',fullfile(meta.path.raw.data,file,[file,'.csv ']), ...
                     fullfile(meta.path.raw.xyz.(maze),[file,'.csv'])]);
        end
    end
    % SYMLINK dat files to processed.ephys
% $$$     system([test,'ln -s ',fullfile(meta.raw.data,file,[file,'.dat']),' ',...
% $$$             fullfile(meta.path.processed.ephys,[file,'.dat'])]);
end

% COPY xml file for ndm processing
system([test,'cp ', fullfile(meta.raw.ephys, meta.xml), ' ', ...
                    fullfile(meta.path.processed.ephys, [sessionName,'.xml'])]);

currentWorkingDirectory = pwd();
cd(meta.path.processed.ephys);

%system(['hsw_rename_files_to_nlx ' meta.path.raw.data]);
system(['merge_hsw ',  sessionName])
system(['process_ndm ',sessionName])


%% event stuff

% FS04-20210326a.evt.beacon
% FS04-20210326a.csv.1
% $$$ for f = reshape(hswfiles,1,[])
% $$$     meta.raw.ephys
% $$$ 
% $$$ % RENAME csv files
% $$$     csvName = {'Take 2021-03-26 08.45.18 PM.csv'};
% $$$     csvDate = struct('year',    '',...
% $$$                      'month',   '',...
% $$$                      'day',     '',...
% $$$                      'hour',    '',...
% $$$                      'minute',  '',...
% $$$                      'second',  '',...
% $$$                      'meridian','');
% $$$     pat = ['Take\s',               ...
% $$$            '(?<year>\d*)[-]',      ...
% $$$            '(?<month>\d*)[-]',     ...
% $$$            '(?<day>\d*)\s',        ...
% $$$            '(?<hour>\d*)[.]',      ...
% $$$            '(?<minute>\d*)[.]',    ...
% $$$            '(?<second>\d*)\s',     ...
% $$$            '(?<meridian>\w*)[.]csv'];
% $$$ 
% $$$     csvDate = regexp(csvName,pat,'names');
% $$$     if strcmp(csvDate.meridian, 'PM')
% $$$         csvDate.hour = num2str(str2num(csvDate.hour)+12);
% $$$     end
% $$$     tmpCsvName = ['Take_',...
% $$$                   csvDate.year,csvDate.month,csvDate.day,'_',...
% $$$                   csvDate.hour,csvDate.minute,csvDate.second,'.csv'];
% $$$     movefile(fullfile(meta.path.raw.xyz,csvName),...
% $$$              fullfile(meta.path.raw.xyz,tmpCsvName));
% $$$ end



% $$$ 'HSW_2021_03_26__20_46_13__31min_07sec__hsamp_64ch_25000sps' 'FS04-20210326.Take001.csv'
% $$$ 'HSW_2021_03_26__21_22_21__42min_14sec__hsamp_64ch_25000sps' 'FS04-20210326.Take002.csv'



% $$$ srs = importdata('/storage/gravio/data/processed/ephys/FS04/FS04-20210318a/FS04-20210318a.srs')
% $$$ 
% $$$ for hsw = srs.textdata'
% $$$     hsw = hsw{1};
% $$$     disp(hsw);
% $$$     % regexp -> date info
% $$$ end


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
m = 1;
%meta.sessionBase = meta.subjectName;
%meta.sessionName = meta.sessionName;
%meta.mazeName    = meta.mazeName;
%meta.trialName = 'all';
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
meta.rotation = 0;
% SET State label collection
meta.stcMode  = 'default';
% SET LFP Theta reference channels
meta.thetaRef = [1:11:64];
meta.thetaRefGeneral = 1;

meta.srs = importdata(fullfile(meta.path.processed.ephys,[meta.sessionName,'.srs']));


% Constructed as cell array with empty strings where no mocap occured.
% WARNING 1 take per record
% EXAMPLE : given 3 records exist in a session
%           and record 1 and 3 have associated csv files
%           the meta.csv field should look like this 
%
% meta.csv = {'Take 2021-03-22 08.24.46 PM.csv',     '', 'Take 2021-03-22 09.24.46 PM.csv'};
%                      record1,                 record2,              record3
%
%meta.csv = {'Take 2021-03-19 07.48.56 PM.csv',''};

meta.csv = repmat({''},[1,length(hswfiles)]);
csvList = dir(fullfile(meta.path.raw.xyz.(maze),['*.csv']));

if ~isempty(csvList),
    for f = 1:numel(hswfiles)
        for c = 1:numel(csvList)
            if strcmp([hswfiles{f},'.csv'],csvList(c).name)
                meta.csv{f} = csvList(c).name;
            end
        end
    end
end


% MOVE csv file from the raw ephys folder to raw xyz folder
%cf(@(csv) movefile(fullfile(meta.path.processed.ephys,csv),meta.path.raw.xyz), meta.csv(~isempty(meta.csv)));


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
% $$$ for maze = meta.mazeList
% $$$     maze = maze{1};

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
%Session = MTASession.validate('FS04-20210319a.vrr.all');




%% 7. Build Trials -----------------------------------------------------------------------------------------------------
QuickTrialSetup(meta);



%% 8. Load Trials ------------------------------------------------------------------------------------------------------
% NOTE the next 3 commands do the same thing
% Meta data struct
Trial = MTATrial.validate(meta);
% Piecing together the session name 
Trial = MTATrial.validate([meta.sessionName,'.',meta.mazeName,'.',meta.trialName]);
% writing out the session name 
%Trial = MTATrial.validate('FS04-20210318a.vrr.all');



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
figure();
plot(ratAC(:,1,1),ratAC(:,1,2),'.')



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
Trial.stc.states{end}.save()
                   

% LABEL speed thresholded periods
ratFilt = filter(copy(rat),'ButFilter',4,1,'low');
%figure();plot(sqrt(sum(diff(ratFilt(:,1,1:2)).^2,3)).*ratFilt.sampleRate./10)
vper = ThreshCross(double(sqrt(sum(diff(ratFilt(:,1,1:2)).^2,3)).*ratFilt.sampleRate./10>5),0.5,10);
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   vper,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'vel','v');
Trial.stc.states{end}.save()


%% 10. Compute ratemaps ------------------------------------------------------------------------------------------------
% SUBSET of xyz within active stae
% p(s)
% 
%
activeState = 'theta&gper&vel'; 

% ROOM frame of reference  (X,Y)
pfsArgs = struct('states',           activeState,                  ... Computational Periods 
                 'binDims',          [20,20],                      ... Physical size of bins in milimeters
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
subplotHandles = gobjects([1,4]);             % INITIALIZE empty GraphicsPlaceholder array (matlab::graphics::objectsystem)
subplotHandles = tight_subplot(2,2,0.1);      % GENERATE 1x2 subplot grid                  (labbox::Graphics)
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
    hfig.CurrentAxes = ax;                    % SET focus on first set of axes
    cla(ax);                                  % CLEAR the contents from axes (ax)          (matlab::graphics)
    ax.Title.String = '';
    plot   (pfz,unit,1,'text',[],false);      % PLOT the place field of unit               (MTA::MTAApfs)
    daspect(ax,[1,1,1]);                      % SET the aspect ratio                       (matlab::graph3d)
    title  (ax,num2str(unit));                % SET the title of the axes (ax)             (matlab::graph2d)
    
    ax = subplotHandles(4);
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
    xlim(ax,[0,size(mspkt,1)]);                % SET the Y axis visible range               (matlab::graph3d)
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







