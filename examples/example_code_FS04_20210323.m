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
% shell $> ./configure -u fabian -s storage2 -d Data -t /storage2/fabian/Code 



%% 2. Data Arangement --------------------------------------------------------------------------------------------------
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
meta.sessionName = 'FS04-20210323a';
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
link_session( meta.sessionName, meta.dPaths);
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
Session = MTASession.validate('FS04-20210323a.vrr.all');

% $$$ Session.load('spk');
% $$$ Session.save();


%% 7. Build Trials -----------------------------------------------------------------------------------------------------
QuickTrialSetup(meta);



%% 8. Load Trials ------------------------------------------------------------------------------------------------------
% NOTE the next 3 commands do the same thing
% Meta data struct
Trial = MTATrial.validate(meta);
% Piecing together the session name 
Trial = MTATrial.validate([meta.sessionName,'.',meta.mazeName,'.',meta.trialName]);
% writing out the session name 
Trial = MTATrial.validate('FS04-20210323a.vrr.all');



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

ratAC   = Trial.load('subject','FS04_AC')
ratAC.data(:,1,3) = ratAC.data(:,1,3) + 270;


%% 9. Label Behaviors --------------------------------------------------------------------------------------------------
% LABEL 

% LABEL theta periods from lfp
label_theta(Session);


beacons = importdata(fullfile(Trial.spath,'BPositions_FS4_20210323-200748','beacons_20210323-200748.txt'));
bPos = importdata(fullfile(Trial.spath,'BPositions_FS4_20210323-200748','position_20210323-200748.txt'));

bTimes = round(interp1(bPos(:,1),bPos(:,8),beacons(:,1)));
diff(bTimes)/ratAC.sampleRate
runId = zeros([size(bTimes)]);
runId(1:2:end) = 1;
runId(2:2:end) = 2;
runId([false;diff(bTimes)/ratAC.sampleRate>60]) = 3;
% $$$ figure();
% $$$ hold('on');
% $$$ plot(ratAC.data(:,1,3)); % height
% $$$ plot(sqrt(sum(diff(ratFilt(:,1,1:2)).^2,3)).*ratFilt.sampleRate./10) % speed
% $$$ Lines(bTimes,[],'r');
posBeacons = unique(beacons(:,end-1:end),'rows');

FS_append_states; % Trial, rat




% $$$ figure();
% $$$ hold('on');
% $$$ plot(ratAC.data(:,1,3)); % height
% $$$ plot(sqrt(sum(diff(ratFilt(:,1,1:2)).^2,3)).*ratFilt.sampleRate./10) % speed
% $$$ Lines(rper.data(:),[],'r');


%% 10. Compute ratemaps ------------------------------------------------------------------------------------------------
% SUBSET of xyz within active stae
% p(s)
% 
%
activeState = 'theta&gper&vel-rear'; 
                       
FS_compute_placefields();                       

% LOAD spk object - holds the spike times and cluster ids ( see MTASpk )
% help MTASpk
spkw = Trial.spk.copy();
spkw.load_spk(Trial);
                       

figpath = fullfile('/storage/gravio/figures/analysis/fabian/',Trial.filebase);

sampleRate = 250;
fet = Trial.load('subject','FS04_AC');
fet.data = sq(fet(:,1,[1,2]));
fet.filter('ButFilter',4,2,'low');
[drz,fc,pfds] = compute_drz(Trial,Trial.spk.map(:,1),pfs,[],[],'Head',...
                       struct('bins',{{linspace(-500,500,200)',              ...
                                       linspace(-1000,1000,200)'}},            ...
                              'nanMaskThreshold', 0,                         ...
                              'methodNanMap',     'linear',                  ...
                              'methodRateMap',    'linear'),                 ...
                       fet,...
                       [],...
                       1.1,...
                       sampleRate ...
);

lfp = Trial.load('lfp',31);
phz = lfp.phase([6,12]);
phz.data = unwrap(phz.data);
phz.resample(sampleRate);
phz.data = mod(phz.data+pi,2*pi)-pi; 

spk = Trial.load('spk',sampleRate,activeState,[]);



create_directory(figpath);


% PLOT Placefields for each Unit
spkScale = 0.2;
hfig = figure(101);
clf(hfig);
subplotHandles = gobjects([1,4]);             % INITIALIZE empty GraphicsPlaceholder array (matlab::graphics::objectsystem)
subplotHandles = tight_subplot(2,2,0.1);      % GENERATE 1x2 subplot grid                  (labbox::Graphics)
for unit = 1:size(Trial.spk.map,1)
    ax = subplotHandles(1);
    hfig.CurrentAxes = ax;                    % SET focus on first set of axes
    cla(ax);                                  % CLEAR the contents from axes (ax)
                                              % (matlab::graphics)
    hold(ax,'on');
    plot   (pfs,unit,1,'colorbar',[],false,'flipAxesFlag',true);      % PLOT the place field of unit               (MTA::MTAApfs)
    %plot(posBeacons(:,1).*1000,posBeacons(:,2).*1000,'om','MarkerSize',10);
    plot(posBeacons(:,2).*1000,posBeacons(:,1).*1000,'om','MarkerSize',10);
    daspect(ax,[1,1,1]);                      % SET the aspect ratio                       (matlab::graph3d)
    title  (ax,num2str(unit));                % SET the title of the axes (ax)             (matlab::graph2d)

    
    ax = subplotHandles(2);
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


    ax = subplotHandles(3);
    cla(ax);
    res = spk(28);    
    plot(ax,[drz(res,unit);drz(res,unit)],[phz(res);phz(res)+2*pi],'.');
    
    %waitforbuttonpress();
    pause(0.1);
    figname = [Trial.filebase,'-vis_unit-',num2str(unit,'%04.f')];
    print(gcf, '-depsc2', fullfile(figpath,[figname,'.eps']));
    print(gcf,'-dpng',    fullfile(figpath,[figname,'.png']));
    
end


% NOTES ----------------------------------------------------------------------------------------------------------------
% graphical handles inherit the hgsetget superclass, 
% which allows you to use the get and set functions.
%
% Try using the get function on a handle from a figure/axes/plot/line 
%  EX: get(hfig)






% PLOT Placefields for each Unit
hfig = figure(102);
clf(hfig);
ax = tight_subplot(2,2,0.1,0.1,0.1)
for unit = 1:size(Trial.spk.map,1)

    mrate = max(max(plot(pfs,unit,1,'colorbar',[],false,'flipAxesFlag',true)))*1.5;
    a = 1;
    axes(ax(a));
    cla(ax(a))
    hold(ax(a),'on');
    plot(pfs,unit,1,'colorbar',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);      % PLOT the place field of unit               (MTA::MTAApfs)
    plot(ax(a),posBeacons(:,2).*1000,posBeacons(:,1).*1000,'om','MarkerSize',10);
    daspect(ax(a),[1,1,1]);                      % SET the aspect ratio                       (matlab::graph3d)
    title  (ax(a),['Unit: ',num2str(unit),' Theta-rear']);                % SET the title of the axes (ax)             (matlab::graph2d)
    
    a = 2;
    axes(ax(a));
    cla(ax(a));
    hold(ax(a),'on');
    plot   (pfa,unit,1,'colorbar',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);      % PLOT the place field of unit               (MTA::MTAApfs)
    plot(posBeacons(:,2).*1000,posBeacons(:,1).*1000,'om','MarkerSize',10);
    daspect(ax(a),[1,1,1]);                      % SET the aspect ratio                       (matlab::graph3d)
    title  (ax(a),['Unit: ',num2str(unit),' B Approach']);                % SET the title of the axes (ax)             (matlab::graph2d)

    a = 4;
    axes(ax(a));
    cla(ax(a));
    hold(ax(a),'on');
    plot   (pfd,unit,1,'colorbar',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);      % PLOT the place field of unit               (MTA::MTAApfs)
    plot(posBeacons(:,2).*1000,posBeacons(:,1).*1000,'om','MarkerSize',10);
    daspect(ax(a),[1,1,1]);                      % SET the aspect ratio                       (matlab::graph3d)
    title  (ax(a),['Unit: ',num2str(unit),' B Departure']);                % SET the title of the axes (ax)             (matlab::graph2d)
    
    a = 3;
    axes(ax(a));    
    cla(ax(a));
    hold(ax(a),'on');
    plot   (pfr,unit,1,'colorbar',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);      % PLOT the place field of unit               (MTA::MTAApfs)
    plot(posBeacons(:,2).*1000,posBeacons(:,1).*1000,'om','MarkerSize',10);
    daspect(ax(a),[1,1,1]);                      % SET the aspect ratio                       (matlab::graph3d)
    title  (ax(a),['Unit: ',num2str(unit),' Rear']);                % SET the title of the axes (ax)             (matlab::graph2d)
    
    %waitforbuttonpress();
    pause(0.1);
    figname = [Trial.filebase,'-vis_beacon_unit-',num2str(unit,'%04.f')];
% $$$     print(gcf, '-depsc2', fullfile(figpath,[figname,'.eps']));
    print(gcf,'-dpng',    fullfile(figpath,[figname,'.png']));
    
end


