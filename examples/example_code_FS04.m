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

% METADATA  ----------------------------------------------------------------------------------
% SET Session info
meta.sessionBase = 'FS04';
meta.sessionName = 'FS04-20210326a';
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


% CONSTRUCT paths
MTA_DATA_PATH = getenv('MTA_DATA_PATH');
meta.path.raw.ephys       = create_directory(fullfile(MTA_DATA_PATH,'raw/ephys/',      meta.sessionBase,meta.sessionName));
meta.path.raw.xyz         = create_directory(fullfile(MTA_DATA_PATH,'raw/xyz/',        meta.sessionBase,meta.sessionName,meta.mazeName));
meta.path.processed.xyz   = create_directory(fullfile(MTA_DATA_PATH,'processed/xyz/',  meta.sessionBase,meta.sessionName,meta.mazeName));
meta.path.processed.ephys = create_directory(fullfile(MTA_DATA_PATH,'processed/ephys/',meta.sessionBase,meta.sessionName));


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
cf(@(csv) movefile(fullfile(meta.path.processed.ephys,csv),meta.path.raw.xyz), meta.csv(~isempty(meta.csv)));

% NOTE : not really a TTLValue
%
meta.TTLValue = convert_srs_to_TTLValue(meta);    
% EXPECTED output
% meta.TTLValue = { ...                                record1
%                   {[46029977]},...                     record sample count
%                   {'FS04-20210322a.take_0001.mat'},... file with extracted rigidbodies
%                   {'FS04'}...                          session's primary subject
% }



%% 4. Preprocess CSV ---------------------------------------------------------------------------------------------------
%
% EXTRACTS rigid bodies from the csv file generated by motive
%
parse_rbo_from_csv(meta);



%% 5. Build Sessions ---------------------------------------------------------------------------------------------------
link_session( meta.sessionName, meta.dPaths)
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
% writing out the session name 
Session = MTASession.validate('FS04-20210322a.vrr.all');

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
Trial = MTATrial.validate('FS04-20210322a.vrr.all');



%% 8. Load Subjects ----------------------------------------------------------------------------------------------------
% problem: getting the modified subject objects passed to the pfs computation
% problem: loading Objects

rat   = Trial.load('subject','FS04');
Arena = Trial.load('subject','Arena');



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
ratRC.label = [ratRC.label,'_RC'];
ratAC.update_filename([Trial.filebase,'.rbo.',rat.name,'_RC.s.mat']);
ratRC.save();


% CONVERT the xyz coordinates to Arena frame of refrences centered on the Arena
ratAC = copy(rat);
% TRANSLATE 
ratAC.data(:,1,1:3) = ratAC(:,'Head',1:3) - [Arena(:,'Arena',1:3);zeros([1,1,3])];
% ROTATE 
o = 0.0349065850398866;
ratAC.data(:,1,1:2) = multiprod([cos(o),-sin(o);sin(o),cos(o)],sq(ratAC.data(:,1,1:2)),[1,2],[2]);
ratAC.label = [ratAC.label,'_AC'];
ratAC.update_filename([Trial.filebase,'.rbo.',rat.name,'_AC.s.mat']);
ratAC.save();




%% 9. Label Behaviors --------------------------------------------------------------------------------------------------
% LABEL 

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



%ArenaConfigA
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

arenaConfB = ThreshCross(double(Arena(:,'Arena',2)<700),0.5,10);
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   arenaConfB,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'ArenaConfB','B');
Trial.stc.states{end}.save(1);

arenaConfC = ThreshCross(double(Arena(:,'Arena',2)<974&Arena(:,'Arena',2)>700),0.5,10);
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   arenaConfC,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'ArenaConfC','C');
Trial.stc.states{end}.save(1);







%% 10. Compute ratemaps ------------------------------------------------------------------------------------------------
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





% CORRECT Yaw for rigid bodies (SIMPLE) ------------------------------------------------------------

% GET head yaw from rat
ang = quaternion2rad(sq(rat(:,1,5:8)));
figure,plot(ang(:,1))

% GET trajectory yaw from rat
trjVect = sq(circshift(ratFilt(:,'Head',[1,2]),-20)-circshift(ratFilt(:,'Head',[1,2]),20)); % vector
trjVect = bsxfun(@rdivide,trjVect,sqrt(sum(trjVect.^2,2)));                                 % basis
trjYaw  = atan2(trjVect(:,1),trjVect(:,2));                                                 % yaw
figure();plot(trjYaw)

% GET XY head speed from rat
ratFilt = copy(rat);
ratFilt.filter('ButFilter',4,1.5,'low');
spdHead = sqrt(sum(sq(ratFilt(:,'Head',[1,2])-circshift(ratFilt(:,'Head',[1,2]),-1)).^2,2)).*ratFilt.sampleRate./10;
figure();plot(spdHead)

figure();hist(circ_dist(ang(spdHead>20,1),trjYaw(spdHead>20,1)),linspace(-pi,pi,100))


