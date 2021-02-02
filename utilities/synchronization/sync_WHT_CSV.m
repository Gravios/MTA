function Session = sync_WHT_CSV(Session,~,xyzSampleRate)
%
% Populate a Session with xyz data synchronized to
% electrophysiological data based on an event file
%
% sync_WHT_CSV assumes the ephys and mocap are activated simultaneously
% but halted asynchronously with ephys preceeding the mocap. 
%
% suports a single subject.
%   %NOTE - multi subject suport will be implemented at a later date
%
% note - the event file must be located in the nlx folder
%
%

ERR.type = 'MTASession:create';
ERR.msg  = ['SyncPeriods is empty, check event file that ' ...
            'the TTLValues corresponding to the recording ' ...
            'trials is equal to: %s'];


%% Setup session synchronization based on the master recording system 

% LOAD Session parameter file
% LOAD single channel of lfp to obtain the exact number of samples
  % NOTE - should be able to compute this from information within
  % parameter file and directory info
% ASSIGN synchronization periods in sesconds as continuous time
% CREATE Sync object for lfp
Par = LoadPar(fullfile(Session.spath, [Session.name '.xml']));            
lfp = LoadBinary(fullfile(Session.spath, [Session.name '.lfp']),1,Par.nChannels,4)';
recordSync = [0,numel(lfp)./Par.lfpSampleRate];
clear('lfp');
Session.lfp = MTADlfp(Session.spath,             ... Path
                      [Session.name '.lfp'],     ... File Name
                      [],                        ... Data
                      Par.lfpSampleRate,         ... Sample Rate
                      MTADepoch([],[],           ... Sync Periods
                                recordSync+[1/Par.lfpSampleRate,1/Par.lfpSampleRate],... Data
                                1,                   ... Sample Rate
                                recordSync,          ... Sync Periods
                                0),                  ... Sync Origin
                      0);%                           Sync Origin
Session.lfp.filename = [Session.name '.lfp']; % NOTE - May be redundant
% SET session master sample rate
Session.sampleRate = Par.SampleRate;


%% Setup session folder in project directory -------------------------

% CREATE directories within current project folder
if ~exist(Session.spath,'dir')
    if ~exist(Session.spath,'dir') && ...
        exist(fullfile(Session.path.project,'nlx',Session.name),'dir') && ...
        exist(fullfile(Session.path.project,'xyz',Session.name),'dir')
        create_directory(Session.spath);        
    else
        e.message    = ['Session: ' Session.name ', cannot be found or does not exist. ' ...
                       'Check paths or see the README for the correct directory structure ' ...
                       'for the MTA toolbox.'];
        e.identifier = 'MTA:utilities:syncViconNlx:Session404';
        error(e);
    end
end


%% Load xyz data and marker model ------------------------------------


markers = {'head_frontL','head_frontR','head_left','head_back','head_right'};

ds = load(fullfile(Session.spath,Session.maze.name,[Session.name,'.pos.data.mat']));
% $$$ figure,
% $$$ plot3(ds.data(10000,:,1),ds.data(10000,:,2),ds.data(10000,:,3),'-+');
% $$$ hold('on');
% $$$ plot3(ds.data(10000,1,1),ds.data(10000,1,2),ds.data(10000,1,3),'or');

xyzData = ds.data;
           
% GENERATE marker model
s.model = MTAModel(markers,'-mar');;


%% Map xyz data to sychronization events -----------------------------


% CONCATENATE start and stop events 
xyzSyncPeriods = [0,ds.timestamps(end)];


% SETUP the xyz synchronization periods (only first and last)
Session.sync = MTADepoch(Session.spath,                      ... Path
                         [Session.filebase '.sync.mat'],     ... FileName
                         xyzSyncPeriods([1,end]),          ... Data
                         1,                                  ... Sample Rate
                         recordSync,                         ... Sync Periods
                         0,                                  ... Sync Origin
                         [],                                 ... Type
                         [],                                 ... Extension
                         [],                                 ... Name
                         'sync');%                               Label
Session.sync.save(1);


%% Initialize MTA data objects ---------------------------------------

% INITIALIZE MTASpk Object - holds all neuronal spiking information
Session.spk = MTASpk;
Session.spk.create(Session);

% UPDATE the synchronization periods of the LFP object
Session.lfp.sync.sync = Session.sync.copy;
Session.lfp.origin =  Session.sync.data(1);

% INITIALIZE MTAStateCollection object, which holds all behavioral sets of periods
Session.stc = MTAStateCollection(Session.spath,Session.filebase,'default',[],[],1);
Session.stc.updateSync(Session.sync);
Session.stc.updateOrigin(0);

% INITIALIZE MTADxyz object to contain vicon data
Session.xyz = MTADxyz(Session.spath,            ... Path
                      Session.filebase,         ... File Name
                      xyzData,                   ... Data
                      xyzSampleRate,          ... Sample Rate
                      MTADepoch([],[],          ... Sync Periods
                                xyzSyncPeriods,    ... Data
                                1,                   ... Sample Rate
                                Session.sync.copy,   ... Sync Periods
                                0),                  ... Sync Sync Origin
                      Session.sync.data(1),     ... XYZ Sync Origin
                      Session.model);%          ... model  
Session.xyz.save;
Session.xyz.clear;

% INITIALIZE MTADang object for inter marker spherical coordinates
Session.ang = MTADang(Session.spath,...
                      Session.filebase,...
                      [],...
                      xyzSampleRate,...
                      Session.xyz.sync,...
                      Session.xyz.origin,...
                      Session.model);

% INITIALIZE MTADufr object to hold unit firing rates computed from MTASpk object
Session.ufr = MTADufr(Session.spath,...
                      Session.filebase);

% INITIALIZE MTADfet object to hold feature data
Session.fet = MTADfet(Session.spath,...
                      [],...
                      [],...
                      [],...
                      Session.sync.copy,...
                      Session.sync.data(1),...
                      []);                  


%% Save the session --------------------------------------------------
Session.save;
Session.spk.clear;


%% Create other useful directories -----------------------------------
create_directory(fullfile(Session.spath,'figures'));
