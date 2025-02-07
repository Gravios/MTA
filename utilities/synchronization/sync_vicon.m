function Session = sync_vicon(Session, viconSampleRate)
%
% Populate a Session with xyz data synchronized to
% electrophysiological data based on an event file
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
Session.sampleRate = viconSampleRate;


%% Setup session folder in project directory -------------------------

% CREATE directories within current project folder
if ~exist(Session.spath,'dir')
    if ~exist(Session.spath,'dir') && ...
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

% CONCATENATE xyz positions files
if isempty(viconSampleRate),
    [xyzData, markers, viconSampleRate] = concatenate_vicon_files(Session);            
    if isempty(viconSampleRate),
        error('MTA:utilities:syncViconNlx:emptySampleRate');
    end
else
    [xyzData, markers] = concatenate_vicon_files(Session);
end
% GENERATE marker model from vsk file (VICON IQ ONLY)
vsk_path = fullfile(Session.spath, Session.maze.name, [Session.name '-' Session.maze.name '.vsk']);
if exist(vsk_path,'file'),
    model = MTAModel(vsk_path,'-vsk');
else
    warning(['VSK file associated with this session was ' ...
             'not found. \nCreating a general marker model.']);
    model = MTAModel(markers,'-mar');
end
Session.model = model;




%% MAP xyz data trials to lfp time periods ---------------------------

viconSyncPeriods = cellfun(@length,xyzData);
viconSyncPeriods = [1,viconSyncPeriods+1;
                    viconSyncPeriods,1];
viconSyncPeriods(:,end) = [];
viconSyncPeriods = viconSyncPeriods';

% ASSIGN xyz data to start/stop periods
% SETUP the xyz synchronization periods (only first and last)
Session.sync = MTADepoch(Session.spath,                      ... Path
                         [Session.filebase '.sync.mat'],     ... FileName
                         viconSyncPeriods([1,end]),          ... Data
                         viconSampleRate,                    ... Sample Rate
                         viconSyncPeriods([1,end]),          ... Sync Periods
                         1,                                  ... Sync Origin
                         [],                                 ... Type
                         [],                                 ... Extension
                         [],                                 ... Name
                         'sync');%                               Label
Session.sync.save(1);


% CONCATENATE all xyz pieces and fill gaps with zeros
nViconTrials = length(xyzData);
xyzLengths = cellfun(@length,xyzData);
xyz = zeros([viconSyncPeriods([end]), size(xyzData{1},2), size(xyzData{1},3)]);
syncXyzStart = viconSyncPeriods(:,1);
for s=1:nViconTrials,
    xyz(syncXyzStart(s):syncXyzStart(s)+xyzLengths(s)-1,:,:) = xyzData{s};
end
xyz = double(xyz);

%% Initialize MTA data objects ---------------------------------------

% INITIALIZE MTADxyz object to contain vicon data
Session.xyz = MTADxyz(Session.spath,            ... Path
                      Session.filebase,         ... File Name
                      xyz,                      ... Data
                      viconSampleRate,          ... Sample Rate
                      MTADepoch([],[],          ... Sync Periods
                                viconSyncPeriods,... Data
                                viconSampleRate, ... Sample Rate
                                Session.sync.copy,   ... Sync Periods
                                1),                  ... Sync Sync Origin
                      Session.sync.data(1),     ... XYZ Sync Origin
                      Session.model);%          ... model  
Session.xyz.save;
Session.xyz.clear;

% INITIALIZE MTADang object for inter marker spherical coordinates
Session.ang = MTADang(Session.spath,...
                      Session.filebase,...
                      [],...
                      viconSampleRate,...
                      Session.xyz.sync,...
                      Session.xyz.origin,...
                      Session.model);

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


%% Create other useful directories -----------------------------------
create_directory(fullfile(Session.spath,'figures'));
