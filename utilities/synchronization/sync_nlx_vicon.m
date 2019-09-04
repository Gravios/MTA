function Session = sync_nlx_vicon(Session,TTLValue,viconSampleRate)
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
vsk_path = fullfile(Session.spath, [Session.name '-' Session.maze.name '.vsk']);
if exist(vsk_path,'file'),
    model = MTAModel(vsk_path,'-vsk');
else
    warning(['VSK file associated with this session was ' ...
             'not found. \nCreating a general marker model.']);
    model = MTAModel(markers,'-mar');
end
Session.model = model;


%% Map xyz data to sychronization events -----------------------------






% LOAD recording events
events = LoadEvents(fullfile(Session.spath, [Session.name '.all.evt']));
% PARSE events
stopEventTTL = '0x0000';
%stopEventTTL = 'TTL Input on AcqSystem1_0 board 0 port 1 value \(0x0000\)\.'
pfirstVStart = find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))),1,'first')-1;
% REMOVE non vicon related events
events.time(1:pfirstVStart)=[];
events.Clu(1:pfirstVStart)=[];
events.description(1:pfirstVStart)=[];
% GET start and stop events
vstarts = events.time(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))));
vstops  = events.time(find(ismember(events.Clu,find(~cellfun(@isempty,regexp(events.Labels,stopEventTTL))))));
if numel(vstops) ~= numel(vstarts),
    s = 1;
    while (numel(vstops) ~= numel(vstarts)) &  (s < numel(vstarts)),
        if (vstops(s) - vstarts(s)) < 0
            vstops(s) = [];
        else
            s = s + 1;
        end
        
    end
    vstops = vstops(1:numel(vstarts));
end
% CONCATENATE start and stop events 
viconPeriods = [vstarts,vstops];


%% MAP xyz data trials to lfp time periods ---------------------------

% ASSIGN xyz data to start/stop periods
viconSyncPeriods=[];
xyzDataInd = [];
for i=1:numel(xyzData),
    if ~isempty(xyzData{i})
        for j=1:size(viconPeriods,1),
            if abs(diff(viconPeriods(j,:))*10-size(xyzData{i},1)/viconSampleRate*10)<0.2,
                % ???Shift time back so index corresponds to
                % viconSyncPeriods(1) :-1/viconSampleRate???
                viconSyncPeriods(end+1,:) = viconPeriods(j,:)-1/viconSampleRate;
                viconPeriods(viconPeriods(j,1)>=viconPeriods(:,1),:) = [];
                xyzDataInd(end+1) = i;
                break;
            end
        end
    end
end
% THROW error if no NLX events match xyzData
assert(~isempty(viconSyncPeriods),ERR.type,ERR.msg,TTLValue);
% SELECT xyzData which match NLX events Pairs
xyzData = xyzData(xyzDataInd);
% SETUP the xyz synchronization periods (only first and last)
Session.sync = MTADepoch(Session.spath,                      ... Path
                         [Session.filebase '.sync.mat'],     ... FileName
                         viconSyncPeriods([1,end]),          ... Data
                         1,                                  ... Sample Rate
                         recordSync,                         ... Sync Periods
                         0,                                  ... Sync Origin
                         [],                                 ... Type
                         [],                                 ... Extension
                         [],                                 ... Name
                         'sync');%                               Label
Session.sync.save(1);


% CONCATENATE all xyz pieces and fill gaps with zeros
nViconTrials = length(xyzData);
xyzLengths = cellfun(@length,xyzData);
xyz = zeros([ceil(diff(viconSyncPeriods([1,end]))*viconSampleRate),size(xyzData{1},2),size(xyzData{1},3)]);
syncXyzStart = round((viconSyncPeriods(1:length(viconSyncPeriods))-viconSyncPeriods(1))*viconSampleRate+1);
for s=1:nViconTrials,
    xyzseg = xyzData{s};
    xyzseg(xyzseg==0)=eps;
    xyz(syncXyzStart(s):syncXyzStart(s)+xyzLengths(s)-1,:,:) = xyzseg;
end
xyz = double(xyz);

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
                      xyz,                      ... Data
                      viconSampleRate,          ... Sample Rate
                      MTADepoch([],[],          ... Sync Periods
                                viconSyncPeriods,    ... Data
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
                      viconSampleRate,...
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
