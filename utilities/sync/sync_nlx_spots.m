function Session = sync_nlx_spots(Session, TTLValue, bhvSampleRate)
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

% >>> ERROR MESSAGES >>>
ERR.type = 'MTASession:create';
ERR.msg  = ['SyncPeriods is empty, check event file that ' ...
            'the TTLValues corresponding to the recording ' ...
            'trials is equal to: %s'];
% <<< ERROR MESSAGES <<<

% >>> Setup session synchronization based on the master recording system >>> --

% LOAD Session parameter file
% LOAD single channel of lfp to obtain the exact number of samples
  % NOTE - should be able to compute this from information within
  % parameter file and directory info
% ASSIGN synchronization periods in sesconds as continuous time
% CREATE Sync object for lfp
Par = LoadPar(fullfile(Session.spath, [Session.name '.xml']));            
lfp = LoadBinary(fullfile(Session.spath, [Session.name '.lfp']),1,Par.nChannels,4)';
lfp_sr = Par.lfpSampleRate;
recordSync = [0,numel(lfp)./Par.lfpSampleRate];
lfp_periods =                                    ...
    MTADepoch([],[],                             ... Sync Periods
              recordSync+[1/lfp_sr,1/lfp_sr],    ... Data
              1,                                 ... Sample Rate
              recordSync,                        ... Sync Periods
              0);%                               ... Sync Origin
clear('lfp');
Session.lfp = MTADlfp(Session.spath,             ... Path
                      [Session.name '.lfp'],     ... File Name
                      [],                        ... Data
                      lfp_sr,                    ... Sample Rate
                      lfp_periods,               ... Sync Periods
                      0);%                       ... Sync Origin
Session.lfp.filename = [Session.name '.lfp']; 
% SET session master sample rate
Session.sampleRate = Par.SampleRate;
Session.par = MTAPar(Session);
Session.load('par');
% <<< Setup session synchronization based on the master recording system <<< --

% >>> CREATE directories within current project folder >>> --------------------
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
% <<< CREATE directories within current project folder <<< --------------------

% >>> Load postion data >>> ---------------------------------------------------
% CONCATENATE xyz positions files
xy = load(fullfile(Session.path.project, Session.name,[Session.name,'.pos']));
xy = cat(2, permute(xy(:,1:2),[1,3,2]), permute(xy(:,3:4),[1,3,2]));
% <<< Load position data <<< --------------------------------------------------

% >>> GENERATE marker model from vsk file (VICON IQ ONLY) >>> -----------------
Session.model = MTAModel( {'head_back', 'head_front'}, '-mar');
% <<< GENERATE marker model from vsk file (VICON IQ ONLY) <<< -----------------

% >>> Map xyz data to sychronization events >>> -------------------------------
% LOAD recording events
events = load_nlx_evt(fullfile(Session.spath, [Session.name '.all.evt']));
% PARSE events
stopEventTTL = 'Stopping Recording';
%stopEventTTL = 'TTL Input on AcqSystem1_0 board 0 port 1 value \(0x0000\)\.'
pfirstVStart = find(events.clu==find(~cellfun(@isempty,regexp(events.labels,TTLValue))),1,'first')-1;
% REMOVE non vicon related events
events.time(1:pfirstVStart)=[];
events.clu(1:pfirstVStart)=[];
events.description(1:pfirstVStart)=[];
% GET start and stop events
vstarts = events.time(                                                      ...
    events.clu == find(                                                     ...
        ~cellfun( @isempty,                                                 ...
                  regexp( events.labels, TTLValue)                          ...
                  )                                                         ...
        )                                                                   ...
);
vstops  = events.time(                                                      ...
    find(                                                                   ...
        ismember(                                                           ...
            events.clu,                                                     ...
            find(                                                           ...
                ~cellfun( @isempty,                                         ...
                          regexp( events.labels, stopEventTTL)              ...
                          )                                                 ...
                )                                                           ...
            )                                                               ...
        )                                                                   ...
);
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
spotsSyncPeriods = [vstarts(1),vstops(end)] + [1,0] .* (1/bhvSampleRate);

% <<< Map xyz data to sychronization events <<< -------------------------------

% >>> MAP xyz data trials to lfp time periods >>> -----------------------------

% SELECT xyzData which match NLX events Pairs
xyzData = xy;
% <<< MAP xyz data trials to lfp time periods <<< -----------------------------

% >>> SETUP the xyz synchronization periods (only first and last) >>> ---------
Session.sync = MTADepoch(Session.spath,                                     ... Path
                         [Session.filebase '.sync.mat'],                    ... FileName
                         spotsSyncPeriods([1,end]),                         ... Data
                         1,                                                 ... Sample Rate
                         recordSync,                                        ... Sync Periods
                         0,                                                 ... Sync Origin
                         [],                                                ... Type
                         [],                                                ... Extension
                         [],                                                ... Name
                         'sync');%                                          ... Label
Session.sync.save(1);
% <<< SETUP the xyz synchronization periods (only first and last) <<< ---------

% >>> CONCATENATE all xyz pieces and fill gaps with zeros >>> -----------------
xy = double(xy);
% <<< CONCATENATE all xyz pieces and fill gaps with zeros <<< -----------------

% >>> INITIALIZE MTA data objects >>> -----------------------------------------

% >>> INITIALIZE MTASpk Object - holds all neuronal spiking information >>> ---
Session.spk = MTASpk;
Session.spk.parent = Session;
Session.spk.create(Session);
% <<< INITIALIZE MTASpk Object - holds all neuronal spiking information <<< ---
% >>> UPDATE the synchronization periods of the LFP object >>> ----------------
Session.lfp.sync.sync = Session.sync.copy;
Session.lfp.origin =  Session.sync.data(1);
% <<< UPDATE the synchronization periods of the LFP object <<< ----------------
% >>> INITIALIZE MTAStateCollection object, sets of behavioral periods >>> ----
Session.stc = MTAStateCollection(Session.spath,Session.filebase,'default',[],[],1);
Session.stc.updateSync(Session.sync);
Session.stc.updateOrigin(0);
% <<< INITIALIZE MTAStateCollection object, sets of behavioral periods <<< ----
% >>> INITIALIZE MTADxyz object to contain vicon data >>> ---------------------
Session.xyz = MTADxyz(Session.spath,            ... Path
                      Session.filebase,         ... File Name
                      xy,                      ... Data
                      bhvSampleRate,          ... Sample Rate
                      MTADepoch([],[],          ... Sync Periods
                                spotsSyncPeriods,    ... Data
                                1,                   ... Sample Rate
                                Session.sync.copy,   ... Sync Periods
                                0),                  ... Sync Sync Origin
                      Session.sync.data(1),     ... XYZ Sync Origin
                      Session.model);%          ... model  
Session.xyz.save;
Session.xyz.clear;
% <<< INITIALIZE MTADxyz object to contain vicon data <<< ---------------------
% >>> INITIALIZE MTADang object for inter marker spherical coordinates >>> ----
Session.ang = MTADang(Session.spath,...
                      Session.filebase,...
                      [],...
                      bhvSampleRate,...
                      Session.xyz.sync,...
                      Session.xyz.origin,...
                      Session.model);
% <<< INITIALIZE MTADang object for inter marker spherical coordinates <<< ----
% >>> INITIALIZE MTADufr object to hold unit firing rates >>> -----------------
Session.ufr = MTADufr(Session.spath,...
                      Session.filebase);
% <<< INITIALIZE MTADufr object to hold unit firing rates <<< -----------------
% >>> INITIALIZE MTADfet object to hold feature data >>> ----------------------
Session.fet = MTADfet(Session.spath,...
                      [],...
                      [],...
                      [],...
                      Session.sync.copy,...
                      Session.sync.data(1),...
                      []);                  
% <<< INITIALIZE MTADfet object to hold feature data <<< ----------------------

% <<< INITIALIZE MTA data objects <<< -----------------------------------------

%% Save the session --------------------------------------------------
Session.save;
Session.spk.clear;


%% Create other useful directories -----------------------------------
create_directory(fullfile(Session.spath,'figures'));
