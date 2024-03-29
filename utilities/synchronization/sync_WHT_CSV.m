function Session = sync_WHT_CSV(Session,recordLengths,xyzSampleRate)
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
%dat = LoadBinary(fullfile(Session.spath, [Session.name '.dat']),1,Par.nChannels,1)';
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
create_directory(Session.spath);        


%% Load xyz data and marker model ------------------------------------

subjectRecords = load_take_files(Session,recordLengths);


% GET first non emtpy record for reference
rec = find(~cellfun(@isempty,subjectRecords));
% ASSERT that the sample rate given matches the sample rate of the record
assert(xyzSampleRate == subjectRecords{rec(1)}(1).sampleRate);
% GET a list of the subject names
subjectNameList = {subjectRecords{rec(1)}.name};
% SET the identity of the primary subject
primarySubject = recordLengths{rec(1)}{3};

rboData = {};
msg_pad = ['[INFO] MTA:utilities:synchronization:sync_WHT_CSV - ' ...
           'subject: %s |  Record: %d | Mode: %s'];
for n = 1:numel(subjectNameList),
    for s = 1:numel(subjectRecords)
        if isempty(subjectRecords{s})
            disp(sprintf(msg_pad,subjectNameList{n},s,'no csv pad zeros'))
            % GENERATE zero padded 
            rboData{n}{s} = zeros([round(recordLengths{s}{1}./Session.sampleRate.*xyzSampleRate), ...
                                    size(subjectRecords{rec(1)}(n).data,2),...
                                    size(subjectRecords{rec(1)}(n).data,3)])+eps;
        else
            % if csv is larger than ephys segment -> truncate csv data
            if size(subjectRecords{s}(n).data,1) >= round(recordLengths{s}{1}./Session.sampleRate.*xyzSampleRate)+1
                disp(sprintf(msg_pad,subjectNameList{n},s,'truncate csv'))
                rboData{n}{s} = subjectRecords{s}(n).data(1:round(recordLengths{s}{1}./Session.sampleRate.*xyzSampleRate),:,:)+ eps;
            else % pad csv data to match ephys segment
                disp(sprintf(msg_pad,subjectNameList{n},s,['pad zeros to end of csv']))
                rboData{n}{s} = cat(1,...
                                subjectRecords{s}(n).data,...
                                zeros([round(recordLengths{s}{1}./Session.sampleRate.*xyzSampleRate)-size(subjectRecords{s}(n).data,1),...
                                       size(subjectRecords{s}(n).data,2),...
                                       size(subjectRecords{s}(n).data,3)]))+ eps;
            end
            %rboData{n}{s} = rboData{n}{s} + eps;
        end
    end
    rboLabels{n} = subjectRecords{rec(1)}(n).rboLabels;
end



cumulativeRecordLength = [0,cumsum(cellfun(@length,rboData{1}))]./xyzSampleRate;

xyzSyncPeriods = [];
%for ind = find(~cellfun(@isempty,subjectRecords))
for ind = 1:numel(subjectRecords)
    xyzSyncPeriods(end+1,:) = cumulativeRecordLength([ind,ind+1])+[1/xyzSampleRate,0];
end
xyzSyncPeriods(end) = xyzSyncPeriods(end) - 1/xyzSampleRate;


%markers = {'head_frontL','head_frontR','head_left','head_back','head_right'};

           
% GENERATE marker model
%s.model = MTAModel(markers,'-mar');


%% Map xyz data to sychronization events -----------------------------


% SETUP the xyz synchronization periods (only first and last)
Session.sync = MTADepoch(Session.spath,                      ... Path
                         [Session.filebase '.sync.mat'],     ... FileName
                         recordSync,                         ... Data
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



xyzData = [];
% INITIALIZE MTADxyz object to contain vicon data
Session.xyz = MTADxyz(Session.spath,           ... Path
                      Session.filebase,        ... File Name
                      xyzData,                 ... Data
                      xyzSampleRate,           ... Sample Rate
                      MTADepoch([],[],         ... Sync Periods
                                xyzSyncPeriods,      ... Sync Data
                                1,                   ... Sample Rate
                                Session.sync.copy,   ... Sync Periods
                                0),                  ... Sync Sync Origin
                      Session.sync.data(1),     ... XYZ Sync Origin %
                      []);%          ... model  
Session.xyz.save;
Session.xyz.clear;


for n = 1:numel(rboData)
    model = MTAModel(rboLabels{n},'-mar');
    data = cat(1,rboData{n}{:});
% INITIALIZE MTADxyz object to contain vicon data
    rbo = MTADrbo(Session.spath,                ... Path
                  Session.filebase,             ... File Name
                  data,                         ... Data
                  xyzSampleRate,                ... Sample Rate
                  MTADepoch([],[],              ... Sync Periods
                            xyzSyncPeriods,... Sync Data
                            1,                            ... Sync Sample Rate
                            Session.sync.copy,            ... Sync Periods
                            0),                           ... Sync Sync Origin
                  Session.sync.data(1),         ... rbo Sync Origin 
                  model,                        ... rbo model  
                  'TimeSeries',                 ... data type
                  'rbo',                        ... extension
                  subjectNameList{n},           ... name
                  subjectNameList{n},           ... label
                  's'                           ... key
    );% rbo
    rbo.save;
    rbo.clear;
    switch subjectNameList{n}
      case primarySubject
        Session.subject = copy(rbo);
        Session.model = Session.subject.model;
      case 'Arena'
        Session.arena = copy(rbo);
    end
end    

% INITIALIZE MTADang object for inter marker spherical coordinates
Session.ang = MTADang(Session.spath,...
                      Session.filebase,...
                      [],...
                      xyzSampleRate,...
                      Session.xyz.sync,...
                      Session.xyz.origin,...
                      []);

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
