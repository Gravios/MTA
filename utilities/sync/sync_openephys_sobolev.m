function Session = sync_openephys_sobolev(Session,TTLValue,mocapSampleRate)
%sync_openephys_sobolev

% required files
% [filebase].positions.traj:  ASCII???, contains a rigidbody model of the head trajcetory 
%                                       1st column is "absolute" time
%                                       columns 2-7 -> X, Z, Y, Yaw, Pitch, Roll
%
% [filebase].messages.events: ASCII???, contains time points relating trajectory time to ephys time
%                                       1st column is ephys samples
%                                       3rd column is "absolute" time

% $$$ Session.name = 'as01-20191118';
% $$$ Session.spath = '/storage/gravio/data/processed/ephys/as01/as01-20191118/';
% $$$ mocapSampleRate = 100;


%% Setup session synchronization based on the master recording system 
% LOAD Session parameter file
% LOAD single channel of lfp to obtain the exact number of samples
  % NOTE - should be able to compute this from information within
  % parameter file and directory info
% ASSIGN synchronization periods in sesconds as continuous time
% CREATE Sync object for lfp
Par = LoadPar(fullfile(Session.spath,[Session.name,'.xml']));
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



%% Load Maze XYZ data and marker model ------------------------------------

% LOAD the xyz data and angles
fid = fopen(fullfile(Session.spath,[Session.name,'.arena.traj']));
data = textscan(fid,                                                            ...
                '%n',                                                           ...
                'Delimiter',               ' ',                                 ...
                'HeaderLines',             1,                                   ...
                'CollectOutput',           false,                               ...
                'ReturnOnError',           false);
fclose(fid);
data = reshape(data{1},8,[])';
data(diff(data(:,1))==0,:) = [];
arenaTrajData = data(:,[2:end]);
clear('data');




%% Load xyz data and marker model ------------------------------------

% LOAD the xyz data and angles
fid = fopen(fullfile(Session.spath,[Session.name,'.positions.traj']));
data = textscan(fid,                                                            ...
                '%n',                                                           ...
                'Delimiter',               ' ',                                 ...
                'HeaderLines',             1,                                   ...
                'CollectOutput',           false,                               ...
                'ReturnOnError',           false);
fclose(fid);
data = reshape(data{1},7,[])';
data(diff(data(:,1))==0,:) = [];
tsTraj = data(:,1);
trajData = data(:,[2:7]);
clear('data');



% GET the events, which map mocap time to ephys time.
buffer    = fileread(fullfile(Session.spath,[Session.name,'.messages.events'])) ;
data      = strsplit(buffer,'\n');
startTime = str2num(data{1}(1:find(data{4}==' ',1,'first')));
data(1:4) = [];
data = regexprep(data,'\w+:','');
data(cellfun(@isempty,data)) = [];
data = cf(@(d) strsplit(d(1:end-1)),data);
data = reshape(cell2mat(cf(@(d) cell2mat(cf(@(e) str2num(e), d)), data)),5,[])';
evtxyz = data(:,3:end);
tsEphySync = (data(:,1)-startTime)./Par.SampleRate; % CONVERT from oepenephys sample rate to seconds
tsTrajSync = data(:,2);
clear('data');

tsTraj = tsTraj-tsTrajSync(1)+tsEphySync(1);
%arenaTrajTS = arenaTrajTS-tsTrajSync(1)+tsEphySync(1);
tsTrajSync = tsTrajSync-tsTrajSync(1)+tsEphySync(1);



% $$$ tsTraj = find(tsTraj>tsTrajSync(1),1,'first') 
% $$$ p = polyfit(tsTrajSync-tsTrajSync(1),tsEphySync,1);
% $$$ tsTrajInterp = polyval(p,tsTraj-tsTrajSync(1));

% $$$ trajData(tsTrajInterp<=0,:) =[];
% $$$ tsTrajInterp(tsTrajInterp<=0) = [];

tsFinal = [tsTraj(1):(1/mocapSampleRate):tsTraj(end)]';

trajDataInterp = interp1(tsTraj,trajData,tsFinal,'linear');
mazeDataInterp = interp1(tsTraj,arenaTrajData,tsFinal,'linear');

% CONVERT to mutli marker rigid body model
markers = {'hcom'};
xyzData = {trajDataInterp(:,[1,3,2]).*1000};
mazeData = {mazeDataInterp(:,[1,3,2]).*1000};


% GENERATE marker model from vsk file (VICON IQ ONLY)
model = MTAModel(markers,'-mar');

Session.model = model;


mocapSyncPeriods = tsFinal([1,end])';


% SETUP the xyz synchronization periods (only first and last)
Session.sync = MTADepoch(Session.spath,                      ... Path
                         [Session.filebase '.sync.mat'],     ... FileName
                         mocapSyncPeriods([1,end]),          ... Data
                         1,                                  ... Sample Rate
                         recordSync,                         ... Sync Periods
                         0,                                  ... Sync Origin
                         [],                                 ... Type
                         [],                                 ... Extension
                         [],                                 ... Name
                         'sync');%                               Label
Session.sync.save(1);


%% NEEDED IN THE FUTUERE MAYBE
% $$$ % CONCATENATE all xyz pieces and fill gaps with zeros
% $$$ nViconTrials = length(xyzData);
% $$$ xyzLengths = cellfun(@length,xyzData);
% $$$ xyz = zeros([ceil(diff(mocapSyncPeriods([1,end]))*viconSampleRate),size(xyzData{1},2),size(xyzData{1},3)]);
% $$$ syncXyzStart = round((mocapSyncPeriods(1:length(mocapSyncPeriods))-mocapSyncPeriods(1))*viconSampleRate+1);
% $$$ for s=1:nViconTrials,
% $$$     xyzseg = xyzData{s};
% $$$     xyzseg(xyzseg==0)=eps;
% $$$     xyz(syncXyzStart(s):syncXyzStart(s)+xyzLengths(s)-1,:,:) = xyzseg;
% $$$ end
xyz = permute(cat(1,xyzData{:}),[1,3,2]);
mazeXYZ = permute(cat(1,mazeData{:}),[1,3,2]);


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
                      mocapSampleRate,          ... Sample Rate
                      MTADepoch([],[],          ... Sync Periods
                                mocapSyncPeriods,    ... Data
                                1,                   ... Sample Rate
                                Session.sync.copy,   ... Sync Periods
                                0),                  ... Sync Sync Origin
                      Session.sync.data(1),     ... XYZ Sync Origin
                      Session.model);%          ... model  
Session.xyz.save;
Session.xyz.clear;

% INITIALIZE MTADxyz object to contain maze position data
% INITIALIZE MTADxyz object to contain vicon data
Session.arena = MTADxyz(Session.spath,            ... Path
                      Session.filebase,         ... File Name
                      mazeXYZ,                      ... Data
                      mocapSampleRate,          ... Sample Rate
                      MTADepoch([],[],          ... Sync Periods
                                mocapSyncPeriods,    ... Data
                                1,                   ... Sample Rate
                                Session.sync.copy,   ... Sync Periods
                                0),                  ... Sync Sync Origin
                      Session.sync.data(1),     ... XYZ Sync Origin
                           [],                       ... model
                           [],                       ... type
                           [],                       ... ext
                           [],                       ... name
                           'maze'                    ... label
);
Session.arena.save();
Session.arena.clear();



% INITIALIZE MTADang object for inter marker spherical coordinates
Session.ang = MTADang(Session.spath,...
                      Session.filebase,...
                      [],...
                      mocapSampleRate,...
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
%create_directory(fullfile(Session.spath,'figures'));
