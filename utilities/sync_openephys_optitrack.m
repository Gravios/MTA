function Session = sync_openephys_optitrack(Session,recSyncPulseChan,camSyncPulseChan)

% find following values for synchronization
% syncPeriods
% mocapSampleRate

% REQVARS for sync_openephys_optitrack
% recSyncPulseChan = 17;
% mocapSampleRate = 18;

ERR.type = 'MTASession:create';            
ERR.msg  = ['SyncPeriods is empty, check recording sync channel: %s'];

    
%recordingSyncPulseChannel
Par = LoadPar(fullfile(Session.spath,[Session.name,'.xml']));

mocapSyncPulse = LoadBinary(fullfile(Session.spath,[Session.name,'.lfp']),recSyncPulseChan,Par.nChannels,4)';

mocapSyncPeriods = ThreshCross(abs((mocapSyncPulse-mean(mocapSyncPulse))./max(mocapSyncPulse(:,1))),...
                          0.5,... threshold
                          0); %   Minumum Interval
mocapSyncPeriods = reshape(round(mean(mocapSyncPeriods,2)),2,[])';


camSyncPulse = LoadBinary(fullfile(Session.spath,[Session.name,'.lfp']),...
                          camSyncPulseChan,Par.nChannels,4)';
camSyncInds = LocalMinima(-camSyncPulse,4,-1e4);



mocapSampleRate = 1/(mean(diff( camSyncInds ))./Par.lfpSampleRate);



if exist('Par','var'),
    %% Load single channel of lfp to check the exact number of samples
    recordSync = [0,numel(mocapSyncPulse)./Par.lfpSampleRate];
    lfpSyncPeriods = MTADepoch([],[],recordSync+[1/Par.lfpSampleRate,1/Par.lfpSampleRate],1,recordSync,0);
    Session.lfp = MTADlfp(Session.spath,[Session.name '.lfp'],[],Par.lfpSampleRate,lfpSyncPeriods,0);    
    Session.lfp.filename = [Session.name '.lfp'];
    Session.sampleRate = Par.SampleRate;
end

if ~exist(Session.spath,'dir')
    if ~exist(Session.spath,'dir') && ...
        exist(fullfile(Session.path.data,'nlx',Session.name),'dir') && ...
        exist(fullfile(Session.path.data,'xyz',Session.name),'dir')
        mkdir(Session.spath);
        
    else
        e.message    = ['Session: ' Session.name ', cannot be found or does not exist. ' ...
                       'Check paths or see the README for the correct directory structure for the MTA toolbox.'];
        e.identifier = 'MTA:utilities:syncViconNlx:SessionNotFound';
        error(e);
    end
end




%% Organize the Sessions trials into a cell array
%% Return the names of the markers
if isempty(mocapSampleRate),
    [xyzData, markers, mocapSampleRate] = concatViconFiles(Session);            
    if isempty(mocapSampleRate),
        error('MTA:utilities:syncViconNlx:emptySampleRate');
    end
else
    [xyzData, markers] = concatViconFiles(Session);            
end

%% Load VSK if possible 
vsk_path = fullfile(Session.spath, [Session.name '-' Session.maze.name '.vsk']);
if exist(vsk_path,'file'),
    model = MTAModel(vsk_path,'-vsk');
else
    warning(['VSK file associated with this session was ' ...
             'not found. \nCreating a general marker model.']);
    model = MTAModel(markers,'-mar');
end
Session.model = model;

%assert(exist([Session.spath Session.name '.all.evt'],'file'))
% Load events
tsRanges = mocapSyncPeriods/Par.lfpSampleRate;


%% Assign xyz data to nlx event epochs
syncPeriods=[];
xyzDataInd = [];
for i=1:numel(xyzData),
    if ~isempty(xyzData{i})
        for j=1:size(tsRanges,1),
            %diff(tsRanges(j,:))-size(xyzData{i},1)/mocapSampleRate
            if abs(diff(tsRanges(j,:))*10-size(xyzData{i},1)/mocapSampleRate*10)<0.2,
                % ???Shift time back so index corresponds to
                % syncPeriods(1) :-1/mocapSampleRate???
                syncPeriods(end+1,:) = tsRanges(j,:)-1/mocapSampleRate;

                tsRanges(tsRanges(j,1)>=tsRanges(:,1),:) = [];
                xyzDataInd(end+1) = i;
                break;
            end
        end
    end
end


%% Error if no NLX events match xyzData
assert(~isempty(syncPeriods),ERR.type,ERR.msg,recSyncPulseChan);

%% Select xyzData which match NLX events Pairs
xyzData = xyzData(xyzDataInd);


%% Setup the Session Synchronization Periods
% syncViconNlx - Sessions are synchronized to the period
% between the start and end  of vicon recording
Session.sync = MTADepoch(Session.spath,[Session.filebase '.sync.mat'],syncPeriods([1,end]),1,recordSync,0,[],[],[],'sync');
Session.sync.save(1);


%% Concatenate all xyz pieces and fill gaps with zeros
nViconTrials = length(xyzData);
xyzLengths = cellfun(@length,xyzData);
xyz = zeros([ceil(diff(syncPeriods([1,end]))*mocapSampleRate),size(xyzData{1},2),size(xyzData{1},3)]);
syncXyzStart = round((syncPeriods(1:length(syncPeriods))-syncPeriods(1))*mocapSampleRate+1);
%syncShift = [0;floor(diff([syncPeriods(1:end-1,2) syncPeriods(2:end,1)],1,2)*mocapSampleRate)];
for s=1:nViconTrials,
% $$$     if syncShift(s),
% $$$         xyz = cat(1,xyz,zeros(syncShift(s),size(xyzData{s},2),size(xyzData{s},3)));
% $$$     end
xyzseg = xyzData{s};
xyzseg(xyzseg==0)=eps;
xyz(syncXyzStart(s):syncXyzStart(s)+xyzLengths(s)-1,:,:) = xyzseg;
% $$$     [size(xyz,1)/mocapSampleRate+syncPeriods.data(1)+1/mocapSampleRate,syncPeriods.data(s)]
% $$$     diff([size(xyz,1)/mocapSampleRate+syncPeriods.data(1)+1/mocapSampleRate,syncPeriods.data(s)])
% $$$     xyz = cat(1,xyz,xyzData{s});
end
%xyz(syncXyzStart(s)+xyzLengths(s):end,:,:)=1;
xyz = double(xyz);


%% MTASpk Object - holds all neuronal spiking information
Session.spk = MTASpk;
Session.spk.create(Session);

%% Update the synchronization periods of the LFP object
Session.lfp.sync.sync = Session.sync.copy;
Session.lfp.origin =  Session.sync.data(1);

%% MTAStateCollection object holds all behavioral sets of periods
Session.stc = MTAStateCollection(Session.spath,Session.filebase,'default',[],[],1);
Session.stc.updateSync(Session.sync);
Session.stc.updateOrigin(0);


syncPeriods = MTADepoch([],[],syncPeriods,1,Session.sync.copy,0);

Session.xyz = MTADxyz(Session.spath,Session.filebase,xyz,mocapSampleRate,...
                      syncPeriods,Session.sync.data(1),Session.model);                  
Session.xyz.save;
Session.xyz.clear;

Session.ang = MTADang(Session.spath,...
                      Session.filebase,...
                      [],...
                      mocapSampleRate,...
                      Session.xyz.sync,...
                      Session.xyz.origin,...
                      Session.model);

Session.ufr = MTADufr(Session.spath,Session.filebase);

Session.fet = MTADfet(Session.spath,...
                      [],...
                      [],...
                      [],...
                      Session.sync.copy,...
                      Session.sync.data(1),...
                      []);                  


Session.save;
Session.spk.clear;

end





% $$$ Par = LoadPar('/storage/javier/Raw_data/ephys/JM11/2016-11-25_18-34-06/processed/2016-11-25_18-34-06.xml');;
% $$$ lfp = LoadBinary('/storage/javier/Raw_data/ephys/JM11/2016-11-25_18-34-06/processed/2016-11-25_18-34-06.lfp',...
% $$$                  145:164,Par.nChannels,4)';    
% $$$ Session.spath = '/storage/javier/data/processed/xyz/JM11-20161126/';
% $$$ Session.maze.name = 'cof';
% $$$ Session.name = 'JM11-20161126';
% $$$ [xyzData, markers, mocapSampleRate] = concatViconFiles(Session);            
% $$$ xyzSampleRate = 199.997752;
% $$$ (diff(ThreshCross(lfp(:,13),1e4,1),1,2)-4)./1250-cellfun(@length,xyzData(1:5))'./xyzSampleRate
