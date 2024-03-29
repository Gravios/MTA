function Session = syncViconNlx(Session,TTLValue,viconSampleRate)
% Session = create(Session,TTLValue)
% Populate a Session with xyz data synchronized to
% electrophysiological data based on an event file
%
% note - the event file must be located in the nlx folder
%        of the MTA directory tree
%
ERR.type = 'MTASession:create';            
ERR.msg  = ['SyncPeriods is empty, check event file that ' ...
            'the TTLValues corresponding to the recording ' ...
            'trials is equal to: %s'];
try
    Par = LoadPar(fullfile(Session.spath, [Session.name '.xml']));            
catch
    warning(['xml parameter file in the nlx folder may ' ...
             'not exist'])
end

%% Because users make me cry 
% $$$ delete(fullfile(Session.spath,['*',Session.maze.name,'.all*']));
% $$$ delete(fullfile(Session.spath,['*.ses.*']));
% $$$ delete(fullfile(Session.spath,['*.trl.*']));


if exist('Par','var'),
    %% Load single channel of lfp to check the exact number of samples
    lfp = LoadBinary(fullfile(Session.spath, [Session.name '.lfp']),1,Par.nChannels,4)';
    recordSync = [0,numel(lfp)./Par.lfpSampleRate];
    lfpSyncPeriods = MTADepoch([],[],recordSync+[1/Par.lfpSampleRate,1/Par.lfpSampleRate],1,recordSync,0);
    clear('lfp');
    Session.lfp = MTADlfp(Session.spath,[Session.name '.lfp'],[],Par.lfpSampleRate,lfpSyncPeriods,0);    
    Session.lfp.filename = [Session.name '.lfp'];
    Session.sampleRate = Par.SampleRate;
end

if ~exist(Session.spath,'dir')
    if ~exist(Session.spath,'dir') && ...
        exist(fullfile(Session.path.project,'nlx',Session.name),'dir') && ...
        exist(fullfile(Session.path.project,'xyz',Session.name),'dir')
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
if isempty(viconSampleRate),
    [xyzData, markers, viconSampleRate] = concatViconFiles(Session);            
    if isempty(viconSampleRate),
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

events = LoadEvents(fullfile(Session.spath, [Session.name '.all.evt']));
stopTTL = '0x0000';

pfirstVStart = find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))),1,'first')-1;
events.time(1:pfirstVStart)=[];
events.Clu(1:pfirstVStart)=[];
events.description(1:pfirstVStart)=[];

vstarts = events.time(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))));
vstops  = events.time(find(ismember(events.Clu,find(~cellfun(@isempty,regexp(events.Labels,stopTTL)))),numel(vstarts),'first'));
tsRanges = [vstarts,vstops];
            
%% MUA HAHAHAHAHAAAaaaaaa 
% but seriously it just finds the events with the TTLValue and then
% the next event which corresponds to the stopTTL
% $$$     tsRanges = [events.time(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue)))),...
% $$$        events.time(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))))+...
% $$$ cell2mat(cellfun(@find,...                   
% $$$ cellfun(@eq,...                   
% $$$ cellfun(@subsref,repmat({events},numel(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))))),1),cellfun(@substruct,...
% $$$ repmat({'.'},numel(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))))),1),...
% $$$ repmat({'Clu'},numel(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))))),1),...
% $$$ repmat({'()'},numel(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))))),1),...
% $$$ mat2cell(cellfun(@colon,mat2cell(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue)))),ones(numel(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))))),1),1),repmat({numel(events.Clu)},numel(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))))),1),'Uniformoutput',false),ones(numel(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))))),1),1),'Uniformoutput',false),'Uniformoutput',false),repmat({find(~cellfun(@isempty,regexp(events.Labels,stopTTL)))},numel(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))))),1),'Uniformoutput',false),repmat({1},numel(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))))),1),repmat({'first'},numel(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))))),1),'uniformoutput',false))-1)];


%% Assign xyz data to nlx event epochs
syncPeriods=[];
xyzDataInd = [];
for i=1:numel(xyzData),
    if ~isempty(xyzData{i})
        for j=1:size(tsRanges,1),
            %diff(tsRanges(j,:))-size(xyzData{i},1)/viconSampleRate
            if abs(diff(tsRanges(j,:))*10-size(xyzData{i},1)/viconSampleRate*10)<0.2,
                % ???Shift time back so index corresponds to
                % syncPeriods(1) :-1/viconSampleRate???
                syncPeriods(end+1,:) = tsRanges(j,:)-1/viconSampleRate;

                tsRanges(tsRanges(j,1)>=tsRanges(:,1),:) = [];
                xyzDataInd(end+1) = i;
                break;
            end
        end
    end
end


%% Error if no NLX events match xyzData
assert(~isempty(syncPeriods),ERR.type,ERR.msg,TTLValue);

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
xyz = zeros([ceil(diff(syncPeriods([1,end]))*viconSampleRate),size(xyzData{1},2),size(xyzData{1},3)]);
syncXyzStart = round((syncPeriods(1:length(syncPeriods))-syncPeriods(1))*viconSampleRate+1);
%syncShift = [0;floor(diff([syncPeriods(1:end-1,2) syncPeriods(2:end,1)],1,2)*viconSampleRate)];
for s=1:nViconTrials,
% $$$     if syncShift(s),
% $$$         xyz = cat(1,xyz,zeros(syncShift(s),size(xyzData{s},2),size(xyzData{s},3)));
% $$$     end
xyzseg = xyzData{s};
xyzseg(xyzseg==0)=eps;
xyz(syncXyzStart(s):syncXyzStart(s)+xyzLengths(s)-1,:,:) = xyzseg;
% $$$     [size(xyz,1)/viconSampleRate+syncPeriods.data(1)+1/viconSampleRate,syncPeriods.data(s)]
% $$$     diff([size(xyz,1)/viconSampleRate+syncPeriods.data(1)+1/viconSampleRate,syncPeriods.data(s)])
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

Session.xyz = MTADxyz(Session.spath,Session.filebase,xyz,viconSampleRate,...
                      syncPeriods,Session.sync.data(1),Session.model);                  
Session.xyz.save;
Session.xyz.clear;

Session.ang = MTADang(Session.spath,...
                      Session.filebase,...
                      [],...
                      viconSampleRate,...
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


