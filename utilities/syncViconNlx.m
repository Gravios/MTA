function Session = syncViconNlx(Session,viconSampleRate,TTLValue)
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

if exist('Par','var'),
    
    lfp = LoadBinary(fullfile(Session.spath, [Session.name '.lfp']),1,Par.nChannels,4)';
    recordSync = [0,numel(lfp)./Par.lfpSampleRate];
    lfpSyncPeriods = MTADepoch([],[],recordSync*Par.lfpSampleRate+[1,0],Par.lfpSampleRate,recordSync,0);
    clear('lfp');
    Session.lfp = MTADlfp(Session.spath,[Session.name '.lfp'],[],Par.lfpSampleRate,lfpSyncPeriods,0);    
    Session.sampleRate = Par.SampleRate;
end

if ~exist(Session.spath,'dir')
    mkdir(Session.spath);
end


%% Organize the Sessions trials into a cell array
%% Return the names of the markers
[xyzData, markers] = concatViconFiles(Session);            

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
tsRanges = [events.time(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue)))),...
       events.time(find(events.Clu==find(~cellfun(@isempty,regexp(events.Labels,TTLValue))))+1)];


%% Assign xyz data to nlx event epochs
syncPeriods=[];
xyzDataInd = [];
for i=1:numel(xyzData),
    if ~isempty(xyzData{i})
        for j=1:size(tsRanges,1),
            if abs(diff(tsRanges(j,:))*10-size(xyzData{i},1)/viconSampleRate*10)<0.2,
                syncPeriods(end+1,:) = tsRanges(j,:);

                tsRanges(tsRanges(j,1)>=tsRanges(:,1),:) = [];
                xyzDataInd(end+1) = i;
                break;
            end
        end
    end
end


assert(~isempty(syncPeriods),ERR.type,ERR.msg,TTLValue);

xyzData = xyzData(xyzDataInd);


Session.sync = MTADepoch(Session.spath,[Session.filebase '.sync.mat'],syncPeriods([1,end]),1,recordSync,0,[],[],[],'sync');
Session.sync.save(1);
Session.lfp.sync.sync = Session.sync.copy;
Session.lfp.origin = round(Session.lfp.sync.sync.data(1)*Par.lfpSampleRate);
Session.stc = MTAStateCollection(Session.spath,Session.filebase,'default',[],[],1);
Session.stc.updateSync(Session.sync);
Session.stc.updateOrigin(0);

nSessions = length(xyzData);
xyz = [];
syncShift = [0;round(diff([syncPeriods(1:end-1,2) syncPeriods(2:end,1)],1,2)*viconSampleRate)];
for s=1:nSessions,
    if syncShift(s),
        xyz = cat(1,xyz,zeros(syncShift(s),size(xyzData{s},2),size(xyzData{s},3)));
    end
    xyz = cat(1,xyz,xyzData{s});
end
xyz = double(xyz);


Session.spk = MTASpk;
Session.spk.create(Session);

syncPeriods = MTADepoch([],[],syncPeriods,1,Session.sync.copy,0);
syncPeriods.resample(viconSampleRate);

Session.xyz = MTADxyz(Session.spath,Session.filebase,xyz,viconSampleRate,...
                      syncPeriods,round(Session.sync.data(1)*viconSampleRate),Session.model);                  
Session.xyz.save;

Session.ang = MTADang(Session.spath,Session.filebase,[],viconSampleRate,...
                      Session.xyz.sync,Session.xyz.origin,Session.model);
Session.ang.save;

Session.ufr = MTADufr(Session.spath,Session.filebase);

Session.save();
Session.spk.clear;

end

