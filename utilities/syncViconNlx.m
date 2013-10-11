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
    Session.lfp = MTADlfp(Session.spath,[Session.name '.lfp'],[],Par.lfpSampleRate);
    Session.sampleRate = Par.SampleRate;
end

if ~exist(Session.spath,'dir')
    mkdir(Session.spath);
end


%% Organize the Sessions trials into a cell array
%% Return the names of the markers
[xyzData, markers] = concatViconFiles(Session);            

%% Load VSK if possible 
vsk_path = fullfile(Session.spath, [Session.name '-' Session.Maze.name '.vsk']);
if exist(vsk_path,'file'),
    Session.Model = MTAModel(vsk_path,'-vsk');
else
    warning(['VSK file associated with this session was ' ...
             'not found. \nCreating a general marker model.']);
    Session.Model = MTAModel(markers,'-mar');
end

%assert(exist([Session.spath Session.name '.all.evt'],'file'))
% Load events
events = textread(fullfile(Session.spath, [Session.name '.all.evt']),'%s','delimiter','\n','emptyvalue',NaN);
% Parse event type
parsedEvents = {};
aborentEvents = {};
for i=1:length(events),
    rStart = regexpi(events{i},'starting recording');
    rStop = regexpi(events{i},'stopping recording');
    eStart = regexpi( events{i},TTLValue);
    eStop = regexpi(events{i},'0x0000');
    ts = events{i}(1:regexpi(events{i},' ','once')-1);
    
    if rStart,
        parsedEvents{end+1} = struct('ts',ts,'event_name','starting_recording','ttl_value','');
    elseif eStart,
        parsedEvents{end+1} = struct('ts',ts,'event_name','vicon_start','ttl_value',TTLValue);
    elseif sum(eStop) && strcmp(parsedEvents{end}.ttl_value, TTLValue),
        parsedEvents{end+1} = struct('ts',ts,'event_name', 'vicon_stop','ttl_value','0x0000');
    elseif rStop,
        parsedEvents{end+1} = struct('ts',ts,'event_name','stopping_recording','ttl_value','');
    else
        aborentEvents{end+1} = events{i};
    end
end

% Create array of time spans between events
tsRanges = [];

for i=1:length(parsedEvents),
    if strcmp(parsedEvents{i}.event_name,'vicon_stop')&&strcmp(parsedEvents{i-1}.event_name,'vicon_start'),
        tsRanges(end+1,:) = [str2double(parsedEvents{i-1}.ts) str2double(parsedEvents{i}.ts)];
    end
end

syncPeriods=[];
xyzDataInd = [];
for i=1:length(xyzData),
    if ~isempty(xyzData{i})
        for j=1:size(tsRanges,1),
            if round(diff(tsRanges(j,:))/100)==round(length(xyzData{i})/viconSampleRate*10),
                syncPeriods(end+1,:) = tsRanges(j,:);
                tsRanges(j,:) = [];
                xyzDataInd(end+1) = i;
                break;
            end
        end
    end
end

assert(~isempty(syncPeriods),ERR.type,ERR.msg,TTLValue);

xyzData = xyzData(xyzDataInd);

Session.sync = MTASync(Session.spath,Session.name,syncPeriods./1000);
Session.stc = MTAStateCollection(Session.spath,Session.filebase,'default');


nSessions = length(xyzData);
xyz = [];
syncShift = [0;round(diff([Session.sync(1:end-1,2) Session.sync(2:end,1)],1,2)*viconSampleRate)];
for s=1:nSessions,
    if syncShift(s),
        xyz = cat(1,xyz,zeros(syncShift(s),size(xyzData{s},2),size(xyzData{s},3)));
    end
    xyz = cat(1,xyz,xyzData{s});
end
xyz = double(xyz);


Session.spk = MTASpk;

Session.xyz = MTADxyz(Session.spath,Session.filebase,xyz,viconSampleRate);
Session.xyz.save;

Session.ang = MTADang(Session.spath,Session.filebase,[],viconSampleRate);
Session.ang.save;

Session.save();

end

