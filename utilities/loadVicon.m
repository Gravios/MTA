function Session = loadVicon(Session,viconSampleRate)

% Load concatinate the xyz data from all c3d.mat files associated with the
% Session
if isempty(viconSampleRate),
    [xyzData, markers, viconSampleRate] = concatViconFiles(Session);            
else
    [xyzData, markers] = concatViconFiles(Session);            
end

% Load VSK to get marker names,colors and connections for creating the
% model
vsk_path = fullfile(Session.spath, Session.maze.name, [Session.name '-' Session.maze.name '.vsk']);
if exist(vsk_path,'file'),
    Session.model = MTAModel(vsk_path,'-vsk');
else
    warning(['VSK file associated with this session was ' ...
             'not found. \nCreating a general marker model.']);
    Session.model = MTAModel(markers,'-mar');
end

% Create the xyzPeriods and concatinate all trials into one xyz array
xyz = cell2mat(xyzData(~cellfun(@isempty,xyzData))');

% Create XYZ data object


syncPeriods = cellfun(@length,xyzData);
syncPeriods(syncPeriods==0)=[];
syncPeriods = [ cumsum([1,syncPeriods(1:end-1)]);cumsum(syncPeriods)]'...
               ./viconSampleRate;
syncPeriods(1) = 0;

Session.sync = MTADepoch(Session.spath,[Session.filebase '.sync.mat'],syncPeriods([1,end]),1,0,0,[],[],[],'sync');


syncPeriods = MTADepoch([],[],syncPeriods,1,Session.sync.copy,0);



Dxyz = MTADxyz(Session.spath,Session.filebase,xyz,viconSampleRate,...
               syncPeriods,0,Session.model);
Dxyz.save;

Session.xyz = Dxyz;

Session.ang = MTADang(Session.spath,Session.filebase,[],viconSampleRate,...
                      Session.xyz.sync,Session.xyz.origin,Session.model);

%% MTAStateCollection object holds all behavioral sets of periods
Session.stc = MTAStateCollection(Session.spath,Session.filebase,'default',[],[],1);
Session.stc.updateSync(Session.sync);
Session.stc.updateOrigin(0);
