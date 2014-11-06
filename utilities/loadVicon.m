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
Dxyz = MTADxyz(Session.spath,Session.filebase,...
               xyz,viconSampleRate);
Dxyz.save;

Session.xyz = Dxyz;