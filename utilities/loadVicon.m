function Session = loadVicon(Session,viconSampleRate)

% Load concatinate the xyz data from all c3d.mat files associated with the
% Session
[xyzData, markers] = concatViconFiles(Session);            

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
xyzDataInd = find(~cellfun(@isempty,xyzData));
xyz = cell2mat(xyzData(xyzDataInd)');

% Create XYZ data object
Dxyz = MTADxyz(Session.spath,Session.filebase,...
               xyz,viconSampleRate);
Dxyz.save;

Session.xyz = Dxyz;