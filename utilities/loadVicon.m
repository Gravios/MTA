function Session = loadVicon(Session)

% Load concatinate the xyz data from all c3d.mat files associated with the
% Session
[xyzData, markers] = concatViconFiles(Session);            

% Load VSK to get marker names,colors and connections for creating the
% model
vsk_path = fullfile(Session.spath, Session.Maze.name, [Session.name '-' Session.Maze.name '.vsk']);
if exist(vsk_path,'file'),
    Session.Model = MTAModel(vsk_path,'-vsk');
else
    warning(['VSK file associated with this session was ' ...
             'not found. \nCreating a general marker model.']);
    Session.Model = MTAModel(markers,'-mar');
end

% Create the xyzPeriods and concatinate all trials into one xyz array
xyzDataInd = find(~cellfun(@isempty,xyzData));
xyz = cell2mat(xyzData(xyzDataInd)');
Session.xyzSegLength = cellfun(@length,xyzData(xyzDataInd));                        
xsl = [0,cumsum(Session.xyzSegLength)];
Session.xyzPeriods = [xsl(1:end-1)+1;xsl(2:end)]';

% Create XYZ data object
Dxyz = MTADxyz(Session.spath,Session.filebase,...
               xyz,viconSampleRate);
Dxyz.save;

Session.xyz = Dxyz;