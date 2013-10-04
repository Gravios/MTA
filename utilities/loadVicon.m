function Session = loadVicon(Session)

[xyzData, markers] = concatViconFiles(Session);            

vsk_path = fullfile(Session.spath, Session.Maze.name, [Session.name '-' Session.Maze.name '.vsk']);
if exist(vsk_path,'file'),
    Session.Model = MTAModel(vsk_path,'-vsk');
else
    warning(['VSK file associated with this session was ' ...
             'not found. \nCreating a general marker model.']);
    Session.Model = MTAModel(markers,'-mar');
end

xyzDataInd = find(~cellfun(@isempty,xyzData));
xyz = cell2mat(xyzData(xyzDataInd)');
Session.xyzSegLength = cellfun(@length,xyzData(xyzDataInd));                        
xsl = [0,cumsum(Session.xyzSegLength)];
Session.xyzPeriods = [xsl(1:end-1)+1;xsl(2:end)]';

Dxyz = MTADxyz(fullfile(Session.spath, Session.Maze.name, [Session.filebase '.xyz.mat']),xyz,'TimeSeries','xyz',viconSampleRate);

Session.xyz = Dxyz;