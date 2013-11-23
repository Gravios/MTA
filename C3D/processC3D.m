function [ output_args ] = processC3D( sessionName, destPath, mazeName, subjectName, username, password, server, sourcePath )
%PROCESSC3D Summary of this function goes here
%   Detailed explanation goes here
%subjectName = 'rat-kenny';
%sessionName = '2010-10-27';
if (nargin < 8)||isempty(sourcePath),
    sourcePath = 'C:\Vicon_Data\Vicon_Down\Top Level 2\';
end
if (nargin < 7)||isempty(server),
    server = '172.25.9.20';
end
if (nargin < 6)||isempty(password),
    password = '';
end
if (nargin < 5)||isempty(username),
    username = 'ricardo';
end
if (nargin < 4)||isempty(subjectName),
    subjectName = 'rs01';
end
if (nargin < 3)||isempty(mazeName),
    mazeName = 'cof';
end
if (nargin < 2)||isempty(destPath),
    system('net use z: /DELETE');
    system(['net use z: \\' server '\' username ' ' password ' /USER:' username ]); 
    destPath = 'z:\data\xyz\';
    if ~exist([destPath '\'  sessionName ],'dir'),
        mkdir([destPath '\'  sessionName '\']);
    end
    if ~exist([destPath '\'  sessionName '\' mazeName],'dir'),
        mkdir([destPath '\'  sessionName '\' mazeName]);
    end
end

if ~exist(sourcePath,'dir'),
    error('Path does not exist')
    return
elseif ~exist([sourcePath subjectName],'dir'),
    error('Subject does not exist')
    return
elseif ~exist([sourcePath  subjectName '\' sessionName '\' mazeName ],'dir'),
    error('Session does not exist')
    return
elseif ~exist(destPath,'dir'),
    error('Destination path does not exist')
    return
elseif ~exist(destPath,'dir'),
    mkdir(destPath);
end
    

% Open a c3dserver to conver from C3D to MAT
itf = c3dserver;
fileList = dir([sourcePath subjectName '\' sessionName '\' mazeName '\' '*.c3d']);
% Convert all available C3D files and save each as a .mat file.
for i = 1:length(fileList),
    openc3d(itf,0,[sourcePath subjectName '\' sessionName '\' mazeName '\' fileList(i).name]);
    [xyzpos,markers] = getAll3dTargets(itf);
    closec3d(itf);
    fileName = [destPath sessionName '\' mazeName '\' sessionName '-' fileList(i).name '.mat' ];
    save(fileName,'xyzpos','markers');
end

% Close Network Drive
system('net use z: /DELETE');

end