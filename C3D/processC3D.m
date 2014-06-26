function [ output_args ] = processC3D( sessionName , varargin )
%PROCESSC3D Summary of this function goes here
%   Detailed explanation goes here
%subjectName = 'rat-kenny';
%sessionName = '2010-10-27';

[destPath, mazeName, subjectName, username, password, server, sourcePath] = DefaultArgs(varargin,{[],'cof','er09','justi_000',[],'172.25.9.20','F:\Vicon_Data\vicon\Top Level 3\'});

if isempty(destPath),
    system('net use z: /DELETE');
    system(['net use z: \\' server '\' username ' ' password ' /USER:' username ]); 
    destPath = 'z:\data\xyz\';
end


if ~exist(fullfile(destPath,sessionName),'dir'),
    mkdir(fullfile(destPath,sessionName));
end
if ~exist(fullfile(destPath,sessionName,mazeName),'dir'),
    mkdir(fullfile(destPath,sessionName,mazeName));
end

% if ~exist(sourcePath,'dir'),
%     error('Path to data does not exist')
% elseif ~exist(fullfile(sourcePath,subjectName),'dir'),
%     error('Path to subject does not exist')
% elseif ~exist(fullfile(sourcePath,subjectName,sessionName,mazeName),'dir'),
%     error('Path to session does not exist')
% end
% 
% if ~exist(destPath,'dir'),
%     error('Path to destination does not exist')
% end
% if ~exist(destPath,'dir'),
%     mkdir(destPath);
% end
%     


% List All c3d files in source directory for conversion
fileList = dir(fullfile(sourcePath,subjectName,sessionName,mazeName,'*.c3d'));


% Open a c3dserver to conver from C3D to MAT
itf = c3dserver;


% Convert all available C3D files and save each as a .mat file.
for i = 1:length(fileList),
    openc3d(itf,0,fullfile(sourcePath,subjectName,sessionName,mazeName,fileList(i).name));
    [xyzpos,markers] = getAll3dTargets(itf);
    closec3d(itf);
    fileName = fullfile(destPath,sessionName,mazeName,[sessionName '-' fileList(i).name '.mat' ]);
    save(fileName,'xyzpos','markers');
end

% Close network drive if open
try
    system('net use z: /DELETE');
end

end