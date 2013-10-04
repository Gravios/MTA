function [ output_args ] = exportC3D( sessionName, subjectName, varargin )
%PROCESSC3D Summary of this function goes here
%   Detailed explanation goes here
%subjectName = 'rat-kenny';
%sessionName = '2010-10-27';
%'C:\Vicon_Data\Vicon_Down\Top Level 2\'
%[ destPath, server, sourcePath ] = DefaultArgs(varargin,{'x:\data\xyz\','172.25.9.20','C:\Vicon_Data\Vicon_Down\Top Level 2\'})
destPath='x:\data\xyz\';
server='172.25.9.20';
sourcePath='C:\Vicon_Data\Vicon_Down\Top Level 2\';

[username, password] = logindlg(['Login:' server],'User');

system('net use x: /DELETE');
system(['net use x: \\' server '\' username ' ' password ' /USER:' username ]); 



if ~exist(sourcePath,'dir'),
    error('Path does not exist')
    return
elseif ~exist([sourcePath subjectName],'dir'),
    error('Subject does not exist')
    return
elseif ~exist([sourcePath  subjectName '\' sessionName],'dir'),
    error('Session does not exist')
    return
elseif ~exist(destPath,'dir'),
    error('Destination path does not exist')
    return
elseif ~exist([destPath  sessionName],'dir');
    mkdir([destPath  sessionName]);
end    


% Open a c3dserver to conver from C3D to MAT
itf = c3dserver;
fileList = dir([sourcePath subjectName '\' sessionName '\' '*.c3d']);
% Convert all available C3D files and save each as a .mat file.
for i = 1:length(fileList),
    openc3d(itf,0,[sourcePath subjectName '\' sessionName '\' fileList(i).name]);
    [xyzpos,markers] = getAll3dTargets(itf);
    closec3d(itf);
    fileName = [destPath sessionName '\' sessionName '-' fileList(i).name '.mat'];
    save(fileName,'xyzpos','markers');
end

% Close Network Drive
system('net use x: /DELETE');

clear

end
