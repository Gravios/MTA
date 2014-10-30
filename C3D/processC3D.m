function [ output_args ] = processC3D( varargin )
%PROCESSC3D Summary of this function goes here
%   Detailed explanation goes here
%subjectName = 'rat-kenny';
%sessionName = '2010-10-27';
%'F:\Vicon_Data\Vicon_Right\Top Level 5\';
%'F:\Vicon_Data\vicon\Top Level 3\'

[sessionName,destPath, mazeName, subjectName, username, password, server, sourcePath] = ...
DefaultArgs(varargin,{[],'F:\data\xyz\',[],'Ed10','justi_000',[],'172.25.9.20','F:\Vicon_Data\testdatabase\Project 2\'});



if isempty(destPath),
    system('net use z: /DELETE');
    system(['net use z: \\' server '\' username ' ' password ' /USER:' username ]); 
    destPath = 'z:\data\xyz\';
end



if isempty(subjectName),
    subjectList= dir(sourcePath);
    subjectList(1:2) = [];
    subjectList(~[subjectList(:).isdir]) = [];
    for k = 1:numel(subjectList),
        if isempty(sessionName),
            sessionList = dir(fullfile(sourcePath,subjectList(k).name));
            sessionList(1:2) = [];
            sessionList(~[sessionList(:).isdir]) = [];
            for i = 1:numel(sessionList),
                if ~exist(fullfile(destPath,sessionList(i).name),'dir'),
                    mkdir(fullfile(destPath,sessionList(i).name));
                end                
                if isempty(mazeName)
                    mazeList = dir(fullfile(sourcePath,subjectList(k).name,sessionList(i).name));
                    mazeList(1:2) = [];
                    mazeList(~[mazeList(:).isdir]) = [];
                    for j = 1:numel(mazeList),
                        if ~exist(fullfile(destPath,sessionName,mazeList(j).name),'dir'),
                            mkdir(fullfile(destPath,sessionName,mazeList(j).name));
                        end
                        %pc3d is a helper function found at the end of this file
                        pc3d(sourcePath,destPath,subjectList(k).name,sessionList(i).name,mazeList(j).name)
                    end
                else
                    if ~exist(fullfile(destPath,sessionList(i).name,mazeName),'dir'),
                        mkdir(fullfile(destPath,sessionList(i).name,mazeName));
                    end
                    pc3d(sourcePath,destPath,subjectList(k).name,sessionList(i).name,mazeName)
                end
                
            end
        end
    end
else
    if isempty(sessionName),
        sessionList = dir(fullfile(sourcePath,subjectName));
        sessionList(1:2) = [];
        sessionList(~[sessionList(:).isdir]) = [];
        for i = 1:numel(sessionList),
            if ~exist(fullfile(destPath,sessionList(i).name),'dir'),
                mkdir(fullfile(destPath,sessionList(i).name));
            end
            if isempty(mazeName)
                mazeList = dir(fullfile(sourcePath,subjectName,sessionList(i).name));
                mazeList(1:2) = [];
                mazeList(~[mazeList(:).isdir]) = [];
                for j = 1:numel(mazeList),
                    if ~exist(fullfile(destPath,sessionList(i).name,mazeList(j).name),'dir'),
                        mkdir(fullfile(destPath,sessionList(i).name,mazeList(j).name));
                    end
                    %pc3d is a helper function found at the end of this file
                    pc3d(sourcePath,destPath,subjectName,sessionList(i).name,mazeList(j).name)
                end
            else
                if ~exist(fullfile(destPath,sessionList(i).name,mazeName),'dir'),
                    mkdir(fullfile(destPath,sessionList(i).name,mazeName));
                end
                pc3d(sourcePath,destPath,subjectName,sessionList(i).name,mazeName)
            end
        end
    else
        if ~exist(fullfile(destPath,sessionName),'dir'),
            mkdir(fullfile(destPath,sessionName));
        end
        if isempty(mazeName)
            mazeList = dir(fullfile(sourcePath,subjectName,sessionName));
            mazeList(1:2) = [];
            mazeList(~[mazeList(:).isdir]) = [];
            for j = 1:numel(mazeList),
                if ~exist(fullfile(destPath,sessionName,mazeList(j).name),'dir'),
                    mkdir(fullfile(destPath,sessionName,mazeList(j).name));
                end
                %pc3d is a helper function found at the end of this file
                pc3d(sourcePath,destPath,subjectName,sessionName,mazeList(j).name)
            end
        else
            if ~exist(fullfile(destPath,sessionName,mazeName),'dir'),
                mkdir(fullfile(destPath,sessionName,mazeName));
            end
            pc3d(sourcePath,destPath,subjectName,sessionName,mazeName)
        end
    end
end





% Close network drive if open
if isempty(destPath),
    try
        system('net use z: /DELETE');
    end
end
end

function pc3d(sourcePath,destPath,subjectName,sessionName,mazeName)
    fileList = dir(fullfile(sourcePath,subjectName,sessionName,mazeName,'*.c3d'));
    % Convert all available C3D files and save each as a .mat file.
    if isempty(fileList),return,end
    itf = c3dserver;    
    for i = 1:length(fileList),
        openc3d(itf,0,fullfile(sourcePath,subjectName,sessionName,mazeName,fileList(i).name));
        [xyzpos,markers] = getAll3dTargets(itf);
        sampleRate = itf.GetVideoFrameRate;
        closec3d(itf);
        fileName = fullfile(destPath,sessionName,mazeName,[sessionName '-' fileList(i).name '.mat' ]);
        save(fileName,'xyzpos','markers','sampleRate');
    end
    copyfile(fullfile(sourcePath,subjectName,sessionName,mazeName,[sessionName '-' mazeName '.vsk' ]),fullfile(destPath,sessionName,mazeName,[sessionName '-' mazeName '.vsk' ]));
    itf.release;
end