function link_session_to_project(sessionName,xyzPath,nlxPath)
%function linkSession(sessionName,xyzPath,nlxPath)
%
%  Creates Session folder in the current MTA project
%  and creates symlinks pointed to folders containing
%  session data.
% 
% TODO make linkSession windows compatible 
% UPDATE requires link shell extension

Session = MTASession([]);
datDir = fullfile(Session.path.data,sessionName);
try,mkdir(datDir);end

if ~isempty(xyzPath),
    % get list of directories along the path xyzPath
    xyz_maze_dirs = dir(fullfile(xyzPath,sessionName));xyz_maze_dirs(1:2)=[];
    if ~strcmp(xyzPath,fullfile(Session.path.data,'xyz')),
        xyzDir = fullfile(Session.path.data,'xyz',sessionName);
        create_directory(xyzDir);
        for maze = 1:numel(xyz_maze_dirs),
            create_directory(fullfile(xyzDir,xyz_maze_dirs(maze).name));
            system(['find ',fullfile(xyzDir,xyzMazeDirs(maze).name),' -type l -delete']);            
            system(['ln -sf ' fullfile(xyzPath,sessionName,xyz_maze_dirs(maze).name) '/* ' fullfile(xyzDir,xyz_maze_dirs(maze).name)]);
        end
    else
        xyzDir = fullfile(xyzPath,sessionName);
    end
end


% Link Nlx data to nlx folder
if exist('nlxPath','var')
    if ~isempty(nlxPath),
        if ~strcmp(nlxPath,fullfile(Session.path.data,'nlx')),
            nlxDir = fullfile(Session.path.data,'nlx',sessionName);
            create_directory(nlxDir);
            system(['find ',nlxDir,' -type l -delete']);
            if numel(dir(fullfile(nlxPath,sessionName)))>2,
                system(['ln -sf ' fullfile(nlxPath,sessionName) '/* ' nlxDir]);
            end
        else
            % if nlxPath is within the MTA data collection
            nlxDir = fullfile(nlxPath,sessionName);
        end
        cd(datDir);
        system(['find ',datDir,' -type l -delete']);                
        system(['ln -s ../nlx/' sessionName '/* ' datDir ]);
    end
end


% For each maze in the xyz data link directory 
if exist('xyzPath','var')    
    if ~isempty(xyzPath),
        for maze = 1:numel(xyz_maze_dirs),
            cd(datDir);
            try,mkdir(xyz_maze_dirs(maze).name);end
            
            try
                system(['ln -sf ../xyz/' sessionName '/' xyz_maze_dirs(maze).name '/' ...
                        sessionName '-' xyz_maze_dirs(maze).name '.vsk ' ...
                        datDir '/' sessionName '-' xyz_maze_dirs(maze).name '.vsk ']);
            end
            
            cd(fullfile(datDir,xyz_maze_dirs(maze).name))
            try,
                system(['find ',datDir,'/',xyzMazeDirs(maze).name,'/ -type l -delete']);
                system(['ln -s ../../xyz/' sessionName '/' xyz_maze_dirs(maze).name '/* ' datDir '/' xyz_maze_dirs(maze).name '/']);
            end
            
        end
    end
end


