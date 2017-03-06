function link_session_to_project(session_name,xyz_path,nlx_path)
%function linkSession(session_name,xyz_path,nlx_path)
%
%  Creates Session folder in the current MTA project
%  and creates symlinks pointed to folders containing
%  session data.
% 
% TODO make linkSession windows compatible 
% UPDATE requires link shell extension

Session = MTASession([]);
datDir = fullfile(Session.path.data,session_name);
try,mkdir(datDir);end

if ~isempty(xyz_path),
    % get list of directories along the path xyz_path
    xyz_maze_dirs = dir(fullfile(xyz_path,session_name));xyz_maze_dirs(1:2)=[];
    if ~strcmp(xyz_path,fullfile(Session.path.data,'xyz')),
        xyzDir = fullfile(Session.path.data,'xyz',session_name);
        try,mkdir(xyzDir);end
        for maze = 1:numel(xyz_maze_dirs),
            try,mkdir(fullfile(xyzDir,xyz_maze_dirs(maze).name));end
            system(['ln -sf ' fullfile(xyz_path,session_name,xyz_maze_dirs(maze).name) '/* ' fullfile(xyzDir,xyz_maze_dirs(maze).name)]);
        end
    else
        xyzDir = fullfile(xyz_path,session_name);
    end
end


% Link Nlx data to nlx folder
if exist('nlx_path','var')
    if ~isempty(nlx_path),
        if ~strcmp(nlx_path,fullfile(Session.path.data,'nlx')),
            nlxDir = fullfile(Session.path.data,'nlx',session_name);
            mkdir(nlxDir);
            if numel(dir(fullfile(nlx_path,session_name)))>2,
                system(['ln -sf ' fullfile(nlx_path,session_name) '/* ' nlxDir]);
            end
        else
            % if nlx_path is within the MTA data collection
            nlxDir = fullfile(nlx_path,session_name);
        end
        cd(datDir);
        system(['ln -sf ../nlx/' session_name '/* ' datDir ])
        
    end
end


% For each maze in the xyz data link directory 
if exist('xyz_path','var')    
    if ~isempty(xyz_path),
        for maze = 1:numel(xyz_maze_dirs),
            cd(datDir);
            try,mkdir(xyz_maze_dirs(maze).name);end
            
            try
                system(['ln -sf ../xyz/' session_name '/' xyz_maze_dirs(maze).name '/' ...
                        session_name '-' xyz_maze_dirs(maze).name '.vsk ' ...
                        datDir '/' session_name '-' xyz_maze_dirs(maze).name '.vsk ']);
            end
            cd(fullfile(datDir,xyz_maze_dirs(maze).name))
            try,
                system(['ln -sf ../../xyz/' session_name '/' ...
                        xyz_maze_dirs(maze).name '/* ' datDir '/' ...
                        xyz_maze_dirs(maze).name '/']);
            end
            
        end
    end
end


