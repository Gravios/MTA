function linkSession(session_name,xyz_path,nlx_path)
%function linkSession(session_name,xyz_path,nlx_path)
    
    Session = MTASession([]);
    datDir = fullfile(Session.path.data,session_name);
    try,mkdir(datDir);end

    xyz_maze_dirs = dir(fullfile(xyz_path,session_name));xyz_maze_dirs(1:2)=[];

    if ~strcmp(xyz_path,fullfile(Session.path.data,'xyz')),
        xyzDir = fullfile(Session.path.data,'xyz',session_name);
        try,mkdir(xyzDir);end
        
        for maze = 1:numel(xyz_maze_dirs),
            try,mkdir(fullfile(xyzDir,xyz_maze_dirs(maze).name));end
            system(['ln -s ' fullfile(xyz_path,session_name,xyz_maze_dirs(maze).name) '/* ' fullfile(xyzDir,xyz_maze_dirs(maze).name)]);
        end
    else
        
        xyzDir = fullfile(xyz_path,session_name);
    end

    
    %% Link Nlx data to nlx folder
    if ~strcmp(nlx_path,fullfile(Session.path.data,'nlx')),
        nlxDir = fullfile(Session.path.data,'nlx',session_name);
        mkdir(nlxDir);
        system(['ln -s ' fullfile(nlx_path,session_name) '/* ' nlxDir]);
    else
        % if nlx_path is within the MTA data collection
        nlxDir = fullfile(nlx_path,session_name);
    end
    
    cd(datDir);
    system(['ln -s ../nlx/' session_name '/* ' datDir ])

    for maze = 1:numel(xyz_maze_dirs),
        try,mkdir(xyz_maze_dirs(maze).name);end
        system(['ln -s ../xyz/' session_name '/' xyz_maze_dirs(maze).name '/' ...
                session_name '-' xyz_maze_dirs(maze).name '.vsk ' ...
                datDir '/']);
        cd(xyz_maze_dirs(maze).name)
        system(['ln -s ../../xyz/' session_name '/' ...
                xyz_maze_dirs(maze).name '/* ' datDir '/' ...
               xyz_maze_dirs(maze).name '/']);
    end

    
    
% $$$     if exist([target_path.xyz session_name],'file')&&exist([target_path.nlx session_name],'file'),
% $$$ 
% $$$        mkdir([ paths.xyz session_name ]);
% $$$        mkdir([ paths.nlx session_name ]);
% $$$        mkdir([ paths.analysis session_name ]);
% $$$ 
% $$$        subsessions = dir([target_path.xyz session_name]);
% $$$        subsessions = subsessions(3:end);
% $$$        for i = 1:length(subsessions),
% $$$            mkdir([ paths.xyz session_name '/' subsessions(i).name]);
% $$$            session_parts = dir([target_path.xyz session_name '/' subsessions(i).name]);
% $$$            session_parts = session_parts(3:end);
% $$$            for j = 1:length(session_parts)
% $$$                system(['ln -s ' target_path.xyz session_name '/' subsessions(i).name '/' session_parts(j).name ' ' paths.xyz session_name '/' subsessions(i).name '/' session_parts(j).name]);
% $$$            end
% $$$        end
% $$$ 
% $$$        nlx_file = dir([target_path.nlx session_name]);
% $$$        nlx_file = nlx_file(3:end);
% $$$        for j = 1:length(nlx_file)
% $$$            system(['ln -s ' target_path.nlx session_name '/' nlx_file(j).name ' ' paths.nlx session_name '/' nlx_file(j).name]);
% $$$        end
% $$$ 
% $$$        analysis_file = dir([target_path.analysis session_name]);
% $$$        analysis_file = analysis_file(3:end);
% $$$        for j = 1:length(analysis_file)
% $$$            system(['ln -s ' target_path.analysis session_name '/' analysis_file(j).name ' ' paths.analysis session_name '/' analysis_file(j).name]);
% $$$        end
% $$$ 
% $$$   end
% $$$ end
