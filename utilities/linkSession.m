function linkSession(session_name,xyz_path,nlx_path)
%function linkSession(session_name,xyz_path,nlx_path)
% TODO make linkSession windows compatible 
    
%need to make cygwin system('C:\cygwin64\bin\bash -c CYGWIN=winsymlinks:native;/bin/ls')    
%system('C:\cygwin64\bin\bash -c ''/bin/ln -s /cygdrive/c/Users/justi_000/Documents/MATLAB/dlmtest.txt  /cygdrive/c/Users/justi_000/Documents/MATLAB/dlmt.txt''')
%system(['C:\cygwin64\bin\bash -c ''export CYGWIN64="winsymlinks:native";/bin/ln -s c:/Users/justi_000/Documents/MATLAB/dlmtest.txt  `c:/Users/justi_000/Documents/MATLAB/dlmt.txt;c:/Users/justi_000/Documents/MATLAB/dlmtest.txt'' -> `c:/Users/justi_000/Documents/MATLAB/dlmt.txt'''''])


%system(['C:\cygwin64\bin\bash -c export CYGWIN="winsymlinks:native";"/bin/ln -s -v c:/Users/justi_000/Documents/MATLAB/dlmtest.txt";`/cygdrive/c/Users/justi_000/Documents/MATLAB/dlmtest.txt -> `/cygdrive/c/Users/justi_000/Documents/MATLAB/dlmt.txt"'])
%system(['C:\cygwin64\bin\bash -c export CYGWIN="winsymlinks:native";"/bin/ln -s -v /cygdrive/c/Users/justi_000/Documents/MATLAB/dlmtest.txt";`/cygdrive/c/Users/justi_000/Documents/MATLAB/dlmtest.txt -> `/cygdrive/c/Users/justi_000/Documents/MATLAB/dlmt.txt"'])
%system(['C:\cygwin64\bin\bash -c export CYGWIN="winsymlinks:native";"/bin/ln -s -v c:/Users/justi_000/Documents/MATLAB/dlmtest.txt";`C:\Users\justi_000\Documents\MATLAB\dlmtest.txt -> `C:\Users\justi_000\Documents\MATLAB\dlmt.txt"'])
%system(['C:\cygwin64\bin\bash -c export CYGWIN="winsymlinks:native";"/bin/ln -s -v c:/Users/justi_000/Documents/MATLAB/dlmtest.txt c:/Users/justi_000/Documents/MATLAB/dlmt.txt";`C:\Users\justi_000\Documents\MATLAB\dlmtest.txt'' -> `C:\Users\justi_000\Documents\MATLAB\dlmt.txt'''])

    Session = MTASession([]);
    datDir = fullfile(Session.path.data,session_name);
    try,mkdir(datDir);end

    if ~isempty(xyz_path),
        xyz_maze_dirs = dir(fullfile(xyz_path,session_name));xyz_maze_dirs(1:2)=[];
   
        if ~strcmp(xyz_path,fullfile(Session.path.data,'xyz')),
            xyzDir = fullfile(Session.path.data,'xyz',session_name);
            try,mkdir(xyzDir);end
            
            for maze = 1:numel(xyz_maze_dirs),
                try,mkdir(fullfile(xyzDir,xyz_maze_dirs(maze).name));end
                system(['ln -sif ' fullfile(xyz_path,session_name,xyz_maze_dirs(maze).name) '/* ' fullfile(xyzDir,xyz_maze_dirs(maze).name)]);
            end
        else
            
            xyzDir = fullfile(xyz_path,session_name);
        end
    end

    
    %% Link Nlx data to nlx folder
    if ~isempty(nlx_path),
        if ~strcmp(nlx_path,fullfile(Session.path.data,'nlx')),
            nlxDir = fullfile(Session.path.data,'nlx',session_name);
            mkdir(nlxDir);
            if numel(dir(fullfile(nlx_path,session_name)))>2,
                system(['ln -sif ' fullfile(nlx_path,session_name) '/* ' nlxDir]);
            end
        else
            % if nlx_path is within the MTA data collection
            nlxDir = fullfile(nlx_path,session_name);
        end
        cd(datDir);
        system(['ln -sif ../nlx/' session_name '/* ' datDir ])
        
    end
    
    
    if ~isempty(xyz_path),
        for maze = 1:numel(xyz_maze_dirs),
            cd(datDir);
            try,mkdir(xyz_maze_dirs(maze).name);end
            
            try
                system(['ln -sif ../xyz/' session_name '/' xyz_maze_dirs(maze).name '/' ...
                    session_name '-' xyz_maze_dirs(maze).name '.vsk ' ...
                    datDir '/' session_name '-' xyz_maze_dirs(maze).name '.vsk ']);
            end
            cd(fullfile(datDir,xyz_maze_dirs(maze).name))
            try,
                system(['ln -sif ../../xyz/' session_name '/' ...
                    xyz_maze_dirs(maze).name '/* ' datDir '/' ...
                    xyz_maze_dirs(maze).name '/']);
            end
            
        end
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
