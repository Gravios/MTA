function linkSession(sessionName,xyzPath,nlxPath)
%function linkSession(sessionName,xyzPath,nlxPath)
%
%  Creates Session folder in the current MTA project
%  and creates symlinks pointed to folders containing
%  session data.
% 
% TODO make linkSession windows compatible 
% UPDATE requires link shell extension
% MOD removed f flag from ln command

Session = MTASession([]);
datDir = fullfile(Session.path.project,sessionName);
try,mkdir(datDir);end

if ~isempty(xyzPath),
% GET directory list along xyzPath
    xyzMazeDirs = dir(fullfile(xyzPath,sessionName));xyzMazeDirs(1:2)=[];
    if ~strcmp(xyzPath,fullfile(Session.path.project,'xyz')),        
        xyzDir = fullfile(Session.path.project,'xyz',sessionName);
        create_directory(xyzDir);
% REMOVE old links
        system(['find ',xyzDir,' -type l -delete']);        
        for maze = 1:numel(xyzMazeDirs),
% CREATE maze directory for xyz data files
            create_directory(fullfile(xyzDir,xyzMazeDirs(maze).name));
            system(['find ',fullfile(xyzDir,xyzMazeDirs(maze).name),' -type l -delete']);
            system(['ln -s ' fullfile(xyzPath,sessionName,xyzMazeDirs(maze).name) '/* ' fullfile(xyzDir,xyzMazeDirs(maze).name)]);
        end
    else
        xyzDir = fullfile(xyzPath,sessionName);
    end
end


% Link Nlx data to nlx folder
if exist('nlxPath','var')
    if ~isempty(nlxPath),
        if ~strcmp(nlxPath,fullfile(Session.path.project,'nlx')),
            nlxDir = fullfile(Session.path.project,'nlx',sessionName);
            create_directory(nlxDir);
            system(['find ',nlxDir,' -type l -delete']);
            if numel(dir(fullfile(nlxPath,sessionName)))>2,
                system(['ln -s ' fullfile(nlxPath,sessionName) '/* ' nlxDir]);
            end
        else
            % if nlxPath is within the MTA data collection
            nlxDir = fullfile(nlxPath,sessionName);
        end
        cd(datDir);
        system(['find ',datDir,' -type l -delete']);        
        system(['ln -s ../nlx/' sessionName '/* ' datDir ])
    end
end


% For each maze in the xyz data link directory 
if exist('xyzPath','var')    
    if ~isempty(xyzPath),
        for maze = 1:numel(xyzMazeDirs),
            cd(datDir);
            try,mkdir(xyzMazeDirs(maze).name);end
            
            try
                system(['ln -sf ../xyz/' sessionName '/' xyzMazeDirs(maze).name '/' ...
                        sessionName '-' xyzMazeDirs(maze).name '.vsk ' ...
                        datDir '/' sessionName '-' xyzMazeDirs(maze).name '.vsk ']);
            end
            cd(fullfile(datDir,xyzMazeDirs(maze).name))
            try,
                system(['find ',datDir,'/',xyzMazeDirs(maze).name,'/ -type l -delete']);
                system(['ln -s ../../xyz/',sessionName,'/',xyzMazeDirs(maze).name,'/* ',datDir,'/',xyzMazeDirs(maze).name,'/']);
            end
            
        end
    end
end


