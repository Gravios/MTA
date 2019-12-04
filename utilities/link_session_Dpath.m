function link_session_Dpath(sessionName,Dpaths)
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
sessionPath = fullfile(Session.path.project,sessionName);
create_directory(sessionPath);

cwd = pwd();
cd(sessionPath);

for field = fieldnames(Dpaths)'
    field = field{1};
    
    create_directory(fullfile(Session.path.project,field));
    create_directory(fullfile(Session.path.project,field,sessionName));
    cd(fullfile(Session.path.project,field,sessionName));

% REMOVE all links if they exist    
    system(['find . -type l -delete'])
% LINK source files recursively
    system(['cp -rs ', fullfile(Dpaths.(field),sessionName),' ',...
                       fullfile(Session.path.project,field)]);    
% LINK interstitial files recursively
    system(['find . -type d -exec mkdir ',fullfile(sessionPath,'{}'),' \;']);
    system(['find . -type l -exec ln -s {} ',fullfile(sessionPath,'{}'),' \;']);
end

cd(cwd);