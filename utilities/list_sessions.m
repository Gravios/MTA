function list_sessions(sessionName)
if ~exist('sessionName','var'),
    sessionName = '';
end

% RETRIVE file list of project directory
ses = MTASession([]);
files = dir(fullfile(ses.path.project,sessionName));

re = '^[a-zA-Z]{1,2}[0-9]{2,4}[-][0-9]{8,8}$'; % session naming convention

% SELECT directories matching session naming convention
sessionList = {files(~cellfun(@isempty,regexp({files.name},re)) & [files(:).isdir]).name};
 
for f = sessionList
    
    fprintf('-%s---------------------------------------------------------------\n',f{1})
    fprintf('|                                                                             |\n')
    fprintf('| Sessions:                                                                   |\n')

    slist = dir(fullfile(ses.path.project,f{1},'*.ses.*'));
    for s = slist'
        fprintf('|    %s                                            |\n',s.name)
    end

    fprintf('|                                                                             |\n')
    fprintf('| Trials:                                                                     |\n')

    tlist = dir(fullfile(ses.path.project,f{1},'*.trl.*'));
    for t = tlist'
        fprintf('|    %s                                            |\n',t.name)
    end


    fprintf('-----------------------------------------------------------------------------\n',f{1})
end


