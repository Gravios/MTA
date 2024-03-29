function listSessions(pattern)
if ~exist('pattern','var'),
    pattern = '.*ses.*';
end

ses = MTASession([]);
ds = dir(fullfile(ses.path.project,pattern));
fname = {ds(cellfun(@numel,{ds.name})==13).name}

for f = fname
    
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


