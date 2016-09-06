ses = MTASession([]);
ds = dir(ses.path.data);
fname = {ds(cellfun(@numel,{ds.name})==13).name}

for f = fname
    
    fprintf('-%s---------------------------------------------------------------\n',f{1})
    fprintf('|                                                                             |\n',f{1})
    fprintf('| Sessions:                                                                   |\n',f{1})

    slist = dir(fullfile(ses.path.data,f{1},'*.ses.*'));
    for s = slist'
        fprintf('|    %s                                            |\n',s.name)
    end

    tlist = dir(fullfile(ses.path.data,f{1},'*.trl.*'));
    for t = tlist'
        fprintf('|    %s                                            |\n',t.name)
    end


    fprintf('-----------------------------------------------------------------------------\n',f{1})
end


