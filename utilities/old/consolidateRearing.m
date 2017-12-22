function consolidateRearing(Session,clean_up_files)

bfccg = [];

files = dir(Session.spath.analysis);
re = ['\.REARING\.'];
rearFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};
for i = 1:length(rearFileList),
    if rearFileList{i},
        ds = load([Session.spath.analysis rearFileList{i}]);
        if i==1, 
            bfccg=ds.bfccg;
        else
            bfccg = cat(1,bfccg,ds.bfccg);
        end
        if clean_up_files,
            system(['rm ' Session.spath.analysis rearFileList{i}]);
        end
    end
end
