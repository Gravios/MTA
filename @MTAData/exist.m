function out = exist(Data)
out = exist(fullfile(Data.path,Data.filename),'file');