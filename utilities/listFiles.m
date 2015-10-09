function listFiles(SessionName,fileExt)
%function listFiles(SessionName,fileExt)
files = dir(fullfile(MTASession().path.data,SessionName));
re = ['\.' fileExt '\.'];% regular expression to match Session naming convention
{files(~cellfun(@isempty,regexp({files.name},re))).name}
