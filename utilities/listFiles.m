function flist = listFiles(SessionName,varargin)
%function listFiles(SessionName,fileExt)
if ~isempty(varargin)
    fileExt = varargin{1};
else
    fileExt = '.*';
end

files = dir(fullfile(MTASession().path.data,SessionName));
re = regexptranslate('escape',fileExt);% regular expression to match Session naming convention
flist = {files(~cellfun(@isempty,regexp({files.name},re))).name};
