function flist = list_files(SessionName,varargin)
%function listFiles(SessionName,fileExt)
if ~isempty(varargin),    fileExt = varargin{1};     varargin(1) = [];    
else,                     fileExt = '.*';            end;   
if ~isempty(varargin),    returnType =  varargin{1}; varargin(1) = [];
else,                     returnType = 'cell';       end;


files = dir(fullfile(MTASession().path.data,SessionName));
re = regexptranslate('escape',fileExt);% regular expression to match Session naming convention
flist = {files(~cellfun(@isempty,regexp({files.name},re))).name};
switch returnType,
  case 'string'
    flist = strjoin(flist,' ');
end

    
