function listSessions(varargin)
files = dir(MTASession().path.data);
re = '^[a-zA-Z]{1,2}[0-9]{2}[-][0-9]{8,8}$';% regular expression to match Session naming convention
{files(~cellfun(@isempty,regexp({files.name},re))).name}
