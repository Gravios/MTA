function Data = load(Data,varargin)
% function Data = load(Data,varargin)
% 
% Loads an MTADepoch object from the location in its path property with the
% filename found in its filename property
    [Session,sync] = DefaultArgs(varargin,{[],[]});

    if ~isempty(Data.filename)
        load(Data.fpath);
    else
        files = dir(Data.path);
        re = '\.sst\.';
        stsFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};                
        if ~isempty(Data.key)
            re = ['\.' Data.key '\.'];
        elseif ~isempty(Data.label)
            re = ['\.' Data.label '\.'];                   
        else
            error('STS file does not exist or key/label is missing')
        end
        Data.filename = stsFileList{~cellfun(@isempty,regexp(stsFileList,re))};
        load(Data.fpath)
    end
    if ~isempty(Session),
        Data.resync(Session);
    end
    if ~isempty(sync),
        Data.resync([],sync);
    end

end
