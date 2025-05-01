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

        if ~isempty(Data.key)
            re = ['.' Data.key '.'];
        elseif ~isempty(Data.label)
            re = ['.' Data.label '.'];
        end
        
        % GET the list of all sst files
        re = regexptranslate('escape', re);
        matchedList = {files(~cellfun(@isempty,regexp({files.name}, re))).name};
        
        if ~isempty(matchedList),
            if numel(matchedList)>1
                nre = [ Data.parent '.' Data.ext ];
                nre = regexptranslate('escape', nre);
                Data.filename = matchedList{~cellfun(@isempty, regexp(matchedList, nre))};
            else
                Data.filename = matchedList{1};
            end
            
        else,
            error(struct('identifier','MTA:MTAData:load:FileNotFound',...
                         'message',   'MTAData object file not found',...
                         'stack',     dbstack));
        end

        load(Data.fpath)
    end
    if ~isempty(Session),
        Data.resync(Session);
    end
    if ~isempty(sync),
        Data.resync([],sync);
    end

end
