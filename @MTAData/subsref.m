function Data = subsref(Data,S)
ni = numel(S);
switch Data.type
  case 'TimeSeries',
    if strcmp(S(1).type,'()')&&ni==1,
        modelRef = (cellfun(@numel,S(1).subs(2:end))>1).*cellfun(@ischar,S(1).subs(2:end))|cellfun(@iscell,S(1).subs(2:end));
        if sum(modelRef)>0&&~isempty(Data.model),
            mri = find(modelRef)+1;
            for  i = mri,
                S(1).subs{i} = Data.model.gmi(S(1).subs{i});
            end
        end
        
        % Convert MTAepoch to nx2 matrix and resample if necessary
        epicDataInd = find(cellfun(@isa,S(1).subs,repmat({'MTADepoch'},1,numel(S(1).subs))));
        for i = epicDataInd,
            epoch = S(1).subs{i}.copy;
            if Data.sampleRate~=epoch.sampleRate,
                epoch.resample(Data);
                epoch.clean;
            end
            S(1).subs{i} = epoch.data;
        end
        
        if numel(S.subs)==0,
            Data = Data.data;
        elseif all(cellfun(@islogical,S.subs)+cellfun(@isnumeric,S.subs)+...
                   cellfun(@strcmp,S.subs,repmat({':'},[1,numel(S.subs)]))),
            if strcmp(S.subs{1},':'),
                Data = builtin('subsref',Data.data,S);
            elseif size(S.subs{1},1)==1&&size(S.subs{1},2)>2,
                Data = builtin('subsref',Data.data,S);
            elseif size(S.subs{1},2)==1,
                Data = builtin('subsref',Data.data,S);
            else
                S.subs = S.subs(~cellfun(@isempty,S.subs));
                Sa = S;
                Sa.subs{1} = ':';
                Data = SelectPeriods(builtin('subsref',Data.data,Sa),S.subs{1},'c');
            end
        else
            S.subs = S.subs(~cellfun(@isempty,S.subs));
            Sa = S;
            Sa.subs{1} = ':';
            Data = SelectPeriods(builtin('subsref',Data.data,Sa),S.subs{1},'c');
        end
        return
    end
  case 'TimePeriods'
    if strcmp(S(1).type,'()')&&ni==1,
        if numel(S.subs)==0,
            Data = Data.data;
        elseif numel(S.subs)==1,
            if  strcmp(S.subs{1},':')
                Data = Data.data(:);
            else
                Data = Data.data(S.subs{1});
            end
        else
            Data = Data.data(S.subs{1},S.subs{2});
        end
        return
    end
end
if strcmp(S(1).type,'.'),
    if isprop(Data,S(1).subs)
        if numel(S)==1,
            Data = Data.(S(1).subs);
        else
            Data = subsref(Data.(S(1).subs),S(2:end));
        end
    elseif ismethod(Data,S(1).subs)
        if numel(S)==1,
            Data = Data.(S(1).subs);
        else
            Data = Data.(S(1).subs)(S(2).subs{:});
        end
    else
        Data = builtin('subsref',Data,S);
    end
end

end
