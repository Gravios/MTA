function Par = subsref(Par,S)

ni = numel(S);


if strcmp(S(1).type,'.'),
    if isprop(Par,S(1).subs)
        if numel(S)==1,
            Par = Par.(S(1).subs);
        else
            Par = subsref(Par.(S(1).subs),S(2:end));
        end
    elseif ismethod(Par,S(1).subs)
        if numel(S)==1,
            Par = Par.(S(1).subs);
        else
            Par = Par.(S(1).subs)(S(2).subs{:});
        end
    elseif isfield(Par.data,S(1).subs)
        Par = Par.data.(S(1).subs);
        S(1) = [];
        if numel(S) > 0,
            Par = builtin('subsref',Par,S);
        end
    else        
        Par = builtin('subsref',Par,S);
    end
end

