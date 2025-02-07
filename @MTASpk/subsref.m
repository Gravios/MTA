
function Data = subsref(Data,S)
% function Data = subsref(Data,S)
% 
% ref type: Data(clu) clu:isnumeric - returns res of units belonging to cluster clu
%
ni = numel(S);
sampleRate = Data.sampleRate;
if strcmp(S(1).type,'()') && ni==1,
    if numel(S.subs)==0,
        Data = Data.res;
    else
        Data = Data.res(ismember(Data.clu,S.subs{1}));
    end
    if numel(S.subs) == 2
        state = S.subs{2};
        Data = Data(within_ranges(Data,state.data));
    end    
    return
end
Data = builtin('subsref',Data,S);

