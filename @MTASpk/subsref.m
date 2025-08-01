function [out,varargout] = subsref(Data,S)
% function Data = subsref(Data,S)
% 
% ref type: Data(clu) clu:isnumeric - returns res of units belonging to cluster clu
%
varargout = {};
ni = numel(S);
sampleRate = Data.sampleRate;
if strcmp(S(1).type,'()') && ni==1,
    if numel(S.subs)==0,
        out = Data.res;
    else
        cid = ismember(Data.clu,S.subs{1});
        out = Data.res(cid);
        if nargout > 1
            varargout{end+1} = Data.fet(cid,:);
        elseif nargout > 2
            varargout{end+1} = Data.spk(cid,:,:);
        end
    end
    if numel(S.subs) == 2
        state = S.subs{2};
        ind = within_ranges(out,state.data);
        out = out(ind);
        if nargout > 1
            varargout{1} = varargout{1}(ind,:)
        elseif nargout > 2
            varargout{2} = varargout{2}(ind,:,:)            
        end
    end    
    return
end
out = builtin('subsref',Data,S);

