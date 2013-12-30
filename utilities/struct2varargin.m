function out = struct2varargin(S)
f = fieldnames(S);
c = struct2cell(S);
out = reshape(cat(1,f',c'),1,numel(c)*2);
end