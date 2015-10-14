function out = struct2varargin(S)
%function out = struct2varargin(S)
%
% Converts struct to parameter name/value pairs .
%
% See also:  DefaultArgs
    f = fieldnames(S);
    c = struct2cell(S);
    out = reshape(cat(1,f',c'),1,numel(c)*2);
end