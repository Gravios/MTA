function varargout = cf(varargin)
%function varargout = cf(varargin)
%Wrapper for cellfun because Uniformity shouldn't be allowed

varargin = cat(2,varargin,{'UniformOutput',false});
if nargout>0
    varargout = cell([1,nargout]);
    [varargout{:}] = cellfun(varargin{:});
else 
    cellfun(varargin{:});
end

    