function varargout = af(varargin)
%function varargout = af(varargin)
%Wrapper for arrayfun because Uniformity shouldn't be allowed

varargin = cat(2,varargin,{'UniformOutput',false});
if nargout>0
    varargout = cell([1,nargout]);
    [varargout{:}] = arrayfun(varargin{:});
else 
    arrayfun(varargin{:});
end

    