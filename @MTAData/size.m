function sdim = size(Data,varargin)
%sdim = size(Data,dimInd)
% Returns the size of the MTAData object's "data" field
%
% Inputs:
%   dimInd - numericArray: indicies of the dimensions. If v 
%                            empty all dimesion sizes are
%                            returned.
%
% Outputs:
%   sdim - numericArray: contains the size of each dimension
%

if ~isempty(cell2mat(varargin)),
    dim = cell2mat(varargin);
    ndim = numel(dim);
    if ndim>1,
        sdim = zeros(1,ndim);
        for i = 1:ndim,
            sdim(i) = size(Data.data,dim(i));
        end
    else
        sdim = builtin('size',Data.data,dim);
    end
else
    sdim = size(Data.data);
end
end

