function out = floor(x,varargin)
%ROUND  Round towards nearest integer.
%   ROUND(X,sig) rounds the elements of X to the nearest integers
%   or significant digit.
%
%   See also FLOOR, CEIL, FIX.
%

if isempty(varargin)
    out = builtin('floor',x);
else
    sig = varargin{1};
    out = builtin('floor',x*10^sig)*10^-sig;
end

