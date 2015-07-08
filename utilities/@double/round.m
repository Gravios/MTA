function out = round(x,varargin)
%ROUND  Round towards nearest integer.
%   ROUND(X,sig) rounds the elements of X to the nearest integers
%   or significant digit.
%
%   See also FLOOR, CEIL, FIX.
%

if isempty(varargin)
    out = builtin('round',x);
else
    sig = varargin{1};
    out = builtin('round',x*10^sig)*10^-sig;
end

