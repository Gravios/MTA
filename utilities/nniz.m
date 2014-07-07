function lind= nniz(var,varargin)
%function lind = nniz(var)
%not nan
%not inf
%not zero
szv = size(var);
if isa(var,'MTAData'),
    lind = var~=0&~var.isnan&~var.isinf;
else
    lind = var~=0&~isnan(var)&~isinf(var);
end
for s = 2:numel(szv),
    lind = sum(lind,s)>0;
end
