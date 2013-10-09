function out = end(Data,k,n)
if isa(Data,'MTAData')||isa(Data,'MTASync')||isa(Data,'MTAStateCollection')
    out=builtin('end',Data.data,k,n);
else
    out=builtin('end',Data,k,n);
end
end