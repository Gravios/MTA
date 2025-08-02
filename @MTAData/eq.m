function out = eq(a,b)
% function out = eq(a,b)
% Equal. Wrapper for MTAData
if  isa(a,'MTAData') &&  isa(b,'MTAData'), out = a.data==b.data; return,end
if ~isa(a,'MTAData') &&  isa(b,'MTAData'), out = a     ==b.data; return,end
if  isa(a,'MTAData') && ~isa(b,'MTAData'), out = a.data==b;      return,end
end
