function out = ne(a,b)
% function out = ne(a,b)
% Not equal. Wapper for MTAData.
if  isa(a,'MTAData') &&  isa(b,'MTAData'), out = a.data~=b.data; return,end
if ~isa(a,'MTAData') &&  isa(b,'MTAData'), out = a     ~=b.data; return,end
if  isa(a,'MTAData') && ~isa(b,'MTAData'), out = a.data~=b;      return,end
end
