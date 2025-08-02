function out = le(a,b)
% function out = le(a,b)    
% Less than or equal. Wrapper for MTAData
if  isa(a,'MTAData') &&  isa(b,'MTAData'), out = a.data<=b.data; return,end
if ~isa(a,'MTAData') &&  isa(b,'MTAData'), out = a     <=b.data; return,end
if  isa(a,'MTAData') && ~isa(b,'MTAData'), out = a.data<=b;      return,end
end
