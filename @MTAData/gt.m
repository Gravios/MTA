function out = gt(a,b)
% function out = gt(a,b)    
% Greater than wrapper for MTAData
if  isa(a,'MTAData') &&  isa(b,'MTAData'), out = a.data>b.data; return,end
if ~isa(a,'MTAData') &&  isa(b,'MTAData'), out = a     >b.data; return,end
if  isa(a,'MTAData') && ~isa(b,'MTAData'), out = a.data>b;      return,end
end
