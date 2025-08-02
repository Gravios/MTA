function out = lt(a,b)
% function out = lt(a,b)
% Less than wrapper 
if  isa(a,'MTAData') &&  isa(b,'MTAData'), out = a.data<b.data; return,end
if ~isa(a,'MTAData') &&  isa(b,'MTAData'), out = a     <b.data; return,end
if  isa(a,'MTAData') && ~isa(b,'MTAData'), out = a.data<b;      return,end
end
