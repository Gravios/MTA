function R = JoinRanges(varargin)
%function R = JoinRanges(varargin)
%
% made for integer ranges of reasonable length
% don't you dare stick a decimal point in there

rc = {};
if nargin==1
rc = varargin{1};
else 
rc = varargin;
end

rmin = min(cellfun(@min,cellfun(@min,rc,'UniformOutput',false)));
rmax = max(cellfun(@max,cellfun(@max,rc,'UniformOutput',false)));
rind = false(rmax-rmin,1);

shift = rmin-1;
for i = 1:numel(rc)
    for j = 1:size(rc{i})        
        rind(rc{i}(j,1)-shift:rc{i}(j,2)-shift) = true;
    end
end
rind = diff([0;rind;0]);
nper = reshape([find(rind==1)+1,find(rind==-1)],[],2);
R = nper + shift-1;



