function R = JoinRanges(R1,R2)
%function R = JoinRanges(R1,R2)
%
% Joins ranges R2 and R1 : R = R1+R2
%


NotR2 = [[-inf, R2(1,1)];[R2(1:end-1,2), R2(2:end,1)]; [R2(end,2), inf]];
NotR1 = [[-inf, R1(1,1)];[R1(1:end-1,2), R1(2:end,1)]; [R1(end,2), inf]];

R = IntersectRanges(NotR1, NotR2)';
R = reshape(R(2:end-1),2,[])';

% $$$ rc = {};
% $$$ if nargin==1
% $$$ rc = varargin{1};
% $$$ else 
% $$$ rc = varargin;
% $$$ end
% $$$ 
% $$$ rmin = min(cellfun(@min,cellfun(@min,rc,'UniformOutput',false)));
% $$$ rmax = max(cellfun(@max,cellfun(@max,rc,'UniformOutput',false)));
% $$$ rind = false(rmax-rmin,1);
% $$$ 
% $$$ shift = rmin-1;
% $$$ for i = 1:numel(rc)
% $$$     for j = 1:size(rc{i})        
% $$$         rind(rc{i}(j,1)-shift:rc{i}(j,2)-shift) = true;
% $$$     end
% $$$ end
% $$$ rind = diff([0;rind;0]);
% $$$ nper = reshape([find(rind==1)+1,find(rind==-1)],[],2);
% $$$ R = nper + shift-1;



