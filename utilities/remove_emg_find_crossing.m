function [val,id] = remove_emg_find_crossing(x,d)
if d>0 % finding from the left to the right side
    id = find(abs(x)<1e-4,1,'first');
else % from the right to the left side
    id = find(abs(x)<1e-4,1,'last');
end
% if none of them get close enough
if isempty(id)
    [~,id] = min(abs(x));
end
val = x(id);
