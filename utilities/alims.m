function alims(varargin)
%function alims(varargin)
%
% accepts a Nx2 where N:={1,2,3} array where each row, following the order {x,y,z}, 
% is a axis limit pair 
%

% Find axes or assign gca as the target axis
ax = [];
vind = [];
vind = cellfun(@isa,varargin, repmat({'matlab.graphics.axis.Axes'},size(varargin)));
if ~any(vind), 
    ax = gca; 
else
    ax = varargin{vind};
end

% Find axes limits or throw a warning and exit
lims = [];
vind = [];
vind = cellfun(@isa,varargin,repmat({'double'},size(varargin)));
if ~any(vind),
    warning('MTA:utilities:alims:NoLimitsProvided');
    return;
else
    lims = varargin{vind};
    assert(size(lims,2)==2,'MTA:utilities:alims:BadDimSize');
    assert(size(lims,1)<=3,'MTA:utilities:alims:BadDimSize');
end

                
for a = 1:size(lims,1)
    if a==1,
        xlim(ax,lims(a,:));
    elseif a==2,
        ylim(ax,lims(a,:));
    elseif a==3,
        zlim(ax,lims(a,:));
    end
end
