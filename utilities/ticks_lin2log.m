function ticks_lin2log(varargin)
% ticks_lin2log(ax,options)
% DefaultArgs
% ax: gca
% options: 'xy'
%
[ax,options] = DefaultArgs(varargin,{gca,'xy'});

if ismember('x',options),
xticks = get(ax,'XTickLabel');
xticks = mat2cell(xticks,ones(1,size(xticks,1)),size(xticks,2));
xticks = cellfun(@str2num,xticks);
xticks = mat2cell(10.^xticks,ones(1,size(xticks,1)),1);
xticks = cellfun(@sprintf,repmat({'%2.2f'},numel(xticks),1),xticks,'uniformoutput',false);
set(ax,'XTickLabel',xticks);
end

if ismember('y',options),
yticks = get(ax,'YTickLabel');
yticks = mat2cell(yticks,ones(1,size(yticks,1)),size(yticks,2));
yticks = cellfun(@str2num,yticks);
yticks = mat2cell(10.^yticks,ones(1,size(yticks,1)),1);
yticks = cellfun(@sprintf,repmat({'%2.2f'},numel(yticks),1),yticks,'uniformoutput',false);
set(ax,'YTickLabel',yticks);
end


if ismember('z',options),
zticks = get(ax,'ZTickLabel');
zticks = mat2cell(zticks,ones(1,size(zticks,1)),size(zticks,2));
zticks = cellfun(@str2num,zticks);
zticks = mat2cell(10.^zticks,ones(1,size(zticks,1)),1);
zticks = cellfun(@sprintf,repmat({'%2.2f'},numel(zticks),1),zticks,'uniformoutput',false);
set(ax,'zTickLabel',zticks);
end