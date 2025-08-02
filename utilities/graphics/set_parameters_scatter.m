function set_parameters_scatter(hax,varargin)

[color,alpha] = DefaultArgs(varargin,{'',0.5});

hax.MarkerFaceAlpha = alpha;
hax.MarkerEdgeAlpha = alpha;

if color,
    hax.MarkerEdgeColor = color;
    hax.MarkerFaceColor = color;
end

