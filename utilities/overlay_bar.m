function overlay_bar(bax,varargin)
[color,alpha] = DefaultArgs(varargin,{'r',.4},true);
bax.FaceColor = color;
bax.FaceAlpha = alpha;
bax.EdgeColor = color;
bax.EdgeAlpha = alpha;