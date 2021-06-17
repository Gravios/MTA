function preserve_figure_layout(hfig)
% function preserve_figure_layout(hfig)
%
% preserve
% 
% see reset_figure_layout.m and set_figure_layout.m

hfig.UserData.figureSettings.Units = hfig.Units;
hfig.UserData.figureSettings.Position = hfig.Position;
hfig.UserData.figureSettings.axes.Units = get(findobj(hfig,'Type','Axes'),'Units');
hfig.UserData.figureSettings.axes.Positions = get(findobj(hfig,'Type','Axes'),'Position');

hfig.UserData.figureSettings.axes.PlotBoxAspectRatioMode ...
    = get(findobj(hfig,'Type','Axes'),'PlotBoxAspectRatioMode');