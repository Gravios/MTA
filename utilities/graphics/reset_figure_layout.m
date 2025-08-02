function reset_figure_layout(hfig)
% function reset_figure_layout(hfig)
%
% restore preserved figure layout
%
% see preserve_figure_layout.m

hfig.Units = hfig.UserData.figureSettings.Units;
hfig.Position = hfig.UserData.figureSettings.Position;

af(@(h,val) set(h,'Units',val{1}), ...
   findobj(hfig,'Type','Axes'),...
   hfig.UserData.figureSettings.axes.Units);

af(@(h,val) set(h,'Position',val{1}), ...
   findobj(hfig,'Type','Axes'),...
   hfig.UserData.figureSettings.axes.Positions);

af(@(h,val) set(h,'PlotBoxAspectRatioMode',val{1}), ...
   findobj(hfig,'Type','Axes'),...
   hfig.UserData.figureSettings.axes.PlotBoxAspectRatioMode);

