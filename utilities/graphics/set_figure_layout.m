function [hfig,figOpts,fax] = set_figure_layout(varargin)
% hfig,
% format
% layout

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('hfig',                          [],                                            ...
                 'format',                        [],                                            ...
                 'layout',                        [],                                            ...
                 'units',                         'centimeters',                                 ...
                 'subplotWidth',                  1.15,                                          ...
                 'subplotHeight',                 1.15,                                          ...
                 'subplotHorizontalPadding',      0.0,                                           ...
                 'subplotVerticalPadding',        0.15                                           ...
);
[hfig,format,layout,units,subplotWidth,subplotHeight,                                            ...
 subplotHorizontalPadding,subplotVerticalPadding] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------       

 
% MAIN --------------------------------------------------------------------------------------------- 
 switch format,
    
  case 'A4'
    switch layout
      case 'portrait'
        figOpts.page.width    = 21.0;
        figOpts.page.height   = 29.7;
      case 'landscape'
        figOpts.page.width    = 29.7;
        figOpts.page.height   = 21.0;
    end
    
    figOpts.page.units                = units;
    figOpts.page.PaperPositionMode    = 'auto';
    figOpts.page.marginLeft           = 2.54;
    figOpts.page.marginTop            = 2.54;
    figOpts.subplot.width             = subplotWidth;
    figOpts.subplot.height            = subplotHeight;
    figOpts.subplot.horizontalPadding = subplotHorizontalPadding;
    figOpts.subplot.verticalPadding   = subplotVerticalPadding;

    figOpts.page.xpos = figOpts.page.marginLeft                                                  ...
                        : (figOpts.subplot.width  + figOpts.subplot.horizontalPadding)           ...
                        : figOpts.page.width;
    figOpts.page.ypos = fliplr((figOpts.page.height - figOpts.page.marginTop)                    ...
                               -(figOpts.subplot.height + figOpts.subplot.verticalPadding)       ...
                               .*floor((figOpts.page.height - figOpts.page.marginTop)            ...
                                       ./(figOpts.subplot.height + figOpts.subplot.verticalPadding))...
                        : figOpts.subplot.height + figOpts.subplot.verticalPadding               ...
                        : figOpts.page.height - figOpts.page.marginTop - figOpts.subplot.height);
end

% SETUP figure
if isempty(hfig), 
    hfig = figure();
end            
clf(hfig);

hfig.Units = figOpts.page.units;
hfig.Position = [1,                               ...
                 1,                               ...
                 figOpts.page.width,                  ...
                 figOpts.page.height];
hfig.PaperPositionMode = figOpts.page.PaperPositionMode;

% SETUP Background axes
fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
xlim([0,hfig.Position(3)]);
ylim([0,hfig.Position(4)]);
 
% END MAIN -----------------------------------------------------------------------------------------