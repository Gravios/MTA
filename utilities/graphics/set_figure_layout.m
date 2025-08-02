function [hfig,figOpts,fax, sax] = set_figure_layout(varargin)
%function [hfig,figOpts,fax, sax] = set_figure_layout(varargin)
% IN:
%     hfig   - Handle: Main figure graphical handle
%     format - String: page format option {'A4'}
%     layout - String: page layout option {'portrait', 'landscape'}
%     units  - String: figure's unit of measure {'normalized','centimeters','inches'...}
%     subplotWidth  - Numeric: 
%     subplotHeight - Numeric: 
%     subplotPaddingHorizontal - Numeric: horizontal spacing of subplots within the same vertical oridnate
%     subplotPaddingHorizontal - Numeric: horizontal spacing of subplots within the same vertical oridnate
%
% OUT:
%     hfig - handle: Main figure graphical handle
%
%     figOpts - struct: Figure layout parameters 
%        figOpts.page.units                
%        figOpts.page.PaperPositionMode    
%        figOpts.page.marginLeft           
%        figOpts.page.marginTop            
%        figOpts.subplot.width             
%        figOpts.subplot.height            
%        figOpts.subplot.horizontalPadding 
%        figOpts.subplot.verticalPadding   
%      
%     fax - handle: Background axes handle
%
%     sax - handle: empty subplot axes handle


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('hfig',                          [],                                            ...
                 'format',                        [],                                            ...
                 'layout',                        [],                                            ...
                 'units',                         'centimeters',                                 ...
                 'subplotWidth',                  1.15,                                          ...
                 'subplotHeight',                 1.15,                                          ...
                 'subplotPaddingHorizontal',      0.0,                                           ...
                 'subplotVerticalPadding',        0.15                                           ...
);
[hfig,format,layout,units,subplotWidth,subplotHeight,                                            ...
 subplotPaddingHorizontal,subplotVerticalPadding] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------       

 
% MAIN --------------------------------------------------------------------------------------------- 

hfig.RendererMode = 'manual';
hfig.Renderer     = 'Painters';
hfig.PaperOrientation = layout;

switch format,
    
  case '1080p'
    units = 'Pixels'
    figOpts.page.width  = 1920;
    figOpts.page.height = 1080;
    figOpts.page.units                = units;
    figOpts.page.PaperPositionMode    = 'auto';
    figOpts.page.marginLeft           = 100;
    figOpts.page.marginTop            = 100;
    figOpts.subplot.units             = units;
    figOpts.subplot.width             = subplotWidth;
    figOpts.subplot.height            = subplotHeight;
    figOpts.subplot.horizontalPadding = subplotPaddingHorizontal;
    figOpts.subplot.verticalPadding   = subplotVerticalPadding;

    
    
  case 'A4'
    hfig.PaperType    = format;
    
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
    figOpts.subplot.units             = units;    
    figOpts.subplot.width             = subplotWidth;
    figOpts.subplot.height            = subplotHeight;
    figOpts.subplot.horizontalPadding = subplotPaddingHorizontal;
    figOpts.subplot.verticalPadding   = subplotVerticalPadding;

end
figOpts.page.xpos = figOpts.page.marginLeft                                                  ...
                    : (figOpts.subplot.width  + figOpts.subplot.horizontalPadding)           ...
                    : figOpts.page.width;
figOpts.page.ypos = fliplr((figOpts.page.height - figOpts.page.marginTop)                    ...
                           -(figOpts.subplot.height + figOpts.subplot.verticalPadding)       ...
                           .*floor((figOpts.page.height - figOpts.page.marginTop)            ...
                                ./(figOpts.subplot.height + figOpts.subplot.verticalPadding))...
                    : figOpts.subplot.height + figOpts.subplot.verticalPadding               ...
                    : figOpts.page.height - figOpts.page.marginTop - figOpts.subplot.height);


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
fax = axes('Position',[0,0,1,1],'Visible','off','Units',units);
xlim([0,hfig.Position(3)]);
ylim([0,hfig.Position(4)]);

% SETUP subplot axes
sax = gobjects([1,0]); 

% END MAIN -----------------------------------------------------------------------------------------