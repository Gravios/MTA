function sax = setup_axes(fig, yind, yoffset, xind, xoffset, varargin)

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('gYoffset' , 0,                                                                 ...
                 'gXoffset' , 0,                                                                 ...
                   'wscale' , 1,                                                                 ...
                   'hscale' , 1                                                                  ...
);
[gXoffset, gYoffset,wscale, hscale ] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

sax = axes('Units','centimeters',                                                                ...
           'Position',[fig.page.xpos(xind) + xoffset + gXoffset,                                 ...
                       fig.page.ypos(yind) + yoffset + gYoffset,                                 ...
                       fig.subplot.width  * wscale,                                              ...
                       fig.subplot.height * hscale],                                             ...
                  'FontSize', 8,                                                                 ...
                  'LineWidth',1);
hold(sax,'on');

% END MAIN -----------------------------------------------------------------------------------------




