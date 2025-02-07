function sax = setup_axes(hfig, yind, yoffset, xind, xoffset, varargin)

% DEFARGS -----------------------------------------------------------------------
defargs = struct('gYoffset' , 0,                                              ...
                 'gXoffset' , 0,                                              ...
                   'hscale' , 1,                                              ...
                   'wscale' , 1,                                              ...
           'axesFuncHandle' , @axes                                           ...
);
[gYoffset, gXoffset,hscale, wscale, axesFuncHandle ] =                            ...
    DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------

% MAIN --------------------------------------------------------------------------

global DELETE_CURRENT_AXES
LOG_PREFIX = 'MTA:utilities:graphics:setup_axes';

if DELETE_CURRENT_AXES & ~isempty(hfig.UserData.sax)
    warning(LOG_PREFIX,'DELETE_CURRENT_AXES is set to True');
    axI = find(hfig.UserData.sax==hfig.CurrentAxes);
    if ~isempty(axI)
        delete(hfig.UserData.sax(axI));
        hfig.UserData.sax(axI) = [];
    end
end

hash = DataHash({yind, yoffset,                                                 ...
                 xind, xoffset,                                                 ...
                 gYoffset, gXoffset,                                            ...
                 hscale, wscale});

sax = findobj(hfig, 'Tag', hash);

if isempty(sax)
    hfig.UserData.sax(end+1) =                                                  ...
        axesFuncHandle('Units','centimeters',                                   ...
             'Position',[hfig.UserData.page.xpos(xind) + xoffset + gXoffset,    ...
                         hfig.UserData.page.ypos(yind) + yoffset + gYoffset,    ...
                         hfig.UserData.subplot.width  * wscale,                 ...
                         hfig.UserData.subplot.height * hscale],                ...
                      'FontSize', 8,                                            ...
                      'LineWidth',1);
    hold(hfig.UserData.sax(end),'on');
    hfig.UserData.sax(end).Tag = hash;
    sax = hfig.UserData.sax(end);
    
else
    axes(sax);
    cla(sax);
end

% END MAIN -----------------------------------------------------------------------------------------





