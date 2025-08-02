function set_figure_layout_mode(hfig,mode)

switch mode
  case 'absolute'
    if ~isempty(hfig.UserData)
        reset_figure_layout(hfig);
    else
        warning('MTA:utilities:graphics:set_figure_layout_mode:NothingDone');
    end
    
  case 'normalized'
    set(findobj(hfig,'Type','Axes'),'PlotBoxAspectRatioMode','auto');        
    set(findobj(hfig,'Type','Axes'),'Units','normalized');
    
  case 'fixedaspect'
    set(findobj(hfig,'Type','Axes'),'PlotBoxAspectRatioMode','manual');    
    set(findobj(hfig,'Type','Axes'),'Units','normalized');    
    
end    

