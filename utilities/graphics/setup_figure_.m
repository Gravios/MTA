function [hfig,sax] = setup_figure_(figName)
    hfig = figure();
    sax = axes();
    sax.Tag = figName;
    hold(sax,'on');
end
