



% FIGURE - Summary of ego-fields: {headBodyAngle X thetaPhase}

% BLOCK - mean ego-field position
xBlockOffset = 0;
yBlockOffset = 0;
ascnPhz = 4;
descPhz = 2;
orientationLabels = 'LLCRR';
[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','portrait',[],4,4,1.5,1.0);
[yind, yOffSet, xind, xOffSet] = deal(1,                         ... yind
                                      0,                         ... yOffSet
                                      1,                         ... xind
                                      0);                  % xOffSet

for a = 1:5,
    [yind, yOffSet, xind, xOffSet] = deal(   a,       0,    1,       0);        
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
                        fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
                        fig.subplot.width,                                ...
                        fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');
    plot(sq(rwpa(~tida,2,descPhz,a)),sq(rwpa(~tida,1,descPhz,a)),'.');        
    plot(mean(rwpa(~tida,2,descPhz,a)),mean(rwpa(~tida,1,descPhz,a)),'*r');        
    xlim([-100,100]);
    ylim([-100,100]);
    daspect([1,1,1]);
    grid('on');
    Lines([],0,'k');
    Lines(0,[],'k');
    set(sax(end),'XTick',[-100,-50,0,50,100]);
    title(['HBA: [',num2str(round(circ_rad2ang(hbaBinEdges(a)))),', ',num2str(round(circ_rad2ang(hbaBinEdges(a+1)))),'] : ',orientationLabels(a)]);
    if a == 1,
        title({'Descending Theta Phase',...
               ['PHZ: [',num2str(circ_rad2ang(binPhzs(descPhz))),', ',num2str(circ_rad2ang(binPhzs(descPhz+1))),']'], ...                       
               ['HBA: [',num2str(circ_rad2ang(hbaBinEdges(a))),', ',num2str(circ_rad2ang(hbaBinEdges(a+1))),'] : ',orientationLabels(a)]});
    end
    if a == 5,
        xlabel({'left  center  right','head FOR (mm)'});
        ylabel({'head FOR (mm)','back  center    front'});
    end
    

    % XY ascending phase                          yind, yOffSet, xind, xOffSet
    [yind, yOffSet, xind, xOffSet] = deal(   a,       0,    2,       0);
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
                        fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
                        fig.subplot.width,                                ...
                        fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');
    plot(  sq(rwpa(~tida,2,ascnPhz,a)),  sq(rwpa(~tida,1,ascnPhz,a)),'.');
    plot(mean(rwpa(~tida,2,ascnPhz,a)),mean(rwpa(~tida,1,ascnPhz,a)),'*r');                
    xlim([-100,100])
    ylim([-100,100])
    daspect([1,1,1]);
    grid('on');
    Lines([],0,'k');
    Lines(0,[],'k');
    set(sax(end),'XTick',[-100,-50,0,50,100]);
    title(['HBA: [',num2str(round(circ_rad2ang(hbaBinEdges(a)))),', ',num2str(round(circ_rad2ang(hbaBinEdges(a+1)))),'] :',orientationLabels(a)]);
    if a == 1,
        title({'Ascending Theta Phase',...
               ['PHZ: [',num2str(circ_rad2ang(binPhzs(ascnPhz))),', ',num2str(circ_rad2ang(binPhzs(ascnPhz+1))),']'], ...
               ['HBA: [',num2str(round(circ_rad2ang(hbaBinEdges(a)))),', ',num2str(round(circ_rad2ang(hbaBinEdges(a+1)))),'] : ',orientationLabels(a)]});                
    end
    
    % PHASE precession vector angle
    [yind, yOffSet, xind, xOffSet] = deal(   a,       0,    3,       0);            
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                        fig.page.ypos(yind)+yOffSet,                      ...
                        fig.subplot.width,                                ...
                        fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    rose(-atan2(diff(sq(rwpa(~tida,2,[2,p],a)),1,2),diff(sq(rwpa(~tida,1,[2,p],a)),1,2)),36);
    daspect([1,1,1]);
    if a==1,
        title({'Angle of phase','precession vector','relative to head direction'});
    end
end

% TITLE of figure
axes(fax)
text(0.5*fig.page.width,...
     0.95*fig.page.height,...
     'Rate Weighted Mean Position in the Egocentric Frame of Reference',...
     'HorizontalAlignment','center',...
     'FontSize', 14);

% END BLOCK - mean ego-field position


% END figure
