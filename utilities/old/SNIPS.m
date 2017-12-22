set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)
hfig = figure(figHnum);clf
set(hfig,'units','centimeters')
set(hfig,'Position',[0,0,(4+2)*numel(mode)+2,4*(3+2)+2])
set(hfig,'PaperPositionMode','auto');

axes('Units','centimeters',...
     'Position',[(4+2)*(m-1)+2,(3+2)*4-(4)+1,4,3])

h=colorbar('EastOutside');
apos = get(gca,'Position');
set(h,'Units','centimeters');
set(h,'Position',[apos(1)+4.1,apos(2),0.3,3])
