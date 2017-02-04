set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)

set(gcf,'Units','centimeters',...
        'PaperPositionMode', 'auto')

hax = axes('Units','centimeters','Position',[]);

print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
