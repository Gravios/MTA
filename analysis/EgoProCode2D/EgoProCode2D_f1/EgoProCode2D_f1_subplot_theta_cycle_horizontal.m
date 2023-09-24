% EgoProCode2D_f1_subplot_theta_cycle_horizontal

[pfig, sax] = setup_figure_(mfilename());

% SUBPLOT <- theta cycle (color=partions)
plot(sax,linspace(0,1),cos(linspace(0,2*pi)),'k','LineWidth',2);
plims = [0.5,              2.26106176905986];
plot(sax,...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)),...
     cos(linspace(plims(1),plims(2))), 'Color',pclr(1,:), 'LineWidth',2);
plims = [2.26106176905986, 4.02212353811972];
plot(sax,...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)),...
     cos(linspace(plims(1),plims(2))), 'Color',pclr(2,:), 'LineWidth',2);
plims = [4.02212353811972, 5.78318530717959];
plot(sax,...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)),...
     cos(linspace(plims(1),plims(2))), 'Color',pclr(3,:),'LineWidth',2);


%FORMAT subplot
ylim(sax,[-1,1]);
xlim(sax,[0,1]);
sax.Visible = 'off';
text(sax,0.05,0.5,'0','FontSize',8);
text(sax,0.5,0.5,'\pi','FontSize',8);
text(sax,0.95,0.5,'2\pi','FontSize',8);

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
close(pfig);
