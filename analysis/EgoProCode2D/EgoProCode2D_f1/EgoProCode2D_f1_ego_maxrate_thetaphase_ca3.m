% EgoProCode2D_f1_ego_maxrate_thetaphase_ca1


[pfig, sax] =  setup_figure_(mfilename());

% SUBPLOT <- theta cycle (color=partions)
plot(sax,-cos(linspace(0,2*pi)),linspace(0,1),'k','LineWidth',2);
plims = [0.5,              2.26106176905986];
plot(sax,...
     -cos(linspace(plims(1),plims(2))),...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(1,:), 'LineWidth',2);
plims = [2.26106176905986, 4.02212353811972];
plot(sax,...
     -cos(linspace(plims(1),plims(2))),...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(2,:), 'LineWidth',2);
plims = [4.02212353811972, 5.78318530717959];
plot(sax,...
     -cos(linspace(plims(1),plims(2))),...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(3,:),'LineWidth',2);

% FORMAT subplot
ylim([0,1]);
sax.Visible = 'off';
text(0.5,0.05,'0','Rotation',90,'FontSize',8);
text(-0.5,0.5,'\pi','Rotation',90,'FontSize',8);
text(0.5,0.95,'2\pi','Rotation',90,'FontSize',8);


savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

