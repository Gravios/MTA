
hfig = figure;

u = 1;
set(hfig,'Name',num2str(u));
while u~=-1,

subplot(251);
hist(mywu(:,u))
xl = xlim
xlim([-0.05,xl(2)])

subplot(252);
imagescnan({xbins,ybins,pfknnmrw(:,:,u)},[],[],0,[0,0,0]),axis xy
title(['walk u ' num2str(Session.spk.map(u,:))])
text(xbins(1)+30,ybins(1)-50,sprintf('%2.1f',max(max(pfknnmrw(:,:,u)))),'Color','w','FontWeight','bold','FontSize',10)

subplot(253);
imagescnan({xbins,ybins,pfknnmrr(:,:,u)},[],[],0,[0,0,0]),axis xy
title(['rear u ' num2str(Session.spk.map(u,:))])
text(xbins(1)+30,ybins(1)-50,sprintf('%2.1f',max(max(pfknnmrr(:,:,u)))),'Color','w','FontWeight','bold','FontSize',10)

subplot(254);
imagescnan({xbins,ybins,pfknnmrg(:,:,u)},[],[],0,[0,0,0]),axis xy
title(['hwalk u ' num2str(Session.spk.map(u,:))])
text(xbins(1)+30,ybins(1)-50,sprintf('%2.1f',max(max(pfknnmrg(:,:,u)))),'Color','w','FontWeight','bold','FontSize',10)

subplot(255);
imagescnan({xbins,ybins,pfknnmrl(:,:,u)},[],[],0,[0,0,0]),axis xy
title(['lwalk u ' num2str(Session.spk.map(u,:))])
text(xbins(1)+30,ybins(1)-50,sprintf('%2.1f',max(max(pfknnmrl(:,:,u)))),'Color','w','FontWeight','bold','FontSize',10)


subplot(256);
bar(tbin,accg(:,u)),axis tight,

subplot(257);
pfw.plot(u),

subplot(258);
pfr.plot(u),

subplot(259);
pfh.plot(u)

subplot(2,5,10);
pfl.plot(u)

%print(gcf,'-dpsc2',[Trial.filebase '.knnpf_100ms_50msol_walk_rear' num2str(u) '.eps']);
    uicontrol(hfig,'String','>','units','normalized','Position',[.97, 0, .03,  1],'Callback',@(hfig,index,indArray,callback)figure_controls_gui(hfig,u,units,'forwardButton_Callback'));
    uicontrol(hfig,'String','<','units','normalized','Position',[  0, 0, .03,  1],'Callback',@(hfig,index,indArray,callback)figure_controls_gui(hfig,u,units,'backwardButton_Callback'));
    uicontrol(hfig,'String','X','units','normalized','Position',[.45, 0, .10,.07],'Callback',@(hfig,index,indArray,callback)figure_controls_gui(hfig,u,units,'exitButton_Callback'));
    uicontrol(hfig,'String','Print','units','normalized','Position',[.55, 0, .10,.07],'Callback',@(hfig,index,indArray,callback)figure_controls_gui(hfig,u,units,'printButton_Callback'));
    while strcmp(num2str(u),get(hfig,'Name'))
        pause(0.2);
    end
    u = str2double(get(hfig,'Name'));
end
 