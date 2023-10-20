
nconfigure_default_args();
EgoProCode2D_load_data();

EgoProCode2D_f2_data_egoHbaPhzPause();
EgoProCode2D_f2_data_egoHbaPhzLoc();

mask = double(sqrt(bsxfun(@plus,egoHbaPhzRmaps_loc.xbins.^2,egoHbaPhzRmaps_loc.ybins'.^2)') < 445);
mask(~mask) = nan;

rat = load_patch_model('rat')

sig = 1-(1-0.05)^(1/numel(unitsEgoCA1))

[hfig,fig,fax,sax] = set_figure_layout(figure(666003),'A4','portrait',[],1.5,1.5,0.1,0.1);

state = {'Locomotion','Pause'};

% loc

exampleUnit.trialIndex = 19;
exampleUnit.close.Xlims = [-200,400];
exampleUnit.close.Ylims = [-400,200];
exampleUnit.id = 50;
exampleUnit.maxRate = 25;
exampleUnit.index = find(unitsEgo{exampleUnit.trialIndex}==exampleUnit.id);
exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];
globalXOffset = -0.8;
globalYOffset = 0;




% SUBPLOTS -- EGO FIELD -- partitioned by theta-phase and head-body-angle
%%%<<<
for phzInd = 1:phzBin.count
    for hbaInd = 1:hbaBin.count
        
        %%%<<< PLOT egoField (phz x hba)
        % ADJUST subplot coordinates
        [yind, yOffSet, xind, xOffSet] = deal(phzBin.count+1-phzInd, 0, 1+hbaInd, 0);
        % CREATE subplot axes
        sax(end+1) = axes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                            fig.page.ypos(yind)+yOffSet+globalYOffset,...
                            fig.subplot.width,                        ...
                            fig.subplot.height],                      ...
                          'FontSize', 8,                                        ...
                          'LineWidth',1);
        hold(sax(end),'on');
        shading(sax(end),'flat');
        set(pcolor(egoHbaPhzRmaps_loc.xpos-diff(egoHbaPhzRmaps_loc.xpos(1:2))/2,...
                   egoHbaPhzRmaps_loc.ypos-diff(egoHbaPhzRmaps_loc.ypos(1:2))/2,...
                   fliplr(rot90(egoHbaPhzRmaps_loc.rmap{exampleUnit.trialIndex}(:,:,exampleUnit.index,phzInd,hbaInd)',-1)).*egoHbaPhzRmaps_loc.mask),'EdgeColor','none');
        axis(sax(end),'xy');
        colormap(sax(end),'jet');
        caxis(sax(end),[0,exampleUnit.maxRate]);
        
        % FORMAT subplot
        sax(end).XTickLabel =[];
        sax(end).YTickLabel =[];
        xlim(sax(end),[-250,250]);
        ylim(sax(end),[-250,250]);
        daspect(sax(end),[1,1,1]);
        box(sax(end),'on');

        Lines([],0,'w');
        Lines(0,[],'w');

        % ANNOTATE 
        subject = struct(rat);
        subject = update_subject_patch(subject, 'head',...
                                       [], false,...
                                       hbaBin.edges,...
                                       hbaBin.centers);
        subject = update_subject_patch(subject, 'body',...
                                       hbaBin.count+1-hbaInd,  true,...
                                       hbaBin.edges,...
                                       hbaBin.centers);
        patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
        patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
        patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
        line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
        line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
        axes(fax);
        rectangle('Position',sax(end).Position, ...
                  'EdgeColor',phzBin.color(phzInd,:), ...
                  'FaceColor','None',...
                  'LineWidth',1);
        if phzInd==3,
            line([sax(end).Position(1),sum(sax(end).Position([1,3]))],...
                  sum(sax(end).Position([2,4])).*[1,1]+0.15,...
                 'LineWidth',2,...
                 'Color',hbaBin.color(hbaInd,:));
            title(sax(end),{hbaBin.label{hbaInd},' '});
        end
        if phzInd==1 && hbaInd==3,
            line(sum(sax(end).Position([1,3])).*[1,1]+0.1,...
                 sax(end).Position(2).*[1,1]+[0,sax(end).Position(4)*0.4],...
                 'LineWidth',2,...
                 'Color',[0,0,0]);
            text(sum(sax(end).Position([1,3]))+0.4,...
                 sax(end).Position(2),...
                 '20 cm',...
                 'Rotation',90);
        end
    end
end
%%%>>>


% SUBPLOTS -- EGO FIELD -- partitioned by theta-phase and head-body-angle
%%%<<<


%% LOCOMOTION STATS -------------------------------------------------------
% SUBPLOT -- LATERAL POS -- left vs Right lateral coordinatats for egoHba
%%%<<<
% ADJUST subplot coordinates
for phzInd = 1:phzBin.count
    [yind, yOffSet, xind, xOffSet] = deal(phzBin.count+1-phzInd, 0, 6, 0);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                        fig.page.ypos(yind)+yOffSet+globalYOffset,...
                        fig.subplot.width,                        ...
                        fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    % PLOT subplot
    xlim(sax(end),[-20,20]);
    ylim(sax(end),[-20,20]);
    Lines([],0,'k');
    Lines(0,[],'k');
    plot(egoHbaPhzLoc.control.meanPos(unitsEgoCA1,phzInd,3,2), ...
         egoHbaPhzLoc.control.meanPos(unitsEgoCA1,phzInd,1,2), ...
         '.',                                           ...
         'MarkerFaceColor',phzBin.color(phzInd,:),      ...
         'MarkerEdgeColor',phzBin.color(phzInd,:));
    % FORMAT subplot
    grid(sax(end),'on');
    xlim(sax(end),[-10,10]);
    ylim(sax(end),[-10,10]);
    sax(end).XTick = [-5,0,5];
    sax(end).YTick = [-5,0,5];
    
    title(sax(end),{'Loc'});

    if phzInd == 2 
        ylabel(sax(end),'cm');
        sax(end).YLabel.Units = 'centimeters';
        sax(end).YLabel.Position = [-0.55,0.74,0]
    end
    if phzInd == 1
        xlabel(sax(end),'cm');
        sax(end).XLabel.Units = 'centimeters';
        sax(end).XLabel.Position = [0.75,-0.4,0];
    else
        sax(end).XTickLabel = {};
    end
    daspect(sax(end),[1,1,1]);
end
%%%>>>


% SUBPLOTS -- LAT POS DISTRIB -- partitioned by theta-phase and head-body-angle
%%%<<<
for phzInd = 1:phzBin.count
    % ADJUST subplot coordinates
    [yind, yOffSet, xind, xOffSet] = deal(phzBin.count+1-phzInd, 0, 7, 1.2);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                        fig.page.ypos(yind)+yOffSet+globalYOffset,...
                        fig.subplot.width,                        ...
                        fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    % PLOT subplot
    % FORMAT subplot
    grid(sax(end),'on');
    % PLOT subplot
    for hbaInd = 1:hbaBin.count
        [ehpcmpKDE,dxi] = ksdensity( egoHbaPhzLoc.control.meanPos(unitsEgoCA1,phzInd,hbaInd,lat) );
        med             = median(    egoHbaPhzLoc.control.meanPos(unitsEgoCA1,phzInd,hbaInd,lat) );
        plot(dxi,ehpcmpKDE,'-','color',hbaBin.color(hbaInd,:))
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hbaBin.color(hbaInd,:));
    end
    % FORMAT subplot
    xlim(sax(end),[-10,10])
    ylim(sax(end),[0,0.16])
    grid(sax(end),'on');
    sax(end).YTick = [0,0.05,0.10,0.15];
    sax(end).YTickLabel = {'0','','0.1',''};
    sax(end).XTick = [-5,0,5];

    if phzInd == 1
        xlabel(sax(end),'cm');
        sax(end).XLabel.Units = 'centimeters';
        sax(end).XLabel.Position = [0.75,-0.4,0];
    end

    if phzInd == 2
        ylabel(sax(end),'Prob')
        sax(end).YLabel.Units = 'centimeters';
        sax(end).YLabel.Position = [-0.55,0.74,0]
    end

end


% SUBPLOTS -- AP POS DISTRIB -- partitioned by theta-phase and head-body-angle
%%%<<<
for phzInd = 1:phzBin.count
    % ADJUST subplot coordinates
    [yind, yOffSet, xind, xOffSet] = deal(phzBin.count+1-phzInd, 0, 8, 1.2);

    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                        fig.page.ypos(yind)+yOffSet+globalYOffset,...
                        fig.subplot.width*1.75,                        ...
                        fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    % PLOT subplot
    % FORMAT subplot
    grid(sax(end),'on');
    % PLOT subplot
    for hbaInd = 1:hbaBin.count
        [ehpcmpKDE,dxi] = ksdensity( egoHbaPhzLoc.control.meanPos( unitsEgoCA1, phzInd, hbaInd, fwd) );
        med             = median(    egoHbaPhzLoc.control.meanPos( unitsEgoCA1, phzInd, hbaInd, fwd) );
        plot(dxi,ehpcmpKDE,'-','color',hbaBin.color(hbaInd,:))
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hbaBin.color(hbaInd,:));
    end
    % FORMAT subplot
    xlim(sax(end),[-15,20]);
    ylim(sax(end),[0,0.16]);
    grid(sax(end),'on');
    sax(end).YTick = [0,0.05,0.10,0.15];
    sax(end).YTickLabel = {};
    sax(end).XTick = [-10,-5,0,5,10,15,20];
    
    
    if phzInd==1,
        xlabel(sax(end),'cm');
        sax(end).XLabel.Units = 'centimeters';
        sax(end).XLabel.Position = [0.75,-0.4,0];
        sax(end).XTickLabel = {'-10','','0','','10','','20'};
    else
        sax(end).XTickLabel = {};
    end

end
%%%>>>


%%%>>>


%% PAUSE -------------------------------------------------------------------


exampleUnit.maxRate = 16;

for phzInd = 1:phzBin.count
    for hbaInd = 1:hbaBin.count
        
        %%%<<< PLOT egoField (phz x hba)
        % ADJUST subplot coordinates
        [yind, yOffSet, xind, xOffSet] = deal(phzBin.count-phzInd+5, 0, 1+hbaInd, 0);        
        % CREATE subplot axes
        sax(end+1) = axes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                            fig.page.ypos(yind)+yOffSet+globalYOffset,...
                            fig.subplot.width,                        ...
                            fig.subplot.height],                      ...
                          'FontSize', 8,                                        ...
                          'LineWidth',1);
        hold(sax(end),'on');
        shading(sax(end),'flat');
        set(pcolor(egoHbaPhzRmaps_pause.xpos-diff(egoHbaPhzRmaps_pause.xpos(1:2))/2,...
                   egoHbaPhzRmaps_pause.ypos-diff(egoHbaPhzRmaps_pause.ypos(1:2))/2,...
                   fliplr(rot90(egoHbaPhzRmaps_pause.rmap{exampleUnit.trialIndex}(:,:,exampleUnit.index,phzInd,hbaInd)',-1)).*egoHbaPhzRmaps_pause.mask),'EdgeColor','none');
        axis(sax(end),'xy');
        colormap(sax(end),'jet');
        caxis(sax(end),[0,exampleUnit.maxRate]);
        
        % FORMAT subplot
        sax(end).XTickLabel =[];
        sax(end).YTickLabel =[];
        xlim(sax(end),[-250,250]);
        ylim(sax(end),[-250,250]);
        daspect(sax(end),[1,1,1]);
        box(sax(end),'on');
        Lines([],0,'w');
        Lines(0,[],'w');
        % ANNOTATE 
        subject = struct(rat);
        subject = update_subject_patch(subject, 'head',...
                                       [], false,...
                                       hbaBin.edges,...
                                       hbaBin.centers);
        subject = update_subject_patch(subject, 'body',...
                                       hbaBin.count+1-hbaInd,  true,...
                                       hbaBin.edges,...
                                       hbaBin.centers);
        patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
        patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
        patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
        line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
        line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
        axes(fax);
        rectangle('Position',sax(end).Position, ...
                  'EdgeColor',phzBin.color(phzInd,:), ...
                  'FaceColor','None',...
                  'LineWidth',1);
        if phzInd==3,
            line([sax(end).Position(1),sum(sax(end).Position([1,3]))],...
                  sum(sax(end).Position([2,4])).*[1,1]+0.1,...
                 'LineWidth',2,...
                 'Color',hbaBin.color(hbaInd,:));
            title(sax(end),{hbaBin.label{hbaInd},' '});
        end
    end
end
%%%>>>    


%globalXOffset = 1;
% SUBPLOT -- LATERAL POS -- left vs Right lateral coordinatats for egoHba
%%%<<<
% ADJUST subplot coordinates
for phzInd = 1:phzBin.count
    [yind, yOffSet, xind, xOffSet] = deal(phzBin.count+1-phzInd+4, 0, 6, 0);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                        fig.page.ypos(yind)+yOffSet+globalYOffset,...
                        fig.subplot.width,                        ...
                        fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    % PLOT subplot
    xlim(sax(end),[-10,10]);
    ylim(sax(end),[-10,10]);
    Lines([],0,'k');
    Lines(0,[],'k');
    plot(egoHbaPhzPause.control.meanPos( unitsEgoCA1, phzInd, 3, lat) ,...
         egoHbaPhzPause.control.meanPos( unitsEgoCA1, phzInd, 1, lat) ,...
         '.',                                                   ...
         'MarkerFaceColor',phzBin.color(phzInd,:),              ...
         'MarkerEdgeColor',phzBin.color(phzInd,:));
    % FORMAT subplot
    grid(sax(end),'on');
    xlim(sax(end),[-10,10]);
    ylim(sax(end),[-10,10]);
    sax(end).XTick = [-5,0,5];
    sax(end).YTick = [-5,0,5];
    title(sax(end),{'Pause'});
    sax(end).YLabel.Units = 'centimeters';
    sax(end).YLabel.Position = [-0.55,0.74,0]
    daspect(sax(end),[1,1,1]);
    if phzInd == 2
        ylabel(sax(end),'cm');
        sax(end).YLabel.Units = 'centimeters';
        sax(end).YLabel.Position = [-0.55,0.74,0]
    end
    if phzInd == 1,
        xlabel(sax(end),'cm');
        sax(end).XLabel.Units = 'centimeters';
        sax(end).XLabel.Position = [0.75,-0.4,0];
    else
        sax(end).XTickLabel = {};
    end
end
%%%>>>



% SUBPLOTS -- LAT POS DISTRIB -- partitioned by theta-phase and head-body-angle
%%%<<<
for phzInd = 1:phzBin.count
    % ADJUST subplot coordinates
    [yind, yOffSet, xind, xOffSet] = deal(phzBin.count+1-phzInd+4, 0, 7, 1.2);

    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                        fig.page.ypos(yind)+yOffSet+globalYOffset,...
                        fig.subplot.width,                        ...
                        fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    % PLOT subplot
    % FORMAT subplot
    grid(sax(end),'on');
    % PLOT subplot
    for hbaInd = 1:hbaBin.count
        [ehpcmpKDE,dxi] = ksdensity( egoHbaPhzPause.control.meanPos( unitsEgoCA1, phzInd, hbaInd, lat) );
        med             = median(    egoHbaPhzPause.control.meanPos( unitsEgoCA1, phzInd, hbaInd, lat) );
        plot(dxi,ehpcmpKDE,'-','color',hbaBin.color(hbaInd,:))
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hbaBin.color(hbaInd,:));
    end
    % FORMAT subplot
    xlim(sax(end),[-10,10])
    ylim(sax(end),[0,0.16])
    grid(sax(end),'on');
    sax(end).YTick = [0,0.05,0.10,0.15];
    sax(end).YTickLabel = {'0','','0.1',''};
    sax(end).XTick = [-5,0,5];

    if phzInd == 2,
        ylabel(sax(end),'Prob');
        sax(end).YLabel.Units = 'centimeters';
        sax(end).YLabel.Position = [-0.55,0.74,0]
    end
    
    if phzInd==1,
        xlabel(sax(end),'cm');
        sax(end).XLabel.Units = 'centimeters';
        sax(end).XLabel.Position = [0.75,-0.4,0];
    else
        sax(end).XTickLabel = {};
    end

end
%%%>>>



% SUBPLOTS -- AP POS DISTRIB -- partitioned by theta-phase and head-body-angle
%%%<<<
for phzInd = 1:phzBin.count
    % ADJUST subplot coordinates
    [yind, yOffSet, xind, xOffSet] = deal(phzBin.count+1-phzInd+4, 0, 8, 1.2);

    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                        fig.page.ypos(yind)+yOffSet+globalYOffset,...
                        fig.subplot.width*1.75,                        ...
                        fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    % PLOT subplot
    % FORMAT subplot
    grid(sax(end),'on');
    % PLOT subplot
    for hbaInd = 1:hbaBin.count
        [ehpcmpKDE,dxi] = ksdensity( egoHbaPhzPause.control.meanPos( unitsEgoCA1, phzInd, hbaInd, fwd) );
        med             = median(    egoHbaPhzPause.control.meanPos( unitsEgoCA1, phzInd, hbaInd, fwd) );
        plot(dxi,ehpcmpKDE,'-','color',hbaBin.color(hbaInd,:))
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hbaBin.color(hbaInd,:));
    end
    % FORMAT subplot
    xlim(sax(end),[-15,20])
    ylim(sax(end),[0,0.16])
    grid(sax(end),'on');
    sax(end).YTick = [0,0.05,0.10,0.15];
    sax(end).YTickLabel = {};
    sax(end).XTick = [-10,-5,0,5,10,15,20];
    
    if phzInd==1,
        xlabel(sax(end),'cm');
        sax(end).XLabel.Units = 'centimeters';
        sax(end).XLabel.Position = [0.75,-0.4,0];
        sax(end).XTickLabel = {'-10','','0','','10','','20'};
    else
        sax(end).XTickLabel = {};
    end

end
%%%>>>



% $$$ 
% $$$ % LEGEND hba colors
% $$$ axes(fax);
% $$$ for hbaInd = 1:hbaBin.count
% $$$     lh(hbaInd) = line(fig.page.xpos(5)+[0,0.5]+1,...
% $$$                       fig.page.ypos(3).*[1,1]-hbaInd*0.2+0.75,...
% $$$                       'Color',hbaBin.color(hbaInd,:),...
% $$$                       'LineWidth',2);
% $$$     text(fig.page.xpos(5)+1.75,...
% $$$                       fig.page.ypos(3)-hbaInd*0.2+0.75,...
% $$$                       hbaBin.label{hbaInd},...
% $$$                       'VerticalAlignment','middle');
% $$$ end




% THETA Phase Vertical
%%%<<<
[yind, yOffSet, xind, xOffSet] = deal(3, 0, 1, 0.6);
subplotWidth = 0.8;
subplotHeight = fig.subplot.height * 3 + fig.subplot.verticalPadding * 2
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              subplotWidth,                             ...
                              subplotHeight],                           ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
% SUBPLOT <- theta cycle (color=partions)
hold(sax(end),'on');
plot(sax(end),-cos(linspace(0,2*pi)),linspace(0,1),'k','LineWidth',2);
plims = [0.5,              2.26106176905986];
plot(sax(end),...
     -cos(linspace(plims(1),plims(2))),...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',phzBin.color(1,:), 'LineWidth',2);
plims = [2.26106176905986, 4.02212353811972];
plot(sax(end),...
     -cos(linspace(plims(1),plims(2))),...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',phzBin.color(2,:), 'LineWidth',2);
plims = [4.02212353811972, 5.78318530717959];
plot(sax(end),...
     -cos(linspace(plims(1),plims(2))),...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',phzBin.color(3,:),'LineWidth',2);

% FORMAT subplot
ylim([0,1]);
sax(end).Visible = 'off';
text(0.5,0.05,'0','Rotation',90,'FontSize',8);
text(-0.5,0.5,'\pi','Rotation',90,'FontSize',8);
text(0.5,0.95,'2\pi','Rotation',90,'FontSize',8);

ylabel(sax(end),'Locomotion')
sax(end).XAxis.Visible = 'off';
sax(end).YAxis.Color = 'k';
sax(end).YTick =[];
sax(end).Color = 'none';


%%%>>>



% THETA Phase Vertical
%%%<<<
[yind, yOffSet, xind, xOffSet] = deal(7, 0, 1, 0.6);
subplotWidth = 0.8;
subplotHeight = fig.subplot.height * 3 + fig.subplot.verticalPadding * 2
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              subplotWidth,                             ...
                              subplotHeight],                           ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
% SUBPLOT <- theta cycle (color=partions)
hold(sax(end),'on');
plot(sax(end),-cos(linspace(0,2*pi)),linspace(0,1),'k','LineWidth',2);
plims = [0.5,              2.26106176905986];
plot(sax(end),...
     -cos(linspace(plims(1),plims(2))),...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',phzBin.color(1,:), 'LineWidth',2);
plims = [2.26106176905986, 4.02212353811972];
plot(sax(end),...
     -cos(linspace(plims(1),plims(2))),...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',phzBin.color(2,:), 'LineWidth',2);
plims = [4.02212353811972, 5.78318530717959];
plot(sax(end),...
     -cos(linspace(plims(1),plims(2))),...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',phzBin.color(3,:),'LineWidth',2);
% FORMAT subplot
ylim([0,1]);
text(0.5,0.05,'0','Rotation',90,'FontSize',8);
text(-0.5,0.5,'\pi','Rotation',90,'FontSize',8);
text(0.5,0.95,'2\pi','Rotation',90,'FontSize',8);
ylabel(sax(end),'Pause');
sax(end).XAxis.Visible = 'off';
sax(end).YAxis.Color = 'k';
sax(end).YTick =[];
sax(end).Color = 'none';
%%%>>>





% $$$ 
% $$$ 
% $$$ for phzInd = 1:phzBin.count
% $$$     for velInd = 1:velBin.count-1
% $$$     [yind, yOffSet, xind, xOffSet] = deal(phzBin.count+1-phzInd+8, 0, velInd, 0);
% $$$     % CREATE subplot axes
% $$$     sax(end+1) = axes('Units','centimeters',                                ...
% $$$                       'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
% $$$                         fig.page.ypos(yind)+yOffSet+globalYOffset,...
% $$$                         fig.subplot.width,                        ...
% $$$                         fig.subplot.height],                      ...
% $$$                       'FontSize', 8,                                        ...
% $$$                       'LineWidth',1);
% $$$     hold(sax(end),'on');
% $$$ plot(pfe{exampleUnit.trialIndex},exampleUnit.id,1,'',[0,exampleUnit.maxRate],'colorMap',@jet,'mazeMaskFlag',false,'flipAxesFlag',true);        
% $$$     plot(pfv{exampleUnit.trialIndex}{phzInd,velInd+1},...
% $$$          exampleUnit.id,...
% $$$          1,...
% $$$          '',...
% $$$          [0,exampleUnit.maxRate],...
% $$$          false,...
% $$$          [],...
% $$$          true,...
% $$$          'colorMap',@jet);
% $$$     sax(end).XTickLabel =[];
% $$$     sax(end).YTickLabel =[];
% $$$     xlim(sax(end),[-250,250])
% $$$     ylim(sax(end),[-250,250])
% $$$     Lines([],0,'w');
% $$$     Lines(0,[],'w');
% $$$     end
% $$$ end




% $$$ 
% $$$ figure,
% $$$ hold('on');
% $$$ for phzInd = 1:phzBin.count
% $$$ plot(1:4,sq(mean(egoHvfPhz.control.meanPos(unitsEgoHvfCA3,phzInd,:,1))),'-+','Color',phzBin.color(phzInd,:));
% $$$ end
% $$$ 
% $$$ figure,
% $$$ hold('on');
% $$$ for phzInd = 1:phzBin.count
% $$$     plot(1:4,sq(mean(egoHvfPhz.control.meanPos(unitsEgoHvfCA1,phzInd,:,1))),'-+','Color',phzBin.color(phzInd,:));
% $$$ end


for phzInd = 1:phzBin.count
    for velInd = 1:velBin.count-1
    [yind, yOffSet, xind, xOffSet] = deal(phzBin.count+1-phzInd+8, 0, velInd, 0);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                        fig.page.ypos(yind)+yOffSet+globalYOffset,...
                        fig.subplot.width,                        ...
                        fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
        shading(sax(end),'flat');
        set(pcolor(xpos-diff(xpos(1:2))/2,...
                   ypos-diff(ypos(1:2))/2,...
                   fliplr(rot90(rmap{exampleUnit.trialIndex}(:,:,exampleUnit.index,phzInd,velInd)',-1)).*mask),'EdgeColor','none');
        axis(sax(end),'xy');
        colormap(sax(end),'jet');
        caxis(sax(end),[0,exampleUnit.maxRate]);
    sax(end).XTickLabel =[];
    sax(end).YTickLabel =[];
    xlim(sax(end),[-250,250])
    ylim(sax(end),[-250,250])
    Lines([],0,'w');
    Lines(0,[],'w');
    end
end




exampleUnit.trialIndex = 20;
exampleUnit.close.Xlims = [-200,400];
exampleUnit.close.Ylims = [-400,200];
exampleUnit.id = 21;
exampleUnit.maxRate = 21;
exampleUnit.index = find(unitsEgo{exampleUnit.trialIndex}==exampleUnit.id);
exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];
uids = unitsEgoCA1;
globalXOffset = 0;
globalYOffset = 0;

figure
sax = gobjects([0,1]);
for velInd = 1:velBin.count-1
    sax(end+1) = subplot(1,3,velInd);
    hold(sax(end),'on');
        shading(sax(end),'flat');
        set(pcolor(xpos-diff(xpos(1:2))/2,...
                   ypos-diff(ypos(1:2))/2,...
                   fliplr(rot90(rmap{exampleUnit.trialIndex}(:,:,1,velInd)',-1)).*mask),'EdgeColor','none');
        axis(sax(end),'xy');
        colormap(sax(end),'jet');
        %caxis(sax(end),[0,exampleUnit.maxRate]);
    sax(end).XTickLabel =[];
    sax(end).YTickLabel =[];
    xlim(sax(end),[-350,350])
    ylim(sax(end),[-350,350])
    Lines([],0,'w');
    Lines(0,[],'w');
end





figure
        set(pcolor(xpos-diff(xpos(1:2))/2,...
                   ypos-diff(ypos(1:2))/2,...
                   fliplr(rot90(rmap{exampleUnit.trialIndex}(:,:,exampleUnit.index,velInd)',-1)).*mask),'EdgeColor','none');



figure,plot(atan2(egoHbaPhzLoc.control.meanPos(unitsEgoCA1,3,3,2),egoHbaPhzLoc.control.meanPos(unitsEgoCA1,3,3,1)),atan2(egoHbaPhzPause.control.meanPos(unitsEgoCA1,3,3,2),egoHbaPhzPause.control.meanPos(unitsEgoCA1,3,3,1)),'.','MarkerSize',10)


ind = false([size(egoHbaPhzLoc.control.meanPos,1),1]);
ind(unitsEgoCA1) = true;
ind(egoHbaPhzLoc.boot.ca1.zscore(:,3,hbaInd,2) < egoHbaPhzLoc.boot.ca1.sig) = false;

figure()
for hbaInd = 1:hbaBin.count
ind = false([size(egoHbaPhzLoc.control.meanPos,1),1]);
ind(unitsEgoCA1) = true;
if hbaInd == 1
    ind(egoHbaPhzLoc.boot.ca1.zscore(:,3,hbaInd,2) > -egoHbaPhzLoc.boot.ca1.sig ...
        &egoHbaPhzPause.boot.ca1.zscore(:,3,hbaInd,2) > -egoHbaPhzPause.boot.ca1.sig) = false;
elseif hbaInd == 2
    %ind(abs(egoHbaPhzLoc.boot.ca1.zscore(:,3,hbaInd,2)) < egoHbaPhzLoc.boot.ca1.sig) = false;
else
    ind(egoHbaPhzLoc.boot.ca1.zscore(:,3,hbaInd,2) < egoHbaPhzLoc.boot.ca1.sig ...
        &egoHbaPhzPause.boot.ca1.zscore(:,3,hbaInd,2) < egoHbaPhzPause.boot.ca1.sig) = false;
end
    subplot(hbaBin.count,1,hbaInd)
    hold('on');
% $$$ plot(atan2(  egoHbaPhzLoc.control.meanPos(unitsEgoCA1,3,hbaInd,2),  egoHbaPhzLoc.control.meanPos(unitsEgoCA1,3,hbaInd,1)),...
% $$$      atan2(egoHbaPhzPause.control.meanPos(unitsEgoCA1,3,hbaInd,2),egoHbaPhzPause.control.meanPos(unitsEgoCA1,3,hbaInd,1)),...
% $$$      '.','Color',hbaBin.color(hbaInd,:),'MarkerSize',30)
plot(atan2(  egoHbaPhzLoc.control.meanPos(ind,3,hbaInd,2),  egoHbaPhzLoc.control.meanPos(ind,3,hbaInd,1)),...
     atan2(egoHbaPhzPause.control.meanPos(ind,3,hbaInd,2),egoHbaPhzPause.control.meanPos(ind,3,hbaInd,1)),...
     '.','Color',hbaBin.color(hbaInd,:),'MarkerSize',30)
grid('on');
xlim([-3,3])
ylim([-3,3])
end

figure,plot(egoHvfPhz.control.meanPos(unitsEgoHvfCA1,3,4,2),  egoHvfPhz.control.meanPos(unitsEgoHvfCA1,3,2,1),'.')

egoHbaPhzLoc.boot.ca1.zscore(:,3,hbaInd,2)> ...
    egoHbaPhzLoc.boot.ca1.sig





figure()
for hbaInd = 1:hbaBin.count
ind = false([size(egoHbaPhzLoc.control.meanPos,1),1]);
ind(unitsEgoCA1) = true;
if hbaInd == 1
    ind(egoHbaPhzLoc.boot.ca1.zscore(:,3,hbaInd,2) > -egoHbaPhzLoc.boot.ca1.sig ...
        &egoHbaPhzPause.boot.ca1.zscore(:,3,hbaInd,2) > -egoHbaPhzPause.boot.ca1.sig) = false;
elseif hbaInd == 2
    %ind(abs(egoHbaPhzLoc.boot.ca1.zscore(:,3,hbaInd,2)) < egoHbaPhzLoc.boot.ca1.sig) = false;
else
    ind(egoHbaPhzLoc.boot.ca1.zscore(:,3,hbaInd,2) < egoHbaPhzLoc.boot.ca1.sig ...
        &egoHbaPhzPause.boot.ca1.zscore(:,3,hbaInd,2) < egoHbaPhzPause.boot.ca1.sig) = false;
end
    subplot(hbaBin.count,1,hbaInd)
    hold('on');
% $$$ plot(atan2(  egoHbaPhzLoc.control.meanPos(unitsEgoCA1,3,hbaInd,2),  egoHbaPhzLoc.control.meanPos(unitsEgoCA1,3,hbaInd,1)),...
% $$$      atan2(egoHbaPhzPause.control.meanPos(unitsEgoCA1,3,hbaInd,2),egoHbaPhzPause.control.meanPos(unitsEgoCA1,3,hbaInd,1)),...
% $$$      '.','Color',hbaBin.color(hbaInd,:),'MarkerSize',30)
plot(egoHbaPhzLoc.control.meanPos(ind,3,hbaInd,2),...
     egoHbaPhzPause.control.meanPos(ind,3,hbaInd,2),...
     '.','Color',hbaBin.color(hbaInd,:),'MarkerSize',30)
plot(mean(egoHbaPhzLoc.control.meanPos(ind,3,hbaInd,2)),...
     mean(egoHbaPhzPause.control.meanPos(ind,3,hbaInd,2)),...
     '.','Color','m','MarkerSize',50)
grid('on');
xlim([-15,15])
ylim([-15,15])
Lines([],0,'k');
Lines(0,[],'k');
end






figure()
for hbaInd = 1:hbaBin.count
ind = false([size(egoHbaPhzLoc.control.meanPos,1),1]);
ind(unitsEgoCA1) = true;
if hbaInd == 1
    ind(egoHbaPhzLoc.boot.ca1.zscore(:,3,hbaInd,2) > -egoHbaPhzLoc.boot.ca1.sig ...
        &egoHbaPhzPause.boot.ca1.zscore(:,3,hbaInd,2) > -egoHbaPhzPause.boot.ca1.sig) = false;
elseif hbaInd == 2
    %ind(abs(egoHbaPhzLoc.boot.ca1.zscore(:,3,hbaInd,2)) < egoHbaPhzLoc.boot.ca1.sig) = false;
else
    ind(egoHbaPhzLoc.boot.ca1.zscore(:,3,hbaInd,2) < egoHbaPhzLoc.boot.ca1.sig ...
        &egoHbaPhzPause.boot.ca1.zscore(:,3,hbaInd,2) < egoHbaPhzPause.boot.ca1.sig) = false;
end
    subplot(hbaBin.count,1,hbaInd)
    hold('on');
% $$$ plot(atan2(  egoHbaPhzLoc.control.meanPos(unitsEgoCA1,3,hbaInd,2),  egoHbaPhzLoc.control.meanPos(unitsEgoCA1,3,hbaInd,1)),...
% $$$      atan2(egoHbaPhzPause.control.meanPos(unitsEgoCA1,3,hbaInd,2),egoHbaPhzPause.control.meanPos(unitsEgoCA1,3,hbaInd,1)),...
% $$$      '.','Color',hbaBin.color(hbaInd,:),'MarkerSize',30)
plot(egoHbaPhzLoc.control.meanPos(ind,3,hbaInd,1),...
     egoHbaPhzPause.control.meanPos(ind,3,hbaInd,1),...
     '.','Color',hbaBin.color(hbaInd,:),'MarkerSize',30)
plot(mean(egoHbaPhzLoc.control.meanPos(ind,3,hbaInd,1)),...
     mean(egoHbaPhzPause.control.meanPos(ind,3,hbaInd,1)),...
     '.','Color','m','MarkerSize',50)
grid('on');
xlim([-15,15])
ylim([-15,15])
Lines([],0,'k');
Lines(0,[],'k');
end
