



% shuffling
% premutation test

%             phz : theta phase
%             hba : head-body-angle ( 0=straight, pi/2=left/right, pi=broken rat)

                    {HBA-L}   {HBA-C}   {HBA-R}
        allofield  egofield  egofield  egofield  plot(HREF-LAT(LxR))  ksden(HREF-LAT(LCR))
{PHZ-A} egofield      X         X         X      plot(HREF-LAT(LxR))  ksden(HREF-LAT(LCR))
{PHZ-T} egofield      X         X         X      plot(HREF-LAT(LxR))  ksden(HREF-LAT(LCR))
(PHZ-D) egofield      X         X         X      plot(HREF-LAT(LxR))  ksden(HREF-LAT(LCR))





global MTA_PROJECT_PATH
partsPath = fullfile(fullfile(MTA_PROJECT_PATH,'analysis','EgoProCode2D','EgoProCode2D_figure_parts'));
overwrite = false;
rat = load_patch_model('rat');

mask = double(sqrt(bsxfun(@plus,egoHbaRmaps.xbins.^2,egoHbaRmaps.ybins'.^2)') < 445);
mask(~mask) = nan;



[hfig,fig,fax,sax] = set_figure_layout(figure(666003),'A4','landscape',[],1.8,1.8,0.2,0.4);

regions = {'ca1','ca3'};

for region = 1%:2
    if region == 1
        % CA1 
        exampleUnit.trialIndex = 20;
        exampleUnit.close.Xlims = [-200,400];
        exampleUnit.close.Ylims = [-400,200];
        exampleUnit.id = 25;
        exampleUnit.maxRate = 18;
        exampleUnit.index = find(unitsEgo{exampleUnit.trialIndex}==exampleUnit.id);
        exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];
        uids = unitsEgoCA1;
        gXoffset = 0;
        gYoffset = 0;
    else
        % CA3
        exampleUnit.trialIndex = 6;
        exampleUnit.close.Xlims = [-200,400];
        exampleUnit.close.Ylims = [-400,200];
        exampleUnit.id = 10;
        exampleUnit.maxRate = 18;
        %exampleUnit.index = find(units{exampleUnit.trialIndex}==exampleUnit.id);
        exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];
        uids = unitsEgoCA3;
        gXoffset = 0;
        gYoffset = -2.4*5;
    end

% SUBPLOT -- PLACE FIELD -- Theta
%%%<<<    
    [yind, yOffSet, xind, xOffSet] = deal(1, 0, 1, 0);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                                  fig.page.ypos(yind)+yOffSet+gYoffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    title(sax(end),{'Place Field'});

    plot(pft{exampleUnit.trialIndex},exampleUnit.id,1,'',[0,exampleUnit.maxRate],'colorMap',@jet);
    sax(end).XTickLabel =[];
    sax(end).YTickLabel =[];
%%%>>>    


% SUBPLOT -- HBA Partition Diagram -- Theta    
%%%<<<
    for d = 1:size(pfs{exampleUnit.trialIndex},2)
        [yind, yOffSet, xind, xOffSet] = deal(1, 0, d+1, 0);
        % CREATE subplot axes
        sax(end+1) = axes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                                      fig.page.ypos(yind)+yOffSet+gYoffset,...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height],                      ...u
                          'FontSize', 8,                                        ...
                          'LineWidth',1);
        hold(sax(end),'on');

        % ADD rat
        shading(sax(end),'flat');
        set(pcolor(egoHba.xpos-diff(egoHba.xpos(1:2))/2,egoHba.ypos-diff(egoHba.ypos(1:2))/2,...
                   fliplr(rot90(egoHba.rmap{exampleUnit.trialIndex}(:,:,exampleUnit.index,d)',-1)).*mask),'EdgeColor','none');
        axis(sax(end),'xy');
        colormap(sax(end),'jet');
        caxis(sax(end),[0,exampleUnit.maxRate]);
        ylim(sax(end),[egoHba.ypos([1,end])+[-1,1].*diff(egoHba.ypos(1:2))/2])
        xlim(sax(end),[egoHba.xpos([1,end])+[-1,1].*diff(egoHba.xpos(1:2))/2])
        Lines([],0,'k');
        Lines(0,[],'k');
        xlim(sax(end),[-250,250]);
        ylim(sax(end),[-250,250]);

        sax(end).XTick =[];
        sax(end).YTick =[];
        sax(end).XTickLabel =[];
        sax(end).YTickLabel =[];

        if d == 2
            title(sax(end),{'EgoFields','Partioned by Head-Body Angle','',hbaBin.label{d}})
        else
            title(sax(end),hbaBin.label{d})
        end
        
        
        daspect(sax(end),[1,1,1]);
        box(sax(end),'on');
        
    end
%%%>>>    

% SUBPLOTS -- EGO FIELD -- partitioned by theta phase
%%%<<<
    for phzInd = 1:phzBin.count
        % ADJUST subplot coordinates
        [yind, yOffSet, xind, xOffSet] = deal(phzInd+1, 0, 1, 0);
        % CREATE subplot axes
        sax(end+1) = axes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                                      fig.page.ypos(yind)+yOffSet+gYoffset,...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height],                      ...
                          'FontSize', 8,                                        ...
                          'LineWidth',1);
        hold(sax(end),'on');

        plot(pfet{exampleUnit.trialIndex}{4-phzInd},                            ...
             exampleUnit.id,1,'',[0,exampleUnit.maxRate],                       ...
             'colorMap',@jet,                                                   ...
             'mazeMaskFlag',false,                                              ...
             'flipAxesFlag',true);

        % FORMAT subplot
        sax(end).XTickLabel =[];
        sax(end).YTickLabel =[];
        xlim(sax(end),[-250,250]);
        ylim(sax(end),[-250,250]);
        daspect(sax(end),[1,1,1]);
        box(sax(end),'on');

        Lines([],0,'w');
        Lines(0,[],'w');

        if phzInd == 2
            ylabel(sax(end),{'EgoFields Partioned by Head-Body Angle',phzBin.label{phzInd}})
        else
            ylabel(sax(end),phzBin.label{phzInd})
        end
        

        %text(ets(exampleRange(1))+0.25,380,'Y');
        axes(fax);
     
    end
%%%>>>



% SUBPLOTS -- EGO FIELD -- partitioned by theta-phase and head-body-angle
%%%<<<
    for phzInd = 1:phzBin.count
        for hbaInd = 1:hbaBin.count
            
            %%%<<< PLOT egoField (phz x hba)
            % ADJUST subplot coordinates
            [yind, yOffSet, xind, xOffSet] = deal(phzInd+1, 0, hbaInd+1, 0);
            % CREATE subplot axes
            sax(end+1) = axes('Units','centimeters',                                ...
                              'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                                          fig.page.ypos(yind)+yOffSet+gYoffset,...
                                          fig.subplot.width,                        ...
                                          fig.subplot.height],                      ...
                              'FontSize', 8,                                        ...
                              'LineWidth',1);
            hold(sax(end),'on');

            plot(pfs{exampleUnit.trialIndex}{4-phzInd,hbaInd},                      ...
                 exampleUnit.id,1,'',[0,exampleUnit.maxRate],                       ...
                 'colorMap',@jet,                                                   ...
                 'mazeMaskFlag',false,                                              ...
                 'flipAxesFlag',true);

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
        end
    end
%%%>>>

    
    
% SUBPLOT -- LATERAL POS -- left vs Right lateral coordinatats for egoHba
%%%<<<
    % ADJUST subplot coordinates
    [yind, yOffSet, xind, xOffSet] = deal(1, 0, 6, 0);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                                  fig.page.ypos(yind)+yOffSet+gYoffset,...
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
    plot(egoHba.control.meanPos(uids,3,2) ,...
         egoHba.control.meanPos(uids,1,2) ,...
         '.',...
         'MarkerFaceColor','b',...
         'MarkerEdgeColor','b');
    % FORMAT subplot
    grid(sax(end),'on');
    xlim(sax(end),[-10,10]);
    ylim(sax(end),[-10,10]);
    sax(end).XTick = [-10,-5,0,5,10];
    sax(end).YTick = [-5,0,5];
    sax(end).XTickLabel = {};
    title(sax(end),{'Lateral'});
    daspect(sax(end),[1,1,1]);
%%%>>>



% $$$ % SUBPLOTS -- LAT POS ECDF -- partitioned by theta-phase and head-body-angle
% $$$     %%%<<< PLOT ecdf ofleft vs Right lateral coordinatats for egoHba
% $$$     % ADJUST subplot coordinates
% $$$     [yind, yOffSet, xind, xOffSet] = deal(1, 0, 8, -1.2);
% $$$     % CREATE subplot axes
% $$$     sax(end+1) = axes('Units','centimeters',                                ...
% $$$                       'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
% $$$                                   fig.page.ypos(yind)+yOffSet+gYoffset,...
% $$$                         fig.subplot.width,                        ...
% $$$                         fig.subplot.height],                      ...
% $$$                       'FontSize', 8,                                        ...
% $$$                       'LineWidth',1);
% $$$     hold(sax(end),'on');
% $$$     % PLOT subplot
% $$$     xlim(sax(end),[-10,10]);
% $$$     ylim(sax(end),[0,1]);
% $$$     Lines([],0.5,'k');
% $$$     Lines(0,[],'k');
% $$$     [F,X] = ecdf((egoHba.control.meanPos(uids,3,2)-egoHba.control.meanPos(uids,2,2)));
% $$$     plot(X,F,'Color','b');
% $$$     [F,X] = ecdf((egoHba.control.meanPos(uids,2,2)-egoHba.control.meanPos(uids,1,2)));
% $$$     plot(X,F,'Color','r');
% $$$     [F,X] = ecdf((egoHba.control.meanPos(uids,3,2)-egoHba.control.meanPos(uids,1,2)));
% $$$     plot(X,F,'Color',[219/255,172/255,52/255]);    
% $$$     % FORMAT subplot
% $$$     grid(sax(end),'on');
% $$$     xlim(sax(end),[-10,10]);
% $$$     ylim(sax(end),[0,1]);
% $$$     sax(end).XTick = [-10,-5,0,5,10];
% $$$     sax(end).YTick = [0,0.5,1];
% $$$     sax(end).XTickLabel = {};
% $$$     title(sax(end),{'cdf','Lateral - Center'});


% SUBPLOT -- LAT POS DISTRIB -- partitioned by head-body-angle
%%%<<<
    % ADJUST subplot coordinates
    [yind, yOffSet, xind, xOffSet] = deal(1, 0, 8, -1.2);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                                  fig.page.ypos(yind)+yOffSet+gYoffset,...
                        fig.subplot.width,                        ...
                        fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    % PLOT subplot
    for hbaInd = 1:hbaBin.count
        [ehpcmpKDE,dxi] = ksdensity(egoHba.control.meanPos(unitsEgoCA1,hbaInd,2));
        plot(dxi,ehpcmpKDE,'-','color',hbaBin.color(hbaInd,:))
        med = median(egoHba.control.meanPos(unitsEgoCA1,hbaInd,2));
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hbaBin.color(hbaInd,:));
    end
    % FORMAT subplot
    xlim(sax(end),[-10,10])
    ylim(sax(end),[0,0.16])
    grid(sax(end),'on');
    sax(end).YTick = [0,0.05,0.10,0.15];
    sax(end).YTickLabel = {'0','','0.1',''};        
    sax(end).XTick = [-10,-5,0,5,10];
    sax(end).XTickLabel = {};
    title(sax(end),{'Lateral', 'Position'});
%%%>>>


% SUBPLOT -- AP POS DISTRIB -- partitioned by head-body-angle   
%%%<<<
    [yind, yOffSet, xind, xOffSet] = deal(1, 0, 9, -1.2);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                                  fig.page.ypos(yind)+yOffSet+gYoffset,...
                        fig.subplot.width*1.75,                        ...
                        fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    % PLOT subplot
    for hbaInd = 1:hbaBin.count
        [ehpcmpKDE,dxi] = ksdensity(egoHba.control.meanPos(unitsEgoCA1,hbaInd,1));
        plot(dxi,ehpcmpKDE,'-','color',hbaBin.color(hbaInd,:));
        med = median(egoHba.control.meanPos(unitsEgoCA1,hbaInd,1));
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hbaBin.color(hbaInd,:));
    end
    % FORMAT subplot
    xlim(sax(end),[-15,25])
    ylim(sax(end),[0,0.16])
    grid(sax(end),'on');
    sax(end).YTick = [0,0.05,0.10,0.15];
    sax(end).YTickLabel = {};        
    sax(end).XTick = [-10,-5,0,5,10,15,20];
    sax(end).XTickLabel = {};
    title(sax(end),{'Anteroposterior', 'Position'});
%%%>>>    


% SUBPLOT -- LAT BOOT ZSCR ECDF -- partitioned by head-body-angle
%%%<<<
    % ADJUST subplot coordinates
    [yind, yOffSet, xind, xOffSet] = deal(1, 0, 11, -1.2);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                        fig.page.ypos(yind)+yOffSet+gYoffset,...
                        fig.subplot.width,                        ...
                        fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    % PLOT subplot
    for hbaInd = 1:hbaBin.count
        [F,X] = ecdf(egoHba.boot.(regions{region}).zscore(:,hbaInd,lat));
        plot(X,F,'Color',hbaBin.color(hbaInd,:));
    end
    
    % FORMAT subplot
    grid(sax(end),'on');
    xlim(sax(end),[-15,15]);
    ylim(sax(end),[0,1]);
    sax(end).XTick = [-10,-5,0,5,10];
    sax(end).YTick = [0,0.5,1];
    sax(end).YTick = [0,0.25,0.5,0.75,1];
    sax(end).YTickLabel = {};
    sax(end).XTickLabel = {};
    title(sax(end),{'Lateral','z-score'});
                    
    
    line(egoHba.perm.(regions{region}).sig.*[1,1],[0,1],'LineStyle','--','Color',[0.25,0.25,0.25]);
    line(-egoHba.perm.(regions{region}).sig.*[1,1],[0,1],'LineStyle','--','Color',[0.25,0.25,0.25]);
%%%>>>


    
% SUBPLOT -- AP BOOT ZSCR ECDF -- partitioned by head-body-angle
%%%<<<
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(1, 0, 12, -1.2);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                        fig.page.ypos(yind)+yOffSet+gYoffset,...
                        fig.subplot.width,                        ...
                        fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    % PLOT subplot
    for hbaInd = 1:hbaBin.count
        [F,X] = ecdf(egoHba.boot.(regions{region}).zscore(:,hbaInd,fwd));
        plot(X,F,'Color',hbaBin.color(hbaInd,:));
    end
    
    % FORMAT subplot
    grid(sax(end),'on');
    xlim(sax(end),[-15,15]);
    ylim(sax(end),[0,1]);
    sax(end).XTick = [-10,-5,0,5,10];
    sax(end).YTick = [0,0.5,1];
    sax(end).YTick = [0,0.25,0.5,0.75,1];
    sax(end).YTickLabel = {};
    sax(end).XTickLabel = {};
    title(sax(end),{'AP','z-score'});
    line(egoHba.perm.(regions{region}).sig.*[1,1],[0,1],'LineStyle','--','Color',[0.25,0.25,0.25]);
    line(-egoHba.perm.(regions{region}).sig.*[1,1],[0,1],'LineStyle','--','Color',[0.25,0.25,0.25]);
%%%>>>


% SUBPLOTS -- LATERAL POS -- left vs Right lateral coordinatats for egoHba
%%%<<<    
    for phzInd = 1:phzBin.count
        % ADJUST subplot coordinates
        [yind, yOffSet, xind, xOffSet] = deal(4-phzInd+1, 0, 6, 0);
        % CREATE subplot axes
        sax(end+1) = axes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                                      fig.page.ypos(yind)+yOffSet+gYoffset,...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height],                      ...
                          'FontSize', 8,                                        ...
                          'LineWidth',1);
        hold(sax(end),'on');
        %sax(end).tag = 
        % PLOT subplot
        xlim(sax(end),[-10,10]);
        ylim(sax(end),[-10,10]);
        Lines([],0,'k');
        Lines(0,[],'k');
        plot(egoHbaPhz.control.meanPos(uids,phzInd,3,2) ,...
             egoHbaPhz.control.meanPos(uids,phzInd,1,2) ,...
             '.',...
             'MarkerFaceColor',phzBin.color(phzInd,:),...
             'MarkerEdgeColor',phzBin.color(phzInd,:));
        % FORMAT subplot
        grid(sax(end),'on');
        xlim(sax(end),[-10,10]);
        ylim(sax(end),[-10,10]);
        sax(end).XTick = [-10,-5,0,5,10];
        sax(end).YTick = [-5,0,5];
        if phzInd ~= 1,    sax(end).XTickLabel = {};                           end % not bottom of column
        if phzInd == 1,    xlabel(sax(end),'R Lat (cm)');                              end % bottom of column
        % $$$         if p == 3,    title(sax(end),{'Lateral Shift','Right VS Left'}); end % top of column
        if phzInd == 1,    ylabel(sax(end),'L Lat (cm)');                              end % middle of column
        daspect(sax(end),[1,1,1]);
    end
%%%>>>

% SUBPLOTS -- LAT POS DISTRIB -- partitioned by theta-phase and head-body-angle
%%%<<<
    for phzInd = 1:phzBin.count
        % ADJUST subplot coordinates
        [yind, yOffSet, xind, xOffSet] = deal(4-phzInd+1, 0, 8, -1.2);
        % CREATE subplot axes
        sax(end+1) = axes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                                      fig.page.ypos(yind)+yOffSet+gYoffset,...
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
            [ehpcmpKDE,dxi] = ksdensity( egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,lat) );
            med             = median(    egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,lat) );
            plot(dxi,ehpcmpKDE,'-','color',hbaBin.color(hbaInd,:))
            [~,xi] = NearestNeighbour(dxi,med);
            line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hbaBin.color(hbaInd,:));
        end
        % FORMAT subplot
        xlim(sax(end),[-10,10])
        ylim(sax(end),[0,0.16])
        grid(sax(end),'on');
        sax(end).YTick = [0,0.05,0.10,0.15];
        sax(end).XTick = [-10,-5,0,5,10];
        if phzInd ~= 1,  sax(end).XTickLabel = {};                           end % not bottom of column
        if phzInd == 1,  xlabel(sax(end),'cm');                              end % bottom of column
        if phzInd == 1,  sax(end).YTickLabel = {'0','','0.1',''};            end % bottom of column
        if phzInd == 2,  ylabel(sax(end),'prob');                            end % middle of column        
    end
%%%>>>


% SUBPLOT -- AP POS DISTRIB -- partitioned by theta-phase and head-body-angle    
%%%<<<
    for phzInd = 1:phzBin.count
        % ADJUST subplot coordinates
        [yind, yOffSet, xind, xOffSet] = deal(4-phzInd+1, 0, 9, -1.2);
        % CREATE subplot axes
        sax(end+1) = axes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                                      fig.page.ypos(yind)+yOffSet+gYoffset,...
                                      fig.subplot.width*1.75,                   ...
                                      fig.subplot.height],                      ...
                          'FontSize', 8,                                        ...
                          'LineWidth',1);
        hold(sax(end),'on');
        % PLOT subplot
        % FORMAT subplot
        grid(sax(end),'on');
        % PLOT subplot
        for hbaInd = 1:hbaBin.count 
            [ehpcmpKDE,dxi] = ksdensity( egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,fwd) );
            med             = median(    egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,fwd) );
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
        
        if phzInd ~= 1,  sax(end).XTickLabel = {};                             end % not bottom of column
        if phzInd == 1,  xlabel(sax(end),'cm');                                end % bottom of column
        if phzInd == 1,  sax(end).XTickLabel = {'-10','','0','','10','','20'}; end % bottom of column
    end
%%%>>>



% $$$ 
% $$$ % SUBPLOTS -- LAT ZCR DISTRIB -- partitioned by theta-phase and head-body-angle
% $$$ %%%<<<
% $$$     for phzInd = 1:phzBin.count
% $$$         % ADJUST subplot coordinates
% $$$         [yind, yOffSet, xind, xOffSet] = deal(4-phzInd+1, 0, 11, -1.2);
% $$$         % CREATE subplot axes
% $$$         sax(end+1) = axes('Units','centimeters',                                ...
% $$$                           'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
% $$$                                       fig.page.ypos(yind)+yOffSet+gYoffset,...
% $$$                                       fig.subplot.width,                        ...
% $$$                                       fig.subplot.height],                      ...
% $$$                           'FontSize', 8,                                        ...
% $$$                           'LineWidth',1);
% $$$         hold(sax(end),'on');
% $$$         % PLOT subplot
% $$$         % FORMAT subplot
% $$$         grid(sax(end),'on');
% $$$         % PLOT subplot
% $$$         for hbaInd = 1:hbaBin.count
% $$$             [ehpcmpKDE,dxi] = ksdensity(egoHbaPhz.boot.(regions{region}).zscore(:,phzInd,hbaInd,2));
% $$$             plot(dxi,ehpcmpKDE,['-',hbaBin.color(hbaInd)])
% $$$             med = median(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,2));
% $$$             [~,xi] = NearestNeighbour(dxi,med);
% $$$             line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hbaBin.color(hbaInd));
% $$$         end
% $$$         % FORMAT subplot
% $$$         xlim(sax(end),[-25,25])
% $$$         ylim(sax(end),[0,0.16])
% $$$         grid(sax(end),'on');
% $$$         sax(end).YTick = [0,0.05,0.10,0.15];
% $$$         sax(end).XTick = [-10,-5,0,5,10];
% $$$         if p ~= 1,  sax(end).XTickLabel = {};                           end % not bottom of column
% $$$         if p == 1,  xlabel(sax(end),'cm');                              end % bottom of column
% $$$         if p == 1,  sax(end).YTickLabel = {'0','','0.1',''};            end % bottom of column
% $$$         if p == 2,  ylabel(sax(end),'prob');                            end % middle of column        
% $$$     end
% $$$ %%%>>>


% $$$ % SUBPLOT -- AP ZCR DISTRIB -- partitioned by theta-phase and head-body-angle    
% $$$ %%%<<<
% $$$     for phzInd = 1:phzBin.count
% $$$         % ADJUST subplot coordinates
% $$$         [yind, yOffSet, xind, xOffSet] = deal(4-phzInd+1, 0, 9, -1.2);
% $$$         % CREATE subplot axes
% $$$         sax(end+1) = axes('Units','centimeters',                                ...
% $$$                           'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
% $$$                                       fig.page.ypos(yind)+yOffSet+gYoffset,...
% $$$                                       fig.subplot.width*1.75,                        ...
% $$$                                       fig.subplot.height],                      ...
% $$$                           'FontSize', 8,                                        ...
% $$$                           'LineWidth',1);
% $$$         hold(sax(end),'on');
% $$$         % PLOT subplot
% $$$         % FORMAT subplot
% $$$         grid(sax(end),'on');
% $$$         % PLOT subplot
% $$$         for hbaInd = 1:hbaBin.count
% $$$             [ehpcmpKDE,dxi] = ksdensity(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,1));
% $$$             plot(dxi,ehpcmpKDE,['-',hbaBin.color(hbaInd)])
% $$$             med = median(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,1));
% $$$             [~,xi] = NearestNeighbour(dxi,med);
% $$$             line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hbaBin.color(hbaInd));
% $$$         end
% $$$         % FORMAT subplot
% $$$         xlim(sax(end),[-15,20])
% $$$         ylim(sax(end),[0,0.16])
% $$$         grid(sax(end),'on');
% $$$         sax(end).YTick = [0,0.05,0.10,0.15];
% $$$         sax(end).YTickLabel = {};        
% $$$         sax(end).XTick = [-10,-5,0,5,10,15,20];
% $$$         
% $$$         if p ~= 1,  sax(end).XTickLabel = {};                             end % not bottom of column
% $$$         if p == 1,  xlabel(sax(end),'cm');                                end % bottom of column
% $$$         if p == 1,  sax(end).XTickLabel = {'-10','','0','','10','','20'}; end % bottom of column
% $$$     end
% $$$ %%%>>>

% SUBPLOT -- LAT BOOT ZSCR ECDF -- partitioned by theta-phase and head-body-angle
%%%<<<
    for phzInd = 1:phzBin.count
        % ADJUST subplot coordinates
        [yind, yOffSet, xind, xOffSet] = deal(4-phzInd+1, 0, 11, -1.2);
        % CREATE subplot axes
        sax(end+1) = axes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                                      fig.page.ypos(yind)+yOffSet+gYoffset,...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height],                      ...
                          'FontSize', 8,                                        ...
                          'LineWidth',1);
        hold(sax(end),'on');
        samples = ':';
        % PLOT subplot
        xlim(sax(end),[-10,10]);
        ylim(sax(end),[0,1]);
        for hbaInd = 1:phzBin.count
            [F,X] = ecdf(egoHbaPhz.boot.(regions{region}).zscore(samples,phzInd,hbaInd,lat));
            plot(X,F,'Color',hbaBin.color(hbaInd,:));
        end
        % FORMAT subplot
        grid(sax(end),'on');
        xlim(sax(end),[-15,15]);
        ylim(sax(end),[0,1]);
        sax(end).XTick = [-10,-5,0,5,10];
        sax(end).XTickLabel = {'-10','','0','','10'};        
        sax(end).YTick = [0,0.25,0.5,0.75,1];
        sax(end).YTickLabel = {'0','','0.5','','1'};
        
        line(egoHbaPhz.perm.(regions{region}).sig.*[1,1],[0,1],'LineStyle','--','Color',[0.25,0.25,0.25]);
        line(-egoHbaPhz.perm.(regions{region}).sig.*[1,1],[0,1],'LineStyle','--','Color',[0.25,0.25,0.25]);
        if phzInd ~= 1,    sax(end).XTickLabel = {};                           end % not bottom of column
        if phzInd == 1,    xlabel(sax(end),'z-score');                              end % bottom of column
        %if p == 3,    title(sax(end),{'cdf','Lateral - Center'}); end % top of column
    end
%%%>>>


    
% SUBPLOT -- AP BOOT ZSCR ECDF -- partitioned by theta-phase and head-body-angle
%%%<<<
    for phzInd = 1:size(pfs{exampleUnit.trialIndex},1)
        % ADJUST subplot coordinates
        [yind, yOffSet, xind, xOffSet] = deal(4-phzInd+1, 0, 12, -1.2);
        % CREATE subplot axes
        sax(end+1) = axes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                                      fig.page.ypos(yind)+yOffSet+gYoffset,...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height],                      ...
                          'FontSize', 8,                                        ...
                          'LineWidth',1);
        hold(sax(end),'on');
        % PLOT subplot
        xlim(sax(end),[-10,10]);
        ylim(sax(end),[0,1]);
        for hbaInd = 1:hbaBin.count
            [F,X] = ecdf(egoHbaPhz.boot.(regions{region}).zscore(:,phzInd,hbaInd,fwd));
            plot(X,F,'Color',hbaBin.color(hbaInd,:));
        end
        
        % FORMAT subplot
        grid(sax(end),'on');
        xlim(sax(end),[-15,15]);
        ylim(sax(end),[0,1]);
        sax(end).XTick = [-10,-5,0,5,10];
        sax(end).YTick = [0,0.5,1];
        sax(end).YTick = [0,0.25,0.5,0.75,1];
        sax(end).YTickLabel = {};
        
        line(egoHbaPhz.perm.(regions{region}).sig.*[1,1],[0,1],'LineStyle','--','Color',[0.25,0.25,0.25]);
        line(-egoHbaPhz.perm.(regions{region}).sig.*[1,1],[0,1],'LineStyle','--','Color',[0.25,0.25,0.25]);
        if phzInd ~= 1,    sax(end).XTickLabel = {};                           end % not bottom of column
        if phzInd == 1,    xlabel(sax(end),'z-score');                              end % bottom of column
        %if p == 3,    title(sax(end),{'cdf','Lateral - Center'}); end % top of column
    end
%%%>>>

% SUBPLOTS -- LATERAL POS -- left vs Right lateral coordinatats for egoHba
%%%<<<    
    phzInd = phzBin.count;
    % ADJUST subplot coordinates
    [yind, yOffSet, xind, xOffSet] = deal(6, 0, 3+bins.hba.count, 0);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                        fig.page.ypos(yind)+yOffSet+gYoffset,...
                        fig.subplot.width*1.25,                        ...
                        fig.subplot.height*1.25],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
    %sax(end).tag = 
    % PLOT subplot
    boxplot( reshape(sqrt(egoHbaPhz.control.meanPos(unitsEgoCA1, phzInd, :, lat).^2 ...
                +(egoHbaPhz.control.meanPos(unitsEgoCA1, phzInd, :, fwd)-2).^2),[],1),...
            reshape(ones([numel(unitsEgoCA1),1])*[1,2,3],[],1),...
            'symbol',     'r.',...
            'plotstyle',  'traditional',...
             'orientation','vertical');
    
    sax(end).XTickLabel = mat2cell(hbaBin.key,1,ones([1,hbaBin.count]));
    xlabel(sax(end),'hb-angle');
    ylabel(sax(end),'radius (cm)');
    ylim(sax(end),[0,20]);
    title(sax(end),{'Ascending Theta Phase','Distance to Head'})
    % FORMAT subplot
    grid(sax(end),'on');
    %ylim(sax(end),[-15,20]);
% $$$     sax(end).XTick = [-10,-5,0,5,10];
% $$$     sax(end).YTick = [-5,0,5];
%%%>>>


% SUBPLOTS -- ASCENDING EGO ANGLE -- 
%%%<<<    
    phzInd = phzBin.count;
    for hbaInd = 1:hbaBin.count
        % ADJUST subplot coordinates
        [yind, yOffSet, xind, xOffSet] = deal(6, 0.5, hbaInd, (hbaInd-1));
        % CREATE subplot axes
        sax(end+1) = polaraxes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,...
                            fig.page.ypos(yind)+yOffSet+gYoffset,          ...
                            fig.subplot.width,                              ...
                            fig.subplot.height],                            ...
                          'FontSize', 8,                                        ...
                               'LineWidth',1);
        polarhistogram(sax(end),...
                             atan2(-egoHbaPhz.control.meanPos(unitsEgoCA1,3, hbaInd, lat), ...
                             egoHbaPhz.control.meanPos(unitsEgoCA1,3, hbaInd, fwd)),...
                       32,... bins
                       'FaceColor',hbaBin.color(hbaInd,:));
        hold(sax(end),'on');
        rticklabels(gca(),{'','','10','','20'})
        thetaticklabels(gca(),{'0','','','90','','','180','','','270','','',})
        if hbaInd==2
            title(sax(end),{'Ascending Theta Phase: Head to EgoField Angle','',hbaBin.label{hbaInd}});
        else
            title(sax(end),...
                  hbaBin.label(hbaInd));
        end
        tids = {3:5,[6,7,27],18:25,29};
        for tid = 1:4
            uinds = ismember(egoCluSessionMap(:,1),tids{tid});
            polarplot(sax(end),...
                      -circ_mean(atan2(egoHbaPhz.control.meanPos(uinds,3,hbaInd,2),...
                                      egoHbaPhz.control.meanPos(uinds,3,hbaInd,1))).*[1,1],...
                      [0,20],...
                      '-',...
                      'LineWidth',2);
        end
    end
    %%%>>>

    
end


nb =4;
mBinHvl.edges = [-40,-5,5,40];%linspace(-40,40,nb);
mBinHvl.centers = mean([mBinHvl.edges(1:end-1);mBinHvl.edges(2:end)]);
mBinHvl.count = numel(mBinHvl.edges)-1;
mBinHba.edges = [-1.2, -0.2, 0.2, 1.2];%linspace(-1.2,1.2,nb);
mBinHba.centers = mean([mBinHba.edges(1:end-1);mBinHba.edges(2:end)]);
mBinHba.count = numel(mBinHba.edges)-1;
mout = zeros([mBinHba.count,mBinHvl.count]);
cout = zeros([mBinHba.count,mBinHvl.count]);
for a = 1:mBinHba.count
    for v = 1:mBinHvl.count
        ind =   WithinRanges(decoded.phz,phzBin.edges([3,4])) ...
              & WithinRanges(-decoded.hvl,mBinHvl.edges([v,v+1])) ...
              & WithinRanges(decoded.hba,mBinHba.edges([a,a+1])) ...              
              & randn(size(decoded.hba))>0;
        mout(a,v) = mean(decoded.clat(ind),'omitnan');
        cout(a,v) = sum(ind);
    end
end



% SUBPLOT  - New Axes - y   yo    x   xo       gyo       gxo     ws    hs
sax(end+1) = setup_axes(fig, 6,   0,   8,   0, gYoffset, gXoffset,  1.25, 1.25);
% PLOT subplot
imagesc(sax(end),            ...
        mBinHba.centers,     ...
        mBinHvl.centers,     ...
        mout',               ...
        [-50,50]);
axis(sax(end),'xy');
axis(sax(end),'tight');
colormap('jet');
sax(end).XTick = mBinHba.centers;
sax(end).YTick = mBinHvl.centers ;   
sax(end).XTickLabel = mat2cell(hbaBin.key,1,ones([1,hbaBin.count]));
sax(end).YTickLabel = mat2cell(hbaBin.key,1,ones([1, ...
                    hbaBin.count]));
xlabel(sax(end),'Lat head speed (cm/s)');
ylabel(sax(end),'HB Angle (rad)');
title(sax(end),{'Mean Decoded','Position'})
cax = colorbar(sax(end));
ylabel(cax,'Mean Lateral Position');
cax.Units = 'centimeters';
cax.Position(1) = cax.Position(1)+1;




% SUBPLOT  - New Axes -      y   yo    x   xo       gyo       gxo     ws    hs
sax(end+1) = setup_axes(fig, 6,   0,  10,   1, gYoffset, gXoffset,  1.25, 1.25);
for hbaInd = 1:hbaBin.count
    ind = WithinRanges(decoded.phz,[4.5,5.5]) ...
          & randn(size(decoded.hba))>0 ...
          & WithinRanges(decoded.hba,hbaBin.edges(hbaInd+[0,1]));      
    [F,xi] = ksdensity(decoded.clt(ind)/10);
    plot(sax(end),xi,F,'-','Color',hbaBin.color(hbaInd,:),'LineWidth',1);
end
xlim(sax(end),[-40,40]);
lgd = legend(mat2cell(bins.hba.key,1,ones(size(bins.hba.key))),...
             'Location','NorthEastOutside');
lgd.Units = 'centimeters';
lgd.Position(1) = lgd.Position(1)+2
xlabel(sax(end),'Lat Head Pos (cm)')
ylabel(sax(end),'Prob')
title(sax(end),{'Mean Decoded','Position by Hba'})
%endfig


figure,
phzInd = 3;
hold('on')
for hbaInd = 1:hbaBin.count
    plot(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,2),egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,1),'.','Color',hbaBin.color(hbaInd))
end

figure,
for hbaInd = 1:hbaBin.count
    for phzInd = 1:phzBin.count
        subplot2(3,3,phzBin.count+1-phzInd,hbaInd);
        histogram(sqrt(egoHbaPhz.control.size(unitsEgoCA1,phzInd,hbaInd))./pi*2,linspace([0,10,20]));;
    end
end

figure,
ind = zeros([size(egoSizeHba,1),1]);
ind(uidsCA3) = 1;
ind = ind & nniz(egoSizeHba);
subplot(211);
histogram(nonzeros(egoSizeHba(ind,3,:))-nonzeros(egoSizeHba(ind,2,:)),linspace(-300,300,20))
subplot(212);
histogram(nonzeros(egoSizeHba(ind,2,:))-nonzeros(egoSizeHba(ind,1,:)),linspace(-300,300,20))


figure,
ind = zeros([size(egoSizeHba,1),1]);
ind(unitsEgoCA1) = 1;
ind = ind & nniz(egoSizeHba);
subplot(211);
histogram(mean(egoSizeHba(ind,3,:),3)-egoSize(ind,3),linspace(-300,300,20))
subplot(212);
histogram(mean(egoSizeHba(ind,2,:),3)-egoSize(ind,2),linspace(-300,300,20))


figure,
hold('on');
plot(egoSize(ind,3),mean(sq(egoSizeHba(ind,3,:)),2),'.');
line([0,800],[0,800]);
plot(egoSize(ind,2),mean(sq(egoSizeHba(ind,2,:)),2),'.g');
line([0,800],[0,400]);

figure,
hold('on');
histogram(egoSize(ind,3)-egoSize(ind,2),linspace(-300,300,20));
histogram(mean(egoSizeHba(ind,3,:),3)-mean(egoSizeHba(ind,2,:),3),linspace(-300,300,20));

figure,
subplot(211);histogram((egoSize(ind,3)-egoSize(ind,2))./(egoSize(ind,3)+egoSize(ind,2)).*100,linspace(-100,100,20));
subplot(212);histogram((median(egoSizeHba(ind,3,:),3)-median(egoSizeHba(ind,2,:),3))./(median(egoSizeHba(ind,3,:),3)+median(egoSizeHba(ind,2,:),3)).*100,linspace(-100,100,20));

[H,P] = ttest2(egoSize(ind,3)-egoSize(ind,2),mean(egoSizeHba(ind,3,:),3)-mean(egoSizeHba(ind,2,:),3));

figure,imagesc(log10(sq(mean(egoMeanRmapRateHba(unitsEgoCA1,:,:),3)))')


figure,
plot(sq(mean(egoMeanRmapRateHba(unitsEgoCA1,1,:),3))./sq(mean(egoMeanRmapRateHba(unitsEgoCA1,2,:),3)),...
     sq(mean(egoMeanRmapRateHba(unitsEgoCA1,3,:),3))./sq(mean(egoMeanRmapRateHba(unitsEgoCA1,2,:),3)),...
        '.')

figure,
plot(egoMeanRmapRate(unitsEgoCA1,2)./egoMeanRmapRate(unitsEgoCA1,3)),...
     sq(mean(egoMeanRmapRateHba(unitsEgoCA1,3,:),3))./sq(mean(egoMeanRmapRateHba(unitsEgoCA1,2,:),3)),...
    '.')




figure
for p = 1:3
    subplot2(3,1,p,1)
    plot(pfsh{exampleUnit.trialIndex}{4-p,1},exampleUnit.id,1,'', ...
         [0,exampleUnit.maxRate],'colorMap',@jet,'mazeMaskFlag',false,'flipAxesFlag',true);
end


%% PLOTING EXAMPLES
lims = {[-250,250],[-250,250]};
figure();
al = 1:5;
numAng = numel(al);
sax = gobjects([0,1]);
for p = 3
    rmap = plot(pfsh{20}{p,1},25,19);
    for a = 1:numAng,
        sax(end+1) = subplot2(1,5,1,a);
        pcolor(pfsh{20}{p,1}.adata.bins{1},...
                    pfsh{20}{p,1}.adata.bins{2},...
                    rmap(:,:,al(a)));
        
        caxis   (sax(end),[0,10]);
        colormap(sax(end),'jet');
        shading (sax(end),'flat');
        axis    (sax(end),'xy');
        xlim    (sax(end),lims{1});
        ylim    (sax(end),lims{2});        
        
        Lines([],0,'k');
        Lines(0,[],'k');
        
        set(sax(end),'XTick',[]);
        set(sax(end),'YTick',[]);        
        
        % ADD subject
% $$$         if p %== 4,
% $$$             subject = struct(rat);
% $$$             subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
% $$$             subject = update_subject_patch(subject,'body', numAng+1-a,  true,hbaBinEdg,hbaBinCtr);
% $$$             patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$             patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$             patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
% $$$             line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
% $$$             line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
% $$$         end
    end
end



%% PLOTING EXAMPLES
lims = {[-300,300],[-300,300]};
figure();
sax = gobjects([0,1]);
trialInd = 20;
unit = 2;
for a = 1:egoHba.hbaBin,
    sax(end+1) = subplot2(1,3,1,a);
    pcolor(egoHba.xbins,...
           egoHba.ybins,...
           egoHba.rmap{trialInd}(:,:,unit,a));
    
    caxis   (sax(end),[0,10]);
    colormap(sax(end),'jet');
    shading (sax(end),'flat');
    axis    (sax(end),'xy');
    xlim    (sax(end),lims{1});
    ylim    (sax(end),lims{2});        
    
    Lines([],0,'k');
    Lines(0,[],'k');
    
    set(sax(end),'XTick',[]);
    set(sax(end),'YTick',[]);        
    
    % ADD subject
% $$$         if p %== 4,
% $$$             subject = struct(rat);
% $$$             subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
% $$$             subject = update_subject_patch(subject,'body', numAng+1-a,  true,hbaBinEdg,hbaBinCtr);
% $$$             patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$             patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$             patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
% $$$             line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
% $$$             line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
% $$$         end
end



