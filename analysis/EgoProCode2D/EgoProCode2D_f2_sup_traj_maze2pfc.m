




[hfig,fig,fax,sax] = set_figure_layout(figure(666013),'A4', 'landscape',[],1.8,1.8,0.2,0.4);
gxoff = 0;
gyoff = 0;
lat = 2;

[scaleW, scaleH] = deal(1,1);
for phzInd = 1:bins.phz.count
    [yind, yoff, xind, xoff] = deal((4-phzInd)+3, 0, 12, -1.2);
    sax(end+1) = setup_axes(fig,yind, yoff, xind, xoff, gyoff, gxoff, scaleW, scaleH);
    
    grid(sax(end),'on');
    for hbaInd = 1:bins.hba.count
        [ehpcmpKDE,dxi] = ksdensity( egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,lat) );
        med             = median(    egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,lat) );
        plot(dxi,ehpcmpKDE,'-','color',bins.hba.color(hbaInd,:))
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',bins.hba.color(hbaInd,:));
    end
    
    ylim(sax(end),[0,0.16]);  sax(end).YTick = [ 0, 0.05, 0.10, 0.15];
    xlim(sax(end),[-15,15]);  sax(end).XTick = [-15,-10,-5,0,5,10,15];
                              sax(end).XTickLabel = {'','-10','','0','','10',''};
    grid(sax(end),'on');
    if phzInd ~= 1,  sax(end).XTickLabel = {};                           end % not bottom
    if phzInd == 1,  xlabel(sax(end),'cm');                              end % bottom
    if phzInd == 2,  ylabel(sax(end),'prob');                            end % middle
    if phzInd == 3,
        title(sax(end),{'Lateral','Mean Field','Position'});   
        lax = legend(gca(),{'L','L','C','C','R','R'},'Location','EastOutside')
        lax.Units = 'centimeters';		 
        lax.Position(1) = lax.Position(1)+2;	         
    end % top
    
end



% CA1 
uids = unitsEgoCA1;
globalXOffset = 0;
globalYOffset = 0;
exampleUnit.trialIndex = 18;
exampleUnit.close.Xlims = [-200,400];
exampleUnit.close.Ylims = [-400,200];

exampleUnit.id = 42;
exampleUnit.maxRate = 30;
exampleUnit.index = find(unitsEgo{exampleUnit.trialIndex}==exampleUnit.id);
exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];
% SUBPLOTS -- EGO FIELD -- partitioned by theta-phase and head-body-angle
    for phzInd = 1:bins.phz.count
        for hbaInd = 1:bins.hba.count
            
            % CREATE subplot axes            
            [yind, yOffSet, xind, xOffSet] = deal((bins.phz.count-phzInd+1)+3, 0, hbaInd+2, 0);
            sax(end+1) = setup_axes(fig,yind, yoff, xind, xoff, gyoff, gxoff, scaleW, scaleH);

            shading(sax(end),'flat');
            set(pcolor(egoHbaPhzRmaps.xpos-diff(egoHbaPhzRmaps.xpos(1:2))/2,...
                       egoHbaPhzRmaps.ypos-diff(egoHbaPhzRmaps.ypos(1:2))/2,...
                       fliplr(rot90(...
                           egoHbaPhzRmaps.rmap {exampleUnit.trialIndex} (:,:,exampleUnit.index,phzInd,hbaInd)',-1)) ...
                       .*egoHbaPhzRmaps.mask),'EdgeColor','none');
            axis(sax(end),'xy');
            colormap(sax(end),'jet');
            caxis(sax(end),[0,exampleUnit.maxRate]);
            
            xlim(sax(end),[-250,250]);  sax(end).XTickLabel =[];
            ylim(sax(end),[-250,250]);  sax(end).YTickLabel =[];
            daspect(sax(end),[1,1,1]);
            box(sax(end),'on');
            Lines([],0,'w');
            Lines(0,[],'w');
            subject = struct(rat);
            subject = update_subject_patch(subject, 'head',...
                                           [], false,...
                                           bins.hba.edges,...
                                           bins.hba.centers);
            subject = update_subject_patch(subject, 'body',...
                                           bins.hba.count+1-hbaInd,  true,...
                                           bins.hba.edges,...
                                           bins.hba.centers);
            patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
            patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
            patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
            line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
            line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
            axes(fax);
        end
    end
%%%>>>



% THETA Phase Vertical
%%%<<<
[yind, yOffSet, xind, xOffSet] = deal(6, 0, 1, 0.6);
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


exampleUnit.trialIndex = 20;
exampleUnit.close.Xlims = [-200,400];
exampleUnit.close.Ylims = [-400,200];
exampleUnit.id = 20;
exampleUnit.maxRate = 40;
exampleUnit.index = find(unitsEgo{exampleUnit.trialIndex}==exampleUnit.id);
exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];
% SUBPLOTS -- EGO FIELD -- partitioned by theta-phase and head-body-angle
    for phzInd = 1:bins.phz.count
        for hbaInd = 1:bins.hba.count
            % CREATE subplot axes            
            [yind, yOffSet, xind, xOffSet] = deal(bins.phz.count-phzInd+1+3, 1, hbaInd+6, 0);
            sax(end+1) = setup_axes(fig,yind, yoff, xind, xoff, gyoff, gxoff, scaleW, scaleH);
            
            shading(sax(end),'flat');
            set(pcolor(egoHbaPhzRmaps.xpos-diff(egoHbaPhzRmaps.xpos(1:2))/2,...
                       egoHbaPhzRmaps.ypos-diff(egoHbaPhzRmaps.ypos(1:2))/2,...
                       fliplr(rot90(...
                           egoHbaPhzRmaps.rmap {exampleUnit.trialIndex} (:,:,exampleUnit.index,phzInd,hbaInd)',-1)) ...
                       .*egoHbaPhzRmaps.mask),'EdgeColor','none');
            axis(sax(end),'xy');
            colormap(sax(end),'jet');
            caxis(sax(end),[0,exampleUnit.maxRate]);
            
            sax(end).XTickLabel =[];
            sax(end).YTickLabel =[];
            xlim(sax(end),[-250,250]);
            ylim(sax(end),[-250,250]);
            daspect(sax(end),[1,1,1]);
            box(sax(end),'on');
            Lines([],0,'w');
            Lines(0,[],'w');
            subject = struct(rat);
            subject = update_subject_patch(subject, 'head',...
                                           [], false,...
                                           bins.hba.edges,...
                                           bins.hba.centers);
            subject = update_subject_patch(subject, 'body',...
                                           bins.hba.count+1-hbaInd,  true,...
                                           bins.hba.edges,...
                                           bins.hba.centers);
            patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
            patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
            patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
            line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
            line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
            axes(fax);
        end
    end
%%%>>>






%% Ego Forward
%% HEAD BODY ANGLE ---------------------------------------------------------------------------------
medD = [];
skwD = [];
stdD = [];
medR = [];
skwR = [];
stdR = [];
medL = [];
skwL = [];
stdL = [];
medC = [];
skwC = [];
stdC = [];
sectors = linspace(-pi,pi,11);
rdists = 5:5:45;
for r = 1:numel(rdists)
    clear('xcomp','ycomp','zcomp','ccomp');
    xcomp.label = 'hba';            xcomp.data = [];    xcomp.edgs = [-1.25,-0.2,0.2,1.25];
    ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5,6 );
    ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-80,80];
    fcomp.data = [];
    mcomp.data = [];
    for t = [1:3,5:8,11],
        %for t = [6]    
        dc = dca{t};
        
        headAngle = sq(dca{t}.xyz(:,'nose',[1,2])-dca{t}.xyz(:,'hcom',[1,2]));
        headAngle = atan2(headAngle(:,2),headAngle(:,1));
        
        mazeAngle =  sq(dca{t}.xyz(:,'hcom',[1,2]));
        mazeAngle = atan2(mazeAngle(:,2),mazeAngle(:,1));
        
        headMazeAngle = circ_dist(headAngle, mazeAngle);
        
        mind =  dc.stcm(:,1)==1                                             ... Theta
                & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5)         ... Walk,Turn,Pause              
                & dc.hvfl(:,1)>-2                                           ... FwdHeadSpeed
                & dc.ucnt>=4 & dc.ucnt<=10                                    ... UnitCount
                & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)<rdists(r)+5 ... DistanceFromMazeCenter
                & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)>rdists(r)-5;

        xcomp.data = cat(1, xcomp.data, dc.hbang(mind,1));
        ycomp.data = cat(1, ycomp.data, dc.phz(mind));
        
        ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)./10);
        %ccomp.data = cat(1, ccomp.data, dc.esax(mind,1)./10);
        mcomp.data = cat(1, mcomp.data, headMazeAngle(mind));
    end
    
    [xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,ccomp);
    zmean   = zcomp.mean;    zmean(zcomp.count<5)  = nan;
    zmedian = zcomp.median;  zmean(zcomp.count<5)  = nan;
    zstd    = zcomp.std;     zstd(zcomp.count<5)   = nan;
    zcount  = zcomp.count;   zcount(zcomp.count<5) = nan;
    sectors = linspace(-pi,pi,11);
    
    sectorsL = linspace(-pi,pi,11);
    for s = 1:numel(sectors)-1
        if r>1
            bind = mcomp.data<sectors(s+1) & mcomp.data>sectors(s) ...
                   & ycomp.data>4.5& ycomp.data<5;
        else
            bind = ycomp.data>4.5& ycomp.data<5;
        end
        

        ind  = bind & xcomp.data<-0.2;        
        medR(r,s) = median(ccomp.data(ind,1));%-rmodel(r,s);
        skwR(r,s) = skew(ccomp.data(ind,1));
        stdR(r,s) = std(ccomp.data(ind,1));
        
        ind = bind & xcomp.data<0.2& xcomp.data>-0.2;
        medC(r,s) = median(ccomp.data(ind,1));%-rmodel(r,s);
        skwC(r,s) = skew(ccomp.data(ind,1));
        stdC(r,s) = std(ccomp.data(ind,1));
        
        ind = bind & xcomp.data>0.2;
        medL(r,s) = median(ccomp.data(ind,1));%-rmodel(r,s);
        skwL(r,s) = skew(ccomp.data(ind,1));
        stdL(r,s) = std(ccomp.data(ind,1));
        
        medD(r,s) = medR(r,s)-medL(r,s);
        skwD(r,s) = skwR(r,s)-skwL(r,s);
        stdD(r,s) = stdR(r,s)-stdL(r,s);
    end
end

sectorc = mean([sectors(2:end);sectors(1:end-1)]);
rdistc = mean([rdists(2:end);rdists(1:end-1)]);



rdiste = [0,rdists+2.5];

%[THETA,RR] = meshgrid(sectors,[rdists,rdists(end)]);
[THETA,RR] = meshgrid(sectors,[rdiste]);
[THETAC,RRC] = meshgrid(sectorc,rdists);

%[A,B] = pol2cart(circ_dist(THETA,diff(sectors([1,2]))),RR);
[X,Y] = pol2cart(THETA,RR);
[Xc,Yc] = pol2cart(THETAC,RRC);


[scaleW, scaleH] = deal(3,3);
[yind, yoff, xind, xoff] = deal(3, 0, 1, 0);
sax(end+1) = setup_axes(fig,yind, yoff, xind, xoff, gyoff, gxoff, scaleW, scaleH);
title('Left HBA')

sax(end).Color='none';
sax(end).XAxis.Visible='off';
sax(end).YAxis.Visible='off';
for t = 1:numel(sectors)-1
    for r = 1:numel(rdists)-1
        patch([X(r,t),X(r,t+1),X(r+1,t+1),X(r+1,t)],...
              [Y(r,t),Y(r,t+1),Y(r+1,t+1),Y(r+1,t)],...
              medL(r,t));
    end
end
colormap(gca,'jet');
caxis([-15,15]);

circle(0, 55, 9, 'k');
circular_arrow(sax(end),5,[0,55],0, 1.5*pi, -1, 'r', [3,5],'vback2');
plot(0, 55,'*k','MarkerSize',5);
text(sax(end), 0, 70, 'CCW Movement','HorizontalAlignment','center');


circle(0, -55, 9, 'k');
circular_arrow(sax(end),5,[0,-55],0, 1.5*pi, 1, 'r', [3,5],'vback2');
plot(0,-55,'*k','MarkerSize',5);
text(sax(end), 0, -70, 'CW Movement','HorizontalAlignment','center');

circle(-55, 0, 9, 'k');
har = annotation('arrow');
har.Parent = sax(end);
har.Position = [-64,0,9,0];
har.HeadLength = 5;
har.HeadWidth = 3;
har.HeadStyle = 'vback2';
har.LineWidth = 1;
har.Color = 'r';
plot(-55,0,'*k','MarkerSize',5);
text(sax(end), -50, 22, {'To','Center'},'HorizontalAlignment','right');

circle(55, 0, 9, 'k');
har1 = annotation('arrow');
har1.Parent = sax(end);
har1.Position = [55,0,9,0];
har1.HeadLength = 5;
har1.HeadWidth = 3;
har1.HeadStyle = 'vback2';
har1.LineWidth = 1;
har1.Color = 'r';
plot(55,0,'*k','MarkerSize',5);
text(sax(end), 50, 22, {'From','Center'},'HorizontalAlignment','left');
xlim(sax(end),[-85,85]);
ylim(sax(end),[-85,85]);




[scaleW, scaleH] = deal(3,3);
[yind, yoff, xind, xoff] = deal(3, 0, 4, 0);
sax(end+1) = setup_axes(fig,yind, yoff, xind, xoff, gyoff, gxoff, scaleW, scaleH);
title('Straight HBA')
sax(end).Color='none';
sax(end).XAxis.Visible='off';
sax(end).YAxis.Visible='off';
for t = 1:numel(sectors)-1
    for r = 1:numel(rdists)-1
        patch([X(r,t),X(r,t+1),X(r+1,t+1),X(r+1,t)],...
              [Y(r,t),Y(r,t+1),Y(r+1,t+1),Y(r+1,t)],...
              medC(r,t));
    end
end
colormap(sax(end),'jet');
caxis(sax(end),[-15,15]);
xlim(sax(end),[-85,85]);
ylim(sax(end),[-85,85]);

[scaleW, scaleH] = deal(3,3);
[yind, yoff, xind, xoff] = deal(3, 0, 7, 0);
sax(end+1) = setup_axes(fig,yind, yoff, xind, xoff, gyoff, gxoff, scaleW, scaleH);
title('Right HBA')
sax(end).Color='none';
sax(end).XAxis.Visible='off';
sax(end).YAxis.Visible='off';
for t = 1:numel(sectors)-1
    for r = 1:numel(rdists)-1
        patch([X(r,t),X(r,t+1),X(r+1,t+1),X(r+1,t)],...
              [Y(r,t),Y(r,t+1),Y(r+1,t+1),Y(r+1,t)],...
              medR(r,t));
    end
end
colormap(sax(end),'jet');
caxis(sax(end),[-15,15]);
xlim(sax(end),[-85,85]);
ylim(sax(end),[-85,85]);



[scaleW, scaleH] = deal(3,3);
[yind, yoff, xind, xoff] = deal(3, 0, 10, 0);
sax(end+1) = setup_axes(fig,yind, yoff, xind, xoff, gyoff, gxoff, scaleW, scaleH);
title({'Decoded','Right Minus Left HBA'})
sax(end).Color='none';
sax(end).XAxis.Visible='off';
sax(end).YAxis.Visible='off';
for t = 1:numel(sectors)-1
    for r = 1:numel(rdists)-1
        patch([X(r,t),X(r,t+1),X(r+1,t+1),X(r+1,t)],...
              [Y(r,t),Y(r,t+1),Y(r+1,t+1),Y(r+1,t)],...
              medR(r,t)-medL(r,t));
    end
end
colormap(sax(end),'jet');
caxis(sax(end),[-15,15]);
xlim(sax(end),[-85,85]);
ylim(sax(end),[-85,85]);
cax = colorbar(sax(end));
cax.Units = 'centimeters';
cax.Position(1) = cax.Position(1)+1;
cax.Position(4) = 3;
cax.Position(2) = cax.Position(2)+1;
cax.Position(3) = 0.25;
ylabel(cax,'cm');

% endfig

% $$$ 
% $$$ 
% $$$ figure();
% $$$ ax = subplot(1,1,1);
% $$$ hold(ax,'on');
% $$$ cmedD = medL; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
% $$$ surf(A,B,cmedD,'edgecolor','none'),view(opt.view);colormap('jet');cax=colorbar();caxis(opt.clim);
% $$$ 
% $$$ daspect([1,1,1])
% $$$ 
% $$$ figure();
% $$$ for phzInd = 1:bins.phz.count
% $$$     subplot(3,1,phzInd);
% $$$     hold('on');
% $$$     for hbaInd = 1:bins.hba.count
% $$$         histogram(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,1),linspace([-20,20,40]),'FaceColor',bins.hba.color(hbaInd,:))
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ figure();
% $$$ for phzInd = 1:bins.phz.count
% $$$     subplot(3,1,phzInd);
% $$$     hold('on');
% $$$     for hbaInd = 1:bins.hba.count
% $$$         histogram(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,2),linspace([-20,20,40]),'FaceColor',bins.hba.color(hbaInd,:))
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure();
% $$$ for phzInd = 1:bins.phz.count
% $$$     subplot(3,1,phzInd);
% $$$     hold('on');
% $$$     for hbaInd = 1:bins.hba.count
% $$$         histogram(egoHbaPhz.control.meanPos(unitsEgoCA3,phzInd,hbaInd,1),linspace([-20,20,40]),'FaceColor',bins.hba.color(hbaInd,:))
% $$$     end
% $$$ end

% $$$ 
% $$$ figure();
% $$$ for phzInd = 1:bins.phz.count
% $$$     subplot(3,1,phzInd);
% $$$     hold('on');
% $$$     for hbaInd = 1:bins.hba.count
% $$$         histogram(egoHbaPhz.control.meanPos(unitsEgoCA,phzInd,hbaInd,2),linspace([-20,20,40]),'FaceColor',bins.hba.color(hbaInd,:))
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ [t,p] = ttest2(egoHbaPhz.control.meanPos(unitsEgoCA3,3,2,2), ...
% $$$                egoHbaPhz.control.meanPos(unitsEgoCA3,3,1,2))
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ [hfig,fig,fax,sax] = set_figure_layout(figure(666014),'A4', 'landscape',[],1.8,1.8,0.2,0.4);
% $$$ gxoff = 0;
% $$$ gyoff = 0;
% $$$ lat = 2;
% $$$ 
% $$$ uids = unitsEgoCA1;
% $$$ globalXOffset = 0;
% $$$ globalYOffset = 0;
% $$$ exampleUnit.trialIndex = 18;
% $$$ exampleUnit.close.Xlims = [-200,400];
% $$$ exampleUnit.close.Ylims = [-400,200];
% $$$ 
% $$$ exampleUnit.id = 42;
% $$$ exampleUnit.maxRate = 30;
% $$$ exampleUnit.index = find(unitsEgo{exampleUnit.trialIndex}==exampleUnit.id);
% $$$ exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];
% $$$ 
% $$$ [scaleW, scaleH] = deal(1,1);
% $$$ for phzInd = 1:bins.phz.count
% $$$     for hbaInd = 1:bins.hba.count
% $$$         
% $$$         % CREATE subplot axes            
% $$$         [yind, yOffSet, xind, xOffSet] = deal((bins.phz.count-phzInd+1)+3, 0, hbaInd+2, 0);
% $$$         sax(end+1) = setup_axes(fig,yind, yoff, xind, xoff, gyoff, gxoff, scaleW, scaleH);
% $$$ 
% $$$         shading(sax(end),'flat');
% $$$         set(pcolor(egoHbaPhzRmaps.xpos-diff(egoHbaPhzRmaps.xpos(1:2))/2,...
% $$$                    egoHbaPhzRmaps.ypos-diff(egoHbaPhzRmaps.ypos(1:2))/2,...
% $$$                    fliplr(rot90(...
% $$$                        egoHbaPhzRmaps.rmap {exampleUnit.trialIndex} (:,:,exampleUnit.index,phzInd,hbaInd)',-1)) ...
% $$$                    .*egoHbaPhzRmaps.mask),'EdgeColor','none');
% $$$         axis(sax(end),'xy');
% $$$         colormap(sax(end),'jet');
% $$$         caxis(sax(end),[0,exampleUnit.maxRate]);
% $$$         
% $$$         xlim(sax(end),[-250,250]);  sax(end).XTickLabel =[];
% $$$         ylim(sax(end),[-250,250]);  sax(end).YTickLabel =[];
% $$$         daspect(sax(end),[1,1,1]);
% $$$         box(sax(end),'on');
% $$$         Lines([],0,'w');
% $$$         Lines(0,[],'w');
% $$$         subject = struct(rat);
% $$$         subject = update_subject_patch(subject, 'head',...
% $$$                                        [], false,...
% $$$                                        bins.hba.edges,...
% $$$                                        bins.hba.centers);
% $$$         subject = update_subject_patch(subject, 'body',...
% $$$                                        bins.hba.count+1-hbaInd,  true,...
% $$$                                        bins.hba.edges,...
% $$$                                        bins.hba.centers);
% $$$         patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$         patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$         patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
% $$$         line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
% $$$         line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
% $$$         axes(fax);
% $$$     end
% $$$ end
% $$$ 
