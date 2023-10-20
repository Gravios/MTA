
[hfig,fig,fax,sax] = set_figure_layout(figure(666003),'A4','portrait',[],1.5,1.5,0.1,0.1);


exampleUnit.trialIndex = 20;
exampleUnit.close.Xlims = [-200,400];
exampleUnit.close.Ylims = [-400,200];
exampleUnit.id = 21;
exampleUnit.maxRate = 16;
exampleUnit.index = find(unitsEgo{exampleUnit.trialIndex}==exampleUnit.id);
exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];
uids = unitsEgoCA1;
globalXOffset = 0;
globalYOffset = 0;



exampleUnit.maxRate = 20;



globalXOffset = -0.8;
globalYOffset = 0;

for phzInd = 1:phzBin.count
    for hvfInd = 1:hvfBin.count-1
        
        %%%<<< PLOT egoField (phz x hvf)
        % ADJUST subplot coordinates
        [yind, yOffSet, xind, xOffSet] = deal(phzBin.count-phzInd+1, phzInd/3-1.5, hvfInd+1, 0);        
        % CREATE subplot axes
        sax(end+1) = axes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                                      fig.page.ypos(yind)+yOffSet+globalYOffset,...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height*1.2],                      ...
                          'FontSize', 8,                                        ...
                          'LineWidth',1);
        hold(sax(end),'on');
        shading(sax(end),'flat');
        set(pcolor(egoHvfPhzRmaps.xpos-diff(egoHvfPhzRmaps.xpos(1:2))/2,...
                   egoHvfPhzRmaps.ypos-diff(egoHvfPhzRmaps.ypos(1:2))/2,...
                   fliplr(rot90(egoHvfPhzRmaps.rmap{exampleUnit.trialIndex}(:,:,exampleUnit.index,phzInd,hvfInd+1)',-1)).*egoHvfPhzRmaps.mask),'EdgeColor','none');
        axis(sax(end),'xy');
        colormap(sax(end),'jet');
        caxis(sax(end),[0,exampleUnit.maxRate]);
        
        % FORMAT subplot
        sax(end).XTickLabel =[];
        sax(end).YTickLabel =[];
        xlim(sax(end),[-250,250]);
        ylim(sax(end),[-300,300]);
        daspect(sax(end),[1,1,1]);
        box(sax(end),'on');
        Lines([],0,'w');
        Lines(0,[],'w');
        % ANNOTATE 
% $$$         subject = struct(rat);
% $$$         subject = update_subject_patch(subject, 'head',...
% $$$                                        [], false,...
% $$$                                        hvfBin.edges,...
% $$$                                        hvfBin.centers);
% $$$         subject = update_subject_patch(subject, 'body',...
% $$$                                        hvfBin.count+1-hvfInd,  true,...
% $$$                                        hvfBin.edges,...
% $$$                                        hvfBin.centers);
% $$$         patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$         patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$         patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
% $$$         line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
% $$$         line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
        axes(fax);
        rectangle('Position',sax(end).Position, ...
                  'EdgeColor',phzBin.color(phzInd,:), ...
                  'FaceColor','None',...
                  'LineWidth',1);
        if phzInd==3,
            line([sax(end).Position(1),sum(sax(end).Position([1,3]))],...
                  sum(sax(end).Position([2,4])).*[1,1]+0.1,...
                 'LineWidth',2,...
                 'Color',hvfBin.color(hvfInd+1,:));
            title(sax(end),{hvfBin.label{hvfInd+1},' '});
        end
    end
end
%%%>>>    

% THETA Phase Vertical
%%%<<<
[yind, yOffSet, xind, xOffSet] = deal(3, phzInd/3-1.5, 1, 0.6);
subplotWidth = 0.8;
subplotHeight = fig.subplot.height *1.2 * 3 + fig.subplot.verticalPadding * 2
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
ylabel(sax(end),'Theta');
sax(end).XAxis.Visible = 'off';
sax(end).YAxis.Color = 'k';
sax(end).YTick =[];
sax(end).Color = 'none';
%%%>>>





%globalXOffset = 1;
% SUBPLOT -- LATERAL POS -- left vs Right lateral coordinatats for egoHba
%%%<<<
% ADJUST subplot coordinates
for phzInd = 1:phzBin.count
    [yind, yOffSet, xind, xOffSet] = deal(phzBin.count+1-phzInd+8, 0, 6, 0);
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
    xlim(sax(end),[-15,15]);
    ylim(sax(end),[-15,15]);
    Lines([],0,'k');
    Lines(0,[],'k');
    plot(egoHvfPhz.control.meanPos( unitsEgoHvfCA1, phzInd, 2, fwd) ,...
         egoHvfPhz.control.meanPos( unitsEgoHvfCA1, phzInd, 4, fwd) ,...
         '.',                                                   ...
         'MarkerFaceColor',phzBin.color(phzInd,:),              ...
         'MarkerEdgeColor',phzBin.color(phzInd,:));
    % FORMAT subplot
    grid(sax(end),'on');
    xlim(sax(end),[-15,15]);
    ylim(sax(end),[-15,15]);
    sax(end).XTick = [-10,0,10];
    sax(end).YTick = [-10,0,10];
    title(sax(end),{'Pause'});
    sax(end).YLabel.Units = 'centimeters';
    sax(end).YLabel.Position = [-0.55,0.74,0];
    daspect(sax(end),[1,1,1]);
    if phzInd == 2
        ylabel(sax(end),'cm');
        sax(end).YLabel.Units = 'centimeters';
        sax(end).YLabel.Position = [-0.55,0.74,0];
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
    [yind, yOffSet, xind, xOffSet] = deal(phzBin.count+1-phzInd+8, 0, 7, 1.2);
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
    for hvfInd = 1:hvfBin.count-1
        [ehpcmpKDE,dxi] = ksdensity( egoHvfPhz.control.meanPos( unitsEgoHvfCA1, phzInd, hvfInd+1, lat) );
        med             = median(    egoHvfPhz.control.meanPos( unitsEgoHvfCA1, phzInd, hvfInd+1, lat) );
        plot(dxi,ehpcmpKDE,'-','color',hvfBin.color(hvfInd+1,:));
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hvfBin.color(hvfInd+1,:));
    end
    % FORMAT subplot
    xlim(sax(end),[-10,10]);
    ylim(sax(end),[0,0.16]);
    grid(sax(end),'on');
    sax(end).YTick = [0,0.05,0.10,0.15];
    sax(end).YTickLabel = {'0','','0.1',''};
    sax(end).XTick = [-5,0,5];
    
    if phzInd == 2,
        ylabel(sax(end),'Prob');
        sax(end).YLabel.Units = 'centimeters';
        sax(end).YLabel.Position = [-0.55,0.74,0];
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
    [yind, yOffSet, xind, xOffSet] = deal(phzBin.count+1-phzInd+8, 0, 8, 1.2);
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
    for hvfInd = 1:hvfBin.count-1
        [ehpcmpKDE,dxi] = ksdensity( egoHvfPhz.control.meanPos( unitsEgoHvfCA1, phzInd, hvfInd+1, fwd) );
        med             = median(    egoHvfPhz.control.meanPos( unitsEgoHvfCA1, phzInd, hvfInd+1, fwd) );
        plot(dxi,ehpcmpKDE,'-','color',hvfBin.color(hvfInd+1,:))
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hvfBin.color(hvfInd+1,:));
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








configure_default_args();
EgoProCode2D_load_data();

% CA1
tind = [3,4,5,17,18,19,20,21,22,23,29];
%tind = [6,7,26,27,30];
sampleRate = 250;

global AP
% compute_ratemaps ---------------------------------------------------------------------------------
AP.compute_ratemaps =                                                                            ...
    struct('get_featureSet',            @fet_xy,                                                 ...
           'sampleRate',                16,                                                      ...
           'pfsArgs',                   struct('states',           'theta-groom-sit-rear',       ...
                                               'binDims',          [50,50],                      ...
                                               'SmoothingWeights', [2.4,2.4],                    ...
                                               'numIter',          1,                            ...
                                               'boundaryLimits',   [-500,500;-500,500],          ...
                                               'halfsample',       false)                        ...
           );
%---------------------------------------------------------------------------------------------------

dca = cf(@(T,U) accumulate_decoding_vars(T,U), Trials(tind),units(tind));


decoded = struct('fwd',[],...
                 'lat',[],...
                 'hvf',[],...
                 'hvl',[],...
                 'hav',[],...
                 'hba',[],...
                 'phz',[]);
for t = [1:3,5:8,11],
    %for t = [1:3,5:8,11],
    mind =    dca{t}.stcm(:,1)==1                                            ...
              & (dca{t}.stcm(:,3)==3|dca{t}.stcm(:,4)==4|dca{t}.stcm(:,5)==5)  ...
              & dca{t}.post>0.005 ...
         ...   & dca{t}.hvfl(:,1)>-2                                             ...
         ...  & abs(dca{t}.hvfl(:,2))>5                                        ...
            & dca{t}.ucnt>=2 & dca{t}.ucnt<8                                 ...
            & sqrt(sum(dca{t}.xyz(:,'hcom',[1,2]).^2,3))<325;
    decoded.fwd = cat(1,decoded.fwd,dca{t}.esax(mind,1));
    decoded.lat = cat(1,decoded.lat,dca{t}.esax(mind,2));%+20*double(t>4));
    decoded.hvf = cat(1,decoded.hvf,dca{t}.hvfl(mind,1));
    decoded.hvl = cat(1,decoded.hvl,dca{t}.hvfl(mind,2));
    decoded.hav = cat(1,decoded.hav,dca{t}.hvang(mind,1));
    decoded.hba = cat(1,decoded.hba,dca{t}.hbang(mind,1));
    decoded.phz = cat(1,decoded.phz,dca{t}.phz(mind,1));
end

ind = WithinRanges(decoded.phz,phzBin.edges(3:4)) ...
      & randn(size(decoded.hba))>0.5 ...
      & abs(decoded.hba)<1.2;
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind)]);
STATS


figure,
ind = WithinRanges(decoded.phz,phzBin.edges([3,4])) ...
      & randn(size(decoded.hba))>0              ...
 ;...     & abs(decoded.hvl)>5;
hist2([decoded.lat(ind),decoded.hba(ind)],linspace(-300,300,24),linspace(-1.2,1.2,24),'xprob')
colormap('jet');
Lines([],0,'w');
Lines(0,[],'w');

figure,
hist2([decoded.lat(ind),decoded.hvl(ind)],linspace(-300,300,24),linspace(-40,40,24),'xprob')
colormap('jet');
Lines([],0,'w');
Lines(0,[],'w');


mBinHvl.edges = linspace(-40,40,4);
mBinHvl.centers = mean([mBinHvl.edges(1:end-1);mBinHvl.edges(2:end)])
mBinHvl.count = numel(mBinHvl.edges)-1;
mBinHba.edges = linspace(-1.2,1.2,4);
mBinHba.centers = mean([mBinHba.edges(1:end-1);mBinHba.edges(2:end)])
mBinHba.count = numel(mBinHba.edges)-1;
mout = zeros([mBinHba.count,mBinHvl.count]);
for a = 1:mBinHba.count
    for v = 1:mBinHvl.count
        ind =   WithinRanges(decoded.phz,phzBin.edges([3,4])) ...
              & WithinRanges(decoded.hvl,mBinHvl.edges([v,v+1])) ...
              & WithinRanges(decoded.hba,mBinHba.edges([a,a+1])) ...              
              & randn(size(decoded.hba))>0;
        mout(a,v) = mean(decoded.lat(ind),'omitnan');
    end
end

figure,imagesc(mBinHba.centers,mBinHvl.centers,mout')
caxis([-60,60])
colormap('jet');



figure,
hist2([decoded.fwd,decoded.hba],linspace(-300,300,24),linspace(-1.2,1.2,24),'xprob')
colormap('jet');
Lines([],0,'w');
Lines(0,[],'w');

figure,
subplot(311);
ind = decoded.hba>0.2;
hist2([decoded.lat(ind),decoded.fwd(ind)],linspace(-400,400,32),linspace(-300,500,32));
caxis([0,4000])
Lines(0,[],'w');
Lines([],0,'w');
subplot(312);
ind = abs(decoded.hba)<0.2;
hist2([decoded.lat(ind),decoded.fwd(ind)],linspace(-400,400,32),linspace(-300,500,32));
caxis([0,4000])
Lines(0,[],'w');
Lines([],0,'w');
subplot(313);
ind = decoded.hba<-0.2;
hist2([decoded.lat(ind),decoded.fwd(ind)],linspace(-400,400,32),linspace(-300,500,32));
caxis([0,4000])
Lines(0,[],'w');
Lines([],0,'w');
colormap('jet');

binPhzs = linspace(0.5,2*pi-0.5,4);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
hbaBinEdg = [-1.2,-0.2,0.2,1.2];
hbaBinCtr = mean([hbaBinEdg(1:end-1);hbaBinEdg(2:end)]);
                          
hbaBin.edges = [-1.2,-0.2,0.2,1.2];
hbaBin.centers = mean([hbaBin.edges(1:end-1);hbaBin.edges(2:end)]);
hbaBin.count = numel(hbaBin.centers);        

phzBin.edges = linspace(0.5,2*pi-0.5,4);
phzBin.centers = (binPhzs(1:end-1)+binPhzs(2:end))./2;
phzBin.count = numel(phzBin.centers);


hbaBnds = {[0.2,1.2],[-0.2,0.2],[-1.2,-0.2]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
havBnds = {[-0.3,-0.018],[-0.018,-0.009],[-0.009,0.009],[0.009,0.018],[0.018,0.3]};
hvlBnds = {[-50,-5],[-5,5],[5,50]};
figure
norm = 'xprob';
xBinEdg = linspace(-300,300,8);
xBinCtr = mean([xBinEdg(1:end-1);xBinEdg(2:end)]);
yBinEdg = linspace(0.5,2*pi-0.5,4);
yBinCtr = mean([yBinEdg(1:end-1);yBinEdg(2:end)]);
clims = [0,0.4];
%clims = 'auto';
out = zeros([7,3,3,3]);
for h = 1:numel(hbaBnds)
    for v = 1:numel(hvlBnds)
        subplot2(numel(hvlBnds),numel(hbaBnds),v,h);
        ind = WithinRanges(decoded.hba,hbaBnds{h}) & ...
              WithinRanges(decoded.hvl,hvlBnds{v}) & ...
              randn(size(decoded.hvl))>0.5;
        out(:,:,h,v) = hist2([decoded.lat(ind), decoded.phz(ind)],xBinEdg,yBinEdg,norm);
        imagesc(xBinCtr, yBinCtr, out(:,:,h,v)');
        colormap('jet');caxis(clims);
        axis('xy');
    end
end


hbaBnds = {[0.2,1.2],[-0.2,0.2],[-1.2,-0.2]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
havBnds = {[-0.3,-0.018],[-0.018,-0.009],[-0.009,0.009],[0.009,0.018],[0.018,0.3]};
hvlBnds = {[-50,-5],[-5,5],[5,50]};
figure
norm = 'xprob';
xBinEdg = linspace(-300,300,8);
xBinCtr = mean([xBinEdg(1:end-1);xBinEdg(2:end)]);
yBinEdg = linspace(0.5,2*pi-0.5,4);
yBinCtr = mean([yBinEdg(1:end-1);yBinEdg(2:end)]);
clims = [0,0.4];
%clims = 'auto';
out = zeros([7,3,3,3]);
for h = 1:numel(hbaBnds)
    for v = 1:numel(hvlBnds)
        subplot2(numel(hvlBnds),numel(hbaBnds),v,h);
        ind = WithinRanges(decoded.hba,hbaBnds{h}) & ...
              WithinRanges(decoded.hvl,hvlBnds{v}) & ...
              randn(size(decoded.hvl))>0.5;
        out(:,:,h,v) = hist2([decoded.fwd(ind), decoded.phz(ind)],xBinEdg,yBinEdg,norm);
        imagesc(xBinCtr, yBinCtr, out(:,:,h,v)');
        colormap('jet');caxis(clims);
        axis('xy');
    end
end

%% CDF
hbaBnds = {[0.2,1.2],[-0.2,0.2],[-1.2,-0.2]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
havBnds = {[-0.3,-0.018],[-0.018,-0.009],[-0.009,0.009],[0.009,0.018],[0.018,0.3]};
hvlBnds = {[-50,-5],[-5,5],[5,50]};
figure
norm = 'xprob';
xBinEdg = linspace(-300,300,8);
xBinCtr = mean([xBinEdg(1:end-1);xBinEdg(2:end)]);
yBinEdg = linspace(0.5,2*pi-0.5,4);
yBinCtr = mean([yBinEdg(1:end-1);yBinEdg(2:end)]);
clims = [0,0.4];
%clims = 'auto';
out = zeros([7,3,3,3]);
for h = 1:numel(hbaBnds)
    for p = 1:3;    
    subplot2(phzBin.count,hbaBin.count,phzBin.count+1-p,h);
    hold('on');
    for v = 1:numel(hvlBnds)
        ind = WithinRanges(decoded.hba,hbaBnds{h}) & ...
              WithinRanges(decoded.hvl,hvlBnds{v}) & ...
              randn(size(decoded.hvl))>0.5;
        pind = WithinRanges(decoded.phz,phzBin.edges(p:p+1));
        cdfplot(decoded.lat(ind&pind));
    end
    xlim([-300,300]);
    end
end



figure,
clims = [0,0.05];
for h = 1:3
    for d = 1:2
        subplot2(2,3,d,h);
        imagesc(xBinCtr,yBinCtr,diff(out(:,:,h,[0:1]+d),[],4)');
        colormap('jet'); caxis(clims); axis('xy');
    end
end


figure,
clims = [-0.1,0.1];
for v = 1:3
    subplot2(3,2,v,1);
    imagesc(xBinCtr,yBinCtr,(out(:,:,1,v)-out(:,:,2,v))');
    colormap('jet'); caxis(clims); axis('xy');
    subplot2(3,2,v,2);
    imagesc(xBinCtr,yBinCtr,(out(:,:,3,v)-out(:,:,2,v))');
    colormap('jet'); caxis(clims); axis('xy');
end


figure,
clims = [-0.1,0.1];
for v = 1:3
    subplot2(2,3,1,v);
    imagesc(xBinCtr,yBinCtr,(out(:,:,v,1)-out(:,:,v,2))');
    colormap('jet'); caxis(clims); axis('xy');
    subplot2(2,3,2,v);
    imagesc(xBinCtr,yBinCtr,(out(:,:,v,3)-out(:,:,v,2))');
    colormap('jet'); caxis(clims); axis('xy');
end

figure
norm = 'xprob';
xBinEdg = linspace(-250,250,8);
yBinEdg = linspace(0.5,2*pi-0.5,4);
clims = [0,0.4];
for v = 1:3
    subplot2(1,3,1,v);
    ind = WithinRanges(decoded.hav,havBnds{v}) & ...
          randn(size(decoded.hba))>0.5;
    hist2([decoded.lat(ind), decoded.phz(ind)],xBinEdg,yBinEdg,norm);
    colormap('jet');caxis(clims);
end



figure
norm = 'xprob';
clims = [0.01,0.4];
sax = tight_subplot(3,1,[0.01,0.01],[0.1,0.1],[0.1,0.1]);
for v = 1:3
    ind = WithinRanges(decoded.hba,hbaBnds{v}) & ...
          randn(size(decoded.hba))>0.5;
    out = hist2([decoded.lat(ind)+25, decoded.phz(ind)],xBinEdg, ...
                yBinEdg,norm);
    for p = 1:3
        axes(sax(4-p))
        hold(sax(p),'on');
        plot(xBinCtr,cumsum(out(:,p))','-+','Color',bclr(v,:));
        Lines([],0.5,'k');
        Lines(0,[],'k');
    end
    ylim([0,1]);
    xlim([-250,250]);
    %colormap('jet');caxis(clims);
    %set(gca,'ColorScale','log');    
end




ind = WithinRanges(decoded.phz,binPhzs(3:4)) & randn(size(decoded.hba))>0.5;
ind = WithinRanges(decoded.phz,binPhzs(1:2)) & randn(size(decoded.hba))>0.5;
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind),decoded.hav(ind)]);
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind),decoded.hvl(ind)]);
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hav(ind)]);

ind = WithinRanges(decoded.phz,phzBin.edges(3:4)) ...
      & randn(size(decoded.hba))>0.5 ...
      & abs(decoded.hba)<1.2;
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind)]);
STATS

glm = fitglm([decoded.hba(ind),decoded.hav(ind)],decoded.lat(ind));
glm = fitglm([decoded.hba(ind),decoded.hvl(ind)],decoded.lat(ind));

% MTAData
% function summarize



%% HEAD BODY ANGLE ---------------------------------------------------------------------------------
clear('xcomp','ycomp','zcomp','ccomp');
xcomp.label = 'hba';            xcomp.data = [];    xcomp.edgs = [-1.25,-0.2,0.2,1.25];
ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5,6 );
ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-80,80];
fcomp.data = [];
mcomp.data = [];
%for t = 1:10
%for t = [2,3,4,5,7,8,9]
for t = [1:10]    
    dc = dca{t};
        
    mang = sq(dca{t}.xyz(:,'hcom',[1,2])-dca{t}.xyz(:,'bcom',[1,2]));
    mbang = atan2(mang(:,2),mang(:,1));
    bang =  sq(dca{t}.xyz(:,'bcom',[1,2]));
    bbang = atan2(bang(:,2),bang(:,1));
    mbbang = circ_dist(bbang,mbang);
    
    mind =  dc.stcm(:,1)==1                                       ... % theta state
           & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5)    ... % states {walk,turn,pause}
           & dc.hvfl(:,1)>2                                       ... % forward movement
           & dc.ucnt>=3 & dc.ucnt<9                               ... % number of coactive units
           & sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))<325;          ... % center of maze
            
    mind(mind==true) = randn([sum(mind),1])>0.5;                  ... % DROP half of the data points
                     
    xcomp.data = cat(1, xcomp.data, -dc.hbang(mind,1));
    ycomp.data = cat(1, ycomp.data, dc.phz(mind));
    %ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+0*double(t<5)-12.5*double(t>=5));    
    ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+25*double(t<4)-12.5*double(t>=4));
    fcomp.data = cat(1, fcomp.data, dc.esax(mind,1));
    mcomp.data = cat(1, mcomp.data, mbbang(mind));
end

[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,ccomp);
zmean   = zcomp.mean;    zmean(zcomp.count<5)  = nan;
zmedian = zcomp.median;  zmean(zcomp.count<5)  = nan;
zstd    = zcomp.std;     zstd(zcomp.count<5)   = nan;
zcount  = zcomp.count;   zcount(zcomp.count<5) = nan;
% $$$ set(figure(),'Units','centimeters','Position',[0,-3,8,28]);
% $$$ subplot(411); imagesc(xcomp.ctrs,ycomp.ctrs,zmean'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,['Mean ',ccomp.label]); colormap('jet');  caxis([ccomp.clim]);
% $$$     ylabel(ycomp.label);    
% $$$ subplot(412); imagesc(xcomp.ctrs,ycomp.ctrs,zstd'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,['Std ',ccomp.label]); 
% $$$     colormap('jet');
% $$$ subplot(413); 
% $$$ imagesc(xcomp.ctrs,ycomp.ctrs,zcount'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,'Count'); colormap('jet'); 
% $$$     xlabel(xcomp.label);
% $$$ subplot(414); 
% $$$     imagesc(xcomp.ctrs,ycomp.ctrs,bsxfun(@rdivide,zcount,sum(zcount,2,'omitnan'))'); axis('xy');    
% $$$     cax = colorbar(); ylabel(cax,'Probability'); colormap('jet'); 
% $$$     xlabel(xcomp.label);

xBlockOffset = 0;
yBlockOffset = 0;
figure
pinds = {ycomp.data > 0.5  & ycomp.data < 2.25,...
         ycomp.data > 2.25  & ycomp.data < 4.25,...
         ycomp.data > 4.25  & ycomp.data < 6};
phaseLabels = {'descending','trough','ascending'};
for p = 1:numel(pinds)
% $$$    [yind, yOffSet, xind, xOffSet] = deal(   nPhz+2-p,       0,    6,     fig.subplot.width/2+1.5);
% $$$     sax(end+1) = axes('Units','centimeters',                                        ...
% $$$                       'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
% $$$                         fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
% $$$                         fig.subplot.width,                                ...
% $$$                         fig.subplot.height],                              ...
% $$$                       'FontSize', 8,                                                ...
% $$$                       'LineWidth',1);
% $$$     hold(sax(end),'on');
    sax = subplot(3,1,4-p)
    hold('on');
    
 
    indL =    pinds{p} & xcomp.data > -1.2 & xcomp.data < -0.2;
    indC =    pinds{p} & xcomp.data > -0.2 & xcomp.data < 0.2;
    indR =    pinds{p} & xcomp.data >  0.2 & xcomp.data < 1.2;    
    cdfplot(ccomp.data(indL)./10);
    cdfplot(ccomp.data(indC)./10);
    cdfplot(ccomp.data(indR)./10);
    xlim([-25,25]);    
    title('');
    %ylim([-25,30]);
    sax(end).XTick = [-20,-10,0,10,20];
    title(phaseLabels{p});
    %sax(end).XTick = [-20,-10,0,10,20];
    if p~=1,
        %sax(end).XTickLabel = [];
        %sax(end).YTickLabel = [];
        sax.XTickLabel = [];
        %sax.YTickLabel = [];
        xlabel('')
        ylabel('')
    end
end
                 

