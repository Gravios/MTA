

configure_default_args();
EgoProCode2D_load_data();


% shuffling
% premutation test






global MTA_PROJECT_PATH
partsPath = fullfile(fullfile(MTA_PROJECT_PATH,'analysis','EgoProCode2D','EgoProCode2D_figure_parts'));
overwrite = false;
rat = load_patch_model('rat');




[hfig,fig,fax,sax] = set_figure_layout(figure(666003),'A4','portrait',[],2,2,0.2,0.4);


for region = 1:2
    if region == 1
% CA1
exampleUnit.trialIndex = 20;
exampleUnit.close.Xlims = [-200,400];
exampleUnit.close.Ylims = [-400,200];
exampleUnit.id = 25;
exampleUnit.maxRate = 18;
exampleUnit.index = find(units{exampleUnit.trialIndex}==exampleUnit.id);
exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];
uids = unitsEgoCA1;
globalXOffset = 0;
globalYOffset = 0;
else
% CA3
exampleUnit.trialIndex = 6;
exampleUnit.close.Xlims = [-200,400];
exampleUnit.close.Ylims = [-400,200];
exampleUnit.id = 10;
exampleUnit.maxRate = 18;
exampleUnit.index = find(units{exampleUnit.trialIndex}==exampleUnit.id);
exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];
uids = unitsEgoCA3;
globalXOffset = 0;
globalYOffset = -2.4*5;
end

[yind, yOffSet, xind, xOffSet] = deal(1, 0, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold(sax(end),'on');

plot(pft{exampleUnit.trialIndex},exampleUnit.id,1,'',[0,exampleUnit.maxRate],'colorMap',@jet);
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];



for d = 1:size(pfs{exampleUnit.trialIndex},2)
[yind, yOffSet, xind, xOffSet] = deal(1, 0, d+1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold(sax(end),'on');

% ADD rat
        subject = struct(rat);
        subject = update_subject_patch(subject,'head', d, true,hbaBinEdg,hbaBinCtr);
        subject = update_subject_patch(subject,'body',[],false,hbaBinEdg,hbaBinCtr);
        patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
        patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
        patch(subject.head.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
        line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
        line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
xlim(sax(end),[-150,150]);
ylim(sax(end),[-150,150]);
sax(end).XTick =[];
sax(end).YTick =[];
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
        
daspect(sax(end),[1,1,1]);
box(sax(end),'on');
 
end


for p = 1:size(pfs{exampleUnit.trialIndex},1)
        
%%%<<< PLOT 
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(p+1, 0, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold(sax(end),'on');

plot(pfet{exampleUnit.trialIndex}{4-p},exampleUnit.id,1,'',[0,exampleUnit.maxRate],'colorMap',@jet,'mazeMaskFlag',false,'flipAxesFlag',true);

% FORMAT subplot
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
xlim(sax(end),[-250,250]);
ylim(sax(end),[-250,250]);
daspect(sax(end),[1,1,1]);
box(sax(end),'on');

Lines([],0,'w');
Lines(0,[],'w');


%text(ets(exampleRange(1))+0.25,380,'Y');
axes(fax);
% $$$ line([sax(end).Position(1)-0.1].*[1,1],                                        ...
% $$$      sax(end).Position(2)+[0,sax(end).Position(4).*(200/diff(ylim(sax(end))))],...
% $$$      'Color','k',                                                              ...
% $$$      'LineWidth',1);

%%%>>>
end



for p = 1:size(pfs{exampleUnit.trialIndex},1)
for d = 1:size(pfs{exampleUnit.trialIndex},2)
        
%%%<<< PLOT 
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(p+1, 0, d+1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold(sax(end),'on');

plot(pfs{exampleUnit.trialIndex}{4-p,d},exampleUnit.id,1,'',[0,exampleUnit.maxRate],'colorMap',@jet,'mazeMaskFlag',false,'flipAxesFlag',true);

% FORMAT subplot
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
xlim(sax(end),[-250,250]);
ylim(sax(end),[-250,250]);
daspect(sax(end),[1,1,1]);
box(sax(end),'on');

Lines([],0,'w');
Lines(0,[],'w');

            subject = struct(rat);
            subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
            subject = update_subject_patch(subject,'body', numel(hbaBinCtr)+1-d,  true,hbaBinEdg,hbaBinCtr);
            patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
            patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
            patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
            line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
            line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);


%text(ets(exampleRange(1))+0.25,380,'Y');
axes(fax);
% $$$ line([sax(end).Position(1)-0.1].*[1,1],                                        ...
% $$$      sax(end).Position(2)+[0,sax(end).Position(4).*(200/diff(ylim(sax(end))))],...
% $$$      'Color','k',                                                              ...
% $$$      'LineWidth',1);

%%%>>>
end
end

for p = 1:size(pfs{exampleUnit.trialIndex},1)
%%%<<< PLOT 
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(4-p+1, 0, 6, -1.2);
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
plot(egoMeanRmapPosHba(uids,p,3,2) ,...
     egoMeanRmapPosHba(uids,p,1,2) ,...
     '.',...
     'MarkerFaceColor',pclr(p,:),...
     'MarkerEdgeColor',pclr(p,:));
% FORMAT subplot
grid(sax(end),'on');
xlim(sax(end),[-10,10]);
ylim(sax(end),[-10,10]);
sax(end).XTick = [-10,-5,0,5,10];
sax(end).YTick = [-5,0,5];
if p ~= 1,    sax(end).XTickLabel = {};                           end % not bottom of column
if p == 1,    xlabel(sax(end),'cm');                              end % bottom of column
if p == 3,    title(sax(end),{'Lateral Offset','Right VS Left'}); end % top of column
if p == 2,    ylabel(sax(end),'cm');                              end % middle of column
daspect(sax(end),[1,1,1]);
end




for p = 1:size(pfs{exampleUnit.trialIndex},1)
%%%<<< PLOT 
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(4-p+1, 0, 8, -2.4);
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
ylim(sax(end),[0,1]);
Lines([],0.5,'k');
Lines(0,[],'k');
[F,X] = ecdf((egoMeanRmapPosHba(uids,p,3,2)-egoMeanRmapPosHba(uids,p,2,2)));
plot(X,F,'Color','b');
[F,X] = ecdf((egoMeanRmapPosHba(uids,p,2,2)-egoMeanRmapPosHba(uids,p,1,2)));
plot(X,F,'Color','r');
[F,X] = ecdf((egoMeanRmapPosHba(uids,p,3,2)-egoMeanRmapPosHba(uids,p,1,2)));
plot(X,F,'Color',[219/255,172/255,52/255]);    
% FORMAT subplot
grid(sax(end),'on');
xlim(sax(end),[-10,10]);
ylim(sax(end),[0,1]);
sax(end).XTick = [-10,-5,0,5,10];
sax(end).YTick = [0,0.5,1];
if p ~= 1,    sax(end).XTickLabel = {};                           end % not bottom of column
if p == 1,    xlabel(sax(end),'cm');                              end % bottom of column
if p == 3,    title(sax(end),{'cdf','Lateral - Center'}); end % top of column
if p == 2,    ylabel(sax(end),'cm');                              end % middle of column
end

end



% $$$ 
% $$$ rmat = [cos(-pi/4),-sin(-pi/4);sin(-pi/4),cos(-pi/4)];
% $$$ rota = multiprod(rmat,[egoMeanRmapPos(uids,p,3,2) - egoMeanRmapPos(uids,p,2,2),egoMeanRmapPos(uids,p,2,2) - egoMeanRmapPos(uids,p,1,2)],[1,2],[2]);
% $$$ figure,
% $$$ subplot(211);    histogram(rota(:,1),linspace(-100,100,20));
% $$$ subplot(212);    histogram(rota(:,2),linspace(-100,100,20));
% $$$ 
% $$$ figure,
% $$$ subplot(211);    histogram(egoMeanRmapPos(uids,p,3,2) - egoMeanRmapPos(uids,p,2,2),linspace(-100,100,20));
% $$$ subplot(212);    histogram(egoMeanRmapPos(uids,p,1,2) - egoMeanRmapPos(uids,p,2,2),linspace(-100,100,20));
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,
% $$$ plot(egoMeanRmapPos(uids,p,3,2),...
% $$$      egoMeanRmapPos(uids,p,1,2),...
% $$$      '.')
% $$$ grid('on');
% $$$ 
% $$$ figure,
% $$$ plot(egoMeanRmapPos(uidsCA3,p,3,2)-egoMeanRmapPos(uidsCA3,p,2,2),...
% $$$      egoMeanRmapPos(uidsCA3,p,1,2)- egoMeanRmapPos(uidsCA3,p,2,2),...
% $$$      '.')
% $$$ grid('on');
% $$$ nn
% $$$ 
% $$$ figure,
% $$$ hold('on');
% $$$ plot(egoMeanRmapPos(uids,p,3,2),...
% $$$      egoMeanRmapPos(uids,p,3,1),...
% $$$      '.r')
% $$$ plot(egoMeanRmapPos(uids,p,1,2),...
% $$$      egoMeanRmapPos(uids,p,1,1),...
% $$$      '.b')
% $$$ xlim([-100,100]);
% $$$ ylim([-100,200]);
% $$$ 
% $$$ figure,
% $$$ hold('on');
% $$$ plot(egoMeanRmapPos(uids,p,3,2)-egoMeanRmapPos(uids,p,2,2),...
% $$$      egoMeanRmapPos(uids,p,3,1)-egoMeanRmapPos(uids,p,2,2),...
% $$$      '.r')
% $$$ plot(egoMeanRmapPos(uids,p,1,2)-egoMeanRmapPos(uids,p,2,2),...
% $$$      egoMeanRmapPos(uids,p,1,1)-egoMeanRmapPos(uids,p,2,2),...
% $$$      '.b')
% $$$ xlim([-150,150]);
% $$$ ylim([-100,200]);
% $$$ 
% $$$ 
% $$$ figure
% $$$ hold('on');
% $$$ histogram(egoMeanRmapPos(uids,p,1,2)-egoMeanRmapPos(uids,p,2,2),...
% $$$           linspace([-100,100,20]));
% $$$ histogram(egoMeanRmapPos(uids,p,3,2)-egoMeanRmapPos(uids,p,2,2),...
% $$$           linspace([-100,100,20]));
% $$$ 


figure,
for x = 1:3
    for y = 1:3
        subplot2(3,3,x,y);
        histogram((egoSizeHba(unitsEgoCA1,x,y)),linspace([0,800,20]));;
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