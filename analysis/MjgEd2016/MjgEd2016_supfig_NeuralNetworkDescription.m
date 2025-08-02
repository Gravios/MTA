
dbstop at 152 in bhv_nn.m

dbstop at 184 in bhv_nn.m

req20160128([],'compute');

set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)
      

% Figure preprocessing
hfig = figure(2016061401);clf
set(hfig,'PaperPositionMode','auto');
set(hfig,'units','centimeters');
set(hfig,'Position',[0,0,18,16])

% Individual output layer plots
w = 4;


hfig = figure(2016061401);
%xlm = [2950,3100];
xlm = round([51549, 55145]./Trial.xyz.sampleRate*12);
w = 4;
% NN example 1
ypos = 0;
axes('Units','centimeters',...
     'Position',[1,1+ypos,w,2]...
)
plot([1:length(d_state)]./12,d_state)
xlim([xlm(1),xlm(2)]./12)
ypos = 0;w = 4;
axes('Units','centimeters',...
     'Position',[1+5.5,1+ypos,w,2]...
)


plotSTC(Stc,1);
xlim([xlm(1),xlm(2)]./12)
ylim([0.5,7.5])
set(gca,'YTick',[1.5:6.5])
set(gca,'YTickLabels',{'walk','rear','turn','pause','groom','sit'})
dbcont

hfig = figure(2016061401);
%xlm = [2950,3100];
xlm = round([51549, 55145]./Trial.xyz.sampleRate*12);
% NN example 2
ypos = 3;w = 4;
axes('Units','centimeters',...
     'Position',[1,1+ypos,w,2]...
)
plot([1:length(d_state)]/12,d_state)
xlim([xlm(1),xlm(2)]./12)
axes('Units','centimeters',...
     'Position',[1+5.5,1+ypos,w,2]...
)

plotSTC(Stc,1);
xlim([xlm(1),xlm(2)]./12)
ylim([0.5,7.5])
set(gca,'YTick',[1.5:6.5])
set(gca,'YTickLabels',{'walk','rear','turn','pause','groom','sit'})
dbcont


hfig = figure(2016061401);
%xlm = [2950,3100];
xlm = round([51549, 55145]./Trial.xyz.sampleRate*12);
% NN example 3
ypos = 6;w = 4;
axes('Units','centimeters',...
     'Position',[1,1+ypos,w,2]...
)
plot([1:length(d_state)]./12,d_state)
xlim([xlm(1),xlm(2)]./12)
axes('Units','centimeters',...
     'Position',[1+5.5,1+ypos,w,2]...
)
plotSTC(Stc,1);
xlim([xlm(1),xlm(2)]./12)
ylim([0.5,7.5])
set(gca,'YTick',[1.5:6.5])
set(gca,'YTickLabels',{'walk','rear','turn','pause','groom','sit'})
dbcont

hfig = figure(2016061401);
xlm = [2950,3100];
xlm = round([51549, 55145]./Trial.xyz.sampleRate*12);
% NN example 4
ypos = 9;w = 4;
axes('Units','centimeters',...
     'Position',[1,1+ypos,w,2]...
)
plot([1:length(d_state)]./12,d_state)
xlim([xlm(1),xlm(2)]./12)
axes('Units','centimeters',...
     'Position',[1+5.5,1+ypos,w,2]...
)
plotSTC(Stc,1);
xlim([xlm(1),xlm(2)]./12)
ylim([0.5,7.5])
set(gca,'YTick',[1.5:6.5])
set(gca,'YTickLabels',{'walk','rear','turn','pause','groom','sit'})
dbcont


dbca

dbstop at 75 in req20160128.m

dbcont

hfig = figure(2016061401);
xlm = round([51549, 55145]./Trial.xyz.sampleRate*12);
% NN composite
ypos = 6;w = 6;
axes('Units','centimeters',...
     'Position',[11.5,1+ypos,w,4]...
)
imagesc([1:length(d_state{end})]./Trial.xyz.sampleRate,1:6,d_state{end}');
xlim([xlm(1),xlm(2)]./12)
axis xy
set(gca,'YTickLabels',{'walk','rear','turn','pause','groom','sit'})


ypos = 3;w = 6;
axes('Units','centimeters',...
     'Position',[11.5,1+ypos,w,1]...
)
plotSTC(stc{end},1,[],[],'brgymk',false);
xlim([xlm(1),xlm(2)]./12)


print(hfig,'-depsc2','-loose',...
      fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_2',...
               ['Supfig2_1_multi_patternnet.eps']))
