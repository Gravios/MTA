

configure_default_args();
MjgER2016_load_data();
unitsEgo = cf(@(T)  T.spk.get_unit_set(T,'egocentric'),  Trials); 
%EgoProCode2D_load_data();


% COMPUTE the size of primary egofields


%%%<<< Load general variables
sampleRate = 250;
%headCenterCorrection = [-25,-8];
pfsState = 'theta-groom-sit-rear';
hbaBinEdges = -1.5:0.6:1.5;
xyz = cf(@(t) preproc_xyz(t,'trb'),             Trials);
      cf(@(x) x.filter('ButFilter',3,30,'low'), xyz);    
      cf(@(x) x.resample(sampleRate),           xyz);
spk = cf(@(t,u) t.load('spk',sampleRate,'gper',u,'deburst'),Trials,units);    
pft = cf(@(t,u)  pfs_2d_theta(t,u),  Trials, units);
%%%>>>


overwrite = false;
%%%<<< (pfe) ego ratemap
pfe = cf(@(t,u,x,s,p)                                          ... Egocentric ratemap .
         compute_ego_ratemap(t,u,x,s,p,'overwrite',overwrite), ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft                                               ... MTAApfs object, theta state placefields 
);
%%%>>>

%%%<<< (pfet) egothp ratemap
pfet = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase.
         compute_egothp_ratemap(t,u,x,s,p,'overwrite',overwrite),...
             Trials,                        ... MTATrial
             units,                         ... Unit subset, placefields away from the maze walls
             xyz,                           ... MTADxyz object, head position
             spk,                           ... MTASpk object, spike time and id collection 
             pft                            ... MTAApfs object, theta state placefields 
);
%%%>>>


%%%<<< (pfs) egohba ratemap
pfs = cf(@(t,u,x,s,p)                       ... Egocentric ratemap | theta phase , head body angle.
compute_egohba_ratemap(t,u,x,s,p,'overwrite',overwrite),   ...
    Trials,                        ... MTATrial
units,                         ... Unit subset, placefields away from the maze walls
xyz,                           ... MTADxyz object, head position
spk,                           ... MTASpk object, spike time and id collection 
pft                            ... MTAApfs object, theta state placefields 
);
% $$$ pfsh = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase , shuffled head body angle.
% $$$          compute_egohba_ratemap_shuffled(t,u,x,s,p,'overwrite',overwrite),   ...
% $$$              Trials,                        ... MTATrial
% $$$              units,                         ... Unit subset, placefields away from the maze walls
% $$$              xyz,                           ... MTADxyz object, head position
% $$$              spk,                           ... MTASpk object, spike time and id collection 
% $$$              pft                            ... MTAApfs object, theta state placefields 
% $$$ );
%%%>>>

tind = [1:30];
%tind = [18];
%%%<<< (pfs) egohvf ratemap
pfv = cf(@(t,u,x,s,p)                       ... Egocentric ratemap | theta phase , head body angle.
         compute_egohvf_ratemap(t,u,x,s,p,'overwrite',overwrite),   ...
             Trials(tind),                        ... MTATrial
             units(tind),                         ... Unit subset, placefields away from the maze walls
             xyz(tind),                           ... MTADxyz object, head position
             spk(tind),                           ... MTASpk object, spike time and id collection 
             pft(tind)                            ... MTAApfs object, theta state placefields 
);

% $$$ pfsh = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase , shuffled head body angle.
% $$$          compute_egohba_ratemap_shuffled(t,u,x,s,p,'overwrite',overwrite),   ...
% $$$              Trials,                        ... MTATrial
% $$$              units,                         ... Unit subset, placefields away from the maze walls
% $$$              xyz,                           ... MTADxyz object, head position
% $$$              spk,                           ... MTASpk object, spike time and id collection 
% $$$              pft                            ... MTAApfs object, theta state placefields 
% $$$ );
%%%>>>


tind = [1:30];
%tind = [18];
%%%<<< (pfs) egohvf ratemap
pfl = cf(@(t,u,x,s,p)                       ... Egocentric ratemap | theta phase , head body angle.
         compute_egohvl_ratemap(t,u,x,s,p,'overwrite',overwrite),   ...
             Trials(tind),                        ... MTATrial
             units(tind),                         ... Unit subset, placefields away from the maze walls
             xyz(tind),                           ... MTADxyz object, head position
             spk(tind),                           ... MTASpk object, spike time and id collection 
             pft(tind)                            ... MTAApfs object, theta state placefields 
);

% $$$ pfsh = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase , shuffled head body angle.
% $$$          compute_egohba_ratemap_shuffled(t,u,x,s,p,'overwrite',overwrite),   ...
% $$$              Trials,                        ... MTATrial
% $$$              units,                         ... Unit subset, placefields away from the maze walls
% $$$              xyz,                           ... MTADxyz object, head position
% $$$              spk,                           ... MTASpk object, spike time and id collection 
% $$$              pft                            ... MTAApfs object, theta state placefields 
% $$$ );
%%%>>>


%%%<<< DIAGNOSTIC FIGURE
% $$$ t = 29;
% $$$ figure();
% $$$ for u = 1:numel(units{t});
% $$$     subplot2(5,2,1,1);
% $$$         plot(pft{t},units{t}(u),1,'colorbar',[],true);
% $$$         title(num2str(u));
% $$$     subplot2(5,2,2,1);
% $$$         plot(pfe{t},units{t}(u),1,'colorbar',[],false);
% $$$         mrate = max(cell2mat(cf(@(p) max(p.data.rateMap(:,p.data.clu==units{t}(u))), pfet{t})));
% $$$     for p = 1:5
% $$$         subplot2(5,2,p,2);
% $$$             plot(pfet{t}{p},units{t}(u),1,'colorbar',[0,mrate],false,'colorMap',@jet);
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ t = 18;
% $$$ figure();
% $$$ for u = 1:numel(units{t});
% $$$     subplot2(5,5,1,1);
% $$$         plot(pft{t},units{t}(u),1,'colorbar',[],true);
% $$$         title(num2str(units{t}(u)));
% $$$     subplot2(5,5,2,1);
% $$$         plot(pfe{t},units{t}(u),1,'colorbar',[],false,'flipAxesFlag',true);
% $$$         mrate = max(cell2mat(cf(@(p) max(p.data.rateMap(:,p.data.clu==units{t}(u))), pfet{t})));
% $$$         xlim([-300,300]);ylim([-300,300]);
% $$$     for p = 1:3
% $$$         subplot2(5,5,p,2);
% $$$             plot(pfet{t}{4-p},units{t}(u),1,'colorbar',[0,mrate],false,'flipAxesFlag',true,'colorMap',@jet);
% $$$             Lines(0,[],'m');
% $$$             Lines([],0,'m');
% $$$             xlim([-300,300]);ylim([-300,300]);            
% $$$     end
% $$$     for p = 1:3
% $$$         for h = 1:3
% $$$         subplot2(5,5,p,h+2);
% $$$             plot(pfs{t}{4-p,h},units{t}(u),1,'colorbar',[0,mrate*1.25],false,'flipAxesFlag',true,'colorMap',@jet);
% $$$             Lines(0,[],'m');
% $$$             Lines([],0,'m');
% $$$             xlim([-300,300]);ylim([-300,300]);
% $$$         end
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end
%%%>>>

%% MAIN FIGURE --------------------------------------------------------------------    
    
% LOAD spikes with behavioral variables
spkv = MjgER2016_load_spikeVars(Trials,units,sessionList,[],[],[],[],overwrite); 

% LOAD patch model 
rat = load_patch_model('rat');

pclr = cool(3);

[hfig,fig,fax,sax] = set_figure_layout(figure(666001),'A4','portrait',[],2,2,1,0.6);

yGlobalOffSet = -6;


%%%<<< First example cell ---------------------------------------------------------

xind = 1;
tid = 18;
uid = 11;
mrate = 18;
%%%<<< PLACEFIELD { x, y | (theta periods) } -------------------------------------------
% ADJUST subplot coordinates 
[yind, yOffSet, xind, xOffSet] = deal(1,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,        ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(pft{tid},uid,1,'text',[0,mrate],'colorMap',@jet);
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];

line(sax(end),...
     [-470,-270],...
     [460].*[1,1],'LineWidth',2,'Color','w');
line(sax(end),...
     [-460].*[1,1],...
     [470,270],'LineWidth',2,'Color','w');
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< SPIKES { thetaPhase(DRZ) | (theta periods) } ------------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,        ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
sind = ismember(spkv.map,[tid,uid],'rows')& spkv.stc(:,1)==1 & spkv.stc(:,2)==2;
out = hist2([spkv.hrz(sind,:),spkv.phz(sind,:)],linspace(-1,1,30),linspace(0,2*pi,30));
imagesc(linspace(-1,1,13),linspace(0,2*pi,13),imgaussfilt(out,2.2)');
colormap(sax(end),'jet')
sax(end).XTick =[-1,0,1];
sax(end).YTickLabel =[];
cax = colorbar(sax(end));
cax.Units = 'centimeters';
drawnow();
pause(0.5);
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
pause(0.5);
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
ylabel(cax,'count');
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< DIAGRAM Theta phase line reference ----------------------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2,0, xind, -0.75);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,        ...
                              fig.subplot.width./4,                             ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
ylim([0,1]);
plot(sax(end),-cos(linspace(0,2*pi)),linspace(0,1),'k','LineWidth',2);
sax(end).Visible = 'off';
text(0.5,0.05,'0','Rotation',90,'FontSize',8);
text(-0.5,0.5,'\pi','Rotation',90,'FontSize',8);
text(0.5,0.95,'2\pi','Rotation',90,'FontSize',8);
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< SPIKES { xEgo, yEgo | theta, Color: thetaPhase} ---------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(3,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,        ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
sind = ismember(spkv.map,[tid,uid],'rows')& spkv.stc(:,1)==1 & spkv.stc(:,2)==2 ;
scatter(spkv.ego(sind,2),spkv.ego(sind,1),2,spkv.phz(sind,:),'Filled');
colormap(sax(end),'hsv');
caxis(sax(end),[0,2*pi]);
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
cax = colorbar();
cax.Units = 'centimeters';
drawnow();
pause(0.5);
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
pause(0.5);
xlim([-300,300]);ylim([-250,350]);
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
ylabel(sax(end),'Forward');
xlabel(sax(end),'Lateral');
ylabel(cax,'Theta Phase');
line(sax(end),...
     [-270,-70],...
     [320].*[1,1],'LineWidth',2,'Color','k');
line(sax(end),...
     [-270].*[1,1],...
     [320,120],'LineWidth',2,'Color','k');
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< EGOFIELD { xEgo, yEgo | theta, Color: rate Hz} ----------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(4,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,        ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(pfe{tid},uid,1,'text',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);
axis(sax(end),'tight');
xlim([-300,300]);ylim([-250,350]);
Lines([],0,'w');
Lines(0,[],'w');
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
% ADD rat model
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false,[],[]);
subject = update_subject_patch(subject,'body',[],false,[],[]);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
line(sax(end),...
     [-270,-70],...
     [320].*[1,1],'LineWidth',2,'Color','w');
line(sax(end),...
     [-270].*[1,1],...
     [320,120],'LineWidth',2,'Color','w');
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< EGOFIELD { xEgo, yEgo | { θ | θ ∈ [4π/3, 2π] }, Color:Hz } ----------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(5,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,        ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(pfet{tid}{3},uid,1,'text',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);
axis(sax(end),'tight');
xlim([-300,300]);ylim([-250,350]);
Lines([],0,'w');
Lines(0,[],'w');
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
% ADD rat model
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false,[],[]);
subject = update_subject_patch(subject,'body',[],false,[],[]);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< EGOFIELD BY PHASE { xEgo, yEgo | { θ | θ ∈ [2π/3,4π/3] }, Color:Hz } ------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(6,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,        ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(pfet{tid}{2},uid,1,'text',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);
axis(sax(end),'tight');
xlim([-300,300]);ylim([-250,350]);
Lines([],0,'w');
Lines(0,[],'w');
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
% ADD rat model
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false,[],[]);
subject = update_subject_patch(subject,'body',[],false,[],[]);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< EGOFIELD BY PHASE { xEgo, yEgo | { θ | θ ∈ [0,2π/3] }, Color:Hz } ---------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(7,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,        ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(pfet{tid}{1},uid,1,'text',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);
axis(sax(end),'tight');
xlim([-300,300]);ylim([-250,350]);
Lines([],0,'w');
Lines(0,[],'w');
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
% ADD rat model
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false,[],[]);
subject = update_subject_patch(subject,'body',[],false,[],[]);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< DIAGRAM Theta phase line reference ----------------------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(7,0, xind, -1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,        ...
                              fig.subplot.width./3,                             ...
                              fig.subplot.height.*3+fig.subplot.verticalPadding*2],...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
ylim([0,1]);
plot(sax(end),-cos(linspace(0,2*pi)),linspace(0,1),'k','LineWidth',2);
plims = [0.5,              2.26106176905986];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(1,:), 'LineWidth',2);
plims = [2.26106176905986, 4.02212353811972];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(2,:), 'LineWidth',2);
plims = [4.02212353811972, 5.78318530717959];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(3,:),'LineWidth',2);
    
sax(end).Visible = 'off';
text(0.5,0.05,'0','Rotation',90,'FontSize',8);
text(-0.5,0.5,'\pi','Rotation',90,'FontSize',8);
text(0.5,0.95,'2\pi','Rotation',90,'FontSize',8);
%%%>>>
%---------------------------------------------------------------------------------------

% END first example cell
%%%>>> 


%%%<<< Second example cell
tid = 7; % Ed10
uid = 57;
xind = 4;
mrate = 26;
%%%<<< PLACEFIELD { x, y | (theta periods) } -------------------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(1,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,        ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(pft{tid},uid,1,'text',[0,mrate],'colorMap',@jet);
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
line(sax(end),...
     [-470,-270],...
     [460].*[1,1],'LineWidth',2,'Color','w');
line(sax(end),...
     [-460].*[1,1],...
     [470,270],'LineWidth',2,'Color','w');
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< SPIKES { thetaPhase(DRZ) | (theta periods) } ------------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,        ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
sind = ismember(spkv.map,[tid,uid],'rows')& spkv.stc(:,1)==1 & spkv.stc(:,2)==2;
out = hist2([spkv.hrz(sind,:),spkv.phz(sind,:)],linspace(-1,1,30),linspace(0,2*pi,30));
imagesc(linspace(-1,1,30),linspace(0,2*pi,30),imgaussfilt(out,2.2)');
cax = colorbar(sax(end));
cax.Units = 'centimeters';
drawnow();
pause(0.5);
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
pause(0.5);
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
ylabel(cax,'count');
axis(sax(end),'xy');
colormap(sax(end),'jet')
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< DIAGRAM Theta phase line reference ----------------------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2,0, xind, -0.75);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width./4,                             ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
ylim([0,1]);
plot(sax(end),-cos(linspace(0,2*pi)),linspace(0,1),'k','LineWidth',2);
sax(end).Visible = 'off';
text(0.5,0.05,'0','Rotation',90,'FontSize',8);
text(-0.5,0.5,'\pi','Rotation',90,'FontSize',8);
text(0.5,0.95,'2\pi','Rotation',90,'FontSize',8);
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< SPIKES { xEgo, yEgo | theta, Color: thetaPhase} ---------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(3,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
sind = ismember(spkv.map,[tid,uid],'rows')& spkv.stc(:,1)==1 & spkv.stc(:,2)==2 ;
scatter(spkv.ego(sind,2),spkv.ego(sind,1),2,spkv.phz(sind,:),'Filled');
colormap(sax(end),'hsv');
caxis(sax(end),[0,2*pi]);
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
cax = colorbar();
cax.Units = 'centimeters';
drawnow();
pause(0.5);
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
pause(0.5);
xlim([-300,300]);ylim([-250,350]);
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
ylabel(cax,'Theta Phase');
line(sax(end),...
     [-270,-70],...
     [320].*[1,1],'LineWidth',2,'Color','k');
line(sax(end),...
     [-270].*[1,1],...
     [320,120],'LineWidth',2,'Color','k');
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< EGOFIELD { xEgo, yEgo | theta, Color: rate Hz} ----------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(4,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(pfe{tid},uid,1,'text',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);
axis(sax(end),'tight');
xlim([-300,300]);ylim([-250,350]);
Lines([],0,'w');
Lines(0,[],'w');
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
% ADD rat model
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false,[],[]);
subject = update_subject_patch(subject,'body',[],false,[],[]);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
line(sax(end),...
     [-270,-70],...
     [320].*[1,1],'LineWidth',2,'Color','w');
line(sax(end),...
     [-270].*[1,1],...
     [320,120],'LineWidth',2,'Color','w');
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< EGOFIELD { xEgo, yEgo | { θ | θ ∈ [4π/3, 2π] }, Color:Hz } ----------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(5,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(pfet{tid}{3},uid,1,'text',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);
axis(sax(end),'tight');
xlim([-300,300]);ylim([-250,350]);
Lines([],0,'w');
Lines(0,[],'w');
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
% ADD rat model
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false,[],[]);
subject = update_subject_patch(subject,'body',[],false,[],[]);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< EGOFIELD BY PHASE { xEgo, yEgo | { θ | θ ∈ [2π/3,4π/3] }, Color:Hz } ------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(6,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(pfet{tid}{2},uid,1,'text',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);
axis(sax(end),'tight');
xlim([-300,300]);ylim([-250,350]);
Lines([],0,'w');
Lines(0,[],'w');
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
% ADD rat model
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false,[],[]);
subject = update_subject_patch(subject,'body',[],false,[],[]);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
%%%>>>
%---------------------------------------------------------------------------------------



%%%<<< EGOFIELD BY PHASE { xEgo, yEgo | { θ | θ ∈ [0,2π/3] }, Color:Hz } ---------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(7,0, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,        ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(pfet{tid}{1},uid,1,'text',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);
axis(sax(end),'tight');
xlim([-300,300]);ylim([-250,350]);
Lines([],0,'w');
Lines(0,[],'w');
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];
% ADD rat model
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false,[],[]);
subject = update_subject_patch(subject,'body',[],false,[],[]);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
%%%>>>
%---------------------------------------------------------------------------------------



%%%<<< DIAGRAM Theta phase line reference ----------------------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(7,0, xind, -1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                             ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                           ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,             ...
                              fig.subplot.width./3,                                  ...
                              fig.subplot.height.*3+fig.subplot.verticalPadding*2],  ...
                  'FontSize', 8,                                                     ...
                  'LineWidth',1);
hold(sax(end),'on');
ylim([0,1]);
plot(sax(end),-cos(linspace(0,2*pi)),linspace(0,1),'k','LineWidth',2);
plims = [0.5,              2.26106176905986];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(1,:), 'LineWidth',2);
plims = [2.26106176905986, 4.02212353811972];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(2,:), 'LineWidth',2);
plims = [4.02212353811972, 5.78318530717959];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(3,:), 'LineWidth',2);
sax(end).Visible = 'off';
text(0.5,0.05,'0','Rotation',90,'FontSize',8);
text(-0.5,0.5,'\pi','Rotation',90,'FontSize',8);
text(0.5,0.95,'2\pi','Rotation',90,'FontSize',8);
%%%>>>
%---------------------------------------------------------------------------------------


% END second example cell
%%%>>> 



%%%<<< COMPUTATION - theta phase resolved egocentric field center of mass --------------
ucounter = 1;
egoMeanRmapPos = [];
egoMaxRmapPos = [];
egoSize = [];
egoMeanRmapRate = [];
egoMaxRmapRate = [];
for t = 1:numel(Trials)
for u = 1:numel(unitsEgo{t})
unit = unitsEgo{t}(u);
for p = 1:3;
binSubsetX = abs(pfet{t}{p}.adata.bins{1})<300;
binSubsetY = abs(pfet{t}{p}.adata.bins{2})<300;
mapPosition = cell([1,2]);
[mapPosition{:}] = ndgrid(pfet{t}{p}.adata.bins{1}(binSubsetX),...
                          pfet{t}{p}.adata.bins{2}(binSubsetY));
mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
rmap = plot(pfet{t}{p},unit,[],[],[],false);
rmap = rmap(binSubsetX,binSubsetY);
nanmap = double(~isnan(rmap));
nanmap(nanmap==0) = nan;
rmap = rmap.*fliplr(nanmap);
rmap(isnan(rmap)) = 0;
rmap(rmap<2) = 0;
% $$$ subplot2(3,3,4-p,a);
% $$$ hold('on');
% $$$ imagescnan({pfet{t}{p}.adata.bins{1}(binSubsetX),...
% $$$                    pfet{t}{p}.adata.bins{2}(binSubsetY),...
% $$$                    rmap'});
nrmap =rmap./sum(rmap(:),'omitnan');
rmapCenter = sq(sum(sum(bsxfun(@times,nrmap,mapPosition),'omitnan'),'omitnan'))';
% $$$ plot(rmapCenter(1),rmapCenter(2),'*m')
% $$$ axis('tight')
% $$$ Lines([],0,'g');
% $$$ Lines(0,[],'g');
egoMeanRmapPos(ucounter,p,:) = rmapCenter./10;
egoSize(ucounter,p) = sum(nniz(nrmap(:)));
egoMeanRmapRate(ucounter,p,:) = mean(nonzeros(rmap));
[~,maxPos] = max(nrmap(:));
egoMaxRmapRate(ucounter,p) = rmap(maxPos);
if ~isempty(maxPos)
[maxX,maxY] = ind2sub(size(nrmap),maxPos);
egoMaxRmapPos(ucounter,p,:) = mapPosition(maxX,maxY,:);
else
egoMaxRmapPos(ucounter,p,:) = nan([1,1,1,2]);
end
end
ucounter = ucounter+1;
end
end

%%%>>>
%---------------------------------------------------------------------------------------


% $$$ uids = [1:164];
uidsCA1 = [1:19,44:149];

uidsCA3 = [20:41,123:127,150:164];

%%%<<< FORWARD egoMeanRmapPos CA1 ------------------------------------------------------
xind = 2;
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(5,0, xind, 1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
[h,L,MX,MED,bw] = violin({(egoMeanRmapPos(uidsCA1,1,1)-2),...
                           egoMeanRmapPos(uidsCA1,2,1)-2,...
                           egoMeanRmapPos(uidsCA1,3,1)-2});
delete(L);
af(@(hndl) set(hndl,'Vertices',fliplr(hndl.Vertices)), h);
af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr,ones([3,1]))');
xlim([-15,15])
ylim([0.5,3.5])
hlns = findobj(sax(end),'Type','line');
for ln = 1:numel(hlns)
    [hlns(ln).XData,hlns(ln).YData] = deal(hlns(ln).YData,hlns(ln).XData);
end
sax(end).XTickLabelMode = 'auto';
sax(end).XTick = [-15,-7.5,0,7.5,15];
sax(end).YTick =[];
grid(sax(end),'on');
title(sax(end),{'Mean Field Position','Forward'})
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< DIAGRAM Theta phase line reference ----------------------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(5,0, xind, 0.25);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width./4,                             ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
ylim([0,1]);
plot(sax(end),-cos(linspace(0,2*pi)),linspace(0,1),'k','LineWidth',2);
plims = [0.5,              2.26106176905986];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(1,:), 'LineWidth',2);
plims = [2.26106176905986, 4.02212353811972];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(2,:), 'LineWidth',2);
plims = [4.02212353811972, 5.78318530717959];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(3,:), 'LineWidth',2);
    
sax(end).Visible = 'off';
text(0.5,0.05,'0','Rotation',90,'FontSize',8);
text(-0.5,0.5,'\pi','Rotation',90,'FontSize',8);
text(0.5,0.95,'2\pi','Rotation',90,'FontSize',8);
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< LATERAL egoMeanRmapPos ----------------------------------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(6,-1, xind, 1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
[h,L,MX,MED,bw] = violin({egoMeanRmapPos(uidsCA1,1,2)+0.8,...
                          egoMeanRmapPos(uidsCA1,2,2)+0.8,...
                          egoMeanRmapPos(uidsCA1,3,2)+0.8})
delete(L);
af(@(hndl) set(hndl,'Vertices',fliplr(hndl.Vertices)), h);
af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr,ones([3,1]))');
xlim([-15,15])
ylim([0.5,3.5])
hlns = findobj(sax(end),'Type','line');
for ln = 1:numel(hlns)
    [hlns(ln).XData,hlns(ln).YData] = deal(hlns(ln).YData,hlns(ln).XData);
end
sax(end).XTickLabelMode = 'auto';
sax(end).XTick = [-15,-7.5,0,7.5,15];
sax(end).YTick =[];
grid(sax(end),'on');
title(sax(end),{'Lateral'})
xlabel(sax(end),'cm')
%%%>>>
%---------------------------------------------------------------------------------------


%%%<<< DIAGRAM Theta phase line reference ----------------------------------------------
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(6,-1, xind, 0.25);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width./4,                             ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
ylim([0,1]);
plot(sax(end),-cos(linspace(0,2*pi)),linspace(0,1),'k','LineWidth',2);
plims = [0.5,              2.26106176905986];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(1,:), 'LineWidth',2);
plims = [2.26106176905986, 4.02212353811972];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(2,:), 'LineWidth',2);
plims = [4.02212353811972, 5.78318530717959];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(3,:), 'LineWidth',2);
    
sax(end).Visible = 'off';
text(0.5,0.05,'0','Rotation',90,'FontSize',8);
text(-0.5,0.5,'\pi','Rotation',90,'FontSize',8);
text(0.5,0.95,'2\pi','Rotation',90,'FontSize',8);
%%%>>>
%---------------------------------------------------------------------------------------




%%%<<< proprotional rates between ascending and descending phase relative to the trough egoMeanRmapPos
[yind, yOffSet, xind, xOffSet] = deal(1,0, xind, 1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(egoMaxRmapRate(uidsCA1,1)./egoMaxRmapRate(uidsCA1,2),...
     egoMaxRmapRate(uidsCA1,3)./egoMaxRmapRate(uidsCA1,2),'.')
xlim([0,1.8]);ylim([0,1.8]);
grid(sax(end),'on');
title(sax(end),'Rate Ratio');
xlabel(sax(end),'Dsc vs Trough');
ylabel(sax(end),'Asc vs Trough');


[yind, yOffSet, xind, xOffSet] = deal(2,-1, xind, 1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot((egoSize(uidsCA1,2))*0.02^2,(egoSize(uidsCA1,1,1))*0.02^2,'.','Color',pclr(1,:));
plot((egoSize(uidsCA1,2))*0.02^2,(egoSize(uidsCA1,3,1))*0.02^2,'.','Color',pclr(3,:));
line([0,0.3],[0,0.3],'Color','k')
grid(sax(end),'on');
title(sax(end),'Field Area')
xlabel(sax(end),'m^2');
ylabel(sax(end),'m^2');


%% CA3
xind = 5;
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(5,0, xind, 1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
[h,L,MX,MED,bw] = violin({egoMeanRmapPos(uidsCA3,1,1)-2,...
                          egoMeanRmapPos(uidsCA3,2,1)-2,...
                          egoMeanRmapPos(uidsCA3,3,1)-2})
delete(L);
af(@(hndl) set(hndl,'Vertices',fliplr(hndl.Vertices)), h);
af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr,ones([3,1]))');
xlim([-15,15])
ylim([0.5,3.5])
hlns = findobj(sax(end),'Type','line');
for ln = 1:numel(hlns)
    [hlns(ln).XData,hlns(ln).YData] = deal(hlns(ln).YData,hlns(ln).XData);
end
sax(end).XTickLabelMode = 'auto';
sax(end).XTick = [-15,-7.5,0,7.5,15];
sax(end).YTick =[];
grid(sax(end),'on');
title(sax(end),{'Mean Field Position','Forward'})

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(5,0, xind, 0.25);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width./4,                             ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
ylim([0,1]);
plot(sax(end),-cos(linspace(0,2*pi)),linspace(0,1),'k','LineWidth',2);
plims = [0.5,              2.26106176905986];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(1,:), 'LineWidth',2);
plims = [2.26106176905986, 4.02212353811972];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(2,:), 'LineWidth',2);
plims = [4.02212353811972, 5.78318530717959];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(3,:),'LineWidth',2);
    
sax(end).Visible = 'off';
text(0.5,0.05,'0','Rotation',90,'FontSize',8);
text(-0.5,0.5,'\pi','Rotation',90,'FontSize',8);
text(0.5,0.95,'2\pi','Rotation',90,'FontSize',8);



% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(6,-1, xind, 1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
[h,L,MX,MED,bw] = violin({egoMeanRmapPos(uidsCA3,1,2)+0.8,...
                          egoMeanRmapPos(uidsCA3,2,2)+0.8,...
                          egoMeanRmapPos(uidsCA3,3,2)+0.8})
delete(L);
af(@(hndl) set(hndl,'Vertices',fliplr(hndl.Vertices)), h);
af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr,ones([3,1]))');
xlim([-15,15])
ylim([0.5,3.5])
hlns = findobj(sax(end),'Type','line');
for ln = 1:numel(hlns)
    [hlns(ln).XData,hlns(ln).YData] = deal(hlns(ln).YData,hlns(ln).XData);
end
sax(end).XTickLabelMode = 'auto';
sax(end).XTick = [-15,-7.5,0,7.5,15];
sax(end).YTick =[];
grid(sax(end),'on');
title(sax(end),{'Lateral'})
xlabel(sax(end),'cm')

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(6,-1, xind, 0.25);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width./4,                             ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
ylim([0,1]);
plot(sax(end),-cos(linspace(0,2*pi)),linspace(0,1),'k','LineWidth',2);
plims = [0.5, 2.26106176905986];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)),...
         'Color',pclr(1,:),...
         'LineWidth',2);
plims = [2.26106176905986, 4.02212353811972];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)),...
         'Color',pclr(2,:),...
         'LineWidth',2);
plims = [4.02212353811972, 5.78318530717959];
    plot(sax(end),...
         -cos(linspace(plims(1),plims(2))),...
         linspace(plims(1)/(2*pi),plims(2)/(2*pi)),...
         'Color',pclr(3,:),'LineWidth',2);
sax(end).Visible = 'off';
text(0.5,0.05,'0','Rotation',90,'FontSize',8);
text(-0.5,0.5,'\pi','Rotation',90,'FontSize',8);
text(0.5,0.95,'2\pi','Rotation',90,'FontSize',8);



[yind, yOffSet, xind, xOffSet] = deal(1,0, xind, 1);
% CREATE subplot axes
figure, sax(end) = axes()
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
% $$$ plot(egoMaxRmapRate(uidsCA3,1)./egoMaxRmapRate(uidsCA3,2),...
% $$$      egoMaxRmapRate(uidsCA3,3)./egoMaxRmapRate(uidsCA3,2),'.')
figure
scatter(egoMaxRmapRate(uidsCA3,1)./egoMaxRmapRate(uidsCA3,2),...
        egoMaxRmapRate(uidsCA3,3)./egoMaxRmapRate(uidsCA3,2),...
        10,...
        log10(max(egoMaxRmapRate(uidsCA3,:),[],2)),...
        'filled' ...
);
colormap('jet');

figure
scatter(egoMaxRmapRate(uidsCA1,1)./egoMaxRmapRate(uidsCA1,2),...
        egoMaxRmapRate(uidsCA1,3)./egoMaxRmapRate(uidsCA1,2),...
        15,...
        log10(max(egoMaxRmapRate(uidsCA1,:),[],2)),...
        'filled' ...
);
colormap('jet');

xlim([0,1.8]);ylim([0,1.8]);
grid(sax(end),'on');
title(sax(end),'Rate Ratio');
xlabel(sax(end),'Dsc vs Trough');
ylabel(sax(end),'Asc vs Trough');


[yind, yOffSet, xind, xOffSet] = deal(2,-1, xind, 1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet+yGlobalOffSet,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot((egoSize(uidsCA3,2))*0.02^2,(egoSize(uidsCA3,1,1))*0.02^2,'.','Color',pclr(1,:));
plot((egoSize(uidsCA3,2))*0.02^2,(egoSize(uidsCA3,3,1))*0.02^2,'.','Color',pclr(3,:));
line([0,0.3],[0,0.3],'Color','k')
grid(sax(end),'on');
title(sax(end),'Field Area')
xlabel(sax(end),'m^2');
ylabel(sax(end),'m^2');

text(fax,...
     fig.page.xpos(4)+2.75, ...
     fig.page.ypos(1)+3+yGlobalOffSet,...
     'CA3');

text(fax,...
     fig.page.xpos(1)+2.75, ...
     fig.page.ypos(1)+3+yGlobalOffSet,...
     'CA1');

xind = 1
[yind, yOffSet, xind, xOffSet] = deal(1,-2, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*5,                                ...
                              fig.subplot.height*3],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');

% $$$ xyz{tid}.model = Trials{tid}.model;
% $$$ ang = create(MTADang,Trials{tid},xyz{tid});
% $$$ for c = 1:numel(xyz{tid}.model.Connections)
% $$$     xyz{tid}.model.Connections{c}.color = [0.25,0.25,0.25];
% $$$ end
% $$$ lfp = Trials{tid}.load('lfp',Trials{tid}.meta.channelGroup.theta);
% $$$ %lfp = Trials{tid}.load('lfp',65);
% $$$ lfp.resample(sampleRate);
% $$$ lfp.data = nunity(lfp.data);
% $$$ phz = load_theta_phase(Trials{tid},sampleRate);
%timeIndex = 2300;
timeIndex = 87000;

tid = 18;
uid = 11;
mrate = 18;
rmap = plot(pft{tid},uid,1,'text',[0,mrate],'colorMap',@jet);
pos=cell([1,2]);
[pos{:}] = ndgrid(pft{tid}.adata.bins{:});
sux = surf(pos{:},rmap*6,rmap);
plotSkeleton(Trials{tid},xyz{tid},timeIndex,'surface',ang,[],[-200,200],{'head_front'});
xlim([50,360]);
ylim([-300,50]);
% $$$ xlim([0,450]);
% $$$ ylim([-300,50]);
sux.CDataMode = 'auto';
colormap(sax(end),'jet')
caxis([0,65])
scatter3(xyz{tid}(mres,'head_front',1),xyz{tid}(mres,'head_front',2),xyz{tid}(mres,'head_front',3),25,'m','filled');
mres = spk{tid}(uid); 
mres  = mres(WithinRanges(mres,timeIndex+[-200,200]));
sax(end).Visible = 'off';

view([-174.260519685344 3.92967460420232]);


xind = 1
[yind, yOffSet, xind, xOffSet] = deal(2,0.5, xind, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*5,                                ...
                              fig.subplot.height*0.5],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(([-200:200])./sampleRate,lfp(timeIndex+[-200:200],1)+3,'k','LineWidth',1)
plot((mres-timeIndex)./sampleRate,lfp(mres,1)+3,'*m');
plot(([-200:200])./sampleRate,phz(timeIndex+[-200:200],1),'g','LineWidth',1)
sax(end).Visible = 'off';




figure,
plot(egoMeanRmapPos(uidsCA1,1,1)-egoMeanRmapPos(uidsCA1,2,1),...
        egoMeanRmapPos(uidsCA1,3,1)-egoMeanRmapPos(uidsCA1,2,1),...
     '.');
figure,
plot(egoMeanRmapPos(uidsCA3,1,1)-egoMeanRmapPos(uidsCA3,2,1),...
        egoMeanRmapPos(uidsCA3,3,1)-egoMeanRmapPos(uidsCA3,2,1),...
     '.');

    
figure();
subplot(311);
    plot(egoMaxRmapRate(uidsCA1,2),egoMaxRmapRate(uidsCA1,1),'.');
    line([0,12.5],[0,25],'Color','r');
    line([0,25],[0,12.5],'Color','r');
    line([0,25],[0,25])    ;
    xlim([0,25]);    
    ylim([0,25]);   
    daspect([1,1,1]);    
    set(gca(),'XTick',[0:5:25]);
    set(gca(),'YTick',[0:5:25]);    
    grid('on');
subplot(312);
    plot(egoMaxRmapRate(uidsCA1,2),egoMaxRmapRate(uidsCA1,3),'.');
    line([0,12.5],[0,25],'Color','r');
    line([0,25],[0,12.5],'Color','r');
    line([0,25],[0,25]);
    xlim([0,25]);    
    ylim([0,25]);    
    daspect([1,1,1]);
    line([0,25],[0,25]);
    set(gca(),'XTick',[0:5:25]);
    set(gca(),'YTick',[0:5:25]);    
    grid('on');    
subplot(313);    
    plot(egoMaxRmapRate(uidsCA1,1),egoMaxRmapRate(uidsCA1,3),'.');
    line([0,12.5],[0,25],'Color','r');
    line([0,25],[0,12.5],'Color','r');
    line([0,25],[0,25]);
    xlim([0,25]);    
    ylim([0,25]);   
    daspect([1,1,1]);   
    set(gca(),'XTick',[0:5:25]);
    set(gca(),'YTick',[0:5:25]);    
    grid('on');
    

figure();
subplot(311);
    plot(egoMeanRmapRate(uidsCA1,2),egoMeanRmapRate(uidsCA1,1),'.');
    line([0,5],[0,10],'Color','r');
    line([0,10],[0,5],'Color','r');
    line([0,10],[0,10])    ;
    xlim([0,10]);    
    ylim([0,10]);   
    daspect([1,1,1]);    
    set(gca(),'XTick',[0:2:10]);
    set(gca(),'YTick',[0:2:10]);    
    grid('on');
subplot(312);
    plot(egoMeanRmapRate(uidsCA1,2),egoMeanRmapRate(uidsCA1,3),'.');
    line([0,5],[0,10],'Color','r');
    line([0,10],[0,5],'Color','r');
    line([0,10],[0,10]);
    xlim([0,10]);    
    ylim([0,10]);    
    daspect([1,1,1]);
    line([0,10],[0,10]);
    set(gca(),'XTick',[0:2:10]);
    set(gca(),'YTick',[0:2:10]);    
    grid('on');    
subplot(313);    
    plot(egoMeanRmapRate(uidsCA1,1),egoMeanRmapRate(uidsCA1,3),'.');
    line([0,5],[0,10],'Color','r');
    line([0,10],[0,5],'Color','r');
    line([0,10],[0,10]);
    xlim([0,10]);    
    ylim([0,10]);   
    daspect([1,1,1]);   
    set(gca(),'XTick',[0:2:10]);
    set(gca(),'YTick',[0:2:10]);    
    grid('on');


% $$$ figure,
% $$$ plot(egoMeanRmapPos(uidsCA1,2,1),...
% $$$     (egoMeanRmapPos(uidsCA1,3,1)+egoMeanRmapPos(uidsCA1,1,1))./2,...
% $$$             '.');
% $$$ line([-20,100],[-20,100]);
% $$$ figure();
% $$$ plot(sqrt(egoSize(uidsCA3,:)'*0.02^2))
% $$$ figure();
% $$$ plot(egoSize(uidsCA3,:)'*0.02^2)
% $$$ line([0,0.3],[0,0.3],'Color','k')




%%%<<< EXAMPLE Allo to ego coordinates


t = 20;
u = 2;
ny = 5;
nx = 3;

hvec = xyz{t}(:,'nose',[1,2])-xyz{t}(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(Trials{t}.meta.correction.headYaw),-sin(Trials{t}.meta.correction.headYaw); ...
                  sin(Trials{t}.meta.correction.headYaw),cos(Trials{t}.meta.correction.headYaw)],...
                 [2,3],...
                 [1,2]);


[mxr,mxp] = pft{t}.maxRate(units{t}(u));        
pfstrj = MTADfet.encapsulate( ...
        Trials{t},...
        multiprod(bsxfun(@minus,...
                         mxp,...
                         sq(xyz{t}(:,'hcom',[1,2]))),...
                  hvec,2,[2,3]),...
        sampleRate,...
        'pfstrj','ptrj','t');

figure
clf();
set(gcf(),'PaperType','A3');
set(gcf(),'PaperOrientation','landscape');
tspn = 3.3e4:4e4;
subplot(121);
    hold('on');
    plot(pft{t},units{t}(u),'colorMap',@bone);
    scatter3(xyz{t}(tspn,'hcom',1),...
             xyz{t}(tspn,'hcom',2),300.*ones(size(tspn)),10,atan2(hvec(tspn,1,2),hvec(tspn,1,1)),'filled')
    plot(mxp(1),mxp(2),'or');
    %plot(xyz{t}(tspn(1),'spine_middle',1),xyz{t}(tspn(1),'spine_middle',2),'^g');    
    %plot(xyz{t}(tspn(end),'spine_middle',1),xyz{t}(tspn(end),'spine_middle',2),'^r');        
    plotSkeleton(Trials{t},xyz{t},tspn(1));        
    plotSkeleton(Trials{t},xyz{t},tspn(round(length(tspn)/2.25)));    
    xlim([-325,350]); % mannually set
    ylim([-450,0]);   % mannually set
    colormap('hsv');    
    caxis([-pi,pi]);        
    title('Allocentric Coordinates')    
    daspect([1,1,1]);    
    grid('on');    
subplot(122);
    hold('on');
    scatter(pfstrj(tspn,2),pfstrj(tspn,1),10,atan2(hvec(tspn,1,2),hvec(tspn,1,1)),'filled')
    %plot(pfstrj(tspn(1),2),pfstrj(tspn(1),1),'^g','MarkerSize',10,'MarkerFaceColor','g')
    %plot(pfstrj(tspn(end),2),pfstrj(tspn(end),1),'^r','MarkerSize',10,'MarkerFaceColor','r')
    %legend({'Trajectory','Start','Stop'});
    xlim([-500,500]);    
    ylim([-500,500]);
    title('Egocentric Coordinates')
    colormap('hsv');
    caxis([-pi,pi]);            
    daspect([1,1,1]);
    colorbar();
    grid('on');


%%%>>>