% EgoProCode2D_f1_parts()

close_figure = true;
close_figure = false;


global MTA_PATH
global MTA_PROJECT_PATH
partsPath = fullfile(fullfile(MTA_PROJECT_PATH,'analysis','EgoProCode2D','EgoProCode2D_figure_parts'));
create_directory(partsPath);


EgoProCode2D_f1_data();




% $$$     tid=21;;
% $$$     uid = 27;
% $$$     u=5;
% $$$     mrate = 10; % Max Rate    
    
% GENERATE subplots    

EgoProCode2D_f1_subplot_allo_to_ego_schematic                     ();

% $$$     % ALLO stuff
EgoProCode2D_f1_subplot_placefieldExample                         ();
EgoProCode2D_f1_subplot_trajectories_allo                         ();
EgoProCode2D_f1_subplot_trajectory_allo_head_direction            ();  

% $$$     % EGO stuff    
EgoProCode2D_f1_subplot_trajectory_ego_head_direction             ();
EgoProCode2D_f1_subplot_trajectories_ego                          ();
EgoProCode2D_f1_subplot_egofieldExample                           ();

EgoProCode2D_f1_subplot_theta_cycle_horizontal                    ();

% EgoProCode2D_f1_subplot_theta_cycle_vertical                      ();

% EgoProCode2D_f1_subplot_traj_allo_theta_descend                   ();
% EgoProCode2D_f1_subplot_traj_allo_theta_trough                    ();
% EgoProCode2D_f1_subplot_traj_allo_theta_ascend                    ();

% EgoProCode2D_f1_subplot_traj_ego_theta_descend                    ();
% EgoProCode2D_f1_subplot_traj_ego_theta_trough                     ();
% EgoProCode2D_f1_subplot_traj_ego_theta_ascend                     ();

EgoProCode2D_f1_subplot_allofieldExample_descending_phase         ();
EgoProCode2D_f1_subplot_allofieldExample_trough_phase             ();
EgoProCode2D_f1_subplot_allofieldExample_ascending_phase          ();

EgoProCode2D_f1_subplot_egofieldExample_descending_phase          ();
EgoProCode2D_f1_subplot_egofieldExample_trough_phase              ();
EgoProCode2D_f1_subplot_egofieldExample_ascending_phase           ();

EgoProCode2D_f1_subplot_forward_field_distrib_every_phase         ();
%EgoProCode2D_f1_subplot_forward_field_distrib_descending_phase    ();
%EgoProCode2D_f1_subplot_forward_field_distrib_trough_phase        ();
%EgoProCode2D_f1_subplot_forward_field_distrib_ascending_phase     ();

EgoProCode2D_f1_subplot_lateral_field_distrib_every_phase         ();
% EgoProCode2D_f1_subplot_lateral_field_distrib_descending_phase    ();
% EgoProCode2D_f1_subplot_lateral_field_distrib_trough_phase        ();    
% EgoProCode2D_f1_subplot_lateral_field_distrib_ascending_phase     ();   

EgoProCode2D_f1_subplot_forward_distrib_stats                     ();
EgoProCode2D_f1_subplot_lateral_distrib_stats                     ();


% $$$ EgoProCode2D_f1_descending_label_ca1                              ();   
% $$$ EgoProCode2D_f1_trough_label_ca1                                  ();   
% $$$ EgoProCode2D_f1_ascending_label_ca1                               ();   

% $$$ EgoProCode2D_f1_descending_label_ca3                              ();   
% $$$ EgoProCode2D_f1_trough_label_ca3                                  ();   
% $$$ EgoProCode2D_f1_ascending_label_ca3                               ();   

% EgoProCode2D_f1_subplot_ca1_rate_ratios                           ();
% EgoProCode2D_f1_subplot_ca1_field_size                            ();
% EgoProCode2D_f1_subplot_ca3_rate_ratios                           ();
% EgoProCode2D_f1_subplot_ca3_field_size                            ();

% EgoProCode2D_f1_subplot_lfp_timeseries                            ();

EgoProCode2D_f1_ego_maxrate_thetaphase_ca1                        ();
EgoProCode2D_f1_ego_maxrate_thetaphase_ca3                        ();

EgoProCode2D_f1_ego_maxrate_stats_ca1                             ();
EgoProCode2D_f1_ego_maxrate_stats_ca3                             ();
EgoProCode2D_f1_ego_field_size_stats_ca1                          ();
EgoProCode2D_f1_ego_field_size_stats_ca3                          ();


%%%<<< PLACEFIELD { x, y | (theta periods) } -------------------------------------------

%function EgoProCode2D_f1_subplot_placefieldExample()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_placefieldExample');
% $$$ 
% $$$     % SUBPLOT <- placefield, scalebars, circle
% $$$     plot(pft{exampleUnit.trialIndex},exampleUnit.id,1,'',[0,mrate],'colorMap',@jet);
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     line(sax,                               ...
% $$$          [-470,-270],                       ...
% $$$          [460].*[1,1],'LineWidth',2,'Color','w');
% $$$     line(sax,                            ...
% $$$          [-460].*[1,1],                  ...
% $$$          [470,270],'LineWidth',2,'Color','w');
% $$$     circle(mxp(1),mxp(2),200,'-r');        
% $$$     
% $$$     savefig(pfig, fullfile(partsPath,[sax.Tag,'.fig']));
% $$$     close(pfig);
%end

%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< TRAJECTORIES { x, y | (theta periods) } -----------------------------------------

%function EgoProCode2D_f1_subplot_trajectories_allo()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_trajectories_allo');
% $$$     
% $$$     % PLOT trajectories
% $$$     plot(sax, ...
% $$$          xyz{exampleUnit.trialIndex}(sper,'hcom',1),...
% $$$          xyz{exampleUnit.trialIndex}(sper,'hcom',2),...
% $$$          '.',...
% $$$          'MarkerFaceColor',[0.25,0.25,0.25],...
% $$$          'MarkerEdgeColor',[0.25,0.25,0.25],...
% $$$          'MarkerSize',1);
% $$$     % PLOT trajectory points at spike times
% $$$     scatter(sax, ...
% $$$             xyz{exampleUnit.trialIndex}(spktemp.res,'hcom',1),...
% $$$             xyz{exampleUnit.trialIndex}(spktemp.res,'hcom',2),...
% $$$             1, ...
% $$$             pclr(spkPhzInd,:));
% $$$     
% $$$     % FORMAT subplot        
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     xlim(sax,exampleUnit.close.Xlims); 
% $$$     ylim(sax,exampleUnit.close.Ylims); 
% $$$     daspect(sax,[1,1,1]);
% $$$     
% $$$     savefig(pfig, fullfile(partsPath,[sax.Tag,'.fig']));
% $$$     close(pfig);
%end

%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< TRAJECTORIES { x, y | (theta descend) } -----------------------------------------
%function EgoProCode2D_f1_subplot_traj_allo_theta_descend()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_traj_allo_theta_descend');
% $$$     
% $$$     % PLOT trajectories
% $$$     plot(sax, ...
% $$$          xyz{exampleUnit.trialIndex}(sper,'hcom',1),...
% $$$          xyz{exampleUnit.trialIndex}(sper,'hcom',2),...
% $$$          '.',...
% $$$          'MarkerFaceColor',[0.25,0.25,0.25],...
% $$$          'MarkerEdgeColor',[0.25,0.25,0.25],...
% $$$          'MarkerSize',1);
% $$$     % PLOT trajectory points at spike times
% $$$     scatter(sax, ...
% $$$             xyz{exampleUnit.trialIndex}(spktemp.res(spkPhzInd==1),'hcom',1),...
% $$$             xyz{exampleUnit.trialIndex}(spktemp.res(spkPhzInd==1),'hcom',2),...
% $$$             1, ...
% $$$             pclr(spkPhzInd(spkPhzInd==1),:));
% $$$     
% $$$     % FORMAT subplot        
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     xlim(sax,exampleUnit.close.Xlims); 
% $$$     ylim(sax,exampleUnit.close.Ylims); 
% $$$     daspect(sax,[1,1,1]);
% $$$     
% $$$     savefig(pfig, fullfile(partsPath,[sax.Tag,'.fig']));
% $$$     close(pfig);
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< TRAJECTORIES { x, y | (theta trough) } ------------------------------------------
%function EgoProCode2D_f1_subplot_traj_allo_theta_trough()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_traj_allo_theta_trough');
% $$$     
% $$$     % PLOT trajectories
% $$$     plot(sax, ...
% $$$          xyz{exampleUnit.trialIndex}(sper,'hcom',1),...
% $$$          xyz{exampleUnit.trialIndex}(sper,'hcom',2),...
% $$$          '.',...
% $$$          'MarkerFaceColor',[0.25,0.25,0.25],...
% $$$          'MarkerEdgeColor',[0.25,0.25,0.25],...
% $$$          'MarkerSize',1);
% $$$     % PLOT trajectory points at spike times
% $$$     scatter(sax, ...
% $$$             xyz{exampleUnit.trialIndex}(spktemp.res(spkPhzInd==2),'hcom',1),...
% $$$             xyz{exampleUnit.trialIndex}(spktemp.res(spkPhzInd==2),'hcom',2),...
% $$$             1, ...
% $$$             pclr(spkPhzInd(spkPhzInd==2),:));
% $$$     
% $$$     % FORMAT subplot        
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     xlim(sax,exampleUnit.close.Xlims); 
% $$$     ylim(sax,exampleUnit.close.Ylims); 
% $$$     daspect(sax,[1,1,1]);
% $$$     
% $$$     savefig(pfig, fullfile(partsPath,[sax.Tag,'.fig']));
% $$$     close(pfig);
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< TRAJECTORIES { x, y | (theta ascend) } ------------------------------------------
%function EgoProCode2D_f1_subplot_traj_allo_theta_ascend()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_traj_allo_theta_ascend');
% $$$     
% $$$     % PLOT trajectories
% $$$     plot(sax, ...
% $$$          xyz{exampleUnit.trialIndex}(sper,'hcom',1),...
% $$$          xyz{exampleUnit.trialIndex}(sper,'hcom',2),...
% $$$          '.',...
% $$$          'MarkerFaceColor',[0.25,0.25,0.25],...
% $$$          'MarkerEdgeColor',[0.25,0.25,0.25],...
% $$$          'MarkerSize',1);
% $$$     % PLOT trajectory points at spike times
% $$$     scatter(sax, ...
% $$$             xyz{exampleUnit.trialIndex}(spktemp.res(spkPhzInd==3),'hcom',1),...
% $$$             xyz{exampleUnit.trialIndex}(spktemp.res(spkPhzInd==3),'hcom',2),...
% $$$             1, ...
% $$$             pclr(spkPhzInd(spkPhzInd==3),:));
% $$$     
% $$$     % FORMAT subplot        
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     xlim(sax,exampleUnit.close.Xlims); 
% $$$     ylim(sax,exampleUnit.close.Ylims); 
% $$$     daspect(sax,[1,1,1]);
% $$$     
% $$$     savefig(pfig, fullfile(partsPath,[sax.Tag,'.fig']));
% $$$     close(pfig);
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< EgoProCode2D_f1_subplot_trajectory_allo_head_direction --------------------------

%function EgoProCode2D_f1_subplot_trajectory_allo_head_direction()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_trajectory_allo_head_direction');
% $$$ 
% $$$     % PLOT allo-ratemap, trajectory (color=thetaPhase), and skeletons
% $$$     plot(pft{exampleUnit.trialIndex},units{exampleUnit.trialIndex}(exampleUnit.index),1,'text',[0,mrate],'colorMap',@bone,'nanColor',[1,1,1]);
% $$$     scatter(xyz{exampleUnit.trialIndex}(exampleUnit.trajectoryTimeSeries,'hcom',1),...
% $$$             xyz{exampleUnit.trialIndex}(exampleUnit.trajectoryTimeSeries,'hcom',2),...
% $$$             4,...
% $$$             atan2(pfstrj(exampleUnit.trajectoryTimeSeries,2),pfstrj(exampleUnit.trajectoryTimeSeries,1)),'filled')
% $$$     circle(mxp(1),mxp(2),200,'-r');
% $$$     plotSkeleton(Trials{exampleUnit.trialIndex},xyz{exampleUnit.trialIndex},exampleUnit.trajectoryTimeSeries(1));        
% $$$     plotSkeleton(Trials{exampleUnit.trialIndex},xyz{exampleUnit.trialIndex},exampleUnit.trajectoryTimeSeries(round(length(exampleUnit.trajectoryTimeSeries)/1.15)));    
% $$$ 
% $$$     % FORMAT subplot
% $$$     xlim(sax,exampleUnit.close.Xlims); 
% $$$     ylim(sax,exampleUnit.close.Ylims); 
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     colormap(sax,'hsv');    
% $$$     caxis([-pi,pi]);        
% $$$     daspect([1,1,1]);    
% $$$     grid('on');    
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end

%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< EgoProCode2D_f1_subplot_trajectory_ego_head_direction() -------------------------
%function EgoProCode2D_f1_subplot_trajectory_ego_head_direction()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_trajectory_ego_head_direction');
% $$$ 
% $$$     % PLOT ego-ratemap, 
% $$$     plot(pfe{exampleUnit.trialIndex},units{exampleUnit.trialIndex}(exampleUnit.index),1,'',[0,mrate],'mazeMaskFlag',false,'colorMap',@bone,'flipAxesFlag',true);
% $$$     patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$     patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$     line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color,'LineWidth',2);
% $$$     line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color,'LineWidth',2);
% $$$     circle(0,0,200,'r-');
% $$$     scatter(pfstrj(exampleUnit.trajectoryTimeSeries,2),                                                        ...
% $$$             pfstrj(exampleUnit.trajectoryTimeSeries,1),                                                        ...
% $$$             4,                                                                                                 ...
% $$$             atan2(pfstrj(exampleUnit.trajectoryTimeSeries,2),pfstrj(exampleUnit.trajectoryTimeSeries,1)),'filled')
% $$$ 
% $$$     % FORMAT subplot
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     xlim    (sax,[-300,300]);    
% $$$     ylim    (sax,[-300,300]);
% $$$     colormap(sax,'hsv');
% $$$     caxis   (sax,[-pi,pi]);            
% $$$     daspect (sax,[1,1,1]);
% $$$     grid    (sax,'on');
% $$$     
% $$$     savefig(pfig, fullfile(partsPath,[sax.Tag, '.fig']));
% $$$     close(pfig);
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< EgoProCode2D_f1_subplot_trajectories_ego ----------------------------------------
%function EgoProCode2D_f1_subplot_trajectories_ego()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_trajectories_ego');
% $$$ 
% $$$     % SUBPLOT <- ego-traj, traj, spktraj(color=thetaPhase), rat model, circle
% $$$     plot(sax, ...
% $$$          pfstrj(sper,2),...
% $$$          pfstrj(sper,1),...
% $$$          '.',...
% $$$          'MarkerFaceColor',[0.25,0.25,0.25],...
% $$$          'MarkerEdgeColor',[0.25,0.25,0.25],...
% $$$          'MarkerSize',1);
% $$$     scatter(sax, ...
% $$$             pfstrj(spktemp.res,2),...
% $$$             pfstrj(spktemp.res,1),...
% $$$             1,...
% $$$             pclr(spkPhzInd,:));
% $$$     patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$     patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$     line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color,'LineWidth',2);
% $$$     line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color,'LineWidth',2);
% $$$     circle(0,0,200,'r-');
% $$$ 
% $$$     % FORMAT subplot
% $$$     xlim(sax,[-300,300]);    
% $$$     ylim(sax,[-300,300]);
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));        
% $$$     close(pfig);        
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< EgoProCode2D_f1_subplot_egofieldExample -----------------------------------------
%function EgoProCode2D_f1_subplot_egofieldExample()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_egofieldExample');
% $$$ 
% $$$     % SUBPLOT <- ego-ratemap, rat model, circle
% $$$     plot(pfe{exampleUnit.trialIndex},exampleUnit.id,1,'',[0,mrate],'colorMap',@jet,'mazeMaskFlag',false,'flipAxesFlag',true);
% $$$     patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$     patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$     line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color,'LineWidth',2);
% $$$     line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color,'LineWidth',2);
% $$$     circle(0,0,200,'r-');
% $$$ 
% $$$     % FORMAT subplot
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     xlim(sax,[-300,300]);
% $$$     ylim(sax,[-300,300]);
% $$$     daspect(sax,[1,1,1]);
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< EgoProCode2D_f1_subplot_traj_ego_theta_descend ----------------------------------
%function  EgoProCode2D_f1_subplot_traj_ego_theta_descend()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_traj_ego_theta_descend');
% $$$ 
% $$$     % SUBPLOT <- ego-traj, traj, spktraj(color=thetaPhase), rat model, circle
% $$$     phzInd = 1;
% $$$     plot(sax, ...
% $$$          pfstrj(sper,2),...
% $$$          pfstrj(sper,1),...
% $$$          '.',...
% $$$          'MarkerFaceColor',[0.25,0.25,0.25],...
% $$$          'MarkerEdgeColor',[0.25,0.25,0.25],...
% $$$          'MarkerSize',1);
% $$$     % overlay spikes    
% $$$     plot(sax, ...
% $$$          pfstrj(spktemp.res(spkPhzInd==phzInd),2),...
% $$$          pfstrj(spktemp.res(spkPhzInd==phzInd),1),...
% $$$          '.',...
% $$$          'MarkerFaceColor',pclr(spkPhzInd(spkPhzInd==phzInd),:),...
% $$$          'MarkerEdgeColor',pclr(spkPhzInd(spkPhzInd==phzInd),:),...
% $$$          'MarkerSize',1);
% $$$     patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$     patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$     line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color,'LineWidth',2);
% $$$     line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color,'LineWidth',2);
% $$$     circle(0,0,200,'r-');
% $$$ 
% $$$     % FORMAT subplot
% $$$     xlim(sax,[-300,300]);    
% $$$     ylim(sax,[-300,300]);
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));        
% $$$     close(pfig);        
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< EgoProCode2D_f1_subplot_traj_ego_theta_trough -----------------------------------
%function  EgoProCode2D_f1_subplot_traj_ego_theta_trough()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_traj_ego_theta_trough');
% $$$ 
% $$$     % SUBPLOT <- ego-traj, traj, spktraj(color=thetaPhase), rat model, circle
% $$$     phzInd = 2;
% $$$     plot(sax, ...
% $$$          pfstrj(sper,2),...
% $$$          pfstrj(sper,1),...
% $$$          '.',...
% $$$          'MarkerFaceColor',[0.25,0.25,0.25],...
% $$$          'MarkerEdgeColor',[0.25,0.25,0.25],...
% $$$          'MarkerSize',1);
% $$$     % overlay spikes    
% $$$     plot(sax, ...
% $$$          pfstrj(spktemp.res(spkPhzInd==phzInd),2),...
% $$$          pfstrj(spktemp.res(spkPhzInd==phzInd),1),...
% $$$          '.',...
% $$$          'MarkerFaceColor',pclr(spkPhzInd(spkPhzInd==phzInd),:),...
% $$$          'MarkerEdgeColor',pclr(spkPhzInd(spkPhzInd==phzInd),:),...
% $$$          'MarkerSize',1);
% $$$     patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$     patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$     line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color,'LineWidth',2);
% $$$     line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color,'LineWidth',2);
% $$$     circle(0,0,200,'r-');
% $$$ 
% $$$     % FORMAT subplot
% $$$     xlim(sax,[-300,300]);    
% $$$     ylim(sax,[-300,300]);
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));        
% $$$     close(pfig);        
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< EgoProCode2D_f1_subplot_traj_ego_theta_ascend -----------------------------------
%function  EgoProCode2D_f1_subplot_traj_ego_theta_ascend()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_traj_ego_theta_ascend');
% $$$ 
% $$$     % SUBPLOT <- ego-traj, traj, spktraj(color=thetaPhase), rat model, circle
% $$$     phzInd = 3;
% $$$     plot(sax, ...
% $$$          pfstrj(sper,2),...
% $$$          pfstrj(sper,1),...
% $$$          '.',...
% $$$          'MarkerFaceColor',[0.25,0.25,0.25],...
% $$$          'MarkerEdgeColor',[0.25,0.25,0.25],...
% $$$          'MarkerSize',1);
% $$$     % overlay spikes
% $$$     plot(sax, ...
% $$$          pfstrj(spktemp.res(spkPhzInd==phzInd),2),...
% $$$          pfstrj(spktemp.res(spkPhzInd==phzInd),1),...
% $$$          '.',...
% $$$          'MarkerFaceColor',pclr(spkPhzInd(spkPhzInd==phzInd),:),...
% $$$          'MarkerEdgeColor',pclr(spkPhzInd(spkPhzInd==phzInd),:),...
% $$$          'MarkerSize',1);
% $$$     patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$     patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$     line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color,'LineWidth',2);
% $$$     line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color,'LineWidth',2);
% $$$     circle(0,0,200,'r-');
% $$$ 
% $$$     % FORMAT subplot
% $$$     xlim(sax,[-300,300]);    
% $$$     ylim(sax,[-300,300]);
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));        
% $$$     close(pfig);        
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< DIAGRAM Theta phase line reference ----------------------------------------------
%function EgoProCode2D_f1_subplot_theta_cycle_horizontal()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_theta_cycle_horizontal');
% $$$     
% $$$     % SUBPLOT <- theta cycle (color=partions)
% $$$     plot(sax,linspace(0,1),cos(linspace(0,2*pi)),'k','LineWidth',2);
% $$$     plims = [0.5,              2.26106176905986];
% $$$     plot(sax,...
% $$$          linspace(plims(1)/(2*pi),plims(2)/(2*pi)),...
% $$$          cos(linspace(plims(1),plims(2))), 'Color',pclr(1,:), 'LineWidth',2);
% $$$     plims = [2.26106176905986, 4.02212353811972];
% $$$     plot(sax,...
% $$$          linspace(plims(1)/(2*pi),plims(2)/(2*pi)),...
% $$$          cos(linspace(plims(1),plims(2))), 'Color',pclr(2,:), 'LineWidth',2);
% $$$     plims = [4.02212353811972, 5.78318530717959];
% $$$     plot(sax,...
% $$$          linspace(plims(1)/(2*pi),plims(2)/(2*pi)),...
% $$$          cos(linspace(plims(1),plims(2))), 'Color',pclr(3,:),'LineWidth',2);
% $$$ 
% $$$     
% $$$     %FORMAT subplot
% $$$     ylim(sax,[-1,1]);
% $$$     xlim(sax,[0,1]);
% $$$     sax.Visible = 'off';
% $$$     text(sax,0.05,0.5,'0','FontSize',8);
% $$$     text(sax,0.5,0.5,'\pi','FontSize',8);
% $$$     text(sax,0.95,0.5,'2\pi','FontSize',8);
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< DIAGRAM Theta phase line reference ----------------------------------------------
%function EgoProCode2D_f1_subplot_theta_cycle_vertical()
% $$$     [pfig, sax] =  setup_figure_('EgoProCode2D_f1_subplot_theta_cycle_vertical');
% $$$ 
% $$$     % SUBPLOT <- theta cycle (color=partions)
% $$$     plot(sax,-cos(linspace(0,2*pi)),linspace(0,1),'k','LineWidth',2);
% $$$     plims = [0.5,              2.26106176905986];
% $$$     plot(sax,...
% $$$          -cos(linspace(plims(1),plims(2))),...
% $$$          linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(1,:), 'LineWidth',2);
% $$$     plims = [2.26106176905986, 4.02212353811972];
% $$$     plot(sax,...
% $$$          -cos(linspace(plims(1),plims(2))),...
% $$$          linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(2,:), 'LineWidth',2);
% $$$     plims = [4.02212353811972, 5.78318530717959];
% $$$     plot(sax,...
% $$$          -cos(linspace(plims(1),plims(2))),...
% $$$          linspace(plims(1)/(2*pi),plims(2)/(2*pi)), 'Color',pclr(3,:),'LineWidth',2);
% $$$ 
% $$$     % FORMAT subplot
% $$$     ylim([0,1]);
% $$$     sax.Visible = 'off';
% $$$     text(0.5,0.05,'0','Rotation',90,'FontSize',8);
% $$$     text(-0.5,0.5,'\pi','Rotation',90,'FontSize',8);
% $$$     text(0.5,0.95,'2\pi','Rotation',90,'FontSize',8);
% $$$     
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< EGOFIELD BY PHASE { xEgo, yEgo | { θ | θ ∈ [0,2π/3] }, Color:Hz } ---------------
%function EgoProCode2D_f1_subplot_egofieldExample_descending_phase()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_egofieldExample_descending_phase');
% $$$ 
% $$$     pfig = figure();sax = axes(); sax.Tag = 'EgoProCode2D_f1_subplot_egofieldExample_descending_phase';hold(sax,'on');
% $$$     
% $$$     plot(pfet{exampleUnit.trialIndex}{1},exampleUnit.id,1,'text',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);
% $$$ 
% $$$     % FORMAT subplot
% $$$     axis(sax,'tight');
% $$$     xlim(sax,[-300,300]);
% $$$     ylim(sax,[-250,350]);
% $$$     Lines([],0,'w');
% $$$     Lines(0,[],'w');
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     
% $$$     % ADD rat model
% $$$     subject = struct(rat);
% $$$     subject = update_subject_patch(subject,'head', [], false,[],[]);
% $$$     subject = update_subject_patch(subject,'body',[],false,[],[]);
% $$$     patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
% $$$     patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
% $$$     line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
% $$$     line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end

%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< EGOFIELD BY PHASE { xEgo, yEgo | { θ | θ ∈ [2π/3,4π/3] }, Color:Hz } ------------
%function EgoProCode2D_f1_subplot_egofieldExample_trough_phase()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_egofieldExample_trough_phase');
% $$$     
% $$$     plot(pfet{exampleUnit.trialIndex}{2},exampleUnit.id,1,'text',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);
% $$$ 
% $$$     % FORMAT subplot
% $$$     axis(sax,'tight');
% $$$     xlim(sax,[-300,300]);
% $$$     ylim(sax,[-250,350]);
% $$$     Lines([],0,'w');
% $$$     Lines(0,[],'w');
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     
% $$$     % ADD rat model to the subplot
% $$$     subject = struct(rat);
% $$$     subject = update_subject_patch(subject,'head', [], false,[],[]);
% $$$     subject = update_subject_patch(subject,'body',[],false,[],[]);
% $$$     patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
% $$$     patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
% $$$     line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
% $$$     line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
    %end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< EGOFIELD BY PHASE { xEgo, yEgo | { θ | θ ∈ [0,2π/3] }, Color:Hz } ---------------

%function EgoProCode2D_f1_subplot_egofieldExample_ascending_phase()
% $$$     [pfig, sax] = setup_figure_( 'EgoProCode2D_f1_subplot_egofieldExample_ascending_phase' );
% $$$ 
% $$$     % PLOT egocentric rate map
% $$$     plot(pfet{exampleUnit.trialIndex}{3},exampleUnit.id,1,'text',[0,mrate],false,'colorMap',@jet,'flipAxesFlag',true);
% $$$ 
% $$$     % FORMAT subplot
% $$$     axis(sax,'tight');
% $$$     xlim(sax,[-300,300]);
% $$$     ylim(sax,[-250,350]);
% $$$     Lines([],0,'w');
% $$$     Lines(0,[],'w');
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     
% $$$     % ADD rat model to subplot
% $$$     subject = struct(rat);
% $$$     subject = update_subject_patch(subject,'head', [], false,[],[]);
% $$$     subject = update_subject_patch(subject,'body',[],false,[],[]);
% $$$     patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
% $$$     patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
% $$$     line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
% $$$     line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end

%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< ALLOFIELD BY PHASE { xEgo, yEgo | { θ | θ ∈ [0,2π/3] }, Color:Hz } ---------------
%function EgoProCode2D_f1_subplot_allofieldExample_descending_phase()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_allofieldExample_descending_phase');
% $$$ 
% $$$     plot(ratemapsAlloThp{exampleUnit.trialIndex}{1},exampleUnit.id,1,'text',[0,mrate],false,'colorMap',@jet);
% $$$ 
% $$$     % FORMAT subplot
% $$$     axis(sax,'tight');
% $$$     
% $$$     %xlim(sax,[-300,300]);
% $$$     %ylim(sax,[-300,300]);
% $$$     %circle(mxp(1),mxp(2),200,'-r');            
% $$$     xlim(sax,[-300,300]+mxp(1));
% $$$     ylim(sax,[-300,300]+mxp(2));
% $$$     
% $$$     %Lines([],mxp(2),'w');
% $$$     %Lines(mxp(1),[],'w');
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end

%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< ALLOFIELD BY PHASE { xEgo, yEgo | { θ | θ ∈ [2π/3,4π/3] }, Color:Hz } ------------

%function EgoProCode2D_f1_subplot_allofieldExample_trough_phase()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_allofieldExample_trough_phase');
% $$$     
% $$$     plot(ratemapsAlloThp{exampleUnit.trialIndex}{2},exampleUnit.id,1,'text',[0,mrate],false,'colorMap',@jet);
% $$$ 
% $$$     % FORMAT subplot
% $$$     axis(sax,'tight');
% $$$     xlim(sax,[-300,300]+mxp(1));
% $$$     ylim(sax,[-300,300]+mxp(2));
% $$$     %Lines([],mxp(2),'w');
% $$$     %Lines(mxp(1),[],'w');
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
% $$$     %end

%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< ALLOFIELD BY PHASE { xEgo, yEgo | { θ | θ ∈ [0,2π/3] }, Color:Hz } ---------------

%function EgoProCode2D_f1_subplot_allofieldExample_ascending_phase()
% $$$     [pfig, sax] = setup_figure_( 'EgoProCode2D_f1_subplot_allofieldExample_ascending_phase' );
% $$$ 
% $$$     % PLOT egocentric rate map
% $$$     plot(ratemapsAlloThp{exampleUnit.trialIndex}{3},exampleUnit.id,1,'text',[0,mrate],false,'colorMap',@jet);
% $$$ 
% $$$     % FORMAT subplot
% $$$     axis(sax,'tight');
% $$$     xlim(sax,[-300,300]+mxp(1));
% $$$     ylim(sax,[-300,300]+mxp(2));
% $$$     %Lines([],mxp(2),'w');
% $$$     %Lines(mxp(1),[],'w');
% $$$     sax.XTickLabel =[];
% $$$     sax.YTickLabel =[];
% $$$     
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end

%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< FORWARD egoMeanRmapPos CA1 Descending -------------------------------------------
%function EgoProCode2D_f1_subplot_forward_field_distrib_descending_phase()
% $$$     [pfig, sax] = setup_figure_( 'EgoProCode2D_f1_subplot_forward_field_distrib_descending_phase' );
% $$$     
% $$$     % PLOT distribution of forward direction of rate map centers
% $$$     hold(sax,'on');
% $$$     [h,L,MX,MED,bw] = violin({(egoMeanRmapPos(uidsCA1,1,1)-2),...
% $$$                                egoMeanRmapPos(uidsCA3,1,1)-2});
% $$$     
% $$$     % FORMAT subplot
% $$$     delete(L);
% $$$     af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr([1,1],:),ones([2,1]))');
% $$$     xlim(sax,[0.5,2.5])
% $$$     %ylim(sax,[-15,15])
% $$$     ylim(sax,[-12.5,17.5]);    
% $$$     sax.YTickLabelMode = 'Manual';
% $$$     sax.YTickLabel = {};
% $$$     sax.YTick = [-15,-7.5,0,7.5,15];
% $$$     %sax.XTickLabel = {'CA1','CA3'};
% $$$     grid(sax,'on');
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< FORWARD egoMeanRmapPos CA1 Trough -----------------------------------------------
%function EgoProCode2D_f1_subplot_forward_field_distrib_trough_phase()
% $$$     [pfig, sax] = setup_figure_( 'EgoProCode2D_f1_subplot_forward_field_distrib_trough_phase');
% $$$     
% $$$     % PLOT distribution of forward direction of rate map centers
% $$$     hold(sax,'on');
% $$$     [h,L,MX,MED,bw] = violin({(egoMeanRmapPos(uidsCA1,2,1)-2),...
% $$$                         egoMeanRmapPos(uidsCA3,2,1)-2});
% $$$     
% $$$     % FORMAT subplot
% $$$     delete(L);
% $$$     af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr([2,2],:),ones([2,1]))');
% $$$     xlim([0.5,2.5])
% $$$     %ylim(sax,[-15,15])
% $$$     ylim(sax,[-12.5,17.5]);    
% $$$     
% $$$     sax.YTickLabelMode = 'Manual';
% $$$     sax.YTickLabel = {};
% $$$     sax.YTick = [-15,-7.5,0,7.5,15];
% $$$     sax.XTickLabel = {};
% $$$     grid(sax,'on');
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< FORWARD egoMeanRmapPos CA1 Ascending -----------------------------------------------

%function EgoProCode2D_f1_subplot_forward_field_distrib_ascending_phase()
% $$$     [pfig, sax] = setup_figure_( 'EgoProCode2D_f1_subplot_forward_field_distrib_ascending_phase');
% $$$     
% $$$     % PLOT distribution of forward direction of rate map centers
% $$$     hold(sax,'on');
% $$$     [h,L,MX,MED,bw] = violin({(egoMeanRmapPos(uidsCA1,3,1)-2),...
% $$$                                egoMeanRmapPos(uidsCA3,3,1)-2});
% $$$     
% $$$     % FORMAT subplot
% $$$     delete(L);
% $$$     af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr([3,3],:),ones([2,1]))');
% $$$     xlim([0.5,2.5])
% $$$     %ylim(sax,[-15,15])
% $$$     ylim(sax,[-12.5,17.5]);    
% $$$     
% $$$     sax.YTickLabelMode = 'Manual';
% $$$     sax.YTickLabel = {};
% $$$     sax.YTick = [-15,-7.5,0,7.5,15];
% $$$     sax.XTickLabel = {};
% $$$     grid(sax,'on');
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end

%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< LATERAL mean field position distribution, descending ----------------------------
%function EgoProCode2D_f1_subplot_lateral_field_distrib_descending_phase()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_lateral_field_distrib_descending_phase');
% $$$     
% $$$     % PLOT distribution of lateral direction of rate map centers
% $$$     [h,L,MX,MED,bw] = violin({egoMeanRmapPos(uidsCA1,1,2)+0.8,...
% $$$                               egoMeanRmapPos(uidsCA3,1,2)+0.8});
% $$$ 
% $$$     % FORMAT subplot
% $$$     delete(L);
% $$$     af(@(hndl) set(hndl,'Vertices',fliplr(hndl.Vertices)), h);
% $$$     af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr([1,1],:),ones([2,1]))');
% $$$     hlns = findobj(sax,'Type','line');
% $$$     for ln = 1:numel(hlns)
% $$$         [hlns(ln).XData,hlns(ln).YData] = deal(hlns(ln).YData,hlns(ln).XData);
% $$$     end
% $$$ 
% $$$     ylim(sax,[0.5,2.5]);
% $$$     xlim(sax,[-30,30]);
% $$$     
% $$$     sax.XTickLabelMode = 'Auto';
% $$$     sax.XTick = [-20,-10,0,10,20];
% $$$     
% $$$     xlabel(sax,'cm');
% $$$     sax.YTickLabel = {};    
% $$$ 
% $$$     grid(sax,'on');
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< LATERAL mean field position distribution, trough --------------------------------
%function EgoProCode2D_f1_subplot_lateral_field_distrib_trough_phase()    
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_lateral_field_distrib_trough_phase');
% $$$ 
% $$$     % PLOT distribution of lateral direction of rate map centers
% $$$     [h,L,MX,MED,bw] = violin({egoMeanRmapPos(uidsCA1,2,2)+0.8,...
% $$$                               egoMeanRmapPos(uidsCA3,2,2)+0.8});
% $$$ 
% $$$     % FORMAT subplot
% $$$     delete(L);
% $$$     af(@(hndl) set(hndl,'Vertices',fliplr(hndl.Vertices)), h);
% $$$     af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr([2,2],:),ones([2,1]))');
% $$$     hlns = findobj(sax,'Type','line');
% $$$     for ln = 1:numel(hlns)
% $$$         [hlns(ln).XData,hlns(ln).YData] = deal(hlns(ln).YData,hlns(ln).XData);
% $$$     end
% $$$     xlim(sax,[-30,30]);
% $$$     ylim(sax,[0.5,2.5]);
% $$$     
% $$$     sax.XTickLabelMode = 'Auto';
% $$$     sax.XTick = [-20,-10,0,10,20];
% $$$     
% $$$     xlabel(sax,'cm');
% $$$     
% $$$     sax.YTickLabel = {};
% $$$     grid(sax,'on');
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< LATERAL mean field position distribution, ascending -----------------------------

%function EgoProCode2D_f1_subplot_lateral_field_distrib_ascending_phase()    
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_lateral_field_distrib_ascending_phase');
% $$$ 
% $$$     % PLOT distribution of lateral direction of rate map centers
% $$$     [h,L,MX,MED,bw] = violin({egoMeanRmapPos(uidsCA1,3,2)+0.8,...
% $$$                               egoMeanRmapPos(uidsCA3,3,2)+0.8});
% $$$ 
% $$$     % FORMAT subplot
% $$$     delete(L);
% $$$     af(@(hndl) set(hndl,'Vertices',fliplr(hndl.Vertices)), h);
% $$$     af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr([3,3],:),ones([2,1]))');
% $$$     hlns = findobj(sax,'Type','line');
% $$$     for ln = 1:numel(hlns)
% $$$         [hlns(ln).XData,hlns(ln).YData] = deal(hlns(ln).YData,hlns(ln).XData);
% $$$     end
% $$$     
% $$$     xlim(sax,[-30,30]);
% $$$     ylim(sax,[0.5,2.5]);
% $$$     
% $$$     sax.XTickLabelMode = 'Auto';
% $$$     sax.XTick = [-20,-10,0,10,20];
% $$$     
% $$$     xlabel(sax,'cm');
% $$$     
% $$$     sax.YTickLabel = {};
% $$$     grid(sax,'on');
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%end

%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< RATIOS - firing rate ratios vs trough CA1 ---------------------------------------
%figure EgoProCode2D_f1_subplot_ca1_rate_ratios()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_ca1_rate_ratios');
% $$$ 
% $$$     % SUBPLOT <- rate ratios
% $$$     plot(egoMaxRmapRate(uidsCA1,1)./egoMaxRmapRate(uidsCA1,2),...
% $$$          egoMaxRmapRate(uidsCA1,3)./egoMaxRmapRate(uidsCA1,2),'.')
% $$$ 
% $$$     % FORMAT subplot
% $$$     xlim(sax,[0,1.8]);
% $$$     ylim(sax,[0,1.8]);
% $$$     grid(sax,'on');
% $$$     title(sax(end),'Rate Ratio CA1');
% $$$     xlabel(sax(end),'Dsc/Trough (A.U.)');
% $$$     ylabel(sax(end),'Asc/Trough (A.U.)');
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
% $$$     %end
%%%>>>    
%--------------------------------------------------------------------------------------- TRANS
%%%<<< FIELD SIZE STATS CA1 ------------------------------------------------------------
%function EgoProCode2D_f1_subplot_ca1_field_size()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_ca1_field_size');
% $$$     
% $$$     % SUBPLOT <- Asce, Desc field size vs Trough field size
% $$$     plot((egoSize(uidsCA1,2))*0.02^2,(egoSize(uidsCA1,1,1))*0.02^2,'.','Color',pclr(1,:));
% $$$     plot((egoSize(uidsCA1,2))*0.02^2,(egoSize(uidsCA1,3,1))*0.02^2,'.','Color',pclr(3,:));
% $$$     
% $$$     % FORMAT subplot
% $$$     line([0,0.3],[0,0.3],'Color','k')
% $$$     grid(sax(end),'on');
% $$$     title(sax(end),'Field Area CA1')
% $$$     xlabel('Trough Field Size (m^2)');
% $$$     ylabel('asce and desc Field Size (m^2)');
% $$$ 
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
% $$$     %end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< RATIOS - firing rate ratios vs trough CA3 ---------------------------------------
%figure EgoProCode2D_f1_subplot_ca3_rate_ratios()    
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_ca3_rate_ratios');
% $$$ 
% $$$     % SUBPLOT <- rate ratios
% $$$     plot(egoMaxRmapRate(uidsCA3,1)./egoMaxRmapRate(uidsCA3,2),...
% $$$          egoMaxRmapRate(uidsCA3,3)./egoMaxRmapRate(uidsCA3,2),'.')
% $$$ 
% $$$     % FORMAT subplot
% $$$     xlim(sax,[0,1.8]);
% $$$     ylim(sax,[0,1.8]);
% $$$     grid(sax,'on');
% $$$     title(sax(end),'Rate Ratio CA3');
% $$$     xlabel(sax(end),'Dsc/Trough (A.U.)');
% $$$     ylabel(sax(end),'Asc/Trough (A.U.)');
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
    %end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< FIELD SIZE STATS CA3 ------------------------------------------------------------
%function EgoProCode2D_f1_subplot_ca3_field_size()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_ca3_field_size');
% $$$     
% $$$     % SUBPLOT <- Asce, Desc field size vs Trough field size
% $$$     plot((egoSize(uidsCA3,2))*0.02^2,(egoSize(uidsCA3,1,1))*0.02^2,'.','Color',pclr(1,:));
% $$$     plot((egoSize(uidsCA3,2))*0.02^2,(egoSize(uidsCA3,3,1))*0.02^2,'.','Color',pclr(3,:));
% $$$     
% $$$     % FORMAT subplot
% $$$     line([0,0.3],[0,0.3],'Color','k')
% $$$     grid(sax(end),'on');
% $$$     title(sax(end),'Field Area CA3')
% $$$     xlabel('Trough Field Size (m^2)');
% $$$     ylabel('asce and desc Field Size (m^2)');
% $$$ 
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
% $$$     %end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS    
%%%<<< LFP timeseries example ----------------------------------------------------------
% function EgoProCode2D_f1_subplot_lfp_timeseries()
    [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_lfp_timeseries');

    timeIndex = 87000;
    lfp = Trials{exampleUnit.trialIndex}.load('lfp',Trials{exampleUnit.trialIndex}.meta.channelGroup.theta);
    lfp.resample(sampleRate);
    lfp.data = nunity(lfp.data);
    phz = load_theta_phase(Trials{exampleUnit.trialIndex},sampleRate);
    mres = spk{exampleUnit.trialIndex}(exampleUnit.id); 
    mres  = mres(WithinRanges(mres,timeIndex+[-200,200]));

    
    plot(([-200:200])./sampleRate,lfp(timeIndex+[-200:200],1)+3,'k','LineWidth',1)
    plot(([-200:200])./sampleRate,phz(timeIndex+[-200:200],1),'g','LineWidth',1)
    scatter((mres-timeIndex)./sampleRate,lfp(mres,1)+3,4,'r','filled');
    sax(end).Visible = 'off';

    savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
    close(pfig);
    %end
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< IMSCALE maxrate stats CA1 -------------------------------------------------------
% function EgoProCode2D_f1_ego_maxrate_stats_CA1()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_ego_maxrate_stats_ca1');
% $$$ 
% $$$     rateRows = log10(egoMaxRmapRate(uidsCA1,:));
% $$$     [mx,mind] = max(log10(egoMaxRmapRate(uidsCA1,:)),[],2);
% $$$     newrows=[];
% $$$     for m = 1:3 
% $$$         [~,sind] = sort(mx(mind==m));
% $$$         myrows = rateRows(mind==m,:);
% $$$         newrows = cat(1, newrows,myrows(sind,:));
% $$$     end
% $$$ 
% $$$     imagesc(sax,newrows')
% $$$     axis(sax,'tight');
% $$$     colormap(sax,'jet');    
% $$$     sax.YTickLabel = {};
% $$$     sax.XTickLabel = {};
% $$$     caxis(sax,[0,1.4]);
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%%%>>>
%--------------------------------------------------------------------------------------- TRANS    
%%%<<< IMSCALE maxrate stats CA3 -------------------------------------------------------
% function EgoProCode2D_f1_ego_maxrate_stats_CA3()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_ego_maxrate_stats_ca3');
% $$$ 
% $$$     rateRows = log10(egoMaxRmapRate(uidsCA3,:));
% $$$     [mx,mind] = max(log10(egoMaxRmapRate(uidsCA3,:)),[],2);
% $$$     newrows=[];
% $$$     for m = 1:3 
% $$$         [~,sind] = sort(mx(mind==m));
% $$$         myrows = rateRows(mind==m,:);
% $$$         newrows = cat(1, newrows,myrows(sind,:));
% $$$     end
% $$$ 
% $$$     imagesc(sax,newrows');
% $$$     axis(sax,'tight');    
% $$$     colormap(sax,'jet');
% $$$     cax = colorbar();
% $$$     ylabel(cax,'log10(Hz)');
% $$$     sax.YTickLabel = {};
% $$$     sax.XTickLabel = {};
% $$$     caxis(sax,[0,1.4]);    
% $$$ 
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< IMSCALE field size stats CA1 ----------------------------------------------------
% function EgoProCode2D_f1_ego_field_size_stats_CA1()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_ego_field_size_stats_ca1');
% $$$ 
% $$$     rateRows = log10(egoSize(uidsCA1,:)*0.02^2);
% $$$     [mx,mind] = max(log10(egoSize(uidsCA1,:)*0.02^2),[],2);
% $$$     newrows=[];
% $$$     for m = 1:3 
% $$$         [~,sind] = sort(mx(mind==m));
% $$$         myrows = rateRows(mind==m,:);
% $$$         newrows = cat(1, newrows,myrows(sind,:));
% $$$     end
% $$$     
% $$$     imagesc(sax,newrows');
% $$$     axis(sax,'tight');
% $$$     colormap(sax,'jet');
% $$$     sax.YTickLabel = {};
% $$$     caxis(sax,[-3,0.6]);        
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< IMSCALE field size stats CA3 ----------------------------------------------------
% function EgoProCode2D_f1_ego_field_size_stats_CA3()
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_ego_field_size_stats_ca3');
% $$$     
% $$$     rateRows = log10(egoSize(uidsCA3,:)*0.02^2);
% $$$     [mx,mind] = max(log10(egoSize(uidsCA3,:)*0.02^2),[],2);
% $$$     newrows=[];
% $$$     for m = 1:3 
% $$$         [~,sind] = sort(mx(mind==m));
% $$$         myrows = rateRows(mind==m,:);
% $$$         newrows = cat(1, newrows,myrows(sind,:));
% $$$     end
% $$$     
% $$$     imagesc(sax,newrows');
% $$$     axis(sax,'tight');
% $$$     colormap(sax,'jet');
% $$$     cax = colorbar();
% $$$     ylabel(cax,'log10(m^2)');
% $$$     sax.YTickLabel = {};
% $$$     caxis(sax,[-3,0.6]);
% $$$     
% $$$     
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%%%>>>
%--------------------------------------------------------------------------------------- TRANS
%%%<<< SCHEMATIC of allo to ego --------------------------------------------------------
% function EgoProCode2D_f1_subplot_allo_to_ego_schematic
    [pfig, sax] = setup_figure_('EgoProCode2D_f1_subplot_allo_to_ego_schematic');
    img = imread(fullfile( MTA_PATH,'analysis','EgoProCode2D','EgoProCode2D_figure_parts',...
                           'EgoProCode2D_f1_subplot_allo_to_ego_schematic.png'));
    image(sax,img);
    axis(sax,'ij');
    axis(sax,'tight');
    axis(sax,'off');
    daspect(sax,[1,1,1]);
    
    savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
    close(pfig);
%%%>>>
%--------------------------------------------------------------------------------------- TRANS    
%%%<<< LABEL descending ca1 ------------------------------------------------------------
% function EgoProCode2D_f1_descending_label_ca1    
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_descending_label_ca1');
% $$$     text(sax,0,0.5,'CA1','FontSize',8);
% $$$     sax.Color = 'w';
% $$$     axis(sax,'off');
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%%%>>>
%--------------------------------------------------------------------------------------- TRANS        
%%%<<< LABEL trough ca1 ----------------------------------------------------------------
% function EgoProCode2D_f1_trough_label_ca1    
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_trough_label_ca1');
% $$$     text(sax,0,0.5,'CA1','FontSize',8);
% $$$     sax.Color = 'w';
% $$$     axis(sax,'off');
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%%%>>>
%--------------------------------------------------------------------------------------- TRANS        
%%%<<< LABEL ascending ca1 -------------------------------------------------------------
% function EgoProCode2D_f1_ascending_label_ca1    
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_ascending_label_ca1');
% $$$     text(sax,0,0.5,'CA1','FontSize',8);
% $$$     sax.Color = 'w';
% $$$     axis(sax,'off');
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%%%>>>
%--------------------------------------------------------------------------------------- TRANS        
%%%<<< LABEL descending ca3 ------------------------------------------------------------
% function EgoProCode2D_f1_descending_label_ca3    
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_descending_label_ca3');
% $$$     text(sax,0,0.5,'CA3','FontSize',8);
% $$$     sax.Color = 'w';
% $$$     axis(sax,'off');
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%%%>>>
%--------------------------------------------------------------------------------------- TRANS        
%%%<<< LABEL trough ca3 ----------------------------------------------------------------
% function EgoProCode2D_f1_trough_label_ca3    
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_trough_label_ca3');
% $$$     text(sax,0,0.5,'CA3','FontSize',8);
% $$$     sax.Color = 'w';
% $$$     axis(sax,'off');
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%%%>>>
%--------------------------------------------------------------------------------------- TRANS        
%%%<<< LABEL ascending ca3 -------------------------------------------------------------
% function EgoProCode2D_f1_ascending_label_ca3    
% $$$     [pfig, sax] = setup_figure_('EgoProCode2D_f1_ascending_label_ca3');
% $$$     text(sax,0,0.5,'CA3','FontSize',8);
% $$$     sax.Color = 'w';
% $$$     axis(sax,'off');
% $$$     savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
% $$$     close(pfig);
%%%>>>
%--------------------------------------------------------------------------------------- TRANS        
    