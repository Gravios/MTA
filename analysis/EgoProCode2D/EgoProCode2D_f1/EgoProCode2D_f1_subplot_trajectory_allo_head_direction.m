%EgoProCode2D_f1_subplot_trajectory_allo_head_direction

[pfig, sax] = setup_figure_(mfilename());

% PLOT allo-ratemap, trajectory (color=thetaPhase), and skeletons
plot(pft{exampleUnit.trialIndex},units{exampleUnit.trialIndex}(exampleUnit.index),1,'text',[0,mrate],'colorMap',@bone,'nanColor',[1,1,1]);
scatter(xyz{exampleUnit.trialIndex}(exampleUnit.trajectoryTimeSeries,'hcom',1),...
        xyz{exampleUnit.trialIndex}(exampleUnit.trajectoryTimeSeries,'hcom',2),...
        4,...
        atan2(pfstrj(exampleUnit.trajectoryTimeSeries,2),pfstrj(exampleUnit.trajectoryTimeSeries,1)),'filled')
circle(mxp(1),mxp(2),200,'-r');
plotSkeleton(Trials{exampleUnit.trialIndex},xyz{exampleUnit.trialIndex},exampleUnit.trajectoryTimeSeries(1));        
plotSkeleton(Trials{exampleUnit.trialIndex},xyz{exampleUnit.trialIndex},exampleUnit.trajectoryTimeSeries(round(length(exampleUnit.trajectoryTimeSeries)/1.15)));    

% FORMAT subplot
xlim(sax,exampleUnit.close.Xlims); 
ylim(sax,exampleUnit.close.Ylims); 
sax.XTickLabel =[];
sax.YTickLabel =[];
colormap(sax,'hsv');    
caxis([-pi,pi]);        
daspect([1,1,1]);    
grid('on');    

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
close(pfig);
