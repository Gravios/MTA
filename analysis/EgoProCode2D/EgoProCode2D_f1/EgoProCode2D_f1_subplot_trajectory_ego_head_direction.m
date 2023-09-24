% EgoProCode2D_f1_subplot_trajectory_ego_head_direction

[pfig, sax] = setup_figure_(mfilename());

% PLOT ego-ratemap, 
plot(pfe{exampleUnit.trialIndex},units{exampleUnit.trialIndex}(exampleUnit.index),1,'',[0,mrate],'mazeMaskFlag',false,'colorMap',@bone,'flipAxesFlag',true);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color,'LineWidth',2);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color,'LineWidth',2);
circle(0,0,200,'r-');
scatter(pfstrj(exampleUnit.trajectoryTimeSeries,2),                                                        ...
        pfstrj(exampleUnit.trajectoryTimeSeries,1),                                                        ...
        4,                                                                                                 ...
        atan2(pfstrj(exampleUnit.trajectoryTimeSeries,2),pfstrj(exampleUnit.trajectoryTimeSeries,1)),'filled')

% FORMAT subplot
sax.XTickLabel =[];
sax.YTickLabel =[];
xlim    (sax,[-300,300]);    
ylim    (sax,[-300,300]);
colormap(sax,'hsv');
caxis   (sax,[-pi,pi]);            
daspect (sax,[1,1,1]);
grid    (sax,'on');

savefig(pfig, fullfile(partsPath,[sax.Tag, '.fig']));
close(pfig);
