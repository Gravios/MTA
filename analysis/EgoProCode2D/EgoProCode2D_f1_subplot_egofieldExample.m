% EgoProCode2D_f1_subplot_egofieldExample

[pfig, sax] = setup_figure_(mfilename());

% SUBPLOT <- ego-ratemap, rat model, circle
plot(pfe{exampleUnit.trialIndex},exampleUnit.id,1,'',[0,exampleUnit.maxRate],'colorMap',@jet,'mazeMaskFlag',false,'flipAxesFlag',true);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color,'LineWidth',2);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color,'LineWidth',2);
circle(0,0,200,'r-');

% FORMAT subplot
sax.XTickLabel =[];
sax.YTickLabel =[];
xlim(sax,[-300,300]);
ylim(sax,[-300,300]);
daspect(sax,[1,1,1]);

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));
close(pfig);
