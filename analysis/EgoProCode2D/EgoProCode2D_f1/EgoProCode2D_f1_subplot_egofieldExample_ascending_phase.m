% EgoProCode2D_f1_subplot_egofieldExample_ascending_phase

[pfig, sax] = setup_figure_( mfilename() );

% PLOT egocentric rate map
plot(pfet{exampleUnit.trialIndex}{3},...
     exampleUnit.id,...
     1, ...
     '',...
     [0,exampleUnit.maxRate],...
     false,...
     'colorMap',@jet,...
     'flipAxesFlag',true);

% FORMAT subplot
axis(sax,'tight');
xlim(sax,exampleUnit.ego.xlims);
ylim(sax,exampleUnit.ego.ylims);
Lines([],0,'w');
Lines(0,[],'w');
sax.XTickLabel =[];
sax.YTickLabel =[];

% ADD rat model to subplot
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false,[],[]);
subject = update_subject_patch(subject,'body',[],false,[],[]);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.25);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

