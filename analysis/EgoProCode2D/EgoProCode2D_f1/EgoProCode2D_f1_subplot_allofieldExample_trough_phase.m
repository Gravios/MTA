% EgoProCode2D_f1_subplot_allofieldExample_trough_phase

[pfig, sax] = setup_figure_( mfilename() );


plot(ratemapsAlloThp{exampleUnit.trialIndex}{2},...
     exampleUnit.id,...
     1,...
     '',...
     [0,exampleUnit.maxRate],...
     false,...
     'colorMap',@jet);


% FORMAT subplot
axis(sax,'tight');
xlim(sax,exampleUnit.allo.xlims+mxp(1));
ylim(sax,exampleUnit.allo.ylims+mxp(2));
sax.XTickLabel =[];
sax.YTickLabel =[];

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end


