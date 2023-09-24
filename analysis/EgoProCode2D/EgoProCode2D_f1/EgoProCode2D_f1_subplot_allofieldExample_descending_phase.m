% EgoProCode2D_f1_subplot_allofieldExample_descending_phase

[pfig, sax] = setup_figure_( mfilename() );

plot(ratemapsAlloThp{exampleUnit.trialIndex}{1},...
     exampleUnit.id,...
     1,...
     '',...
     [0,exampleUnit.maxRate],...
     false,...
     'colorMap',@jet);

% FORMAT subplot
axis(sax,'tight');

%xlim(sax,[-300,300]);
%ylim(sax,[-300,300]);
%circle(mxp(1),mxp(2),200,'-r');            
xlim(sax,exampleUnit.allo.xlims+mxp(1));
ylim(sax,exampleUnit.allo.ylims+mxp(2));

ylabel('Allocentric');

%Lines([],mxp(2),'w');
%Lines(mxp(1),[],'w');
sax.XTickLabel =[];
sax.YTickLabel =[];


savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

