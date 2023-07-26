
%EgoProCode2D_f1_subplot_placefieldExample


[pfig, sax] = setup_figure_(mfilename());

% SUBPLOT <- placefield, scalebars, circle
plot(pft{exampleUnit.trialIndex},exampleUnit.id,1,'',[0,exampleUnit.maxRate],'colorMap',@jet);
sax.XTickLabel =[];
sax.YTickLabel =[];
line(sax,                               ...
     [-470,-270],                       ...
     [460].*[1,1],'LineWidth',2,'Color','w');
line(sax,                            ...
     [-460].*[1,1],                  ...
     [470,270],'LineWidth',2,'Color','w');
circle(mxp(1),mxp(2),200,'-r');        

savefig(pfig, fullfile(partsPath,[sax.Tag,'.fig']));
close(pfig);
