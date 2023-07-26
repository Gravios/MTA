% EgoProCode2D_f1_descending_label_ca1

[pfig, sax] = setup_figure_( mfilename() );
text(sax,0,0.5,'CA1','FontSize',8);
sax.Color = 'w';
axis(sax,'off');
savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

