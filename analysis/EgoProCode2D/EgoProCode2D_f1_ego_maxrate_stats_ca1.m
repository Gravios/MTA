% EgoProCode2D_f1_ego_maxrate_stats_ca1

[pfig, sax] = setup_figure_( mfilename() );

imagesc(sax,EgoProCode2D_f1_ego_maxrate_stats_ca1_rows')
axis(sax,'tight');
colormap(sax,'jet');    
sax.YTickLabel = {};
sax.XTickLabel = {};
caxis(sax,[0,1.5]);

%title(sax,'NEW TITLE');

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

