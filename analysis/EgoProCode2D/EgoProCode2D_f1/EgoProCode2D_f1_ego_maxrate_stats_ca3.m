% EgoProCode2D_f1_ego_maxrate_stats_ca3

[pfig, sax] = setup_figure_( mfilename() );


imagesc(sax,EgoProCode2D_f1_ego_maxrate_stats_CA3_rows');
axis(sax,'tight');    
colormap(sax,'jet');
cax = colorbar();
ylabel(cax,'log10(Hz)');
sax.YTickLabel = {};
sax.XTickLabel = {};
caxis(sax,[0,1.5]);    

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

