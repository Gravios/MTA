% EgoProCode2D_f1_subplot_allo_to_ego_schematic

[pfig, sax] = setup_figure_( mfilename() );
img = imread(fullfile( MTA_PATH,'analysis','EgoProCode2D','EgoProCode2D_figure_parts',...
                       [ mfilename(),'.png'] ));
image(sax,img);
axis(sax,'ij');
axis(sax,'tight');
axis(sax,'off');
daspect(sax,[1,1,1]);

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

