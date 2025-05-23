% EgoProCode2D_f1_ego_field_size_stats_ca3

[pfig, sax] = setup_figure_( mfilename() );


rateRows = log10(egoSize(uidsCA3,:) .* 4 ./10000 );
[mx,mind] = max(log10(egoSize(uidsCA3,:)*0.02^2),[],2);
newrows=[];
EgoProCode2D_f1_ego_field_size_stats_ca3_rowlengths(m) = [];
for m = 1:3 
    [~,sind] = sort(mx(mind==m));
    myrows = rateRows(mind==m,:);
    newrows = cat(1, newrows,myrows(sind,:));
    EgoProCode2D_f1_ego_field_size_stats_ca3_rowlengths(m) = size(myrows,1);
end

imagesc(sax,newrows');
axis(sax,'tight');
colormap(sax,'jet');
cax = colorbar();
ylabel(cax,'log10(m^2)');
sax.YTickLabel = {};
caxis(sax,[-3,-0.1]);



savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

