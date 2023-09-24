% EgoProCode2D_f1_subplot_lateral_distrib_stats

[pfig, ax] = setup_figure_( mfilename() );

tag = ax.Tag;
delete(ax);
ax1 = axes();
ax1.Tag = tag;
axPos = ax1.Position;
im1 = imagesc(ax1,egoMeanPosMeanLateral);
colormap(ax1,'jet');
cax1 = colorbar(ax1,'Location','eastoutside');
mask = triu(ones(size(egoMeanPosPvalLateral))) > 0;
im1.AlphaData = ones(size(egoMeanPosPvalLateral)).*mask;
ax1.Position = axPos;
ax1.Color = 'none';
ax1.Visible = 'off';
ax1.Tag = tag;
caxis(ax1,[-10,10]);
ylabel(cax1,'cm');

ax2 = axes();
im2 = imagesc(ax2,egoMeanPosPvalLateral);
mask = tril(ones(size(egoMeanPosPvalLateral))) > 0;
im2.AlphaData = ones(size(egoMeanPosPvalLateral)).*mask;
colormap(ax2,'summer');
cax2 = colorbar(ax2,'Location','southoutside');
ax2.Color = 'none';
ax2.Visible = 'off';
ax2.Position = axPos;
ax2.Tag = tag;
caxis(ax2,[-30,-1.301]);
ylabel(cax2,'P-val');

ax3 = axes();
im3 = imagesc(ax3,tril(ones(size(egoMeanPosPvalLateral))).*isnan(egoMeanPosPvalLateral) .*0.8);
mask = tril(ones(size(egoMeanPosPvalLateral))).*isnan(egoMeanPosPvalLateral) .*0.8 ;
im3.AlphaData = ones(size(egoMeanPosPvalLateral)).*mask;
colormap(ax3,'bone');
ax3.Color = 'none';
ax3.Visible = 'off';
ax3.Position = axPos;
ax3.Tag = tag;
caxis(ax3,[0,1]);

for l = 1:5,
    line([0.5,6.5],[0.5,0.5]+l,'Color','k','LineWidth',1);
    line([0.5,0.5]+l,[0.5,6.5],'Color','k','LineWidth',1);    
end

%ylabel(ax1,'Mean Lat Diff');

daspect(ax1,[1,1,1]);
daspect(ax2,[1,1,1]);
daspect(ax3,[1,1,1]);

%title(ax1,'Lateral Shift Difference');

savefig(pfig, fullfile(partsPath, [ax1.Tag,'.fig']));

if close_figure
    close(pfig);
end
