% EgoProCode2D_f1_subplot_forward_distrib_stats

[pfig, ax] = setup_figure_( mfilename() );

tag = ax.Tag;
delete(ax);
ax1 = axes();
axPos = ax1.Position;
im1 = imagesc(ax1,egoMeanPosMeanForward);
colormap(ax1,'jet');
cax1 = colorbar(ax1,'Location','eastoutside');
mask = triu(ones(size(egoMeanPosPvalForward))) > 0;
im1.AlphaData = ones(size(egoMeanPosPvalForward)).*mask;
ax1.Position = axPos;
ax1.Color = 'none';
ax1.Visible = 'off';
ax1.Tag = tag;
caxis(ax1,[-10,10]);
ylabel(cax1,'cm');

ax2 = axes();
im2 = imagesc(ax2,egoMeanPosPvalForward);
mask = tril(ones(size(egoMeanPosPvalForward))) > 0;
im2.AlphaData = ones(size(egoMeanPosPvalForward)).*mask;
colormap(ax2,'summer');
%cax2 = colorbar(ax2,'Location','southoutside');
ax2.Color = 'none';
ax2.Visible = 'off';
ax2.Position = axPos;
ax2.Tag = tag;
caxis(ax2,[-30,-1.301]);
%ylabel(cax2,'P-val');


ax3 = axes();
im3 = imagesc(ax3,tril(ones(size(egoMeanPosPvalForward))).*isnan(egoMeanPosPvalForward) .*0.8);
mask = tril(ones(size(egoMeanPosPvalForward))).*isnan(egoMeanPosPvalForward) .*0.8 ;
im3.AlphaData = ones(size(egoMeanPosPvalForward)).*mask;
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


daspect(ax1,[1,1,1]);
daspect(ax2,[1,1,1]);
daspect(ax3,[1,1,1]);

ylabel(ax1,'Mean AP Diff');


drawnow();
pause(1);

%title(ax1,'Anteroposterior Shift Difference');

savefig(pfig, fullfile(partsPath, [ax1.Tag,'.fig']));

if close_figure
    close(pfig);
end
