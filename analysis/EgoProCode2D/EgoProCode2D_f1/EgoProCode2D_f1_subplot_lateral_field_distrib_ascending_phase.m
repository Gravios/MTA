% EgoProCode2D_f1_subplot_lateral_field_distrib_ascending_phase

[pfig, sax] = setup_figure_( mfilename() );

% PLOT distribution of lateral direction of rate map centers
[h,L,MX,MED,bw] = violin({egoMeanRmapPos(uidsCA3,3,2),...
                    egoMeanRmapPos(uidsCA1,3,2)});

% FORMAT subplot
delete(L);
af(@(hndl) set(hndl,'Vertices',fliplr(hndl.Vertices)), h);
af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr([3,3],:),ones([2,1]))');
hlns = findobj(sax,'Type','line');
for ln = 1:numel(hlns)
    [hlns(ln).XData,hlns(ln).YData] = deal(hlns(ln).YData,hlns(ln).XData);
end

xlim(sax,[-30,30]);
ylim(sax,[0.5,2.5]);

sax.XTickLabelMode = 'Auto';
sax.XTick = [-20,-10,0,10,20];

xlabel(sax,'cm');

sax.YTickLabel = {};
grid(sax,'on');

patch([-30,  -30,  30,  30],...
      [ 0.5, 1.5, 1.5, 0.5],...
      [  -1,  -1,  -1,  -1],...
      [ 0.5, 0.5, 0.5],...
      'FaceAlpha',0.2,...
      'EdgeAlpha',0);

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

