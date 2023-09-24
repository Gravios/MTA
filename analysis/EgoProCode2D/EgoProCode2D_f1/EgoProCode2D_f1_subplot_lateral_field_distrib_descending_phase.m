% EgoProCode2D_f1_subplot_lateral_field_distrib_descending_phase

[pfig, sax] = setup_figure_( mfilename() );

% PLOT distribution of lateral direction of rate map centers
[h,L,MX,MED,bw] = violin({egoMeanRmapPos(uidsCA3,1,2),...
                    egoMeanRmapPos(uidsCA1,1,2)});

% FORMAT subplot
delete(L);
af(@(hndl) set(hndl,'Vertices',fliplr(hndl.Vertices)), h);
af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr([1,1],:),ones([2,1]))');
hlns = findobj(sax,'Type','line');
for ln = 1:numel(hlns)
    [hlns(ln).XData,hlns(ln).YData] = deal(hlns(ln).YData,hlns(ln).XData);
end

ylim(sax,[0.5,2.5]);
xlim(sax,[-30,30]);

sax.XTickLabelMode = 'Auto';
sax.XTick = [-20,-10,0,10,20];

% xticks at 1 and 2


xlabel(sax,'cm');
ylabel(sax,'Lat');
sax.YTickLabel = {};    

patch([-30,  -30,  30,  30],...
      [ 0.5, 1.5, 1.5, 0.5],...
      [  -1,  -1,  -1,  -1],...
      [ 0.5, 0.5, 0.5],...
      'FaceAlpha',0.2,...
      'EdgeAlpha',0);

text( sax,...
      -20,...
      1,...
      'CA3',...
      'HorizontalAlignment','center');
text( sax,...
      -20,...
      2,...
      'CA1',...
      'HorizontalAlignment','center');
      


grid(sax,'on');

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

