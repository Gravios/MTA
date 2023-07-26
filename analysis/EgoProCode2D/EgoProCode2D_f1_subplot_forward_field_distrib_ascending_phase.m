% EgoProCode2D_f1_subplot_forward_field_distrib_ascending_phase

[pfig, sax] = setup_figure_( mfilename() );

% PLOT distribution of forward direction of rate map centers
hold(sax,'on');
[h,L,MX,MED,bw] = violin({(egoMeanRmapPos(uidsCA1,3,1)),...
                    egoMeanRmapPos(uidsCA3,3,1)});

% FORMAT subplot
delete(L);
af(@(hndl,c) set(hndl,'FaceColor',c{1}), h,mat2cell(pclr([3,3],:),ones([2,1]))');
xlim([0.5,2.5])
ylim(sax,[-25,35]);    

sax.YAxisLocation = 'right';
sax.YTickLabelMode = 'Auto';
sax.YTick = [-20,-10,0,10,20,30];

sax.XTickLabel = {};
grid(sax,'on');

ylabel(sax,'cm');

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

