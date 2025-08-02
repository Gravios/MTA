% EgoProCode2D_f1_subplot_ca1_rate_ratios

[pfig, sax] = setup_figure_( mfilename() );

% SUBPLOT <- rate ratios
plot(egoMaxRmapRate(uidsCA1,1)./egoMaxRmapRate(uidsCA1,2),...
     egoMaxRmapRate(uidsCA1,3)./egoMaxRmapRate(uidsCA1,2),'.')

% FORMAT subplot
xlim(sax,[0,1.8]);
ylim(sax,[0,1.8]);
grid(sax,'on');
title(sax(end),'Rate Ratio CA1');
xlabel(sax(end),'Dsc/Trough (A.U.)');
ylabel(sax(end),'Asc/Trough (A.U.)');

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

