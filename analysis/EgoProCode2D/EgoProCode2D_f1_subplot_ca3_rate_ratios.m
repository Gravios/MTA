% EgoProCode2D_f1_subplot_ca3_rate_ratios

[pfig, sax] = setup_figure_( mfilename() );

% SUBPLOT <- rate ratios
plot(egoMaxRmapRate(uidsCA3,1)./egoMaxRmapRate(uidsCA3,2),...
     egoMaxRmapRate(uidsCA3,3)./egoMaxRmapRate(uidsCA3,2),'.')

% FORMAT subplot
xlim(sax,[0,1.8]);
ylim(sax,[0,1.8]);
grid(sax,'on');
title(sax(end),'Rate Ratio CA3');
xlabel(sax(end),'Dsc/Trough (A.U.)');
ylabel(sax(end),'Asc/Trough (A.U.)');

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end

