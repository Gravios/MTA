% EgoProCode2D_f1_subplot_ca3_field_size

[pfig, sax] = setup_figure_( mfilename() );

% SUBPLOT <- Asce, Desc field size vs Trough field size
plot((egoSize(uidsCA3,2))*0.02^2,(egoSize(uidsCA3,1,1))*0.02^2,'.','Color',pclr(1,:));
plot((egoSize(uidsCA3,2))*0.02^2,(egoSize(uidsCA3,3,1))*0.02^2,'.','Color',pclr(3,:));

% FORMAT subplot
line([0,0.3],[0,0.3],'Color','k')
grid(sax(end),'on');
title(sax(end),'Field Area CA3')
xlabel('Trough Field Size (m^2)');
ylabel('asce and desc Field Size (m^2)');

savefig(pfig, fullfile(partsPath, [sax.Tag,'.fig']));

if close_figure
    close(pfig);
end


