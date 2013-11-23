function pf = qpf(Trial,display,overwrite)

pf{1} = MTAPlaceField(Trial,[],{'walk'},overwrite);
pf{2} = MTAPlaceField(Trial,[],{'rear'},overwrite);
pf{3} = MTAPlaceField(Trial,[],{'theta'},overwrite);
pf{4} = MTAPlaceField(Trial,[],{'theta'},overwrite,[],'all','pfcrz');
pf{5} = MTAPlaceField(Trial,[],{'bturn'},overwrite);
pf{6} = MTAPlaceField(Trial,[],{'hturn'},overwrite);

if display
    figure,
    unit = 1;
    while unit~=-1
        subplot2(2,3,2,3)
        pf{6}.plot(unit,1)
        title(pf{6}.stateLabel)
        subplot2(2,3,1,3)
        pf{5}.plot(unit,1)
        title(pf{5}.stateLabel)
        subplot2(2,3,1,1)
        pf{4}.plot(unit,1)
        title(pf{4}.stateLabel)
        subplot2(2,3,1,2)
        pf{2}.plot(unit,1)
        title(pf{2}.stateLabel)
        subplot2(2,3,2,1)
        pf{3}.plot(unit,1)
        title(pf{3}.stateLabel)
        subplot2(2,3,2,2)
        pf{1}.plot(unit,1)
        title([pf{1}.stateLabel ' ' num2str(unit)])
        unit=figure_controls(gcf,unit);
    end
end