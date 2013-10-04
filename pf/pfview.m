function pfview(Trial,Pfs,varargin)
[name,display_mode] = DefaultArgs(varargin,{'all','report'});


figure(102)
set(gcf,'CurrentCharacter','l');
unit = 1;
while 1,
    clf
    if ~iscell(Pfs),
        Pfs.plot(unit)
        title([Pfs.stateLabel ' ' num2str(unit)])
        numClu = length(Pfs.rateMap);
    else
        for pf = 1:length(Pfs),
            %% Rate Map
            subplotfit(pf,length(Pfs));
            Pfs{pf}.plot(unit)
            title([Pfs{pf}.stateLabel ' ' num2str(unit)])
        end
        numClu = length(Pfs{1}.rateMap);
    end


    %% START - Controls
    switch display_mode
      case 'report'
        %% ReportFig
        set(gcf,'CurrentCharacter','n');
        reportfig(102,[Trial.filebase '.pfs.' name],[],['unit: ' num2str(unit)],[],0);
        unit = unit+1;
        if unit > numClu, return,end
      case 'display'
        %% manual controls
        unit = figure_controls(gcf,unit);
        if unit==-1;return,end
    end
    %% END - Controls
    
    if ~iscell(Pfs)
        if unit>size(Pfs.cluMap,1), unit = 1; end
    else 
        if unit>size(Pfs{1}.cluMap,1), unit = 1; end
    end
end
