function gen_ccg_figs(Trial,varargin)
[ccg_name,trialName,mazeName,display_mode] = DefaultArgs(varargin,{'rear_par10','all','cof','report'});


if ~isa(Trial,'MTATrial'),
    Trial = MTATrial(Trial,{},trialName,[],[],mazeName);
end


load([Trial.spath.analysis Trial.filebase '.ccg.' ccg_name '.mat']);

fccg = Bccg.filter(gausswin(7));

numClu = size(Bccg.cluMap,1);

[accg,atbin] = autoccg(Trial.name);
clear('Session');



unit = 1;
figure(2993),
set(gcf,'CurrentCharacter','l');
unit = 1;
while 1,
    clf
    if Bccg.partitions>1,
        subplot2(1,3,1,2);imagesc(Bccg.tbin,0:Bccg.partitions:1000,sq(fccg(:,unit,1,:,1))'),axis xy,
        title('CCG centered on rear onset')

        subplot2(1,3,1,3);imagesc(Bccg.tbin,0:Bccg.partitions:1000,sq(fccg(:,unit,2,:,1))'),axis xy
        title('CCG centered on rear offset')
    else
        subplot2(1,3,1,2);bar(Bccg.tbin,sq(fccg(:,unit,1,:,1))),axis xy,axis tight
        title('CCG centered on rear onset')

        subplot2(1,3,1,3);bar(Bccg.tbin,sq(fccg(:,unit,2,:,1))),axis tight,axis xy
        title('CCG centered on rear offset')
    end

    subplot2(1,3,1,1);
    bar(atbin,accg(:,unit))
    title(['unit: ' num2str(unit)]);
    axis tight


    %% START - Controls
    switch display_mode

      case 'report'
        %% ReportFig
        set(gcf,'CurrentCharacter','n');
        reportfig(2993,[Trial.filebase '.ccg.' ccg_name],[],['unit: ' num2str(unit)],[],0);
        unit = unit+1;
        if unit > numClu, return,end

      case 'display'
        %% manual controls
        unit = figure_controls(unit);
        if unit==-1;return,end

    end
    %% END - Controls

end

