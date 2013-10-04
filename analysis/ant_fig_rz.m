function ant_fig_rz(Trial,varargin)
[event_type,trialName,mazeName,display_mode] = DefaultArgs(varargin,{'onset','all','cof','report'});

if ~isa(Trial,'MTATrial'),
    Trial = MTATrial(Trial,{},trialName,[],[],mazeName);
end

load([Trial.spath.analysis Trial.filebase '.ccg.rear_par1.mat']);

fccg = Bccg.filter(gausswin(5));

load([Trial.spath.analysis Trial.filebase '.ccg.rear_par10.mat']);

dccg = Bccg.filter(gausswin(5));

numClu = size(Bccg.cluMap,1);

[accg,atbin] = autoccg(Trial.name);
clear('Session');

pf_search = MTAPlaceField([]);
pf_search.type = 'xy';
pf_search.mazeName = 'cof';
pf_search.trialName = 'all';
pf_search.trackingMarker = 'head_front';
pf_search.stateLabel = 'theta';
pf_search.spk_shuffle = 'n';
pf_search.pos_shuffle = 0;
pf_search.numBSiterations = 1;
pf_search.nbins = 50;
pf_search.smooth = 0.03;

Trial = Trial.load_Pfs();

unit = 1;
figure(2993),
set(gcf,'CurrentCharacter','l');
unit = 1;
x = 3;
y = 4;

while 1,
    clf

    subplot2(y,x,[1,2],1); 
    pf_search.type = 'xy';
    pf_search.stateLabel = 'head.theta';
    Pfs = Trial.getPfs(pf_search);
    Pfs.plot(unit); hold on,
    plot(Trial.xyz(Trial.Bhv.getState('rear').state(:,1),Trial.Model.gmi('head_front'),1),...
         Trial.xyz(Trial.Bhv.getState('rear').state(:,1),Trial.Model.gmi('head_front'),2),...
         '.','Color','m')
    title('Place Field: head.theta - xy')

    subplot2(y,x,[1,2],2); 
    pf_search.type = 'pfcrz';
    pf_search.stateLabel = 'theta';
    Pfs = Trial.getPfs(pf_search);
    Pfs.plot(unit);
    title('Place Field: theta - rz')

    subplot2(y,x,3,2);bar(Bccg.tbin,sq(fccg(:,unit,1,:,1))),axis xy,axis tight
    title('CCG centered on rear onset')

    subplot2(y,x,4,2);bar(Bccg.tbin,sq(fccg(:,unit,2,:,1))),axis tight,axis xy
    title('CCG centered on rear offset')

    subplot2(y,x,3,3);imagesc(Bccg.tbin,0:Bccg.partitions:1000,sq(dccg(:,unit,1,:,1))'),axis xy,
    title('CCG centered on rear onset')

    subplot2(y,x,4,3);imagesc(Bccg.tbin,0:Bccg.partitions:1000,sq(dccg(:,unit,2,:,1))'),axis xy
    title('CCG centered on rear offset')

    % Place Field: Rearing 
    subplot2(y,x,[1,2],3); 
    pf_search.type = 'xy';
    pf_search.stateLabel = 'pre_peri_rear_onset';
    Pf1 = Trial.getPfs(pf_search);
    pf_search.stateLabel = 'post_peri_rear_onset';
    Pf2 = Trial.getPfs(pf_search);
    ppf(Pf1.xbin,Pf1.ybin,log10((Pf2.rateMap{unit}+1)./(Pf1.rateMap{unit}+1)))
    colorbar
    title({'log10 ratio between ratemaps around rearing onset','(0ms to 750ms) and (-750ms to 0ms)'});

    subplot2(y,x,[3,4],1);
    bar(atbin,accg(:,unit))
    title(['u: ' num2str(unit) ' e: ' num2str(Pfs.cluMap(unit,2)) ' c: ' num2str(Pfs.cluMap(unit,3))]);
    axis tight

    %% START - Controls
    switch display_mode

      case 'report'
        %% ReportFig
        set(gcf,'CurrentCharacter','n');
        reportfig(2993,[Trial.filebase '.ccg.ant_fig_rz'],[],['unit: ' num2str(unit)],[],0);
        print(2993,'-depsc2','-r200',['/data/homes/gravio/mrep/ant_figs/' ...
                            Trial.name '\' Trial.name '.ant_fig_rz.' ...
                            num2str(unit) '.eps'])
        print(2993,'-dpng',['/data/homes/gravio/mrep/ant_figs/' ...
                            Trial.name '\' Trial.name '.ant_fig_rz.' ...
                            num2str(unit) '.png'])

        unit = unit+1;
        if unit > numClu, return,end

      case 'display'
        %% manual controls
        unit = figure_controls(unit);
        if unit==-1;return,end

    end
    %% END - Controls

end
