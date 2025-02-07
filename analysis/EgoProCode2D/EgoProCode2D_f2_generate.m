% shuffling
% premutation test
%             phz : theta phase
%             hba : head-body-angle ( 0=straight, pi/2=left/right, pi=broken rat)
%                     {HBA-L}   {HBA-C}   {HBA-R}
%         allofield  egofield  egofield  egofield  plot(HREF-LAT(LxR))  ksden(HREF-LAT(LCR))
% {PHZ-A} egofield      X         X         X      plot(HREF-LAT(LxR))  ksden(HREF-LAT(LCR))
% {PHZ-T} egofield      X         X         X      plot(HREF-LAT(LxR))  ksden(HREF-LAT(LCR))
% (PHZ-D) egofield      X         X         X      plot(HREF-LAT(LxR))  ksden(HREF-LAT(LCR))



%%%<<< LOADDATA -----------------------------------------------------------------
%%%<<< ASSETS and paths----------------------------------------------------------
global MTA_PROJECT_PATH
EgoProCode2D_load_data();
EgoProCode2D_f2_data_egoHba();
EgoProCode2D_f2_data_egoHbaPhz();
EgoProCode2D_f2_data_decoding();
EgoProCode2D_f2_load_subject_model();

partsPath = fullfile( ...
    fullfile( MTA_PROJECT_PATH,                                               ...
              'analysis',                                                     ...
              'EgoProCode2D',                                                 ...
              'EgoProCode2D_figure_parts'));

overwrite = false;
rat = load_patch_model('rat');
mask = double(sqrt(bsxfun(@plus,egoHbaRmaps.xbins.^2,egoHbaRmaps.ybins'.^2)') < 445);
mask(~mask) = nan;
regions = {'ca1','ca3'};
region = 1;
pc = bins.phz.count;
hc = bins.hba.count;
vc = bins.hav.count;
%%%>>>---------------------------------------------------------------------------
%%%<<< CA1 example unit ---------------------------------------------------------

exampleUnit.trialIndex = 20;
exampleUnit.close.Xlims = [-200,400];
exampleUnit.close.Ylims = [-400,200];
exampleUnit.id = 25;
exampleUnit.maxRate = 18;
exampleUnit.colormap = @jet;
exampleUnit.index = find(unitsEgo{exampleUnit.trialIndex}==exampleUnit.id);
exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];
uids = unitsEgoCA1;
gXoffset = 0;
gYoffset = 0;

%%%>>>---------------------------------------------------------------------------
%%%<<< CA3 example unit ---------------------------------------------------------

% $$$         % CA3
% $$$         exampleUnit.trialIndex = 6;
% $$$         exampleUnit.close.Xlims = [-200,400];
% $$$         exampleUnit.close.Ylims = [-400,200];
% $$$         exampleUnit.id = 10;
% $$$         exampleUnit.maxRate = 18;
% $$$         %exampleUnit.index = find(units{exampleUnit.trialIndex}==exampleUnit.id);
% $$$         exampleUnit.trajectoryTimeSeries = [274*250]:[278.5*250];
% $$$         uids = unitsEgoCA3;
% $$$         gXoffset = 0;
% $$$         gYoffset = -2.4*5;

%%%>>>---------------------------------------------------------------------------
%%%<<< DATA THETA Subject -------------------------------------------------------

% $$$ decoded = EgoProCode2D_comp_decoded(dca, 'theta', 1, mfun, beta);
% $$$ ffun = @(x,xdata) x.*sin(xdata);
% $$$ ind = within_ranges(decoded.phz,[4,5.5]);
% $$$ x = lsqcurvefit(ffun, 1, decoded.ang(ind)', decoded.hba(ind)');
%figure,plot(decoded.ang(ind), decoded.hba(ind),'.')

%%%>>>
%%%<<< DATA THETA CONDITIONAL COPULA --------------------------------------------

decoded = EgoProCode2D_comp_decoded(dca, 'theta', 'all', mfun, beta);

[decoded.hbaUni,decoded.hbaUx] = MakeUniformDistr(decoded.hba);
[decoded.havUni,decoded.havUx] = MakeUniformDistr(decoded.hav);

mBinHav.edges = linspace(-0.4, 0.4, 10);
mBinHav.centers = mean([mBinHav.edges(1:end-1); mBinHav.edges(2:end)]);
mBinHav.count = numel(mBinHav.edges)-1;
mBinHba.edges = linspace(-3, 3, 10);
mBinHba.centers = mean([mBinHba.edges(1:end-1); mBinHba.edges(2:end)]);
mBinHba.count = numel(mBinHba.edges)-1;
mout  = zeros([mBinHba.count, mBinHav.count, pc]);
sout  = zeros([mBinHba.count, mBinHav.count, pc]);
tout  = zeros([mBinHba.count, mBinHav.count, pc]);
cout  = zeros([mBinHba.count, mBinHav.count, pc]);
fout  = zeros([mBinHba.count, mBinHav.count, pc]);
fouts = zeros([mBinHba.count, mBinHav.count, pc]);
rout  = zeros([mBinHba.count, mBinHav.count, pc]);
aout  = zeros([mBinHba.count, mBinHav.count, pc]);
for phzI = 1 : pc
    for hbaI = 1 : mBinHba.count
        for havI = 1 : mBinHav.count
            ind =  WithinRanges(decoded.phz, bins.phz.edges([phzI, phzI+1])) ...
                 & WithinRanges(-decoded.havUni, mBinHav.edges([havI, havI+1])) ...
                 & WithinRanges(decoded.hbaUni,  mBinHba.edges([hbaI, hbaI+1])) ...
                 & randn(size(decoded.hba))>0;
            [aout(hbaI, havI, phzI),rout(hbaI, havI, phzI)] = ...
                cart2pol(mean(decoded.fwd(ind),'omitnan'),    ...
                         mean(decoded.clat(ind),'omitnan'));
            mout(hbaI,  havI, phzI) = mean(decoded.clat(ind),'omitnan');
            sout(hbaI,  havI, phzI) =  std(decoded.clat(ind),'omitnan');
            fout(hbaI,  havI, phzI) = mean( decoded.fwd(ind),'omitnan');
            fouts(hbaI, havI, phzI) = std( decoded.fwd(ind),'omitnan');
            cout(hbaI,  havI, phzI) = sum(ind);
        end
    end
end
dmask = ones([mBinHba.count, mBinHav.count]);

%%%>>>
%%%<<< DATA THETA ---------------------------------------------------------------

decoded = EgoProCode2D_comp_decoded(dca, 'theta', 'CA1', mfun, beta);
mBinHav.edges = linspace(-0.3, 0.3, 10);
mBinHav.centers = mean([mBinHav.edges(1:end-1); mBinHav.edges(2:end)]);
mBinHav.count = numel(mBinHav.edges)-1;
mBinHba.edges = linspace(-1.2, 1.2, 10);
mBinHba.centers = mean([mBinHba.edges(1:end-1); mBinHba.edges(2:end)]);
mBinHba.count = numel(mBinHba.edges)-1;
mout  = zeros([mBinHba.count, mBinHav.count, pc]);
sout  = zeros([mBinHba.count, mBinHav.count, pc]);
tout  = zeros([mBinHba.count, mBinHav.count, pc]);
cout  = zeros([mBinHba.count, mBinHav.count, pc]);
fout  = zeros([mBinHba.count, mBinHav.count, pc]);
fouts = zeros([mBinHba.count, mBinHav.count, pc]);
rout  = zeros([mBinHba.count, mBinHav.count, pc]);
aout  = zeros([mBinHba.count, mBinHav.count, pc]);
for phzI = 1 : pc
    for hbaI = 1 : mBinHba.count
        for havI = 1 : mBinHav.count
            ind =  WithinRanges(decoded.phz, bins.phz.edges([phzI, phzI+1])) ...
                 & WithinRanges(decoded.hav, mBinHav.edges([havI, havI+1])) ...
                 & WithinRanges(decoded.hba,  mBinHba.edges([hbaI, hbaI+1])) ...
                 & randn(size(decoded.hba))>0;
            [aout(hbaI, havI, phzI),rout(hbaI, havI, phzI)] = ...
                cart2pol(mean(decoded.fwd(ind),'omitnan'),    ...
                         mean(decoded.clat(ind),'omitnan'));
            mout(hbaI,  havI, phzI) = mean(decoded.clat(ind),'omitnan');
            sout(hbaI,  havI, phzI) =  std(decoded.clat(ind),'omitnan');
            fout(hbaI,  havI, phzI) = mean( decoded.fwd(ind),'omitnan');
            fouts(hbaI, havI, phzI) = std( decoded.fwd(ind),'omitnan');
            cout(hbaI,  havI, phzI) = sum(ind);
        end
    end
end
dmask = ones([mBinHba.count, mBinHav.count]);
dmask(cout(:,:,1) < 32) = nan;


decoded_t = EgoProCode2D_comp_decoded(dca, 'theta', 'all', mfun, beta);
mout_t  = zeros([mBinHba.count, mBinHav.count]);
cout_t  = zeros([mBinHba.count, mBinHav.count]);
rout_t  = zeros([mBinHba.count, mBinHav.count]);
aout_t  = zeros([mBinHba.count, mBinHav.count]);
for a = 1:mBinHba.count
    for v = 1:mBinHav.count
        ind =     WithinRanges( decoded_t.hav, mBinHav.edges([v,v+1])) ...
                & WithinRanges( decoded_t.hba, mBinHba.edges([a,a+1])) ...              
                & randn(size(decoded_t.hba))>0;
        mout_t(a,v) = mean(decoded_t.clat(ind),'omitnan');
        [aout_t(a,v), rout_t(a,v)] = ...
            cart2pol(mean(decoded_t.fwd(ind), 'omitnan'),...
                     mean(decoded_t.clat(ind),'omitnan'));
        cout_t(a,v) = sum(ind);
    end
end
tmask = ones([mBinHba.count, mBinHav.count]);
tmask(cout_t < 32) = nan;

%%%>>>---------------------------------------------------------------------------
%%%>>>---------------------------------------------------------------------------


%%% FIGURE **********************************************************************
% BEGINFIG **********************************************************************

%%%<<< SETUP FIGURE -------------------------------------------------------------

% CONFIG 
hfig = figure(666002);
figureFormat      = 'A4';
figureOrientation = 'landscape';
figureUnits       = 'centimeters';
subplotHeight     = 1.4;%cm
subplotWidth      = 1.4;%cm
subplotPadVert    = 0.1;%cm
subplotPadHorz    = 0.1;%cm
% SETUP 
setup_figure(hfig,                  ... 
             figureFormat,          ... 
             figureOrientation,     ... 
             figureUnits,           ...
             subplotWidth,          ...
             subplotHeight,         ...
             subplotPadHorz,        ...
             subplotPadVert         ...
             );

%%%>>>---------------------------------------------------------------------------

%%%<<< ALLO/EGO EXAMPLE FIELD ---------------------------------------------------
% DESCRIPTION 
% <description>
% Egocentric fields, data partitioned by head-body-angle (HBA, degrees): Left
% green,(-70, -10); Center: blue (-10, 10); Right: red, (10, 70), and
% theta-phase (TP, degrees): Descending: cyan (30, 130); Trough: purple, (130,
% 230); Ascending: magenta, (230, 330).
% </description>
%%%<<< PLACE FIELD --THETA-- allocentric ratemap --------------------------------

% DESCRIPTION 
% <description>
% The firing rate of an example neuron conditioned on allo-centric position
% (room coordinates), and periods of attentive-prone behavior.</description>
% COORDINATES 
[yindex, yoffset] = deal( 1,  0 );
[xindex, xoffset] = deal( 1,  0 );
[   gyo, gxo    ] = deal( 0,  0 );
[yscale, xscale ] = deal( 1,  1 );
% SETUP 
sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% DATA 
etrial = exampleUnit.trialIndex;
eunit  = exampleUnit.id;
emax   = exampleUnit.maxRate;
ecmap  = exampleUnit.colormap;
% PLOT 
plot(pft{etrial}, eunit, 1, '', [0, emax], 'colorMap', ecmap);
% FORMAT 
title(sax(end), {'Place', 'Field'});
sax(end).XTickLabel =[];
sax(end).YTickLabel =[];

%--------------------------------------------------------------------------------

%%%>>>    
%%%<<< EGO FIELD   --THETA-- egocentric ratemap ---------------------------------
% DESCRIPTION 
% <description>
% The firing rate of an example neuron conditioned on ego-centric position
% (head coordinates), and periods of attentive-prone behavior.<\description>
% COORDINATES 
[yindex, yoffset] = deal( 1,  0 );
[xindex, xoffset] = deal( 2,  0 );
[   gyo, gxo    ] = deal( 0,  0 );
[yscale, xscale ] = deal( 1,  1 );
% SETUP 
sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% DATA 
etrial = exampleUnit.trialIndex;
eunit  = exampleUnit.id;
emax   = exampleUnit.maxRate;
ecmap  = exampleUnit.colormap;
% PLOT - egofield
plot(pfe{etrial}, eunit, 1, '', ...
     [0, emax],                 ...
     'colorMap',     ecmap,     ...
     'mazeMaskFlag', false,     ...
     'flipAxesFlag', true);
% FORMAT 
sax.XTickLabel =[];
sax.YTickLabel =[];
xlim   (sax, [-250,250]);
ylim   (sax, [-250,250]);
daspect(sax, [1,1,1]);
title  (sax, {'Ego', 'Field'});
% ANNOTATE - rat, circle
subject = struct(rat);
subject = update_subject_patch(subject, 'head',              ...
                               [], false,                    ...
                               bins.hba.edges,               ...
                               bins.hba.centers);
subject = update_subject_patch(subject, 'body',              ...
                               1,  true,                     ...
                               bins.hba.edges([1,end]),      ...
                               [0]);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
patch(subject.body.overlay.vert{:}, [0.75,0.50,0.50], 'FaceAlpha', 0.3, 'EdgeColor','w');
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
circle(0,0,200,'r-');

%--------------------------------------------------------------------------------
%%%>>>
%%%<<< EGO FIELD   --THETA-- HBA Partitioned egofields --------------------------
% DESCRIPTION 
% <description>
% Egocentric ratemaps partitioned by HBA.</description>
for hbaI = 1:bins.hba.count
% COORDINATES 
    [yindex, yoffset] = deal(      1,  0  );
    [xindex, xoffset] = deal( 2+hbaI,  0.15  );
    [   gyo, gxo    ] = deal(      0,  0    );
    [yscale, xscale ] = deal(      1,  1    );
% SETUP         
    sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% DATA 
    etrial = exampleUnit.trialIndex;
    eunit  = exampleUnit.index;
    ecmap = exampleUnit.colormap;    
    xbins  = egoHba.xpos-diff(egoHba.xpos(1:2))/2;
    ybins  = egoHba.ypos-diff(egoHba.ypos(1:2))/2;
    %mask = double(sqrt(bsxfun(@plus,egoHbaRmaps.xbins.^2,egoHbaRmaps.ybins'.^2)') < 445);
    %mask(~mask) = nan;
    %ratemap = ...
    %    fliplr(...
    %        rot90(egoHba.rmap{etrial}(:,:,eunit,hbaI)',-1)).*mask;
% PLOT 
    %set(pcolor(xbins, ybins,ratemap),'EdgeColor','none');
% FORMAT 
    shading (sax, 'flat');
    colormap(sax, func2str(ecmap));
    caxis   (sax, [0, exampleUnit.maxRate]);
    daspect (sax, [1,1,1]);
    title   (sax, hbaBin.label{hbaI});
    box     (sax, 'on');
    axis    (sax, 'xy');
    ylim    (sax, [ybins([1,end])]);
    ylim    (sax, [xbins([1,end])]);
    xlim    (sax, [-250,250]);
    ylim    (sax, [-250,250]);
    sax.XTick =[];
    sax.YTick =[];
    sax.XTickLabel =[];
    sax.YTickLabel =[];
% ANNOTATE 
    Lines([],0,'w');
    Lines(0,[],'w');
    subject = struct(rat);
    subject = update_subject_patch(subject, 'head',              ...
                                   [], false,                    ...
                                   bins.hba.edges,               ...
                                   bins.hba.centers);
    subject = update_subject_patch(subject, 'body',              ...
                                   bins.hba.count+1-hbaI,  true, ...
                                   bins.hba.edges,               ...
                                   bins.hba.centers);
    patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
    patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
    patch(subject.body.overlay.vert{:}, [0.75,0.50,0.50], 'FaceAlpha', 0.3,'EdgeColor','w');
    line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
    line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
end

%--------------------------------------------------------------------------------
%%%>>>    
%%%<<< THETA PHASE --THETA-- sine wave  -----------------------------------------

% DESCRIPTION 
% <description>
% Single theta cycle overlayed with the bins used in the computation of theta
% phase partitioned egocentric rate maps.
% </description>
% COORDINATES 
yscale = (fig.subplot.height*3+fig.subplot.verticalPadding*2)./fig.subplot.height;
[yindex, yoffset] = deal(      4,  -0.15 );
[xindex, xoffset] = deal(      1,  0.35 );
[   gyo, gxo    ] = deal(      0,  0   );
[yscale, xscale ] = deal( yscale,  0.5 );    
% SETUP 
sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% PLOT 
plot(sax, -cos(linspace(0,2*pi)), linspace(0,1), 'k', 'LineWidth', 2);
% Desc
plims = [0.5, 2.26106176905986];
plot(sax,                                        ...
     -cos(linspace(plims(1),plims(2))),          ...
     linspace(plims(1)/(2*pi), plims(2)/(2*pi)), ...
     'Color', bins.phz.color( 1, :),             ...
     'LineWidth', 2);
% Trgh
plims = [2.26106176905986, 4.02212353811972];
plot(sax,                                        ...
     -cos(linspace(plims(1), plims(2))),         ...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)),  ...
     'Color', bins.phz.color( 2, :),             ...
     'LineWidth', 2);
% Ascn 
plims = [4.02212353811972, 5.78318530717959];
plot(sax,                                        ...
     -cos(linspace(plims(1),plims(2))),          ...
     linspace(plims(1)/(2*pi),plims(2)/(2*pi)),  ...
     'Color', bins.phz.color( 3, :),             ...
     'LineWidth', 2);
% FORMAT 
ylim([0,1]);
textOpts = {'Rotation', 90, ...
            'FontSize', 8};
text(sax, 0.5, 0.05,   '0', textOpts);
text(sax,-0.5, 0.5,  '180', textOpts);
text(sax, 0.5, 0.95, '360', textOpts);
% $$$ text( 0.5, 0.05,    '0', textOpts);
% $$$ text(-0.5, 0.5,   '\pi', textOpts);
% $$$ text( 0.5, 0.95, '2\pi', textOpts);
sax.XAxis.Visible = 'off';
sax.YAxis.Visible = 'off';
sax.YTick =[];
sax.Color = 'none';
axes(hfig.UserData.fax);
text(sax.Position(1)-0.25,              ...
     sum(sax.Position([2,4]).*[1, 0.5]),...
     'Theta Phase',                     ...
     'Color', 'k',                      ...
     'Rotation', 90,                    ...
     'HorizontalAlignment', 'center',   ...
     'VerticalAlignment', 'middle',     ...
     'FontSize', 10);
%--------------------------------------------------------------------------------

%%%>>>
%%%<<< EGO FIELD   --THETA-- PHZ partitioned egofields --------------------------

% DESCRIPTION 
% <description>
% Egocentric ratemaps partitioned by TP.
% </description>
for phzI = 1:bins.phz.count
% COORDINATES 
    [yindex, yoffset] = deal( pc+2-phzI,  -0.15  );
    [xindex, xoffset] = deal(         2,  0  );
    [   gyo, gxo    ] = deal(         0,  0  );
    [yscale, xscale ] = deal(         1,  1  );
% SETUP         
    sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% DATA 
    etrial = exampleUnit.trialIndex;
    eunit = exampleUnit.id;
    emax = exampleUnit.maxRate;
    ecmap = exampleUnit.colormap;
% PLOT 
    plot(pfet{etrial}{phzI}, eunit, 1 ,'' , [0, emax],  ...
         'colorMap',     ecmap,                         ...
         'mazeMaskFlag', false,                         ...
         'flipAxesFlag', true);
% FORMAT 
    daspect(sax, [1,1,1]);
    box    (sax,'on');
    xlim   (sax, [-250,250]);
    ylim   (sax, [-250,250]);
    sax(end).XTickLabel =[];
    sax(end).YTickLabel =[];
    Lines([],0,'w');
    Lines(0,[],'w');
end

%--------------------------------------------------------------------------------

%%%>>>
%%%<<< EGO FIELD   --THETA-- partitioned by theta-phase and head-body-angle -----

% DESCRIPTION
% <description>
% Egocentric ratemaps partitioned by HBA and TP.
% </description>
for phzI = 1:bins.phz.count
    for hbaI = 1:bins.hba.count
% COORDINATES 
        [yindex, yoffset] = deal( pc+2-phzI,  -0.15  );
        [xindex, xoffset] = deal(  hbaI + 2,  0.15  );
        [   gyo, gxo    ] = deal(         0,  0  );
        [yscale, xscale ] = deal(         1,  1  );
% SETUP         
        sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% DATA 
        etrial = exampleUnit.trialIndex;
        eunit = exampleUnit.id;
        emax = exampleUnit.maxRate;
        ecmap = exampleUnit.colormap;
% PLOT 
        plot(pfs{exampleUnit.trialIndex}{phzI, hbaI},             ...
             exampleUnit.id, 1, '', [0, exampleUnit.maxRate],     ...
             'colorMap',     ecmap,                               ...
             'mazeMaskFlag', false,                               ...
             'flipAxesFlag', true);
% FORMAT 
        sax.XTickLabel =[];
        sax.YTickLabel =[];
        daspect (sax, [1,1,1]);        
        xlim    (sax, [-250,250]);
        ylim    (sax, [-250,250]);
        box     (sax, 'on');
        caxis   (sax, [0,emax]);
        colormap(sax, func2str(ecmap));
% ANNOTATE 
        Lines([],0,'w');
        Lines(0,[],'w');
        axes(hfig.UserData.fax);
        if phzI == 3
            line([sax.Position(1),sum(sax.Position([1,3]))],...
                 sum(sax.Position([2,4])).*[1,1]+0.1,...
                 'LineWidth',2,...
                 'Color',bins.hba.color(hbaI,:));
        end
        if hbaI == 1
            line([sax.Position(1)-[0.05,0.05]],...
                 [sax.Position([2]),sum(sax.Position([2,4]))],...
                 'LineWidth',2,...
                 'Color',bins.phz.color(phzI,:));
        end
        if phzI==1 && hbaI==3,
            line(sum(sax.Position([1,3])).*[1,1]+0.1,...
                 sax.Position(2).*[1,1]+[0,sax.Position(4)*0.4],...
                 'LineWidth',2,...
                 'Color',[0,0,0]);
            text(sum(sax.Position([1,3]))+0.4,...
                 sax.Position(2),...
                 '20 cm',...
                 'Rotation',90);
            cax = colorbar(sax, 'SouthOutside');
            cax.Units = 'Centimeters';
            cax.Position(2) = cax.Position(2)-1;
            cax.Position(1) = sax.Position(1);
            cax.Position(3) = sax.Position(3);
            cax.Ticks = [0, 7.5, 15];
            cax.TickLabels ={'0', '7.5','   15Hz'};
        end
    end
end

%--------------------------------------------------------------------------------

%%%>>>
%%%>>>---------------------------------------------------------------------------

%%%<<< STATS-BLOCK --------------------------------------------------------------
% START block
% DESCRIPTION 
% <description>
% </description>
%%%<<< (NULL) LAT POS --THETA-- left vs Right lateral coordinates for egoHba ----
% $$$ % DESCRIPTION 
% $$$ % <description>
% $$$ % The lateral, egocentric coordinate (cm) of each field for the left vs right
% $$$ % HBA bins.
% $$$ % </description>
% $$$ % COORDINATES 
% $$$ [yindex, yoffset] = deal(  1,  0  );
% $$$ [xindex, xoffset] = deal(  7,  0  );
% $$$ [   gyo, gxo    ] = deal(  0,  0  );
% $$$ [yscale, xscale ] = deal(  1,  1  );
% $$$ % SETUP 
% $$$ sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% $$$ % PLOT 
% $$$ plot(egoHba.control.meanPos(uids,3,2) ,...
% $$$      egoHba.control.meanPos(uids,1,2) ,...
% $$$      '.',...
% $$$      'MarkerFaceColor','b',...
% $$$      'MarkerEdgeColor','b');
% $$$ % FORMAT subplot
% $$$ grid(sax(end),'on');
% $$$ xlim(sax(end),[-10,10]);
% $$$ ylim(sax(end),[-10,10]);
% $$$ sax(end).XTick = [-10,-5,0,5,10];
% $$$ sax(end).YTick = [-5,0,5];
% $$$ sax(end).XTickLabel = {};
% $$$ title(sax, {'Lateral','Position'}, 'FontWeight', 'Normal');    
% $$$ daspect(sax(end),[1,1,1]);
% $$$ % ANNOTATE 
% $$$ Lines([],0,'k');
% $$$ Lines(0,[],'k');
%%%>>>---------------------------------------------------------------------------
%%%<<< LAT POS     --THETA-- left vs Right lateral coordinates for egoHbaPhz ----
% DESCRIPTION 
% <description>
% The lateral, egocentric coordinate (cm) of each egocentric field for the
% left vs right HBA bins, split by theta phase bins.
% </description>
for phzI = 1:pc
% COORDINATES 
    [yindex, yoffset] = deal( pc+2-phzI, -0.15  );
    [xindex, xoffset] = deal(         7,  0  );
    [   gyo, gxo    ] = deal(         0,  0  );
    [yscale, xscale ] = deal(         1,  1  );
% SETUP 
    sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% PLOT 
    xlim(sax,[-20, 20]);
    ylim(sax,[-20, 20]);
    Lines([],0,'k');
    Lines(0,[],'k');
    plot(egoHbaPhz.control.meanPos(uids,phzI,3,2) ,...
         egoHbaPhz.control.meanPos(uids,phzI,1,2) ,...
         '.',...
         'MarkerFaceColor',phzBin.color(phzI,:),...
         'MarkerEdgeColor',phzBin.color(phzI,:));
% FORMAT 
    grid (sax, 'on');
    xlim (sax, [-10,10]);
    ylim (sax, [-10,10]);    
    sax.XTick = [-5,0,5];
    sax.YTick = [-5,0,5];
    title(sax, '');
    if phzI == 2 
        ylabel(sax,'HBA L (cm)');
        sax.YLabel.Units = 'centimeters';
        sax.YLabel.Position = [-0.55,0.74,0];
    end
    if phzI == 1
        xlabel(sax,'HBA R (cm)');
        sax.XLabel.Units = 'centimeters';
        sax.XLabel.Position = [0.75,-0.4,0];
    else
        sax.XTickLabel = {};
    end
% $$$     if phzI == 3
% $$$         title(sax, {'Mean Lat','Field Pos'}, 'FontWeight', 'Normal');
% $$$     end
    daspect(sax, [1,1,1]);
end
%%%>>>---------------------------------------------------------------------------
%%%<<< LAT POS DSB --THETA-- partitioned by head-body-angle ---------------------
% DESCRIPTION 
% <description>
% The distribution of the lateral, egocentric coordinate (cm) of each egocentric
% field for each HBA bin: Left (green), Center (blue), Right (red).
% </description>
% COORDINATES 
[yindex, yoffset] = deal(  1,  0 );
[xindex, xoffset] = deal(  9, -0.5    );
[   gyo, gxo    ] = deal(  0,  0    );
[yscale, xscale ] = deal(  1,  1   );
% SETUP 
sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% PLOT 
for hbaI = 1:bins.hba.count
    [ehpcmpKDE,dxi] = ksdensity(egoHba.control.meanPos(unitsEgoCA1,hbaI,2));
    plot(dxi,ehpcmpKDE,'-','color',bins.hba.color(hbaI,:))
    med = median(egoHba.control.meanPos(unitsEgoCA1,hbaI,2));
    [~,xi] = NearestNeighbour(dxi,med);
    line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',bins.hba.color(hbaI,:));
end
% FORMAT 
xlim  (sax, [-10,10])
ylim  (sax, [0,0.15])
grid  (sax, 'on');
title (sax, {'Lateral', 'Position'});
sax.YTick = [0,0.05,0.10,0.15];
sax.YTickLabel = {'0','','0.1',''};        
sax.XTick = [-10,-5,0,5,10];
sax.XTickLabel = {};

%--------------------------------------------------------------------------------
%%%>>>
%%%<<< LAT POS DSB --THETA-- partitioned by theta-phase and head-body-angle------
% DESCRIPTION 
% <description>
% The distribution of the lateral, egocentric coordinate (cm) of each egocentric
% field for each HBA bin: Left (green), Center (blue), Right (red). Top:
% ascending TP, Middle: trough TP, Bottom: descending TP.
% </description>
% COORDINATES 

for phzI = 1:pc
% COORDINATES 
    [yindex, yoffset] = deal(  pc-phzI+2, -0.15 );
    [xindex, xoffset] = deal(  9, -0.5    );
    [   gyo, gxo    ] = deal(  0,  0    );
    [yscale, xscale ] = deal(  1,  1   );
% SETUP     
    sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% PLOT 
    for hbaI = 1:bins.hba.count
        meanLatPos = egoHbaPhz.control.meanPos(unitsEgoCA1, phzI, hbaI, lat);
        [ehpcmpKDE,dxi] = ksdensity(meanLatPos);
        med             = median(meanLatPos);
        plot(dxi,ehpcmpKDE,'-','color',bins.hba.color(hbaI,:));
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',bins.hba.color(hbaI,:));
    end
% FORMAT 
    xlim(sax,[-10, 10   ]);
    ylim(sax,[  0,  0.15]);
    grid(sax,'on');
    sax.YTick = [0,0.05,0.10,0.15];
    sax.YTickLabel = {'0','','0.1',''};
    sax.XTick      = [-5, 0, 5];
    if phzI == 1
        xlabel(sax,'cm');
        sax.XLabel.Units = 'centimeters';
        sax.XLabel.Position = [0.75,-0.4,0];
        sax.YTickLabel = {'0', '', '0.1', ''};
    else
        sax.XTickLabel = {};
    end
    if phzI == 2
        ylabel(sax,'(Probability/cm)')
        sax.YLabel.Units = 'centimeters';
        sax.YLabel.Position = [-0.55,0.74,0];
    end
end

%--------------------------------------------------------------------------------
%%%>>>
%%%<<< AP POS DSB  --THETA-- partitioned by head-body-angle ---------------------
% DESCRIPTION 
% <description>
% The distribution of the anteroposterior, egocentric coordinate (cm) of each
% egocentric field for each HBA bin: Left (green), Center (blue), Right (red).
% </description>
% COORDINATES 
    [yindex, yoffset] = deal(  1,  0    );
    [xindex, xoffset] = deal( 10, -0.35 );
    [   gyo, gxo    ] = deal(  0,  0    );
    [yscale, xscale ] = deal(  1,  1.75 );
% SETUP     
    sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% PLOT 
    for hbaI = 1:bins.hba.count
        [ehpcmpKDE,dxi] = ksdensity(egoHba.control.meanPos(unitsEgoCA1,hbaI,1));
        plot(dxi,ehpcmpKDE,'-','color',bins.hba.color(hbaI,:));
        med = median(egoHba.control.meanPos(unitsEgoCA1,hbaI,1));
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',bins.hba.color(hbaI,:));
    end
% FORMAT 
    xlim(sax, [-15,25]);
    ylim(sax, [0,0.15]);
    grid(sax, 'on');
    sax.YTick = [0, 0.05, 0.10, 0.15];
    sax.YTickLabel = {};
    sax.XTick = [-10, -5, 0, 5, 10, 15, 20];
    sax.XTickLabel = {};
    title(sax,                     ...
          {'Mean Field', 'AP Pos'},...
          'FontWeight','Normal',   ...
          'FontSize', 8);

%--------------------------------------------------------------------------------
%%%>>>    
%%%<<< AP POS DSB  --THETA-- partitioned by theta-phase and head-body-angle -----
% DESCRIPTION 
% <description>
% The distribution of the anteroposterior, egocentric coordinate (cm) of each egocentric
% field for each HBA bin: Left (green), Center (blue), Right (red). Top:
% ascending TP, Middle: trough TP, Bottom: descending TP.
% </description>
for phzI = 1:pc
% COORDINATES 
    [yindex, yoffset] = deal( pc-phzI+2, -0.15 );
    [xindex, xoffset] = deal(        10, -0.35 );
    [   gyo, gxo    ] = deal(         0,  0    );
    [yscale, xscale ] = deal(         1,  1.75 );
% SETUP 
    sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% PLOT 
    for hbaI = 1:bins.hba.count
        meanFwdPos = egoHbaPhz.control.meanPos( unitsEgoCA1, phzI, hbaI, fwd);
        [ehpcmpKDE,dxi] = ksdensity(meanFwdPos);
        med             = median(meanFwdPos);
        plot(dxi, ehpcmpKDE, '-', 'color', bins.hba.color(hbaI, :));
        [~,xi] = NearestNeighbour(dxi, med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',bins.hba.color(hbaI, :));
    end
% FORMAT 
    xlim(sax,[-15,20]);
    ylim(sax,[0,0.15]);
    grid(sax,'on');
    sax.YTick      = [0, 0.05, 0.10];
    sax.YTickLabel = {};
    sax.XTick      = [-10, -5, 0, 5, 10, 15, 20];
    sax.XTickLabel = {};    
    switch phzI
      case 1
        xlabel(sax,'    AP Pos (cm)');
        sax.XLabel.Units = 'centimeters';
        sax.XLabel.Position = [0.75, -0.4, 0];
        sax.XTickLabel = {'-10', '', '0', '', '10', '', '20'};
    end
end

%--------------------------------------------------------------------------------
%%%>>>
%%%<<< LAT ZCR CDF --THETA-- partitioned by head-body-angle ---------------------

% DESCRIPTION 
% <description>
% The distribution of the z-scores (Null: randomized HBA) of lateral, egocentric
% coordinate of each egocentric field over HBA bins.
% </description>
% COORDINATES 
    [yindex, yoffset] = deal(  1,  0    );
    [xindex, xoffset] = deal( 12, -0.2  );
    [   gyo, gxo    ] = deal(  0,  0    );
    [yscale, xscale ] = deal(  1,  1    );
% SETUP     
    sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% PLOT 
    for hbaI = 1:bins.hba.count
        [F,X] = ecdf(egoHba.boot.(regions{region}).zscore(:,hbaI,lat));
        plot(X,F,'Color',bins.hba.color(hbaI,:));
    end
% FORMAT 
    grid(sax(end),'on');
    xlim(sax(end),[-15,15]);
    ylim(sax(end),[0,1]);
    sax(end).XTick = [-10,-5,0,5,10];
    sax(end).XTickLabel = {};    
    sax(end).YTick = [0, 0.25,0.5,0.75,1];
    sax(end).YTickLabel = {'', '', '0.5', '', '1'};
    title(sax(end),{'Lateral','z-score'});
% ANNOTATE 
    line(egoHba.perm.(regions{region}).sig.*[1,1],[0,1],...
         'LineStyle', '--', 'Color', [0.25, 0.25, 0.25]);
    line(-egoHba.perm.(regions{region}).sig.*[1,1],[0,1],...
         'LineStyle', '--', 'Color', [0.25, 0.25, 0.25]);

%--------------------------------------------------------------------------------

%%%>>>
%%%<<< LAT ZCR CDF --THETA-- partitioned by theta-phase and head-body-angle -----

% DESCRIPTION 
% <description>
% The distribution of the z-scores (Null: randomized HBA) of lateral,
% egocentric coordinate of each egocentric field for each HBA bin: Left
% (green), Center (blue), Right (red). Top: ascending TP, Middle: trough TP,
% Bottom: descending TP.
% </description>
    for phzI = 1:bins.phz.count
% COORDINATES 
        [yindex, yoffset] = deal( pc-phzI+2, -0.15 );
        [xindex, xoffset] = deal(        12, -0.2  );
        [   gyo, gxo    ] = deal(         0,  0    );
        [yscale, xscale ] = deal(         1,  1    );
% SETUP     
        sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% PLOT 
        samples = ':';
        for hbaI = 1:bins.phz.count
            [F,X] = ecdf(egoHbaPhz.boot.(regions{region}).zscore(samples,phzI,hbaI,lat));
            plot(X,F,'Color',bins.hba.color(hbaI,:));
        end
        % FORMAT subplot
        grid(sax, 'on');
        xlim(sax, [-15, 15]);
        ylim(sax, [0, 1]);
        sax.XTick = [-10, -5, 0, 5, 10];
        sax.XTickLabel = {};        
        sax.YTick = [0, 0.25, 0.5, 0.75, 1];
        sax.YTickLabel = {'', '', '0.5', '', '1'};
        if phzI == 1
            xlabel(sax(end),'z-score');
            sax.XTickLabel = {'-10', '', '0', '', '10'};
        end 
% ANNOTATE 
        line(egoHbaPhz.perm.(regions{region}).sig .* [1, 1], [0, 1], ...
             'LineStyle', '--', 'Color', [0.25, 0.25, 0.25]);
        line(-egoHbaPhz.perm.(regions{region}).sig .* [1, 1], [0, 1],...
             'LineStyle', '--', 'Color', [0.25, 0.25, 0.25]);
    end

%--------------------------------------------------------------------------------    

%%%>>>
%%%<<< AP ZCR CDF  --THETA-- partitioned by head-body-angle ---------------------

% DESCRIPTION 
% <description>
% The distribution of the z-scores (Null: randomized HBA) of anteroposterior,
% egocentric coordinate of each egocentric field over HBA bins.
% </description>
% COORDINATES 
    [yindex, yoffset] = deal(  1,  0   );
    [xindex, xoffset] = deal( 13, -0.2 );
    [   gyo, gxo    ] = deal(  0,  0   );
    [yscale, xscale ] = deal(  1,  1   );
% SETUP     
    sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% PLOT 
    for hbaI = 1:bins.hba.count
        [F,X] = ecdf(egoHba.boot.(regions{region}).zscore(:,hbaI,fwd));
        plot(X,F,'Color',bins.hba.color(hbaI,:));
    end
% FORMAT 
    grid(sax, 'on');
    xlim(sax, [-15, 15]);
    ylim(sax, [0, 1]);
    sax.XTick = [-10 ,-5 ,0 ,5 ,10 ];
    sax.XTickLabel = {};
    sax.YTick = [0, 0.25, 0.5, 0.75, 1];
    sax.YTickLabel = {};
    title(sax, {'AP','z-score'});
    line(egoHba.perm.(regions{region}).sig.*[1, 1], [0, 1],...
         'LineStyle','--','Color',[0.25,0.25,0.25]);
    line(-egoHba.perm.(regions{region}).sig.*[1, 1], [0, 1],...
         'LineStyle', '--', 'Color', [0.25, 0.25, 0.25]);

%--------------------------------------------------------------------------------

%%%>>>
%%%<<< AP ZCR CDF  --THETA-- partitioned by theta-phase and head-body-angle -----

% DESCRIPTION 
% <description>
% The distribution of the z-scores (Null: randomized HBA) of anteroposterior,
% egocentric coordinate of each egocentric field for each HBA bin: Left (green),
% Center (blue), Right (red). Top: ascending TP, Middle: trough TP, Bottom:
% descending TP.
% </description>
    for phzI = 1:bins.phz.count
% COORDINATES 
        [yindex, yoffset] = deal( pc-phzI+2, -0.15 );
        [xindex, xoffset] = deal(        13, -0.2  );
        [   gyo, gxo    ] = deal(         0,  0    );
        [yscale, xscale ] = deal(         1,  1    );
% SETUP     
        sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% PLOT 
        samples = ':';
        for hbaI = 1:bins.phz.count
            [F,X] = ecdf(egoHbaPhz.boot.(regions{region}).zscore(samples,phzI,hbaI,fwd));
            plot(X,F,'Color',bins.hba.color(hbaI,:));
        end
        % FORMAT subplot
        grid(sax, 'on');
        xlim(sax, [-15, 15]);
        ylim(sax, [0, 1]);
        sax.XTick = [-10, -5, 0, 5, 10];
        sax.XTickLabel = {};        
        sax.YTick = [0, 0.25, 0.5, 0.75, 1];
        sax.YTickLabel = {};
        %sax.YTickLabel = {'0', '', '0.5', '', '1'};
        if phzI == 1
            xlabel(sax(end),'z-score');
            sax.XTickLabel = {'-10', '', '0', '', '10'};
        end 
% ANNOTATE 
        line(egoHbaPhz.perm.(regions{region}).sig .* [1, 1], [0, 1], ...
             'LineStyle', '--', 'Color', [0.25, 0.25, 0.25]);
        line(-egoHbaPhz.perm.(regions{region}).sig .* [1, 1], [0, 1],...
             'LineStyle', '--', 'Color', [0.25, 0.25, 0.25]);
    end

%--------------------------------------------------------------------------------    

%%%>>>
%%%<<< SIGNIFICANT POLAR  -------------------------------------------------------
hfig = figure();

%%%<<< SETUP FIGURE -------------------------------------------------------------

% CONFIG 
hfig = figure(666002);
figureFormat      = 'A4';
figureOrientation = 'landscape';
figureUnits       = 'centimeters';
subplotHeight     = 5;%cm
subplotWidth      = 5;%cm
subplotPadVert    = 1;%cm
subplotPadHorz    = 2;%cm
% SETUP 
setup_figure(hfig,                  ... 
             figureFormat,          ... 
             figureOrientation,     ... 
             figureUnits,           ...
             subplotWidth,          ...
             subplotHeight,         ...
             subplotPadHorz,        ...
             subplotPadVert         ...
             );

%%%>>>---------------------------------------------------------------------------
sigranges = [egoHbaPhz.perm.ca1.sig, inf;...
            -egoHbaPhz.perm.ca1.sig,egoHbaPhz.perm.ca1.sig;...
            egoHbaPhz.perm.ca1.sig,inf];
phzI = bins.phz.count;
for hbaI = 1:bins.hba.count
% DESCRIPTION 
% <description>
% The firing rate of an example neuron conditioned on allo-centric position
% (room coordinates), and periods of attentive-prone behavior.</description>
% COORDINATES 
    [yindex, yoffset] = deal(    1,  0 );
    [xindex, xoffset] = deal( hbaI,  0 );
    [   gyo, gxo    ] = deal(    0,  0 );
    [yscale, xscale ] = deal(    1,  1 );
% SETUP 
    sax = setup_axes(hfig,                                               ...
                     yindex, yoffset,                                    ...
                     xindex, xoffset,                                    ...
                     gyo, gxo,                                           ...
                     yscale, xscale,                                     ...
                     @polaraxes);
    posLat = -egoHbaPhz.control.meanPos(unitsEgoCA1, phzI, hbaI, lat);
    posFwd = egoHbaPhz.control.meanPos(unitsEgoCA1, phzI, hbaI, fwd);
    sigInd = WithinRanges(egoHbaPhz.perm.ca1.zscore(:,phzI,hbaI,lat), sigranges(hbaI,:));
    polarhistogram(sax(end),...
                   atan2( posLat(sigInd), posFwd(sigInd)), ...
                   16,... bins
                   'FaceColor',bins.hba.color(hbaI,:));
    rlim(           sax(end), [0,20]);
    hold(           sax(end), 'on');
    rticklabels(    sax(end), {'','','10','','20'});
    thetaticklabels(sax(end), {'0','','','90','','','180','','','270','','',});
    if hbaI==2
        title( sax(end), ...
               {'Ascending Theta Phase: Head to EgoField Angle','',bins.hba.label{hbaI}});
    else
        title(sax(end),...
              bins.hba.label(hbaI));
    end
    tids = {3:5,[6,7,27],18:25,29};
    for tid = 1:4
        uinds = ismember(egoCluSessionMap(:,1),tids{tid});
        -circ_mean(atan2(egoHbaPhz.control.meanPos(uinds,3,hbaI,2),...
                                  egoHbaPhz.control.meanPos(uinds,3,hbaI,1)))
        polarplot(sax(end),...
                  -circ_mean(atan2(egoHbaPhz.control.meanPos(uinds,3,hbaI,2),...
                                  egoHbaPhz.control.meanPos(uinds,3,hbaI,1))).*[1,1],...
                  [0,20],...
                  '-',...
                  'LineWidth',2);
    end
end

%%%>>>---------------------------------------------------------------------------
%%%>>>---------------------------------------------------------------------------

%%%<<< EGO DECODING -------------------------------------------------------------
% DESCRIPTION 
% <description>
% The mean polar cooridinates (angle: left column, radius, right column) of
% the decoded position in the head frame of reference, resolved by theta
% phase (top: ascending, middle: trough, bottom: descending).
% </description>
%%%<<< (NULL) MEAN angle  --THETA-- hba vs hvl -----------------------------------------

% $$$ % COORDINATES 
% $$$     [yindex, yoffset] = deal(  1,  0   );
% $$$     [xindex, xoffset] = deal( 15, -0.4 );
% $$$     [   gyo, gxo    ] = deal(  0,  0   );
% $$$     [yscale, xscale ] = deal(  1,  1   );
% $$$ % SETUP 
% $$$     sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% $$$ % DATA 
% $$$     xbins = mBinHba.centers;
% $$$     ybins = mBinHav.centers;
% $$$     maskedMean = aout_t' .* double(cout_t>200)';
% $$$     maskedMean(maskedMean==0) = nan;
% $$$ % PLOT 
% $$$     imagescnan({xbins, ybins, -maskedMean},...
% $$$                [-pi,pi],...
% $$$                'circular',...
% $$$                false,...
% $$$                [0.8, 0.8, 0.8], 1, 1,...
% $$$                @(x) circshift(hsv(x),round(x/3)));
% $$$     %imagesc(sax, xbins, ybins, -maskedMean);
% $$$ % FORMAT 
% $$$     axis(sax,'ij');
% $$$     axis(sax,'tight');
% $$$     sax.XTick = mBinHba.centers;
% $$$     sax.XTickLabel = {};
% $$$     sax.YTick = [];
% $$$     sax.YTickLabelRotation = 90;
% $$$     sax.YTickLabel = {};
% $$$     title(sax, {'L   C   R',' '},  'FontWeight', 'Normal');
% $$$     axes(hfig.UserData.fax);
% $$$     for hbaI = 1:bins.hba.count
% $$$         xcoords = sax.Position(1) + [hbaI-1, hbaI] *  (sax.Position(3)./3);
% $$$         ycoords = sum(sax.Position([2,4])).*[1,1]+0.15;
% $$$         line(xcoords, ycoords,            ...
% $$$              'LineWidth', 2,              ...
% $$$              'Color',bins.hba.color(hbaI,:));
% $$$     end
% $$$     sax.YTickLabel = round(mBinHav.centers(2:3:end),2);    

%%%>>>
%%%<<< mean angle  --THETA-- hba vs hvl by theta phase --------------------------

for phzI = 1:pc
% COORDINATES 
    [yindex, yoffset] = deal( pc-phzI+2, -0.15 );
    [xindex, xoffset] = deal(        15, -0.4  );
    [   gyo, gxo    ] = deal(         0,  0    );
    [yscale, xscale ] = deal(         1,  1    );
% SETUP 
    sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% DATA 
    xbins = mBinHba.centers;
    ybins = mBinHav.centers;
    maskedMean = aout(:, :, phzI)' .* dmask';
% PLOT 
    imagescnan({xbins, ybins, -maskedMean},...
               [-pi,pi],...
               'circular',...
               false,...
               [0.8, 0.8, 0.8], 1, 1,...
               @(x) circshift(hsv(x),round(x/3)));
% FORMAT 
    axis    (sax, 'ij');
    axis    (sax, 'tight');
    colormap(sax, flipud(circshift(hsv(100),round(100/3))));
    caxis   (sax, [-pi, pi]);
    sax.XTick = mBinHba.centers(2:3:end);
    sax.XTickLabel = {};
    sax.YTick = mBinHav.centers(2:3:end);
    sax.YTickLabelRotation = 90;
    sax.YTickLabel = {};
    switch phzI
      case 1
        cax = colorbar(sax, 'SouthOutside');
        cax.Units = 'Centimeters';
        cax.Position(2) = cax.Position(2) - 1.5;
        cax.Position(1) = sax.Position(1);
        cax.Position(3) = sax.Position(3);
        cax.Ticks = [-pi, 0, pi];
        cax.TickLabels ={'-\pi', '0','\pi'};
        sax.XTick      = [mBinHba.edges(1), mBinHba.centers(5), mBinHba.edges(9)];
        sax.XTickLabel = [mBinHba.edges(1), mBinHba.centers(5), mBinHba.edges(end)];
      case 2
        ylabel(sax, {'Head Angular Vel. (HAV) Direction',' '});
        sax.XTick = mBinHba.centers(2:3:end);
        sax.XTickLabel = {};
        sax.YTick      = mBinHav.centers(2:3:end);
        sax.YTickLabel = mat2cell(bins.hav.key, 1, ones([1,bins.hav.count]));
        sax.YTickLabelRotation = 90;
      case 3
        ylabel(sax,{'Head-Body','Ang Vel'});
        axes (hfig.UserData.fax);
        for hbaI = 1:bins.hba.count
            xcoords = sax.Position(1) + [hbaI-1, hbaI] *  (sax.Position(3)./3)
            ycoords = sum(sax.Position([2,4])) .* [1, 1] + 0.15;
            line(xcoords, ycoords,            ...
                 'LineWidth', 2,              ...
                 'Color',bins.hba.color(hbaI,:));
        end%for
    end%switch
end%phzI




%%%------------------------------------------------------------------------------

%%%>>>
%%%<<< MEAN radius --THETA-- hba vs hvl -----------------------------------------

% COORDINATES 
    [yindex, yoffset] = deal(  1,  0   );
    [xindex, xoffset] = deal( 16, -0.4 );
    [   gyo, gxo    ] = deal(  0,  0   );
    [yscale, xscale ] = deal(  1,  1   );
% SETUP 
    sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% DATA 
    xbins = mBinHba.centers;
    ybins = mBinHav.centers;
    maskedMean = rout_t' .* double(cout_t>200)';
    % PLOT 
    imagescnan({xbins, ybins, maskedMean},...
               [0,10],...
               'linear',...
               false,...
               [0.8, 0.8, 0.8], 1, 1,...
               @copper);
% FORMAT 
    axis(sax,'ij');
    axis(sax,'tight');
    sax.XTick = mBinHba.centers;
    sax.XTickLabel = {};
    sax.YTick = mBinHav.centers;
    sax.YTickLabelRotation = 90;
    sax.YTickLabel = {};
    title(sax, {'L   C   R',' '},  'FontWeight', 'Normal');
    axes(hfig.UserData.fax);
    for hbaI = 1:bins.hba.count
        xcoords = sax.Position(1) + [hbaI-1, hbaI] *  (sax.Position(3)./3)
        ycoords = sum(sax.Position([2,4])).*[1,1]+0.15;
        line(xcoords, ycoords,            ...
             'LineWidth', 2,              ...
             'Color',bins.hba.color(hbaI,:));
    end

%%%>>>
%%%<<< MEAN radius --THETA-- hba vs hvl by theta phase --------------------------

for phzI = 1:pc
% COORDINATES 
    [yindex, yoffset] = deal( pc-phzI+2, -0.15 );
    [xindex, xoffset] = deal(        16, -0.4  );
    [   gyo, gxo    ] = deal(         0,  0    );
    [yscale, xscale ] = deal(         1,  1    );
% SETUP 
    sax = setup_axes(hfig, yindex, yoffset, xindex, xoffset, gyo, gxo, yscale, xscale);
% DATA 
    xbins = mBinHba.centers;
    ybins = mBinHav.centers;
    maskedMean = rout(:, :, phzI)' .* dmask';
% PLOT 
    imagescnan({xbins, ybins, maskedMean},...
               [0,10],...
               'linear',...
               false,...
               [0.8, 0.8, 0.8],1,1,...
               @copper);
% FORMAT 
    axis    (sax, 'ij');
    axis    (sax, 'tight');
    colormap(sax, 'copper');
    caxis   (sax, [0,10]);
    ylabel  (sax, '');
    xlabel  (sax, '');
    sax.XTick = mBinHba.centers;
    sax.XTickLabel = {};
    sax.YTick = mBinHav.centers;        
    sax.YTickLabelRotation = 90;
    sax.YTickLabel = {};
    switch phzI
      case 1        
        cax = colorbar(sax, 'SouthOutside');
        cax.Units = 'Centimeters';
        cax.Position(2) = cax.Position(2)-1.5;
        cax.Position(1) = sax.Position(1);
        cax.Position(3) = sax.Position(3);
        cax.Ticks = [0.25, 7.5];
        cax.TickLabels ={'0', '   7.5cm'};
% $$$       case 3
% $$$         title( sax,{'L   C   R',' '},  'FontWeight', 'Normal');
% $$$         axes (hfig.UserData.fax);
% $$$         for hbaI = 1:bins.hba.count
% $$$             xcoords = [sax.Position(1)+(hbaI-1)*(sax.Position(3)./3), ...
% $$$                        sax.Position(1) + sax.Position(3)./3 * hbaI];
% $$$             ycoords = sum(sax.Position([2,4])).*[1,1]+0.15;
% $$$             line(xcoords, ycoords,            ...
% $$$                  'LineWidth', 2,              ...
% $$$                  'Color',bins.hba.color(hbaI,:));
% $$$         end
    end%switch phzI
end%for phzI

%%%------------------------------------------------------------------------------
%%%>>

% END -- LOC ----- LOC ----- LOC ----- LOC ----- LOC ----- LOC ----- LOC --------
%%%>>>
%%%>>> --------------------------------------------------------------------------

% ENDFIG ************************************************************************
%%% END FIGURE ******************************************************************



%%%<<< SCRATCH ------------------------------------------------------------------

phzI = 3;
lateral = 2;


figure();
hbaI = 1;
plot(  egoHbaPhzLoc.control.meanPos(unitsEgoCA1, phzI, hbaI, lateral),...
     egoHbaPhzPause.control.meanPos(unitsEgoCA1, phzI, hbaI, lateral),'.')

figure();
imagesc(rout(:,:,3)'),
axis('xy'),
colormap('jet');
colorbar();
xlabel('hba');
ylabel('hav');


figure()
for phzI = 1:bins.phz.count
    xco = sq(egoHbaPhz.control.meanPos(uids,phzI,3,1));
    yco = sq(egoHbaPhz.control.meanPos(uids,phzI,3,2));
    [gThetaR,gRadiusR] = cart2pol(xco,yco);
    xco = sq(egoHbaPhz.control.meanPos(uids,phzI,1,1));
    yco = sq(egoHbaPhz.control.meanPos(uids,phzI,1,2));
    [gThetaL,gRadiusL] = cart2pol(xco,yco);
    subplot(3,1,phzI);
    plot(gThetaR,...
         gthetal,...
         '.',...
         'MarkerFaceColor',phzBin.color(phzI,:),...
         'MarkerEdgeColor',phzBin.color(phzI,:));
    Lines([],0,'k');Lines(0,[],'k');
end


%%%<<< HBA SPECTROGRAPH ---------------------------------------------------------

fhba = hba.copy();
%fhba.data = fhba.data(1:100000,:);
defspec = struct('nFFT',2^9,'Fs',fhba.sampleRate,...
                 'WinLength',2^8,'nOverlap',2^8*.875,...
                 'FreqRange',[1,20]);
[ys, fs, ts ] = fet_spec(Trial, fhba, 'mtchglong', true, 'defspec', defspec);
rhm = fet_rhm(Trial);
defspec = struct('nFFT',2^9,'Fs',rhm.sampleRate,...
                 'WinLength',2^8,'nOverlap',2^8*.875,...
                 'FreqRange',[1,20]);
[rys, fs, ts ] = fet_spec(Trial, rhm, 'mtchglong', true, 'defspec', defspec);
figure();
subplot(411);
    imagesc(ts,fs,log10(ys.data'));
    axis('xy');
    caxis([-4,-1.5]);
    colormap('jet');
subplot(412);
    imagesc(ts,fs,log10(rys.data'));
    axis('xy');
    caxis([-4,-1.5]);
    colormap('jet');
subplot(413);
    hold('on');
    plot([1:size(hvang,1)]./hvang.sampleRate,-hvang.data);
    plot([1:size(hba,1)]./hba.sampleRate,hba.data);
subplot(414);
    plotSTC(Trial.stc, 1, 'text',...
            {'lpause','lloc','hpause','hloc','rear','groom','sit'},...
            'bcggrmy');
linkx();

%%%>>> --------------------------------------------------------------------------

figure();
for hbaInd = 1:3,
%%%<<< HBA elips ---------------------------------------------

% Generate or import your data
% Example: random data points
x = egoHbaPhz.control.meanPos(unitsEgoCA1,3,hbaInd,2);
y = egoHbaPhz.control.meanPos(unitsEgoCA1,3,hbaInd,1);

% Combine into a data matrix
data = [x, y];

% Calculate the mean of the data
mu = mean(data);

% Calculate the covariance matrix
cov_matrix = cov(data);

% Compute the eigenvalues and eigenvectors of the covariance matrix
[eigvec, eigval] = eig(cov_matrix);

% Desired confidence level (e.g., 95%)
p = 0.65;
chisquare_val = chi2inv(p, 2);

% Calculate the ellipse parameters
theta_grid = linspace(0, 2*pi, 100);
phi = atan2(eigvec(2,1), eigvec(1,1));
X0 = mu(1);
Y0 = mu(2);
a = sqrt(chisquare_val * eigval(1,1));
b = sqrt(chisquare_val * eigval(2,2));

% Define the ellipse in x and y coordinates
ellipse_x_r = a * cos(theta_grid);
ellipse_y_r = b * sin(theta_grid);

% Define the rotation matrix
R = [cos(phi) sin(phi); -sin(phi) cos(phi)];

% Rotate the ellipse to the correct orientation
r_ellipse = [ellipse_x_r; ellipse_y_r]' * R;

% Plot the data points
hold on;
plot(x, y, '.','Color', bins.hba.color(hbaInd,:),'MarkerSize',20);

% Plot the ellipse
hold on;
plot(r_ellipse(:,1) + X0, r_ellipse(:,2) + Y0, '-', 'LineWidth', 2, 'Color', bins.hba.color(hbaInd,:));
hold off;

% Set the plot labels
xlabel('X-axis');
ylabel('Y-axis');
title(sprintf('Ellipse Containing %.1f%% of Data Points', p*100));
axis equal;
grid on;

%%%>>>
end


lspk = Trials{20}.spk.copy();
lspk.create(Trials{20},1250);
pind = [Trials{20}.stc{'t&w'}];

phz = load_theta_phase(Trials{20},1250);

mxyz = xyz{20}.copy();
mxyz.resample(phz);
headYawCorrection = Trials{20}.meta.correction.headYaw;
headCenterCorrection = Trials{20}.meta.correction.headCenter;


hvec = mxyz(:,'nose',[1,2])-mxyz(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(headYawCorrection),-sin(headYawCorrection);...
                  sin(headYawCorrection), cos(headYawCorrection)],...
                 [2,3],...
                 [1,2]);

unit = 119;
[mxr,mxp] = pft{20}.maxRate(unit);
pfsCenterHR = MTADfet.encapsulate(Trial,                                            ...
                                      [bsxfun(@plus,                                    ...
                                              multiprod(bsxfun(@minus,                  ...
                                                               mxp,                     ...
                                                               sq(mxyz(:,'hcom',[1,2]))),...
                                                         hvec,2,[2,3]),                 ...
                                              headCenterCorrection)],                   ...
                                       mxyz.sampleRate,                                      ...
                                      'egocentric_placefield',                          ...
                                      'egopfs',                                         ...
                                      'p'                                               ...
                                      );
[hfang,hfrad] = cart2pol(pfsCenterHR(:,1),pfsCenterHR(:,2));
figure,
for phzInd = 1:3
    lres = lspk(unit);
    %sind =   (abs(hfang(lres))<0.3|hfrad(lres)<100));
    sind =   (abs(pfsCenterHR(lres,2))<100);
switch phzInd
  case 1
    lres = lres(phz(lres,1)>4&sind);
  case 2
    lres = lres(phz(lres,1)>2&phz(lres,1)<4& sind);
  case 3
    lres = lres(phz(lres,1)<2 & sind);
end
[mccg,t]= CCG(lres,ones(size(lres)), 10, 50, 1250, [1], 'hz');
subplot(3,1,phzInd);
bar(t,mccg);
ylim([0,10]);
Lines([],5,'k');
end

%%%>>> --------------------------------------------------------------------------

mxyz = xyz{20};
mxyz.resample(hba);

hav = copy(hba);
hav.filter('ButFilter',4,2,'low');
hav.data = [0;diff(hba.data)];

figure,
plot(mxyz(:,1)


% ANGULAR velocity head
hvang = filter(copy(mxyz),'ButFilter',4,2,'low');
xycoor = cat(2,...
             hvang(:,'spine_upper',[1,2])-hvang(:,'bcom',[1,2]),...
             hvang(:,'nose',[1,2])-hvang(:,'hcom',[1,2]));
hvang.data = cart2pol(xycoor(:,:,1), xycoor(:,:,2));
% Positive: CCW (Left)     Negative: CW (Right)
hvang.data = circ_dist(circshift(hvang.data(:,2),-10),...
                          circshift(hvang.data(:,2),+10));



%%%<<< PARTS --------------------------------------------------------------------

% SUBPLOTS -- LATERAL POS -- left vs Right lateral coordinatats for egoHba
phzI = bins.phz.count;
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(6, 0, 3+bins.hba.count, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                         ...
                  'Position',[fig.page.xpos(xind)+xOffSet+gXoffset,              ...
                              fig.page.ypos(yind)+yOffSet+gYoffset,              ...
                              fig.subplot.width*1.25,                            ...
                              fig.subplot.height*1.25],                          ...
                  'FontSize', 8,                                                 ...
                  'LineWidth',1);
hold(sax(end),'on');
% PLOT subplot
boxplot( reshape(sqrt(egoHbaPhz.control.meanPos(unitsEgoCA1, phzI, :, lat).^2   ...
          +(egoHbaPhz.control.meanPos(unitsEgoCA1, phzI, :, fwd)-2).^2),[],1),...
         reshape(ones([numel(unitsEgoCA1),1])*[1,2,3],[],1),...
        'symbol',     'r.',...
        'plotstyle',  'traditional',...
         'orientation','vertical');
    
sax(end).XTickLabel = mat2cell(bins.hba.key,1,ones([1,bins.hba.count]));
label(sax(end),'hb-angle');
ylabel(sax(end),'radius (cm)');
ylim(sax(end),[0,20]);
title(sax(end),{'Ascending Theta Phase','Distance to Head'})
% FORMAT subplot
grid(sax(end),'on');
% ylim(sax(end),[-15,20]);
% $$$     sax(end).XTick = [-10,-5,0,5,10];
% $$$     sax(end).YTick = [-5,0,5];



phzI = 3;
unitsEgoCA1Index = false([size(egoCluSessionMap,1),1]);
unitsEgoCA1Index(unitsEgoCA1) = true;
tids = {[3:5], [6,7,27], [18:25], [29]};
for tid = 1:4
    uinds = ismember(egoCluSessionMap(:,1),tids{tid});
    mean(reshape(sqrt(egoHbaPhz.control.meanPos(uinds, phzI, 1, LAT).^2          ...
                    +(egoHbaPhz.control.meanPos(uinds, phzI, 1, FWD)-2).^2),     ...
                 [],1),...
         'omitnan')
end
    


% COLLATE the circular and radial means into the excel format
%     - rows: Samples
%
%     - col A: subject name
%     - col B: sample count (number of cells)
%     - col C: circular mean descending phase 
%     - col D: circular mean trough phase
%     - col E: circular mean ascending phase
%     - col F: radial mean descending phase
%     - col G: radial mean trough phase
%     - col H: radial mean ascending phase


decoded = EgoProCode2D_comp_decoded(dca,'theta','CA1',mfun,beta);
mout = zeros([bins.hba.count,bins.hvl.count]);
sout = zeros([bins.hba.count,bins.hvl.count]);
cout = zeros([bins.hba.count,bins.hvl.count]);
for a = 1:bins.hba.count
    for v = 1:bins.hvl.count
        ind =   WithinRanges( decoded.phz, bins.phz.edges([3,4])) ...
              & WithinRanges(-decoded.hvl, bins.hvl.edges([v,v+1])) ...
              & WithinRanges( decoded.hba, bins.hba.edges([a,a+1])) ...              
              & randn(size(decoded.hba))>0;
        mout(a,v) = mean(decoded.clat(ind),'omitnan');
        sout(a,v) = std(decoded.clat(ind),'omitnan');
        cout(a,v) = sum(ind);
    end
end

SEM = sout./sqrt(cout);               % Standard Error
tscr_l = arrayfun(@(x) tinv([0.025],x), cout-1)
tscr_h = arrayfun(@(x) tinv([0.975],x), cout-1)
CI_l = mout + tscr_l*SEM;
CI_h = mout + tscr_h*SEM;

decoded = EgoProCode2D_comp_decoded(dca,'theta','all',mfun,beta);
%mBinHvl.edges = linspace(-35,35,8);
mBinHvl.edges = linspace(-45,45,8);
%mBinHvl.edges = linspace(-0.5,0.5,8);
mBinHvl.centers = mean([mBinHvl.edges(1:end-1);mBinHvl.edges(2:end)]);
mBinHvl.count = numel(mBinHvl.edges)-1;
mBinHba.edges = linspace(-1.2,1.2,8);
mBinHba.centers = mean([mBinHba.edges(1:end-1);mBinHba.edges(2:end)]);
mBinHba.count = numel(mBinHba.edges)-1;
mout = zeros([mBinHba.count,mBinHvl.count]);
sout = zeros([mBinHba.count,mBinHvl.count]);
tout = zeros([mBinHba.count,mBinHvl.count]);
cout = zeros([mBinHba.count,mBinHvl.count]);
for a = 1:mBinHba.count
    for v = 1:mBinHvl.count
        ind =   WithinRanges(decoded.phz,phzBin.edges([3,4])) ...
              & WithinRanges(decoded.hvl,mBinHvl.edges([v,v+1])) ...
              & WithinRanges(decoded.hba,mBinHba.edges([a,a+1])) ...              
              & randn(size(decoded.hba))>0;
        mout(a,v) = circ_mean(atan2(decoded.clat(ind),decoded.fwd(ind)));
        sout(a,v) = std(decoded.fwd(ind)/10,'omitnan');
        cout(a,v) = sum(ind);
    end
end
mask = ones([mBinHba.count,mBinHvl.count]);
mask(cout<100) = nan;

figure
imagesc(mBinHba.centers,     ...
        mBinHvl.centers,     ...
        mout',               ...
        [-2,2]);
colormap(gca(),'jet')
axis(gca(),'xy');


figure
imagesc(mBinHba.centers,     ...
        mBinHvl.centers,     ...
        mout',               ...
        [-2,2]);
colormap(gca(),'jet')
axis(gca(),'xy');




mBinHvf.edges = linspace(-5,80,8);
mBinHvf.centers = mean([mBinHvf.edges(1:end-1);mBinHvf.edges(2:end)]);
mBinHvf.count = numel(mBinHvf.edges)-1;

mBinHvl.edges = linspace(-45,45,8);
mBinHvl.centers = mean([mBinHvl.edges(1:end-1);mBinHvl.edges(2:end)]);
mBinHvl.count = numel(mBinHvl.edges)-1;

mBinHba.edges = linspace(-1.2,1.2,8);
mBinHba.centers = mean([mBinHba.edges(1:end-1);mBinHba.edges(2:end)]);
mBinHba.count = numel(mBinHba.edges)-1;
mout = zeros([mBinHba.count,mBinHvl.count]);
sout = zeros([mBinHba.count,mBinHvl.count]);
tout = zeros([mBinHba.count,mBinHvl.count]);
cout = zeros([mBinHba.count,mBinHvl.count]);
for a = 1:mBinHba.count
    for v = 1:mBinHvf.count
        ind =   WithinRanges(decoded.phz,phzBin.edges([3,4])) ...
              & WithinRanges(decoded.hvf,mBinHvf.edges([v,v+1])) ...
              & WithinRanges(decoded.hba,mBinHba.edges([a,a+1])) ...              
              & randn(size(decoded.hba))>0;
        mout(a,v) = mean(decoded.fwd(ind)/10,'omitnan');
        sout(a,v) = std(decoded.fwd(ind)/10,'omitnan');
        cout(a,v) = sum(ind);
    end
end
mask = ones([mBinHba.count,mBinHvl.count]);
mask(cout<100) = nan;

figure
imagesc(mBinHba.centers,     ...
        mBinHvf.centers,     ...
        mout',               ...
        [-2,12]);
colormap(gca(),'jet')
axis(gca(),'xy');



decoded = EgoProCode2D_comp_decoded(dca,'theta','all',mfun,beta);

mBinHvf.edges = linspace(-5,80,8);
mBinHvf.centers = mean([mBinHvf.edges(1:end-1);mBinHvf.edges(2:end)]);
mBinHvf.count = numel(mBinHvf.edges)-1;

mBinHvl.edges = linspace(-45,45,8);
mBinHvl.centers = mean([mBinHvl.edges(1:end-1);mBinHvl.edges(2:end)]);
mBinHvl.count = numel(mBinHvl.edges)-1;

mBinHba.edges = linspace(-1.2,1.2,8);
mBinHba.centers = mean([mBinHba.edges(1:end-1);mBinHba.edges(2:end)]);
mBinHba.count = numel(mBinHba.edges)-1;
mout = zeros([mBinHba.count,mBinHvl.count]);
sout = zeros([mBinHba.count,mBinHvl.count]);
tout = zeros([mBinHba.count,mBinHvl.count]);
cout = zeros([mBinHba.count,mBinHvl.count]);
for a = 1:mBinHba.count
    for v = 1:mBinHvf.count
        ind =   WithinRanges(decoded.phz,phzBin.edges([3,4])) ...
              & WithinRanges(decoded.hvf,mBinHvf.edges([v,v+1])) ...
              & WithinRanges(decoded.hba,mBinHba.edges([a,a+1])) ...              
              & randn(size(decoded.hba))>0;
        mout(a,v) = mean(decoded.fwd(ind)/10,'omitnan');
        sout(a,v) = std(decoded.fwd(ind)/10,'omitnan');
        cout(a,v) = sum(ind);
    end
end
mask = ones([mBinHba.count,mBinHvl.count]);
mask(cout<100) = nan;

figure
imagesc(mBinHba.centers,     ...
        mBinHvf.centers,     ...
        mout',               ...
        [-2,12]);
colormap(gca(),'jet')
axis(gca(),'xy');

figure,hist2([decoded.hba,decoded.hvf],linspace(-1.2,1.2,10),linspace(-5,80,10));

% SUBPLOT  - New Axes - y   yo    x   xo       gyo       gxo     ws    hs
sax(end+1) = setup_axes(fig, 6,   0,   8,   0, gYoffset, gXoffset,  1.25, 1.25);
% PLOT subplot
imagesc(sax(end),            ...
        mBinHba.centers,     ...
        mBinHvl.centers,     ...
        mout',               ...
        [-50,50]);
axis(sax(end),'xy');
axis(sax(end),'tight');
colormap('jet');
sax(end).XTick = mBinHba.centers;
sax(end).YTick = mBinHvl.centers ;   
sax(end).XTickLabel = mat2cell(bins.hba.key,1,ones([1,bins.hba.count]));
sax(end).YTickLabel = mat2cell(bins.hba.key,1,ones([1, ...
                    bins.hba.count]));
xlabel(sax(end),'Head Movement');
ylabel(sax(end),'HB Angle (rad)');
sax(end).YTickLabelRotation = 90;
title(sax(end),{'Mean Decoded','Position'})
cax = colorbar(sax(end));
ylabel(cax,'Mean Lateral Position');
cax.Units = 'centimeters';
cax.Position(1) = cax.Position(1)+1;




% SUBPLOT  - New Axes -      y   yo    x   xo       gyo       gxo     ws    hs
sax(end+1) = setup_axes(fig, 6,   0,  10,   1, gYoffset, gXoffset,  1.25, 1.25);
for hbaI = 1:bins.hba.count
    ind = WithinRanges(decoded.phz,[4.5,5.5]) ...
          & randn(size(decoded.hba))>0 ...
          & WithinRanges(decoded.hba,bins.hba.edges(hbaI+[0,1]));      
    [F,xi] = ksdensity(decoded.clat(ind)/10);
    plot(sax(end),xi,F,'-', 'Color', bins.hba.color(hbaI,:),'LineWidth',1);
end
xlim(sax(end),[-40,40]);
lgd = legend(mat2cell(bins.hba.key,1,ones(size(bins.hba.key))),...
             'Location','NorthEastOutside');
lgd.Units = 'centimeters';
lgd.Position(1) = lgd.Position(1)+2
xlabel(sax(end),'Lat Head Pos (cm)')
ylabel(sax(end),'Prob')
title(sax(end),{'Mean Decoded','Position by Hba'})
grid(sax(end),'on');

% $$$ 
% $$$ for t = 1:11
% $$$     ind = WithinRanges(decoded.phz,[4.5,5.5]) ...
% $$$           & randn(size(decoded.hba))>0 ...
% $$$           & WithinRanges(decoded.hba,bins.hba.edges(hbaI+[0,1])),...
% $$$ 


%%%<<< SG2 (theta) MEAN LAT hba vs hvl -------------------------------------------

% SUBPLOT  -   New Axes   -  y   yo   x   xo  gyo  gxo   ws    hs
sax(end+1) = setup_axes(fig, 9,  0,  12,  3,  gyo, gxo,  1.25, 1.25);
% PLOT subplot
imagesc(sax(end), bins.hba.centers, bins.hvl.centers, mout', [-80,80]);
axis(sax(end),'xy');
axis(sax(end),'tight');
colormap(sax(end),'jet');
sax(end).XTick = bins.hba.centers;
sax(end).YTick = bins.hvl.centers;
sax(end).XTickLabel = mat2cell(bins.hba.key,1,ones([1,bins.hba.count]));
sax(end).YTickLabel = mat2cell(bins.hvl.key,1,ones([1,bins.hvl.count]));
xlabel(sax(end),{'Head-Body','Angle'});
ylabel(sax(end),'Head Movement');
sax(end).YTickLabelRotation = 90;
title(sax(end),{'Mean Decoded','Lateral Position'});
cax = colorbar(sax(end));
ylabel(cax,{'cm'});
cax.Units = 'centimeters';
cax.Position(1) = cax.Position(1)+1;
%%%-----------------------------------------------------------------------------

%%%>>>

normalization ='';
figure;
mstates = {'theta','pause','loc'};
for sts = 1:numel(mstates)
decoded = EgoProCode2D_comp_decoded(dca,mstates{sts},'all',mfun,beta);
subplot2(3,3,1,sts);
hist2([decoded.hba,decoded.hav],linspace([-1.2,1.2,30]),linspace(-0.5,0.5,30),normalization);
xlabel('hba');
ylabel('hav');
title(mstates{sts});
subplot2(3,3,2,sts);
hist2([decoded.hba,decoded.hvl],linspace([-1.2,1.2,30]),linspace(-30,30,30),normalization);
xlabel('hba');
ylabel('hvl');
subplot2(3,3,3,sts);
hist2([decoded.hav,decoded.hvl],linspace([-.5,.5,30]),linspace(-30,30,30),normalization);
xlabel('hav');
ylabel('hvl');
end

figure
%%%<<< BINNING ----

for sts = 1:3;
    decoded = EgoProCode2D_comp_decoded(dca, mstates{sts}, 'all', mfun, beta);
% $$$     nmBinHvl.edges = linspace(-45, 45, 8);
% $$$     mBinHvl.centers = mean([mBinHvl.edges(1:end-1);mBinHvl.edges(2:end)]);
% $$$     mBinHvl.count = numel(mBinHvl.edges)-1;
    mBinHva.edges = linspace(-0.3, 0.3, 10);
    mBinHva.centers = mean([mBinHva.edges(1:end-1);mBinHva.edges(2:end)]);
    mBinHva.count = numel(mBinHva.edges)-1;
    mBinHba.edges = linspace(-1.2, 1.2, 10);
    mBinHba.centers = mean([mBinHba.edges(1:end-1);mBinHba.edges(2:end)]);
    mBinHba.count = numel(mBinHba.edges)-1;
    mout  = zeros([mBinHba.count, mBinHva.count, 3]);
    sout  = zeros([mBinHba.count, mBinHva.count, 3]);
    tout  = zeros([mBinHba.count, mBinHva.count, 3]);
    cout  = zeros([mBinHba.count, mBinHva.count, 3]);
    fout  = zeros([mBinHba.count, mBinHva.count, 3]);
    fouts = zeros([mBinHba.count, mBinHva.count, 3]);
    rout  = zeros([mBinHba.count, mBinHva.count, 3]);
    aout  = zeros([mBinHba.count, mBinHva.count, 3]);
    rmout = zeros([mBinHba.count, mBinHva.count, 3]);
    amout = zeros([mBinHba.count, mBinHva.count, 3]);
    for a = 1:mBinHba.count
        for v = 1:mBinHva.count
            ind =    WithinRanges(decoded.phz,bins.phz.edges([3,4])) ...
                   & WithinRanges(-decoded.hav,mBinHva.edges([v,v+1])) ...
                   & WithinRanges(decoded.hba,mBinHba.edges([a,a+1])) ...
                   & randn(size(decoded.hba))>0;
            %$ mean(F),mean(L) -> Polar -> theta, radius
            [aout(a,v,sts),rout(a,v,sts)] = ...
                cart2pol(mean(decoded.fwd(ind),'omitnan'),...
                         mean(decoded.clat(ind),'omitnan'));
            % F, L -> Polar -> mean(Theta), mean(radius)
            [tttt,rrrr] = cart2pol(decoded.fwd(ind), decoded.clat(ind));
            amout(a,v,sts) = circ_mean(tttt);
            rmout(a,v,sts) = mean(rrrr);
            mout(a,v,sts) = mean(decoded.clat(ind),'omitnan');
            sout(a,v,sts) = std(decoded.clat(ind),'omitnan');
            fout(a,v,sts) = mean(decoded.fwd(ind),'omitnan');
            fouts(a,v,sts) = std(decoded.fwd(ind),'omitnan');
            %[~,tout(a,v,sts)] = ttest(decoded.clat(ind));
            cout(a,v,sts) = sum(ind);
        end
    end
    
    mask = ones([mBinHba.count, mBinHva.count]);
    mask(cout(:,:,sts) < 100) = nan;
    subplot2(6,3,1,sts);
    imagescnan({mBinHba.centers, ...
                mBinHva.centers, ...
                rout(:,:,sts)' .* mask'}, ...
               [0,8],           ...
               'linear',true,    ...
               'colorMap',@jet);
    title([mstates{sts}, ' RADIUS']);
    xlabel('HBA');
    ylabel('HAV');
    axis('xy');
    subplot2(6,3,2,sts);    
    imagescnan({mBinHba.centers, ...
                mBinHva.centers, ...
                aout(:,:,sts)' .* mask'}, ...
               [-pi,pi],           ...
               'linear', true,   ...
               'colorMap',@hsv);
    title('ANGLE');
    xlabel('HBA');
    ylabel('HAV');
    axis('xy');
    subplot2(6,3,3,sts);
    imagescnan({mBinHba.centers, ...
                mBinHva.centers, ...
                rmout(:,:,sts)' .* mask'}, ...
               [0,20],           ...
               'linear',true,    ...
               'colorMap',@jet);
    title([mstates{sts}, ' RADIUS']);
    xlabel('HBA');
    ylabel('HAV');
    axis('xy');
    subplot2(6,3,4,sts);    
    imagescnan({mBinHba.centers, ...
                mBinHva.centers, ...
                amout(:,:,sts)' .* mask'}, ...
               [-pi,pi],           ...
               'linear', true,   ...
               'colorMap',@hsv);
    title('ANGLE');
    xlabel('HBA');
    ylabel('HAV');
    axis('xy');
    subplot2(6,3,5,sts);
    plot(rmout(:),rout(:),'.');
    subplot2(6,3,6,sts);
    plot(amout(:),aout(:),'.');
end

%%%>>>


figure();
subplot(211);
imagescnan({mBinHba.centers, ...
            mBinHva.centers, ...
            mout' .* mask'}, ...
           [-7,7],           ...
           'linear',true,    ...
           'colorMap',@jet);
title([state, ' Lat']);
xlabel('HBA');
ylabel('HAV');
axis('xy');
subplot(212);
imagescnan({mBinHba.centers, ...
            mBinHva.centers, ...
            fout' .* mask'}, ...
           [0,10],           ...
           'linear', true,   ...
           'colorMap',@jet);
title('Fwd');
xlabel('HBA');
ylabel('HAV');
axis('xy');


figure
hist2([decoded.dst,decoded.hav],linspace([0,400,30]),linspace(-0.5,0.5,30),'xprob');

figure
hist2([decoded.hba,decoded.hvl],linspace([-1.2,1.2,30]),linspace(-30,30,30),'xprob');

figure
hist2([decoded.hav,decoded.hvl],linspace([-.5,.5,30]),linspace(-30,30,30),'xprob');

%endfig


figure,
phzI = 3;
hold('on')
for hbaI = 1:bins.hba.count
    plot(egoHbaPhz.control.meanPos(unitsEgoCA1,phzI,hbaI,2),egoHbaPhz.control.meanPos(unitsEgoCA1,phzI,hbaI,1),'.','Color',bins.hba.color(hbaI))
end

figure,
for hbaI = 1:bins.hba.count
    for phzI = 1:bins.phz.count
        subplot2(3,3,bins.phz.count+1-phzI,hbaI);
        histogram(sqrt(egoHbaPhz.control.size(unitsEgoCA1,phzI,hbaI))./pi*2,linspace([0,10,20]));;
    end
end

figure,
ind = zeros([size(egoSizeHba,1),1]);
ind(uidsCA3) = 1;
ind = ind & nniz(egoSizeHba);
subplot(211);
histogram(nonzeros(egoSizeHba(ind,3,:))-nonzeros(egoSizeHba(ind,2,:)),linspace(-300,300,20))
subplot(212);
histogram(nonzeros(egoSizeHba(ind,2,:))-nonzeros(egoSizeHba(ind,1,:)),linspace(-300,300,20))


figure,
ind = zeros([size(egoSizeHba,1),1]);
ind(unitsEgoCA1) = 1;
ind = ind & nniz(egoSizeHba);
subplot(211);
histogram(mean(egoSizeHba(ind,3,:),3)-egoSize(ind,3),linspace(-300,300,20))
subplot(212);
histogram(mean(egoSizeHba(ind,2,:),3)-egoSize(ind,2),linspace(-300,300,20))


figure,
hold('on');
plot(egoSize(ind,3),mean(sq(egoSizeHba(ind,3,:)),2),'.');
line([0,800],[0,800]);
plot(egoSize(ind,2),mean(sq(egoSizeHba(ind,2,:)),2),'.g');
line([0,800],[0,400]);

figure,
hold('on');
histogram(egoSize(ind,3)-egoSize(ind,2),linspace(-300,300,20));
histogram(mean(egoSizeHba(ind,3,:),3)-mean(egoSizeHba(ind,2,:),3),linspace(-300,300,20));

figure,
subplot(211);histogram((egoSize(ind,3)-egoSize(ind,2))./(egoSize(ind,3)+egoSize(ind,2)).*100,linspace(-100,100,20));
subplot(212);histogram((median(egoSizeHba(ind,3,:),3)-median(egoSizeHba(ind,2,:),3))./(median(egoSizeHba(ind,3,:),3)+median(egoSizeHba(ind,2,:),3)).*100,linspace(-100,100,20));

[H,P] = ttest2(egoSize(ind,3)-egoSize(ind,2),mean(egoSizeHba(ind,3,:),3)-mean(egoSizeHba(ind,2,:),3));

figure,imagesc(log10(sq(mean(egoMeanRmapRateHba(unitsEgoCA1,:,:),3)))')


figure,
plot(sq(mean(egoMeanRmapRateHba(unitsEgoCA1,1,:),3))./sq(mean(egoMeanRmapRateHba(unitsEgoCA1,2,:),3)),...
     sq(mean(egoMeanRmapRateHba(unitsEgoCA1,3,:),3))./sq(mean(egoMeanRmapRateHba(unitsEgoCA1,2,:),3)),...
        '.')

figure,
plot(egoMeanRmapRate(unitsEgoCA1,2)./egoMeanRmapRate(unitsEgoCA1,3)),...
     sq(mean(egoMeanRmapRateHba(unitsEgoCA1,3,:),3))./sq(mean(egoMeanRmapRateHba(unitsEgoCA1,2,:),3)),...
    '.')




figure
for p = 1:3
    subplot2(3,1,p,1)
    plot(pfsh{exampleUnit.trialIndex}{4-p,1},exampleUnit.id,1,'', ...
         [0,exampleUnit.maxRate],'colorMap',@jet,'mazeMaskFlag',false,'flipAxesFlag',true);
end


%% PLOTING EXAMPLES
lims = {[-250,250],[-250,250]};
figure();
al = 1:5;
numAng = numel(al);
sax = gobjects([0,1]);
for p = 3
    rmap = plot(pfsh{20}{p,1},25,19);
    for a = 1:numAng,
        sax(end+1) = subplot2(1,5,1,a);
        pcolor(pfsh{20}{p,1}.adata.bins{1},...
                    pfsh{20}{p,1}.adata.bins{2},...
                    rmap(:,:,al(a)));
        
        caxis   (sax(end),[0,10]);
        colormap(sax(end),'jet');
        shading (sax(end),'flat');
        axis    (sax(end),'xy');
        xlim    (sax(end),lims{1});
        ylim    (sax(end),lims{2});        
        
        Lines([],0,'k');
        Lines(0,[],'k');
        
        set(sax(end),'XTick',[]);
        set(sax(end),'YTick',[]);        
        
        % ADD subject
% $$$         if p %== 4,
% $$$             subject = struct(rat);
% $$$             subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
% $$$             subject = update_subject_patch(subject,'body', numAng+1-a,  true,hbaBinEdg,hbaBinCtr);
% $$$             patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$             patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$             patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
% $$$             line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
% $$$             line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
% $$$         end
    end
end



%% PLOTING EXAMPLES
lims = {[-300,300],[-300,300]};
figure();
sax = gobjects([0,1]);
trialInd = 20;
unit = 2;
for a = 1:egoHba.hbaBin,
    sax(end+1) = subplot2(1,3,1,a);
    pcolor(egoHba.xbins,...
           egoHba.ybins,...
           egoHba.rmap{trialInd}(:,:,unit,a));
    
    caxis   (sax(end),[0,10]);
    colormap(sax(end),'jet');
    shading (sax(end),'flat');
    axis    (sax(end),'xy');
    xlim    (sax(end),lims{1});
    ylim    (sax(end),lims{2});        
    
    Lines([],0,'k');
    Lines(0,[],'k');
    
    set(sax(end),'XTick',[]);
    set(sax(end),'YTick',[]);        
    
    % ADD subject
% $$$         if p %== 4,
% $$$             subject = struct(rat);
% $$$             subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
% $$$             subject = update_subject_patch(subject,'body', numAng+1-a,  true,hbaBinEdg,hbaBinCtr);
% $$$             patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$             patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$             patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
% $$$             line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
% $$$             line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
% $$$         end
end

%%%>>>



%[p,h,stats]=signtest(egoHbaPhzLoc.control.size(unitsEgoCA1,1,2),egoHbaPhzLoc.control.size(unitsEgoCA1,2,2))
%[p,h,stats]=signtest(egoHbaPhzLoc.control.size(unitsEgoCA1,2,2),egoHbaPhzLoc.control.size(unitsEgoCA1,3,2))

% $$$ mean(egoHbaPhzLoc.control.size(unitsEgoCA1,2,2)*4- ...
% $$$      egoHbaPhzLoc.control.size(unitsEgoCA1,1,2)*4)
% $$$ mean(egoHbaPhzLoc.control.size(unitsEgoCA1,2,1)*4- ...
% $$$      egoHbaPhzLoc.control.size(unitsEgoCA1,2,2)*4)
% $$$ mean(egoHbaPhzLoc.control.size(unitsEgoCA1,2,1)*4- ...
% $$$      egoHbaPhzLoc.control.size(unitsEgoCA1,1,1)*4)

%for hba1 = bins.
[p,h,stats] = signtest(egoHbaPhzLoc.control.size(unitsEgoCA1,3,1)*4, ...
                       egoHbaPhzLoc.control.size(unitsEgoCA1,1,1)*4)

mean(egoHbaPhzLoc.control.size(unitsEgoCA1,3,2)*4-egoHbaPhzLoc.control.size(unitsEgoCA1,2,2)*4)
% 20.6754385964912
mean(egoHbaPhzLoc.control.size(unitsEgoCA1,2,2)-egoHbaPhzLoc.control.size(unitsEgoCA1,1,2))
% 13.1140350877193
[p,h,stats]=signtest(egoHbaPhz.control.size(unitsEgoCA1,1,2),egoHbaPhz.control.size(unitsEgoCA1,2,2))
figure()
for phzI = 1:bins.phz.count
        subplot(1,3,phzI);
        boxplot(sq(bsxfun(...
            @minus,...
            egoHbaPhz.control.size(unitsEgoCA1,phzI,:), ...
            egoPhz.control.size(unitsEgoCA1,phzI))).*16/10000);
        ylim([-0.350,0.1]);
    end
dfs = (egoHbaPhz.control.size(unitsEgoCA1,phzI,2) ...
         - egoPhz.control.size(unitsEgoCA1,phzI));
signtest(dfs(dfs>-200))


phzI = 1;
figure,
for hbaI = 1:bins.hba.count
    subplot2(2,3,1,hbaI);
    hold('on');
    plot(sqrt(egoPhz.control.size(unitsEgoCA1,phzI)*4/pi), ...
     sqrt(egoHbaPhz.control.size(unitsEgoCA1,phzI,hbaI)*4/pi),'.');
line([0,14],[0,14]);
xlabel('CTRL Radius (cm)');
ylabel('HBA Radius (cm)');
title(['Ego-Field Radius HBA - ' bins.hba.label{hbaI}])
set(gca(),'FontSize',24)
subplot2(2,3,2,hbaI);
histogram(sqrt(egoPhz.control.size(unitsEgoCA1,phzI)*4/pi) ...
          -sqrt(egoHbaPhz.control.size(unitsEgoCA1,phzI,hbaI)*4/pi),... 
          linspace(-5,15,50));
title(['Ego-Field Radius Difference CTRL vs HBA - ' bins.hba.label{hbaI}]);
set(gca(),'FontSize',24)
end



figure,plot(sq(sqrt(egoHbaPhz.control.size(unitsEgoCA1,:,1)*4/pi))')

figure,
    plot(reshape(sqrt(egoHbaPhz.control.size(unitsEgoCA1,1,:)*4/pi),[],1), ...
         reshape(sqrt(egoHbaPhz.control.size(unitsEgoCA1,3,:)*4/pi),[],1),...
         '.','MarkerSize',20);
line([0,14],[0,14]);

% high
% hba
% phz,phz+1


figure,plot(decoded.phz())
dclt = decoded.clt(WithinRanges(decoded.phz,bins.phz.edges(end-1:end)));
figure,plot(dclt,circshift(dclt,-1),'.');
figure,plot(dclt(1:6:end),circshift(dclt(1:6:end),-1),'.');
dphz = load_theta_phase(Trials{20},250);


figure,plot(dca{7}.ind,dphz(dca{7}.ind),'-');

dtc = ThreshCross(dphz.data>bins.phz.centers(3),0.5,1);

figure,
plot(dphz.data)
Lines(dtc(:,1),[],'k');

dcai =ismember(dca{7}.ind,dtc(:,1));
% ~13k theta cycles

s = [Trials{7}.stc{'loc&t',1}];
sum(diff(s.data,1,2))*8


figure
for hbaI = 1:3
    for lag = 1:2
        subplot2(3,3,lag,hbaI);
            ind = logical(dca{7}.stcm(:,1)==1) ...
                & logical((dca{7}.stcm(:,3)==3) | (dca{7}.stcm(:,4)==4)| (dca{7}.stcm(:,5)==5) ) ...
 ...                & logical((dca{7}.stcm(:,6)==6) | (dca{7}.stcm(:,7)==7) ) ...
                & sqrt(sum(dca{7}.xyz(:,'hcom',[1,2]).^2,3)) < 300 ...
...                    & WithinRanges( dca{7}.hbang, bins.hba.edges(hbaI+[0,1])) ...
                & WithinRanges( dca{7}.hbang, bins.hba.edges([1,end])) ...                 
                & dcai ;
            tind = false(size(dca{7}.ind));
            tind(ind) = (circshift( dca{7}.ind(ind),-1) - dca{7}.ind(ind)) < 40 ...
                      & (circshift( dca{7}.ind(ind),-2) - dca{7}.ind(ind)) < 80;
            plot(atan2(dca{7}.ecom(tind,2),dca{7}.ecom(tind,1)),...
                 atan2(dca{7}.ecom(circshift(find(tind),-lag),2),dca{7}.ecom(circshift(find(tind),-lag),1)),...
                 '.',...
                 'MarkerSize',20);
            daspect([1,1,1]);grid('on');%xlim([-350,350]);ylim([-350,350]);
            title(['Lag ' num2str(lag) ' vs start, corr:' num2str(corrcoef(atan2(dca{7}.ecom(tind,2),dca{7}.ecom(tind,1)),atan2(dca{7}.ecom(circshift(find(tind),-lag),2),dca{7}.ecom(circshift(find(tind),-lag),1))))])
            set(gca(),'FontSize',24)
            Lines(0,[],'r');
            Lines([],0,'r');
    end
    subplot2(3,3,3,hbaI);    
    plot(dca{7}.ecom(circshift(find(tind), -1),2) - dca{7}.ecom(tind,2),...
         dca{7}.ecom(circshift(find(tind), -2),2) - dca{7}.ecom(tind,2),...
         '.',...
         'MarkerSize',20);
    title('Lag 1 - start vs lag 2 - start')
    daspect([1,1,1]);grid('on');xlim([-350,350]);ylim([-350,350]);set(gca(),'FontSize',24)
end



phzI = 1;
figure,
for hbaI = 1:bins.hba.count
    subplot2(2,3,1,hbaI);
    hold('on');
    plot(sqrt(egoHbaPhz.control.size(unitsEgoCA1,1,hbaI)*4/pi), ...
         sqrt(egoHbaPhz.control.size(unitsEgoCA1,3,hbaI)*4/pi),'.','MarkerSize',20);
line([0,14],[0,14]);
xlabel('Desc Radius (cm)');
ylabel('Asce Radius (cm)');
    title(['Ego-Field Radius Desc vs Asc HBA - ' bins.hba.label{hbaI}])
        set(gca(),'FontSize',24)
        xlim([0,18]);ylim([0,18]);
    subplot2(2,3,2,hbaI);
    hold('on');
    plot(sqrt(egoHbaPhz.control.size(unitsEgoCA1,2,hbaI)*4/pi), ...
         sqrt(egoHbaPhz.control.size(unitsEgoCA1,3,hbaI)*4/pi),'.','MarkerSize',20);
line([0,14],[0,14]);
xlabel('Trgh Radius (cm)');
ylabel('Asce Radius (cm)');
title(['Ego-Field Radius Trgh vs Asc HBA - ' bins.hba.label{hbaI}])
set(gca(),'FontSize',24)
xlim([0,18]);ylim([0,18]);
end

ny = bins.phz.count;
nx = bins.hba.count;
for hbaI = 1:bins.hba.count
    for phzI = 1:bins.phz.count
        subplot2(ny,nx,4-phzI,hbaI);
            plot(log10(egoHbaPhz.control.size(unitsEgoCA1, phzI, hbaI).*16),...
                 log10(egoPhz.control.size(unitsEgoCA1, phzI).*16),'.');
            line(log10([eps,2000]),log10([eps,2000])) 
            line(log10([eps,3448]),log10([eps,3448]))
            %line(log10([eps,2000]),log10([eps,333]))
            daspect([1,1,1])
            xlim([2, 3.55])
            ylim([2, 3.55])
    end
end


figure,
ny = bins.phz.count;
nx = bins.hba.count;
for hbaI = 1:bins.hba.count
    for phzI = 1:bins.phz.count
        subplot2(ny,nx,4-phzI,hbaI);
            plot(log10(egoHbaPhz.control.size(unitsEgoCA1, phzI, hbaI).*16),...
                 log10(egoPhz.control.size(unitsEgoCA1, phzI).*16),'.');
            line(log10([eps,3448/2]),log10([eps,3448]))
            line(log10([eps,3448]),log10([eps,3448]))
            line(log10([eps,3448]),log10([eps,3448/2]))
            %line(log10([eps,2000]),log10([eps,333]))
            daspect([1,1,1])
            xlim([2, 3.55])
            ylim([2, 3.55])
    end
end




u = 21;
egoHbaPhz.control.size(ismember(egoCluSessionMap, [20,u], 'rows'), 3, 1).*16
egoHbaPhz.control.size(ismember(egoCluSessionMap, [20,u], 'rows'), 3, 2).*16
egoHbaPhz.control.size(ismember(egoCluSessionMap, [20,u], 'rows'), 3, 3).*16
egoPhz.control.size(ismember(egoCluSessionMap, [20,u], 'rows'), 3, 1).*16

rmap = pfet{20}{3}.data.rateMap(:,pfet{20}{3}.data.clu==25);
sum(rmap>2)*16
figure,imagesc(reshape(rmap,[41,41])')
imagesc(reshape(rmap,[41,41])')


figure()
plot(pfet{20}{phzI}, 25, 1 ,'' , [0, 2],            ...
     'colorMap',     ecmap,                         ...
     'mazeMaskFlag', false,                         ...
     'flipAxesFlag', true);

figure,
ny = bins.phz.count;
nx = bins.hba.count;
for hbaI = 1:bins.hba.count
    for phzI = 1:bins.phz.count
        subplot2(ny,nx,4-phzI,hbaI);
        histogram(egoHbaPhz.control.size(unitsEgoCA1, phzI, hbaI) ...
                  - egoPhz.control.size(unitsEgoCA1,  phzI),...
                  linspace(-40,40,25))
    end
end



etrial = 20;
eunit = 1;
figure();
for phzI = 1:3;
subplot2(3,4,4-phzI,1)
xbins  = egoPhzRmaps.xpos-diff(egoPhzRmaps.xpos(1:2))/2;
ybins  = egoPhzRmaps.ypos-diff(egoPhzRmaps.ypos(1:2))/2;
ratemap = ...
    fliplr(...
        rot90(egoPhzRmaps.rmap{etrial}(:,:,eunit,phzI)',-1)).*mask;
set(pcolor(xbins, ybins,ratemap),'EdgeColor','none'); 
    caxis([0,2])
    for hbaI = 1:3
        subplot2(3,4,4-phzI,hbaI+1)
    ratemap = ...
    fliplr(...
        rot90(egoHbaPhzRmaps.rmap{etrial}(:,:,eunit,phzI,hbaI)',-1)).*mask;
    set(pcolor(xbins, ybins,ratemap),'EdgeColor','none');
    caxis([0,2])
    colormap('jet');
end
end