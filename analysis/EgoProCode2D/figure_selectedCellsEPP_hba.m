
configure_default_args();
EgoProCode2D_load_data();

% binPhzc
binPhzs = linspace(0.5,2*pi-0.5,4);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
hbaBinEdg = [-1.2,-0.2,0.2,1.2];
hbaBinCtr = mean([hbaBinEdg(1:end-1);hbaBinEdg(2:end)]);


if ~exist('bfrm','var'),bfrm = cf(@(t,u) compute_bhv_ratemaps(t,u),Trials,unitsEgo);end;

lims = {[-300,300],[-300,300]};

% FIGURE - Examples of individual units
%[ 0-120, 120-240, 240-360 ]
ascnPhz = 3;
descPhz = 1;


orientationLabels = 'LCR';

% LOAD patch model 
rat = load_patch_model('rat');

% GENERAL unit info
%    UNIT_ID, SESSION
%    ACCG
%    theta placefield
%    theta behaviorfield (restricted to placefield)
%    ?? HRZ Vs THETA state ratemaps
%     
% MAIN figure
%    BLOCK 3x3 HBAxPHZ

mazeMaskFlag = true;                                      
pfsPlotMode  = 1;
reportMethod ='text';
pfsColorMap  = @jet;
displayModelFlag = true;



%%%<<< Load general variables
sampleRate = 250;
%headCenterCorrection = [-25,-8];
pfsState = 'theta-groom-sit-rear';
hbaBinEdges = -1.5:0.6:1.5;
xyz = cf(@(t) preproc_xyz(t,'trb'),             Trials);
      cf(@(x) x.filter('ButFilter',3,30,'low'), xyz);    
      cf(@(x) x.resample(sampleRate),           xyz);
spk = cf(@(t,u) t.load('spk',sampleRate,'gper',u,'deburst'),Trials,unitsEgo);    
pft = cf(@(t,u)  pfs_2d_theta(t,u),  Trials, unitsEgo);
%%%>>>


overwrite = false;
%%%<<< (pfe) ego ratemap
pfe = cf(@(t,u,x,s,p)                                          ... Egocentric ratemap .
         compute_ego_ratemap(t,u,x,s,p,'overwrite',overwrite), ...
             Trials,                                           ... MTATrial
             unitsEgo,                                         ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft                                               ... MTAApfs object, theta state placefields 
);
%%%>>>

%%%<<< (pfet) egothp ratemap
pfet = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase.
         compute_egothp_ratemap(t,u,x,s,p,'overwrite',overwrite),...
             Trials,                        ... MTATrial
             unitsEgo,                      ... Unit subset, placefields away from the maze walls
             xyz,                           ... MTADxyz object, head position
             spk,                           ... MTASpk object, spike time and id collection 
             pft                            ... MTAApfs object, theta state placefields 
);
%%%>>>


%%%<<< (pfs) egohba ratemap
pfs = cf(@(t,u,x,s,p)                       ... Egocentric ratemap | theta phase , head body angle.
         compute_egohba_ratemap(t,u,x,s,p,'overwrite',overwrite),   ...
             Trials,                        ... MTATrial
             unitsEgo,                      ... Unit subset, placefields away from the maze walls
             xyz,                           ... MTADxyz object, head position
             spk,                           ... MTASpk object, spike time and id collection 
             pft                            ... MTAApfs object, theta state placefields 
);
% $$$ pfsh = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase , shuffled head body angle.
% $$$          compute_egohba_ratemap_shuffled(t,u,x,s,p,'overwrite',overwrite),   ...
% $$$              Trials,                        ... MTATrial
% $$$              units,                         ... Unit subset, placefields away from the maze walls
% $$$              xyz,                           ... MTADxyz object, head position
% $$$              spk,                           ... MTASpk object, spike time and id collection 
% $$$              pft                            ... MTAApfs object, theta state placefields 
% $$$ );
%%%>>>



t = 20;
u = 3;

pfo = pfs; % hba <- req20200630.m
% $$$ pfo = pfl; % hvl <- req20200630.m
% $$$ pfo = pfv; % hvf <- req20200630.m
% $$$ pfo = pfr; % roll<- req20200630.m

% change incase first is empty
nBins = numel(hbaBinCtr);
nPhz = numel(binPhzc);


for  t = 1:numel(Trials),
    %for  t = 1:numel(Trials),
    Trial = Trials{t};    
    Trial.load('nq');
    figDir = fullfile('/storage/share/Projects/EgoProCode2D/eppAllUnits',Trial.filebase);
    create_directory(figDir);
    [accg,tbins] = autoccg(Trials{t});

% HBA setup
    headBodyAng = [xyz{t}(:,'spine_upper',[1,2])-xyz{t}(:,'bcom',[1,2]),...
                   xyz{t}(:,'nose',[1,2])-xyz{t}(:,'hcom',[1,2])];
    headBodyAng = sq(bsxfun(@rdivide,headBodyAng,sqrt(sum(headBodyAng.^2,3))));
    headBodyAng = cart2pol(headBodyAng(:,:,1),headBodyAng(:,:,2));
    headBodyAng = circ_dist(headBodyAng(:,2),headBodyAng(:,1));
    headBodyAng = MTADfet.encapsulate(Trials{t},                                    ...
                                      -(headBodyAng+Trial.meta.correction.headBody),...
                                      sampleRate,'hba','hba','h');
    hbaBinInd = discretize(headBodyAng.data, hbaBinEdg);

% $$$ 
% $$$ % HVF setup
% $$$ hvfl = fet_href_HXY(Trials{t},sampleRate,false,'trb',4);
% $$$ hvfBinEdg = -10:20:90;
% $$$ hvfBinCtr = mean([hvfBinEdg(1:end-1);hvfBinEdg(2:end)]);
% $$$ hvfBinInd = discretize(hvfl(:,1), hvfBinEdg);
% $$$ 
% $$$ % HRL setup
% $$$ hang = transform_origin(Trial,filter(copy(xyz{t}),'ButFilter',4,20),'hcom','nose',{'hcom','head_right'});
% $$$ hroll = MTADfet.encapsulate(Trial,...
% $$$                                   -(hang.roll+hrlCorrection(t)),...
% $$$                                   xyz{t}.sampleRate,...
% $$$                                   'hrl','hrl','r');
% $$$ hrlBinEdg = -0.6:0.24:0.6;
% $$$ hrlBinCtr = mean([hrlBinEdg(1:end-1);hrlBinEdg(2:end)]);
% $$$ hrlBinInd = discretize(hroll(:,1), hrlBinEdg);
% $$$ 


varBinEdg = hbaBinEdg;
varBinCtr = hbaBinCtr;
varBinInd = hbaBinInd;

% $$$ varBinEdg = hvfBinEdg;
% $$$ varBinCtr = hvfBinCtr;
% $$$ varBinInd = hvfBinInd;
% $$$ 
% $$$ 
% $$$ varBinEdg = hrlBinEdg;
% $$$ varBinCtr = hrlBinCtr;
% $$$ varBinInd = hrlBinInd;


hvec = xyz{t}(:,'nose',[1,2])-xyz{t}(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(Trial.meta.correction.headYaw),-sin(Trial.meta.correction.headYaw);...
                  sin(Trial.meta.correction.headYaw),cos(Trial.meta.correction.headYaw)],...
                 [2,3],...
                 [1,2]);


ind = [Trials{t}.stc{'w+p+n&t',sampleRate}];
%ind = [Trials{t}.stc{'theta-groom-sit-rear',sampleRate}];

ind.cast('TimeSeries');
ind.resample(xyz{t});


if isempty(unitsEgo{t}),
    continue,
end
for u = 1:numel(unitsEgo{t}),

[mxr,mxp] = pft{t}.maxRate(unitsEgo{t}(u));        
pfstrj = MTADfet.encapsulate( ...
        Trials{t},...
        multiprod(bsxfun(@minus,...
                         mxp,...
                         sq(xyz{t}(:,'hcom',[1,2]))),...
                  hvec,2,[2,3]),...
        sampleRate,...
        'pfstrj','ptrj','t');

    
% SET figure layout
%[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','landscape',[],2.5,2.5,0.5,0.5);
[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','portrait',[],2.5,2.5,0.2,0.2);
[yind, yOffSet, xind, xOffSet] = deal(1,                         ... yind
                                      0,                         ... yOffSet
                                      1,                         ... xind
                                      0);                          % xOffSet

xBlockOffset = 0;
yBlockOffset = 0;

unit = unitsEgo{t}(u);

%rateMinMax = max(cell2mat(cf(@(p) prctile(p.data.rateMap(:,p.data.clu==unit,1),99.5), pfo{t})));
rateMinMax = max(cell2mat(cf(@(p) prctile(p.data.rateMap(:,p.data.clu==unit,1),99.9), {pfe{t}})));

% INFORMATION -------------------------------------------------------------------------------------
[yind, yOffSet, xind, xOffSet] = deal(   1,       -1,    1,       0);        
FigInfo =                                                                                       ...
    uicontrol('Parent',hfig,                                                                    ...
              'Style','text',                                                                   ...
              'String',{['Unit: ',num2str(unit)],                                               ...
                         Trial.filebase,                                                        ...
                         ['stcMode: ',Trial.stc.mode],                                          ...
                         ['eDist:   ',num2str(Trial.nq.eDist(unit))],                           ...
                         ['Refrac:  ',num2str(log10(Trial.nq.Refrac(unit)))],                   ...
                         ['SNR:     ',num2str(Trial.nq.SNR(unit))],                             ...
                         ['AmpSym:  ',num2str(Trial.nq.AmpSym(unit))],                          ...
                         ['SpkWidthR:  ',num2str(Trial.nq.SpkWidthR(unit))]                     ...
              },                                                                                ...
              'Units','centimeters',                                                            ...
              'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,                             ...
                          fig.page.ypos(yind+yBlockOffset)+yOffSet,                             ...
                          fig.subplot.width*2,                                                  ...
                          fig.subplot.height*2],                                                  ...
              'FontSize', 8);
% -------------------------------------------------------------------------------------------------

% THETA PLACE FIELD -------------------------------------------------------------------------------
[yind, yOffSet, xind, xOffSet] = deal(   3,       1,    2,       fig.subplot.width/2);        
sax(end+1) = axes('Units','centimeters',                                                        ...
                  'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,                         ...
                              fig.page.ypos(yind+yBlockOffset)+yOffSet,                         ...
                              fig.subplot.width,                                                ...
                              fig.subplot.height],                                              ...
                  'FontSize', 8,                                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(pft{t},                  ... % MTAApfs object
     unit,                    ... % unit id
     pfsPlotMode,             ...
     reportMethod,            ... 
     rateMinMax,              ...
     mazeMaskFlag,            ...
     'colorMap',pfsColorMap);
sax(end).XTick = [];
sax(end).YTick = [];
% -------------------------------------------------------------------------------------------------

% THETA BEHAVIOR FIELD-----------------------------------------------------------------------------
[yind, yOffSet, xind, xOffSet] = deal(   3,       1,    3,       fig.subplot.width/2);        
sax(end+1) = axes('Units','centimeters',                                                        ...
                  'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,                         ...
                              fig.page.ypos(yind+yBlockOffset)+yOffSet,                         ...
                              fig.subplot.width,                                                ...
                              fig.subplot.height],                                              ...
                  'FontSize', 8,                                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(bfrm{t},                 ... % MTAApfs object
     unit,             ... % unit id
     pfsPlotMode,             ...
     reportMethod,            ... 
     rateMinMax,              ...
     ~mazeMaskFlag,           ...
     'colorMap',pfsColorMap);
sax(end).XTick = [];
sax(end).YTick = [];
xlabel('head pitch');
ylabel('body pitch');
axis(sax(end),'tight')
% -------------------------------------------------------------------------------------------------

% THETA accg --------------------------------------------------------------------------------------
[yind, yOffSet, xind, xOffSet] = deal(   3,       1,    4,       fig.subplot.width/2);        
sax(end+1) = axes('Units','centimeters',                                                        ...
                  'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,                         ...
                              fig.page.ypos(yind+yBlockOffset)+yOffSet,                         ...
                              fig.subplot.width,                                                ...
                              fig.subplot.height],                                              ...
                  'FontSize', 8,                                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
bar(tbins,accg(:,unit));axis tight;
sax(end).XTick = [];
sax(end).YTick = [];
axes(fax);
line([sax(end).Position(1),sax(end).Position(1)+sax(end).Position(3)/6],...
      sax(end).Position(2).*[1,1]-0.25,'Color','k','LineWidth',2);
text([sax(end).Position(1)+sax(end).Position(3)/5],...
     sax(end).Position(2)-0.25,...
     '10ms');
% -------------------------------------------------------------------------------------------------

% THETA ego ---------------------------------------------------------------------------------------
[yind, yOffSet, xind, xOffSet] = deal(   3,       1,    1,       fig.subplot.width/2);        
sax(end+1) = axes('Units','centimeters',                                                        ...
                  'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,                         ...
                              fig.page.ypos(yind+yBlockOffset)+yOffSet,                         ...
                              fig.subplot.width,                                                ...
                              fig.subplot.height],                                              ...
                  'FontSize', 8,                                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plot(pfe{t},                  ... % MTAApfs object
     unit,                    ... % unit id
     pfsPlotMode,             ...
     reportMethod,            ... 
     rateMinMax,              ...
     false,            ...
     'colorMap',pfsColorMap,'flipAxesFlag',true);
xlim    (sax(end),lims{1});
ylim    (sax(end),lims{2});        
sax(end).XTick = [];
sax(end).YTick = [];
% -------------------------------------------------------------------------------------------------

% THETA ego phz -----------------------------------------------------------------------------------
for p = 1:nPhz,
[yind, yOffSet, xind, xOffSet] = deal(   p+5,       1,    1,       fig.subplot.width/2);        
sax(end+1) = axes('Units','centimeters',                                                        ...
                  'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,                         ...
                              fig.page.ypos(yind+yBlockOffset)+yOffSet,                         ...
                              fig.subplot.width,                                                ...
                              fig.subplot.height],                                              ...
                  'FontSize', 8,                                                                ...
                  'LineWidth',1);

hold(sax(end),'on');

pcolor(pfet{t}{p}.adata.bins{1},...
       pfet{t}{p}.adata.bins{2},...
       plot(pfet{t}{nPhz+1-p},unit,[],[],[],false));
        
caxis   (sax(end),[0,rateMinMax]);
colormap(sax(end),'jet');
shading (sax(end),'flat');
axis    (sax(end),'xy');
xlim    (sax(end),lims{1});
ylim    (sax(end),lims{2});        

Lines([],0,'k');
Lines(0,[],'k');
sax(end).XTick = [];
sax(end).YTick = [];
end
% -------------------------------------------------------------------------------------------------

% 
for a = 1:nBins,
    [yind, yOffSet, xind, xOffSet] = deal(   5,       1,    a+2,       0);        
    sax(end+1) = axes('Units','centimeters',                                                        ...
                      'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,                         ...
                                  fig.page.ypos(yind+yBlockOffset)+yOffSet,                         ...
                                  fig.subplot.width,                                                ...
                                  fig.subplot.height],                                              ...
                      'FontSize', 8,                                                                ...
                      'LineWidth',1);
    
    hold(sax(end),'on');
    
    indR = ind.data.*varBinInd==a;
    %hist2([pfstrj(indR,2),pfstrj(indR,1)],linspace([-300,300,100]),linspace([-300,300,100]));
    %scatter(pfstrj(indR,2),pfstrj(indR,1),1,headBodyAng(indR),'filled');
    plot(pfstrj(indR,2),pfstrj(indR,1),'.','MarkerSize',1);
    %colormap('jet');
    %caxis([-1,1]);
    xlim(sax(end),lims{1});
    ylim(sax(end),lims{2});        
end

% MAINBLOCK - egofields----------------------------------------------------------------------------

if displayModelFlag,
    xBlockOffset = 2;
    yBlockOffset = 0;
    for a = 1:nBins
        [yind, yOffSet, xind, xOffSet] = deal(   nBins+6,       1,    a,       0);                
        sax(end+1) = axes('Units','centimeters',                                        ...
                          'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
                            fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
                            fig.subplot.width,                                ...
                            fig.subplot.height],                              ...
                          'FontSize', 8,                                                ...
                          'LineWidth',1);
        hold(sax(end),'on');
        subject = struct(rat);
        subject = update_subject_patch(subject,'head', a, true,varBinEdg,varBinCtr);
        subject = update_subject_patch(subject,'body',[],false,varBinEdg,varBinCtr);
        patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
        patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
        patch(subject.head.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
        line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
        line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
        daspect(sax(end),[1,1,1]);
        xlim([-100,100]);
        ylim([-100,100]);    
        sax(end).XTick = [];
        sax(end).YTick = [];
    end
end% if displayModelFlag

xBlockOffset = 2;
yBlockOffset = 0;



for p = 1:nPhz,

    for a = 1:nBins,
        rmap = plot(pfo{t}{p,a},unit,1,[],[],false);        
        [yind, yOffSet, xind, xOffSet] = deal(   nPhz+6-p,       1,    a,       0);                
        sax(end+1) = axes('Units','centimeters',                                        ...
                          'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
                                      fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
                                      fig.subplot.width,                                ...
                                      fig.subplot.height],                              ...
                          'FontSize', 8,                                                ...
                          'LineWidth',1);
        hold(sax(end),'on');
        
        pcolor(pfo{t}{p,a}.adata.bins{1},...
               pfo{t}{p,a}.adata.bins{2},...
               rmap(:,:));
        
        caxis   (sax(end),[0,rateMinMax]);
        colormap(sax(end),'jet');
        shading (sax(end),'flat');
        axis    (sax(end),'xy');
        xlim    (sax(end),lims{1});
        ylim    (sax(end),lims{2});        
        
        Lines([],0,'k');
        Lines(0,[],'k');
        
        if a == 1,    
            ylabel({round(circ_rad2ang(binPhzc(p))),' '});
            if p==round(nPhz/2),
                ylabel({'Theta Phase',round(circ_rad2ang(binPhzc(p))),' '});
            end
        end
        if p == nPhz && a == round(nBins/2),
            %title([pfs{t}{end}.tag,'_unit_',num2str(unit)]);
        end
        
        set(sax(end),'XTick',[]);
        set(sax(end),'YTick',[]);        

        
        % ADD subject
        if displayModelFlag,
            subject = struct(rat);
            subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
            subject = update_subject_patch(subject,'body', nBins+1-a,  true,hbaBinEdg,hbaBinCtr);
            patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
            patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
            patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
            line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
            line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
        else
            if p == 1,
                xlabel([num2str(varBinCtr(a)),' mm/s']);
            end
        end;%if displayModelFlag
    end;%for a
end;%for p

axes(fax);
line(sax(end).Position(1)+sax(end).Position(3).*[1,1]+0.25,               ...
     [sax(end).Position(2),                                               ...
      sax(end).Position(2)+sax(end).Position(4)/5],                       ...
     'LineWidth',2,                                                       ...
     'Color','k');

text(sax(end).Position(1)+sax(end).Position(3)+0.5,                       ...
     sax(end).Position(2),                                                ...
     '10cm',                                                              ...
     'Rotation',90,                                                       ...
     'FontSize',8);

line(sax(end-nPhz*nBins+1).Position(1).*[1,1]-0.25,                    ...
     [sax(end-nPhz*nBins+1).Position(2),                               ...
      sax(end-nPhz+1).Position(2)+sax(end-nPhz).Position(4)],         ...
     'LineWidth',2,                                                       ...
     'Color','k');

figName = [pfo{t}{end}.tag,'_unit_',num2str(unit,'%04.f'),'.png'];
print(hfig,                                                               ...% figure handle
      '-dpng',                                                            ...% image format
      fullfile(figDir,figName));                                             % image path

% -------------------------------------------------------------------------------------------------
        
end% for u : units
end% for t : trials

%% COMPUTE rate weighted field center 


% $$$ figure();
ucounter = 1;
egoMeanRmapPosHba = [];
egoMaxRmapPosHba = [];
egoSizeHba = [];
egoMeanRmapRateHba = [];
egoMaxRmapRateHba = [];

for t = 1:numel(Trials)
for u = 1:numel(unitsEgo{t})
unit = unitsEgo{t}(u);
for p = 1:3;
for a = 1:3,
binSubsetX = abs(pfs{t}{p,a}.adata.bins{1})<300;
binSubsetY = abs(pfs{t}{p,a}.adata.bins{2})<300;
mapPosition = cell([1,2]);
[mapPosition{:}] = ndgrid(pfs{t}{p,a}.adata.bins{1}(binSubsetX),...
                          pfs{t}{p,a}.adata.bins{2}(binSubsetY));
mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
rmap = plot(pfs{t}{p,a},unit,[],[],[],false);
rmap = rmap(binSubsetX,binSubsetY);
nanmap = double(~isnan(rmap));
nanmap(nanmap==0) = nan;
rmap = rmap.*fliplr(nanmap);
rmap(rmap<2) = 0;
% $$$ subplot2(3,3,4-p,a);
% $$$ hold('on');
% $$$ imagescnan({pfs{t}{p,a}.adata.bins{1}(binSubsetX),...
% $$$                    pfs{t}{p,a}.adata.bins{2}(binSubsetY),...
% $$$                    rmap'});
nrmap =rmap./sum(rmap(:),'omitnan');
rmapCenter = sq(sum(sum(bsxfun(@times,nrmap,mapPosition),'omitnan'),'omitnan'))';
egoSizeHba(ucounter,p,a) = sum(nniz(nrmap(:)));
egoMeanRmapRateHba(ucounter,p,a,:) = mean(rmap(nniz(rmap(:))),'omitnan');

% $$$ plot(rmapCenter(1),rmapCenter(2),'*m')
% $$$ axis('tight')
% $$$ Lines([],0,'g');
% $$$ Lines(0,[],'g');
egoMeanRmapPosHba(ucounter,p,a,:) = rmapCenter./10;%+ [-2,1.6]*ismember(t,[3:5,18:26,30]);

[~,maxPos] = max(nrmap(:));
if ~isempty(maxPos)
[maxX,maxY] = ind2sub(size(nrmap),maxPos);
egoMaxRmapPosHba(ucounter,p,a,:) = mapPosition(maxX,maxY,:);
else
egoMaxRmapPosHba(ucounter,p,a,:) = nan([1,1,1,2]);
end
end
end
ucounter = ucounter+1;
end
end

uids = uidsCA1;
uids = uidsCA3;

figure,
hold('on');
plot(egoMeanRmapPos(uids,p,1,2)-egoMeanRmapPos(uids,p,2,2),...
     egoMeanRmapPos(uids,p,1,1)-egoMeanRmapPos(uids,p,2,1),['.',hclr(1)])
plot(egoMeanRmapPos(uids,p,3,2)-egoMeanRmapPos(uids,p,2,2),...
     egoMeanRmapPos(uids,p,3,1)-egoMeanRmapPos(uids,p,2,1),['.',hclr(3)])

figure,
hold('on');
histogram(egoMeanRmapPos(uids,p,1,2)-egoMeanRmapPos(uids,p,2,2),-150:25:150)
histogram(egoMeanRmapPos(uids,p,3,2)-egoMeanRmapPos(uids,p,2,2),-150:25:150)

figure,
hold('on');
histogram(egoMeanRmapPos(uids,p,1,2)-egoMeanRmapPos(uids,p,3,2),-150:25:150)

figure();
hold('on');
histogram(egoMeanRmapPos(uids,p,1,2)-egoMeanRmapPos(uids,p,2,2),-150:25:150,'Normalization','probability')
histogram(egoMeanRmapPos(uids,p,3,2)-egoMeanRmapPos(uids,p,2,2),-150:25:150,'Normalization','probability')

figure,
    plot(egoMeanRmapPos(uids,p,1,2),...
         egoMeanRmapPos(uids,p,3,2),...
         '.');
    line([-150,150],[150,-150],'Color','k')
    xlim([-150,150]);ylim([-150,150]);
    Lines([],0,'k');Lines(0,[],'k');    
    xlabel('Theta240 Turned Left');
    ylabel('Theta240 Turned Right');

figure,
    plot(egoMeanRmapPos(uids,p,1,2),...
         egoMeanRmapPos(uids,p,2,2),...
         '.');
    Lines([],0,'k');Lines(0,[],'k');

figure();
    hold('on');
    xlim([-150,150]);ylim([-150,150]);    
    line([-150,150],[150,-150],'Color','k')
    Lines([],0,'k');Lines(0,[],'k');    
    plot(egoMeanRmapPos(uids,p,1,2)-egoMeanRmapPos(uids,p,2,2),...
         egoMeanRmapPos(uids,p,3,2)-egoMeanRmapPos(uids,p,2,2),...
         '.');
    
p = 3
corr(egoMeanRmapPos(uids,p,1,2)-egoMeanRmapPos(uids,p,2,2),...
     egoMeanRmapPos(uids,p,3,2)-egoMeanRmapPos(uids,p,2,2))




figure,
hold('on')
cdfplot(egoMeanRmapPos(uids,p,1,2)-egoMeanRmapPos(uids,p,2,2));
cdfplot(egoMeanRmapPos(uids,p,3,2)-egoMeanRmapPos(uids,p,2,2));
cdfplot(egoMeanRmapPos(uids,p,3,2)-egoMeanRmapPos(uids,p,1,2));
xlim([-150,150]);

hclr = 'bgr';
fclr = 'ckm';

figure();
for p = 1:3;
    subplot2(3,4,1,p);
    hold('on');
    for a = 1:3
        plot(egoMeanRmapPos(uids,p,a,2),egoMeanRmapPos(uids,p,a,1),['.',hclr(a)])
        xlim([-150,150]);
        ylim([-150,250]);
    end

    subplot2(3,4,2,p);

    hold('on');
    hclr = 'bgr';
    for a = 1:3
        histogram(egoMeanRmapPos(uids,p,a,2),linspace(-300,300,50),'FaceColor',hclr(a))
        xlim([-150,150]);
    end

    subplot2(3,4,3,p);
    hold('on');
    for a = 1:3
        lh = cdfplot(egoMeanRmapPos(uids,p,a,2));
        lh.Color = hclr(a);
        %histogram(egoMeanRmapPos(:,3,a,2),linspace(-300,300,50),'FaceColor',hclr(a))
        xlim([-150,150]);
    end

    subplot2(3,4,1,4);
    hold('on');
    lh = cdfplot(nonzeros(egoMeanRmapPos(uids,p,:,1)));
    lh.Color = fclr(p);
    %histogram(egoMeanRmapPos(:,3,a,2),linspace(-300,300,50),'FaceColor',hclr(a))
    xlim([-150,250]);
end


[H,P,CI,STATS] = ttest2(egoMeanRmapPos(uids,3,1,2),egoMeanRmapPos(uids,3,2,2))
[H,P,CI,STATS] = ttest2(egoMeanRmapPos(uids,3,2,2),egoMeanRmapPos(uids,3,3,2))
[H,P,CI,STATS] = ttest2(egoMeanRmapPos(uids,3,1,2),egoMeanRmapPos(uids,3,3,2))


figure
hold('on');
[h,L,MX,MED,bw] = violin({egoMeanRmapPos(uidsCA1,3,3,2)-egoMeanRmapPos(uidsCA1,3,2,2), ...
                          egoMeanRmapPos(uidsCA1,3,1,2)-egoMeanRmapPos(uidsCA1,3,2,2), ...
                          egoMeanRmapPos(uidsCA1,3,1,2)-egoMeanRmapPos(uidsCA1,3,3,2) ...
                   });


figure
hold('on');
[h,L,MX,MED,bw] = violin({egoMeanRmapPos(uidsCA1,3,3,2), ...
                          egoMeanRmapPos(uidsCA1,3,2,2) ...
                          egoMeanRmapPos(uidsCA1,3,1,2) ...
                   });
