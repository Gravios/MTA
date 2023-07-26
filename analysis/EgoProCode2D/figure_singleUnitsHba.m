
% $$$ configure_default_args();
% $$$ EgoProCode2D_load_data();

% binPhzc
binPhzs = linspace(0.5,2*pi-0.5,4);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
hbaBinEdg = [-1.2,-0.2,0.2,1.2];
hbaBinCtr = mean([hbaBinEdg(1:end-1);hbaBinEdg(2:end)]);

if ~exist('bfrm','var'),bfrm = cf(@(t,u) compute_bhv_ratemaps(t,u),Trials,unitsEgo);end;

lims = {[-300,300],[-300,300]};
lims = {[-200,200],[-200,200]};

orientationLabels = 'LCR';

% LOAD patch model 
rat = load_patch_model('rat');

mazeMaskFlag = true;                                      
pfsPlotMode  = 1;
reportMethod ='text';
pfsColorMap  = @jet;
displayModelFlag = true;

pfo = pfs;

% change incase first is empty
nBins = numel(hbaBinCtr);
nPhz = numel(binPhzc);



    
% SET figure layout
%[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','landscape',[],2.5,2.5,0.5,0.5);
[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','portrait',[],1.65,1.65,0.1,0.1);
[yind, yOffSet, xind, xOffSet] = deal(1,                         ... yind
                                      0,                         ... yOffSet
                                      1,                         ... xind
                                      0);                          % xOffSet

t = 20;
u = 3;

%for  t = 1:numel(Trials),
Trial = Trials{t};    
Trial.load('nq');


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

varBinEdg = hbaBinEdg;
varBinCtr = hbaBinCtr;
varBinInd = hbaBinInd;

hvec = xyz{t}(:,'nose',[1,2])-xyz{t}(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(Trial.meta.correction.headYaw),-sin(Trial.meta.correction.headYaw);...
                  sin(Trial.meta.correction.headYaw),cos(Trial.meta.correction.headYaw)],...
                 [2,3],...
                 [1,2]);

ind = [Trials{t}.stc{'w+p+n&t',sampleRate}];

ind.cast('TimeSeries');
ind.resample(xyz{t});


[mxr,mxp] = pft{t}.maxRate(unitsEgo{t}(u));        
pfstrj = MTADfet.encapsulate( ...
        Trials{t},...
        multiprod(bsxfun(@minus,...
                         mxp,...
                         sq(xyz{t}(:,'hcom',[1,2]))),...
                  hvec,2,[2,3]),...
        sampleRate,...
        'pfstrj','ptrj','t');

xBlockOffset = 0;
yBlockOffset = 0;

unit = unitsEgo{t}(u);

rateMinMax = max(max(cell2mat(cf(@(p) prctile(p.data.rateMap(:,p.data.clu==unit,1),99.9), pfo{t}))));

% THETA PLACE FIELD -------------------------------------------------------------------------------
[yind, yOffSet, xind, xOffSet] = deal(   1,       0,    1,       fig.subplot.width/2);        
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


% THETA ego phz -----------------------------------------------------------------------------------
for p = 1:nPhz,
[yind, yOffSet, xind, xOffSet] = deal(   p+1,       0,    1,       fig.subplot.width/2);        
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

for a = 1:nBins,
    [yind, yOffSet, xind, xOffSet] = deal(   1,       0,    a+1,       fig.subplot.width/2);        
    sax(end+1) = axes('Units','centimeters',                                                        ...
                      'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,                         ...
                                  fig.page.ypos(yind+yBlockOffset)+yOffSet,                         ...
                                  fig.subplot.width,                                                ...
                                  fig.subplot.height],                                              ...
                      'FontSize', 8,                                                                ...
                      'LineWidth',1);
    
    hold(sax(end),'on');
    indR = ind.data.*varBinInd==a;
    plot(pfstrj(indR,2),pfstrj(indR,1),'.','MarkerSize',1);
    xlim(sax(end),lims{1});
    ylim(sax(end),lims{2});        
    sax(end).XTick = [];
    sax(end).YTick = [];
end

% MAINBLOCK - egofields----------------------------------------------------------------------------

% PLOT EGO fields vs hba
for p = 1:nPhz,
    for a = 1:nBins,
        rmap = plot(pfo{t}{p,a},unit,1,[],[],false);        
        [yind, yOffSet, xind, xOffSet] = deal(   nPhz+2-p,       0,    a+1,     fig.subplot.width/2);
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
        
        Lines([],0,'w');
        Lines(0,[],'w');
        
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


% -------------------------------------------------------------------------------------------------

t = 6;
u = 7;
lims = {[-200,200],[-200,200]};
%for  t = 1:numel(Trials),
Trial = Trials{t};    
Trial.load('nq');
xyz{t} = preproc_xyz(Trial,sampleRate);

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

varBinEdg = hbaBinEdg;
varBinCtr = hbaBinCtr;
varBinInd = hbaBinInd;

hvec = xyz{t}(:,'nose',[1,2])-xyz{t}(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(Trial.meta.correction.headYaw),-sin(Trial.meta.correction.headYaw);...
                  sin(Trial.meta.correction.headYaw),cos(Trial.meta.correction.headYaw)],...
                 [2,3],...
                 [1,2]);

ind = [Trials{t}.stc{'w+p+n&t',sampleRate}];

ind.cast('TimeSeries');
ind.resample(xyz{t});


[mxr,mxp] = pft{t}.maxRate(unitsEgo{t}(u));        
pfstrj = MTADfet.encapsulate( ...
        Trials{t},...
        multiprod(bsxfun(@minus,...
                         mxp,...
                         sq(xyz{t}(:,'hcom',[1,2]))),...
                  hvec,2,[2,3]),...
        sampleRate,...
        'pfstrj','ptrj','t');

xBlockOffset = 0;
yBlockOffset = 5;

unit = unitsEgo{t}(u);

rateMinMax = max(max(cell2mat(cf(@(p) prctile(p.data.rateMap(:,p.data.clu==unit,1),99.9), pfo{t}))));

% THETA PLACE FIELD -------------------------------------------------------------------------------
[yind, yOffSet, xind, xOffSet] = deal(   1,       0,    1,       fig.subplot.width/2);        
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


% THETA ego phz -----------------------------------------------------------------------------------
for p = 1:nPhz,
[yind, yOffSet, xind, xOffSet] = deal(   p+1,       0,    1,       fig.subplot.width/2);        
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

for a = 1:nBins,
    [yind, yOffSet, xind, xOffSet] = deal(   1,       0,    a+1,       fig.subplot.width/2);        
    sax(end+1) = axes('Units','centimeters',                                                        ...
                      'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,                         ...
                                  fig.page.ypos(yind+yBlockOffset)+yOffSet,                         ...
                                  fig.subplot.width,                                                ...
                                  fig.subplot.height],                                              ...
                      'FontSize', 8,                                                                ...
                      'LineWidth',1);
    
    hold(sax(end),'on');
    indR = ind.data.*varBinInd==a;
    plot(pfstrj(indR,2),pfstrj(indR,1),'.','MarkerSize',1);
    xlim(sax(end),lims{1});
    ylim(sax(end),lims{2});        
    sax(end).XTick = [];
    sax(end).YTick = [];
end

% MAINBLOCK - egofields----------------------------------------------------------------------------

% PLOT EGO fields vs hba
for p = 1:nPhz,
    for a = 1:nBins,
        rmap = plot(pfo{t}{p,a},unit,1,[],[],false);        
        [yind, yOffSet, xind, xOffSet] = deal(   nPhz+2-p,       0,    a+1,     fig.subplot.width/2);
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
        
        Lines([],0,'w');
        Lines(0,[],'w');
        
        set(sax(end),'XTick',[]);
        set(sax(end),'YTick',[]);        
        
        % ADD subject
        if displayModelFlag,
            subject = struct(rat);
            subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
            subject = update_subject_patch(subject,'body', nBins+1-a,  true,hbaBinEdg,hbaBinCtr);
            patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
            patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
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


% -------------------------------------------------------------------------------------------------




ucounter = 1;
egoMeanRmapPos = [];
egoMaxRmapPos = [];
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
nrmap =rmap./sum(rmap(:),'omitnan');
rmapCenter = sq(sum(sum(bsxfun(@times,nrmap,mapPosition),'omitnan'),'omitnan'))';
egoMeanRmapPos(ucounter,p,a,:) = rmapCenter;
[~,maxPos] = max(nrmap(:));
if ~isempty(maxPos)
[maxX,maxY] = ind2sub(size(nrmap),maxPos);
egoMaxRmapPos(ucounter,p,a,:) = mapPosition(maxX,maxY,:);
else
egoMaxRmapPos(ucounter,p,a,:) = nan([1,1,1,2]);
end
end
end
ucounter = ucounter+1;
end
end


uids = [1:164];
uidsCA1 = [1:19,44:149];
uids = [1:19,44:149];
uidsER06 = [1:19];
uidsJG05 = [44:149];

uids = [20:41,123:127];
uidsCA3 = [20:41,123:127];
uidsED10 = [20:41];
uidsER06 = [123:127];
uidsJG05 = [150:164];
uids = [20:41,123:127,150:164];

hclr = 'bgr';


xBlockOffset =0;
yBlockOffset =0;

for p = 1:3
    [yind, yOffSet, xind, xOffSet] = deal(   nPhz+2-p,       0,    5,     fig.subplot.width/2+1);
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
                        fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
                        fig.subplot.width,                                ...
                        fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');
    plot(egoMeanRmapPos(uidsCA1,p,1,2)-egoMeanRmapPos(uidsCA1,p,2,2),...
         egoMeanRmapPos(uidsCA1,p,3,2)-egoMeanRmapPos(uidsCA1,p,2,2),...
         '.');
    scatter(mean(egoMeanRmapPos(uidsER06,p,1,2)-egoMeanRmapPos(uidsER06,p,2,2)),...
            mean(egoMeanRmapPos(uidsER06,p,3,2)-egoMeanRmapPos(uidsER06,p,2,2)),...
            10,...
            'c',...
            'Filled');
    scatter(mean(egoMeanRmapPos(uidsJG05,p,1,2)-egoMeanRmapPos(uidsJG05,p,2,2)),...
            mean(egoMeanRmapPos(uidsJG05,p,3,2)-egoMeanRmapPos(uidsJG05,p,2,2)),...
            10,...
            'r',...
            'Filled');
    line([-100,100],[100,-100],'Color','k')
    xlim([-100,100]);ylim([-100,100]);
    Lines([],0,'k');Lines(0,[],'k');    
    if p~=1,
        sax(end).XTick = [];
        sax(end).YTick = [];
    end
end


xBlockOffset =0;
yBlockOffset =5;
uidsCA3 = [20:41,123:127];
uidsED10 = [20:41];
uidsER06 = [123:127];
for p = 1:3
    [yind, yOffSet, xind, xOffSet] = deal(   nPhz+2-p,       0,    5,     fig.subplot.width/2+1);
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
                        fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
                        fig.subplot.width,                                ...
                        fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');
    plot(egoMeanRmapPos(uidsCA3,p,1,2)-egoMeanRmapPos(uidsCA3,p,2,2),...
         egoMeanRmapPos(uidsCA3,p,3,2)-egoMeanRmapPos(uidsCA3,p,2,2),...
         '.');
    scatter(mean(egoMeanRmapPos(uidsER06,p,1,2)-egoMeanRmapPos(uidsER06,p,2,2)),...
            mean(egoMeanRmapPos(uidsER06,p,3,2)-egoMeanRmapPos(uidsER06,p,2,2)),...
            10,...
            'c',...
            'Filled');
    scatter(mean(egoMeanRmapPos(uidsED10,p,1,2)-egoMeanRmapPos(uidsED10,p,2,2)),...
            mean(egoMeanRmapPos(uidsED10,p,3,2)-egoMeanRmapPos(uidsED10,p,2,2)),...
            10,...
            'g',...
            'Filled');
    line([-100,100],[100,-100],'Color','k')
    xlim([-100,100]);ylim([-100,100]);
    Lines([],0,'k');Lines(0,[],'k');    
    if p~=1,
        sax(end).XTick = [];
        sax(end).YTick = [];
    end
end



tind = [3,4,5,17,18,19,20,21,22,23,29];
%tind = [6,7,26,27,30];
sampleRate = 250;

global AP
% compute_ratemaps ---------------------------------------------------------------------------------
AP.compute_ratemaps =                                                                            ...
    struct('get_featureSet',            @fet_xy,                                                 ...
           'sampleRate',                16,                                                      ...
           'pfsArgs',                   struct('states',           'theta-groom-sit-rear',       ...
                                               'binDims',          [50,50],                      ...
                                               'SmoothingWeights', [2.4,2.4],                    ...
                                               'numIter',          1,                            ...
                                               'boundaryLimits',   [-500,500;-500,500],          ...
                                               'halfsample',       false)                        ...
           );
%---------------------------------------------------------------------------------------------------

dca = cf(@(T,U) accumulate_decoding_vars(T,U), Trials(tind),units(tind));

dca{7} = accumulate_decoding_vars(Trials{20},units{20});

%% HEAD BODY ANGLE ---------------------------------------------------------------------------------
clear('xcomp','ycomp','zcomp','ccomp');
xcomp.label = 'hba';            xcomp.data = [];    xcomp.edgs = [-1.25,-0.2,0.2,1.25];
ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5,6 );
ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-80,80];
fcomp.data = [];
mcomp.data = [];
%for t = 1:10
%for t = [2,3,4,5,7,8,9]
for t = [1:10]    
    dc = dca{t};
        
    mang = sq(dca{t}.xyz(:,'hcom',[1,2])-dca{t}.xyz(:,'bcom',[1,2]));
    mbang = atan2(mang(:,2),mang(:,1));
    bang =  sq(dca{t}.xyz(:,'bcom',[1,2]));
    bbang = atan2(bang(:,2),bang(:,1));
    mbbang = circ_dist(bbang,mbang);
    
    mind =  dc.stcm(:,1)==1                                       ... % theta state
           & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5)    ... % states {walk,turn,pause}
           & dc.hvfl(:,1)>2                                       ... % forward movement
           & dc.ucnt>=3 & dc.ucnt<9                               ... % number of coactive units
           & sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))<325;          ... % center of maze
            
    mind(mind==true) = randn([sum(mind),1])>0.5;                  ... % DROP half of the data points
                     
    xcomp.data = cat(1, xcomp.data, -dc.hbang(mind,1));
    ycomp.data = cat(1, ycomp.data, dc.phz(mind));
    %ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+0*double(t<5)-12.5*double(t>=5));    
    ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+25*double(t<4)-12.5*double(t>=4));
    fcomp.data = cat(1, fcomp.data, dc.esax(mind,1));
    mcomp.data = cat(1, mcomp.data, mbbang(mind));
end

[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,ccomp);
zmean   = zcomp.mean;    zmean(zcomp.count<5)  = nan;
zmedian = zcomp.median;  zmean(zcomp.count<5)  = nan;
zstd    = zcomp.std;     zstd(zcomp.count<5)   = nan;
zcount  = zcomp.count;   zcount(zcomp.count<5) = nan;
% $$$ set(figure(),'Units','centimeters','Position',[0,-3,8,28]);
% $$$ subplot(411); imagesc(xcomp.ctrs,ycomp.ctrs,zmean'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,['Mean ',ccomp.label]); colormap('jet');  caxis([ccomp.clim]);
% $$$     ylabel(ycomp.label);    
% $$$ subplot(412); imagesc(xcomp.ctrs,ycomp.ctrs,zstd'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,['Std ',ccomp.label]); 
% $$$     colormap('jet');
% $$$ subplot(413); 
% $$$ imagesc(xcomp.ctrs,ycomp.ctrs,zcount'); axis('xy');
% $$$     cax = colorbar(); ylabel(cax,'Count'); colormap('jet'); 
% $$$     xlabel(xcomp.label);
% $$$ subplot(414); 
% $$$     imagesc(xcomp.ctrs,ycomp.ctrs,bsxfun(@rdivide,zcount,sum(zcount,2,'omitnan'))'); axis('xy');    
% $$$     cax = colorbar(); ylabel(cax,'Probability'); colormap('jet'); 
% $$$     xlabel(xcomp.label);

xBlockOffset = 0;
yBlockOffset = 0;
figure
pinds = {ycomp.data > 0.5  & ycomp.data < 2.25,...
         ycomp.data > 2.25  & ycomp.data < 4.25,...
         ycomp.data > 4.25  & ycomp.data < 6};
for p = 1:numel(pinds)
% $$$    [yind, yOffSet, xind, xOffSet] = deal(   nPhz+2-p,       0,    6,     fig.subplot.width/2+1.5);
% $$$     sax(end+1) = axes('Units','centimeters',                                        ...
% $$$                       'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
% $$$                         fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
% $$$                         fig.subplot.width,                                ...
% $$$                         fig.subplot.height],                              ...
% $$$                       'FontSize', 8,                                                ...
% $$$                       'LineWidth',1);
% $$$     hold(sax(end),'on');
    sax = subplot(3,1,4-p)
    hold('on');
    
 
    indL =    pinds{p} & xcomp.data > -1.2 & xcomp.data < -0.2;
    indC =    pinds{p} & xcomp.data > -0.2 & xcomp.data < 0.2;
    indR =    pinds{p} & xcomp.data >  0.2 & xcomp.data < 1.2;    
    cdfplot(ccomp.data(indL)./10);
    cdfplot(ccomp.data(indC)./10);
    cdfplot(ccomp.data(indR)./10);
    xlim([-25,25]);    
    title('');
    %ylim([-25,30]);
    sax(end).XTick = [-20,-10,0,10,20];
    %sax(end).XTick = [-20,-10,0,10,20];
    if p~=1,
        %sax(end).XTickLabel = [];
        %sax(end).YTickLabel = [];
        sax.XTickLabel = [];
        sax.YTickLabel = [];
        xlabel('')
        ylabel('')
    end
end


meanDecPos = [];
for iter = 1:1000
    clear('xcomp','ycomp','zcomp','ccomp');
    xcomp.label = 'hba';            xcomp.data = [];    xcomp.edgs = [-1.25,-0.2,0.2,1.25];
    ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5,6 );
    ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-80,80];
    fcomp.data = [];
    for t = 1:10
        
        mang = sq(dca{t}.xyz(:,'hcom',[1,2])-dca{t}.xyz(:,'bcom',[1,2]));
        mbang = atan2(mang(:,2),mang(:,1));
        bang =  sq(dca{t}.xyz(:,'bcom',[1,2]));
        bbang = atan2(bang(:,2),bang(:,1));
        mbbang = circ_dist(bbang,mbang);

        dc = dca{t};
        mind =  dc.stcm(:,1)==1                                      ...
                & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5)  ...               
                & dc.hvfl(:,1)>2 ...
                & dc.ucnt>=3 & dc.ucnt<10 ...
                & sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))<325;
        if iter>=2
            mind(mind==true) = randn([sum(mind),1])>0.25;
        else
            mind(mind==true) = randn([sum(mind),1])>0;
        end
        xcomp.data = cat(1, xcomp.data, -dc.hbang(mind,1));
        ycomp.data = cat(1, ycomp.data, dc.phz(mind));
        %ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+0*double(t<5)-12.5*double(t>=5));    
        ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+25*double(t<4)-12.5*double(t>=4));
        fcomp.data = cat(1, fcomp.data, dc.esax(mind,1));
        mcomp.data = cat(1, mcomp.data, mbbang(mind,1));
    end

    [xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,ccomp);
    zmean   = zcomp.mean;    zmean(zcomp.count<5)  = nan;
    zmedian = zcomp.median;  zmean(zcomp.count<5)  = nan;
    zstd    = zcomp.std;     zstd(zcomp.count<5)   = nan;
    zcount  = zcomp.count;   zcount(zcomp.count<5) = nan;

    pinds = {ycomp.data > 0.5  & ycomp.data < 2.25,...
             ycomp.data > 2.25  & ycomp.data < 4.25,...
             ycomp.data > 4.25  & ycomp.data < 6};

    for p = 1:numel(pinds)
        indL =    pinds{p} & xcomp.data > -1.2 & xcomp.data < -0.2;
        indC =    pinds{p} & xcomp.data > -0.2 & xcomp.data < 0.2;
        indR =    pinds{p} & xcomp.data >  0.2 & xcomp.data < 1.2;    
        meanDecPos(iter,p,1,:) = [mean(ccomp.data(indL)),mean(fcomp.data(indL))-22];
        meanDecPos(iter,p,2,:) = [mean(ccomp.data(indC)),mean(fcomp.data(indC))-22];
        meanDecPos(iter,p,3,:) = [mean(ccomp.data(indR)),mean(fcomp.data(indR))-22];
    end
    disp(['[INFO] iteration: ', num2str(iter),' of 1000']);
end
save(fullfile(Trials{20}.path.project,'analysis','EgoProCode2D_hba_decode.mat'), 'meanDecPos');



figure();
%plot(sq(meanDecPos(:,:,1,2)),'g-+')
histogram(meanDecPos(1,1,1,2)-meanDecPos(2:end,1,1,2),linspace(-4,4,100))


hvfl = fet_href_HXY(Trials{20},sampleRate,[],'trb');

figure();
hold('on');
for a = 1:3
plot(sq(meanDecPos(:,1,a,1)),sq(meanDecPos(:,1,a,2)),'g-+')
plot(sq(meanDecPos(:,2,a,1)),sq(meanDecPos(:,2,a,2)),'b-+')
plot(sq(meanDecPos(:,3,a,1)),sq(meanDecPos(:,3,a,2)),'r-+')
end
xlim([-60,60]);
ylim([-30,60]);






%% HEAD BODY ANGLE ---------------------------------------------------------------------------------
medD = [];
skwD = [];
stdD = [];
medR = [];
skwR = [];
stdR = [];
medL = [];
skwL = [];
stdL = [];
medC = [];
skwC = [];
stdC = [];

rdists = 0:5:30;
for r = 1:numel(rdists)
clear('xcomp','ycomp','zcomp','ccomp');
xcomp.label = 'hba';            xcomp.data = [];    xcomp.edgs = [-1.25,-0.2,0.2,1.25];
ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5,6 );
ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-80,80];
fcomp.data = [];
mcomp.data = [];
%for t = 1:10
for t = [2,3,4,5,7,8,9]
%for t = [6]    
    dc = dca{t};
        
    mang = sq(dca{t}.xyz(:,'nose',[1,2])-dca{t}.xyz(:,'hcom',[1,2]));
    mbang = atan2(mang(:,2),mang(:,1));
    bang =  sq(dca{t}.xyz(:,'hcom',[1,2]));
    bbang = atan2(bang(:,2),bang(:,1));
    mbbang = circ_dist(bbang,mbang);
    
    mind =  dc.stcm(:,1)==1                                      ...
           & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5)  ...               
           & dc.hvfl(:,1)>2 ...
           & dc.ucnt>=2 & dc.ucnt<10 ...
           & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)<rdists(r)+10 ...
           & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)>rdists(r);
% $$$                        & sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))>200;

    mind(mind==true) = randn([sum(mind),1])>0.5;
    xcomp.data = cat(1, xcomp.data, -dc.hbang(mind,1));
    ycomp.data = cat(1, ycomp.data, dc.phz(mind));
    %ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+0*double(t<5)-12.5*double(t>=5));    
    %ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+25*double(t<4)-12.5*double(t>=4));
    ccomp.data = cat(1, ccomp.data, dc.esax(mind,2));
    fcomp.data = cat(1, fcomp.data, dc.esax(mind,1));
    mcomp.data = cat(1, mcomp.data, mbbang(mind));
end

[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,ccomp);
zmean   = zcomp.mean;    zmean(zcomp.count<5)  = nan;
zmedian = zcomp.median;  zmean(zcomp.count<5)  = nan;
zstd    = zcomp.std;     zstd(zcomp.count<5)   = nan;
zcount  = zcomp.count;   zcount(zcomp.count<5) = nan;
% $$$ 
% $$$ 
% $$$ figure
% $$$ normType = 'probability';
% $$$ %normType = 'count';
% $$$ sectors = linspace(-pi,pi,13);
% $$$ medD = [];
% $$$ skwD = [];
% $$$ medR = [];
% $$$ skwR = [];
% $$$ medL = [];
% $$$ skwL = [];
% $$$ medC = [];
% $$$ skwC = [];
% $$$ for s = 1:numel(sectors)-1
% $$$     subplot(4,4,s);
% $$$     hold('on');
% $$$ bind = mcomp.data<sectors(s+1) & mcomp.data>sectors(s) & ycomp.data>4& ycomp.data<6;
% $$$ ind = bind & xcomp.data>0.2;
% $$$ histogram(ccomp.data(ind,1),linspace(-400,400,50),'Normalization',normType);
% $$$ medR(s) = median(ccomp.data(ind,1));
% $$$ skwR(s) = skew(ccomp.data(ind,1));
% $$$ ind = bind & xcomp.data<0.2& xcomp.data>-0.2;
% $$$ histogram(ccomp.data(ind,1),linspace(-400,400,50),'Normalization',normType);
% $$$ medC(s) = median(ccomp.data(ind,1));
% $$$ skwC(s) = skew(ccomp.data(ind,1));
% $$$ 
% $$$ ind = bind & xcomp.data<-0.2;
% $$$ histogram(ccomp.data(ind,1),linspace(-400,400,50),'Normalization',normType);title('Towards Center');
% $$$ medL(s) = median(ccomp.data(ind,1));
% $$$ skwL(s) = skew(ccomp.data(ind,1));
% $$$ medD(s) = medR(s)-medL(s);
% $$$ skwD(s) = skwR(s);
% $$$ title([num2str(sectors(s)),'-',num2str(sectors(s+1))]);
% $$$ end

sectors = linspace(-pi,pi,13);
for s = 1:numel(sectors)-1
    bind = mcomp.data<sectors(s+1) & mcomp.data>sectors(s) & ycomp.data>4& ycomp.data<6;
    ind = bind & xcomp.data>0.2;
    medR(r,s) = median(ccomp.data(ind,1)./10);%-rmodel(r,s);
    skwR(r,s) = skew(ccomp.data(ind,1)./10);
    stdR(r,s) = std(ccomp.data(ind,1)./10);    
    ind = bind & xcomp.data<0.2& xcomp.data>-0.2;
    medC(r,s) = median(ccomp.data(ind,1)./10);%-rmodel(r,s);
    skwC(r,s) = skew(ccomp.data(ind,1)./10);
    stdC(r,s) = std(ccomp.data(ind,1)./10);    
    ind = bind & xcomp.data<-0.2;
    medL(r,s) = median(ccomp.data(ind,1)./10);%-rmodel(r,s);
    skwL(r,s) = skew(ccomp.data(ind,1)./10);
    stdL(r,s) = std(ccomp.data(ind,1)./10);    
    medD(r,s) = medR(r,s)-medL(r,s);
    skwD(r,s) = skwR(r,s)-skwL(r,s);
    stdD(r,s) = stdR(r,s)-stdL(r,s);
end
end

sectorc = mean([sectors(2:end);sectors(1:end-1)]);

% $$$ figure,imagesc(medR')


[THETA,RR] = meshgrid(sectors,[rdists,35]);
% $$$ 
% $$$ THETA = cat(2,THETA(:,end),THETA);
% $$$ THETA = cat(1,THETA(end,:),THETA);
% $$$ 
% $$$ RR = cat(2,RR(:,end),RR);
% $$$ RR = cat(1,RR(end,:),RR);


[A,B] = pol2cart(THETA,RR);

% SUPFIG - Non-Uniformity of decoded lateral position: subject-arena orientation and distance to maze center
figure
cmedD = medL; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,1);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-7.5,7.5]);
ylabel(cax,'cm');
title({'Leftward','Median'});
cmedD = stdL; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,1);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([0,20]);
ylabel(cax,'Std')
cmedD = skwL; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,1);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-2,2]);
ylabel(cax,'Skew')

cmedD = medC; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,2);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-7.5,7.5]);
ylabel(cax,'cm');
title({'Centered','Median'})
cmedD = stdC; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,2);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([0,20]);
ylabel(cax,'Std');
cmedD = skwC; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,2);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-2,2]);
title('Skew')

cmedD = medR; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,3);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-7.5,7.5])
ylabel(cax,'cm');
title({'Rightward','Median'});
cmedD = stdR; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,3);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([0,20]);
ylabel(cax,'Std');
cmedD = skwR; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,3);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-2,2]);
ylabel(cax,'Skew');
       
cmedD = medD; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,4);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-7.5,7.5]);
title({'Rightward-Leftward','\DeltaMedian'})
ylabel(cax,'cm');
cmedD = stdD; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,4);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([0,5]);
ylabel(cax,'\DeltaStd');
cmedD = skwD; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,4);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-1,1]);
ylabel(cax,'\DeltaSkew');

ForAllSubplots('daspect([1,1,1]);');
ForAllSubplots('ylim([-40,40]);');
ForAllSubplots('xlim([-40,40]);');
fax = axes('Position',[0,0,1,1],'Visible','off');
text(fax,0.1,0.9,{'Non-Uniformity of decoded lateral position ',...
                  'depedent on subject-arena orientation and distance to maze center'});
set(gcf(),'PaperOrientation','landscape');


rdistc = rdists+50;

figure,plot(mean(medR))
figure,plot(medR')
figure,
plot(rmodel')

mvar = [];
ms = [100:10:10000];
for m = 1:numel(ms)
rmodel = [repmat(sin(sectorc-pi/4)./ms(m),[numel(rdistc),1]).*repmat(rdistc'.^2,[1,numel(sectorc)])];
mvar(m) = var((medL(:)-rmodel(:)));
end
figure,plot(mvar)
[~,m] = min(mvar);
figure,plot([repmat(sin(sectorc)./ms(m),[numel(rdistc),1]).*repmat(rdistc'.^2,[1,numel(sectorc)])]');

figure,plot(rmodel')
figure,
subplot(321);plot(medR');
subplot(322);histogram((medR-medC)',linspace(-100,100,20));
subplot(323);plot(medL');
subplot(324);histogram((medL-medC)',linspace(-100,100,20));
subplot(325);histogram((medL-medR)',linspace(-100,100,20));

figure,plot(medR')
figure,plot([repmat(sin(sectorc)./ms(m),[numel(rdistc),1]).*repmat(rdistc'.^2,[1,numel(sectorc)])]');
figure,plot((medL-rmodel)')
figure,plot((medR-rmodel)')

figure,imagesc(medD')
figure,hist(medD(:))
figure;hold('on');
plot(sectorc,medR-rmodel)
plot(sectorc,medC-rmodel)
plot(sectorc,medL-rmodel)

figure;plot(sectorc,medD);

figure,
subplot(211);polar(sectorc([1:end,1]),medD([1:end,1]))
subplot(212);polar(sectorc([1:end,1]),skwD([1:end,1]))

normType = 'probability';
normType = 'count';
figure,
subplot(2,2,1);hold('on'); bind = mcomp.data<pi/8 & mcomp.data>-pi/8 & ycomp.data>4& ycomp.data<6;
ind = bind & xcomp.data>0.2;histogram(ccomp.data(ind,1),linspace(-600,600,100),'Normalization',normType);ind = bind & xcomp.data<0.2& xcomp.data>-0.2;histogram(ccomp.data(ind,1),linspace(-600,600,100),'Normalization',normType);ind = bind & xcomp.data<-0.2;histogram(ccomp.data(ind,1),linspace(-600,600,100),'Normalization',normType);title('Towards Center');
subplot(2,2,2);hold('on'); bind = mcomp.data>7*pi/8 | mcomp.data<-7*pi/8& ycomp.data>4& ycomp.data<6;
ind = bind & xcomp.data>0.2;histogram(ccomp.data(ind,1),linspace(-600,600,100),'Normalization',normType);ind = bind & xcomp.data<0.2& xcomp.data>-0.2;histogram(ccomp.data(ind,1),linspace(-600,600,100),'Normalization',normType);ind = bind & xcomp.data<-0.2;histogram(ccomp.data(ind,1),linspace(-600,600,100),'Normalization',normType);title('Away Center');
subplot(2,2,3); hold('on'); bind = mcomp.data<5*pi/8 & mcomp.data>3*pi/8& ycomp.data>4& ycomp.data<6;
ind = bind & xcomp.data>0.2;histogram(ccomp.data(ind,1),linspace(-600,600,100),'Normalization',normType);ind = bind & xcomp.data<0.2& xcomp.data>-0.2;histogram(ccomp.data(ind,1),linspace(-600,600,100),'Normalization',normType);ind = bind & xcomp.data<-0.2;histogram(ccomp.data(ind,1),linspace(-600,600,100),'Normalization',normType);title('CW Center');
subplot(2,2,4); hold('on'); bind = mcomp.data>-5*pi/8 & mcomp.data<-3*pi/8& ycomp.data>4& ycomp.data<6;
ind = bind & xcomp.data>0.2;histogram(ccomp.data(ind,1),linspace(-600,600,100),'Normalization',normType);ind = bind & xcomp.data<0.2& xcomp.data>-0.2;histogram(ccomp.data(ind,1),linspace(-600,600,100),'Normalization',normType);ind = bind & xcomp.data<-0.2;histogram(ccomp.data(ind,1),linspace(-600,600,100),'Normalization',normType);title('CCW Center');

figure,
subplot(2,2,1);
hold('on');
bind = mcomp.data<pi/4 & mcomp.data>-pi/4& ycomp.data>4& ycomp.data<6;
ind = bind & xcomp.data>0.2;
cdfplot(ccomp.data(ind,1))
ind = bind & xcomp.data<0.2& xcomp.data>-0.2;
cdfplot(ccomp.data(ind,1))
ind = bind & xcomp.data<-0.2;
cdfplot(ccomp.data(ind,1))
xlim([-400,400]); 
title('towards center')
subplot(2,2,2);
hold('on');
bind = mcomp.data>3*pi/4 & mcomp.data>-3*pi/4& ycomp.data>4& ycomp.data<6;
ind = bind & xcomp.data>0.2;
cdfplot(ccomp.data(ind,1))
ind = bind & xcomp.data<0.2& xcomp.data>-0.2;
cdfplot(ccomp.data(ind,1))
ind = bind & xcomp.data<-0.2;
cdfplot(ccomp.data(ind,1))
xlim([-400,400]); 
title('away center')
subplot(2,2,4);
hold('on');
bind = mcomp.data>-3*pi/4 & mcomp.data<-pi/4& ycomp.data>4& ycomp.data<6;
ind = bind & xcomp.data>0.2;
cdfplot(ccomp.data(ind,1))
ind = bind & xcomp.data<0.2& xcomp.data>-0.2;
cdfplot(ccomp.data(ind,1))
ind = bind & xcomp.data<-0.2;
cdfplot(ccomp.data(ind,1))
xlim([-400,400]); 
title('CCW')
subplot(2,2,3);
hold('on');
bind = mcomp.data<3*pi/4 & mcomp.data>pi/4& ycomp.data>4& ycomp.data<6;
ind = bind & xcomp.data>0.2;
cdfplot(ccomp.data(ind,1)); skew(ccomp.data(ind,1))
ind = bind & xcomp.data<0.2& xcomp.data>-0.2;
cdfplot(ccomp.data(ind,1)); skew(ccomp.data(ind,1))
ind = bind & xcomp.data<-0.2;
cdfplot(ccomp.data(ind,1)); skew(ccomp.data(ind,1))
xlim([-400,400]); 
title('CW')


figure,
hold('on');
bind = mcomp.data<pi/4 & mcomp.data>-pi/4& ycomp.data>4& ycomp.data<6;
ind = bind & xcomp.data>0.2;
cdfplot(ccomp.data(ind,1))
ind = bind & xcomp.data<0.2& xcomp.data>-0.2;
cdfplot(ccomp.data(ind,1))
ind = bind & xcomp.data<-0.2;
cdfplot(ccomp.data(ind,1))
xlim([-400,400]); 
title('CW')

figure,
hold('on');
% $$$ bind = dca{tid}.stcm(:,1)==1 ...
% $$$        & (dca{tid}.stcm(:,3)==3 | dca{tid}.stcm(:,4)==4| dca{tid}.stcm(:,5)==5) ...
% $$$        & dca{tid}.ucnt>2 & dca{tid}.phz>4 & dca{tid}.phz<6 & (mbbang>-pi/2&mbbang<pi/2) & dca{tid}.hdist>350; 
% $$$ bind = dca{tid}.stcm(:,1)==1 ...
% $$$        & (dca{tid}.stcm(:,3)==3 | dca{tid}.stcm(:,4)==4| dca{tid}.stcm(:,5)==5) ...
% $$$        & dca{tid}.ucnt>2 & dca{tid}.phz>4 & dca{tid}.phz<6 & (mbbang<-pi/2|mbbang>pi/2) & dca{tid}.hdist<300; 
bind = dca{tid}.stcm(:,1)==1                                                    ...
       & (dca{tid}.stcm(:,3)==3 | dca{tid}.stcm(:,4)==4| dca{tid}.stcm(:,5)==5) ...
       & dca{tid}.ucnt>2 & dca{tid}.phz>4 & dca{tid}.phz<6 & (mbbang<pi/4 & mbbang>-pi/4) & dca{tid}.hdist<350; 
% $$$ bind = dca{tid}.stcm(:,1)==1                                                    ...
% $$$        & (dca{tid}.stcm(:,3)==3 | dca{tid}.stcm(:,4)==4| dca{tid}.stcm(:,5)==5) ...
% $$$        & dca{tid}.ucnt>2 & dca{tid}.phz>4 & dca{tid}.phz<6 & (mbbang<pi/4 & mbbang>-pi/4) & dca{tid}.hdist<300; 
ind = bind & dca{tid}.hbang.data>0.2;
histogram(dca{tid}.esax(ind,2),linspace(-300,300,100));
ind = bind & dca{tid}.hbang.data<0.2& dca{tid}.hbang.data>-0.2;
histogram(dca{tid}.esax(ind,2),linspace(-300,300,100));
ind = bind & dca{tid}.hbang.data<-0.2;
histogram(dca{tid}.esax(ind,2),linspace(-300,300,100));


tid = 7;

figure,
hold('on');
bind = dca{tid}.stcm(:,1)==1 ...
       & (dca{tid}.stcm(:,3)==3 | dca{tid}.stcm(:,4)==4| dca{tid}.stcm(:,5)==5) ...
       & dca{tid}.ucnt>2 & dca{tid}.phz>4 & dca{tid}.phz<6 & dca{tid}.hdist<325; 
ind = bind & dca{tid}.hbang.data>0.2;
histogram(mbang(ind,1),linspace(-pi,pi,100));
ind = bind & dca{tid}.hbang.data<0.2& dca{tid}.hbang.data>-0.2;
histogram(mbang(ind,1),linspace(-pi,pi,100));
ind = bind & dca{tid}.hbang.data<-0.2;
histogram(mbang(ind,1),linspace(-pi,pi,100));


figure,
hold('on');
plot(hvfl.data(:,2));
plot(dlat.data(:,1));
plot(dlat.data(:,2).*100);
Lines([],0,'k');



for p = 1:numel(pinds)
subplot2(3,3,2,p);
    hold('on');
    indL =    pinds{p} & xcomp.data > -1.2 & xcomp.data < -0.2;
    indC =    pinds{p} & xcomp.data > -0.2 & xcomp.data < 0.2;
    indR =    pinds{p} & xcomp.data >  0.2 & xcomp.data < 1.2;    
    
    histogram(ccomp.data(indL),linspace(-400,400,50),'Normalization','probability');
    histogram(ccomp.data(indC),linspace(-400,400,50),'Normalization','probability');
    histogram(ccomp.data(indR),linspace(-400,400,50),'Normalization','probability');
    
    xlim([-300,300]);
    ylim([0,0.1]);
    
    Lines(mean(ccomp.data(indL)),[],'b');    
    Lines(mean(ccomp.data(indC)),[],'g');
    Lines(mean(ccomp.data(indR)),[],'r');
    
    indA = indR|indC|indL;
    [B,BINT,R,RINT,RSTATS] = regress(ccomp.data(indA),[ones([sum(indA),1]),xcomp.data(indA)]);
% RSTATS R-square statistic,  F statistic and p value, error variance
subplot2(3,3,1,p);
    P = polyfit(xcomp.data(indA),ccomp.data(indA),1);
    hold('on');
    hist2([ccomp.data(indA),xcomp.data(indA)],linspace(-300,300,30),linspace(-1.25,1.25,12),'xprob');
    axis('tight');
    line(polyval(P,[-1.25,1.25]),[-1.25,1.25],'Color','m');
subplot2(3,3,3,p);
    hold('on');
    cdfplot(ccomp.data(indL))
    cdfplot(ccomp.data(indC))
    cdfplot(ccomp.data(indR))
    xlim([-300,300]);
RSTATS    
end
