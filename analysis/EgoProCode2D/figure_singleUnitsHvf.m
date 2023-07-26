% $$$ configure_default_args();
% $$$ EgoProCode2D_load_data();

% binPhzc
binPhzs = linspace(0.5,2*pi-0.5,4);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
hbvBinEdg = [-25,-5,5,25,80];
hbvBinCtr = mean([hbvBinEdg(1:end-1);hbvBinEdg(2:end)]);

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
displayModelFlag = false;

pfo = pfv;

% change incase first is empty
nBins = numel(hbvBinCtr);
nPhz = numel(binPhzc);


% SET figure layout
%[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','landscape',[],2.5,2.5,0.5,0.5);
[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','portrait',[],1.65,1.65,0.1,0.1);
[yind, yOffSet, xind, xOffSet] = deal(1,                         ... yind
                                      0,                         ... yOffSet
                                      1,                         ... xind
                                      0);                          % xOffSet

t = 18;
u = 2;

%for  t = 1:numel(Trials),
Trial = Trials{t};    
Trial.load('nq');


% HVF setup
hvfl = fet_href_HXY(Trials{t},sampleRate,false,'trb',4);
hvfBinInd = discretize(hvfl(:,1), hvfBinEdg);


varBinEdg = hvfBinEdg;
varBinCtr = hvfBinCtr;
varBinInd = hvfBinInd;

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


[mxr,mxp] = pft{t}.maxRate(unitsEgoHvf{t}(u));        
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

unit = unitsEgoHvf{t}(u);

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

for a = 1:nBins-1,
    [yind, yOffSet, xind, xOffSet] = deal(   1,       0,    a+1,       fig.subplot.width/2);        
    sax(end+1) = axes('Units','centimeters',                                                        ...
                      'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,                         ...
                                  fig.page.ypos(yind+yBlockOffset)+yOffSet,                         ...
                                  fig.subplot.width,                                                ...
                                  fig.subplot.height],                                              ...
                      'FontSize', 8,                                                                ...
                      'LineWidth',1);
    
    hold(sax(end),'on');
    indR = ind.data.*varBinInd==(a+1);
    plot(pfstrj(indR,2),pfstrj(indR,1),'.','MarkerSize',1);
    xlim(sax(end),lims{1});
    ylim(sax(end),lims{2});        
    sax(end).XTick = [];
    sax(end).YTick = [];
end

% MAINBLOCK - egofields----------------------------------------------------------------------------

% PLOT EGO fields vs hba
for p = 1:nPhz,
    for a = 1:nBins-1,
        rmap = plot(pfo{t}{p,a+1},unit,1,[],[],false);        
        [yind, yOffSet, xind, xOffSet] = deal(   nPhz+2-p,       0,    a+1,     fig.subplot.width/2);
        sax(end+1) = axes('Units','centimeters',                                        ...
                          'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
                                      fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
                                      fig.subplot.width,                                ...
                                      fig.subplot.height],                              ...
                          'FontSize', 8,                                                ...
                          'LineWidth',1);
        hold(sax(end),'on');
        
        pcolor(pfo{t}{p,a+1}.adata.bins{1},...
               pfo{t}{p,a+1}.adata.bins{2},...
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
            subject = update_subject_patch(subject,'head',[], false,[-1.2,1.2],0);
            subject = update_subject_patch(subject,'body', 1,  true,[-1.2,1.2],0);
            patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
            patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
            patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
            line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
            line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
        else
% $$$             if p == 1,
% $$$                 xlabel([num2str(varBinCtr(a)),' mm/s']);
% $$$             end
        end;%if displayModelFlag
    end;%for a
end;%for p
 

% -------------------------------------------------------------------------------------------------

t = 6;
u = 3;
lims = {[-200,200],[-200,200]};
%for  t = 1:numel(Trials),
Trial = Trials{t};    
Trial.load('nq');
%xyz{t} = preproc_xyz(Trial,sampleRate);

% HVF setup
hvfl = fet_href_HXY(Trials{t},sampleRate,false,'trb',4);
hvfBinInd = discretize(hvfl(:,1), hvfBinEdg);


varBinEdg = hvfBinEdg;
varBinCtr = hvfBinCtr;
varBinInd = hvfBinInd;

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
ind.resample(hvfl);


[mxr,mxp] = pft{t}.maxRate(unitsEgoHvf{t}(u));        
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

unit = unitsEgoHvf{t}(u);

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

for a = 1:nBins-1,
    [yind, yOffSet, xind, xOffSet] = deal(   1,       0,    a+1,       fig.subplot.width/2);        
    sax(end+1) = axes('Units','centimeters',                                                        ...
                      'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,                         ...
                                  fig.page.ypos(yind+yBlockOffset)+yOffSet,                         ...
                                  fig.subplot.width,                                                ...
                                  fig.subplot.height],                                              ...
                      'FontSize', 8,                                                                ...
                      'LineWidth',1);
    
    hold(sax(end),'on');
    indR = ind.data.*varBinInd==(a+1);
    plot(pfstrj(indR,2),pfstrj(indR,1),'.','MarkerSize',1);
    xlim(sax(end),lims{1});
    ylim(sax(end),lims{2});        
    sax(end).XTick = [];
    sax(end).YTick = [];
end

% MAINBLOCK - egofields----------------------------------------------------------------------------

% PLOT EGO fields vs hba
for p = 1:nPhz,
    for a = 1:nBins-1,
        rmap = plot(pfo{t}{p,a+1},unit,1,[],[],false);        
        [yind, yOffSet, xind, xOffSet] = deal(   nPhz+2-p,       0,    a+1,     fig.subplot.width/2);
        sax(end+1) = axes('Units','centimeters',                                        ...
                          'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
                                      fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
                                      fig.subplot.width,                                ...
                                      fig.subplot.height],                              ...
                          'FontSize', 8,                                                ...
                          'LineWidth',1);
        hold(sax(end),'on');
        
        pcolor(pfo{t}{p,a+1}.adata.bins{1},...
               pfo{t}{p,a+1}.adata.bins{2},...
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
                xlabel([num2str(varBinCtr(a+1)),' mm/s']);
            end
        end;%if displayModelFlag
    end;%for a
end;%for p





% er01-20110719
unitsEgoHvf{1} = [15,42,99];
% er01-20110721
unitsEgoHvf{2} = [23,75,86,99];
% ER06-20130612
unitsEgoHvf{3} = [38,61,80,151,158];
% ER06-20130613
unitsEgoHvf{4} = [23,33,38,49,51,69,94,107,117,126,140,145,154,166,175];
% ER06-20130614
unitsEgoHvf{5} = [15,34,35,38,40,44,86,90,99,112,121,159];
% Ed10-20140816
unitsEgoHvf{6} = [4,7,10,13,18,25,35,37,49,89];
% Ed10-20140817
unitsEgoHvf{7} = [1,10,17,23,33,38,57,63,64,66,73,82,99,104,105,108];% 2
% jg04-20120128
unitsEgoHvf{8} = [10];
% jg04-20120129
unitsEgoHvf{9} = [39];
% jg04-20120130
unitsEgoHvf{10} = [];
% jg04-20120131
unitsEgoHvf{11} = [];
% jg04-20120201
unitsEgoHvf{12} = [];
% jg04-20120210
unitsEgoHvf{13} = [];
% jg04-20120211
unitsEgoHvf{14} = [10,17,20];
% jg04-20120212
unitsEgoHvf{15} = [];
% jg04-20120213
unitsEgoHvf{16} = [5];
% jg05-20120309
unitsEgoHvf{17} = [70];
% jg05-20120310
unitsEgoHvf{18} = [9,11,17,18,19,20,21,22,24,25,29,33,34,35,42,...
                  44,49,52,56,60,61,72,74,75,78,81,82,84];
% jg05-20120311
unitsEgoHvf{19} = [10,13,15,27,63,66,67,97,143,150,152,153,160];
% jg05-20120312
unitsEgoHvf{20} = [20,21,25,31,35,41,44,52,59,61,72,79,80,81,85,86,89,103,...
                   105,110,111,119,129,139,141,144,151];
% jg05-20120315
unitsEgoHvf{21} = [6,22,24,25,37,43,61,63,68,77,97];
% jg05-20120316
unitsEgoHvf{22} = [8,11,13,18,19,22,30,41,42,43,48,50,56,58,61,65];
% jg05-20120317
unitsEgoHvf{23} = [9,12,18,29,31,36,40,50,51,54,56,63,69,70,72];
% jg05-20120323
unitsEgoHvf{24} = [22,26,32];
% jg05-20120324
unitsEgoHvf{25} = [10,29];
% ER06-20130624
unitsEgoHvf{26} = [82,140,175,218];
% Ed10-20140815
unitsEgoHvf{27} = [6,10,39,51,92,98,100,102];
% er01-20110722
unitsEgoHvf{28} = [15,45,74,90];
% FS03-2020122
unitsEgoHvf{29} = [11,12,24,27,33,43,63,64,65,69,70,71,76,116,120,126];
% jg05-20120329
unitsEgoHvf{30} = [20,23,29,30,55,56,63,82,83,84,97,102,107];
                   
                   



% $$$ figure();

ucounter = 1;
egoMeanRmapPos = [];
egoMaxRmapPos = [];
for t = 1:numel(Trials)
for u = 1:numel(unitsEgoHvf{t})
unit = unitsEgoHvf{t}(u);
for p = 1:3;
for a = 2:4,
binSubsetX = abs(pfv{t}{p,a}.adata.bins{1})<300;
binSubsetY = abs(pfv{t}{p,a}.adata.bins{2})<300;
mapPosition = cell([1,2]);
[mapPosition{:}] = ndgrid(pfv{t}{p,a}.adata.bins{1}(binSubsetX),...
                          pfv{t}{p,a}.adata.bins{2}(binSubsetY));
mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
rmap = plot(pfv{t}{p,a},unit,[],[],[],false);
rmap = rmap(binSubsetX,binSubsetY);
nanmap = double(~isnan(rmap));
nanmap(nanmap==0) = nan;
%rmap = rmap.*fliplr(nanmap);
rmap(rmap<3) = 0;
% $$$ subplot2(3,3,4-p,a);
% $$$ hold('on');
% $$$ imagescnan({pfs{t}{p,a}.adata.bins{1}(binSubsetX),...
% $$$                    pfs{t}{p,a}.adata.bins{2}(binSubsetY),...
% $$$                    rmap'});
nrmap =rmap./sum(rmap(:),'omitnan');
rmapCenter = sq(sum(sum(bsxfun(@times,nrmap,mapPosition),'omitnan'),'omitnan'))';
% $$$ plot(rmapCenter(1),rmapCenter(2),'*m')
% $$$ axis('tight')
% $$$ Lines([],0,'g');
% $$$ Lines(0,[],'g');
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

unitsEgoHvfMap =[];
for c = 1:numel(unitsEgoHvf)
    unitsEgoHvfMap = [unitsEgoHvfMap;ones([numel(unitsEgoHvf{c}),1]).*c];
end

uids = ismember(unitsEgoHvfMap,[3,4,5,17:25,29]);
uidsCA1 = ismember(unitsEgoHvfMap,[3,4,5,17:25]);
uidsER06 = ismember(unitsEgoHvfMap,[3,4,5]);
uidsJG05 = ismember(unitsEgoHvfMap,[17:25]);
uidsCA3 = ismember(unitsEgoHvfMap,[1,2,6,7,15,16,26,30]);

hclr = 'ckm';
figure,violin

figure,
e1 = egoMeanRmapPos(uidsCA1,1,3,1);
e2 = egoMeanRmapPos(uidsCA1,1,4,1);
    plot(sign(e1).*sqrt(abs(e1)+1),...
         sign(e2).*sqrt(abs(e2)+1),...
         '.');
line([-15,15],[-15,15]);
figure,plot(sq(egoMeanRmapPos(uidsCA1,1,[2:4],1))');

figure,
subplot(211);histogram(egoMeanRmapPos(uidsCA1,1,3,1)-egoMeanRmapPos(uidsCA1,1,2,1),linspace(-150,150,30))
subplot(212);histogram(egoMeanRmapPos(uidsCA1,1,4,1)-egoMeanRmapPos(uidsCA1,1,3,1),linspace(-150,150,30))



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
    %figure,hold('on');
    plot(egoMeanRmapPos(uidsCA1,p,3,1)-egoMeanRmapPos(uidsCA1,p,2,1),...
         egoMeanRmapPos(uidsCA1,p,4,1)-egoMeanRmapPos(uidsCA1,p,2,1),...
         '.');
% $$$     scatter(mean(egoMeanRmapPos(uidsER06,p,3,1)-egoMeanRmapPos(uidsER06,p,2,1)),...
% $$$             mean(egoMeanRmapPos(uidsER06,p,4,1)-egoMeanRmapPos(uidsER06,p,2,1)),...
% $$$             10,...
% $$$             'c',...
% $$$             'Filled');
% $$$     scatter(mean(egoMeanRmapPos(uidsJG05,p,3,1)-egoMeanRmapPos(uidsJG05,p,2,1)),...
% $$$             mean(egoMeanRmapPos(uidsJG05,p,4,1)-egoMeanRmapPos(uidsJG05,p,2,1)),...
% $$$             10,...
% $$$             'r',...
% $$$             'Filled');
    %line([-150,150],[150,-150],'Color','k')
    xlim([-150,150]);ylim([-150,150]);
    Lines([],0,'k');Lines(0,[],'k');    
    if p~=1,
        sax(end).XTick = [];
        sax(end).YTick = [];
    end
end



xBlockOffset =0;
yBlockOffset =5;
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
    plot(egoMeanRmapPos(uidsCA3,p,3,1)-egoMeanRmapPos(uidsCA3,p,2,1),...
         egoMeanRmapPos(uidsCA3,p,4,1)-egoMeanRmapPos(uidsCA3,p,2,1),...
         '.');
% $$$     scatter(mean(egoMeanRmapPos(uidsER06,p,3,1)-egoMeanRmapPos(uidsER06,p,2,1)),...
% $$$             mean(egoMeanRmapPos(uidsER06,p,4,1)-egoMeanRmapPos(uidsER06,p,2,1)),...
% $$$             10,...
% $$$             'c',...
% $$$             'Filled');
% $$$     scatter(mean(egoMeanRmapPos(uidsED10,p,1,2)-egoMeanRmapPos(uidsED10,p,2,2)),...
% $$$             mean(egoMeanRmapPos(uidsED10,p,3,2)-egoMeanRmapPos(uidsED10,p,2,2)),...
% $$$             10,...
% $$$             'g',...
% $$$             'Filled');
    %line([-100,100],[100,-100],'Color','k')
    xlim([-150,150]);ylim([-150,150]);
    Lines([],0,'k');Lines(0,[],'k');    
    if p~=1,
        sax(end).XTick = [];
        sax(end).YTick = [];
    end
end




tind = [3,4,5,17,18,19,20,21,22,23];
tind = [6,7,26,27,30];
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


%% HEAD VELOCITY FORWARD ---------------------------------------------------------------------------
clear('xcomp','ycomp','zcomp','ccomp');
xcomp.label = 'hvf';            xcomp.data = [];    xcomp.edgs = [-25,-5,5,25,80]; 
ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5, 4);
ccomp.label = 'ego fwd';        ccomp.data = [];    ccomp.clim = [-40,40];
%for t = 4:10%:numel(tind)-1
%for t = 4:10%:numel(tind)-1    
for t = 1:4%:numel(tind)-1        
    dc = dca{t};
    mind = dc.stcm(:,1)==1                                      ...
           & (dc.stcm(:,3)==3|dc.stcm(:,5)==5)  ...
           & dc.ucnt>=2                                         ...
           & sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))<325;
    mind(mind==true) = randn([sum(mind),1])>0.25;
    xcomp.data = cat(1,xcomp.data,dc.hvfl(mind,1));    
    ycomp.data = cat(1,ycomp.data,dc.phz(mind));
    ccomp.data = cat(1,ccomp.data,dc.esax(mind,1)-25);
end
[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,ccomp);
zmean   = zcomp.mean;    zmean(zcomp.count<5)  = nan;
zmedian = zcomp.median;  zmean(zcomp.count<5)  = nan;
zstd    = zcomp.std;     zstd(zcomp.count<5)   = nan;
zcount  = zcomp.count;   zcount(zcomp.count<5) = nan;
set(figure(),'Units','centimeters','Position',[0,-3,8,28]);
subplot(511); imagesc(xcomp.ctrs,ycomp.ctrs,zmean'); axis('xy');
    cax = colorbar(); ylabel(cax,['Mean ',ccomp.label]); colormap('jet');  caxis([ccomp.clim]);
subplot(512); imagesc(xcomp.ctrs,ycomp.ctrs,zmedian'); axis('xy');
    cax = colorbar(); ylabel(cax,['Median ',ccomp.label]); colormap('jet');  caxis([ccomp.clim]);
subplot(513); imagesc(xcomp.ctrs,ycomp.ctrs,zstd'); axis('xy');
    cax = colorbar(); ylabel(cax,['Std ',ccomp.label]); ylabel(ycomp.label);
    colormap('jet');
subplot(514); 
imagesc(xcomp.ctrs,ycomp.ctrs,zcount'); axis('xy');
    cax = colorbar(); ylabel(cax,'Count'); colormap('jet'); 
    xlabel(xcomp.label);
subplot(515); 
    imagesc(xcomp.ctrs,ycomp.ctrs,bsxfun(@rdivide,zcount,sum(zcount,2,'omitnan'))'); axis('xy');    
    cax = colorbar(); ylabel(cax,'Probability'); colormap('jet'); 
    xlabel(xcomp.label);


figure,
pinds = {ycomp.data > 0.5  & ycomp.data < 2.25,...
         ycomp.data > 2.25  & ycomp.data < 4.25,...
         ycomp.data > 4.25  & ycomp.data < 6};
for p = 1:numel(pinds)
subplot2(3,3,2,p);
    hold('on');
    indL =    pinds{p} & xcomp.data > -5 & xcomp.data < 5;
    indC =    pinds{p} & xcomp.data >  5 & xcomp.data < 25;
    indR =    pinds{p} & xcomp.data < 40 & xcomp.data > 25;    
    indZ =    pinds{p} & xcomp.data < 100 & xcomp.data > 40;
    
    histogram(ccomp.data(indL),linspace(-400,400,50),'Normalization','probability');
    histogram(ccomp.data(indC),linspace(-400,400,50),'Normalization','probability');
    histogram(ccomp.data(indR),linspace(-400,400,50),'Normalization','probability');
    histogram(ccomp.data(indZ),linspace(-400,400,50),'Normalization','probability');
    
    Lines(mean(ccomp.data(indL)),[],'b');    
    Lines(mean(ccomp.data(indC)),[],'g');
    Lines(mean(ccomp.data(indR)),[],'r');
    Lines(mean(ccomp.data(indZ)),[],'m');    
    xlim([-300,300]);
    
    indA = indR|indC|indZ;%indL|
    [B,BINT,R,RINT,RSTATS] = regress(ccomp.data(indA),[ones([sum(indA),1]),xcomp.data(indA)]);
    
% RSTATS R-square statistic,  F statistic and p value, error variance
subplot2(3,3,1,p);
    P = polyfit(xcomp.data(indA),ccomp.data(indA),1);
    hold('on');
    hist2([ccomp.data(indA),xcomp.data(indA)],linspace(-300,300,30),linspace(-10,100,8),'xprob');
    line(polyval(P,[-10,100]),[-10,100],'Color','m');
    % $$$ figure,hist2([ccomp.data(indA),ycomp.data(indA)],linspace(-pi,pi,11),linspace(-1,1,11),'xprob');
    % $$$ figure,plot(ccomp.data(indA),ycomp.data(indA),'.');
% $$$     [h,p,ci,tstats] = ttest2(ccomp.data(indR),ccomp.data(indC))
% $$$     [h,p,ci,tstats] = ttest2(ccomp.data(indL),ccomp.data(indC))
% $$$     [h,p,ci,tstats] = ttest2(ccomp.data(indL),ccomp.data(indR))
subplot2(3,3,3,p);
    hold('on');
    cdfplot(ccomp.data(indL))
    cdfplot(ccomp.data(indC))
    cdfplot(ccomp.data(indR))
    cdfplot(ccomp.data(indZ))    
    xlim([-300,300]);
RSTATS    
end
% sessions with retrospective coding
% ec16-00267278 





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
xcomp.label = 'hba';            xcomp.data = [];    xcomp.edgs = [-5,5,25,80];
ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5,6 );
ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-100,100];
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
           & (dc.stcm(:,3)==3|dc.stcm(:,5)==5)  ...%           
           & abs(dc.hvfl(:,2))<50 ...
           & dc.ucnt>=3 & dc.ucnt<10 ...
           & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)<rdists(r)+10 ...
           & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)>rdists(r);
% $$$                        & sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))>200;

    mind(mind==true) = randn([sum(mind),1])>0.5;
    xcomp.data = cat(1, xcomp.data, dc.hvfl(mind,1));
    ycomp.data = cat(1, ycomp.data, dc.phz(mind));
    %ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+0*double(t<5)-12.5*double(t>=5));    
    %ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+25*double(t<4)-12.5*double(t>=4));
    ccomp.data = cat(1, ccomp.data, dc.esax(mind,1)-25);
    fcomp.data = cat(1, fcomp.data, dc.esax(mind,2));
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
    bind = mcomp.data<sectors(s+1) & mcomp.data>sectors(s) & ycomp.data>0.25& ycomp.data<2.25;
    %bind = mcomp.data<sectors(s+1) & mcomp.data>sectors(s) & ycomp.data>4& ycomp.data<6;    
    ind = bind & xcomp.data<5;
    medR(r,s) = median(ccomp.data(ind,1)./10);%-rmodel(r,s);
    skwR(r,s) = skew(ccomp.data(ind,1)./10);
    stdR(r,s) = std(ccomp.data(ind,1)./10);    
    ind = bind & xcomp.data<25& xcomp.data>5;
    medC(r,s) = median(ccomp.data(ind,1)./10);%-rmodel(r,s);
    skwC(r,s) = skew(ccomp.data(ind,1)./10);
    stdC(r,s) = std(ccomp.data(ind,1)./10);    
    ind = bind & xcomp.data>25;
    medL(r,s) = median(ccomp.data(ind,1)./10);%-rmodel(r,s);
    skwL(r,s) = skew(ccomp.data(ind,1)./10);
    stdL(r,s) = std(ccomp.data(ind,1)./10);    
    medD(r,s) = medR(r,s)-medL(r,s);
    skwD(r,s) = skwR(r,s)-skwL(r,s);
    stdD(r,s) = stdR(r,s)-stdL(r,s);
end
end

% $$$ figure,plot([medR(:,1);medR(:,7)])
% $$$ hold('on');
% $$$ plot([medR(:,end);medR(:,6)])
% $$$ 
% $$$ figure,histogram(medD(:))
% $$$  figure,histogram(medL(:)-medC(:))

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
subplot2(4,4,2,1);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-10,10]);
ylabel(cax,'cm');
title({'Leftward','Median'});
cmedD = stdL; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,1);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([0,20]);
ylabel(cax,'Std')
cmedD = skwL; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,1);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-2,2]);
ylabel(cax,'Skew')

cmedD = medC; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,2);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-10,10]);
ylabel(cax,'cm');
title({'Centered','Median'})
cmedD = stdC; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,2);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([0,20]);
ylabel(cax,'Std');
cmedD = skwC; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,2);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-2,2]);
title('Skew')

cmedD = medR; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,3);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-10,10])
ylabel(cax,'cm');
title({'Rightward','Median'});
cmedD = stdR; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,3,3);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([0,20]);
ylabel(cax,'Std');
cmedD = skwR; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,4,3);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-2,2]);
ylabel(cax,'Skew');
       
cmedD = medD; cmedD = cat(2,cmedD(:,end),cmedD);cmedD = cat(1,cmedD(end,:),cmedD);
subplot2(4,4,2,4);surf(A,B,cmedD,'edgecolor','none'),view([0,90]);colormap('jet');cax=colorbar();caxis([-10,10]);
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
