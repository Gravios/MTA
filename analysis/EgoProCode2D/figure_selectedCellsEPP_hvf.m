

cf(@(T) T.load('nq'), Trials);

binPhzs = linspace(0.5,2*pi-0.5,4);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
binHvfs = [-20,-5,5,25,80];
binHvfc = (binHvfs(1:end-1)+binHvfs(2:end))./2;

hvfBinEdg = binHvfs;
hvfBinCtr = binHvfc;

%pfo = cell([1,30]);
%pfo{18} = pfv{1}; % hba <- req20200630.m
pfo = pfv; % hba <- req20200630.m
% $$$ pfo = pfl; % hvl <- req20200630.m
% $$$ pfo = pfv; % hvf <- req20200630.m
% $$$ pfo = pfr; % roll<- req20200630.m

% change incase first is empty
nBins = numel(hvfBinCtr);
nPhz = numel(binPhzc);


for  t = 1:numel(Trials),
    %for  t = 1:numel(Trials),
    Trial = Trials{t};    
    Trial.load('nq');
    figDir = fullfile('/storage/share/Projects/EgoProCode2D/eppHvf',Trial.filebase);
    create_directory(figDir);
    
    [accg,tbins] = autoccg(Trial);



% HVF setup
hvfl = fet_href_HXY(Trials{t},sampleRate,false,'trb',4);
hvfBinInd = discretize(hvfl(:,1), hvfBinEdg);


varBinEdg = hvfBinEdg;
varBinCtr = hvfBinCtr;
varBinInd = hvfBinInd;

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


ind = [Trials{t}.stc{'w+p&t',sampleRate}];
%ind = [Trials{t}.stc{'theta-groom-sit-rear',sampleRate}];

ind.cast('TimeSeries');
ind.resample(xyz{t});



if isempty(units{t}),
    continue,
end
for u = 1:numel(units{t}),

[mxr,mxp] = pft{t}.maxRate(units{t}(u));        
pfstrj = MTADfet.encapsulate( ...
        Trials{t},...
        multiprod(bsxfun(@minus,...
                         mxp,...
                         sq(xyz{t}(:,'hcom',[1,2]))),...
                  hvec,2,[2,3]),...
        sampleRate,...
        'pfstrj','ptrj','t');

    
% SET figure layout
[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','landscape',[],2.5,2.5,0.5,0.5);
%[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','portrait',[],2.5,2.5,0.2,0.2);
[yind, yOffSet, xind, xOffSet] = deal(1,                         ... yind
                                      0,                         ... yOffSet
                                      1,                         ... xind
                                      0);                          % xOffSet

xBlockOffset = 0;
yBlockOffset = 0;

unit = units{t}(u);


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
[yind, yOffSet, xind, xOffSet] = deal(   1,       1,    3,       fig.subplot.width/2);        
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
[yind, yOffSet, xind, xOffSet] = deal(   1,       1,    4,       fig.subplot.width/2);        
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
[yind, yOffSet, xind, xOffSet] = deal(   1,      1,    5,       fig.subplot.width/2);        
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
[yind, yOffSet, xind, xOffSet] = deal(   1,       1,    6,       fig.subplot.width/2);        
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
[yind, yOffSet, xind, xOffSet] = deal(   p+2,       1,    1,       fig.subplot.width/2);        
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
    [yind, yOffSet, xind, xOffSet] = deal(   2,       0,    a+2,       0);        
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
% $$$ 
% $$$ if displayModelFlag,
% $$$     xBlockOffset = 2;
% $$$     yBlockOffset = 0;
% $$$     for a = 1:nBins
% $$$         [yind, yOffSet, xind, xOffSet] = deal(   nBins+6,       1,    a,       0);                
% $$$         sax(end+1) = axes('Units','centimeters',                                        ...
% $$$                           'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
% $$$                             fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
% $$$                             fig.subplot.width,                                ...
% $$$                             fig.subplot.height],                              ...
% $$$                           'FontSize', 8,                                                ...
% $$$                           'LineWidth',1);
% $$$         hold(sax(end),'on');
% $$$         subject = struct(rat);
% $$$         subject = update_subject_patch(subject,'head', a, true,varBinEdg,varBinCtr);
% $$$         subject = update_subject_patch(subject,'body',[],false,varBinEdg,varBinCtr);
% $$$         patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$         patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$         patch(subject.head.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
% $$$         line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
% $$$         line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
% $$$         daspect(sax(end),[1,1,1]);
% $$$         xlim([-100,100]);
% $$$         ylim([-100,100]);    
% $$$         sax(end).XTick = [];
% $$$         sax(end).YTick = [];
% $$$     end
% $$$ end% if displayModelFlag
% $$$ 
xBlockOffset = 2;
yBlockOffset = 0;


for p = 1:nPhz,
    for a = 1:nBins,
        if isempty(pfo{t}{p,a}), continue; end;
        rmap = plot(pfo{t}{p,a},unit,1,[],[],false);        
        [yind, yOffSet, xind, xOffSet] = deal(   nPhz+3-p,       0,    a,       0);                
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
% $$$ 
% $$$         
% $$$         % ADD subject
% $$$         if displayModelFlag,
% $$$             subject = struct(rat);
% $$$             subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
% $$$             subject = update_subject_patch(subject,'body', nBins+1-a,  true,hbaBinEdg,hbaBinCtr);
% $$$             patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$             patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$             patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
% $$$             line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
% $$$             line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
% $$$         else
            if p == 1,
                xlabel([num2str(varBinCtr(a)),' mm/s']);
            end
% $$$         end;%if displayModelFlag
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
rmap(rmap<2) = 0;
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

hclr = 'ckm';
figure,
hold('on')
for a = 1:3
plot(egoMeanRmapPos(:,1,a+1,1),egoMeanRmapPos(:,1,a+1,2),['.',hclr(a)]);
end

figure
violin(sq(egoMeanRmapPos(:,1,2:4,1)))

figure
violin(sq(egoMeanRmapPos(:,1,[3,4],1))-sq(egoMeanRmapPos(:,1,[2,2],1)))

[H,P,CI,STATS] = ttest2(sq(egoMeanRmapPos(:,1,[2],1))-sq(egoMeanRmapPos(:,1,[3],1)),...
       sq(egoMeanRmapPos(:,1,[4],1))-sq(egoMeanRmapPos(:,1,[3],1)));

[H,P,CI,STATS] = ttest(sq(egoMeanRmapPos(:,1,[2],1)));
[H,P,CI,STATS] = ttest(sq(egoMeanRmapPos(:,1,[3],1)));
[H,P,CI,STATS] = ttest(sq(egoMeanRmapPos(:,1,[4],1)));

[H,P,CI,STATS] = ttest(sq(egoMeanRmapPos(:,1,[3],1))-sq(egoMeanRmapPos(:,1,[2],1)));

figure,histogram(sq(egoMeanRmapPos(:,1,[4],1))-sq(egoMeanRmapPos(:,1,[2],1)),linspace(-200,150,30))

unitsEgoHvfMap =[];
for c = 1:numel(unitsEgoHvf)
    unitsEgoHvfMap = [unitsEgoHvfMap;ones([numel(unitsEgoHvf{c}),1]).*c];
end

uids = ismember(unitsEgoHvfMap,[3,4,5,17:25,29]);
uids = ismember(unitsEgoHvfMap,[3,4,5,17:25]);
uids = ismember(unitsEgoHvfMap,[3,4,5]);
uids = ismember(unitsEgoHvfMap,[17:25]);
uids = ismember(unitsEgoHvfMap,[29]);


[H,P,CI,STATS] = ttest(sq(egoMeanRmapPos(uids,1,[2],1)));
[H,P,CI,STATS] = ttest(sq(egoMeanRmapPos(uids,1,[3],1)));
[H,P,CI,STATS] = ttest(sq(egoMeanRmapPos(uids,1,[4],1)));


[H,P,CI,STATS] = ttest2(sq(egoMeanRmapPos(uids,3,[2],1)),...
                        sq(egoMeanRmapPos(uids,3,[3],1)))


[H,P,CI,STATS] = ttest2(sq(egoMeanRmapPos(uids,1,[2],1)),...
                        sq(egoMeanRmapPos(uids,1,[4],1)))



figure,
for p = 1:3
    subplot(3,1,4-p);
    hold('on');
    histogram(sq(egoMeanRmapPos(uids,p,[2],1)-25),linspace(-150,250,30));
    histogram(sq(egoMeanRmapPos(uids,p,[3],1)-25),linspace(-150,250,30));
    histogram(sq(egoMeanRmapPos(uids,p,[4],1)-25),linspace(-150,250,30));
end

figure
for p = 1:3
    subplot(3,1,p);
    violin(sq(egoMeanRmapPos(uids,4-p,2:4,1))-25);
end







