

cf(@(T) T.load('nq'), Trials);

binPhzs = linspace(0.5,2*pi-0.5,4);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
binHvls = [-60,-5,5,60];
binHvlc = (binHvls(1:end-1)+binHvls(2:end))./2;

hvlBinEdg = binHvls;
hvlBinCtr = binHvlc;
lims = {[-300,300],[-300,300]};
%pfo = cell([1,30]);
%pfo{18} = pfv{1}; % hba <- req20200630.m
pfo = pfl; % hba <- req20200630.m
% $$$ pfo = pfl; % hvl <- req20200630.m
% $$$ pfo = pfv; % hvf <- req20200630.m
% $$$ pfo = pfr; % roll<- req20200630.m
mazeMaskFlag = true;                                      
pfsPlotMode  = 1;
reportMethod ='text';
pfsColorMap  = @jet;
displayModelFlag = true;
% change incase first is empty
nBins = numel(hvlBinCtr);
nPhz = numel(binPhzc);
if ~exist('bfrm','var'),bfrm = cf(@(t,u) compute_bhv_ratemaps(t,u),Trials,units);end;

for  t = 20:numel(Trials),
    %for  t = 1:numel(Trials),
    Trial = Trials{t};    
    Trial.load('nq');
    figDir = fullfile('/storage/share/Projects/EgoProCode2D/eppHvl',Trial.filebase);
    create_directory(figDir);
    
    [accg,tbins] = autoccg(Trial);



% HVF setup
hvfl = fet_href_HXY(Trials{t},sampleRate,false,'trb',4);
hvlBinInd = discretize(hvfl(:,2), hvlBinEdg);


varBinEdg = hvlBinEdg;
varBinCtr = hvlBinCtr;
varBinInd = hvlBinInd;


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



