% req20200630_figure_singleUnitEPP
%
%    for each unit 
%        - plot theta spatial ratemap
%        - plot theta-pfc posture ratemap
%        - plot grid
%            - theta phase
%            - head-body-angle 
%        
%


MjgER2016_defargs();

if ~exist('bfrm','var')
    bfrm = cf(@(t,u)   compute_bhv_ratemaps(t,u),                 Trials, units);
end

lims = {[-250,250],[-250,250]};

% FIGURE - Examples of individual units

hbaBinEdg = -1.5:0.6:1.5;
hbaBinCtr = mean([hbaBinEdg(1:end-1);hbaBinEdg(2:end)]);

ascnPhz = 4;
descPhz = 2;

%[ 0-120, 120-240, 240-360 ]

orientationLabels = 'LLCRR';

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
%    BLOCK 5x5 HBAxPHZ


mazeMaskFlag = true;                                      
pfsPlotMode  = 1;
reportMethod ='text';
pfsColorMap  = @jet;

t = 20;
u = 2;

for  t = 1:numel(Trials),
Trial = Trials{t};    
Trial.load('nq');
figDir = fullfile('/storage/share/Projects/BehaviorPlaceCode/egocentricPP/singleUnitEPP',Trial.filebase);
create_directory(figDir);
[accg,tbins] = autoccg(Trial);
if isempty(units{t}),
    continue,
end
for u = 1:numel(units{t}),

% SET figure layout
[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','landscape',[],2.5,2.5,0.5,0.5);
[yind, yOffSet, xind, xOffSet] = deal(1,                         ... yind
                                      0,                         ... yOffSet
                                      1,                         ... xind
                                      0);                          % xOffSet

xBlockOffset = 0;
yBlockOffset = 0;

unit = units{t}(u);

rateMinMax = max(cell2mat(cf(@(p) prctile(p.data.rateMap(:,p.data.clu==unit,1),99.5), pfs{t})));


% INFORMATION -------------------------------------------------------------------------------------
[yind, yOffSet, xind, xOffSet] = deal(   1,       1,    1,       0);        
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
                          fig.subplot.height],                                                  ...
              'FontSize', 8);
% -------------------------------------------------------------------------------------------------

% THETA PLACE FIELD -------------------------------------------------------------------------------
[yind, yOffSet, xind, xOffSet] = deal(   3,       1,    1,       fig.subplot.width/2);        
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
[yind, yOffSet, xind, xOffSet] = deal(   4,       1,    1,       fig.subplot.width/2);        
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
%sax(end).XTick = [];
%sax(end).YTick = [];
xlabel('head pitch');
ylabel('body pitch');
% -------------------------------------------------------------------------------------------------

% THETA accg --------------------------------------------------------------------------------------
[yind, yOffSet, xind, xOffSet] = deal(   2,       1,    1,       fig.subplot.width/2);        
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



% MAINBLOCK - egofields----------------------------------------------------------------------------
al = 1:pfs{t}{1}.adata.binSizes(end);
xBlockOffset = 2;
yBlockOffset = 0;
for a = 1:numel(al)
    [yind, yOffSet, xind, xOffSet] = deal(   numel(al)+1,       1,    a,       0);                
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
                                  fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
                                  fig.subplot.width,                                ...
                                  fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');
    subject = struct(rat);
    subject = update_subject_patch(subject,'head', a, true,hbaBinEdg,hbaBinCtr);
    subject = update_subject_patch(subject,'body',[],false,hbaBinEdg,hbaBinCtr);
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

xBlockOffset = 2;
yBlockOffset = 0;


numPhz = numel(pfs{t});
numAng = numel(al);
for p = 1:numPhz,
    rmap = plot(pfs{t}{p},unit);
    for a = 1:numAng,
        [yind, yOffSet, xind, xOffSet] = deal(   numPhz+1-p,       1,    a,       0);                
        sax(end+1) = axes('Units','centimeters',                                        ...
                          'Position',[fig.page.xpos(xind+xBlockOffset)+xOffSet,         ...
                                      fig.page.ypos(yind+yBlockOffset)+yOffSet,         ...
                                      fig.subplot.width,                                ...
                                      fig.subplot.height],                              ...
                          'FontSize', 8,                                                ...
                          'LineWidth',1);
        hold(sax(end),'on');
        
        pcolor(pfs{t}{p}.adata.bins{1},...
                       pfs{t}{p}.adata.bins{2},...
                       rmap(:,:,al(a)));
        
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
            if p==round(numPhz/2),
                ylabel({'Theta Phase',round(circ_rad2ang(binPhzc(p))),' '});
            end
        end
        if p == numPhz && a == round(numAng/2),
            title([pfs{t}{end}.tag,'_unit_',num2str(unit)]);
        end
        
        set(sax(end),'XTick',[]);
        set(sax(end),'YTick',[]);        
        
        % ADD subject
        if p %== 4,
            subject = struct(rat);
            subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
            subject = update_subject_patch(subject,'body', numAng+1-a,  true,hbaBinEdg,hbaBinCtr);
            patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
            patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
            patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
            line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
            line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
        end
    end
end

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

line(sax(end-numPhz*numAng+1).Position(1).*[1,1]-0.25,                    ...
     [sax(end-numPhz*numAng+1).Position(2),                               ...
      sax(end-numPhz+1).Position(2)+sax(end-numPhz).Position(4)],         ...
     'LineWidth',2,                                                       ...
     'Color','k');

print(hfig,                                                               ...% figure handle
      '-dpng',                                                            ...% image format
      fullfile(figDir,[pfs{t}{end}.tag,'_unit_',num2str(unit),'.png']));     % image path

% -------------------------------------------------------------------------------------------------
        
end% for u : units
end% for t : trials

