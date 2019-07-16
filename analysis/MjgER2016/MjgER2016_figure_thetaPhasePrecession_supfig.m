% REAR ,               ,  7 ,  3.2
% REAR , high @ peak   , 15 ,  6.9
% REAR , high @ trough ,  0 ,  0.0

% HIGH ,               ,  5 ,  2.3
% HIGH , rear @ trough ,  3 ,  1.4
% HIGH , rear @ peak   ,  3 ,  1.4
% HIGH , rear rate mod ,  0 ,  0.0
% HIGH , low @ peak    , 31 , 14.4
% HIGH , Low @ trough  ,  1 ,  0.5
% HIGH , Low rate mod  , 10 ,  4.6

% LOW  ,               , 15 ,  6.9
% LOW  , rear @ peak   ,  0 ,  0.0
% LOW  , rear @ trough ,  0 ,  0.0
% LOW  , rear rate mod ,  0 ,  0.0
% LOW  , high @ peak   , 15 ,  6.9
% LOW  , high @ trough , 16 ,  7.4
% LOW  , high rate mod , 28 , 13.0

% REAR & HIGH,         , 16 ,  7.4         
% REAR & LOW ,         ,  2 ,  0.9
% LOW & HIGH ,         , 35 , 16.2    
% REAR & HIGH & LOW,   , 16 ,  7.4
                       
% TOTAL ,              , 216,
% GHZ VS PHZ ratemap examples ----------------------------------------------------------------------

statesIndsGPE = [2,3,4];

mrMaxPhzRate = [0,35];


tppsUnits = {};
tppsNames = {};
tppsTitles  = {};
tppsStsInds = {};
% REAR 
tppsNames{end+1}  = 'r';
tppsTitles{end+1} = 'rear';
tppsStsInds{end+1} = [1,2];
tppsUnits{end+1} = { [ 3], [62,130];... 62hbadedge
                     [ 4], [100];... 100hbadedge
                     [ 5], [55];...
                     [20], [17];...
                     [21], [40];...
                     [23], [58];...
                   };

% REAR , high @ peak
tppsNames{end+1}  = 'r-hp';
tppsTitles{end+1} = 'rear , high @ peak';
tppsStsInds{end+1} = [1,2];                
tppsUnits{end+1} = { [ 3], [59,145];... tp
                     [ 4], [26];... tp 26rear,bedg
                     [ 5], [4,73];... tt
                     [10], [5];...
                     [18], [39,46];...
                     [19],[117,136];... rear and high walk have bad estimates need central values
                     [20],[133];...
                     [21], [13,14];...
                     [22], [24];... poor edged estimates in rear
                     [23], [18];...
                   };

% REAR & HIGH
tppsNames{end+1}  = 'rh';
tppsTitles{end+1} = 'rear & high';
tppsStsInds{end+1} = [1,2];
tppsUnits{end+1} = { [ 3], [33,38];...
                     [ 4], [168,184,197];... 168nostruct 184rearbadedge
                     [ 5], [15,71,159];...
                     [ 9], [20,23];...
                     [11], [16];...
                     [18], [21,50];...                
                     [20], [83];...
                     [21], [27,52];...
                   };

% REAR & LOW
tppsNames{end+1}  = 'rl';
tppsTitles{end+1} = 'rear & low';
tppsStsInds{end+1} = [1,3];
tppsUnits{end+1} = { [5], [129];...
                     [20], [85];...
                   };

% LOW & HIGH
tppsNames{end+1}  = 'lh';
tppsTitles{end+1} = 'low & high';
tppsStsInds{end+1} = [2,3];
cstates = [2,3];                
tppsUnits{end+1} = { [3], [114,139,165];...
                [5], [72,99];...
               [17], [16,19,69];...
               [18], [17,23,35,42,60,84];...    
               [19], [ 3,10,27,28,50,63,160];...
               [20], [25,28,44,51,60,118,140];...
               [21], [ 6,32,77];...
               [22], [13,30,48];...
               [23], [50];...
              };

% LOW & HIGH with recession
tppsNames{end+1}  = 'lh-recession';
tppsTitles{end+1} = 'low & high with recession';
tppsStsInds{end+1} = [2,3];
tppsUnits{end+1} = {[20], [44];...
                   };

% REAR & HIGH & LOW
tppsNames{end+1}  = 'rhl';
tppsTitles{end+1} = 'rear & high & low';
tppsStsInds{end+1} = [1,2,3];
tppsUnits{end+1} = { [ 3], [50,111];...
                     [ 4], [69,141];...
                     [ 5], [16,163];...
                     [18], [44,52,91];...
                     [19], [13,97,141];...
                     [20], [110];...
                     [21], [43,64];...
                     [22], [38];...
                   };

% REAR , HIGH & LOW flat
tppsNames{end+1}  = 'r-hl';
tppsTitles{end+1} = 'rear , high & low flat';
tppsStsInds{end+1} = [1,2,3];
tppsUnits{end+1} = {[17], [14,107];...
                   };
               
               
% HIGH ----------------------------------------
tppsNames{end+1}  = 'h';
tppsTitles{end+1} = 'high';
tppsStsInds{end+1} = [2,3];
tppsUnits{end+1} = { [ 4], [90];...
                     [ 9], [24];...
                     [18], [25];... low has recession
                     [19], [15];...
                     [22], [23];...
                   };
               
% HIGH , rear @ trough 
tppsNames{end+1}  = 'h-rt';
tppsTitles{end+1} = 'high';
tppsStsInds{end+1} = [1,2];
tppsUnits{end+1} = { [3], [159];...
                     [19], [150];... low @ peak ??
                     [23], [36];...
                   };

% HIGH , rear @ peak
tppsNames{end+1}  = 'h-rp';
tppsTitles{end+1} = 'high';
tppsStsInds{end+1} = [1,2];
tppsUnits{end+1} = { [4], [126];...
               [20], [20];...
               [21], [59];...
               };
% HIGH , low @ peak
tppsNames{end+1}  = 'h-lp';
tppsTitles{end+1} = 'high , low @ peak';
tppsStsInds{end+1} = [2,3];
tppsUnits{end+1} = { [3],  [23,80,158,173];...
                     [5],  [31,67,81];...
                     [9],  [23];...
                     [10], [28];...
                     [17], [58,61];... 
                     [18], [19,24,57,74,75];...
                     [19], [21,45,95,140,150];...                
                     [20], [20,72,79];...
                     [21], [22,61];...                
                     [22], [8,11,37,61,65];... bexp 61
                   };

% HIGH , Low rate mod
tppsNames{end+1}  = 'h-l-rm';
tppsTitles{end+1} = 'high , low rate mod';
tppsStsInds{end+1} = [2,3];
tppsUnits{end+1} = { [3],  [135];...
                     [11], [14,21];...
                     [18], [13];...
                     [19], [172];...
                     [20], [111,134];...
                     [21], [2,59];...
                     [22], [43]
                   };

% HIGH , Low @ trough 
tppsNames{end+1}  = 'h-lt';
tppsTitles{end+1} = 'high , low @ trough';
tppsStsInds{end+1} = [2,3];
tppsUnits{end+1} = { [3], [80]};

% HIGH , LOW RECESSION
tppsNames{end+1}  = 'h-l-recession';
tppsTitles{end+1} = 'high , low recession';
tppsStsInds{end+1} = [2,3];
tppsUnits{end+1} = {[17], [58];...
                   };%
                

% LOW ----------------------------------------
tppsNames{end+1}  = 'l';
tppsTitles{end+1} = 'low';
tppsStsInds{end+1} = [2,3];
tppsUnits{end+1} = { [ 3], [61,131];... 61badedge
                     [ 4], [94];...
                     [ 5], [90];...
                     [19], [19,108,114,115,158,167];...
                     [20], [89,103,104];...
                     [23], [56,66];...
                   };
               
% LOW , high @ peak
tppsNames{end+1}  = 'l-hp';
tppsTitles{end+1} = 'low , high @ peak';
tppsStsInds{end+1} = [2,3];
tppsUnits{end+1} = { [ 4], [33,38];...
                     [ 5], [35,56,93];...
                     [17], [ 8,32,87];...
                     [18], [11,18];...
                     [19], [32,66];...
                     [20], [80,141];...
                     [22], [15];...
                   };%


% LOW , high @ trough
tppsNames{end+1}  = 'l-ht';
tppsTitles{end+1} = 'low , high @ trough';
tppsStsInds{end+1} = [2,3];
tppsUnits{end+1} = { [ 5], [40];... 
                     [17], [20];...                 
                     [18], [39,46,49];... 
                     [19], [94,143,152,153];...
                     [20], [35,139];... 
                     [22], [50,58];... 
                     [23], [48,66,69];...
                   };

% WEAK LOW , high @ trough
tppsNames{end+1}  = 'wl-ht';
tppsTitles{end+1} = 'weak low , high @ trough';
tppsStsInds{end+1} = [2,3];
tppsUnits{end+1} = {[21], [6]};% low     tt
                

% LOW , high rate mod
tppsNames{end+1}  = 'l-h-rm';
tppsTitles{end+1} = 'low , high rate mod';
tppsStsInds{end+1} = [2,3];
tppsUnits{end+1} = { [ 3], [151];...
                     [ 5], [74,93,112];...
                     [12], [8];... 
                     [17], [30,44,98];...                
                     [18], [20];...
                     [19], [53,96];...
                     [20], [21,31,59,81,119,129,151];...
                     [21], [63,65,97];...
                     [22], [56];...
                     [23], [9,29,51,54,70,72];... 29r@tp
              };
                 

% LOW , RECESSION
tppsNames{end+1}  = 'l-recession';
tppsTitles{end+1} = 'low , recession';
tppsStsInds{end+1} = [2,3];
tppsUnits{end+1} = { [5], [86];...
                     [18], [60];...
                     [20], [141];...
                     [21], [32];...
                     [23], [48];...
                   };


% MESS 
% $$$ tppsUnits{end+1} = {[18], [16];...
% $$$               };%
                
                 

% drop index: bad occ in high or low bad
%             bad rate est in rear
%             bad
%             messy field
%             bad edge est
%             low rate
dind = [85,52,91,82,85,89,110,152,205,246,248,308,367];
rind = [44,46,49,59,62,106,112,157,158,186,187,188,191,200,201,202,205,217,220,230,243,253,256,...
        258,259,262,272,275,286,287,317,330,339,348,357,364,371,374,377,398,402];
bind = [25,34,40,57,64,65,69,77,87,117,161,199,204,208,209,213,268,285,308,350,361,376,377,381,...
        383,388];
mind = [18,37,48,121,153,156,196,210,301,305,324,361,375];
eind = [191,196,202,372,406];
lind = [27,72,76,305];
               





% FIGURE 
% ∀ u ∈ tppsUnits{end+1}
%     PLOT place field
%     PLOT behavior field
%     PLOT rates asfo theta phase and drz
%     PLOT drz corrected rates asfo theta phase
% PLOT overlay subset onto population for rate shift vs phz shift

% SET
% figure opts

for g = 1:numel(tppsTitles),

[hfig,fig,fax,sax] = set_figure_layout(figure(666002),'A4','portrait',[],1.1,1.1,0);

expUnitsCount = sum(cellfun(@numel,tppsUnits{g}(:,2)));
yind = 0;
ucnt = 1;
xshift = 0;
for tind = 1:size(tppsUnits{g},1)
% SET indexing -------------------------------------------------------------------------------------
    t = tppsUnits{g}{tind,1}(1);
    for uid = 1:numel(tppsUnits{g}{tind,2});
        u = tppsUnits{g}{tind,2}(uid);
        uind = find(ismember(cluSessionMapSubset,[t,u],'rows'));
        if ucnt==21        
            yind = 0;
        end
        if ucnt>20
            xshift = 8;
        else
            xshift = 0;
        end
        xind = 1+xshift;
        yind = yind+1;                                    
        
% PLOT place field ---------------------------------------------------------------------------------
        sax(end+1) = axes('Units','centimeters',                                             ...
                          'Position',[fig.page.xpos(xind),                                   ...
                                      fig.page.ypos(yind)-fig.subplot.height/2,              ...
                                      fig.subplot.width,                                     ...
                                      fig.subplot.height],                                   ...
                          'FontSize', 8,                                                     ...
                          'LineWidth',1);
        hold(sax(end),'on');
        plot(pfts{t},u,'mean',false,[],true,0.5,false,interpParPfs,@bone,[],nanColor);
        text(-480,-380,num2str(cond_round(pfts{t}.maxRate(u,true,'interpPar',interpParPfs))),...
             'FontSize',8,'Color',[1,1,1]);
        sax(end).XTickLabels = {};
        sax(end).YTickLabels = {};
        ylabel(sprintf('%i-%i',t,u));
        % plot state placefield centers 
        pfssCntrMrate = zeros([numel(statesPfssInd),1]);
        pfssCntrPos = zeros([numel(statesPfssInd),2]);
        for s = 1:numel(statesPfssInd),
            [pfssCntrMrate(s),pfssCntrPos(s,:)] = maxRate(pfss{t}{statesPfssInd(s)},...
                                                          u,true,'interpPar',interpParPfs);
            if s~=1 & any(ismember(round(pfssCntrPos(1:s-1,:)),round(pfssCntrPos(s,:)),'rows')),
                shift = 10;
            else
                shift = 0;
            end
            
            if pfssCntrMrate(s)>2,
                hax = scatter(pfssCntrPos(s,1)+shift,pfssCntrPos(s,2),15,sclr(s),'filled');
            end
            
        end

        axis(sax(end),'tight');
        axes(fax);
        rectangle('Position',sax(end).Position,'LineWidth',1);

        
        
% PLOT behavior field ------------------------------------------------------------------------------
        xind = xind+1;
        sax(end+1) = axes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind),                      ...
                                      fig.page.ypos(yind)-fig.subplot.height/2, ...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height],                      ...
                          'FontSize', 8,                                        ...
                          'LineWidth',1);
        plot(pfbs{t},u,'mean',[],[],false,0.5,false,...
             interpParDfs,@bone,reshape(validDims,pfbs{t}.adata.binSizes'),nanColor);    
        text(-1.7,-0.45,num2str(cond_round(pfbs{t}.maxRate(u,false,'mask',validDims))),...
             'FontSize',8,'Color',[1,1,1]);
        %text(-1.7,-0.45,num2str(cond_round(maxPfsRate)),'FontSize',10,'Color',[1,1,1]);
        sax(end).XTickLabels = {};
        sax(end).YTickLabels = {};
        xlim([-1.8,0.75]);
        ylim([-0.75,1.8]);
        axes(fax);
        rectangle('Position',sax(end).Position,'LineWidth',1);

        
        
% PLOT rates asfo theta phase and drz --------------------------------------------------------------
        mrateHZTPD = max(cell2mat(cf(@(p) p{t}.maxRate(u,false), pftHZTPD(statesIndsGPE(:)))));
        for s = 1:numel(statesIndsGPE),
            % SET horizontal offset    
            xind = 2+s+xshift;    
            % CREATE subplot axes
            sax(end+1) = axes('Units','centimeters',                                             ...
                              'Position',[fig.page.xpos(xind)+fig.subplot.width/2,               ...
                                          fig.page.ypos(yind)-fig.subplot.height/2,              ...
                                          fig.subplot.width,                                     ...
                                          fig.subplot.height/2],                                 ...
                              'FontSize', 8,                                                     ...
                              'LineWidth',1);
            
            plot(pftHZTPD{statesIndsGPE(s)}{t},u,'mean','',[0,mrateHZTPD],false,[],[],[],@copper,[],nanColor);
            text(-.9,-1.8,num2str(cond_round(mrateHZTPD)),...
                 'FontSize',8,'Color',[1,1,1]);
            [hl,hp] = boundedline(mprPos(:,uind,s)',...
                        pftHZTPD{1}{1}.adata.bins{2}',...
                        0,...
                        ...% mprPosBSStd(:,uind,s)'.*1.96./10,...
                        ['-',sclr(s)],...
                        'alpha',...                        
                        'transparency',0.3,...
                        'orientation','horiz');
            hl.LineWidth = 1;
            hp.EdgeColor = 'none';
            Lines(0,[],'w');
            
            sax(end+1) = axes('Units','centimeters',                                             ...
                              'Position',[fig.page.xpos(xind)+fig.subplot.width/2,               ...
                                          fig.page.ypos(yind),                                   ...
                                          fig.subplot.width,                                     ...
                                          fig.subplot.height/2],                                 ...
                              'FontSize', 8,                                                     ...
                              'LineWidth',1);
            plot(pftHZTPD{statesIndsGPE(s)}{t},u,'mean','',[0,mrateHZTPD],false,[],[],[],@copper,[],nanColor);
            [hl,hp] = boundedline(mprPos(:,uind,s)',...
                                  pftHZTPD{1}{1}.adata.bins{2}',...
                                  0,...
                                    ... %mprPosBSStd(:,uind,s)'.*1.96./10,...
                                 ['-',sclr(s)],...
                                 'alpha',...
                        'transparency',0.3,...
                        'orientation','horiz');
            hl.LineWidth = 1;            
            hp.EdgeColor = 'none';
            Lines(0,[],'w');            
            if ucnt == 1,
                title(stateLabels(statesIndsGPE(s)));            
                sax(end  ).YTickLabels = {};
                sax(end-1).YTickLabels = {};                
                sax(end-1).XTickLabels = {};
                sax(end-1).XTickLabels = {};   
            elseif ucnt == expUnitsCount,
                sax(end  ).YTickLabels = {};
                sax(end-1).YTickLabels = {};
                if s~=1,
                    sax(end-1).XTickLabels = {};
                    sax(end-1).XTickLabels = {};   
                end
            else
                sax(end  ).YTickLabels = {};
                sax(end-1).YTickLabels = {};                
                sax(end  ).XTickLabels = {};
                sax(end-1).XTickLabels = {};   
            end
            
            sax(end).Color = 'none';
            sax(end-1).Color = 'none';            

            axes(fax);
            rectangle('Position',[sax(end-1).Position.*[1,1,1,2]],'LineWidth',1);
            
        end%for s
        
% PLOT drz corrected rates asfo theta phase --------------------------------------------------------
        sax(end+1) = axes('Units','centimeters',                                             ...
                          'Position',[fig.page.xpos(xind+1)+fig.subplot.width/2,             ...
                                      fig.page.ypos(yind)-fig.subplot.height/2,              ...
                                      fig.subplot.width,                                     ...
                                      fig.subplot.height],                                   ...
                          'FontSize', 8,                                                     ...
                          'LineWidth',1);
        hold('on');
        
        %for s = 1:numel(statesIndsGPE),
        for s = tppsStsInds{g},
            [hl,hp] = boundedline(repmat(maxPhzRate(:,uind,s),[2,1])',...
                                  pbins',...
                                  repmat(maxPhzRateStd(:,uind,s),[2,1])',...
                                  ['-',sclr(s)],...
                                  'alpha',...
                                  'transparency',0.3,...
                                  'orientation','horiz');
            hl.LineWidth = 1;
            hp.EdgeColor = 'none';
            [hl,hp] = boundedline(repmat(maxPhzRate(:,uind,s),[2,1])',...
                                  pbins',...
                                  repmat(maxPhzRateStd(:,uind,s),[2,1])'.*1.96./10,...
                                  ['-',sclr(s)],...
                                  'transparency',1,...
                                  'orientation','horiz');
            hl.LineWidth = 1;
            hp.EdgeColor = 'none';
        end
        %plot(mrMaxPhzRate(2)*0.8,360,markers(ucnt),'MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10);
        sax(end).YAxisLocation = 'right';
        if ucnt == expUnitsCount,
            xlabel('Firing Rate');
            ylabel('Theta Phase');
        else
            sax(end).XTickLabels = {};
        end
        
        sax(end).YTick = [0,180,360,540];
        sax(end).YTickLabels = [0,180,360,540];
        ylim(pbins([1,end]));
        %xlim(mrMaxPhzRate);
        axis('tight');
        xlim([0,max(xlim)+2]);

        
        ucnt = 1 + ucnt;
    end%for uid
end%for tind

% $$$ axes(fax);
% $$$ line(fax,fig.page.xpos(1)+[0,2/2],fig.page.ypos(yind).*[1,1]-0.15-fig.subplot.height/2,'Color','k','LineWidth',2);
% $$$ text(fax,fig.page.xpos(1),fig.page.ypos(yind)-0.35-fig.subplot.height/2,'50 cm','FontSize',8,'Color',[0,0,0]);
% $$$ line(fax,fig.page.xpos(2)+[0,2/2],fig.page.ypos(yind).*[1,1]-0.15-fig.subplot.height/2,'Color','k','LineWidth',2);
% $$$ text(fax,fig.page.xpos(2),fig.page.ypos(yind)-0.35-fig.subplot.height/2,'1 rad','FontSize',8,'Color',[0,0,0]);

if ucnt>15&ucnt<21
    xind = 9;
    yind = 4;
else
    xind = xshift+1;    
    yind = yind+4;
end



% PLOT overlay subset onto population for rate shift vs phz shift  ---------------------------------
sax(end+1) = axes('Units','centimeters',                                            ...
                  'Position',[fig.page.xpos(xind)+fig.subplot.width,                ...
                              fig.page.ypos(yind)-fig.subplot.height,               ...
                              fig.subplot.width*4,                                  ...
                              fig.subplot.height.*4],                               ...
                  'FontSize', 8,                                                    ...
                  'LineWidth',1);
% POPULATION 
%smarkers = '^sv';
smarkers = 'ooo';
for s = 1:3
    ind = dsi(:,1)==s & validPosOccRlx & sesInds(unitSubset) & mprMaxRateSrt(:,2) > 4 & validUnits;
    scatter(mprMeanRateSrt(ind,1)'-mprMeanRateSrt(ind,2)',            ...
            -circ_dist(mprMeanPhzAngSrt(ind,1)',                       ...
                      mprMeanPhzAngSrt(ind,2)'),                      ...
            10,'b',smarkers(s),'filled');
end    
Lines([],0,'k');
Lines(0,[],'k');    
colormap(hsv(3));
caxis([0,2*pi]);
xlim([0,15]);
%Lines(4.5,[],'k');
xlabel('Firing Rate Difference (S_1 - S_2)')
ylabel('Phase Preference Difference (S_2 - S_1)')

% PLOT rate diff vs phase diff
hold('on');
ucnt = 1;
% EXAMPLES 
for tind = 1:size(tppsUnits{g},1)
    t = tppsUnits{g}{tind,1}(1);
    dispOpts = {['.'],                              ...
                'MarkerFaceColor','m',                            ...
                'MarkerEdgeColor','m',                            ...
                'MarkerSize',10};
    for uid = 1:numel(tppsUnits{g}{tind,2});
        u = tppsUnits{g}{tind,2}(uid);
        ind = ismember(cluSessionMapSubset,[t,u],'rows');% example unit index        
        plot(mprMeanRateSrt(ind,1)'-mprMeanRateSrt(ind,2)',       ...
             -circ_dist(mprMeanPhzAngSrt(ind,1)',                  ...
                       mprMeanPhzAngSrt(ind,2)'),                 ...
         dispOpts{:});
    end
    ucnt = ucnt+1;
end

axes(fax);
text(fax,fig.page.width/2,fig.page.height-fig.page.marginTop,tppsTitles{g},'FontSize',12);

% $$$ print(gcf,'-depsc2',...
% $$$       ['/storage/share/Projects/BehaviorPlaceCode/thetaPhasePrecession/',...
% $$$        ['unit_examples_DRZCTPP_',tppsNames{g} ,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/thetaPhasePrecession/',...
       ['unit_examples_DRZCTPP_',tppsNames{g} ,'.png']]);


%print tppsNames{g}
end%for g







cfig = figure();
while true,     
    figure(hfig);    
    pause(0.1);
    hax = gca(hfig);
    pause(0.1);    
    clf(cfig);
    cax = copyobj(hax,cfig); 
    set(cax,'Position',[0,0,10,5]);
    cax = copyobj(hax,cfig);     
    set(cax,'Position',[0,5,10,5]);    
    drawnow();
    figure(hfig);        
    waitforbuttonpress(); 
end