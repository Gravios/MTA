

qnt = load(fullfile(MTA_PROJECT_PATH,'analysis','MjgER2016_req20191104_qntlav.mat'));

!rm /storage/gravio/data/project/general/analysis/MjgER2016_req20191104_jpdf_2d_PS_6c0239fd1f24d4bd226d1a9d44277e2c.mat

% CA1 
trialIndex = 1:7;
clear('sBySes');
sBySes(1) = load('/storage/gravio/data/project/general/analysis/MjgER2016_req20191104_jpdf_2d_PS_11c7acec34049037cf62dc1324558890.mat');
%sBySes(end+1) = load('/storage/gravio/data/project/general/analysis/MjgER2016_req20191104_jpdf_2d_PS_053796b66bf5c99fd9a69b898770bfae.mat');
sBySes(end+1) = load('/storage/gravio/data/project/general/analysis/MjgER2016_req20191104_jpdf_2d_PS_466e97cdf155c3f48319dd13a5b485e7.mat');
fwdOffset = -25;
latOffset = 8;
phzAsc = 2;
phzDsc = 6;
phzShift = 0;
cblen = dblen+20*double(ismember(dtind,[3,4,5]));


% CA3 
trialIndex = 1:4;
% hba X hvl
clear('sBySes');
sBySes(1) = load(['/storage/gravio/data/project/general/analysis/', ...
                  'MjgER2016_req20191104_jpdf_2d_PS_ca3_11c7acec34049037cf62dc1324558890.mat']);
% hba X hvf
sBySes(2) = load(['/storage/gravio/data/project/general/analysis/',...
                  'MjgER2016_req20191104_jpdf_2d_PS_ca3_053796b66bf5c99fd9a69b898770bfae.mat']);
% $$$ sBySes(2) = load(['/storage/gravio/data/project/general/analysis/',...
% $$$                   'MjgER2016_req20191104_jpdf_2d_PS_ca3_37b2060b443d93adcab47ed0681c0059.mat']);
fwdOffset = -25;
latOffset = 0;
phzAsc = 5;
phzDsc = 2;
phzShift = 3;



%%%<<< Compute population mean-EPP over behavioral variables 
ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[9:11])),2)            ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<320;


vLbl = {};
vEds = {};
vCtr = {};
vInd = {};
vBinLbl = {};
varBinCnt = 3;

vLbl{end+1} = 'hba';
vEds{end+1} = [-1.2,-0.3,0.3,1.2];  % length: varBinCnt
vCtr{end+1} = mean([vEds{end}(2:end); vEds{end}(1:end-1)]);
vInd{end+1} = discretize(-(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]))),vEds{end});
vBinLbl{end+1} = {'Left','Center','Right'};

vLbl{end+1} = 'hvl';
vEds{end+1} = [-80,-5,5,80]; % length: varBinCnt
vCtr{end+1} = mean([vEds{end}(2:end); vEds{end}(1:end-1)]);
vInd{end+1} = discretize(dfhrvl,vEds{end});
vBinLbl{end+1} = {'CCW','0','CW'};

vLbl{end+1} = 'hvf';
vEds{end+1} = [-10,10,25,80]; % length: varBinCnt
vCtr{end+1} = mean([vEds{end}(2:end); vEds{end}(1:end-1)]);
vInd{end+1} = discretize(dfhrvf,vEds{end});
vBinLbl{end+1} = {'rest','slow','fast'};

vLbl{end+1} = 'hbd';
%vEds{end+1} = [100,125,145,165];
vEds{end+1} = [100,130,150,175];
vCtr{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vInd{end+1} = discretize(cblen,vEds{end});
vBinLbl{end+1} = {'short','medium','long'};



% COMPUTE JPDF(FERROR, PHZ | HBA, HRVL)
tdphz = dphz(ind);

out =  zeros([numel(ferrorBinCenters{e}),numel(phzBinCenters),varBinCnt, varBinCnt, 2,numel(vLbl),numel(vLbl)]);
for v = 1;%:numel(vLbl),
    for b = 1:numel(vLbl),
        if v==b, continue, end        
% SETUP conditioned vars
        edx = numel(vCtr{v});
        edy = numel(vCtr{b});
        txi = vInd{v}(ind);
        tyi = vInd{b}(ind);

        for e = 1:2;
% GET dependent var { ego-centric phase precession: forward, lateral }
            tferr = ferr{e}(ind);
            if e==2, tferr = tferr+latOffset; 
            else,    tferr = tferr+fwdOffset; 
            end
            
            
% COMPUTE conditional means
            for x = 1:edx
                for y = 1:edy,            
                    indb = tyi==y & txi==x;
                    out(:,:,x,y,e,v,b) = hist2([tferr(indb),tdphz(indb)],  ferrorBinEdges{e},  phzBins);
                end
            end
        end        
    end
end


mout =  zeros([1,numel(phzBinCenters),varBinCnt, varBinCnt, 2,numel(vLbl),numel(vLbl)]);
for v = 1,%:numel(vLbl),
    for b = 1:numel(vLbl),
        if v==b, continue, end
        for y = 1:edy,    
            for x = 1:edx
                for p = 1:numel(phzBinCenters),
                   for e = 1:2,
                        try
                            [g, index] = unique(cumsum(out(:,p,x,y,e,v,b)./sum(out(:,p,x,y,e,v,b))),'last'); 
                            index(g==0|g==1) = [];
                            g(g==0|g==1) = [];
                            mout(:,p,x,y,e,v,b) = interp1(g,ferrorBinCenters{e}(index),0.5);
                        end
                    end
                end
            end
        end
    end 
end
%%%>>>





%%%<<< STARTFIG ------------------------------------------------------------------------------


[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','portrait',[],1,1,0.0,0.0);

%       hba,hvl  hvl,hvf  hvf,hbp  hba,hvf
vGrp = { [1,2],   [1,3]};
%vGrp = { [1,2],   [1,4]};
cmap = 'bcr';
cmap = 'bgr';    

xl = [-150,150];
xls = xl;
bxpltOpts.plotstyle = 'traditional';

%%%<<< Forward PP HBAxHVL 
% Forward PP HBAxHVL JPDF
for g = 1:2;
e = 1;
v = vGrp{g}(1);
b = vGrp{g}(2);
[yind, yOffSet, xind, xOffSet] = deal(y,                         ...
                                      0,                         ...
                                      1,                         ...
                                      (g-1)*1.5);

for x = 1:edx,
    for y = 1:edy,                
        [yind, yOffSet, xind, xOffSet] = deal(y,                                        ...
                                              0,           ...
                                              1+6*(g-1)+(x-1)*2,                         ...
                                              xOffSet);
        sax(end+1) = axes('Units','centimeters',                                        ...
                          'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                            fig.page.ypos(yind)+yOffSet,                      ...
                            fig.subplot.width*2 ,                               ...
                            fig.subplot.height],                            ...
                          'FontSize', 8,                                                ...
                          'LineWidth',1);
        hold(sax(end),'on');
        imagesc(ferrorBinEdges{e},                                           ...
                [phzBins,phzBins+2*pi],                                      ...
                imgaussfilt(circshift(repmat(out(:,:,x,4-y,e,v,b),[1,2]),phzShift,2),[3,0.2])');
        axis('tight');  axis('xy');  ylim([0,2.*pi]);
        xlim(xl);  
        colormap(sax(end),'jet');
% $$$         Lines([],pi,'k','--',1);
        Lines(0,[],'k','--',1);
        sax(end).XTickLabel = {};
        sax(end).YTickLabel = {};
        if y == 2 && x == 1 && e == 1 && g==1,
            ylabel('Forward PP');
        end
        if x==2&&y==1&&e==1,
            title(sax(end),[vLbl{vGrp{g}(1)},' vs ',vLbl{vGrp{g}(2)}]);
        end
        if x==3,
            sax(end).YAxisLocation = 'right';
            yh = ylabel(sax(end),vBinLbl{b}(4-y),'Rotation',0,'Color',cmap(4-y));
            %yh = ylabel(sax(end),vBinLbl{b}(4-y),'Rotation',0,'Color',cmap(4-y));            
            yh.Position(1:2) = [yh.Position(1)+65,4];
            %uistack(yh,'top');
        end
        plot(max(xlim()).*0.5*cos(linspace(-pi,pi,100)),linspace([ylim(),100]),'w','LineWidth',2)
    end
end



[yind, yOffSet, xind, xOffSet] = deal(yind+1,yOffSet-0.25,1+6*(g-1)+(x-1)*2,xOffSet);
% Forward PP HBAxHVL Ascending Phase session distribution
for x = 1:edx;
    [yind, yOffSet, xind, xOffSet] = deal(yind, yOffSet, 1+6*(g-1)+(x-1)*2, xOffSet);    
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width*2 ,                             ...
                                  fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    
    if x==3, blbls = vBinLbl{vGrp{g}(2)};
    else,    blbls = {'','',''};
    end

    hold(sax(end),'on');
    ylim([0.5,3.5]);    
    %Lines(0,[],'k',[],1);    
    grid(sax(end),'on');
    
    boxplot(reshape(sq(mean(sBySes(g).smJpdf(2,phzOrder(phzAsc),x,:,e,trialIndex,:),ndims(sBySes(g).smJpdf)))'+fwdOffset,[],1),...
            reshape(ones([numel(trialIndex),1]).*[1:3],[],1),...
            'ori','horizontal',...
            'Colors','bgr',...
            'plotstyle', bxpltOpts.plotstyle,...
            'Labels', blbls);
    sax(end).Position = [fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width*2 ,                             ...
                                  fig.subplot.height];
    sax(end).YAxisLocation = 'right';
    sax(end).XTickLabel = '';        
    xlim([xls]);
    %Lines(0,[],'k',[],1);    
    %ylabel({'\theta @',['290',char(176)]});        

    if x == 1 && g == 1,
        axes(fax);        
        text(sax(end).Position([1])-0.5,... x
             sum(sax(end).Position([2,4]).*[1,0.5]),... y
             {'Phase',['290',char(176)]},... t
             'Rotation',90, ...
             'VerticalAlignment','middle',...
             'HorizontalAlignment','center'); 
    end

end


x = 1;
[yind, yOffSet, xind, xOffSet] = deal(yind+1,yOffSet,1+6*(g-1)+(x-1)*2,xOffSet);
e = 1
% Forward PP HBAxHVL Descending Phase session distribution
for x = 1:edx;
    [yind, yOffSet, xind, xOffSet] = deal(yind, yOffSet, 1+6*(g-1)+(x-1)*2, xOffSet);    
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width*2 ,                             ...
                                  fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    
    if x==3, blbls = vBinLbl{vGrp{g}(2)};
    else,    blbls = {'','',''};
    end
    
    hold(sax(end),'on');
    ylim([0.5,3.5]);    
    %Lines(0,[],'k',[],1);    
    grid(sax(end),'on');    
    
    boxplot(reshape(sq(mean(sBySes(g).smJpdf(2,phzOrder(phzDsc),x,:,e,trialIndex,:),ndims(sBySes(g).smJpdf)))'+fwdOffset,[],1),...
            reshape(ones([numel(trialIndex),1]).*[1:3],[],1),...
            'ori','horizontal',...
            'Colors','bgr',...
            'plotstyle', bxpltOpts.plotstyle,...
            'Labels', blbls);
    sax(end).Position = [fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width*2 ,                             ...
                                  fig.subplot.height];
    sax(end).YAxisLocation = 'right';
    sax(end).XTickLabel = '';    
    xlim(xls);
    %Lines(0,[],'k',[],1);    
    %ylabel({'\theta @',['115',char(176)]});

    if x == 1 && g == 1,    
        axes(fax);
        text(sax(end).Position([1])-0.5,... x
             sum(sax(end).Position([2,4]).*[1,0.5]),... y
             {'Theta',['115',char(176)]},... t
             'Rotation',90, ...
             'VerticalAlignment','middle',...
             'HorizontalAlignment','center'); 
    end
end
%%%>>>

%%%<<< Lateral PP HBAxHVL 
% Lateral PP HBAxHVL JPDF
[yind, yOffSet, xind, xOffSet] = deal(yind,yOffSet-0.25,1+6*(g-1)+(x-1)*2,xOffSet);
yio = yind;
e = 2;
for x = 1:edx
    for y = 1:edy,                
        [yind, yOffSet, xind, xOffSet] = deal(y+yio,                                        ...
                                              yOffSet,           ...
                                              1+6*(g-1)+(x-1)*2,                         ...
                                              xOffSet);
        sax(end+1) = axes('Units','centimeters',                                        ...
                          'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                            fig.page.ypos(yind)+yOffSet,                      ...
                            fig.subplot.width*2 ,                               ...
                            fig.subplot.height],                            ...
                          'FontSize', 8,                                                ...
                          'LineWidth',1);
        hold(sax(end),'on');
        imagesc(ferrorBinEdges{e},                                           ...
                [phzBins,phzBins+2*pi],                                      ...
                imgaussfilt(circshift(repmat(out(:,:,x,4-y,e,v,b),[1,2]),phzShift,2),[3,0.2])');
        axis('tight');  axis('xy');  ylim([0,2.*pi]);
        xlim(xl);  
        colormap(sax(end),'jet');
% $$$         Lines([],pi,'k','--',1);
        Lines(0,[],'k','--',1);

% $$$         hln = Lines([],pi,'k',[],1);
% $$$         hln.LineWidth = 1;            
% $$$         hln = Lines(0,[],'k',[],1);
% $$$         hln.LineWidth = 1;            
        sax(end).XTickLabel = {};
        sax(end).YTickLabel = {};
        if y == 2 && x == 1 && e == 2 && g==1,
            ylabel('Lateral PP');
        end
        if x==3,
            sax(end).YAxisLocation = 'right';
            yh = ylabel(sax(end),vBinLbl{b}(4-y),'Rotation',0,'Color',cmap(4-y));
            %yh = ylabel(sax(end),vBinLbl{b}(4-y),'Rotation',0,'Color',cmap(4-y));            
            yh.Position(1:2) = [yh.Position(1)+65,4];
            %uistack(yh,'top');
        end
        plot(max(xlim()).*0.5*cos(linspace(-pi,pi,100)),linspace([ylim(),100]),'w','LineWidth',2)        
    end
end


% Lateral PP HBAxHVL Ascending Phase session distribution
[yind, yOffSet, xind, xOffSet] = deal(yind+1,yOffSet-0.25,1+6*(g-1)+(x-1)*2,xOffSet);
for x = 1:edx;
    [yind, yOffSet, xind, xOffSet] = deal(yind, yOffSet, 1+6*(g-1)+(x-1)*2, xOffSet);    
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width*2 ,                             ...
                                  fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    
    if x==3, blbls = vBinLbl{vGrp{g}(2)};
    else,    blbls = {'','',''};
    end
    
    hold(sax(end),'on');
    ylim([0.5,3.5]);    
    %Lines(0,[],'k',[],1);    
    grid(sax(end),'on');    
    
    boxplot(reshape(sq(mean(sBySes(g).smJpdf(2,phzOrder(phzAsc),x,:,e,trialIndex,:),ndims(sBySes(g).smJpdf)))'+latOffset,[],1),...
            reshape(ones([numel(trialIndex),1]).*[1:3],[],1),...
            'ori','horizontal',...
            'Colors','bgr',...
            'plotstyle', bxpltOpts.plotstyle,...
            'Labels', blbls);
    sax(end).Position = [fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width*2 ,                             ...
                                  fig.subplot.height];
    sax(end).YAxisLocation = 'right';
    sax(end).XTickLabel = '';        
    xlim([xls]);
    %ylabel({'\theta @',['290',char(176)]});    
    if x == 1 && g == 1,
        axes(fax);            
        text(sax(end).Position([1])-0.5,... x
             sum(sax(end).Position([2,4]).*[1,0.5]),... y
             {'Phase',['290',char(176)]},... t
             'Rotation',90, ...
             'VerticalAlignment','middle',...
             'HorizontalAlignment','center'); 
    end
end

% Lateral PP HBAxHVL Descending Phase session distribution
[yind, yOffSet, xind, xOffSet] = deal(yind+1,yOffSet,1+6*(g-1)+(x-1)*2,xOffSet);
for x = 1:edx;
    [yind, yOffSet, xind, xOffSet] = deal(yind, yOffSet, 1+6*(g-1)+(x-1)*2, xOffSet);    
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width*2 ,                             ...
                                  fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    
    if x==3, blbls = vBinLbl{vGrp{g}(2)};
    else,    blbls = {'','',''};
    end
    
    hold(sax(end),'on');
    ylim([0.5,3.5]);
    %Lines(0,[],'k',[],1);
    grid(sax(end),'on');
    boxplot(reshape(sq(mean(sBySes(g).smJpdf(2,phzOrder(phzDsc),x,:,e,trialIndex,:),ndims(sBySes(g).smJpdf)))'+latOffset,[],1),...
            reshape(ones([numel(trialIndex),1]).*[1:3],[],1),...
            'ori','horizontal',...
            'Colors','bgr',...
            'plotstyle', bxpltOpts.plotstyle,...
            'Labels', blbls);
    sax(end).Position = [fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width*2 ,                             ...
                                  fig.subplot.height];
    sax(end).YAxisLocation = 'right';
    sax(end).XTickLabel = '';    
    xlim(xls);
    %ylabel({'\theta @',['115',char(176)]});
    if x == 1 && g == 1,
        axes(fax);            
        text(sax(end).Position([1])-0.5,... x
             sum(sax(end).Position([2,4]).*[1,0.5]),... y
             {'Theta',['115',char(176)]},... t
             'Rotation',90, ...
             'VerticalAlignment','middle',...
             'HorizontalAlignment','center'); 
    end
% $$$     hold(sax(end),'on');
% $$$     plot(sq(mout(1,phzOrder(6),x,:,e,v,b)),[1,2,3],'*m');
end
%%%>>>

%%%<<< Egocentric asc-des vector HBAxHVL
[yind, yOffSet, xind, xOffSet] = deal(yind+3,yOffSet-0.25,1+6*(g-1)+(x-1)*2,xOffSet);

for x = 1:edx,
    [yind, yOffSet, xind, xOffSet] = deal(yind, ...
                                          yOffSet,...
                                          1+6*(g-1)+(x-1)*2, ...
                                          xOffSet);
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width*2,   ...
                                  fig.subplot.height*3],...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');
    sax(end).UserData = [vLbl{v},'_',vLbl{b}];

    for i = 1:varBinCnt,
        plot(sq(mout(1,phzOrder(fliplr([phzAsc,phzDsc])),x,i,2,v,b))+latOffset,...
             sq(mout(1,phzOrder(fliplr([phzAsc,phzDsc])),x,i,1,v,b)),[cmap(i),'-'],'LineWidth',1.5)
        plot(sq(mout(1,phzOrder(fliplr(phzDsc)),  x,i,2,v,b))+latOffset,...
             sq(mout(1,phzOrder(fliplr(phzDsc))  ,x,i,1,v,b)),[cmap(i),'o'],'LineWidth',2)
    end
    xlim([xls]);
    ylim([-60,75]);
% $$$     Lines(-offset,[],'k',[],1);
% $$$     Lines([],0,'k',[],1);
% $$$     Lines(offset,[],'k',[],1);
% $$$     Lines(0,[],'k',[],1);
    sax(end).XTick = [-75,0,75];
    sax(end).XTickLabel = {-7.5,0,7.5};    
    sax(end).YTick = [-50,0,50];
    sax(end).YTickLabel = {-5, 0, 5};
    sax(end).YAxisLocation = 'right';
    grid(sax(end),'on');    
    
    axes(fax);
    if x < 3,
        sax(end).YTickLabel = '';        
    end
    
    text(sum(sax(end).Position([1,3]).*[1,0.5]),                                    ... x
         sax(end).Position([2])-1,                                                  ... y
         vBinLbl{1}{x},                                                             ... text
         'HorizontalAlignment', 'center');
end

end



af(@(h) set(h,'LineWidth',1),sax);
axes(fax);
hrct = af(@(h) rectangle('Position',h.Position,'LineWidth',1),sax);


% $$$ axes(sax(end));
% $$$ yOffSet = -0.5;
% $$$ xOffSet = fig.subplot.width/2;
% $$$ for s = 1:numel(stsLbls),
% $$$     [yind, yOffSet, xind, xOffSet] = deal(6,yOffSet,s,xOffSet);
% $$$     sax(end+1) = axes('Units','centimeters',                                        ...
% $$$                       'Position',[fig.page.xpos(xind)+xOffSet,                      ...
% $$$                                   fig.page.ypos(yind)+yOffSet,                      ...
% $$$                                   fig.subplot.width ,                               ...
% $$$                                   fig.subplot.height],                              ...
% $$$                       'FontSize', 8,                                                ...
% $$$                       'LineWidth',1);
% $$$     hold(sax(end),'on');
% $$$     imagesc(sq(mean(qnt.qntlav( 2, 3, 2:end-1,1:5,s,1,:),7,'omitnan'))');
% $$$     sax(end).XTickLabel = {};
% $$$     sax(end).YTickLabel = {};
% $$$     colormap(sax(end),'jet');caxis([-30,120]);title(qnt.stsLbls{s});axis('tight');axis('xy');
% $$$ end
% $$$ 
% $$$ yOffSet = -0.5;
% $$$ xOffSet = fig.subplot.width/2;
% $$$ for s = 1:numel(stsLbls),
% $$$     [yind, yOffSet, xind, xOffSet] = deal(7,yOffSet,s,xOffSet);
% $$$     sax(end+1) = axes('Units','centimeters',                                        ...
% $$$                       'Position',[fig.page.xpos(xind)+xOffSet,                      ...
% $$$                                   fig.page.ypos(yind)+yOffSet,                      ...
% $$$                                   fig.subplot.width ,                               ...
% $$$                                   fig.subplot.height],                              ...
% $$$                       'FontSize', 8,                                                ...
% $$$                       'LineWidth',1);
% $$$     hold(sax(end),'on');
% $$$     imagesc(sq(mean(qnt.qntlav( 2, 3, 2:end-1,1:5,s,2,:),7,'omitnan'))');    
% $$$     sax(end).XTickLabel = {};
% $$$     sax(end).YTickLabel = {};
% $$$     colormap(sax(end),'jet');caxis([-70,70]);title(qnt.stsLbls{s});axis('tight');axis('xy');
% $$$ end


%%%>>>



