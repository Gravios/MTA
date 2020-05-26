

qnt = load(fullfile(MTA_PROJECT_PATH,'analysis','MjgER2016_req20191104_qntlav.mat'));

sBySesHbaHvl = load('/storage/gravio/data/project/general/analysis/MjgER2016_req20191104_jpdf_2d_PS_11c7acec34049037cf62dc1324558890.mat');

fwdOffset = 25;
latOffset = 8;

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
%vEds{end+1} = linspace(0,1.2,6);
vCtr{end+1} = mean([vEds{end}(2:end); vEds{end}(1:end-1)]);
vInd{end+1} = discretize(-(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]))),vEds{end});
%vBinLbl{end+1} = {'L','CL','C','CR','R'};
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
            if e==2, tferr = tferr+latOffset; end
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

[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','portrait',[],1.3,1.3,0.0,0.0);

%       hba,hvl  hvl,hvf  hvf,hbp  hba,hvf
vGrp = { [1,2],   [3,2],   [3,1]};
vGrp = { [1,2],   [1,3]};
cmap = 'bcr';
cmap = 'bgr';    
for g = 1:numel(vGrp);
    v = vGrp{g}(1); 
    b = vGrp{g}(2); 
    % PLOT JPDF(FERROR, PHZ | HBA, HRVL)
    [yind, yOffset, xind, xOffset] = deal(1,0, 1, 0);
    for e = 1:2,
        for x = 1:edx
            for y = 1:edy,                
                [yind, yOffSet, xind, xOffSet] = deal(y,                                        ...
                                                      y*fig.subplot.height/2-(g-1)*4,           ...
                                                      x+3*double(e==2),                         ...
                                                      fig.subplot.width/3+double(e==2)./4-1);
                sax(end+1) = axes('Units','centimeters',                                        ...
                                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                                    fig.page.ypos(yind)+yOffSet,                      ...
                                    fig.subplot.width ,                               ...
                                    fig.subplot.height/2],                            ...
                                  'FontSize', 8,                                                ...
                                  'LineWidth',1);
                hold(sax(end),'on');
                imagesc(ferrorBinEdges{e},                                           ...
                        [phzBins,phzBins+2*pi],                                      ...
                        imgaussfilt(repmat(out(:,:,x,4-y,e,v,b),[1,2]),[3,0.2])');
                axis('tight');  axis('xy');  ylim([0,2.*pi]);
                if e == 1,  xlim([-250,250]);  else  xlim([-250,250]);  end
                colormap(sax(end),'jet');
                Lines([],pi,'k');
                Lines(0,[],'k');
                sax(end).XTickLabel = {};
                sax(end).YTickLabel = {};
                if y == 1 && x == 2 && e == 1 && g==1,
                    title('Forward');
                elseif y == 1 && x == 2 && e == 2 && g==1,
                    title('Lateral');
                end
                if y==3,
                    xlabel(sax(end),vBinLbl{v}(x));
                end
                if x==1&&y==2&&e==1,
                    ylabel(sax(end),[vLbl{vGrp{g}(1)},' vs ',vLbl{vGrp{g}(2)}]);
                end
                if x==3 & e==2,
                    sax(end).YAxisLocation = 'right';
                    yh = ylabel(sax(end),vBinLbl{b}(4-y),'Rotation',0,'Color',cmap(4-y));
                    yh.Position(1:2) = [yh.Position(1)+220,6];
                    %uistack(yh,'top');
                end
            end
        end
    end


    [yind, yOffSet, xind, xOffSet] = deal(y, ...
                                          y*fig.subplot.height/2-(g-1)*4,...
                                          x+4*double(e==2), ...
                                          fig.subplot.width/2+double(e==2));
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                        fig.page.ypos(yind)+yOffSet,                      ...
                        fig.subplot.width*3+fig.subplot.horizontalPadding.*2.5,   ...
                        fig.subplot.height/2*4+2*fig.subplot.verticalPadding],...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');
    sax(end).UserData = [vLbl{v},'_',vLbl{b}];

    ofs = 75;
    offset = [-ofs:ofs:ofs];
    for j = 1:varBinCnt,
        for i = 1:varBinCnt,
            plot(sq(mout(1,phzOrder(fliplr([2,6])),j,i,2,v,b))+offset(j),...
                 sq(mout(1,phzOrder(fliplr([2,6])),j,i,1,v,b)),[cmap(i),'-'],'LineWidth',1.5)
            plot(sq(mout(1,phzOrder(fliplr(6)),  j,i,2,v,b))+offset(j),...
                 sq(mout(1,phzOrder(fliplr(6))  ,j,i,1,v,b)),[cmap(i),'o'],'LineWidth',2)
        end
    end
    xlim([offset(1)-diff(offset(1:2)),offset(end)+diff(offset(1:2))]);
    ylim([-80,150])
    Lines(-offset,[],'k');
    Lines([],0,'k');
    Lines(offset,[],'k');
    Lines(0,[],'k');
    sax(end).XTick = offset;
    sax(end).XTickLabel = vBinLbl{v};
    sax(end).YAxisLocation = 'right';
end


axes(findobj(sax,'UserData','hba_hvl'));
for x = 1:3;
    boxplot(reshape(sq(mean(sBySesHbaHvl.smJpdf(2,phzOrder([2]),x,:,2,1:7,:),ndims(sBySesHbaHvl.smJpdf)))'+8+offset(x),[],1),...
            reshape(ones([7,1]).*[1:3],[],1),...
            'ori','horizontal',...
            'positions',[110,125,140]',...
            'Colors','bgr',...
            'plotstyle', 'compact'...
            );
end

for x = 1:3;
    boxplot(reshape(sq(mean(sBySesHbaHvl.smJpdf(2,phzOrder([6]),x,:,2,1:7,:),ndims(sBySesHbaHvl.smJpdf)))'+8+offset(x),[],1),...
            reshape(ones([7,1]).*[1:3],[],1),...
            'ori','horizontal',...
            'positions',[-70,-55,-40]',...
            'Colors','bgr',...
            'plotstyle', 'compact'...
            );
end
xlim([offset(1)-diff(offset(1:2)),offset(end)+diff(offset(1:2))]);
ylim([-80,150])



axes(sax(end));
yOffSet = -0.5;
xOffSet = fig.subplot.width/2;
for s = 1:numel(stsLbls),
    [yind, yOffSet, xind, xOffSet] = deal(6,yOffSet,s,xOffSet);
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width ,                               ...
                                  fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');
    imagesc(sq(mean(qnt.qntlav( 2, 3, 2:end-1,1:5,s,1,:),7,'omitnan'))');
    sax(end).XTickLabel = {};
    sax(end).YTickLabel = {};
    colormap(sax(end),'jet');caxis([-30,120]);title(qnt.stsLbls{s});axis('tight');axis('xy');
end

yOffSet = -0.5;
xOffSet = fig.subplot.width/2;
for s = 1:numel(stsLbls),
    [yind, yOffSet, xind, xOffSet] = deal(7,yOffSet,s,xOffSet);
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width ,                               ...
                                  fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');
    imagesc(sq(mean(qnt.qntlav( 2, 3, 2:end-1,1:5,s,2,:),7,'omitnan'))');    
    sax(end).XTickLabel = {};
    sax(end).YTickLabel = {};
    colormap(sax(end),'jet');caxis([-70,70]);title(qnt.stsLbls{s});axis('tight');axis('xy');
end


%%%>>>



