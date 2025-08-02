% req20210705
% 
% XYHB 30sps 300ms

global MTA_PROJECT_PATH

configure_default_args();

MjgER2016_load_data();

sampleRate = 30;

states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};

%%%<<< MAIN FIGURE decoding ------------------------------------------------------------------------
t = 20; 
exyz = preproc_xyz(Trials{t},'trb');
exyz.resample(dc{t==tind}.sampleRate);
efet =  fet_HB_pitchB(Trials{t},dc{t==tind}.sampleRate);
ets = [1:size(exyz,1)]./sampleRate;
maskcirc = create_tensor_mask(dc{1}.pfs.adata.bins(1:2));
maskcirc = sq(mask(:,:,14,14));
maskbhv = sq(mask(11,11,:,:));

eind =   dc{t==tind}.ucnt >= 2;

%exampleRange = [3332,3417];
exampleRange = [8950,10100];
%exampleRange = [27770,28474];
exampleTimestamps = (dc{t==tind}.ind-1)/dc{t==tind}.sampleRate;
%eTS = [290,345];
%eTS = [530,600];
%eTS = [1000,1045];
eInds = find(WithinRanges(dc{t==tind}.ind,exampleRange));
exampleSubsetRange = [9101,9200];
exampleSubsetRange = [9650,9950];

eIndsSub = find( WithinRanges( dc{t==tind}.ind, exampleSubsetRange));

exampleTimestampsSubset = (dc{t==tind}.ind(eIndsSub)-1)./dc{t==tind}.sampleRate;

enanmaskSubset = double( dc{t==tind}.ucnt(eIndsSub) >= 3     ...
                 & dc{t==tind}.stcm(eIndsSub,1)==1);
enanmaskSubset(~enanmaskSubset)=nan;

enanmask = double( dc{t==tind}.ucnt(eInds) >= 3     ...
                 & dc{t==tind}.stcm(eInds,1)==1);
enanmask(~enanmask)=nan;


[hfig,fig,fax,sax] = set_figure_layout(figure(666007),'A4','landscape',[],2,2,0.1,0.1);

globalXOffset = 0;
globalYOffset = 1;
%%%<<< PLOT timeseries of xcoord vs decoded xcoord

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(1, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width.*5,                     ...
                              fig.subplot.height.*0.5],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(exampleRange(1):exampleRange(2)),...
     exyz(exampleRange,'hcom',1),...
     '-k','LineWidth',1);
plot(exampleTimestamps(eInds),...
     sq(dc{t==tind}.lax(eInds,1)).*enanmask,...
     '-r','LineWidth',1);
plot(exampleTimestampsSubset,...
     sq(dc{t==tind}.lax(eIndsSub,1)).*enanmaskSubset,...
     '-c','LineWidth',1);
sax(end).XTickLabels = [];
sax(end).YTick = [];
xlim(ets([exampleRange(1),exampleRange(2)]));
ylim(sax(end),[-500,500]);
box(sax(end),'on');
text(ets(exampleRange(1))+0.25,380,'X');
axes(fax);
line([sax(end).Position(1)-0.1].*[1,1],                                        ...
     sax(end).Position(2)+[0,sax(end).Position(4).*(200/diff(ylim(sax(end))))],...
     'Color','k',                                                              ...
     'LineWidth',1);
%%%>>>

%%%<<< PLOT timeseries of ycoord vs decoded ycoord
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2, 2, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width.*5,                     ...
                              fig.subplot.height.*0.5],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(exampleRange(1):exampleRange(2)),...
     exyz(exampleRange,'hcom',2),...
     '-k','LineWidth',1);
plot(exampleTimestamps(eInds),...
     sq(dc{t==tind}.lax(eInds,2)).*enanmask,...
     '-r','LineWidth',1);
plot(exampleTimestampsSubset,...
     sq(dc{t==tind}.lax(eIndsSub,2)).*enanmaskSubset,...
     '-c','LineWidth',1);
sax(end).XTickLabels = [];
sax(end).YTick = [];
xlim(sax(end),ets([exampleRange(1),exampleRange(2)]));
ylim(sax(end),[-500,500]);
box(sax(end),'on');
text(ets(exampleRange(1))+0.25,380,'Y');
axes(fax);
line([sax(end).Position(1)-0.1].*[1,1],                                        ...
     sax(end).Position(2)+[0,sax(end).Position(4).*(200/diff(ylim(sax(end))))],...
     'Color','k',                                                              ...
     'LineWidth',1);

%%%>>>

%%%<<< PLOT timeseries of head pitch vs decoded head pitch
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(3, 3, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width.*5,                     ...
                              fig.subplot.height.*0.5],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(exampleRange(1):exampleRange(2)),...
     efet(exampleRange,1),...
     '-k','LineWidth',1);
plot(exampleTimestamps(eInds),...
     sq(dc{t==tind}.lax(eInds,3)).*enanmask,...
     '-r','LineWidth',1);
plot(exampleTimestampsSubset,...
     sq(dc{t==tind}.lax(eIndsSub,3)).*enanmaskSubset,...
     '-c','LineWidth',1);
sax(end).XTickLabels = [];
sax(end).YTick = [];
xlim(ets([exampleRange(1),exampleRange(2)]));
ylim([-1.6,0.5]);
box(sax(end),'on');
text(ets(exampleRange(1))+0.25,0.15,'HP');

axes(fax);
line([sax(end).Position(1)-0.1].*[1,1],                                      ...
     sax(end).Position(2)+[0,sax(end).Position(4).*(0.5/diff(ylim(sax(end))))],...
     'Color','k',                                                            ...
     'LineWidth',1);

%%%>>>

%%%<<< PLOT timeseries of body pitch vs decoded body pitch

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(4, 4, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width.*5,                     ...
                              fig.subplot.height.*0.5],                 ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(exampleRange(1):exampleRange(2)),...
     efet(exampleRange,2),...
     '-k','LineWidth',1);
plot(exampleTimestamps(eInds),...
     sq(dc{t==tind}.lax(eInds,4)).*enanmask,...
     '-r','LineWidth',1);
plot(exampleTimestampsSubset,...
     sq(dc{t==tind}.lax(eIndsSub,4)).*enanmaskSubset,...
     '-c','LineWidth',1);

sax(end).XTickLabels = [];
sax(end).YTick = [];
xlim(ets([exampleRange(1),exampleRange(2)]));
ylim([-0.5,1.6]);
box(sax(end),'on');
text(ets(exampleRange(1))+0.25,1.25,'BP');

axes(fax);
line([sax(end).Position(1)-0.1].*[1,1],                                        ...
     sax(end).Position(2)+[0,sax(end).Position(4).*(0.5/diff(ylim(sax(end))))],  ...
     'Color','k',                                                              ...
     'LineWidth',1);


%%%>>>

%%%<<< PLOT timeseries of states

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(5, 4.5, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width*5,                      ...
                              fig.subplot.height.*0.75],                ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plotSTC(Trials{t}.stc,1,'text',states([1,fliplr(2:6)]),'kbbggr');
xlim(ets([exampleRange(1),exampleRange(2)]));
box(sax(end),'on');
xlabel(sax(end),'Seconds');

%%%>>>


%%%<<< PLOT example trajectory within physical space

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(1, 0, 6, 0.2);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,        ...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,        ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
imagesc(dc{t==tind}.pfs.adata.bins{1:2},~maskcirc');
colormap(sax(end),'gray')
caxis(sax(end),[-2,1])

circle(0,0,480,'-k');
plot(exyz(exampleSubsetRange,'hcom',1),...
     exyz(exampleSubsetRange,'hcom',2),...
     '-k','LineWidth',1);
plot(sq(dc{t==tind}.lax(eIndsSub,1)).*enanmaskSubset,...
     sq(dc{t==tind}.lax(eIndsSub,2)).*enanmaskSubset,...
     '-c','LineWidth',1);
xlim(sax(end),[-500,500]);
ylim(sax(end),[-500,500]);
line(sax(end),...
     [280,480],...
     [480,480],...
     'Color','k',...
     'LineWidth',1);
line(sax(end),...
     [480,480],...
     [280,480],...
     'Color','k',...
     'LineWidth',1);
sax(end).XTick = [];
sax(end).YTick = [];
box(sax(end),'on');

%%%>>>

%%%<<< PLOT example trajectory within behavior space
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2, 0, 6, 0.2);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,        ...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,        ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
imagesc(dc{1}.pfs.adata.bins{3:4},~maskbhv');
colormap(sax(end),'gray')
caxis(sax(end),[-2,1])

axis('xy');
plot(efet(exampleSubsetRange,1),...
     efet(exampleSubsetRange,2),...
     '-k','LineWidth',1);
plot(sq(dc{t==tind}.lax(eIndsSub,3)).*enanmaskSubset,...
     sq(dc{t==tind}.lax(eIndsSub,4)).*enanmaskSubset,...
     '-c','LineWidth',1);
xlim(sax(end),[-1.6,0.5]);
ylim(sax(end),[-0.5,1.6]);
line(sax(end),...
     [-0.15,0.45],...
     [1.55,1.55],...
     'Color','k',...
     'LineWidth',1);
line(sax(end),...
     [0.45,0.45],...
     [1.05,1.55],...
     'Color','k',...
     'LineWidth',1);
sax(end).XTick = [];
sax(end).YTick = [];
box(sax(end),'on');

%%%>>>


%%%<<< PLOT HB error conditioned on XY

% ADJUST subplot coordinates

[yind, yOffSet, xind, xOffSet] = deal(1, 0, 7, 0.5);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,                      ...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');

ind =   nniz(dca.posi(:,1:2)) ...
      & dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2|dca.stcm(:,3)==3|dca.stcm(:,4)==4|dca.stcm(:,5)==5|dca.stcm(:,6)==6) ...
        & dca.ucnt >= 2;
ind(ind) = logical(randn([sum(ind),1])>0.65);
out = accumarray(dca.posi(ind,1:2),sqrt(sum(dca.ecom(ind,3:4).^2,2)),[20,20],@mean);
eCount = accumarray(dca.posi(ind,1:2),ones([sum(ind),1]),[20,20],@sum);
out(eCount<100) = nan;
circle(0,0,480,'-k');

set(pcolor(dca.pfs(1).adata.bins{1}-25,...
           dca.pfs(1).adata.bins{2}-25,out'),'EdgeColor','none');
axis('xy');
xlim([-500,500]);
ylim([-500,500]);
line(sax(end),...
     [280,480],...
     [480,480],...
     'Color','k',...
     'LineWidth',1);
line(sax(end),...
     [480,480],...
     [280,480],...
     'Color','k',...
     'LineWidth',1);
xlabel(sax(end),'Seconds');
colormap(sax(end),'jet');
sax(end).XTick = [];
sax(end).YTick = [];
box(sax(end),'on');
cax = colorbar(sax(end),'eastoutside');
cax.Units = 'Centimeters';
cax.Position(1) = fig.page.xpos(xind)+xOffSet+globalXOffset+0.1 + fig.subplot.width;
cax.Position(2) = fig.page.ypos(yind)+yOffSet+globalYOffset;
cax.Label.String = 'rad';
%%%>>>

%%%<<< PLOT XY error conditioned on HB

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2, 0, 7, 0.5);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,                      ...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');

ind =   nniz(dca.feti(:,1:2)) ...
      & dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2|dca.stcm(:,3)==3|dca.stcm(:,4)==4|dca.stcm(:,5)==5|dca.stcm(:,6)==6) ...
      & dca.ucnt >= 2;
ind(ind) = logical(randn([sum(ind),1])>0.65);
out = accumarray(dca.feti(ind,1:2),sqrt(sum(bsxfun(@plus,dca.elax(ind,1:2),[-35,0]).^2,2)),[28,28],@mean);
eCount = accumarray(dca.feti(ind,1:2),ones([sum(ind),1]),[28,28],@sum);
out(eCount<100) = nan;

set(pcolor(dca.pfs(1).adata.bins{3},...
           dca.pfs(1).adata.bins{4},out'/10),'EdgeColor','none');
axis('xy');
xlim([-1.6,0.5]);
ylim([-0.5,1.6]);
line(sax(end),...
     [-0.15,0.45],...
     [1.55,1.55],...
     'Color','k',...
     'LineWidth',1);
line(sax(end),...
     [0.45,0.45],...
     [1.05,1.55],...
     'Color','k',...
     'LineWidth',1);
colormap(sax(end),'jet');
sax(end).XTick = [];
sax(end).YTick = [];
box(sax(end),'on');
cax = colorbar(sax(end),'eastoutside');
cax.Units = 'Centimeters';
cax.Position(1) = fig.page.xpos(xind)+xOffSet+globalXOffset+0.1 + fig.subplot.width;
cax.Position(2) = fig.page.ypos(yind)+yOffSet+globalYOffset;
cax.Label.String = 'cm';
%%%>>>

%%%<<< PLOT XY error as function of included units
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2, 0, 9, 0.1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,                      ...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,                      ...
                              fig.subplot.width,                                ...
                              fig.subplot.height],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');

ind = dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2 | dca.stcm(:,3)==3 | dca.stcm(:,4)==4 | dca.stcm(:,5)==5 | dca.stcm(:,6)==6);
hist2([dca.ucnt(ind),sqrt(sum(bsxfun(@plus,dca.elax(ind,1:2), [-35,0]).^2,2))/10],...
      linspace(0.5,20.5,21),...
      linspace(0,80,41),'');
sax(end).YTick = [0,20,40,60];
axis(sax(end),'tight');
colormap(sax(end),'jet');
xlabel(sax(end),'Units');
%%%>>>

%%%<<< PLOT HB error as function of included units
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(1, 0, 9, 0.1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                         ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,                       ...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,                       ...
                              fig.subplot.width,                                 ...
                              fig.subplot.height],                               ...
                  'FontSize', 8,                                                 ...
                  'LineWidth',1);
hold(sax(end),'on');

ind = dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2 | dca.stcm(:,3)==3 | dca.stcm(:,4)==4 | dca.stcm(:,5)==5 | dca.stcm(:,6)==6);
hist2([dca.ucnt(ind),sqrt(sum(dca.elax(ind,3:4).^2,2))],...
      linspace(0.5,20.5,21),...
      linspace(0,1.5,41),...
      '');
sax(end).YTick = [0,0.5,1,1.5];
sax(end).XTickLabel = [];
axis(sax(end),'tight');
colormap(sax(end),'jet');
%%%>>>



%%%<<< PLOT JPDF of Uniform dist of Pos and Bhv

% discretize the uniform distributions of position and behavior
ind =   dca.stcm(:,1)==1 ...
      & dca.stcm(:,7)~=7 ...
      & dca.stcm(:,8)~=8;
      
nEdges = 11;
nBins = nEdges-1;
[uPosErr,xPosErr] = MakeUniformDistr(sqrt(sum(bsxfun(@plus,dca.elax(ind,1:2), [-35,0]).^2,2))/10);
[uBhvErr,xBhvErr] = MakeUniformDistr(sqrt(sum(dca.elax(ind,3:4).^2,2)));
uPosErrInd = discretize(uPosErr,linspace([xPosErr([1,end])',nEdges]));
uBhvErrInd = discretize(uBhvErr,linspace([xBhvErr([1,end])',nEdges]));

nXTicks = 2;
nYTicks = 2;
xticks = round(xPosErr(round([1:size(xPosErr,1)/(nBins/nXTicks):size(xPosErr,1),size(xPosErr,1)])),0);
yticks = round(xBhvErr(round([1:size(xBhvErr,1)/(nBins/nYTicks):size(xBhvErr,1),size(xBhvErr,1)])),2);

% JPDF of uniform Pos and Bhv error
nind = nniz(uPosErrInd) & nniz(uBhvErrInd);
outE = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],ones([sum(nind),1]),[nBins,nBins],@sum);
outE = outE./sum(outE(:),'omitnan');

% Mean Number of coactive units
ucntErr = dca.ucnt(ind);
nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(ucntErr);
outC = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],...
                 ucntErr(nind),...
                 [nBins,nBins],...
                 @mean);

% Mean Pos Info
uposInfoM = dca.PosInfoM(ind);
nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(uposInfoM);
outP = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],uposInfoM(nind),[nBins,nBins],@mean);

% Mean Bhv Info
ubhvInfoM = dca.BhvInfoM(ind);
nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(ubhvInfoM);
outB = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],ubhvInfoM(nind),[nBins,nBins],@mean);

%%%<<< PLOT JPDF of uniform marginals of bhv and spc error
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(1,0, 10, 1.5);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                         ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,         ...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,         ...
                              fig.subplot.width,                                 ...
                              fig.subplot.height],                               ...
                  'FontSize', 8,                                                 ...
                  'LineWidth',1);
hold(sax(end),'on');
imagesc(outE')
sax(end).XTick = 0:nXTicks:nBins;
sax(end).YTick = 0:nYTicks:nBins;
sax(end).XTickLabel = [];
sax(end).YTickLabel = yticks;
box(sax(end),'on');
ylabel(sax(end),'Error (rad)');
axis('xy');
tpos = sax(end).Position;
colormap(sax(end),'jet');
cax = colorbar(sax(end));
sax(end).Position = tpos;
cax.Units = 'centimeters';
drawnow();
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
ylabel(cax,{'Probability'});
axis('tight');
%%%>>>

%%%<<< PLOT conditional mean included units
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2,0, 10, 1.5);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                         ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,                       ...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,                       ...
                              fig.subplot.width,                                 ...
                              fig.subplot.height],                               ...
                  'FontSize', 8,                                                 ...
                  'LineWidth',1);
hold(sax(end),'on');
imagesc(outC')
set(sax(end),'XTick',0:nXTicks:nBins);
set(sax(end),'YTick',0:nYTicks:nBins);
set(sax(end),'XTickLabel',[]);
set(sax(end),'YTickLabel',yticks);
box(sax(end),'on');
axis('xy');
tpos = sax(end).Position;
colormap(sax(end),'jet');
cax = colorbar(sax(end));
sax(end).Position = tpos;
cax.Units = 'centimeters';
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
ylabel(cax,{'Average','Unit Count'});
axis('tight');
%%%>>>

%%%<<< PLOT geometric conditional mean of bhv and spc info
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(3,0, 10, 1.5);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                         ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,                       ...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,                       ...
                              fig.subplot.width,                                 ...
                              fig.subplot.height],                               ...
                  'FontSize', 8,                                                 ...
                  'LineWidth',1);
hold(sax(end),'on');
imagesc(sqrt(outP.*outB)')
set(sax(end),'XTick',0:nXTicks:nBins);
set(sax(end),'YTick',0:nYTicks:nBins);
set(sax(end),'XTickLabel',xticks);
set(sax(end),'YTickLabel',yticks);
box(sax(end),'on');
xlabel('Error (cm)');

axis('xy');
tpos = sax(end).Position;
colormap(sax(end),'jet');
cax = colorbar(sax(end));
sax(end).Position = tpos;
cax.Units = 'centimeters';
ylabel(cax,{'Geometric Mean','(SI,BI)'})
drawnow();
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
axis('tight');

%%%>>>

%%%>>>

%%%<<< Save/LOAD intermediate figure
savefig(hfig,'/storage/share/Projects/BehaviorPlaceCode/decoding/decoding_pos_bhv_part1.fig');

%%%>>>


%%%<<< PLOT 
% Fwd phase precession theta -> states
% HP vs phase, theta
% BP vs phase, theta
%%%>>>

%%%<<< LOAD intermediate figure
hfig = open('/storage/share/Projects/BehaviorPlaceCode/decoding/decoding_pos_bhv_part2.fig');
[~,fig,fax,sax] = set_figure_layout(figure(666007),'A4','landscape',[],2,2,0.1,0.1);
hfig.Position(3:4) = [fig.page.width,fig.page.height];
close('666007');
%%%>>>

%%%<<< PLOT
% rear onset decoding error head and body pitch
% rear offset decoding error head and body pitch

% INSERT decoded value into full session length matrix

stsTransitions = {{'rear','hbhv'},...
                  {'hbhv','rear'},...
                  {'lbhv','rear'},...
                  {'hbhv','lbhv'},...
                  {'lbhv','hbhv'}};
stsTransitionLabels = {'High -> Rear','Rear -> High','Rear -> Low','Low -> High','High -> Low'};



for t = 1:28
    Trial = Trials{t};
    hper = [Trial.stc{'hloc+hpause'}];
    hper.fillgaps();
    hper.label = 'hbhv';
    hper.key   = 'u';
    hper.update_filename(Trial);
    hper.update_path(Trial);
    hper.update_hash();
    hper.save(1);

    lper = [Trial.stc{'lloc+lpause'}];
    lper.fillgaps();
    lper.label = 'lbhv';
    lper.key   = 'd';
    lper.update_filename(Trial);
    lper.update_path(Trial);
    lper.update_hash();
    lper.save(1);
end
cf(@(T) T.load('stc','msnn_ppsvd_raux'), Trials);



sTransitions = {'rear','hloc','hpause','lloc','lpause'};
sTransMatrix = zeros([numel(sTransitions),numel(sTransitions)]);
for t = 1:12
    Trial = Trials{tind(t)};
    disp(['[INFO] ',Trial.filebase]);
    
    pch = fet_HB_pitchB(Trial,30);
    
    for stst1 = 1:numel(sTransitions)
        for stst2 = 1:numel(sTransitions)        
            if stst1==stst2
                continue
            end
            stsTrans = Trial.stc.filter(30,{sTransitions{stst1},...
                                {'select_boarder_states',sTransitions{stst2},0.1}});
            sTransMatrix(stst1,stst2) = sTransMatrix(stst1,stst2)+numel(stsTrans{1});
        end
    end
end
sTransMatrix = sTransMatrix./sum(sTransMatrix(:));

% $$$ figure();
% $$$ sax = gobjects([0,1]);
% $$$ sax(end+1) = axes()
% $$$ imagesc(1:numel(sTransitions),1:numel(sTransitions),round(sTransMatrix/sum(sTransMatrix(:)),3));
% $$$ sax(end).XTick = [1:5];
% $$$ sax(end).YTick = [1:5];
% $$$ sax(end).XTickLabel = sTransitions;
% $$$ sax(end).YTickLabel = sTransitions;
% $$$ colormap(sax(end),'jet');
% $$$ 
% $$$ G = graph([1 1], [2 3]);
% $$$ e = G.Edges
% $$$ G = addedge(G,2,3)
% $$$ G = addnode(G,4)
% $$$ plot(G)
% $$$ 
% $$$ figure
% $$$ G=digraph(floor(sTransMatrix',2));
% $$$ ax = axes();
% $$$ gx = plot(ax,G);
% $$$ gx.MarkerSize =10;
% $$$ %ax.NodeCData = [1,0,0,0,1,0,0,1,0,0,0,1,0,0,1];
% $$$ gx.NodeCData = [5,3.15,3.15,1,1];
% $$$ %ax.NodeCData = 'rggbb';
% $$$ gx.NodeLabel = sTransitions;
% $$$ colormap(ax,'jet');
% $$$ caxis(ax,[0.8,5.5]);
% $$$ gx.LineWidth = 2;
% $$$ gx.XData = [2,3,2.5,1,1.5];
% $$$ gx.YData = [3,1.75,1,1.75,1];


esegs = {}; rsegs = {};
for t = 1:12
    Trial = Trials{tind(t)};
    disp(['[INFO] ',Trial.filebase]);
    
    pch = fet_HB_pitchB(Trial,30);
    
    rpitch = {};
    rpitch{1} = copy(pch);
    rpitch{1}.data = pch(:,1);
    rpitch{2} = copy(pch);
    rpitch{2}.data = pch(:,2);
    
    ind = dca.sessionId==tind(t) & dca.ucnt>2 & dca.BhvInfoM > 0.2;
    
    epitch{1} = copy(pch);
    epitch{1}.data = nan([size(pch,1),1]);
    %epitch{1}.data(dca.ind(ind)) = dca.lax(ind,3);
    epitch{1}.data(dca.ind(ind)) = dca.lom(ind,3);    
    
    epitch{2} = copy(pch);
    epitch{2}.data = nan([size(pch,1),1]);
    %epitch{2}.data(dca.ind(ind)) = dca.lax(ind,4);
    epitch{2}.data(dca.ind(ind)) = dca.lom(ind,4);
    
    for stst = 1:numel(stsTransitions)
        stsTrans = Trial.stc.filter(30,{stsTransitions{stst}{1},...
                                        {'duration',0.75},...
                                        {'select_boarder_states',stsTransitions{stst}{2},0.1}});
        esegs{t,stst} = cf(@(e)  e.segs(stsTrans{1}-31,61), epitch);
        rsegs{t,stst} = cf(@(r)  r.segs(stsTrans{1}-31,61), rpitch);
        
    end
end



% $$$ csegs = {}; 
% $$$ for t = 1:12
% $$$     Trial = Trials{tind(t)};
% $$$     disp(['[INFO] ',Trial.filebase]);
% $$$     
% $$$     ind = dca.sessionId==tind(t);
% $$$     
% $$$     pch = fet_HB_pitchB(Trial,30);
% $$$     
% $$$     eucount{1} = copy(pch);
% $$$     eucount{1}.data = nan([size(pch,1),1]);
% $$$     eucount{1}.data(dca.ind(ind)) = dca.ucnt(ind);
% $$$     
% $$$     for stst = 1:numel(stsTransitions)
% $$$         stsTrans = Trial.stc.filter(30,{stsTransitions{stst}{1},...
% $$$                                         {'duration',0.75},...
% $$$                                         {'select_boarder_states',stsTransitions{stst}{2},0.1}});
% $$$         csegs{t,stst} = cf(@(e)  e.segs(stsTrans{1}-31,61), eucount);
% $$$         
% $$$     end
% $$$ end
% $$$ 
% $$$ tbins = ([1:size(allcsegs,1)]-round(size(allcsegs,1))/2)/sampleRate;
% $$$ cbins = linspace(0.5,30.5,31);
% $$$ 
% $$$ figure();
% $$$ for stst = 1:5;
% $$$     allcsegs = [];
% $$$ for t = 1:12
% $$$     allcsegs = cat(2,allcsegs,csegs{t,stst}{1});
% $$$ end
% $$$ subplot(1,5,stst);
% $$$ out = [];
% $$$ for tm = 1:size(allcsegs,1)
% $$$     out(:,tm) = histcounts(allcsegs(tm,:),cbins);
% $$$ end
% $$$ sum(out)
% $$$ imagesc(tbins,cbins,out);
% $$$ hold('on')
% $$$ plot(tbins,mean(allcsegs,2,'omitnan'),'m')
% $$$ title(stsTransitionLabels{stst})
% $$$ axis('xy')
% $$$ end



hfig = open('/storage/share/Projects/BehaviorPlaceCode/decoding/decoding_pos_bhv_part2.fig');
[~,fig,fax,sax] = set_figure_layout(figure(666007),'A4','landscape',[],2,2,0.1,0.1);
hfig.Position(3:4) = [fig.page.width,fig.page.height];
close('666007');



ststPitchType = [1,2,3,4,5;
                 2,2,2,1,1];
xOffsets = [-0.1,-0.1,-0.1,1.2,1.1];

%hfig = figure();

sax = gobjects([1,0]);
plotCount = 0;
% $$$ for stst = 1:numel(stsTransitions)
% $$$ for pitchType = 1:2
for ststp = ststPitchType
stst = ststp(1);
pitchType = ststp(2);


%t = 1; 
allrsegs = [];
allesegs = [];
for t = 1:12
    allrsegs = cat(2,allrsegs,rsegs{t,stst}{pitchType});
    allesegs = cat(2,allesegs,esegs{t,stst}{pitchType});
end
if pitchType == 2,
    xlimits = [-0.25,1.25];
    ylimits = [-0.25,1.25];
else
    xlimits = [-1.6,0.25];
    ylimits = [-1.6,0.25];
% $$$     xlimits = prctile(allrsegs(nniz(allrsegs(:))),[5,95])+[-0.15,0.15];
% $$$     ylimits = prctile(allesegs(nniz(allesegs(:))),[5,95])+[-0.15,0.15];
end


% ADDPLOT : invisible plot which overlays onto the data scatterplot
sax(end+1) = axes('Units','Centimeters','Position',[fig.page.marginLeft+4*(stst-1)+xOffsets(stst)+1,8.75,3,3]);

    rax = sax(end);
    hold('on');
    for seg = 1:2:size(allrsegs,2),
        scatter(allrsegs(1:2:end,seg),...
                allesegs(1:2:end,seg),...
                1,...
                linspace(-1,1,ceil(size(allrsegs,1)./2)),...
                'Filled');
    end
    xlim(sax(end),xlimits);
    ylim(sax(end),ylimits);
    sax(end).XTickLabel = [];
    sax(end).YTickLabel = [];
    colormap(sax(end),'cool');
    daspect(sax(end),[1,1,1]);
    set(line(rax,[-2,2],[-2,2]),'Color','k','LineWidth',2);
    grid(rax,'on');
    box(rax,'on');
    % SERIOUSLY **** matlab
    drawnow(); pause(1);
    % SERIOUSLY 
% ADDPLOT : invisible plot which overlays onto the data scatterplot
sax(end+1) = axes('Units','centimeters','Position',rax.Position,'Visible','off');
    lax = sax(end);
    hold('on');
    scatter(median(allrsegs(1:2:end,:),2,'omitnan'),...
            median(allesegs(1:2:end,:),2,'omitnan'),...
            15,...
            linspace(-1,1,round(size(allrsegs,1)/2)),...
            'Filled','MarkerEdgeColor','k');
    plot(lax,median(allrsegs(1:2:end,:),2,'omitnan'),...
         median(allesegs(1:2:end,:),2,'omitnan')+std(allesegs(1:2:end,:),[],2,'omitnan'),...
         '-k',...
         'LineWidth',2);
    plot(lax,median(allrsegs(1:2:end,:),2,'omitnan'),...
         median(allesegs(1:2:end,:),2,'omitnan')-std(allesegs(1:2:end,:),[],2,'omitnan'),...
         '-k',...
         'LineWidth',2);
% $$$ plot(lax,median(allrsegs,2,'omitnan')+std(allrsegs,[],2,'omitnan'),...
% $$$      median(allesegs,2,'omitnan'),...
% $$$      '-k',...
% $$$      'LineWidth',2);
% $$$ plot(lax,median(allrsegs,2,'omitnan')-std(allrsegs,[],2,'omitnan'),...
% $$$      median(allesegs,2,'omitnan'),...
% $$$      '-k',...
% $$$      'LineWidth',2);
    set(sax(end),'Visible','off');
    colormap(lax,'jet');
    xlim(sax(end),xlimits);
    ylim(sax(end),ylimits);
    daspect(sax(end),[1,1,1]);
    % SERIOUSLY **** matlab
    drawnow(); pause(1);
    % SERIOUSLY 
% ADDPLOT : marginal of the observed pitch during state transition
sax(end+1) = axes('Units','centimeters','Position',rax.Position+[0,-0.5,0,-rax.Position(4)+0.5]);
    set(histogram(allrsegs(:),linspace([xlimits,100])),'EdgeColor','none');
    sax(end).YTickLabel = [];
    xlim(sax(end),xlimits)
    xlabel(sax(end),'Observed');
    sax(end).Position(3) = rax.Position(3)./rax.PlotBoxAspectRatio(2);
    sax(end).Position(1) = sax(end).Position(1)+(rax.Position(3)-rax.Position(3)./rax.PlotBoxAspectRatio(2))/2;
    sax(end).YTick = [];
    grid(sax(end),'on');
% ADDPLOT : marginal of the decoded pitch during state transition
sax(end+1) = axes('Units','centimeters','Position',rax.Position ./ [1,1,1,1] + [0,0,-rax.Position(3)+0.5,0]);
    sax(end).Position(1) = sax(end-1).Position(1)-0.5;
    set(histogram(allesegs(:),linspace([ylimits,100]),'Orientation','horizontal'),'EdgeColor','none');
    sax(end).XTickLabel = [];
    sax(end).XTick = [];
    if stst==2 | stst==3 | stst==5
        sax(end).YTickLabel = [];
    end
    ylim(sax(end),ylimits)
    if stst == 1
        ylabel(sax(end),{'Body Pitch','Decoded'});
    elseif stst == 4
        ylabel(sax(end),{'Head Pitch','Decoded'});
    end
    rax.Position(1) = rax.Position(1);%+0.175;
    lax.Position(1) = rax.Position(1);
    lax.PlotBoxAspectRatio = rax.PlotBoxAspectRatio;
    lax.OuterPosition = rax.OuterPosition;
    grid(sax(end),'on');
    title(rax,stsTransitionLabels{stst})
% ADDPLOT : bounded line plot of decoded and observed transitions
sax(end+1) = axes('Units','centimeters','Position',rax.Position+[0,-3,0,-rax.Position(4)+1.5]);
    sax(end).Position(3) = rax.Position(3)./rax.PlotBoxAspectRatio(2);
    sax(end).Position(1) = sax(end).Position(1)+(rax.Position(3)-rax.Position(3)./rax.PlotBoxAspectRatio(2))/2;
    grid(sax(end),'on');
    
    if pitchType==1,
        cc = 'kb';
        ylim(sax(end),[-1.5,0.3]);            
    else
        cc='kr';
        ylim(sax(end),[-0.2,1.5]);            
    end
    boundedline(tbins,median(allrsegs,2,'omitnan'),std(allrsegs,[],2,'omitnan'),cc(1),'alpha');
    boundedline(tbins,median(allesegs,2,'omitnan'),std(allesegs,[],2,'omitnan'),cc(2),'alpha');        
    grid(sax(end),'on');
    if stst==1 
        if pitchType == 1
            ylabel(sax(end),'Head-Body Pitch');
        else
            ylabel(sax(end),'Body Pitch');        
        end
    end
    if stst==numel(stsTransitions)
        lh = legend({'','Observed','','Decoded'},'Location','NorthEastOutside');
        lh.Units = 'centimeters';
        lh.Position(1) = sum(rax.Position([1,3]))+0.1;
    end
    if stst == 1
        xlabel(sax(end),'Seconds');
    end
    
end
cax = gobjects([1,0]);
cax(end+1) = colorbar(rax,'EastOutside');
cax(end).Units = 'centimeters';
cax(end+1) = colorbar(lax,'EastOutside');
cax(end).Units = 'centimeters';
cax(1).Position(1) = cax(1).Position(1) + 0.1;
cax(2).Position(1) = cax(1).Position(1) + 0.1;
cax(1).TickLabels = [];
cax(1).Position(1) = sum(rax.Position([1,3]))-(rax.Position(3)-rax.Position(3)./rax.PlotBoxAspectRatio(2))/2+0.1;
cax(2).Position(1) = sum(rax.Position([1,3]))-(rax.Position(3)-rax.Position(3)./rax.PlotBoxAspectRatio(2))/2+0.5;
ylabel(cax(2),'Seconds')



% ADDPLOT : state transition matrix
[yind, yOffSet, xind, xOffSet] = deal(3,-0.5, 7, -1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                         ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,                       ...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,                       ...
                              fig.subplot.width+3,...
                              fig.subplot.height], ...
                  'FontSize', 8,                                                 ...
                  'LineWidth',1);
hold(sax(end),'on');
imagesc(1:numel(sTransitions),1:numel(sTransitions),round(sTransMatrix/sum(sTransMatrix(:)),3));
sax(end).XTick = [1:5];
sax(end).YTick = [1:5];
sax(end).XTickLabel = sTransitions;
sax(end).YTickLabel = sTransitions;
axis(sax(end),'ij');
colormap(sax(end),'jet');
cax = colorbar(sax(end));
cax.Units = 'centimeters';
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
drawnow();
cax.Position(1) = sum(sax(end).Position([1,3]))+0.1;
xlim(sax(end),[0.5,5.5]);
ylim(sax(end),[0.5,5.5]);



savefig(hfig,'/storage/share/Projects/BehaviorPlaceCode/decoding/decoding_pos_bhv_final.fig');



% SUPFIG : All state transitions observed vs decoded state space  ----------------------------------
[hfig,fig,fax,sax] = set_figure_layout(figure(666008),'A4','landscape',[],1.25,1.25,0.1,0.1);
%hfig = figure();
set(hfig,'Units','Centimeters');
sax = gobjects([1,0]);
for stst = 1:numel(stsTransitions)
for pitchType = 1:2
%t = 1; 
allrsegs = [];
allesegs = [];
for t = 1:12
    allrsegs = cat(2,allrsegs,rsegs{t,stst}{pitchType});
    allesegs = cat(2,allesegs,esegs{t,stst}{pitchType});
end
if pitchType == 2,
    xlimits = [-0.25,1.25];
    ylimits = [-0.25,1.25];
else
    xlimits = prctile(allrsegs(nniz(allrsegs(:))),[5,95])+[-0.15,0.15];
    ylimits = prctile(allesegs(nniz(allesegs(:))),[5,95])+[-0.15,0.15];
end

% ADDPLOT : invisible plot which overlays onto the data scatterplot
sax(end+1) = axes('Units','Centimeters','Position',[3+4*(stst-1),3+4*(pitchType-1),3,3]);
    rax = sax(end);
    hold('on');
    for seg = 1:2:size(allrsegs,2),
        scatter(allrsegs(1:2:end,seg),...
                allesegs(1:2:end,seg),...
                1,...
                linspace(-1,1,ceil(size(allrsegs,1)./2)),...
                'Filled');
    end
    xlim(sax(end),xlimits);
    ylim(sax(end),ylimits);
    sax(end).XTickLabel = [];
    sax(end).YTickLabel = [];
    colormap(sax(end),'cool');
    daspect(sax(end),[1,1,1]);
    set(line(rax,[-2,2],[-2,2]),'Color','k','LineWidth',2);
    grid(rax,'on');
    box(rax,'on');
    % SERIOUSLY **** matlab
    drawnow(); pause(1);
    % SERIOUSLY 
% ADDPLOT : invisible plot which overlays onto the data scatterplot
sax(end+1) = axes('Units','centimeters','Position',rax.Position,'Visible','off');
    lax = sax(end);
    hold('on');
    scatter(median(allrsegs(1:2:end,:),2,'omitnan'),...
            median(allesegs(1:2:end,:),2,'omitnan'),...
            15,...
            linspace(-1,1,round(size(allrsegs,1)/2)),...
            'Filled','MarkerEdgeColor','k');
    plot(lax,median(allrsegs(1:2:end,:),2,'omitnan'),...
         median(allesegs(1:2:end,:),2,'omitnan')+std(allesegs(1:2:end,:),[],2,'omitnan'),...
         '-k',...
         'LineWidth',2);
    plot(lax,median(allrsegs(1:2:end,:),2,'omitnan'),...
         median(allesegs(1:2:end,:),2,'omitnan')-std(allesegs(1:2:end,:),[],2,'omitnan'),...
         '-k',...
         'LineWidth',2);
% $$$ plot(lax,median(allrsegs,2,'omitnan')+std(allrsegs,[],2,'omitnan'),...
% $$$      median(allesegs,2,'omitnan'),...
% $$$      '-k',...
% $$$      'LineWidth',2);
% $$$ plot(lax,median(allrsegs,2,'omitnan')-std(allrsegs,[],2,'omitnan'),...
% $$$      median(allesegs,2,'omitnan'),...
% $$$      '-k',...
% $$$      'LineWidth',2);
    set(sax(end),'Visible','off');
    colormap(lax,'jet');
    xlim(sax(end),xlimits);
    ylim(sax(end),ylimits);
    daspect(sax(end),[1,1,1]);
    % SERIOUSLY **** matlab
    drawnow(); pause(1);
    % SERIOUSLY 
% ADDPLOT : marginal of the observed pitch during state transition
sax(end+1) = axes('Units','centimeters','Position',rax.Position+[0,-0.5,0,-rax.Position(4)+0.5]);
    set(histogram(allrsegs(:),linspace([xlimits,100])),'EdgeColor','none');
    sax(end).YTickLabel = [];
    xlim(sax(end),xlimits)
    if pitchType == 1,
        xlabel(sax(end),'Observed');
    end
    sax(end).Position(3) = rax.Position(3)./rax.PlotBoxAspectRatio(2);
    sax(end).Position(1) = sax(end).Position(1)+(rax.Position(3)-rax.Position(3)./rax.PlotBoxAspectRatio(2))/2;
    sax(end).YTick = [];
    grid(sax(end),'on');
% ADDPLOT : marginal of the decoded pitch during state transition
sax(end+1) = axes('Units','centimeters','Position',rax.Position ./ [1,1,1,1] + [0,0,-rax.Position(3)+0.5,0]);
    sax(end).Position(1) = sax(end-1).Position(1)-0.5;
    set(histogram(allesegs(:),linspace([ylimits,100]),'Orientation','horizontal'),'EdgeColor','none');
    sax(end).XTickLabel = [];
    sax(end).XTick = [];
    ylim(sax(end),ylimits)
    if stst == 1 && pitchType == 1
        ylabel(sax(end),{'Head Pitch','Decoded'});
    elseif stst == 1 && pitchType == 2
        ylabel(sax(end),{'Body Pitch','Decoded'});
    end
    rax.Position(1) = rax.Position(1);%+0.175;
    lax.Position(1) = rax.Position(1);
    lax.PlotBoxAspectRatio = rax.PlotBoxAspectRatio;
    lax.OuterPosition = rax.OuterPosition;
    grid(sax(end),'on');
    if pitchType == 2
        title(rax,stsTransitionLabels{stst})
    end
    
end
end
cax = gobjects([1,0]);
cax(end+1) = colorbar(rax,'EastOutside');
cax(end).Units = 'centimeters';
cax(end+1) = colorbar(lax,'EastOutside');
cax(end).Units = 'centimeters';
cax(1).Position(1) = cax(1).Position(1) + 0.1;
cax(2).Position(1) = cax(1).Position(1) + 0.1;
cax(1).TickLabels = [];
cax(1).Position(1) = sum(rax.Position([1,3]))-(rax.Position(3)-rax.Position(3)./rax.PlotBoxAspectRatio(2))/2+0.1;
cax(2).Position(1) = sum(rax.Position([1,3]))-(rax.Position(3)-rax.Position(3)./rax.PlotBoxAspectRatio(2))/2+0.5;
ylabel(cax(2),'Seconds')


% END SUPFIG ---------------------------------------------------------------------------------------
[hfig,fig,fax,sax] = set_figure_layout(figure(666011),'A4','landscape',[],1.25,1.25,0.1,0.1);
for stst = 1:numel(stsTransitions)
    sax(end+1) = subplot2(ny,nx,1,stst);    
    hold(sax(end),'on');    
    for pitchType = 1:2
        allrsegs = [];
        allesegs = [];
        for t = 1:12
            allrsegs = cat(2,allrsegs,rsegs{t,stst}{pitchType});
            allesegs = cat(2,allesegs,esegs{t,stst}{pitchType});
        end
        if pitchType==1,
            cc = 'b';
        else
            cc='r';
        end
        tbins = ([1:size(allrsegs,1)]-round(size(allrsegs,1))/2)/sampleRate;
        plot(tbins,median(allesegs-allrsegs,2,'omitnan'),cc,'LineWidth',2);
    end
    title(stsTransitionLabels{stst});
    grid(sax(end),'on');
    Lines([],0,'k');
    ylim(sax(end),[-0.3,0.3]);
    ylabel(sax(end),'decoded - observed');
    if stst ==numel(stsTransitions),
        lh = legend('HP','BP');
    end
end

[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','landscape',[],1.25,1.25,0.1,0.1);
ny = 2;
sax = gobjects([0,1]);
for stst = 1:numel(stsTransitions)
    for pitchType = 1:2
        sax(end+1) = subplot2(ny,nx,pitchType,stst);    
        sax(end).Units = 'centimeters';
        hold(sax(end),'on');    
        allrsegs = [];
        allesegs = [];
        for t = 1:12
            allrsegs = cat(2,allrsegs,rsegs{t,stst}{pitchType});
            allesegs = cat(2,allesegs,esegs{t,stst}{pitchType});
        end
        tbins = ([1:size(allrsegs,1)]-round(size(allrsegs,1))/2)/sampleRate;        
        if pitchType==1,
            cc = 'kb';
            ylim(sax(end),[-1.5,0.3]);            
        else
            cc='kr';
            ylim(sax(end),[-0.2,1.5]);            
        end
        boundedline(tbins,median(allrsegs,2,'omitnan'),std(allrsegs,[],2,'omitnan'),cc(1),'alpha');
        boundedline(tbins,median(allesegs,2,'omitnan'),std(allesegs,[],2,'omitnan'),cc(2),'alpha');        
        grid(sax(end),'on');
        if stst==1 
            if pitchType == 1
                ylabel(sax(end),'Head-Body Pitch');
            else
                ylabel(sax(end),'Body Pitch');        
            end
        end
        if pitchType == 1
            title(stsTransitionLabels{stst});        
        end
        if stst==numel(stsTransitions)
            lh = legend({'','Observed','','Decoded'},'Location','NorthEastOutside');
            lh.Units = 'centimeters';
            lh.Position(1) = sum(sax(end).Position([1,3]));
        end
    end    
end


% SUPFIG : jpdf of state transition segments with mean over lay for observed and decoded  head and body pitch

[hfig,fig,fax,sax] = set_figure_layout(figure(666009),'A4','landscape',[],1.25,1.25,0.1,0.1);
sax = gobjects([1,0]);
nx = numel(stsTransitions);
ny = 7;
for stst = 1:numel(stsTransitions)
    for pitchType = 1:2
        allrsegs = [];
        allesegs = [];
        for t = 1:12
            allrsegs = cat(2,allrsegs,rsegs{t,stst}{pitchType});
            allesegs = cat(2,allesegs,esegs{t,stst}{pitchType});
        end
        tbins = ([1:size(allrsegs,1)]-round(size(allrsegs,1))/2)/sampleRate;
        if pitchType == 2,
            pbins = linspace(-0.2,1.25,30);
        else
            pbins = linspace(-1.8,0.5,30);
        end
        ebins = linspace(-1,1,30);
        
    sax(end+1) = subplot2(ny,nx,pitchType*3+pitchType-1,stst);
        sax(end).Units = 'centimeters';    
        out = [];
        for tm = 1:size(allesegs,1)
            out(:,tm) = histcounts(allesegs(tm,:)-allrsegs(tm,:),ebins);
        end
        imagesc(tbins,ebins,out);
        axis('xy');
        hold('on');
        plot(tbins,median(allesegs-allrsegs,2,'omitnan'),'r','LineWidth',2);
        title(stsTransitionLabels{stst});        
        if stst == 1
        if pitchType == 1
            ylabel({'Head Pitch','Error'});
        else
            ylabel({'Body Pitch','Error'});
        end
        end
        sax(end).Position(2) = sax(end).Position(2)-2*double(pitchType==1);        
        
    sax(end+1) = subplot2(ny,nx,pitchType*3-1+pitchType-1,stst);
        sax(end).Units = 'centimeters';
        out = [];
        for tm = 1:size(allrsegs,1)
            out(:,tm) = histcounts(allrsegs(tm,:),pbins);
        end
        imagesc(tbins,pbins,out);
        axis('xy');
        hold('on');
        plot(tbins,median(allrsegs,2,'omitnan'),'k','LineWidth',2);
        plot(tbins,median(allesegs,2,'omitnan'),'r','LineWidth',2);
        if stst == 1        
        if pitchType == 1
            ylabel({'Head Pitch','Observed'});
        else
            ylabel({'Body Pitch','Observed'});
        end
        end
        sax(end).Position(2) = sax(end).Position(2)-2*double(pitchType==1);
        
    sax(end+1) = subplot2(ny,nx,pitchType*3-2+pitchType-1,stst);
        sax(end).Units = 'centimeters';    
        out = [];
        for tm = 1:size(allesegs,1)
            out(:,tm) = histcounts(allesegs(tm,:),pbins);
        end
        imagesc(tbins,pbins,out);
        axis('xy');
        hold('on');
        plot(tbins,median(allrsegs,2,'omitnan'),'k','LineWidth',2);
        plot(tbins,median(allesegs,2,'omitnan'),'r','LineWidth',2);
        title(stsTransitionLabels{stst});        
        if stst == 1
        if pitchType == 1
            ylabel({'Head Pitch','Observed'});
        else
            ylabel({'Body Pitch','Observed'});
        end
        end
        sax(end).Position(2) = sax(end).Position(2)-2*double(pitchType==1);        
        
        
    end
end


% ACCUMULATE transitions

% high to low decoding error head and body pitch
% low to high decoding error head and body pitch
%%%>>>

% ENDFIG





%%%>>>
