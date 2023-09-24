
configure_default_args();
EgoProCode2D_load_data();

% CA1
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


decoded = struct('fwd',[],...
                 'lat',[],...
                 'hvf',[],...
                 'hvl',[],...
                 'hav',[],...
                 'hba',[],...
                 'phz',[]);
for t = [1:3,5:8,11],
    %for t = [1:3,5:8,11],
    mind =    dca{t}.stcm(:,1)==1                                            ...
              & (dca{t}.stcm(:,3)==3|dca{t}.stcm(:,4)==4|dca{t}.stcm(:,5)==5)  ...
              & dca{t}.post>0.005 ...
         ...   & dca{t}.hvfl(:,1)>-2                                             ...
         ...  & abs(dca{t}.hvfl(:,2))>5                                        ...
            & dca{t}.ucnt>=2 & dca{t}.ucnt<8                                 ...
            & sqrt(sum(dca{t}.xyz(:,'hcom',[1,2]).^2,3))<325;
    decoded.fwd = cat(1,decoded.fwd,dca{t}.esax(mind,1));
    decoded.lat = cat(1,decoded.lat,dca{t}.esax(mind,2));%+20*double(t>4));
    decoded.hvf = cat(1,decoded.hvf,dca{t}.hvfl(mind,1));
    decoded.hvl = cat(1,decoded.hvl,dca{t}.hvfl(mind,2));
    decoded.hav = cat(1,decoded.hav,dca{t}.hvang(mind,1));
    decoded.hba = cat(1,decoded.hba,dca{t}.hbang(mind,1));
    decoded.phz = cat(1,decoded.phz,dca{t}.phz(mind,1));
end

ind = WithinRanges(decoded.phz,phzBin.edges(3:4)) ...
      & randn(size(decoded.hba))>0.5 ...
      & abs(decoded.hba)<1.2;
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind)]);
STATS


figure,
ind = WithinRanges(decoded.phz,phzBin.edges([3,4])) ...
      & randn(size(decoded.hba))>0              ...
 ;...     & abs(decoded.hvl)>5;
hist2([decoded.lat(ind),decoded.hba(ind)],linspace(-300,300,24),linspace(-1.2,1.2,24),'xprob')
colormap('jet');
Lines([],0,'w');
Lines(0,[],'w');

figure,
hist2([decoded.lat(ind),decoded.hvl(ind)],linspace(-300,300,24),linspace(-40,40,24),'xprob')
colormap('jet');
Lines([],0,'w');
Lines(0,[],'w');


mBinHvl.edges = linspace(-40,40,4);
mBinHvl.centers = mean([mBinHvl.edges(1:end-1);mBinHvl.edges(2:end)])
mBinHvl.count = numel(mBinHvl.edges)-1;
mBinHba.edges = linspace(-1.2,1.2,4);
mBinHba.centers = mean([mBinHba.edges(1:end-1);mBinHba.edges(2:end)])
mBinHba.count = numel(mBinHba.edges)-1;
mout = zeros([mBinHba.count,mBinHvl.count]);
for a = 1:mBinHba.count
    for v = 1:mBinHvl.count
        ind =   WithinRanges(decoded.phz,phzBin.edges([3,4])) ...
              & WithinRanges(decoded.hvl,mBinHvl.edges([v,v+1])) ...
              & WithinRanges(decoded.hba,mBinHba.edges([a,a+1])) ...              
              & randn(size(decoded.hba))>0;
        mout(a,v) = mean(decoded.lat(ind),'omitnan');
    end
end

figure,imagesc(mBinHba.centers,mBinHvl.centers,mout')
caxis([-60,60])
colormap('jet');



figure,
hist2([decoded.fwd,decoded.hba],linspace(-300,300,24),linspace(-1.2,1.2,24),'xprob')
colormap('jet');
Lines([],0,'w');
Lines(0,[],'w');

figure,
subplot(311);
ind = decoded.hba>0.2;
hist2([decoded.lat(ind),decoded.fwd(ind)],linspace(-400,400,32),linspace(-300,500,32));
caxis([0,4000])
Lines(0,[],'w');
Lines([],0,'w');
subplot(312);
ind = abs(decoded.hba)<0.2;
hist2([decoded.lat(ind),decoded.fwd(ind)],linspace(-400,400,32),linspace(-300,500,32));
caxis([0,4000])
Lines(0,[],'w');
Lines([],0,'w');
subplot(313);
ind = decoded.hba<-0.2;
hist2([decoded.lat(ind),decoded.fwd(ind)],linspace(-400,400,32),linspace(-300,500,32));
caxis([0,4000])
Lines(0,[],'w');
Lines([],0,'w');
colormap('jet');

binPhzs = linspace(0.5,2*pi-0.5,4);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
hbaBinEdg = [-1.2,-0.2,0.2,1.2];
hbaBinCtr = mean([hbaBinEdg(1:end-1);hbaBinEdg(2:end)]);
                          
hbaBin.edges = [-1.2,-0.2,0.2,1.2];
hbaBin.centers = mean([hbaBin.edges(1:end-1);hbaBin.edges(2:end)]);
hbaBin.count = numel(hbaBin.centers);        

phzBin.edges = linspace(0.5,2*pi-0.5,4);
phzBin.centers = (binPhzs(1:end-1)+binPhzs(2:end))./2;
phzBin.count = numel(phzBin.centers);


hbaBnds = {[0.2,1.2],[-0.2,0.2],[-1.2,-0.2]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
havBnds = {[-0.3,-0.018],[-0.018,-0.009],[-0.009,0.009],[0.009,0.018],[0.018,0.3]};
hvlBnds = {[-50,-5],[-5,5],[5,50]};
figure
norm = 'xprob';
xBinEdg = linspace(-300,300,8);
xBinCtr = mean([xBinEdg(1:end-1);xBinEdg(2:end)]);
yBinEdg = linspace(0.5,2*pi-0.5,4);
yBinCtr = mean([yBinEdg(1:end-1);yBinEdg(2:end)]);
clims = [0,0.4];
%clims = 'auto';
out = zeros([7,3,3,3]);
for h = 1:numel(hbaBnds)
    for v = 1:numel(hvlBnds)
        subplot2(numel(hvlBnds),numel(hbaBnds),v,h);
        ind = WithinRanges(decoded.hba,hbaBnds{h}) & ...
              WithinRanges(decoded.hvl,hvlBnds{v}) & ...
              randn(size(decoded.hvl))>0.5;
        out(:,:,h,v) = hist2([decoded.lat(ind), decoded.phz(ind)],xBinEdg,yBinEdg,norm);
        imagesc(xBinCtr, yBinCtr, out(:,:,h,v)');
        colormap('jet');caxis(clims);
        axis('xy');
    end
end


hbaBnds = {[0.2,1.2],[-0.2,0.2],[-1.2,-0.2]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
havBnds = {[-0.3,-0.018],[-0.018,-0.009],[-0.009,0.009],[0.009,0.018],[0.018,0.3]};
hvlBnds = {[-50,-5],[-5,5],[5,50]};
figure
norm = 'xprob';
xBinEdg = linspace(-300,300,8);
xBinCtr = mean([xBinEdg(1:end-1);xBinEdg(2:end)]);
yBinEdg = linspace(0.5,2*pi-0.5,4);
yBinCtr = mean([yBinEdg(1:end-1);yBinEdg(2:end)]);
clims = [0,0.4];
%clims = 'auto';
out = zeros([7,3,3,3]);
for h = 1:numel(hbaBnds)
    for v = 1:numel(hvlBnds)
        subplot2(numel(hvlBnds),numel(hbaBnds),v,h);
        ind = WithinRanges(decoded.hba,hbaBnds{h}) & ...
              WithinRanges(decoded.hvl,hvlBnds{v}) & ...
              randn(size(decoded.hvl))>0.5;
        out(:,:,h,v) = hist2([decoded.fwd(ind), decoded.phz(ind)],xBinEdg,yBinEdg,norm);
        imagesc(xBinCtr, yBinCtr, out(:,:,h,v)');
        colormap('jet');caxis(clims);
        axis('xy');
    end
end

%% CDF
hbaBnds = {[0.2,1.2],[-0.2,0.2],[-1.2,-0.2]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
%havBnds = {[-0.3,-0.015],[-0.015,0.015],[0.015,0.3]};
havBnds = {[-0.3,-0.018],[-0.018,-0.009],[-0.009,0.009],[0.009,0.018],[0.018,0.3]};
hvlBnds = {[-50,-5],[-5,5],[5,50]};
figure
norm = 'xprob';
xBinEdg = linspace(-300,300,8);
xBinCtr = mean([xBinEdg(1:end-1);xBinEdg(2:end)]);
yBinEdg = linspace(0.5,2*pi-0.5,4);
yBinCtr = mean([yBinEdg(1:end-1);yBinEdg(2:end)]);
clims = [0,0.4];
%clims = 'auto';
out = zeros([7,3,3,3]);
for h = 1:numel(hbaBnds)
    for p = 1:3;    
    subplot2(phzBin.count,hbaBin.count,phzBin.count+1-p,h);
    hold('on');
    for v = 1:numel(hvlBnds)
        ind = WithinRanges(decoded.hba,hbaBnds{h}) & ...
              WithinRanges(decoded.hvl,hvlBnds{v}) & ...
              randn(size(decoded.hvl))>0.5;
        pind = WithinRanges(decoded.phz,phzBin.edges(p:p+1));
        cdfplot(decoded.lat(ind&pind));
    end
    xlim([-300,300]);
    end
end



figure,
clims = [0,0.05];
for h = 1:3
    for d = 1:2
        subplot2(2,3,d,h);
        imagesc(xBinCtr,yBinCtr,diff(out(:,:,h,[0:1]+d),[],4)');
        colormap('jet'); caxis(clims); axis('xy');
    end
end


figure,
clims = [-0.1,0.1];
for v = 1:3
    subplot2(3,2,v,1);
    imagesc(xBinCtr,yBinCtr,(out(:,:,1,v)-out(:,:,2,v))');
    colormap('jet'); caxis(clims); axis('xy');
    subplot2(3,2,v,2);
    imagesc(xBinCtr,yBinCtr,(out(:,:,3,v)-out(:,:,2,v))');
    colormap('jet'); caxis(clims); axis('xy');
end


figure,
clims = [-0.1,0.1];
for v = 1:3
    subplot2(2,3,1,v);
    imagesc(xBinCtr,yBinCtr,(out(:,:,v,1)-out(:,:,v,2))');
    colormap('jet'); caxis(clims); axis('xy');
    subplot2(2,3,2,v);
    imagesc(xBinCtr,yBinCtr,(out(:,:,v,3)-out(:,:,v,2))');
    colormap('jet'); caxis(clims); axis('xy');
end

figure
norm = 'xprob';
xBinEdg = linspace(-250,250,8);
yBinEdg = linspace(0.5,2*pi-0.5,4);
clims = [0,0.4];
for v = 1:3
    subplot2(1,3,1,v);
    ind = WithinRanges(decoded.hav,havBnds{v}) & ...
          randn(size(decoded.hba))>0.5;
    hist2([decoded.lat(ind), decoded.phz(ind)],xBinEdg,yBinEdg,norm);
    colormap('jet');caxis(clims);
end



figure
norm = 'xprob';
clims = [0.01,0.4];
sax = tight_subplot(3,1,[0.01,0.01],[0.1,0.1],[0.1,0.1]);
for v = 1:3
    ind = WithinRanges(decoded.hba,hbaBnds{v}) & ...
          randn(size(decoded.hba))>0.5;
    out = hist2([decoded.lat(ind)+25, decoded.phz(ind)],xBinEdg, ...
                yBinEdg,norm);
    for p = 1:3
        axes(sax(4-p))
        hold(sax(p),'on');
        plot(xBinCtr,cumsum(out(:,p))','-+','Color',bclr(v,:));
        Lines([],0.5,'k');
        Lines(0,[],'k');
    end
    ylim([0,1]);
    xlim([-250,250]);
    %colormap('jet');caxis(clims);
    %set(gca,'ColorScale','log');    
end




ind = WithinRanges(decoded.phz,binPhzs(3:4)) & randn(size(decoded.hba))>0.5;
ind = WithinRanges(decoded.phz,binPhzs(1:2)) & randn(size(decoded.hba))>0.5;
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind),decoded.hav(ind)]);
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind),decoded.hvl(ind)]);
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hav(ind)]);

ind = WithinRanges(decoded.phz,phzBin.edges(3:4)) ...
      & randn(size(decoded.hba))>0.5 ...
      & abs(decoded.hba)<1.2;
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind)]);
STATS

glm = fitglm([decoded.hba(ind),decoded.hav(ind)],decoded.lat(ind));
glm = fitglm([decoded.hba(ind),decoded.hvl(ind)],decoded.lat(ind));

% MTAData
% function summarize



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
phaseLabels = {'descending','trough','ascending'};
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
    title(phaseLabels{p});
    %sax(end).XTick = [-20,-10,0,10,20];
    if p~=1,
        %sax(end).XTickLabel = [];
        %sax(end).YTickLabel = [];
        sax.XTickLabel = [];
        %sax.YTickLabel = [];
        xlabel('')
        ylabel('')
    end
end
                 

