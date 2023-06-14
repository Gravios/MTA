% req20220310
%     Tags:  decoding egocentric hba hvfl
%     Status: Active
%     Type: Analysis
%     Author: Justin Graboski
%     Final_Forms: NA
%     Project: General
%     Description: computing the stats on 


configure_default_args();

MjgER2016_load_data();

tind = [3,4,5,17,18,19,20,21,22,23,29];
% $$$ tind = [6,7,26,27,30];
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



% $$$ roll = cf(@(T) fet_roll(T,sampleRate), Trials(tind));

%% HEAD BODY ANGLE ---------------------------------------------------------------------------------
clear('xcomp','ycomp','zcomp','ccomp');
xcomp.label = 'hba';            xcomp.data = [];    xcomp.edgs = [-1.2,-0.2,0.2,1.2];
ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.28, 2*pi-0.5,6 );
ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-80,80];
%ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-1,1];
acomp.label = 'hvang';          acomp.data = [];   
hcomp.label = 'hvang';          hcomp.data = [];   
fcomp.label = 'egofwd';         fcomp.data = [];   
rcomp.label = 'hvang';          rcomp.data = [];   
%for t = 5:10
for t = 1:10
    dc = dca{t};
    mind =  dc.stcm(:,1)==1                                      ...
           & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5)  ...               
           & dc.hvfl(:,1)>2 ...
           & dc.post > 0.032 ...
           & dc.ucnt>=2 & dc.ucnt<10 ...
           & sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))<325;
    mind(mind==true) = randn([sum(mind),1])>0.75;
    xcomp.data = cat(1, xcomp.data, -dc.hbang(mind,1));
    ycomp.data = cat(1, ycomp.data, dc.phz(mind));
    ccomp.data = cat(1, ccomp.data, dc.esax(mind,2));    % center corrections 
    %ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)+25*double(t<4)-12.5*double(t>=4));
    fcomp.data = cat(1, fcomp.data, dc.esax(mind,1)-25);
    acomp.data = cat(1,acomp.data,dc.hvang.data(mind));
    hcomp.data = cat(1,hcomp.data,dc.hvfl(mind,2));    
end

% SELECT states: theta and ( walk or turn or pause )
% $$$ ind = dc.stcm(:,1)==1 & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5);
% $$$ figure,
% $$$ out = hist2([dc.hvfl(ind,2),dc.hvfl(ind,1)],linspace(-60,60,30),linspace(-40,100,50));
% $$$ set(gca(),'ColorScale','log');
% $$$ imagesc(log10(out'))
% $$$ axis('xy');
% $$$ %caxis([0,1000])



[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,ccomp);
zmean   = zcomp.mean;    zmean(zcomp.count<5)  = nan;
zmedian = zcomp.median;  zmean(zcomp.count<5)  = nan;
zstd    = zcomp.std;     zstd(zcomp.count<5)   = nan;
zcount  = zcomp.count;   zcount(zcomp.count<5) = nan;
set(figure(),'Units','centimeters','Position',[0,-3,8,28]);
subplot(411); imagesc(xcomp.ctrs,ycomp.ctrs,zmean'); axis('xy');
    cax = colorbar(); ylabel(cax,['Mean ',ccomp.label]); colormap('jet');  caxis([ccomp.clim]);
    ylabel(ycomp.label);    
subplot(412); imagesc(xcomp.ctrs,ycomp.ctrs,zstd'); axis('xy');
    cax = colorbar(); ylabel(cax,['Std ',ccomp.label]); 
    colormap('jet');
subplot(413); 
imagesc(xcomp.ctrs,ycomp.ctrs,zcount'); axis('xy');
    cax = colorbar(); ylabel(cax,'Count'); colormap('jet'); 
    xlabel(xcomp.label);
subplot(414); 
    imagesc(xcomp.ctrs,ycomp.ctrs,bsxfun(@rdivide,zcount,sum(zcount,2,'omitnan'))'); axis('xy');    
    cax = colorbar(); ylabel(cax,'Probability'); colormap('jet'); 
    xlabel(xcomp.label);


figure,
pinds = {ycomp.data > 0.28  & ycomp.data < 2.18,...
         ycomp.data > 2.18  & ycomp.data < 4,...
         ycomp.data > 4  & ycomp.data < 6};
for p = 1:numel(pinds)
subplot2(3,3,2,p);
    hold('on');
    indL =    pinds{p} & xcomp.data > -1.2 & xcomp.data < -0.2;
    indC =    pinds{p} & xcomp.data > -0.2 & xcomp.data < 0.2;
    indR =    pinds{p} & xcomp.data >  0.2 & xcomp.data < 1.2;    
    
    histogram(ccomp.data(indL),linspace(-400,400,50),'Normalization','probability','FaceColor','g');
    histogram(ccomp.data(indC),linspace(-400,400,50),'Normalization','probability','FaceColor','b');
    histogram(ccomp.data(indR),linspace(-400,400,50),'Normalization','probability','FaceColor','r');
    
    xlim([-300,300]);
    ylim([0,0.1]);
    
    Lines(mean(ccomp.data(indL)),[],'g');    
    Lines(mean(ccomp.data(indC)),[],'b');
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

% PERMUTATION TEST
distPerms = zeros([1000,3]);
for p = 1:numel(pinds)
    indL =    pinds{p} & xcomp.data > -1.2 & xcomp.data < -0.2;
    indC =    pinds{p} & xcomp.data > -0.2 & xcomp.data < 0.2;
    indR =    pinds{p} & xcomp.data >  0.2 & xcomp.data < 1.2;
    
    sinds =  [find(indL);find(indR)];
    hsCnt = floor(numel(sinds)/2);

    for n = 1:10000
        sinds = sinds(randperm(numel(sinds)));
        distPerms(n,p) = mean(ccomp.data(sinds(1:hsCnt)))-mean(ccomp.data(sinds(hsCnt:end)));
    end
end    

figure,
hold('on');
histogram(distPerms(:,1)/10,linspace(-4,4,50),'Normalization','probability','FaceColor','g');
histogram(distPerms(:,2)/10,linspace(-4,4,50),'Normalization','probability','FaceColor','b');
histogram(distPerms(:,3)/10,linspace(-4,4,50),'Normalization','probability','FaceColor','r');
Lines(mean(ccomp.data(indL))/10,[],'g');
Lines(mean(ccomp.data(indC))/10,[],'b');
Lines(mean(ccomp.data(indR))/10,[],'r');
xlim([-6,6]);    
xlabel('cm');
ylabel('normalized count');
title({'Permutation Test','Distribution of Permuted Differences',...
       'between leftward, centered, and rightward',...
       'head-body postures n = 10000'});




% sessions with retrospective coding
% ec16-00267278 SOF
% ec16-00228240 LIN

subplot(412);
    hold('on');
    histogram(ccomp.data(indL),linspace(-400,400,50),'Normalization','probability');
    histogram(ccomp.data(indC),linspace(-400,400,50),'Normalization','probability');
    histogram(ccomp.data(indR),linspace(-400,400,50),'Normalization','probability');
    Lines(mean(ccomp.data(indC)),[],'g');
    Lines(mean(ccomp.data(indR)),[],'r');
    Lines(mean(ccomp.data(indL)),[],'b');
    xlim([-300,300]);
    
    indA = indR|indC|indL;
    [B,BINT,R,RINT,RSTATS] = regress(ccomp.data(indA),[ones([sum(indA),1]),xcomp.data(indA)]);
%[B,BINT,R,RINT,RSTATS] = regress(ccomp.data(indA),[ones([sum(indA),1]),xcomp.data(indA),acomp.data(indA),hcomp.data(indA)]);
    
% RSTATS R-square statistic,  F statistic and p value, error variance
subplot(411);
    P = polyfit(xcomp.data(indA),ccomp.data(indA),1);
    hold('on');
    hist2([ccomp.data(indA),xcomp.data(indA)],linspace(-300,300,30),linspace(-1.25,1.25,12),'xprob');
    axis('tight');
    line(polyval(P,[-1.25,1.25]),[-1.25,1.25],'Color','m');
    % $$$ figure,hist2([ccomp.data(indA),ycomp.data(indA)],linspace(-pi,pi,11),linspace(-1,1,11),'xprob');
    % $$$ figure,plot(ccomp.data(indA),ycomp.data(indA),'.');
    [h,p,ci,tstats] = ttest2(ccomp.data(indR),ccomp.data(indC))
    [h,p,ci,tstats] = ttest2(ccomp.data(indL),ccomp.data(indC))
    [h,p,ci,tstats] = ttest2(ccomp.data(indL),ccomp.data(indR))
subplot(413);
    hold('on');
    cdfplot(ccomp.data(indL))
    cdfplot(ccomp.data(indC))
    cdfplot(ccomp.data(indR))
    xlim([-300,300]);
subplot(414);
    hold('on');
    cdfplot(ccomp.data(indA&randn([size(indA),1])>0.5))
    cdfplot(ccomp.data(indA&randn([size(indA),1])>0.5))    
    cdfplot(ccomp.data(indA&randn([size(indA),1])>0.5))    
    xlim([-300,300]);
RSTATS    

figure,
pind = ycomp.data > 4.25   & ycomp.data < 5.75;
indR = pind  &  xcomp.data > 0.2 & xcomp.data < 1.25;
indC = pind  &  xcomp.data > -0.2 & xcomp.data < 0.2;
indL = pind  &  xcomp.data < -0.2 & xcomp.data > -1.25;
subplot(311); 
    out = hist2([ccomp.data(indL),fcomp.data(indL)],linspace(-300,300,50),linspace(-300,300,50));
    imagesc(linspace(-300,300,30),linspace(-300,300,30),imgaussfilt(out',1.2));axis('xy');
    Lines([],0,'w');
    Lines(0,[],'w');
    set(gca,'ColorScale','log');
    caxis([1,26]);    
subplot(312); 
    out = hist2([ccomp.data(indC),fcomp.data(indC)],linspace(-300,300,50),linspace(-300,300,50));
    imagesc(linspace(-300,300,30),linspace(-300,300,30),imgaussfilt(out',1.2));axis('xy');
    Lines([],0,'w');
    Lines(0,[],'w');
    set(gca,'ColorScale','log');
    caxis([1,26]);    
subplot(313); 
    out = hist2([ccomp.data(indR),fcomp.data(indR)],linspace(-300,300,50),linspace(-300,300,50));
    imagesc(linspace(-300,300,30),linspace(-300,300,30),imgaussfilt(out',1.2));axis('xy');
    Lines([],0,'w');
    Lines(0,[],'w');
    set(gca,'ColorScale','log');
    caxis([1,26]);
colormap('jet');





figure,
hold('on');
indL =    xcomp.data > 1   & xcomp.data < 2.5 ...
       & ycomp.data > 0.2 & ycomp.data < 0.8;
histogram(ccomp.data(indL),linspace(-400,400,50),'Normalization','probability');
indC =    xcomp.data > 1  & xcomp.data < 2.5 ...
       & ycomp.data > -0.2 & ycomp.data < 0.2;
histogram(ccomp.data(indC),linspace(-400,400,50),'Normalization','probability');
indR =    xcomp.data > 2   & xcomp.data < 2.5 ...
       & ycomp.data < -0.2 & ycomp.data > -0.8;
histogram(ccomp.data(indR),linspace(-400,400,50),'Normalization','probability');

[h,p,ci,tstats] = ttest2(ccomp.data(indR),ccomp.data(indC))
[h,p,ci,tstats] = ttest2(ccomp.data(indL),ccomp.data(indC))
[h,p,ci,tstats] = ttest2(ccomp.data(indL),ccomp.data(indR))

figure,
hold('on');
cdfplot(ccomp.data(indL))
cdfplot(ccomp.data(indC))
cdfplot(ccomp.data(indR))
xlim([-300,300]);



%% HEAD VELOCITY FORWARD ---------------------------------------------------------------------------
clear('xcomp','ycomp','zcomp','ccomp');
xcomp.label = 'hvf';            xcomp.data = [];    xcomp.edgs = [-25,-5,5,25,80]; 
ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5, 4);
ccomp.label = 'ego fwd';        ccomp.data = [];    ccomp.clim = [-50,50];
%for t = 4:10%:numel(tind)-1
%for t = 4:10%:numel(tind)-1    
for t = 1:4%:numel(tind)-1        
    dc = dca{t};
    mind = dc.stcm(:,1)==1                                      ...
           & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5)  ...
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



t = 7;
ind = find(dca{t}.stcm(:,1)==1,1,'first');
ind = [1:2000]+33000;
ind = ind(dca{t}.ucnt(ind)>1&dca{t}.stcm(ind,1)==1);
figure,
hold('on');
plot(dca{t}.sax(ind,1),...
     dca{t}.sax(ind,2),'.')
plot(dca{t}.com(ind,1),...
     dca{t}.com(ind,2),'.g')
plot(dca{t}.xyz(ind,'hcom',1),...
     dca{t}.xyz(ind,'hcom',2),'.r')
plot(dca{t}.xyz(ind(1),'hcom',1),...
     dca{t}.xyz(ind(1),'hcom',2),'*m')
xlim([-450,450]);
ylim([-450,450]);
circle(0,0,400);
circle(0,0,325);

[h,p,ci,tstats] = ttest2(ccomp.data(indR),ccomp.data(indC))
[h,p,ci,tstats] = ttest2(ccomp.data(indL),ccomp.data(indC))
[h,p,ci,tstats] = ttest2(ccomp.data(indL),ccomp.data(indR))

[B,BINT,R,RINT,RSTATS] = regress(ccomp.data(indA),[ones([sum(indA),1]),ycomp.data(indA)]);
[B,BINT,R,RINT,RSTATS] = regress(ccomp.data(indA),[ones([sum(indA),1]),ycomp.data(indA),acomp.data(indA)]);
[B,BINT,R,RINT,RSTATS] = regress(ccomp.data(indA),[ones([sum(indA),1]),ycomp.data(indA),acomp.data(indA),hcomp.data(indA)]);


figure,hist2([ccomp.data(indA),ycomp.data(indA)],linspace(-300,300,10),linspace(-10,80,10),'xprob');




%% HEAD VELOCITY LATERAL ---------------------------------------------------------------------------
clear('xcomp','ycomp','zcomp','ccomp');
xcomp.label = 'hvl';            xcomp.data = [];    xcomp.edgs = [-60,-10,-2,2,10,60]; 
ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5, 4);
ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-50,50];
for t = 1:4%:numel(tind)-1
    dc = dca{t};
    mind = dc.stcm(:,1)==1                                      ...
           & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5)  ...
           & dc.post>0.005                                       ...
           & dc.ucnt>=4                                        ...
           & sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))<325;
    mind(mind==true) = randn([sum(mind),1])>0.25;
    xcomp.data = cat(1,xcomp.data,dc.hvfl(mind,2));    
    ycomp.data = cat(1,ycomp.data,dc.phz(mind));
    ccomp.data = cat(1,ccomp.data,dc.esax(mind,2));    
    %ccomp.data = cat(1,ccomp.data,dc.esax(mind,2)-12.5*double(ismember(t,4:10))+25*double(ismember(t,1:3)));
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
pinds = ycomp.data > 4  & ycomp.data < 6;
%pinds = ycomp.data > 0.5  & ycomp.data < 2.26;
subplot(312);
    hold('on');
    indL =    pinds & -60 < xcomp.data & xcomp.data < -5;
    indC =    pinds &  -5 < xcomp.data & xcomp.data <  5;
    indR =    pinds &   5 < xcomp.data & xcomp.data < 60;    

    histogram(ccomp.data(indL),linspace(-400,400,50),'Normalization','probability');
    histogram(ccomp.data(indC),linspace(-400,400,50),'Normalization','probability');
    histogram(ccomp.data(indR),linspace(-400,400,50),'Normalization','probability');
    
    Lines(mean(ccomp.data(indL)),[],'b');    
    Lines(mean(ccomp.data(indC)),[],'g');
    Lines(mean(ccomp.data(indR)),[],'r');
    xlim([-300,300]);
    
    indA = indR|indC|indL;
    [B,BINT,R,RINT,RSTATS] = regress(ccomp.data(indA),[ones([sum(indA),1]),xcomp.data(indA)]);
    
% RSTATS R-square statistic,  F statistic and p value, error variance
subplot(311);
    P = polyfit(xcomp.data(indA),ccomp.data(indA),1);
    hold('on');
    hist2([ccomp.data(indA),xcomp.data(indA)],linspace(-300,300,30),linspace(-60,60,8),'xprob');
    line(polyval(P,[-60,60]),[-60,60],'Color','m');
    % $$$ figure,hist2([ccomp.data(indA),ycomp.data(indA)],linspace(-pi,pi,11),linspace(-1,1,11),'xprob');
    % $$$ figure,plot(ccomp.data(indA),ycomp.data(indA),'.');
    [h,p,ci,tstats] = ttest2(ccomp.data(indR),ccomp.data(indC))
    [h,p,ci,tstats] = ttest2(ccomp.data(indL),ccomp.data(indC))
    [h,p,ci,tstats] = ttest2(ccomp.data(indL),ccomp.data(indR))
subplot(313);
    hold('on');
    cdfplot(ccomp.data(indL))
    cdfplot(ccomp.data(indC))
    cdfplot(ccomp.data(indR))
    xlim([-300,300]);
RSTATS    





% TEMPORAL shift analysis
% find point in trajectory closest to decoded point
% report time shift and distance

shiftvals = -round(sampleRate*2):round(sampleRate/50):round(sampleRate*2);
% SET time bins
slidingTimeWindowSize = round(sampleRate*1);
timevals = 1:size(xyz,1);
timevals(~nniz(xyz)) = [];
timevals = timevals(1:slidingTimeWindowSize:numel(timevals));
timeBinInds  = discretize([1:size(xyz,1)]',timevals);
badTimeBinInds = find(diff(timevals)>slidingTimeWindowSize)+1;
timeBinInds(ismember(timeBinInds,badTimeBinInds)) = 0;
% SET phase bins
phasevals = linspace(-pi,pi,13);
phaseBinInds = discretize(phz.data,phasevals);
% SET speed bins
speedvals = linspace(0,1.8,7);
speedBinInds = discretize(vxy(:,2),speedvals);

eposTPRErrorXY = nan([numel(timevals)-1,numel(shiftvals),numel(phasevals)-1,numel(speedvals)-1]);
%eposTPRErrorHP = nan([numel(timevals)-1,numel(shiftvals),numel(phasevals)-1,numel(speedvals)-1]);
%eposTPRErrorBP = nan([numel(timevals)-1,numel(shiftvals),numel(phasevals)-1,numel(speedvals)-1]);

eposTPRErrorXYTrj = nan([numel(timevals)-1,numel(shiftvals),numel(phasevals)-1,numel(speedvals)-1,2]);

meanTPRSpeed = nan([numel(timevals)-1,numel(phasevals)-1,numel(speedvals)-1]);
meanTPRPostMax = nan([numel(timevals)-1,numel(phasevals)-1,numel(speedvals)-1]);
meanTPRPhase = nan([numel(timevals)-1,numel(phasevals)-1,numel(speedvals)-1]);

for time = 1:numel(timevals)-2,
    %tic
    disp([num2str(time),' of ',num2str(numel(timevals)-3)])

    timeWindow = time-2:time+2;
    
    timeWindowIndex = ismember(timeBinInds,timeWindow);

    txyz = sq(xyz(timeWindowIndex,'nose',[1,2]));
    tvxy = vxy(timeWindowIndex,2);
    tbvec = bvec(timeWindowIndex,:,:);    
    tposteriorMax = posteriorMax(timeWindowIndex);
    tphz = phz(timeWindowIndex);    
    
    tind =   unitInclusion(timeWindowIndex)>=3  ...
             & stcm(timeWindowIndex,1)==1       ...
             & stcm(timeWindowIndex,2)~=2       ...
             & timeBinInds((timeWindowIndex))==time...
             & tposteriorMax>0.002;

    tempPhaseBinInds = phaseBinInds(timeWindowIndex);
    tempSpeedBinInds = speedBinInds(timeWindowIndex);

    eposThetaPhaseRes = nan([size(tind,1),2]);
    eposThetaPhaseRes(timeBinInds(timeWindowIndex)==time,:) = posEstCom(timeBinInds==time,1:2);
    %eposThetaPhaseRes(timeBinInds(timeWindowIndex)==time,:) = posEstMax(timeBinInds==time,1:2);
% $$$     eposThetaPhaseRes = nan([size(tind,1),size(posEstCom,2)]);
% $$$     eposThetaPhaseRes(timeBinInds(timeWindowIndex)==time,:) = posEstCom(timeBinInds==time,1:2);

    for speed = 1:numel(speedvals)-1,   
        for phase = 1:numel(phasevals)-1,   
            ind = tind  &  tempPhaseBinInds==phase  &  tempSpeedBinInds==speed;
            
            if sum(ind)>3
                meanTPRSpeed(time,phase,speed) = mean(tvxy(ind));
                meanTPRPostMax(time,phase,speed) = mean(tposteriorMax(ind));
                meanTPRPhase(time,phase,speed) = circ_mean(tphz(ind));
                tepos = eposThetaPhaseRes;
                tepos(~ind,:) = nan;


                                    
                for shift = 1:numel(shiftvals)
                    eposTPRErrorXYTrj(time,shift,phase,speed,:) = ...
                       mean(multiprod(circshift(tepos(:,1:2),shiftvals(shift))-txyz,tbvec,2,[2,3]),'omitnan');
                    
                    eposTPRErrorXY(time,shift,phase,speed) = ...
                        mean(sqrt(sum((txyz-circshift(tepos,shiftvals(shift))).^2,2)),'omitnan');

% $$$                     %                    eposTPRErrorHP(time,shift,phase,speed) = mean(sqrt(sum((fet(:,1) ...
% $$$                                        -circshift(eposThetaPhaseRes(:,3),shiftvals(shift))).^2,2)),'omitnan');
% $$$                     %eposTPRErrorBP(time,shift,phase,speed) = mean(sqrt(sum((fet(:,2)...
% $$$                                        -circshift(eposThetaPhaseRes(:,4),shiftvals(shift))).^2,2)),'omitnan');
                end
            end
        end
    end
    %toc
end



% $$$ figure();
% $$$ for speed = 1:numel(speedvals)-1,
% $$$ subplot2(2,numel(speedvals)-1,1,speed);    
% $$$     imagesc(phasevals,shiftvals/sampleRate,eposTPRErrorXY(:,:,speed));axis('xy');
% $$$ caxis([100,300])
% $$$ subplot2(2,numel(speedvals)-1,2,speed);
% $$$     imagesc(phasevals,shiftvals/sampleRate0,1./(bsxfun(@rdivide,eposTPRErrorXY(:,:,speed),min( ...
% $$$         eposTPRErrorXY(:,:,speed)))+eps));axis('xy');
% $$$     title(num2str(mean(10.^speedvals(speed:speed+1))))
% $$$ end



derror = eposTPRErrorXY;
%derror = abs(eposTPRErrorXYTrj(:,1));
% $$$ derror = eposTPRErrorBP;
% $$$ derror = eposTPRErrorHP;

mind = [];
minv = [];
%for time = 1:numel(timevals)-1,
for t = 1:time
    for speed = 1:numel(speedvals)-1,
        %[~,mind(t,speed,:)] = max(RectFilter(RectFilter(1./(bsxfun(@rdivide,sq(derror(t,:,:,speed)),min(sq(derror(t,:,:,speed))))+eps)',3,1)',3,1));
        [minv(t,speed,:),mind(t,speed,:)] = min(RectFilter(RectFilter(sq(derror(t,:,:,speed))',3,1)',3,1));
    end
end
mind(mind==0) = nan;
minv(minv==0) = nan;


shiftTimeBins = shiftvals(2:end)./sampleRate;

figure;
for s = 1:numel(speedvals)-1,   
    for p = 1:numel(phasevals)-1,   
        subplot2(numel(phasevals)-1,numel(speedvals)-1,p,s);

        ind = mind(:,s,p)~=1&mind(:,s,p)~=201;

        bar(shiftTimeBins,histc(shiftTimeBins(mind(ind,s,p)),shiftTimeBins),'histc');
        %hist2([shiftTimeBins(mind(ind,s,p))',minv(ind,s,p)],linspace(-1,2,30),linspace(40,800,6));
% $$$         hist2([shiftTimeBins(mind(ind,s,p))',minv(ind,s,p)],linspace(-1,2,30),linspace(40,800,6));
        Lines(median(shiftTimeBins(mind(ind,s,p)),'omitnan'),[],'g');
        Lines(0,[],'m');

% $$$         scatter(log10(eposTPRErrorXY([ind;false;false],50,p,s)),log10(minv(ind,s,p)),5, ...
% $$$                 shiftTimeBins(mind(ind,s,p)),'filled');
% $$$         xlim([1,2.7]);
% $$$         ylim([1,2.7]);
        
% $$$         scatter(log10([1;1;eposTPRErrorXY([ind;false;false],50,p,s,1)]),...
% $$$                 [1;1;log10(minv(ind,s,p))],5, ...
% $$$                 shiftTimeBins([1;150;mind(ind,s,p)]),'filled');
% $$$         xlim([1,2.7]);
% $$$         ylim([1,2.7]);
% $$$         scatter(([0;0;eposTPRErrorXYTrj([ind;false;false],50,p,s,1)]),...
% $$$                 [1;1;log10(minv(ind,s,p))],5, ...
% $$$                 shiftTimeBins([1;150;mind(ind,s,p)]),'filled');
% $$$         xlim([-400,400]);
% $$$         ylim([1,2.7]);
        
        xlim(shiftvals([1,end])./sampleRate);
        ylim([0,30]);
        if s==1,ylabel(num2str(mean(phasevals([p:p+1]))));end
        if p==(numel(phasevals)-1),xlabel(num2str(speedvals([s:s+1])));end
            
    end
end



