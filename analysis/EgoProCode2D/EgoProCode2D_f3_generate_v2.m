

configure_default_args();
EgoProCode2D_load_data();



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


rat = load_patch_model('rat');
subject = struct(rat);
subject = update_subject_patch(subject,'head', [], false);
subject = update_subject_patch(subject,'body',[],false);
              

% CA1
tind = [3,4,5,17,18,19,20,21,22,23,29];
%tind = [6,7,26,27,30];
sampleRate = 250;
halfSpikeWindow = 0.015;

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

dca = cf(@(T,U) accumulate_decoding_vars( T, U, sampleRate, halfSpikeWindow), Trials(tind),units(tind));                                                                                                                                                                                                                                                    


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
         ...     & (dca{t}.stcm(:,5)==5)  ...
             & (dca{t}.stcm(:,3)==3|dca{t}.stcm(:,4)==4|dca{t}.stcm(:,5)==5)  ...              
              & dca{t}.hvfl(:,1)>-2                                             ...
         ...     & abs(dca{t}.hvfl(:,2))>5                                        ...
            & dca{t}.ucnt>=1 & dca{t}.ucnt<=8                                 ...
            & sqrt(sum(dca{t}.xyz(:,'hcom',[1,2]).^2,3))<325;
    decoded.fwd = cat(1,decoded.fwd,dca{t}.esax(mind,1));
    decoded.lat = cat(1,decoded.lat,dca{t}.esax(mind,2));%+20*double(t>4));
    decoded.hvf = cat(1,decoded.hvf,dca{t}.hvfl(mind,1));
    decoded.hvl = cat(1,decoded.hvl,dca{t}.hvfl(mind,2));
    decoded.hav = cat(1,decoded.hav,dca{t}.hvang(mind,1));
    decoded.hba = cat(1,decoded.hba,dca{t}.hbang(mind,1));
    decoded.phz = cat(1,decoded.phz,dca{t}.phz(mind,1));
end

ind = WithinRanges(decoded.phz,phzBin.edges(2:3)) ...
      & randn(size(decoded.hba))>0.5 ...
      & abs(decoded.hba)<1.2;
[B,BINT,R,RINT,STATS] = regress(decoded.lat(ind),[ones([sum(ind),1]),decoded.hba(ind)]);
STATS

figure,hist2([decoded.hvl,decoded.hba],linspace([-80,80,24]),linspace([-1.5,1.5,24]));


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




mBinHvl.edges = linspace(-60,60,10);
mBinHvl.centers = mean([mBinHvl.edges(1:end-1);mBinHvl.edges(2:end)]);
mBinHvl.count = numel(mBinHvl.edges)-1;
mBinHba.edges = linspace(-1.2,1.2,10);
mBinHba.centers = mean([mBinHba.edges(1:end-1);mBinHba.edges(2:end)]);
mBinHba.count = numel(mBinHba.edges)-1;
mout = zeros([mBinHba.count,mBinHvl.count]);
sout = zeros([mBinHba.count,mBinHvl.count]);
tout = zeros([mBinHba.count,mBinHvl.count]);
cout = zeros([mBinHba.count,mBinHvl.count]);
for a = 1:mBinHba.count
    for v = 1:mBinHvl.count
        ind =   WithinRanges(decoded.phz,phzBin.edges([3,4])) ...
              & WithinRanges(decoded.hvl,mBinHvl.edges([v,v+1])) ...
              & WithinRanges(decoded.hba,mBinHba.edges([a,a+1])) ...              
              & randn(size(decoded.hba))>0;
        mout(a,v) = mean(decoded.lat(ind),'omitnan');
        sout(a,v) = std(decoded.lat(ind),'omitnan');
        [~,tout(a,v)] = ttest(decoded.lat(ind));
        cout(a,v) = sum(ind);
    end
end

mask = ones([mBinHba.count,mBinHvl.count]);
mask(cout<100) = nan;
%mask(tout>abs(norminv(1-(1-0.05)^(1/sum(~isnan(mask(:))))))) = nan;
mask(tout>(1-(1-0.05)^(1/sum(~isnan(mask(:)))))) = nan;
%mask(cout<100) = nan;

figure();
%imagescnan({mBinHba.centers,mBinHvl.centers,(mout.*mask)'},[-80,80],'linear',true,'colorMap',@jet);
subplot(211);
imagescnan({mBinHba.centers,mBinHvl.centers,(mout.*mask)'},[-80,80],'linear',true,'colorMap',@jet);
axis('xy');
subplot(212);
imagescnan({mBinHba.centers,mBinHvl.centers,(sout.*mask)'},[50,150],'linear',true,'colorMap',@jet);
axis('xy');



figure,
hist2([decoded.fwd,decoded.hba],linspace(-300,300,24),linspace(-1.2,1.2,24),'xprob')
colormap('jet');
Lines([],0,'w');
Lines(0,[],'w');




figure,
subplot(131);
ind = decoded.hba>0.2;
hist2([decoded.lat(ind),decoded.fwd(ind)],linspace(-250,250,32),linspace(-200,300,32));
set(gca(),'ColorScale','log')
caxis([100,2000])
Lines(0,[],'w');
Lines([],0,'w');
d = 1
subject = struct(rat);
subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
subject = update_subject_patch(subject,'body', hbaBin.count+1-d,  true,hbaBin.edges,hbaBin.centers);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);    

subplot(132);
ind = abs(decoded.hba)<0.2;
hist2([decoded.lat(ind),decoded.fwd(ind)],linspace(-250,250,32),linspace(-200,300,32));
set(gca(),'ColorScale','log')
caxis([100,2000])
Lines(0,[],'w');
Lines([],0,'w');
d = 2
subject = struct(rat);
subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
subject = update_subject_patch(subject,'body', hbaBin.count+1-d,  true,hbaBin.edges,hbaBin.centers);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);    

subplot(133);
ind = decoded.hba<-0.2;
hist2([decoded.lat(ind),decoded.fwd(ind)],linspace(-250,250,32),linspace(-200,300,32));
set(gca(),'ColorScale','log')
caxis([100,2000])
Lines(0,[],'w');
Lines([],0,'w');
colormap('jet');
d = 3
subject = struct(rat);
subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
subject = update_subject_patch(subject,'body', hbaBin.count+1-d,  true,hbaBin.edges,hbaBin.centers);
patch(subject.body.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.head.patch.vert{:},   [0.75,0.75,0.75],'FaceAlpha',0.3);
patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);    




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
