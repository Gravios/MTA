function interneuron = accumulate_interneuron_stats(Trial,units,thetaRef,phzCorrection)
% It got crazy twoards the end
unitsPyr = Trial.spk.get_unit_set(Trial,'pyramidal');
unitsInt = Trial.spk.get_unit_set(Trial,'interneurons');
units = [unitsPyr,unitsInt];

phz  = load_theta_phase(Trial,Trial.lfp.sampleRate,thetaRef,phzCorrection);
spkt = Trial.load('spk',Trial.lfp.sampleRate,'theta-sit-groom',units,'');
phzLim = [0,2*pi];
phzBin = linspace([phzLim,36]);
phzCtr = mean([phzBin(1:end-1);phzBin(2:end)]);

% COMPUTE ccg of interneuron and SWR
spks = Trial.load('spk',Trial.lfp.sampleRate,'sit',units,'');
ripples = [Trial.stc{'R&sit',Trial.lfp.sampleRate}];

rccg = [];
for u = 1:numel(units),
    uRes = spks(units(u));
    [tccg,tbin] = CCG([uRes;mean(ripples.data,2)],...
                      [ones([size(uRes,1),1]);2*ones([size(ripples,1),1])],...
                      16,40,spks.sampleRate,[1,2],'scale');
    rccg(:,u) = tccg(:,2,1);
end

fccg = RectFilter(bsxfun(@minus,rccg,mean(rccg(tbin<-200|tbin>200,:))),3,3);

gvect = 1./(50.*sqrt(2*pi)).*exp(-0.5.*(tbin.^2./50.^2));
gvect = gvect./sqrt(sum(gvect.^2));
gvectL = 1./(25.*sqrt(2*pi)).*exp(-0.5.*((tbin+50).^2./25.^2));
gvectL = gvectL./sqrt(sum(gvectL.^2));
gvectR = 1./(25.*sqrt(2*pi)).*exp(-0.5.*((tbin-50).^2./25.^2));
gvectR = gvectR./sqrt(sum(gvectR.^2));


clear('interneuron');
for u = 1:numel(units)
    interneuron(u).tpp = circ_mean(phz(spkt(units(u))));
    interneuron(u).tpr = circ_r(phz(spkt(units(u))));
    if interneuron(u).tpp<0
        interneuron(u).tpp = interneuron(u).tpp+2*pi;
    end
end

for u = 1:numel(units)    
interneuron(u).tailLeft   = mean(fccg(tbin>-500&tbin<-200,u));
interneuron(u).centerLeft = sum(fccg(:,u).*gvectL');
interneuron(u).center     = sum(fccg(:,u).*gvect');
interneuron(u).centerRight= sum(fccg(:,u).*gvectR');
interneuron(u).tailRight  = mean(fccg(tbin<500&tbin>200,u));
end


spkr = Trial.load('spk',Trial.lfp.sampleRate,'R&s',unitsInts{20},'');
[tccg,tbin] = CCG(spkr.res,spkr.clu,1,40,spkr.sampleRate,unique(spkr.clu),'scale');
[tccg,tbin] = CCG(spks.res,spks.clu,16,40,spks.sampleRate,unique(spks.clu),'scale');
[tccg,tbin] = CCG(spkt.res,spkt.clu,2,40,spkt.sampleRate,unique(spkt.clu),'scale');

figure,
sp = tight_subplot(14,14,0,0);
for u1 = 1:17
    for u2 = u1:17
        axes(sp((u1-1)*14+u2));
        bar(tbin,tccg(:,u1,u2));
        Lines(0,[],'r');
    end
end


% unitsPyr 
%      5    17    20    21    22    23    24    25    28    31    32    33    34
%     35    37    41    44    48    51    52    59    60    61    63    68    72
%     73    74    79    80    81    83    85    86    89    90    93    96    98
%    102   103   104   105   109   110   111   113   116   118   119   120   128
%    129   132   133   134   136   138   139   140   141   142   144   145   149
%    151   159   179   181

%unitsInt
%3     7     8    15    16    43    45    50    76    77    92   106   124   184
cunits = unique(spkt.clu);

figure,
for u = unitsPyr
    clf();
sp = tight_subplot(4,4,0.01,0);
unitsIntInd = find(ismember(cunits,unitsInt));
punit = find(cunits==u);
for u1 = 1:numel(unitsIntInd)
    axes(sp(u1));
    bar(tbin,tccg(:,punit,unitsIntInd(u1)));
    ylabel(num2str(u));
    %Lines(0,[],'r',1);
end
waitforbuttonpress();
end



xyzH = preproc_xyz(Trial,'trb');
xyzH.data = sq(xyzH(:,{'hcom','nose'},[1,2]));
xyzH.resample(Trial.lfp.sampleRate);
spkt = Trial.load('spk',Trial.lfp.sampleRate,'walk&theta',unitsInt,'');
pft = pfs_2d_theta(Trial,units);

hvec = copy(xyzH);
hvec.data = xyzH(:,2,:)-xyzH(:,1,:);
hvec.data = sq(bsxfun(@rdivide,hvec.data,sqrt(sum(hvec.data.^2,3))));
hvec.data = cat(3,hvec.data,sq(hvec.data)*[0,-1;1,0]);
hRot = 0.264;
hvec.data = multiprod(hvec.data,                  ...
                 [cos(hRot),-sin(hRot); ...
                  sin(hRot), cos(hRot)],...
                 [2,3],                 ...
                 [1,2]);

shvc = hvec(spkt.res,:,:);
spos = sq(xyzH(spkt.res,1,:));


targetPos = [100,100];
sposH = multiprod(bsxfun(@minus,spos,targetPos),shvc,2,[2,3]);
targetDst = sqrt(sum(sposH.^2,2));
sposHLim = [-200,200];
sposHBin = linspace([sposHLim,11]);
sposHCtr = mean([sposHBin(1:end-1);sposHBin(2:end)]);
sposHInd = discretize(sposH(:,1),sposHBin);
dind = sposH(:,2)<100;
tccg = [];
for p = 1:numel(sposHCtr),
rind = dind &sposHInd==p;
[tccg(:,:,:,p),tbin] = CCG(spkt.res(rind),spkt.clu(rind),8,30,spkt.sampleRate,unique(spkt.clu),'hz');
end
out = RectFilter(sq(tccg(:,6,11,:)),3,3);
out = bsxfun(@rdivide,out,sum(out));

for f = 13%1:14
figure,
sp = tight_subplot(4,4,0.01,0);
for p = 1:14,
    axes(sp(p));
    out = RectFilter(sq(tccg(:,p,f,:)),3,3);
    %out = bsxfun(@rdivide,out,sum(out));
    imagesc(tbin,sposHCtr,out');axis('xy');colormap('jet');
end
end


ppos = sq(xyzH(Trial.stc{'x&t'},1,:));
pposH = multiprod(bsxfun(@minus,ppos,targetPos),hvec(Trial.stc{'x&t'},:,:),2,[2,3]);
sposH = multiprod(bsxfun(@minus,spos,targetPos),shvc,2,[2,3]);
[ratemap,Bins] = PlotPF(Trial,sposH,pposH(




%unitsInt
%3     7     8    15    16    43    45    50    76    77    92   106   124   184


figure
sp = tight_subplot(14,14,0.01,0);

for u1 = 1:numel(unitsIntInd)
    for u2 = u1:numel(unitsIntInd)        
        axes(sp((u1-1)*14+u2));   
        if u1==u2,
            plot(pft,unitsInt(u1),1,'text','colorMap',@jet);
        else
            %continue
            out = RectFilter(sq(tccg(:,u1,u2,:)),3,3);
            %out = bsxfun(@rdivide,out,sum(out));
            imagesc(tbin,sposHCtr,out');axis('xy');colormap('jet');
        end
    end
end


% ego
spkt = Trial.load('spk',16,'walk&theta',units,'');
xyz = preproc_xyz(Trial,'trb');
xyz.data = sq(xyz(:,{'hcom','nose'},[1,2]));
xyz.resample(16);

hvec = copy(xyz);
hvec.data = xyz(:,2,:)-xyz(:,1,:);
hvec.data = sq(bsxfun(@rdivide,hvec.data,sqrt(sum(hvec.data.^2,3))));
hvec.data = cat(3,hvec.data,sq(hvec.data)*[0,-1;1,0]);
hRot = 0%0.264;
hvec.data = multiprod(hvec.data,                  ...
                 [cos(hRot),-sin(hRot); ...
                  sin(hRot), cos(hRot)],...
                 [2,3],                 ...
                 [1,2]);


hang = atan2(hvec(:,1,1),hvec(:,1,2));

x = -400:50:400;
y = -400:50:400;
ratemap  = [];
for xind = 1:numel(x)
    for yind = 1:numel(y)
        shvc = hvec(spkt.res(spkt.clu==106),:,:);
        spos = sq(xyz(spkt.res(spkt.clu==106),1,:));
        targetPos = [x(xind),y(yind)];
        sposH = multiprod(bsxfun(@minus,targetPos,spos),shvc,2,[2,3]);
        ppos = sq(xyz(Trial.stc{'x&t'},1,:));
        pposH = multiprod(bsxfun(@minus,targetPos,ppos),hvec(Trial.stc{'x&t'},:,:),2,[2,3]);
        sposH = multiprod(bsxfun(@minus,targetPos,spos),shvc,2,[2,3]);
        [ratemap(:,xind,yind),Bins] = PlotPF(Trial,sposH,pposH,[20,20],[2.3,2.3],'xy',[-200,200;-200,200],16);
    end
end

figure();
for xind = (2:2:16)
    for yind = (2:2:16)
        subplot2(8,8,yind/2,xind/2);
        imagesc(reshape(ratemap(:,xind,yind)',[40,40]));
        caxis([0,40]);
        colormap('jet');
    end
end


figure();
for xind = 1:30
    subplot(5,6,xind);
    imagesc(sq(ratemap(sub2ind([30,30],xind,15),:,:)));
    colorbar();
end

    figure,imagesc(imagesc(sq(ratemap(sub2ind([30,30],20,15),:,:))');


Trial = Trials{20};
unitsPyr = Trial.spk.get_unit_set(Trial,'pyramidal');
unitsInt = Trial.spk.get_unit_set(Trial,'interneurons');
units = [unitsPyr,unitsInt];

    
    
xyz = preproc_xyz(Trial,'trb');
xyz.data = sq(xyz(:,{'hcom','nose'},[1,2]));
xyz.resample(16);

hvec = copy(xyz);
hvec.data = xyz(:,2,:)-xyz(:,1,:);
hvec.data = sq(bsxfun(@rdivide,hvec.data,sqrt(sum(hvec.data.^2,3))));
hvec.data = cat(3,hvec.data,sq(hvec.data)*[0,-1;1,0]);
hRot = 0%0.264;
hvec.data = multiprod(hvec.data,                  ...
                 [cos(hRot),-sin(hRot); ...
                  sin(hRot), cos(hRot)],...
                 [2,3],                 ...
                 [1,2]);


hang = copy(hvec);
hang.data = atan2(hvec(:,1,2),hvec(:,1,1));


ind = Trial.stc{'hbhv&t'};
spkt = Trial.load('spk',16,'hbhv&theta',units,'');

indP = Trial.stc{'lbhv&t'};
spktP = Trial.load('spk',16,'lbhv&theta',units,'');

ind = Trial.stc{'walk&t'};
spkt = Trial.load('spk',16,'walk&theta',units,'');

indP = Trial.stc{'pause&t'};
spktP = Trial.load('spk',16,'pause&theta',units,'');

ind = Trial.stc{'hloc&t'};
spkt = Trial.load('spk',16,'hloc&theta',units,'');

indP = Trial.stc{'lloc&t'};
spktP = Trial.load('spk',16,'lloc&theta',units,'');

ind = Trial.stc{'lloc&t'};
spkt = Trial.load('spk',16,'lloc&theta',units,'');

indP = Trial.stc{'lpause&t'};
spktP = Trial.load('spk',16,'lpause&theta',units,'');


ind = Trial.stc{'walk+pause&t'};
spkt = Trial.load('spk',16,'walk+pause&theta',units,'');

indP = Trial.stc{'rear&t'};
spktP = Trial.load('spk',16,'rear&theta',units,'');


mask = create_tensor_mask(Bins([1,3]));

%4     6     7     8    27    28    43    59    71    86    99   100   101    102
%3     7     8    15    16    43    45    50    76    77    92   106   124    184
hangInd = discretize(hang.data,linspace(-pi,pi,5));

u = 106;
figure
ratemap  = [];
sind = Trial.stc{'walk+pause&t'};
sind.cast('TimeSeries');
sind.resample(16);
shvc = hvec(spkt.res(spkt.clu==u),:,:);
sposH = [xyz(spkt.res(spkt.clu==u),1,1),hang(spkt.res(spkt.clu==u)),xyz(spkt.res(spkt.clu==u),1,2)];
pposH = [xyz(ind,1,1),hang(ind,1),xyz(ind,1,2)];
[ratemap,Bins] = PlotPFCirc(Trial,sposH,pposH,[20,pi/2,20],[3.3,0.8,3.3],'xy',[-500,500;-pi,pi;-500,500],16);
rmap = reshape(ratemap,cellfun(@numel,Bins));
clim = [0,prctile(rmap(:),99)];
for p = 1:numel(Bins{2})
    subplot2(2,numel(Bins{2})+2,1,p);
    imagesc(Bins{1},Bins{3},(sq(rmap(:,p,:)).*mask)');
    axis('xy');
    colormap(gca(),'jet');
    %caxis(clim);
    caxis([0,50]);    
    if p == 1
        ylabel(num2str(u));
    end
    hold('on');
    plot(xyz(sind.data&hangInd==p,1,1),xyz(sind&hangInd==p,1,2),'.')
    xlim([-500,500]);
    ylim([-500,500]);
end
phzRateMeanCpx = rmap.*repmat(exp(i.*permute(Bins{2},[2,1,3])),[numel(Bins{1}),1,numel(Bins{3})]);
phzRateMeanCpxMean = angle(sq(bsxfun(@rdivide,sum(phzRateMeanCpx,2,'omitnan'),sum(rmap,2,'omitnan'))));
phzRateMeanCpxRl = abs(sq(bsxfun(@rdivide,sum(phzRateMeanCpx,2,'omitnan'),sum(rmap,2,'omitnan'))));
subplot2(2,numel(Bins{2})+2,1,p+1);
imagesc(Bins{1},Bins{3},(phzRateMeanCpxMean.*mask)');
axis('xy');
colormap(gca(),'hsv');
caxis([-pi,pi]);
subplot2(2,numel(Bins{2})+2,1,p+2);
imagesc(Bins{1},Bins{3},(phzRateMeanCpxRl.*mask)');
axis('xy');
colormap(gca(),'jet');
caxis([0,0.3]);
 
sind = Trial.stc{'rear&t'};
sind.cast('TimeSeries');
sind.resample(16);
 
ratemap  = [];
shvc = hvec(spktP.res(spktP.clu==u),:,:);
sposH = [xyz(spktP.res(spktP.clu==u),1,1),hang(spktP.res(spktP.clu==u)),xyz(spktP.res(spktP.clu==u),1,2)];
pposH = [xyz(indP,1,1),hang(indP,1),xyz(indP,1,2)];
[ratemap,Bins] = PlotPFCirc(Trial,sposH,pposH,[20,pi/2,20],[3.3,0.8,3.3],'xy',[-500,500;-pi,pi;-500,500],16);
rmap = reshape(ratemap,cellfun(@numel,Bins));
for p = 1:numel(Bins{2})
    subplot2(2,numel(Bins{2})+2,2,p);
    imagesc(Bins{1},Bins{3},(sq(rmap(:,p,:)).*mask)');
    axis('xy');
    colormap(gca(),'jet');
    %caxis(clim);
    caxis([0,50]);
    hold('on');
    plot(xyz(sind.data&hangInd==p,1,1),xyz(sind&hangInd==p,1,2),'.')
    xlim([-500,500]);
    ylim([-500,500]);
end
phzRateMeanCpx = rmap.*repmat(exp(i.*permute(Bins{2},[2,1,3])),[numel(Bins{1}),1,numel(Bins{3})]);
phzRateMeanCpxMean = angle(sq(bsxfun(@rdivide,sum(phzRateMeanCpx,2,'omitnan'),sum(rmap,2,'omitnan'))));
phzRateMeanCpxRl = abs(sq(bsxfun(@rdivide,sum(phzRateMeanCpx,2,'omitnan'),sum(rmap,2,'omitnan'))));
subplot2(2,numel(Bins{2})+2,2,p+1);
imagesc(Bins{1},Bins{3},(phzRateMeanCpxMean.*mask)');
axis('xy');
colormap(gca(),'hsv');
caxis([-pi,pi]);
subplot2(2,numel(Bins{2})+2,2,p+2);
imagesc(Bins{1},Bins{3},(phzRateMeanCpxRl.*mask)');
axis('xy');
colormap(gca(),'jet');
caxis([0,0.3]);



sind = Trial.stc{'walk+pause&t'};
sind.cast('TimeSeries');
sind.resample(16);
figure
hangInd = discretize(hang.data,linspace(-pi,pi,5));
for a = 1:4,
    subplot(1,4,a);
    plot(xyz(sind.data&hangInd==a,1,1),xyz(sind&hangInd==a,1,2),'.')
    xlim([-500,500]);
    ylim([-500,500]);
end

