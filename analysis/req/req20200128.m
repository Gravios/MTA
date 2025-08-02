;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.

Trial = MTATrial.validate('jg05-20120312.cof.all');


dg = LoadBinary(fullfile(Trial.spath,[Trial.name,'.dat']),82:89,96)';

periods = round(Trial.sync.data([1,end])*Trial.sampleRate);
dg = LoadBinary(fullfile(Trial.spath,[Trial.name,'.dat']),84:91,96,4,[],[],periods)';
dglfp = Trial.load('lfp',82:89);
ind = round((241*60+50-Trial.sync.data(1))*Trial.lfp.sampleRate):round((241*60+50.1-Trial.sync.data(1))*Trial.lfp.sampleRate)+3000;
dgcsd = circshift(dglfp(ind,:),-1,2)+circshift(dglfp(ind,:),1,2)-2*dglfp(ind,:);
fdgcsd = dgcsd;%RectFilter(dgcsd(:,2:6),3,21);

figure();
subplot(311);
PlotTraces(dgcsd(:,2:6),[],[],2);
subplot(312);
PlotTraces(fdgcsd,[],[],2);
subplot(313);
PlotTraces(diff(fdgcsd),[],[],2);
linkaxes();
round((241*60+50.1-Trial.sync.data(1))*Trial.lfp.sampleRate);



start = round((240*60+53.350-Trial.sync.data(1))*Trial.sampleRate)
(240*60+53.568-Trial.sync.data(1))*Trial.sampleRate
(240*60+53.350-Trial.sync.data(1))*Trial.sampleRate
start = round((241*60+50-Trial.sync.data(1))*Trial.sampleRate);

dgcsd = circshift(dgsp(ind,:),-1,2)+circshift(dgsp(ind,:),1,2)-2*dgsp(ind,:);

fdgcsd = zeros([size(dg,1),6]);
for i = 1:floor(size(dg,1)/2^20)
    disp(['i: ',num2str(i)]);
    tic
    ind = (i-1)*2^20+1:(i-1)*2^20+2^20;
    dgcsd = circshift(dg(ind,:),-1,2)+circshift(dg(ind,:),1,2)-2*dg(ind,:);
    fdgcsd(ind,:) = RectFilter(dgcsd(:,2:7),21,21);
    toc
end
i = i+1;
ind = (i-1)*2^20+1:size(dg,1);
dgcsd = circshift(dg(ind,:),-1,2)+circshift(dg(ind,:),1,2)-2*dg(ind,:);
fdgcsd(ind,:) = RectFilter(dgcsd(:,2:7),21,21);
 
% PATCH boundaries
for i = 2:floor(size(dg,1)/2^20)
    disp(['i: ',num2str(i)]);
    tic
    ind = (i-1)*2^20+1-2^10:(i-1)*2^20+2^10;
    dgcsd = circshift(dg(ind,:),-1,2)+circshift(dg(ind,:),1,2)-2*dg(ind,:);
    tfdgcsd = RectFilter(dgcsd(:,2:7),21,21);
    sind = false(size(ind));
    sind(size(ind,2)/4:size(ind,2)*3/4) = true;
    fdgcsd(ind(sind),:) = tfdgcsd(sind,:);
    toc
end

lmins = cell([1,floor(size(dg,1)/2^20)]);
lvals = cell([1,floor(size(dg,1)/2^20)]);
for i = 1:floor(size(dg,1)/2^20)
    disp(['i: ',num2str(i)]);    
    tic
    ind = (i-1)*2^20+1:(i-1)*2^20+2^20;    
    [lmins{i},lvals{i} ] = LocalMinima(-sq(sqrt(sum(bsxfun(@times,GetSegs(diff(fdgcsd(ind,3)),1:numel(ind), ...
                                                      256),gausswin(256)).^2))),256,inf);
    toc
end

tlmins = lmins;
tlvals = lvals;
% $$$ lmins = tlmins;
% $$$ lvals = tlvals;

lmins = cf(@(l,i)  l+(i-1)*2^20,  lmins,num2cell(1:numel(lmins)));
lmins = cat(1,lmins{:});
lmins = lmins+128;
lvals = cat(1,lvals{:});
lvals = -lvals;

save(fullfile(Trial.spath,[Trial.filebase,'.dgspk.mat']),'lmins','lvals','-v7.3');

figure,
hold('on');
hist(log10(lvals),1000);
mlvals = mean(log10(lvals));
slvals = std(log10(lvals));
Lines(mlvals,[],'g');
Lines(mlvals+3*slvals,[],'r');
Lines(mlvals+4*slvals,[],'r');
Lines(mlvals-3*slvals,[],'r');

% local events in the dentate gyrus molecular layer
% $$$ lmins = lmins./Trial.sampleRate;


% DETECTION END

% CLASSIFICATION START
dfdgcsd = diff(fdgcsd);
dgsEvts = GetSegs(dfdgcsd,lmins(lvals>300)-256,512);
[~,si] = sort(lvals(lvals>300));

figure,
for c = 1:6,
    subplot(1,6,c);
    imagesc(dgsEvts(:,si,c)');
end
linkaxes();

figure,
imagesc(sq(mean(dgsEvts,2))');



figure,
plot(dfdgcsd(1:2^20,3));
Lines(lmins(1:2100),[],'k');
Lines(lmins(62),[],'r');



% QUICK analysis

rper = Trial.stc{'R'};
rper = Trial.stc{'R&s'};

[tccg,t] = CCG([lmins(lvals>350);mean(rper.data,2)./rper.sampleRate.*Trial.sampleRate],...
               [ones([sum(lvals>350),1]);2*ones([size(rper,1),1])],...
               163,50,Trial.sampleRate,[1,2],'count');

figure
subplot(311);bar(t,tccg(:,2,1));xlim([t([1,end])]);
subplot(312);bar(t,tccg(:,1,1));xlim([t([1,end])]);
subplot(313);bar(t,tccg(:,2,2));xlim([t([1,end])]);


MjgER2016_load_data();

sampleRate = 250;   % Hz
spikeWindow = 0.05; % ms [0.75,0.05,0.025,0.02,0.015]
mode = 'xy'; % alt vals: 'xy'
posteriorMaxThreshold = 0.001;

% LOAD Trial
trialIndex = 20;
Trial = Trials{trialIndex}; 
% SELECT units
unitSubset = units{trialIndex};
% LOAD state collection
stc = Trial.load('stc','msnn_ppsvd_raux');
% LOAD subject position object
xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
fxyz = filter(copy(xyz),'ButFilter',3,20,'low');
% COMPUTE polar coordinates of horizontal position
mazeCenterDist = sqrt(sum(xyz(:,'hcom',[1,2]).^2,3));
mazeCenterAng = circ_dist(atan2(xyz(:,'hcom',2),xyz(:,'hcom',1)),...
                atan2(diff(xyz(:,{'hcom','nose'},2),1,2),diff(xyz(:,{'hcom','nose'},1),1,2)));

states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit','turn','walk'};
stcm = stc2mat(stc,xyz,states);

sampleRate = 250;
spikeWindow = 0.02;
%spikeWindow = 0.05;
overwrite = false;
[posEstCom,posEstMax,posEstSax,posteriorMax] = bhv_decode(Trial,sampleRate,[],unitSubset,mode,[],[],spikeWindow,overwrite);
ds = load('/storage/gravio/data/project/general/jg05-20120312/jg05-20120312.cof.all.decoded_xy_efa6ce4029c67c67ad2b8aea6ee3b84c.mat');
ds = load('/storage/gravio/data/project/general/jg05-20120312/jg05-20120312.cof.all.decoded_xy_256fb4ae77860ccb808b7b5d2a476be8.mat');
ds = load('/storage/gravio/data/project/general/jg05-20120312/jg05-20120312.cof.all.decoded_xy_ec622ef3590c960961960e07d73911a8.mat');
posEstCom = ds.posEstCom;
posEstMax = ds.posEstMax;
posteriorMax = ds.posteriorMax;

hvec = fxyz(:,'hcom',[1,2])-fxyz(:,'nose',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);

decError = multiprod(posEstMax(:,[1,2])-sq(xyz(:,'hcom',[1,2])),hvec(:,:,:),2,[2,3]);
decError = multiprod(posEstCom(:,[1,2])-sq(xyz(:,'hcom',[1,2])),hvec(:,:,:),2,[2,3]);
decError(:,2) = decError(:,2)+16;
decError(:,1) = decError(:,1)-45;
ind = stcm(:,1)==1;
for j = 1:135,
    decError(:,1) = decError(:,1)+1;
    mde(j) = median(sqrt(sum(decError(ind,:).^2,2)),'omitnan');
end
figure,plot(mde)
decError(:,1) = decError(:,1)-(135-LocalMinima(mde,30,[],1));


figure,hist(sqrt(sum(decError(ind,:).^2,2)),linspace(0,1000,1000));

ind = round(lmins(lvals>350)./Trial.sampleRate.*xyz.sampleRate);
ind = ind(stcm(ind,8)==8);

derr = GetSegs(decError,ind-1000,2000);
figure,imagescnan(sqrt(sum(derr.^2,3))');

figure();
hold('on');
plot(nanmean(sqrt(sum(derr.^2,3)),2));
plot(nanmean(sqrt(sum(derr.^2,3)),2)+nanstd(sqrt(sum(derr.^2,3)),[],2),'r');


% $$$ ind = round(lmins(lvals>350)./Trial.sampleRate.*xyz.sampleRate);
ind = round(lmins./Trial.sampleRate.*xyz.sampleRate);
slvals = lvals;
slvals = slvals(stcm(ind,8)==8&stcm(ind,1)~=1&posteriorMax(ind)>0.0025);
ind = ind(stcm(ind,8)==8&stcm(ind,1)~=1&posteriorMax(ind)>0.0025);

% $$$ ind = ind(stcm(ind,8)==8&stcm(ind,1)~=1&[false;diff(lmins(lvals>350))>100]);

figure
hist2([slvals,sqrt(sum(decError(ind,:).^2,2))],linspace(0,600,50),linspace(0,800,100));
caxis([0,100])

nbins = 200;

figure,
subplot(311);
hist(sqrt(sum(decError(stcm(:,8)==8&stcm(:,1)~=1,:).^2,2)),linspace(0,1000,nbins))
subplot(312);
hist(sqrt(sum(decError(ind,:).^2,2)),linspace(0,1000,nbins))
subplot(313);
hist(sqrt(sum(decError(ind+250,:).^2,2)),linspace(0,1000,nbins))


%ind = round((242*60+45-Trial.sync(1)).*250);
start = round((242*60+45-Trial.sync(1)).*250);
start = round((241*60+52-Trial.sync(1)).*250);
start = round((244*60+50-Trial.sync(1)).*250);
start = round((245*60+12-Trial.sync(1)).*250);
start = 4172*250

start = 4499*250
start = 4545*250;
ind = start:start+3*250;
ind = ind(posteriorMax(ind)>0.001);

posest = posEstMax;
posest = posEstCom;

figure();
hold('on');
plot(xyz(ind,'hcom',1),xyz(ind,'hcom',2),'.')
plot(posest(ind,1),posest(ind,2),'.r');
plot(posest(ind(1),1),posest(ind(1),2),'og');
plot(posest(start+25,1),posest(start+75,2),'oc');
xlim([-500,500]);
ylim([-500,500]);

sqrt(sum(decError(stcm(:,8)==8&stcm(:,1)~=1,:).^2,2))


pyr = Trial.load('lfp',70);
rad = Trial.load('lfp',77);
lmo = Trial.load('lfp',83);
dgm = Trial.load('lfp',87);
dgc = Trial.load('lfp',89);
hil = Trial.load('lfp',96);


phzPyr = phase(copy(pyr),[5,12]);
phzRad = phase(copy(rad),[5,12]);
phzLmo = phase(copy(lmo),[5,12]);


vxy = vel(filter(copy(xyz),'ButFilter',4,5,'low'),{'spine_lower','hcom'},[1,2,3]);
vxy = fet_href_HXY(Trial,250,[],'trb');

ind = 4499*1250;
ind = 4172*1250;
ind = 4545*1250;

start = 242*60+25-10060;

w = 100;
start = 1400;

start = 495;
ind = start*1250;
ind = ind:ind+w*1250;
vind = start*250;
vind = vind:vind+w*250;



figure,
hold('on');
plot([0:w*1250]'./1250,circ_dist(phzRad(ind),phzLmo(ind)).*10,'r');
plot([0:w*1250]'./1250,circ_dist(phzRad(ind),phzPyr(ind)).*10,'b');
plot([0:w*1250]'./1250,circ_dist(phzPyr(ind),phzLmo(ind)).*10,'c');

plot([0:w*250]'./250,(vxy(vind,1))./8,'k');
plot([0:w*250]'./250,(vxy(vind,2))./8,'g');

plot([1:w*250]'./250,diff(vxy(vind,1)).*4,'k');
plot([1:w*250]'./250,diff(vxy(vind,2)).*4,'g');



dphzRL = interp1([1:size(phzRad,1)]./1250,circ_dist(phzLmo.data,phzPyr.data),[1:size(vxy,1)]/250);

fet = copy(vxy);
fet.data = cat(2,fet.data,dphzRL');

defspec = struct('nFFT',2^10,'Fs',250,...
                 'WinLength',2^9,'nOverlap',2^9*.875,...
                 'FreqRange',[1,10]);
[ysz,fsz,tsz,phi] = fet_spec(Trial,fet,'mtchglong','defspec',defspec,'flagCrossSpec',true);

figure,
subplot(411);
imagesc(tsz,fsz,log10(ysz(:,:,2,2))');caxis([-4,0]);
subplot(412);
imagesc(tsz,fsz,ysz(:,:,1,3)');caxis([0.5,1]);
subplot(413);
imagesc(tsz,fsz,ysz(:,:,2,3)');caxis([0.5,1]);
subplot(414);
imagesc(tsz,fsz,phi(:,:,2,3)');colormap(gca,'hsv');caxis([-pi,pi]);
linkaxes();

vxy(vind,1),
help mtcdglong

defspec = struct('nFFT',2^11,'Fs',pyr.sampleRate,...
                 'WinLength',2^10,'nOverlap',2^10*.875,...
                 'FreqRange',[1,50]);
[ys,fs,ts] = fet_spec(Trial,pyr,'defspec',defspec);
[ysr,fs,ts] = fet_spec(Trial,rad,'defspec',defspec);
[ysl,fs,ts] = fet_spec(Trial,lmo,'defspec',defspec);

defspec = struct('nFFT',2^8,'Fs',pyr.sampleRate,...
                 'WinLength',2^7,'nOverlap',2^7*.875,...
                 'FreqRange',[25,200]);
[gys,gfs,gts] = fet_spec(Trial,pyr,'defspec',defspec);
[gysr,gfsr,gtsr] = fet_spec(Trial,rad,'defspec',defspec);
[gysl,gfsl,gtsl] = fet_spec(Trial,lmo,'defspec',defspec);

figure();
sp = tight_subplot(6,1,[0.01,0.01,],0.1);
i = 1;
axes(sp(i));i=i+1;
imagesc(ts,fs,log10(ys.data)');
caxis([1,3])
colormap('jet');
axes(sp(i));i=i+1;
imagesc(ts,fs,log10(ysr.data)');
caxis([1,3])
colormap('jet');
axes(sp(i));i=i+1;
imagesc(ts,fs,log10(ysl.data)');
caxis([1,3])
colormap('jet');
axes(sp(i));i=i+1;
imagesc(gts,gfs,log10(gys.data)');
caxis([2,3.5]);
axes(sp(i));i=i+1;
imagesc(gtsr,gfsr,log10(gysr.data)');
caxis([2,3.5]);
axes(sp(i));i=i+1;
imagesc(gtsl,gfsl,log10(gysl.data)');
caxis([2,3.5]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');

pyr = LoadBinary(fullfile(Trial.spath,[Trial.name,'.dat']),70,96,[],[],[],ind([1,end]))';


[ys,fs,ts,phi,fstat] = mtchglong(WhitenSignal(lm,[],1),...
                        2^12,...
                        32552,...
                        2^11,...
                        2^11.*0.875,[],[],[],[15,200]);


figure();
sp = tight_subplot(11,2,0,0);
sp = reshape(sp',[2,11])';
for i = 1:11,
    axes(sp(i,1));
    imagesc(ts,fs,sq(circ_dist(phi(:,i*2,12,:),circshift(phi(:,i*2,12,:),-1,4)))');
    %imagesc(ts,fs,sq(phi(:,i*2,12,:))');
    colormap('hsv'); 
    caxis([-pi/2,pi/2]);
    axes(sp(i,2));    
    tys = [];
    for c = 1:32,
    tys(:,c) =     sq(log10(ys(:,i*2,c,c)));
    end
    imagesc(ts,fs,tys');
    colormap('jet');
    title(num2str(fs(i*2)))
end
linkaxes