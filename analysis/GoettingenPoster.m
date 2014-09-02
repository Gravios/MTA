%% Goettingen Poster 2014
%  

%% Introduction
% Previous work: showed place cells specific for the rearing state.


%% Methods
MTAstartup('cin','bach')
Trial = MTATrial('jg05-20120310');
xyz = Trial.xyz.copy;xyz.load(Trial);
ang = Trial.ang.copy;ang.load(Trial);
% Equipment/Subjects
%     Vicon - Arena with Cameras
%     /data/homes/gravio/Dropbox/figures/PaperFigs/fig1c.png
%
%     Rat - Rat on maze with markers
%
%
% Walk Segmentation
%     1. JPDF mean marker speed versus marker speed variance
%        1.1 JDPF - All data 
%        1.2 JPDF - 
%     2. Examples of walking subtypes
%        2.1 High walk,

[rhm,fs,ts] = fet_rhm(Trial,[],'Swspectral');

s='c';
gper = find(diff(Trial.stc{s}.data,1,2)>140);
perind = gper(1);
perind = (Trial.stc{s}(perind,1)+20):(Trial.stc{s}(perind,1)+100);
figure,
subplot2(3,2,[1,2],1),hold on
for i= 5:8;
    plot3(xyz(perind,i,1),xyz(perind,i,2),xyz(perind,i,3))
end

plotSkeleton(xyz,perind(end));
zlim([0,200]);

subplot2(3,2,3,1)
imagesc(1:(diff(Trial.stc{s,rhm.sampleRate}(gper(1),:))+4),fs,log10(rhm((Trial.stc{s,rhm.sampleRate}(gper(1),1)-2): ...
                  (Trial.stc{s,rhm.sampleRate}(gper(1),2)+2),:))'),axis xy,,caxis([-6,-4])


%        2.2 Low walk
s='p';
gper = find(diff(Trial.stc{s}.data,1,2)>140);
perind = gper(4);
perind = (Trial.stc{s}(perind,1)+20):(Trial.stc{s}(perind,1)+100);
subplot2(3,2,[1,2],2),hold on
for i= 5:8;
    plot3(xyz(perind,i,1),xyz(perind,i,2),xyz(perind,i,3))
end
plotSkeleton(xyz,perind(end));
zlim([0,200]);

subplot2(3,2,3,2)
imagesc(1:(diff(Trial.stc{s,rhm.sampleRate}(gper(4),:))+4),fs,log10(rhm((Trial.stc{s,rhm.sampleRate}(gper(4),1)-2): ...
                  (Trial.stc{s,rhm.sampleRate}(gper(4),2)+2),:))'),axis xy,,caxis([-6,-4])




%% Example Skeletons of Behaviors
figure,
s='c';
gper = find(diff(Trial.stc{s}.data,1,2)>140);
perind = gper(1);
perind = (Trial.stc{s}(perind,1)+20):(Trial.stc{s}(perind,1)+100);
subplot2(1,3,[1],1),hold on
for i= 5:8;
    plot3(xyz(perind,i,1),xyz(perind,i,2),xyz(perind,i,3))
end

plotSkeleton(xyz,perind(end));
zlim([0,300]);
title('High Walk')
set(gca,'YTickLabelMode','manual');
set(gca,'YTickLabel',{});
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',{});
set(gca,'ZTickLabelMode','manual');
set(gca,'ZTickLabel',{});



s='p';
gper = find(diff(Trial.stc{s}.data,1,2)>140);
perind = gper(4);
perind = (Trial.stc{s}(perind,1)+20):(Trial.stc{s}(perind,1)+100);
subplot2(1,3,[1],2),hold on
for i= 5:8;
    plot3(xyz(perind,i,1),xyz(perind,i,2),xyz(perind,i,3))
end
plotSkeleton(xyz,perind(end));
zlim([0,300]);
title('Low Walk')
set(gca,'YTickLabelMode','manual');
set(gca,'YTickLabel',{});
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',{});
set(gca,'ZTickLabelMode','manual');
set(gca,'ZTickLabel',{});


s='r';
gper = find(diff(Trial.stc{s}.data,1,2)>140);
perind = gper(21);
perind = (Trial.stc{s}(perind,1)-40):(Trial.stc{s}(perind,1)+80);
subplot2(1,3,[1],3),hold on
for i= 5:8;
    plot3(xyz(perind,i,1),xyz(perind,i,2),xyz(perind,i,3))
end
plotSkeleton(xyz,perind(end));
zlim([0,300]);
title('Rearing')
set(gca,'YTickLabelMode','manual');
set(gca,'YTickLabel',{});
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',{});
set(gca,'ZTickLabelMode','manual');
set(gca,'ZTickLabel',{});



figure,

% updateSkeleton
% $$$ for kk=1:numel(sticks),
% $$$     set(sticks{kk},'XData',[xyz.data(idx,markerConnections(kk,1),1),xyz.data(idx,markerConnections(kk,2),1)],...
% $$$                    'YData',[xyz.data(idx,markerConnections(kk,1),2),xyz.data(idx,markerConnections(kk,2),2)],...
% $$$                    'ZData',[xyz.data(idx,markerConnections(kk,1),3),xyz.data(idx,markerConnections(kk,2),3)])
% $$$ end
% $$$ 
% $$$ 
% $$$ for kk=1:xyz.model.N,
% $$$     set(markers{kk},'XData',[xyz.data(perind(end),kk,1),xyz.data(perind(end),kk,1)],...
% $$$                     'YData',[xyz.data(perind(end),kk,2),xyz.data(perind(end),kk,2)],...
% $$$                     'ZData',[xyz.data(perind(end),kk,3),xyz.data(perind(end),kk,3)])
% $$$ end
% $$$ 


%% BHV Heirarchical segmentation
wper = Trial.stc{'w'}.cast('TimeSeries');
rper = Trial.stc{'r'}.cast('TimeSeries');


%% Rearing 
xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.5,xyz.sampleRate));
ang = Trial.ang.copy;
ang.create(Trial,xyz);

bang = ang(:,1,7,2);
bang = log10(xyz(:,7,3)-xyz(:,1,3));
hang = ang(:,3,4,2);

figure,
hist2([bang(nniz(bang)&nniz(hang)),...
       hang(nniz(bang)&nniz(hang))],...
      linspace(1,2.5,60),linspace(-.8,1.5,60)),
caxis([0,4000])
set(gca,'YTickLabelMode','manual');
set(gca,'YTickLabel',{});
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',{});
saveas(gcf,fullfile('/gpfs01/sirota/homes/gravio/','figures','GoettingenPoster','JPDF_hbZdiff_angSMSU_all.png'),'png');


figure,
hist2([bang(nniz(bang)&nniz(hang)&~rper.data),...
       hang(nniz(bang)&nniz(hang)&~rper.data)],...
      linspace(1,2.5,60),linspace(-.8,1.5,60)),
caxis([0,4000])
%linspace(-.1,1.5,60)
set(gca,'YTickLabelMode','manual');
set(gca,'YTickLabel',{});
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',{});
saveas(gcf,fullfile('/gpfs01/sirota/homes/gravio/','figures','GoettingenPoster','JPDF_hbZdiff_angSMSU_all_wo_rear.png'),'png');



figure,
hist2([bang(nniz(bang)&nniz(hang)&rper.data),...
       hang(nniz(bang)&nniz(hang)&rper.data)],...
      linspace(1,2.5,60),linspace(-.8,1.5,60)),
%linspace(-.1,1.5,60)
caxis([0,4000])
set(gca,'YTickLabelMode','manual');
set(gca,'YTickLabel',{});
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',{});
saveas(gcf,fullfile('/gpfs01/sirota/homes/gravio/','figures','GoettingenPoster','JPDF_hbZdiff_angSMSU_rear_only.png'),'png');


%% walk, turning, immobility
xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.1,xyz.sampleRate));
vel = xyz.vel([1:4,5,7],[1,2]);
nembed = 90;
N=10000;
shift=2000;
X = GetSegs(vel(:,end),(1+shift):(N+shift),nembed,0)'/sqrt(N);
[U,S,V] = svd(X);
Xa = zeros(vel.size([1,2]));
for i = 1:vel.size(2),
Xa(:,i) = GetSegs(vel(:,i),1:vel.size(1),nembed,0)'/sqrt(N)*V(:,1);
end
mXa = abs(median(Xa,2));
vXa = var(Xa,[],2);
vXa(vXa<.0001) = .0001;

% JPDF_mVel_vVel_all_wo_rear.png
figure,hist2([log10(abs(mXa(nniz(mXa)&nniz(vel)&~rper.data,1))),...
              log10(vXa(nniz(mXa)&nniz(vel)&~rper.data))],...
             linspace(-2.5,.75,60),linspace(-4,.5,60));
caxis([0,800]);
set(gca,'YTickLabelMode','manual');
set(gca,'YTickLabel',{});
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',{});
saveas(gcf,fullfile('/gpfs01/sirota/homes/gravio/','figures','GoettingenPoster','JPDF_mVel_vVel_all_wo_rear.png'),'png');

% JPDF_mVel_vVel_all_wo_walk.png
figure,hist2([log10(abs(mXa(nniz(mXa)&nniz(vel)&~wper.data&~rper.data,1))),...
              log10(vXa(nniz(mXa)&nniz(vel)&~wper.data&~rper.data))],...
             linspace(-2.5,.75,60),linspace(-4,.5,60));
caxis([0,800]);
set(gca,'YTickLabelMode','manual');
set(gca,'YTickLabel',{});
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',{});
saveas(gcf,fullfile('/gpfs01/sirota/homes/gravio/','figures','GoettingenPoster','JPDF_mVel_vVel_all_wo_walk.png'),'png');

% JPDF_mVel_vVel_walk_only.png'
figure,hist2([log10(abs(mXa(nniz(mXa)&nniz(vel)&wper.data,1))),...
              log10(vXa(nniz(mXa)&nniz(vel)&wper.data))],...
             linspace(-2.5,.75,60),linspace(-4,.5,60));
caxis([0,800]);
set(gca,'YTickLabelMode','manual');
set(gca,'YTickLabel',{});
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',{});
saveas(gcf,fullfile('/gpfs01/sirota/homes/gravio/','figures','GoettingenPoster','JPDF_mVel_vVel_walk_only.png'),'png');



%% Time Series Example of sniffing sub-behavior
MTAstartup
Trial = MTATrial('jg05-20120310');
[rhm,fs,ts] = fet_rhm(Trial,[],'Sspectral');
xyz = Trial.xyz.copy;xyz.load(Trial);xyz.filter(gtwin(.05,xyz.sampleRate));
ang = Trial.ang.copy;ang.create(Trial,xyz);
tx = [0:(xyz.size(1)-1)]/xyz.sampleRate;
tx = MTADxyz('data',tx,'sampleRate',xyz.sampleRate);

th = tx(Trial.stc{'hang'});
tl = tx(Trial.stc{'lang'});

hper = Trial.stc{'hang'};hper.resample(1);
lper = Trial.stc{'lang'};lper.resample(1);
rper = Trial.stc{'rear'};rper.resample(1);

sp = [];
figure,
sp(1) = subplot(2,1,1);

%area([0,0;0,300;300,0;300,300],[-2,2;2,-2],'facecolor',[.2,.5,1])

plot(tx.data,ang(:,5,7,2),'k'),hold on,
Lines([],-.47,'r')
ylabel('Head Pitch (rad)')

t = hper.data(1,:)';
psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
ptx = tx(psts);
ptx = [ptx(1),ptx,ptx(end)];
fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'b');
t = lper.data(1,:)',
psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
ptx = tx(psts);
ptx = [ptx(1),ptx,ptx(end)];
fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'g');
t = rper.data(1,:)',
psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
ptx = tx(psts);
ptx = [ptx(1),ptx,ptx(end)];
fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'r');
legend('Head Pitch','Angle Threshold','High Walk','Low Walk','Rearing');

for t = hper.data(1:20,:)',
    psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
    ptx = tx(psts);
    ptx = [ptx(1),ptx,ptx(end)];
    fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'b');
end


for t = lper.data(1:20,:)',
    psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
    ptx = tx(psts);
    ptx = [ptx(1),ptx,ptx(end)];
    fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'g');
end


for t = rper.data(1:20,:)',
    psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
    ptx = tx(psts);
    ptx = [ptx(1),ptx,ptx(end)];
    fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'r');
end


%boundedline(tx.data(1:10000),hper.data(1:10000),2*hper.data(1:10000),'b','alpha');
%boundedline(tx.data(1:10000),lper.data(1:10000),2*lper.data(1:10000),'g','alpha');

p
sp(2) = subplot(2,1,2);
imagesc([0:(rhm.size(1)-1)]/rhm.sampleRate,fs,log10(rhm.data)');
axis xy;
caxis([-5,-2.8]);
linkaxes(sp,'x');
ylabel('RHM Frequencey Hz')
xlabel('Time (s)')


%% Place fields

%MTAstartup('cin','cin');
%Trial = MTATrial('Ed10-20140812');
MTAstartup;
%Trial = MTATrial('jg05-20120317');
Trial = MTATrial('jg05-20120310');
%Trial = MTATrial('jg05-20120309');
%states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta'};
states = {'theta','rear&theta','walk&theta','hswalk&theta','lswalk&theta'};
%states = {'hswalk&theta','lswalk&theta'};
nsts = numel(states);
units = select_units(Trial,18,'pyr');

pfs ={};
for i = 1:numel(states),
pfs{i} = MTAAknnpfs(Trial,units,states{i},false,'numIter',1,'ufrShufBlockSize',0,'binDims',[10,10],'distThreshold',125,'nNearestNeighbors',110);
end

% $$$ pfs = {};
% $$$ for s = 1:nsts,
% $$$     pfs{s} =  MTAAknnpfs(Trial,units,states{s},false,'numIter',1000, ...
% $$$                          'ufrShufBlockSize',0.5,'binDims',[30,30],'distThreshold',70);
% $$$ end

[accg,tbins] = autoccg(Trial);

figure(28384)
for u = pfs{1}.data.clu,
    mrate =zeros([1,5]);
    for s = 1:nsts,
        mrate(s) = max(pfs{s}.data.rateMap(:,pfs{s}.data.clu==u));
    end
    mrate = max(mrate);
    if mrate<5,continue,end
    subplot(nsts+1,1,1);
    bar(tbins,accg(:,u));axis tight;
    set(gca,'YTickLabelMode','manual');
    set(gca,'YTickLabel',{});
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
    %set(gca,'Position',get(gca,'Position').*[0,1,1,1]+[0,0,.1475,0]);
    set(gca,'OuterPosition',get(gca,'OuterPosition').*[0,1,1.1,1.2]);
    for s = 1:nsts,
        subplot(nsts+1,1,s+1);
        pfs{s}.plot(u,[],[],[0,mrate]);
        set(gca,'YTickLabelMode','manual');
        set(gca,'YTickLabel',{});
        set(gca,'XTickLabelMode','manual');
        set(gca,'XTickLabel',{});
        %set(gca,'Position',get(gca,'Position').*[0,1,1,1]+[0,0,.1475,0]);
        set(gca,'OuterPosition',get(gca,'OuterPosition').*[0,1,1.1,1.2]);
        %if s == 1,         title([Trial.filebase '-' num2str(u)]);end
    end
    set(gcf,'paperposition',[0,0,3.4,21])
    saveas(gcf,fullfile('/gpfs01/sirota/homes/gravio/',...
                        'figures','GoettingenPoster',...
                        ['pfsShmm_' Trial.filebase '-' num2str(u) '.png']),'png');
end


%% Segmentation of walking

%MTAstartup('cin','cin');
%Trial = MTATrial('Ed10-20140813');

MTAstartup;
Trial = MTATrial('jg05-20120309');
%Trial = MTATrial('jg05-20120310');
%Trial = MTATrial('er01-20110719');
[rhm,fs,ts] = fet_rhm(Trial,[],'Sspectral');
xyz = Trial.xyz.copy;xyz.load(Trial);xyz.filter(gtwin(.05,xyz.sampleRate));
ang = Trial.ang.copy;ang.create(Trial,xyz);

nrhm = rhm.copy;
nrhm.data = log10(nrhm.data);
nrhm.data(nrhm.data<-9) = nan;

% $$$ mean_rhm =  nanmean(nrhm(nniz(nrhm),:));
% $$$ std_rhm =   nanstd(nrhm(nniz(nrhm),:));
% $$$ nrhm.data = bsxfun(@ldivide,bsxfun(@minus,nrhm.data,mean_rhm),std_rhm);

%figure,imagesc(ts,fs,nrhm.data'),axis xy,
%figure,plot(median(nrhm(:,fs>5&fs<14),2))

rhmpow =median(nrhm(:,fs>5&fs<14),2);
rhmpow = MTADlfp('data',rhmpow,'sampleRate',nrhm.sampleRate);
ang.resample(rhmpow);
xyz.resample(rhmpow);

%hper = Trial.stc{'hang'};
%lper = Trial.stc{'lang'};
%rper = Trial.stc{'rear'};

rhmlims = linspace(-6,-2,40);
anglims = linspace(-1.4,1.4,40);
figure,hist2([ang(Trial.stc{'h'},5,7,2),rhmpow(Trial.stc{'h'})],anglims,rhmlims);
figure,hist2([ang(Trial.stc{'l'},5,7,2),rhmpow(Trial.stc{'l'})],anglims,rhmlims);
figure,hist2([ang(Trial.stc{'w'},5,7,2),rhmpow(Trial.stc{'w'})],anglims,rhmlims);
set(gca,'YTickLabelMode','manual');
set(gca,'YTickLabel',{});
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',{});

bper = Trial.stc{'w',rhmpow.sampleRate}.copy;
bper.cast('TimeSeries');
xaind = nniz(ang(:,5,7,2))&nniz(rhmpow.data)&bper(1:ang.size(1));
fet = [ang(xaind,5,7,2),rhmpow(xaind)];
[bsts,bhmm,bdcd] = gausshmm(fet,2);

fasts = zeros([ang.size(1),1]);
fasts(xaind) = bsts;

figure,
plot(fasts),
Lines(Trial.stc{'w'}(:,1),[],'m');
Lines(Trial.stc{'w'}(:,2),[],'m');
ylim([-1,3])
figure,hist(ang(fasts==1,5,7,2),100)
figure,hist(ang(fasts==2,5,7,2),100)


g = 2;l = 1;
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   fasts==g,...
                   rhmpow.sampleRate,...
                   Trial.xyz.sync.copy,...
                   Trial.xyz.origin,...
                   'hswalk','h','TimeSeries');
Trial.stc.states{Trial.stc.gsi('h')}.cast('TimePeriods');

Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   fasts==l,...
                   rhmpow.sampleRate,...
                   Trial.xyz.sync.copy,...
                   Trial.xyz.origin,...
                   'lswalk','l','TimeSeries');
Trial.stc.states{Trial.stc.gsi('l')}.cast('TimePeriods');

    
%% Results

% Network Dynamics 
%     1. mean PSD with confidence intervals
%         1.1 High Walk
%         1.2 Low Walk
%     2. mean change in theta power along depth profile 
%         2.1 High walk onset
%         2.2 High walk offset
%         2.3 Low walk onset
%         2.4 Low walk offset
%     3. mean change in low gama power along depth profile 
%         3.1 High walk onset
%         3.2 High walk offset
%         3.3 Low walk onset
%         3.4 Low walk offset
%     4. mean change in high gama power along depth profile 
%         4.1 High walk onset
%         4.2 High walk offset
%         4.3 Low walk onset
%         4.4 Low walk offset
%
% Behavior dependent Place Field expression
%     1. Place Field - Walk
%     2. Place Field - High Walk
%     3. Place Field - Low Walk





%% Conclusions


