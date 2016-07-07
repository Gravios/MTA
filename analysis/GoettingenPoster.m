
%% Goettingen Poster 2014
%  

%% Introduction
% Previous work: showed place cells specific for the rearing state.

nn
%% Methods
Trial = MTATrial('jg05-20120310');
xyz = Trial.load('xyz');
ang = Trial.load('ang');
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

[rhm,fs,ts] = fet_rhm(Trial,[],'wcsd');

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
s='h';
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



s='l';
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




Trial = MTATrial('jg05-20120310');
xyz = Trial.load('xyz');

%31600- 33700 turn -> walk -> rear

ftit = 'Turning';
ind = 14500;
perind = (ind-40):(ind+80);
figure,hold on
for i= 1:4;
    plot3(xyz(perind,i,1),xyz(perind,i,2),xyz(perind,i,3))
end
plotSkeleton(xyz,perind(end));
zlim([0,300]);
title(ftit)
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
xyz = Trial.load('xyz');
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


%% head marker error
Trial = MTATrial('jg05-20120310');

xyz = Trial.load('xyz');
ang = Trial.ang.copy;
ang.create(Trial,xyz);

vel = xyz.vel([5:9]);
vel.data(vel<0.01) = 0.01;
vel.data = log10(vel.data);

m = {'head_back','head_front'};
vbins = -1.5:.05:2;
[~,mb(:,1)] = histc(vel(:,m{1}),vbins);
[~,mb(:,2)] = histc(vel(:,m{2}),vbins);

mb = MTADxyz('data',mb,'sampleRate',xyz.sampleRate);

ind = Trial.stc{'a'}.cast('TimeSeries')&nniz(mb);
%ind = ':';


A = accumarray(mb(ind,:),ang(ind,m{1},m{2},3),repmat(numel(vbins),[1,2]),@mean);
B = accumarray(mb(ind,:),ang(ind,m{1},m{2},3),repmat(numel(vbins),[1,2]),@std);
S = accumarray(mb(ind,:),ind(ind),repmat(numel(vbins),[1,size(mb,2)]),@sum);
A(S<100)=nan;
B(S<100)=nan;

figure,imagescnan({vbins,vbins,A'},[],[],true,[0,0,0]),axis xy
figure,imagescnan({vbins,vbins,B'},[],[],true,[0,0,0]),axis xy

figure,
hist(ang(Trial.stc{'a'},'head_back','head_front',3),20:.02:37)
figure,
hist(ang(Trial.stc{'a'},'head_right','head_front',3),20:.02:37)
figure,
hist(ang(Trial.stc{'a'},'head_right','head_left',3),20:.02:37)
figure,
hist(ang(Trial.stc{'a'},'head_back','head_right',3),10:.02:37)


%% Time Series Example of sniffing sub-behavior
Trial = MTATrial('jg05-20120310');
[rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong',true);
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,40,'low');
ang = create(MTADang,Trial,xyz);
tx = [0:(xyz.size(1)-1)]/xyz.sampleRate;
tx = MTADxyz('data',tx,'sampleRate',xyz.sampleRate);

th = tx(Trial.stc{'hwalk'});
tl = tx(Trial.stc{'lwalk'});

hper = Trial.stc{'hwalk'};hper.resample(1);
lper = Trial.stc{'lwalk'};lper.resample(1);
rper = Trial.stc{'rear'};rper.resample(1);

sp = [];
hfig = figure,
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

for t = hper.data(1:50,:)',
    psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
    ptx = tx(1,psts);
    ptx = [ptx(1),ptx,ptx(end)];
    fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'b');
end


for t = lper.data(1:50,:)',
    psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
    ptx = tx(1,psts);
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


sp(2) = subplot(2,1,2);
imagesc([0:(rhm.size(1)-1)]/rhm.sampleRate,fs,log10(rhm.data)');
axis xy;
caxis([-5,-2.8]);
linkaxes(sp,'x');
ylabel('RHM Frequencey Hz');
xlabel('Time (s)');
colormap('jet');

%% Place fields

%MTAstartup('cin','cin');
%Trial = MTATrial('Ed10-20140812');
MTAstartup;
%Trial = MTATrial('jg05-20120317');
Trial = MTATrial('jg05-20120310');
%Trial = MTATrial('jg05-20120309');
states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta'};
%states = {'theta','rear&theta','walk&theta','hswalk&theta','lswalk&theta'};
states = {'hswalk&theta','lswalk&theta'};
nsts = numel(states);
units = select_units(Trial,18,'pyr');

%ow = false;
ow = true;
pfs ={};
for i = 1:numel(states),
pfs{i} = MTAAknnpfs(Trial,units,states{i},ow,'numIter',1,'ufrShufBlockSize',0,'binDims',[30,30],'distThreshold',125,'nNearestNeighbors',110);
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

%MTAstartup('cin','cin');
%Trial = MTATrial('Ed10-20140812');
MTAstartup;
%Trial = MTATrial('jg05-20120317');
Trial = MTATrial('jg05-20120310');
%Trial = MTATrial('jg05-20120309');

states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta';'theta','rear&theta','walk&theta','hswalk&theta','lswalk&theta'};
%states = {'hswalk&theta','lswalk&theta'};
nsts = size(states,2);



units = select_units(Trial,18,'pyr');

pfs ={};
for c = 1:size(states,1),
for i = 1:size(states,2),
%pfs{i,c} = MTAAknnpfs(Trial,units,states{c,i},false,'numIter',1,'ufrShufBlockSize',0,'binDims',[10,10],'distThreshold',125,'nNearestNeighbors',110);
    pfs{i,c} = MTAAknnpfs(Trial,units,states{c,i},false,'numIter',1,'ufrShufBlockSize',0,'binDims',[30,30],'distThreshold',125,'nNearestNeighbors',110);
end
end

% $$$ pfs = {};
% $$$ for s = 1:nsts,
% $$$     pfs{s} =  MTAAknnpfs(Trial,units,states{s},false,'numIter',1000, ...
% $$$                          'ufrShufBlockSize',0.5,'binDims',[30,30],'distThreshold',70);
% $$$ end

[accg,tbins] = autoccg(Trial);

u=units(1);
while u~=-1,
for c = 1:2,

    mrate =zeros([1,5]);
    for s = 1:nsts,
        mrate(s) = max(pfs{s}.data.rateMap(:,pfs{s}.data.clu==u));
    end
    mrate = max(mrate);
    if mrate<5,continue,end
    subplot2(nsts+1,2,1,c);
    bar(tbins,accg(:,u));axis tight;
    set(gca,'YTickLabelMode','manual');
    set(gca,'YTickLabel',{});
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
    %set(gca,'Position',get(gca,'Position').*[0,1,1,1]+[0,0,.1475,0]);
    %set(gca,'OuterPosition',get(gca,'OuterPosition').*[0,1,1.1,1.2]);
    for s = 1:nsts,
        subplot2(nsts+1,2,s+1,c);
        pfs{s,c}.plot(u,[],[],[0,mrate]);
        set(gca,'YTickLabelMode','manual');
        set(gca,'YTickLabel',{});
        set(gca,'XTickLabelMode','manual');
        set(gca,'XTickLabel',{});
        %set(gnca,'Position',get(gca,'Position').*[0,1,1,1]+[0,0,.1475,0]);
        %set(gca,'OuterPosition',get(gca,'OuterPosition').*[0,1,1.1,1.2]);
        %if s == 1,         title([Trial.filebase '-' num2str(u)]);end
    end
    %set(gcf,'paperposition',[0,0,3.4,21])

    %saveas(gcf,fullfile('/gpfs01/sirota/homes/gravio/',...
    %                    'figures','GoettingenPoster',...
    %                    ['pfsShmm_' Trial.filebase '-' num2str(u) '.png']),'png');
end
u = figure_controls(gcf,u,units);
end



figure,
subplot(121),
hist2([ang(Trial.stc{'w'},5,7,2),rhmpow(Trial.stc{'w'})],anglims,rhmlims);
title('JPFD for Walking Periods')
xlabel('head pitch');;
ylabel('rhm(6-13) pow');;
title(
subplot(122),
hist2([ang(Trial.stc{'r'}+[-2,0],5,7,2),rhmpow(Trial.stc{'r'}+[-2,0])],anglims,rhmlims);
title('JPFD for Rearing Periods and 2 seconds before rear')
xlabel('head pitch');;
ylabel('rhm(6-13) pow');;


%% CCGs and phases
Trial = MTATrial('jg05-20120317');

%Calculate bhv trans ccgs
states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta';'theta','rear&theta','walk&theta','hswalk&theta','lswalk&theta'};
nsts = size(states,2);
units = select_units(Trial,18,'pyr');
pfs ={};
for c = 1:size(states,1),
for i = 1:size(states,2),
%pfs{i,c} = MTAAknnpfs(Trial,units,states{c,i},false,'numIter',1,'ufrShufBlockSize',0,'binDims',[10,10],'distThreshold',125,'nNearestNeighbors',110);
    pfs{i,c} = MTAAknnpfs(Trial,units,states{c,i},false,'numIter',1,'ufrShufBlockSize',0,'binDims',[10,10],'distThreshold',125,'nNearestNeighbors',110);
end
end

cstates = {'theta','rear','walk','hswalk','lswalk'};
%cstates = {'theta','rear','walk','hang','lang'};

Sccg = {};
Spkb = {};
for s = 1:nsts,
Sccg{s} = gen_bhv_ccg(Trial,cstates{s},1);
Spkb{s} = Trial.spk.copy;
Spkb{s}.create(Trial,Trial.xyz.sampleRate,cstates{s},[],'deburst');
end

xyz = Trial.load('xyz');
lfp = Trial.lfp.copy;
lfp.create(Trial,[61,75,82,88]);
lfp.resample(xyz);
phs = lfp.phase;



hfig = figure(933848);
auto = true;
ny = 4;
set(hfig,'paperposition',[0,0,12,4])
u = units(1);
while u~=-1

    onm = max([Sccg{1}.ccg(:,u,1),Sccg{2}.ccg(:,u,1),Sccg{3}.ccg(:,u,1),Sccg{4}.ccg(:,u,1),Sccg{5}.ccg(:,u,1)]);
    ofm = max([Sccg{1}.ccg(:,u,2),Sccg{2}.ccg(:,u,2),Sccg{3}.ccg(:,u,2),Sccg{4}.ccg(:,u,2),Sccg{5}.ccg(:,u,2)]);
    mylim = [0,max([onm,ofm])];
    mrate =zeros([1,5]);
    for s = 1:nsts,
        mrate(s) = max(pfs{s}.data.rateMap(:,pfs{s}.data.clu==u));
    end
    mrate = max(mrate);
    %if mrate<5,continue,end

    
    for s = 1:nsts
        subplot2(ny,nsts,[1,2],s);
        pfs{s,1}.plot(u,[],[],[0,mrate]);
        set(gca,'YTickLabelMode','manual');
        set(gca,'YTickLabel',{});
        set(gca,'XTickLabelMode','manual');
        set(gca,'XTickLabel',{});

        %title(num2str(u));
        subplot2(ny,nsts,3,s);
        Sccg{s}.plot(u,1,gausswin(11)./sum(gausswin(11)));axis tight
        ylim(mylim);
        if s~=1,
                    set(gca,'YTickLabelMode','manual');
        set(gca,'YTickLabel',{});
        set(gca,'XTickLabelMode','manual');
        set(gca,'XTickLabel',{});
        
        end
        subplot2(ny,nsts,4,s);
        Sccg{s}.plot(u,2,gausswin(11)./sum(gausswin(11)));axis tight
        ylim(mylim);
        if s~=1,
                    set(gca,'YTickLabelMode','manual');
        set(gca,'YTickLabel',{});
        
        end 
       
        
% $$$         subplot2(ny,nsts,[5,6],s);
% $$$         try,circ_plot(phs(Spkb{s}(u),1),'hist',[],30,true,true),end
% $$$ try,circ_plot(circ_dist(phs(Spkb{s}(u),4),phs(Spkb{s}(u),2)),'hist',[],30,true,true),end
% $$$         subplot2(ny,nsts,[7,8],s);
% $$$         try,circ_plot(phs(Spkb{s}(u),2),'hist',[],30,true,true),end
% $$$         subplot2(ny,nsts,[9,10],s);
% $$$         try,circ_plot(phs(Spkb{s}(u),3),'hist',[],30,true,true),end
% $$$         subplot2(ny,nsts,[11,12],s);
% $$$         try,circ_plot(circ_dist(phs(Spkb{s}(u),4),phs(Spkb{s}(u),3)),'hist',[],30,true,true),end
    end

    saveas(gcf,fullfile('/gpfs01/sirota/homes/gravio/',...
                        'figures','GoettingenPoster',...
                        ['pfsCCG_' Trial.filebase '-' num2str(u) '.png']),'png');

    u = figure_controls(hfig,u,units,auto);
end
 


%% Mean RHM PSD on bhv onset
figure

mrrhm = GetSegs(rhm.data,Trial.stc{'r',rhm.sampleRate}(:,1)-24,48,nan);
%figure,
subplot(1,3,1)
imagesc(linspace(-24/rhm.sampleRate,24/rhm.sampleRate,48),fs,sq(nanmean(mrrhm,2))'),axis xy
title('mean RHM PSD triggered on rearing onset')
ylabel('frequencey')
xlabel('time (s)')
caxis([0,1.1])

mrrhm = GetSegs(rhm.data,Trial.stc{'h',rhm.sampleRate}(:,1)-24,48,nan);
subplot(1,3,2)
%figure,
imagesc(linspace(-24/rhm.sampleRate,24/rhm.sampleRate,48),fs,sq(nanmean(mrrhm,2))'),axis xy
title('mean RHM PSD triggered on high walk onset')
ylabel('frequencey')
xlabel('time (s)')
caxis([0,1.1])

mrrhm = GetSegs(rhm.data,Trial.stc{'l',rhm.sampleRate}(:,1)-24,48,nan);
%figure,
subplot(1,3,3)
imagesc(linspace(-24/rhm.sampleRate,24/rhm.sampleRate,48),fs,sq(nanmean(mrrhm,2))'),axis xy
title('mean RHM PSD triggered on low walk onset')
ylabel('frequencey')
xlabel('time (s)')
caxis([0,1.1])



%% Walks segmentation into low and high walk
nrhm = rhm.copy;
%nrhm.data = log10(nrhm.data);
%nrhm.data(nrhm.data<-9) = nan;

rhmpow =median(nrhm(:,fs>5&fs<14),2);
rhmpow = MTADlfp('data',rhmpow,'sampleRate',nrhm.sampleRate);

rang = Trial.ang.copy;rang.load(Trial);
rang.resample(rhmpow);



figure,hist2([rang(Trial.stc{'w'},5,7,2),rhmpow(Trial.stc{'w'})], ...
           linspace(-1.2,1.2,30),linspace(-.5,2,30)),caxis([0,45])
        set(gca,'YTickLabelMode','manual');
        set(gca,'YTickLabel',{});
        set(gca,'XTickLabelMode','manual');
        set(gca,'XTickLabel',{});

saveas(gcf,fullfile('/gpfs01/sirota/homes/gravio/',...
                    'figures','GoettingenPoster',...
                    ['walk_seg_walk' Trial.filebase '.png']),'png');

figure,hist2([rang(Trial.stc{'h'},5,7,2),rhmpow(Trial.stc{'h'})],linspace(-1.2,1.2,30),linspace(-0.5,2,30)),caxis([0,45])
        set(gca,'YTickLabelMode','manual');
        set(gca,'YTickLabel',{});
        set(gca,'XTickLabelMode','manual');
        set(gca,'XTickLabel',{});
saveas(gcf,fullfile('/gpfs01/sirota/homes/gravio/',...
                    'figures','GoettingenPoster',...
                    ['walk_seg_hwalk' Trial.filebase '.png']),'png');

figure,hist2([rang(Trial.stc{'l'},5,7,2),rhmpow(Trial.stc{'l'})],linspace(-1.2,1.2,30),linspace(-0.5,2,30)),caxis([0,45])
        set(gca,'YTickLabelMode','manual');
        set(gca,'YTickLabel',{});
        set(gca,'XTickLabelMode','manual');
        set(gca,'XTickLabel',{});
saveas(gcf,fullfile('/gpfs01/sirota/homes/gravio/',...
                    'figures','GoettingenPoster',...
                    ['walk_seg_lwalk' Trial.filebase '.png']),'png');



%%Final components of  poster

%% good rear ccg temporal dynamics


