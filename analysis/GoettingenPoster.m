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



%

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


