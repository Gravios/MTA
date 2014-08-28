%% Goettingen Poster 2014
%  

%% Introduction



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

figure,hist(ang(Trial.stc{'c'},'head_back','head_front',2),100)
figure,hist(ang(Trial.stc{'p'},'head_back','head_front',2),100)


s='c';
gper = find(diff(Trial.stc{s}.data,1,2)>100);
perind = gper(1);
perind = Trial.stc{s}(perind,1):Trial.stc{s}(perind,2)-20;
figure,hold on,
for i= 5:8;
    plot3(xyz(perind,i,1),xyz(perind,i,2),xyz(perind,i,3))
end
plotSkeleton(xyz,perind(end));
zlm = zlim;

%        2.2 Low walk


s='p';
gper = find(diff(Trial.stc{s}.data,1,2)>100);
perind = gper(3);
perind = Trial.stc{s}(perind,1):Trial.stc{s}(perind,2);
figure,hold on,
for i= 5:8;
    plot3(xyz(perind,i,1),xyz(perind,i,2),xyz(perind,i,3))
end
plotSkeleton(xyz,perind(end));
zlim(zlm);

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


