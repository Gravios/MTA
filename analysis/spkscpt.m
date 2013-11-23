
state = 'rear'
bspx = s.spkxyz(state);


% $$$ unit = 17;
% $$$ marker = 7;
% $$$ 
% $$$ figure
% $$$ sp1 = subplot(121)
% $$$ plot3(spx{unit}(:,marker,1),spx{unit}(:,marker,2),spx{unit}(:,marker,3),'.')
% $$$ sp2 = subplot(122)
% $$$ plot3(rpx{unit}(:,marker,1),rpx{unit}(:,marker,2),rpx{unit}(:,marker,3),'.')


spkvec = cell(length(bspx),1);
for unit = 17%1:length(bspx)
for i =2:size(rpx{unit},1)
traj = 1;
while sqrt(sum((rpx{unit}(i,s.Model.gmi(s.trackingMarker),:)-rpx{unit}(i-1,s.Model.gmi(s.trackingMarker),:)).^2,3))<50
spkvec{unit,traj}(i,:) = sq(rpx{unit}(i,s.Model.gmi(s.trackingMarker),:));
end
traj = traj+1;
end
end
