function pXY(Trial,key)
xyz = Trial.load('xyz');
figure,hold on
plot(xyz(:,Trial.trackingMarker,1),xyz(:,Trial.trackingMarker,2),'.');
alims(gca,Trial.maze.visible_volume);
try,Lines(Trial.stc{key}(:),[],'r');end