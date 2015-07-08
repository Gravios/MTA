function pXY(Trial,key)
xyz = Trial.load('xyz');
figure,hold on
plot(xyz(:,Trial.trackingMarker,1),xyz(:,Trial.trackingMarker,2),'.');
try,Lines(Trial.stc{key}(:),[],'r');end