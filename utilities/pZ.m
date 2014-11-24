function pZ(Trial,key)
xyz = Trial.load('xyz');
figure,hold on
plot(xyz(:,Trial.trackingMarker,3));
try,Lines(Trial.stc{key}(:),[],'r');end