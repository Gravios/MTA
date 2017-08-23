function pXYV(Trial,key)
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');
vxy = xyz.vel(Trial.trackingMarker,[1,2]);
plot(vxy(:,1));
try,Lines(Trial.stc{key}(:),[],'r');end