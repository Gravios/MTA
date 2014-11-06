function pZ(Trial,key)
xyz = Trial.load('xyz');
figure,hold on
plot(xyz(:,'head_front',3));
try,Lines(Trial.stc{key}(:),[],'r');end