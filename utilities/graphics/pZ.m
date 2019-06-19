function pZ(Trial,varargin)
[hfig,key,marker] = DefaultArgs(varargin,struct('hfig',gcf(),'key','r','marker','head_back'),'--struct');
xyz = Trial.load('xyz');
figure(hfig),hold on
plot(xyz(:,marker,3));
try,Lines(Trial.stc{key,xyz.sampleRate}(:),[],'r');end