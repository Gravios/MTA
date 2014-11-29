function pB(Trial)
xyz = Trial.load('xyz');
figure,hold on
plot(xyz(:,Trial.trackingMarker,3));
Lines((xyz.sync.data(:,1)-xyz.sync.data(1))*xyz.sampleRate,[],'r');
Lines((xyz.sync.data(:,2)-xyz.sync.data(1))*xyz.sampleRate,[],'k');