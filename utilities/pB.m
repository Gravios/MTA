function pB(Trial)
xyz = Trial.load('xyz');
figure,hold on
plot(xyz(:,Trial.trackingMarker,3));
Lines((s.xyz.sync.data(:,1)-s.xyz.sync.data(1))*s.xyz.sampleRate,[],'r')
Lines((s.xyz.sync.data(:,2)-s.xyz.sync.data(1))*s.xyz.sampleRate,[],'k')