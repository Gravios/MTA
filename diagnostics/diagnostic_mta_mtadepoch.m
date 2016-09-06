function diagnostic_mta_mtadepoch(Trial)

Trial = MTATrial('jg05-20120310');
Stc = Trial.stc.copy;
xyz = Trial.load('xyz');
xyz.resample(30);


rper = Stc{'r',30};
rper.clean;


figure,
plot(xyz(:,7,3))
Lines(rper(:),[],'r');

sum(diff(bsxfun(@minus,rper.data,[1,0]),1,2))

size(xyz(rper,7,3))