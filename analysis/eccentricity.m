
t = MTATrial('jg05-20120317');
t.ufr.create(t,t.xyz);
vel = MTADxyz([],[],t.vel(7),t.xyz.sampleRate);

units =[18,24,29,32,37,41,42];

figure,hist(t.ufr(:,units==29),100)
u = 29;
figure,hold on
uind = t.ufr(:,u==units)>10&t.ufr(:,u==units)<20&[0;vel(:)]>2;
plot3(t.xyz(uind,7,1),t.xyz(uind,7,2),t.xyz(uind,7,3),'.r')
pfs.plot(u)


figure,
ts = 20;
spkind = t.spk(18);
sind = reshape(repmat(spkind,1,ts)+repmat(1:ts,size(spkind,1),1),[],1);
sind(sind>t.xyz.size(1))=t.xyz.size(1);
plot3(t.xyz(sind,7,1),t.xyz(sind,7,2),t.xyz(sind,7,3),'.')
hold on

figure,
plot3(t.xyz(1:10:end,7,1),t.xyz(1:10:end,7,2),t.xyz(1:10:end,7,3),'.')
xlim([-500,500]),ylim([-500,500])

