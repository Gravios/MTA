


Trial = MTATrial('Ed10-20140815');
ds = load(fullfile(Trial.path.data,'Ed10-20140815.Breath_classification.lfp.66.mat'));
ds.insp_periods = ds.insp_periods - round(Trial.sync.data(1)*Trial.lfp.sampleRate);

% drop periods before vicon
ind = ds.insp_periods(:,1)<0;
ds.insp_periods(ind,:) = [];
ds.features(ind,:) = [];
ds.clus(ind,:) = [];
ds.breath_trains(ind,:) = [];

% downsample timestamps
ds.insp_periods = round(ds.insp_periods./Trial.lfp.sampleRate.*Trial.xyz.sampleRate);

xyz = Trial.load('xyz');

% drop periods after vicon
ind = ds.insp_periods(:,2)>xyz.size(1);
ds.insp_periods(ind,:) = [];
ds.features(ind,:) = [];
ds.clus(ind,:) = [];
ds.breath_trains(ind,:) = [];



[at,aind] =  SelectPeriods(ds.insp_periods(:,1),Trial.stc{'a'}.data,'d',1,0);
[wt,wind] =  SelectPeriods(ds.insp_periods(:,1),Trial.stc{'w'}.data,'d',1,0);
[rt,rind] =  SelectPeriods(ds.insp_periods(:,1),Trial.stc{'r'}.data,'d',1,0);
[ht,hind] =  SelectPeriods(ds.insp_periods(:,1),Trial.stc{'h'}.data,'d',1,0);
[lt,lind] =  SelectPeriods(ds.insp_periods(:,1),Trial.stc{'l'}.data,'d',1,0);

f = [2,1];
figure,
plot(ds.features(:,f(1)),ds.features(:,f(2)),'.')
hold on
plot(ds.features(wind,f(1)),ds.features(wind,f(2)),'.r')
plot(ds.features(rind,f(1)),ds.features(rind,f(2)),'.g')

xyz.filter(gtwin(.75,xyz.sampleRate));
vel = xyz.vel('spine_lower',[1,2]);
vel.data = log10(vel.data);
ang = Trial.ang.copy;
ang.create(Trial,xyz);




fet = 1;
figure
plot(xyz(ds.insp_periods(:,1),7,3),ds.features(:,fet),'.')


figure
hist2(log10([xyz(ds.insp_periods(aind,1),7,3),ds.features(aind,fet)]+10),1.2:.01:2.5,2:.01:3)
caxis([0,100])

plot(ang(ds.insp_periods(:,1),5,7,2),ds.features(:,fet),'.')

nind = nniz(vel);

figure,
plot(vel(ds.insp_periods(:,1)),ds.features(:,3),'.')
hold on
% $$$ plot(vel(ds.insp_periods(wind,1)),log10(ds.features(wind,1)),'.r')
% $$$ plot(vel(ds.insp_periods(rind,1)),log10(ds.features(rind,1)),'.g')
plot(vel(ds.insp_periods(hind,1)),ds.features(hind,3),'.r')
plot(vel(ds.insp_periods(lind,1)),ds.features(lind,3),'.g')

figure
hist2([vel(ds.insp_periods(nind,1)),ds.features(nind,2)],100,100)


rhm = fet_rhm(Trial,[],'default');

lfp = Trial.lfp.copy;
lfp.load(Trial,66);
lfp.resample(xyz);

figure,
sp(1) = subplot(211);
plot(rhm.data);
sp(2) = subplot(212);
plot(lfp.data);
linkaxes(sp,'x');




%% NEW Session JM11-20161126
SessionName = 'JM11-20161126'; % what was once known as filebase
xyz_path = '/storage/javier/data/processed/xyz';
MazeName = 'cof';
TrialName = 'all'; 
overwrite = true;
xyzSampleRate = 199.997752;

linkSession(SessionName,xyz_path,[]);
s = MTASession(SessionName,MazeName,overwrite,[],'vicon','box',xyzSampleRate);
xyz = s.load('xyz');
xyz.data(:,:,1) = xyz(:,:,1)-85;
xyz.data(:,:,2) = xyz(:,:,2)-95;
xyz.save;

s = MTASession(SessionName);

pZ(s);
pXY(s);
PlotSessionErrors(s);

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
vxy = fxyz.vel([1],[1,2]);

[ys,fs,ts] = fet_rhm(s,xyz,'mtchglong',true);

figure,
sp =[];
sp(1)=subplot(211);
imagesc(ts,fs,log10(ys.data)');
axis xy
caxis([-6,-4])
sp(2)=subplot(212)
plot((1:size(vxy,1))./vxy.sampleRate,vxy.data)
linkaxes(sp,'x')