% Find head only features of head scratching 

Trial = MTATrial('jg05-20120317');

xyz = Trial.load('xyz');

ang = create(MTADang,Trial,xyz);

hang = Trial.transform_origin(xyz,'head_back','head_front',{'head_left','head_right'});
name = 'head roll'; label = 'hroll'; key = 'r';
zv = MTADfet.encapsulate(Trial,...
                         hang.roll,...
                         xyz.sampleRate,...
                         name,label,key);
zv.filter('ButFilter',5,50);
zv.data = diff(zv.data);
                     
dspec = struct('nFFT',2^8,'Fs',zv.sampleRate,...
               'WinLength',2^7,'nOverlap',2^7*.875,...
               'FreqRange',[1,30]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',true,[],dspec,true);
figure,imagesc(ts,fs,log10(ys.data)'),axis xy, colormap jet
Lines(Trial.stc{'m',1}(:),[],'m');


edx = linspace(60,140,100);
edy = linspace(-pi/2,pi/2,100);
figure
ind = Trial.stc{'m'};
hist2([ang(ind,1,4,3).*cos(ang(ind,1,4,2)),hroll(ind)],edx,edy);


edx = linspace(-pi/2,pi/2,100);
edy = linspace(-pi/2,pi/2,100);
figure
ind = Trial.stc{'m'};
hist2([ang(ind,5,7,2),hroll(ind)],edx,edy);


figure,plot(hroll(:))
Lines(Trial.stc{'m'},[],'m');

[S,U,V] = svd(cov(ys(Trial.stc{'a'},:)));


figure,plot(log10(mean(ys(:,fs>10&fs<16),2)./mean(ys(:,(fs>4&fs<8)|(fs>16&fs<20)),2)))
Lines(Trial.stc{'m',ys.sampleRate}(:),[],'m');

spow = ys.copy;
spow.data = log10(mean(ys(:,fs>10&fs<16),2)./mean(ys(:,(fs>4&fs<8)|(fs>16&fs<20)),2));
spow.resample(xyz);

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1);
fvxy = fxyz.vel(5,[1,2]);
fvxy.data(fvxy.data<1e-3)=1e-3;
fvxy.data = log10(fvxy.data);

edx = linspace(-2.9,2,100);
edy = linspace(-2,2,100);
figure
ind = Trial.stc{'m'};
hist2([fvxy(ind),spow(ind)],edx,edy);

edx = linspace(-pi/2,pi/2,100);
edy = linspace(-2,2,100);
figure
ind = Trial.stc{'m'};
hist2([hroll(ind),spow(ind)],edx,edy);


fang = create(MTADang,Trial,fxyz);



name = 'head roll'; label = 'hroll'; key = 'r';
zv = MTADfet.encapsulate(Trial,...
                         diff(fang(:,5,7,2)),...
                         xyz.sampleRate,...
                         name,label,key);
zv.filter('ButFilter',5,50);
zv.data = diff(zv.data);
                     
dspec = struct('nFFT',2^8,'Fs',zv.sampleRate,...
               'WinLength',2^7,'nOverlap',2^7*.875,...
               'FreqRange',[1,30]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',true,[],dspec,true);



figure,plot(fang(:,5,7,2))
Lines(Trial.stc{'m'}(:),[],'m');




