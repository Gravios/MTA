
Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');



fang = create(MTADang,Trial,fxyz);



[rhm,] = fet_spec(Trial,fet,varargin);

[mode,wsig,sampleRate,defspec,overwrite] = ...
     DefaultArgs(varargin,{'raw',true,120,...
                     struct('nFFT',2^9,'Fs',fet.sampleRate,...
                            'WinLength',2^7,'nOverlap',2^7*.875,...
                            'FreqRange',[1,40]),...
                     true});

                 
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,.1,'high');                 
figure,plot(diff(fxyaz(:,[1,2,3,4],3)))
Lines(Trial.stc{'n'}(:),[],'g');Lines(Trial.stc{'w'}(:),[],'c');

fvfz = xyz.copy;
fvfz.data = diff(fxyz(:,[1,2,3,4],3));
fvfz.filter('ButFilter',3,2,'low');                 
figure,plot(fvfz.data)
Lines(Trial.stc{'n'}(:),[],'g');Lines(Trial.stc{'w'}(:),[],'c');

figure,plot(sqrt(sum([circshift(xyz(:,1,[1,2]),-30)-circshift(xyz(:,1,[1,2]),30)].^2,3)))
Lines(Trial.stc{'n'}(:),[],'g');Lines(Trial.stc{'w'}(:),[],'c');

xyz = Trial.load('xyz');
vfx = xyz.copy;
vfx.filter('ButFilter',3,2.4,'low');
vfx = xyz.vel(1,[1,2]);

shift = 30;
for shift = 5:5:80;
vcs = xyz.copy;
vcs.data = sqrt(sum([circshift(xyz(:,1,[1,2]),-shift)-circshift(xyz(:,1,[1,2]),shift)].^2,3));
end


%figure,plot(log10(vcs.data))
%figure,plot(log10(vfx.data))
figure,
subplot(221)
ind = Trial.stc{'a-w-n'};
hist2(log10([vfx(ind),vcs(ind)]),linspace(-1,2,100),linspace(-3,2.5,100))
caxis([0,200])
subplot(222)
ind = Trial.stc{'w'};
hist2(log10([vfx(ind),vcs(ind)]),linspace(-1,2,100),linspace(-3,2.5,100))
caxis([0,200])
subplot(223)
ind = Trial.stc{'a'};
hist2(log10([vfx(ind),vcs(ind)]),linspace(-1,2,100),linspace(-3,2.5,100))
caxis([0,200])
subplot(224)
ind = Trial.stc{'n'};
hist2(log10([vfx(ind),vcs(ind)]),linspace(-1,2,100),linspace(-3,2.5,100))
caxis([0,200])



hvar = vcs;
eds = linspace(-3,3,100);
figure,hold on
ind = Trial.stc{'a-w-n'};
hs = bar(eds,histc(log10(hvar(ind)),eds),'histc');
hs.FaceColor = 'c';
hs.FaceAlpha = .4;

ind = Trial.stc{'w'};
hs = bar(eds,histc(log10(hvar(ind)),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;


