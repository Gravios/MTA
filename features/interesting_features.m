Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,10,'low');
ang = create(MTADang,Trial,xyz);



figure,hold on
plot(circ_dist(ang(:,1,3,1),ang(:,2,4,1)))
plot(circ_dist(ang(:,2,4,1),ang(:,3,7,1)))
Lines(Trial.stc{'w'}(:),[],'b');
Lines(Trial.stc{'n'}(:)-1,[],'g');


name = 'spine and head wag'; label = 'shwag'; key = 'w';
zv = MTADfet.encapsulate(Trial,...
                         [circ_dist(ang(:,1,3,1),ang(:,2,4,1)),...
                          circ_dist(ang(:,2,4,1),ang(:,3,7,1))],...
                         xyz.sampleRate,...
                         name,label,key);


dspec = struct('nFFT',2^8,'Fs',zv.sampleRate,...
               'WinLength',2^7,'nOverlap',2^7*.875,...
               'FreqRange',[1,30]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',true,[],dspec,true);


figure,imagesc(ts,fs,ys(:,:,1,2)'),axis xy,colormap jet

figure,imagesc(ts,fs,log10(ys(:,:,2,2))'),axis xy,colormap jet


eds = linspace(-9,1,100);
figure,hold on
ind = Trial.stc{'a-p'};
hs = bar(eds,histc(log10(ys(ind,2,2,2)),eds),'histc');hs.FaceColor = 'c';hs.FaceAlpha = .5;
ind = Trial.stc{'p'};
hs = bar(eds,histc(log10(ys(ind,2,2,2)),eds),'histc');hs.FaceColor = 'r';hs.FaceAlpha = .5;

name = 'zv and wag pw'; label = 'zvwp'; key = 'w';
zfet = MTADfet.encapsulate(Trial,...
                         log10(nanmedian(ys(:,1:fs<7,1,1),2)),...
                         ys.sampleRate,...
                         name,label,key);
zfet.data(~nniz(zfet),:)=-20;
zfet.resample(bfet);



%exp in normalizing distances
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');

xyz = Trial.load('xyz');
xyz.addMarker('bcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
ang = create(MTADang,Trial,xyz);
NBINS = 100;
VEL_HISTOGRAM_BOUNDARIES = linspace(-3,2,NBINS);
vxy = xyz.copy;
vxy.filter('ButFilter',3,2.4,'low');
vxy = vxy.vel({'bcom','hcom'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);
ind = Trial.stc{'a-m'}.cast('TimeSeries');
[~,ind_v_b] = histc(vxy(:,'bcom'),VEL_HISTOGRAM_BOUNDARIES);
[~,ind_v_h] = histc(vxy(:,'hcom'),VEL_HISTOGRAM_BOUNDARIES);
mind = nniz([ind_v_b,ind_v_h])&ang(:,3,4,2)<.2&ind.data;
mz = accumarray([ind_v_b(mind),ind_v_h(mind)],...
                ang(mind,4,7,2),...
                [NBINS,NBINS],...
                @nanmean);
vz = accumarray([ind_v_b(mind),ind_v_h(mind)],...
                ang(mind,4,7,2),...
                [NBINS,NBINS],...
                @nanstd);


% $$$ figure,subplot(1,2,1),imagesc(mz'),subplot(1,2,2),imagesc(vz')
% $$$ figure,hist(mz(nniz(mz(:))&vz(:)<10),100)


Trial = MTATrial('Ed03-20140624');
Trial.load('stc','hand_labeled_rev2_alt');

xyz = Trial.load('xyz');
xyz.addMarker('bcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
ang = create(MTADang,Trial,xyz);
NBINS = 100;
VEL_HISTOGRAM_BOUNDARIES = linspace(-3,2,NBINS);
vxy = xyz.copy;
vxy.filter('ButFilter',3,2.4,'low');
vxy = vxy.vel({'bcom','hcom'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);
ind = Trial.stc{'a-m'}.cast('TimeSeries');
[~,ind_v_b] = histc(vxy(:,'bcom'),VEL_HISTOGRAM_BOUNDARIES);
[~,ind_v_h] = histc(vxy(:,'hcom'),VEL_HISTOGRAM_BOUNDARIES);
mind = nniz([ind_v_b,ind_v_h])&ang(:,3,4,2)<.2&ind.data;
mzo = accumarray([ind_v_b(mind),ind_v_h(mind)],...
                ang(mind,4,7,2),...
                [NBINS,NBINS],...
                @nanmean);
vzo = accumarray([ind_v_b(mind),ind_v_h(mind)],...
                ang(mind,4,7,2),...
                [NBINS,NBINS],...
                @nanstd);


figure,subplot(1,2,1),imagesc(mz'),axis xy,subplot(1,2,2), ...
    imagesc(mzo'),axis xy
ind = vz(:)<.2&vzo(:)<.2&nniz(mzo(:))&nniz(mz(:));
figure,hist(mz(ind)-mzo(ind),100)




eds = linspace(-pi/2,pi/2,100);
figure,hold on 

Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');

xyz = Trial.load('xyz');
xyz.addMarker('bcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
ang = create(MTADang,Trial,xyz);

ind = Trial.stc{'w'};
hs = bar(eds,histc(ang(ind,4,7,2),eds),'histc'); hs.FaceColor = 'c'; hs.FaceAlpha = .5;
Trial = MTATrial('Ed03-20140624');
Trial.load('stc','hand_labeled_rev2_alt');

xyz = Trial.load('xyz');
xyz.addMarker('bcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
ang = create(MTADang,Trial,xyz);

ind = Trial.stc{'w'};
hs = bar(eds,histc(ang(ind,4,7,2),eds),'histc'); hs.FaceColor = 'r'; hs.FaceAlpha = .5;















figure,plot()
figure,plot(ang(:,'bcom','hcom',2))


fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1,'low');
fang = create(MTADang,Trial,fxyz);

name = 'head and body com'; label = 'hbcom'; key = 'c';
zv = MTADfet.encapsulate(Trial,...
                         [ang(:,'bcom','hcom',2)-fang(:,'bcom','hcom',2)],...
                         xyz.sampleRate,...
                         name,label,key);


dspec = struct('nFFT',2^9,'Fs',zv.sampleRate,...
               'WinLength',2^8,'nOverlap',2^8*.875,...
               'FreqRange',[1,30]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',true,[],dspec,true);


figure,imagesc(ts,fs,ys(:,:,1,2)'),axis xy,colormap jet

figure,imagesc(ts,fs,log10(ys(:,:,1,1))'),axis xy,colormap jet
figure,imagesc(ts,fs,log10(ys(:,:,2,2))'),axis xy,colormap jet

eds = linspace(-.2,1.5,100);
figure,hold on
ind = Trial.stc{'s'};
hs = bar(eds,histc(ang(ind,1,3,2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'p'};
hs = bar(eds,histc(ang(ind,1,3,2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='y';
ind = Trial.stc{'w'};
hs = bar(eds,histc(ang(ind,1,3,2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='b';
ind = Trial.stc{'m'};
hs = bar(eds,histc(ang(ind,1,3,2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='m';
ind = Trial.stc{'n'};
hs = bar(eds,histc(ang(ind,1,3,2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='g';
ind = Trial.stc{'r'};
hs = bar(eds,histc(ang(ind,1,3,2),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';




tfet = MTADxyz('data',bsxfun(@rdivide,ang_b_vel,10.^fang(:,3,4,2)),'sampleRate',fang.sampleRate);
tfet.filter('ButFilter',3,2,'low');
tfet.data = mean(tfet.data,2);
tfet.data(tfet.data<1e-6) = 1e-6;
tfet.data = log10(tfet.data);

% $$$ figure,plot(tfet);Lines(Trial.stc{'w'}(:),[],'b');Lines(Trial.stc{'n'}(:)+1,[],'g');



eds = linspace(-6,1,100);
figure,hold on
ind = Trial.stc{'s'};
hs = bar(eds,histc(tfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'p'};
hs = bar(eds,histc(tfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='y';
ind = Trial.stc{'w'};
hs = bar(eds,histc(tfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='b';
ind = Trial.stc{'m'};
hs = bar(eds,histc(tfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='m';
ind = Trial.stc{'n'};
hs = bar(eds,histc(tfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='g';
ind = Trial.stc{'r'};
hs = bar(eds,histc(tfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';



Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
figure,plot(diff(fxyz(:,1,3)));Lines(Trial.stc{'w'}(:),[],'b');Lines(Trial.stc{'n'}(:)+1,[],'g');




fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1.5,'low');

bfet = Trial.xyz.copy;
tsh = 20;
afet = Trial.xyz.copy;
afet.data = circshift(fxyz.data,-tsh)-circshift(fxyz.data,tsh);
afet.data = reshape(afet.data,[],3);
afet.data = permute(bsxfun(@dot,permute(reshape(repmat(fxyz(:,4,:)-fxyz(:,1,:),[1,fxyz.size(2),1]),[],3),[2,1]),permute(afet.data,[2,1])),[2,1]);
afet.data = reshape(afet.data,[],fxyz.size(2));
bfet.data = mean(afet(:,1:4),2)./1000.*log10(var(afet(:,1:4),[],2));
bfet.filter('ButFilter',3,2.4,'low');

lbfet = bfet.copy;

lbfet.data = abs(log10(bfet.data(:,1)+10))

figure,plot(lbfet.data)
Lines(Trial.stc{'w'}(:),[],'b');Lines(Trial.stc{'n'}(:)+1,[],'g');


eds = linspace(0,3,100);
figure,hold on
ind = Trial.stc{'s'};
hs = bar(eds,histc(lbfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';
ind = Trial.stc{'p'};
hs = bar(eds,histc(lbfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='y';
ind = Trial.stc{'w'};
hs = bar(eds,histc(lbfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='b';
ind = Trial.stc{'m'};
hs = bar(eds,histc(lbfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='m';
ind = Trial.stc{'n'};
hs = bar(eds,histc(lbfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='g';
ind = Trial.stc{'r'};
hs = bar(eds,histc(lbfet(ind),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';

mwp = [];
for per = Trial.stc{'w'}.data'
    mwp(end+ =maxlbfet(per(1):per(