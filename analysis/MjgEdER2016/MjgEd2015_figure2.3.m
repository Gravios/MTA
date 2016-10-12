%Spline Iterpolation of spine

Trial = MTATrial('Ed03-20140624');
xyz = Trial.load('xyz');
ssp = fet_spline_spine(Trial);
sxyz = Trial.load('xyz','seh');

hfig = figure(2016061403),clf
ind = 2000;
plot3(sxyz(ind,:,1),sxyz(ind,:,2),sxyz(ind,:,3),'.m')
set(get(gca,'Children'),'MarkerSize',20)
plotSkeleton(Trial,xyz,ind);
hold on
plot3(ssp(ind,:,1),ssp(ind,:,2),ssp(ind,:,3),'r')

view([-117,18])
       
print(hfig,'-depsc2','-loose',...
      fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_3',...
               ['fig3_A_',num2str(s),'_sk_',slist{s},'.eps']))



sfxyz = sxyz.copy;
sfxyz.filter('ButFilter',3,.5,'low');

rb = sxyz.model.rb({'head_back','head_left','head_front','head_right'});
% find the center of mass of the model
hcom = sfxyz.com(rb);
% add coordinates of the model's center of mass to the xyz object
sxyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);

sang = create(MTADang,Trial,sxyz);




fxyz = sxyz.copy;
fxyz.filter('ButFilter',3,2,'low');

rb = xyz.model.rb({'head_back','head_left','head_front','head_right'});
% find the center of mass of the model
hcom = fxyz.com(rb);
% add coordinates of the model's center of mass to the xyz object
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);

ang = create(MTADang,Trial,xyz);



ncp = Trial.load('lfp',2);
ncp.resample(sxyz);

fssp = ssp.copy;
fssp.filter('ButFilter',3,55,'low');
dssp = MTADxyz('data',sum(sqrt(sum([fssp.data-circshift(fssp.data,-1,2)].^2,3)),2),'sampleRate',ssp.sampleRate);
dssp.data = [0;diff(dssp.data)];
dssp.filter('ButFilter',3,[2,20],'bandpass');
dssp.data = [0;diff(dssp.data)];
dssp.filter('ButFilter',3,[2,20],'bandpass');
dssp.name = 'rythmic body motion';


srhm = MTADxyz('data',sang(:,5,11,3),'sampleRate',sang.sampleRate);
srhm.filter('ButFilter',3,[2,50],'bandpass');
srhm.data = [0;diff(srhm.data)];
srhm.name = 'rhm trb-spline corrected';

spt = MTADxyz('data',sang(:,5,7,2),'sampleRate',sang.sampleRate);
spt.data = [0;diff(spt.data)];
spt.filter('ButFilter',3,[2,50],'bandpass');
spt.name = 'rythmic head pitch';


rhm = MTADxyz('data',ang(:,5,10,3),'sampleRate',ang.sampleRate);
rhm.filter('ButFilter',3,[2,20],'bandpass');
rhm.data = [0;diff(rhm.data)];


% $$$ figure,clf,hold on,
% $$$ plot(dssp.data)
% $$$ plot(ncp.data/5000)


x = dssp.copy;
x.data = [ncp.data/5000,x.data,srhm.data,rhm.data];

[ys,fs,ts] = fet_spec(Trial,x,'mtchglong',1);

sp = [];
n = 4;
figure,
sp(end+1)=subplot(n,1,1);
imagesc(ts,fs,log10(ys(:,:,1,1))');
colormap jet
caxis([-5,-3.2])
axis xy
for j = 2:n,
sp(end+1)=subplot(n,1,j);
imagesc(ts,fs,ys(:,:,1,j)');
colormap jet
caxis([0.4,.8])
axis xy
end
linkaxes(sp,'xy');



sp = [];
figure,
for j = 1:4,
sp(end+1)=subplot(4,1,j);
imagesc(ts,fs,log10(ys(:,:,j,j))');
colormap jet
caxis([-5,-3])
axis xy
end
linkaxes(sp,'xy');

