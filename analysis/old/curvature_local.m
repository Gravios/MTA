
Trial = MTATrial('jg05-20120317');

fxyz = Trial.load('xyz');
fxyz.filter('ButFilter',3,1.5);
fang = create(MTADang,Trial,fxyz);



xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,50);

xyz.data = bsxfun(@minus,xyz.data,xyz(:,5,:));
dx = [zeros([1,size(xyz,[2,3])]);diff(xyz.data)];
ddx = [diff(dx);zeros([1,size(xyz,[2,3])])];

kappa = sqrt((ddx(:,:,3).*dx(:,:,2)-ddx(:,:,2).*dx(:,:,3)).^2+...
             (ddx(:,:,1).*dx(:,:,3)-ddx(:,:,3).*dx(:,:,1)).^2+...
             (ddx(:,:,2).*dx(:,:,1)-ddx(:,:,1).*dx(:,:,2)).^2)./...
        ((dx(:,:,1).^2+dx(:,:,2).^2+dx(:,:,3).^2).^(3/2));



ka = sq(nanmean(GetSegs(log10(1+1./kappa),1:length(kappa),60,nan)));
%kv = sq(nanvar(GetSegs(log10(1+1./kappa),1:length(kappa),60,nan)));

figure,hold on
plot(circshift(ka(:,[1,2,4,5]),30))
%plot(circshift(kv(:,1),40))
Lines(Trial.stc{'w'}(:),[],'c');
Lines(Trial.stc{'r'}(:),[],'r');
Lines(Trial.stc{'n'}(:),[],'g');


kat = Trial.xyz.copy;
kat.data = log10(circshift(ka(:,1),30));

eds = linspace(-3,1,100);
figure,hold on
ind = Trial.stc{'a'};
hs = bar(eds,histc(kat(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'k';
ind = Trial.stc{'r'};
hs = bar(eds,histc(kat(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'm';
ind = Trial.stc{'p'};
hs = bar(eds,histc(kat(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'y';
ind = Trial.stc{'w'};
hs = bar(eds,histc(kat(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{'n'};
hs = bar(eds,histc(kat(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';



figure,plot(kk(:,1:3))
Lines(Trial.stc{'w'}(:),[],'c');

figure,plot(ka(:,1))
Lines(Trial.stc{'w'}(:),[],'c');

figure,plot(kv)


figure,plot(nanmean(ka(:,1:5),2)./(1+nanvar(ka(:,1:5),[],2)));
Lines(Trial.stc{'w'}(:),[],'c');


fv = fxyz.vel(1,[1,2]);
fv.data(fv.data<1e-3) = 1e-3;
fv.data = log10(fv.data);

%za = fv.copy;
at = fv.copy;
za = zav.copy;
edx = linspace(-2.9,2,100);
%edx = linspace(-1.5,0.5,100);
edy = linspace(-5,-1,100);
%edy = linspace(-2.9,2,100);
figure,k = 1;
subplot(2,2,k);k = k+1;
ind = Trial.stc{'a'};
hist2([at(ind),za(ind)],edx,edy);
subplot(2,2,k);k = k+1;
ind = Trial.stc{'n'};
hist2([at(ind),za(ind)],edx,edy);
subplot(2,2,k);k = k+1;
ind = Trial.stc{'p'};
hist2([at(ind),za(ind)],edx,edy);
subplot(2,2,k);k = k+1;
ind = Trial.stc{'w'};
hist2([at(ind),za(ind)],edx,edy);
ForAllSubplots('caxis([0,100]),grid on')




xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,50);


xs = circshift(xyz.segs(1:xyz.size(1),60,nan),30,2);
xs = bsxfun(@minus,xs,xs(:,:,5,:));
ad = sqrt(sum(sq(xs(end,:,:,:)-xs(1,:,:,:)).^2,3));
sd = sq(sum(sqrt(sum(diff(xs).^2,4))));

figure,plot(ad(:,7)./sd(:,7))
Lines(Trial.stc{'m'}(:),[],'m');
Lines(Trial.stc{'w'}(:),[],'c');
Lines(Trial.stc{'n'}(:),[],'g');

figure,plot(ad(:,7)./sd(:,7)-ad(:,5)./sd(:,5))

figure,plot(diff(ad(:,7)./sd(:,7)))

name = 'head traj sinuosity'; label = 'hts'; key = 'h';
zv = MTADfet.encapsulate(Trial,...
                         ad(:,7)./sd(:,7),...
                         xyz.sampleRate,...
                         name,label,key);


dspec = struct('nFFT',2^10,'Fs',zv.sampleRate,...
               'WinLength',2^9,'nOverlap',2^9*.875,...
               'FreqRange',[.1,5]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',true,[],dspec,true);

figure,imagesc(ts,fs,log10(ys(:,:))'),axis xy,colormap jet

figure,plot(ys(:,10))
Lines(Trial.stc{'w',ys.sampleRate}(:),[],'c');

ufr = Trial.ufr.copy;
ufr.create(Trial,ys,'a',1:90,1);

ysf = ys.copy;
[~,ysf.data] = max(ys.data,[],2);
ysf.data = fs(ysf.data);

figure
for unit = 1:90;
    plot(log10(ys(:,20))+randn([ufr.size(1),1])/2,ufr(:,unit)+ randn([ufr.size(1),1])/2,'.')
    %plot(ysf.data+randn([ufr.size(1),1])/2,ufr(:,unit)+ randn([ufr.size(1),1])/2,'.')
pause(.2)
end




[rhm,fs] = fet_rhm(Trial,[],'mtchglong',true);
rhm.data = median(rhm(:,fs>6&fs<14),2);
% Add to feature matrix
rhm.data(rhm.data<1e-8)=1e-8;

ufr = Trial.ufr.copy;
ufr.create(Trial,rhm,'a',1:90,1);

ysf = ys.copy;
[~,ysf.data] = max(ys.data,[],2);
ysf.data = fs(ysf.data);

figure
noise = ufr.copy;
noise.data = randn([ufr.size(1),1])/2;
ind = Trial.stc{'a'};
for unit = 1:90;
    plot(log10(rhm(ind))+noise(ind),ufr(ind,unit)+ noise(ind),'.')
    %plot(ysf.data+randn([ufr.size(1),1])/2,ufr(:,unit)+ randn([ufr.size(1),1])/2,'.')
pause(.8)
end

