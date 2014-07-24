
Trial = MTATrial('jg05-20120311');

fwin = gtwin(1.25,Trial.xyz.sampleRate);
swin = gtwin(.1,Trial.xyz.sampleRate);

xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gausswin(5)./sum(gausswin(5)));

rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
rbb = Trial.xyz.model.rb({'spine_lower','pelvis_root'});

hcom = xyz.com(rb);
xyz.addMarker('hcom',[.7,0,.7],{{'head_back','head_front',[0,0,1]}},hcom);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0, ...
                    1]}},permute(Filter0(fwin,hcom),[1,3,2]));
xyz.addMarker('fbcom',[.7,1,.7],{{'head_back','head_front',[0,0, ...
                    1]}},ButFilter(xyz(:,7,:),3,[2]./(Trial.ang.sampleRate/2),'low'));
xyz.addMarker('bbcom',[.7,1,.7],{{'head_back','head_front',[0,0, ...
                    1]}},ButFilter(xyz(:,5,:),3,[2]./(Trial.ang.sampleRate/2),'low'));
xyz.addMarker('chcom',[.7,1,.7],{{'head_back','head_front',[0,0, ...
                    1]}},ButFilter(hcom,3,[2]./(Trial.ang.sampleRate/2),'low'));






%% Stupid geometry stuff
xyz.data = xyz.data-repmat(xyz(:,'head_back',:),[1,xyz.size(2),1]);



%fline = xyz(:,{'bbcom','fbcom'},:);
fxy = xyz(:,'head_front',:);
rxy = xyz(:,'head_right',:);
lxy = xyz(:,'head_left',:);
fnorm = cross(xyz(:,{'bbcom'},:),xyz(:,{'fbcom'},:));

rnx = rxy.*repmat(dot(fnorm,rxy,3)./sum(rxy.^2,3),[1,1,3]);
lnx = lxy.*repmat(dot(fnorm,lxy,3)./sum(lxy.^2,3),[1,1,3]);
fnx = fxy.*repmat(dot(fnorm,fxy,3)./sum(fxy.^2,3),[1,1,3]);

fnpr = rnx.*fnx;
fnpl = lnx.*fnx;


angR = atan2(sqrt(sum(cross(fnpr,fxy).^2,3)), dot(fnpr,fxy,3));
angL = atan2(sqrt(sum(cross(fnpl,fxy).^2,3)), dot(fnpl,fxy,3));
repmat(dot(fnp(:,1,:),fxy,3)./sum(fxy.^2,3),[1,1,3]);


rxy.*(dot(fnorm.*rxy))./sum(rxy.^2,3))

%B*(A*B/|B|)

rplane = xyz(:,{'head_front','head_right'},:);

(sqrt(sum(sq(rplane(1,1,:)).^2))*sqrt(sum(pp.^2)))

pp = (sq(rplane(1,1,:))*sq(rplane(1,2,:))')*sq(fnorm(1,:,:));
dot(sq(rplane(1,1,:)),pp)/(sqrt(sum(sq(rplane(1,1,:)).^2))*sqrt(sum(pp.^2)))

%% Newest vector stuff
% $$$ xyz.addMarker('fbcom',[.7,1,.7],{{'head_back','head_front',[0,0, ...
% $$$                     1]}},ButFilter(xyz(:,7,:),3,[2]./(Trial.ang.sampleRate/2),'low'));
% $$$ xyz.data = xyz.data-repmat(xyz(:,'head_back',:),[1,xyz.size(2),1]);

fd = (circshift(xyz(:,'fbcom',:),-1)-circshift(xyz(:,'fbcom',:),1))/2;
%fd = (circshift(xyz(:,'head_front',:),-1)-circshift(xyz(:,'head_front',:),1))/2;
td = fd+xyz(:,'head_front',:);
[~,~,sd] = cart2sph(td(:,1,1),td(:,1,2),td(:,1,3));

[ys,fs,ts,phi,fst] = mtchglong(WhitenSignal([sd,hrol],[],1),2^8,Trial.ang.sampleRate,2^7,2^6,[],'linear',[],[1,40]);
ts = ts+diff(ts(1:2));

figure
for i = 1:4
sp(i) = subplot(4,1,i); 
%imagesc(ts,fs,log10(ys(:,:,i,i)'));axis xy
imagesc(ts,fs,(fsta(:,:,i)'));axis xy
end
linkaxes(sp,'xy')


figure
for i = 1:4
sp(i) = subplot(4,1,i); 
imagesc(ts,fs,log10(ys(:,:,i,i)'));axis xy
end
linkaxes(sp,'xy')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555555555555
Trial = MTATrial('Ed05-20140528');
%Trial = MTATrial('jg05-20120317');
%stc_mode = 'qda_filtf1p5';
stc_mode = 'auto_wbhr';
Trial.stc.updateMode(stc_mode);Trial.stc.load;
xyz = Trial.xyz.copy;
xyz.load(Trial);

%xyz.filter(gtwin(.1,xyz.sampleRate));
rb = Trial.xyz.model.rb({'head_back','head_left','head_right'});
hcom = xyz.com(rb);

xyz.addMarker('fbcom',[.7,1,.7],{{'head_back','head_front',[0,0, ...
                    1]}},ButFilter(xyz(:,7,:),3,[2]./(Trial.ang.sampleRate/2),'low'));
xyz.addMarker('bbcom',[.7,1,.7],{{'head_back','head_front',[0,0, ...
                    1]}},ButFilter(xyz(:,5,:),3,[2]./(Trial.ang.sampleRate/2),'low'));
xyz.addMarker('chcom',[.7,1,.7],{{'head_back','head_front',[0,0, ...
                    1]}},ButFilter(hcom,3,[2]./(Trial.ang.sampleRate/2),'low'));


%xyz.data = xyz.data-repmat(xyz(:,'bbcom',:),[1,xyz.size(2),1]);
txyz = xyz.copy;
txyz.data = txyz.data-repmat(txyz(:,'head_back',:),[1,txyz.size(2),1]);
fd = (circshift(txyz(:,'chcom',:),-1)-circshift(txyz(:,'chcom',:),1))/2;
td = fd+txyz(:,'head_front',:);
[~,~,sd] = cart2sph(td(:,1,1),td(:,1,2),td(:,1,3));
sd = Filter0(gausswin(7)./sum(gausswin(7)),sd);

sd = [sq(xyz(:,'head_front',:)),sd];

ang = Trial.ang.copy;
ang.create(Trial,Trial.xyz);
pang = sq(ButFilter(ang(:,5,7,2),3,[1,20]./(Trial.ang.sampleRate/2),'bandpass'));

%rhm = fet_rhm(Trial);
%ncp = fet_ncp(Trial,'chans',[1,2]);



wang = [sd,rhm,pang,ncp(:,2)];
wang = WhitenSignal(wang,[],1);

[ys,fs,ts,phi,fsta] = mtchglong(wang,2^9,Trial.ang.sampleRate,2^7,2^7*.875,3,'linear',[],[1,20]);


figure
for i = 1:7
sp(i) = subplot(7,1,i); 
imagesc(ts,fs,log10(ys(:,:,i,i)'));axis xy
end
linkaxes(sp,'xy')



figure
sp(1) = subplot(6,1,1); 
imagesc(ts,fs,log10(ys(:,:,1,1)'));axis xy,caxis([-7,-3.4])
sp(2) = subplot(6,1,2); 
imagesc(ts,fs,nys(:,:,1,1)');      axis xy,caxis([-.5,3])
sp(3) = subplot(6,1,3); 
imagesc(ts,fs,log10(ys(:,:,2,2)'));axis xy,caxis([-7,-3])
sp(4) = subplot(6,1,4); 
imagesc(ts,fs,nys(:,:,2,2)');axis xy,caxis([-.50,3])
sp(5) = subplot(6,1,5); 
imagesc(ts,fs,log10(ys(:,:,4,4)'));axis xy,caxis([0,4.75])
sp(6) = subplot(6,1,6); 
imagesc(ts,fs,nys(:,:,4,4)');      axis xy,caxis([0,3])
linkaxes(sp,'xy')





caxis([-8,-4.2])

ang = Trial.ang.copy;
ang.create(Trial,xyz);

bang = ButFilter(ang(:,'head_back','fhcom',3),3,[1,30]./(Trial.ang.sampleRate/2),'bandpass');
nang = ButFilter(ang(:,'head_back','chcom',3),3,[1,30]./(Trial.ang.sampleRate/2),'bandpass');
fang = ButFilter(ang(:,'head_front','fhcom',3),3,[1,30]./(Trial.ang.sampleRate/2),'bandpass');
rang = ButFilter(ang(:,'head_front','chcom',3),3,[1,30]./(Trial.ang.sampleRate/2),'bandpass');
lang = ButFilter(ang(:,'head_left','fhcom',3),3,[1,30]./(Trial.ang.sampleRate/2),'bandpass');
pang = ButFilter(ang(:,'head_back','bhcom',3),3,[1,30]./(Trial.ang.sampleRate/2),'bandpass');
sang = ButFilter(ang(:,'spine_upper','head_back',3),3,[1,30]./(Trial.ang.sampleRate/2),'bandpass');

figure,hold on,plot(bang),plot(nang,'r')

hang = Trial.transformOrigin;
hrol = hang.roll;
hind = nniz(hrol);
hrol(hind) = ButFilter(hrol(hind),3,[1,30]./(Trial.ang.sampleRate/2),'bandpass');

ind = nniz(nang)&nniz(rang)&nniz(lang)&nniz(hang.roll)&nniz(sang);
wang = nan(size(hrol),5);
wang(ind,:) = WhitenSignal([nang(ind),rang(ind),lang(ind),hrol(ind),sang(ind)]);


[ys,fs,ts,phi,fst] = mtchglong(wang,2^8,Trial.ang.sampleRate,2^7,2^6,[],'linear',[],[1,30]);
ts = ts+diff(ts(1:2));



figure
for i = 1:5
sp(i) = subplot(6,1,i);
%imagesc(ts,fs,nys(:,:,i,i)');
imagesc(ts,fs,log10(ys(:,:,i,i)'));
Lines([],10,'k');
%caxis([-2,2])
caxis([-6,-2])
axis xy
end
sp(6) = subplot(6,1,6);
imagesc(ts,fs,ys(:,:,1,5)'+ys(:,:,4,5)'+ys(:,:,4,3)'),axis xy
Lines(Trial.stc{'w',1}(:),[],'w');
caxis([1.4,3])
linkaxes(sp,'xy')


yst = log10(ys(:,:,3,3));
yst(yst<-8) = nan;
yst(~nniz(yst),:) = nan;

slims = prctile(yst,[2,99]);
smi = nan(size(yst,2),size(yst,2));
for i = 1:size(yst,2),
for j = 1:size(yst,2),
    v1e = linspace(slims(1,i),slims(2,i),16);
    v2e = linspace(slims(1,j),slims(2,j),16);
    v1m = histc(yst(:,i),v1e);
    v2m = histc(yst(:,j),v2e);
    v1m = v1m(1:end-1)'/sum(v1m);
    v2m = v2m(1:end-1)'/sum(v2m);
    v12jpdf = hist2(yst(nniz(yst),[i,j]),v1e,v2e);
    v12jpdf = v12jpdf/nansum(v12jpdf(:));
    smi(i,j) = nansum(nansum(v12jpdf.*log(v12jpdf./(v1m'*v2m))));
end    
end

 

lys = ys./repmat(nanmean(ys,2),[1,size(ys,2),1,1]);
lys = log10(lys);
%figure,imagesc(ts,fs,lys(:,:,2,2)'),axis xy
nys = (lys-repmat(nanmean(lys(nniz(lys),:,:,:)),[size(lys,1),1,1,1]))./repmat(nanstd(lys(nniz(lys),:,:,:)),[size(lys,1),1,1,1]);

%lys = ys./repmat(mean(ys,2),[1,size(ys,2),1,1]);
%lys = log10(lys);
%figure,imagesc(ts,fs,lys(:,:,2,2)'),axis xy
%nys = (ys-repmat(nanmean(ys(nniz(ys),:,:,:)),[size(ys,1),1,1,1]))./repmat(nanstd(ys(nniz(ys),:,:,:)),[size(ys,1),1,1,1]);
%nys = log10(nys+1);
%figure,imagesc(ts,fs,nys(:,:,2,2)'),axis xy

figure
for i = 1:5
sp(i) = subplot(6,1,i);
imagesc(ts,fs,nys(:,:,i,i)');
%imagesc(ts,fs,log10(10.^nys(:,:,i,i)'./repmat(mean(10.^nys(:,:,i,i),2),1,size(nys,2))'))
caxis([1,3])
%caxis([-3,3])
%imagesc(ts,fs,log10(ys(:,:,i,i))');
%caxis([-5,-2])
axis xy
end
sp(6) = subplot(6,1,6);
imagesc(ts,fs,ys(:,:,1,5)'+ys(:,:,4,5)'+ys(:,:,4,2)'),axis xy
caxis([1.4,3])
Lines(Trial.stc{'w',1}(:),[],'w');
linkaxes(sp,'xy')

nnys = (10.^nys(:,:,i,i)./repmat(mean(10.^nys(:,:,i,i),2),1,size(nys,2)));

figure,
sp(1) = subplot(311);imagesc(ts,fs,log10(nnys(:,:)')),axis xy,caxis([-2,1])
sp(2) = subplot(312);imagesc(ts,fs,log10(ys(:,:,2,2)')),axis xy,caxis([-6-2])
sp(3) = subplot(313);,plot(ts,(dscore))
linkaxes(sp,'x')

[U,S,V] = svd(cov(nys(nniz(nys),:,2,2)));

figure,imagesc(fs,fs,V')

tscore = nys(:,:,2,2)*V;
figure,hist(tscore(:,2),100)
figure,hist(Filter0(gausswin(7)./sum(gausswin(7)),tscore(:,[2])),100)
figure,plot(ts,Filter0(gausswin(5)./sum(gausswin(5)),tscore(:,[2,4])))
figure,hist2(clip(tscore(:,[2,3]),-.01,.01),50,50)


dscore = nanmedian(nnys(:,fs>6&15>fs,2,2),2);%./nanmean(nys(:,fs<4|17<fs,2,2),2);

figure,plot(ts,log10(dscore))
