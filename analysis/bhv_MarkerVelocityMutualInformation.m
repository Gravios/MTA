MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');

Trial = MTATrial('jg05-20120317')

Trial.load('xyz');
Trial.xyz.filter(gausswin(11)./sum(gausswin(11)));

v = Trial.vel;


edges = linspace(-2,2,16);
filts = 20;
filtn = 5;
sbound = -20:20;
ixy = repmat({zeros(size(v,2),size(v,2))},[numel(sbound)],1);



v = Trial.vel;
v = Filter0(gausswin(11)./sum(gausswin(11)),v);

vind = v(:,1)~=0;
nind = numel(vind);


if matlabpool('size')==0,matlabpool open 12,end

for m = 1:size(v,2)
for o = 1:size(v,2)
parfor s = 1:numel(sbound)
[out,xb,yb,p]=hist2([log10(v(vind,m)),circshift(log10(v(vind,o)),sbound(s))],edges,edges);
pxy = out./nind;
px = histc(log10(v(vind,m)),xb);
px = px(1:end-1)/nind;
py = histc(circshift(log10(v(vind,o)),sbound(s)),yb);
py = py(1:end-1)/nind;
ixy{s}(m,o) = nansum(nansum(pxy.*log2(pxy./(px*py'))));
end
end
end


ixy = permute(reshape(cell2mat(ixy),size(v,2),[],size(v,2)),[2,1,3]);

[mixy,sixy] = max(ixy);
mixy = sq(mixy);
sixy = sq(sixy)-ceil(numel(sbound)/2);


subplot2(1,2,1,1),imagesc(mixy)
subplot2(1,2,1,2),imagesc(sixy./Trial.xyz.sampleRate*1000)





%% Vel Seg

MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');

Trial = MTATrial('jg05-20120317')

figure,
hold on,

v = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate);
v.filter(gausswin(121)./sum(gausswin(121)));
plot(log10(v(:,1)))

v = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate);
v.filter(gausswin(241)./sum(gausswin(241)));
plot(log10(v(:,1)),'r')

v = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate);
v.filter(gausswin(481)./sum(gausswin(481)));
plot(log10(v(:,1)),'g')


figure % Spine Lower
v = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate);
v.filter(gausswin(121)./sum(gausswin(121)));
subplot(4,1,1),hist(log10(v(v(:,1)~=0,1)),1000)

v = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate);
v.filter(gausswin(241)./sum(gausswin(241)));
subplot(4,1,2),hist(log10(v(v(:,1)~=0,1)),1000)

v = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate);
v.filter(gausswin(481)./sum(gausswin(481)));
subplot(4,1,3),hist(log10(v(v(:,1)~=0,1)),1000)

v = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate);
v.filter(gausswin(961)./sum(gausswin(961)));
subplot(4,1,4),hist(log10(v(v(:,1)~=0,1)),1000)


figure % head front


for m = 1:8,
v = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate);
v.filter(gausswin(11)./sum(gausswin(11)));
subplot2(5,8,1,m),hist(log10(v(v(:,m)~=0,m)),128);

v = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate);
v.filter(gausswin(121)./sum(gausswin(121)));
subplot2(5,8,2,m),hist(log10(v(v(:,m)~=0,m)),128);

v = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate);
v.filter(gausswin(241)./sum(gausswin(241)));
subplot2(5,8,3,m),hist(log10(v(v(:,m)~=0,m)),128);

v = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate);
v.filter(gausswin(481)./sum(gausswin(481)));
subplot2(5,8,4,m),hist(log10(v(v(:,m)~=0,m)),128);

v = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate);
v.filter(gausswin(961)./sum(gausswin(961)));
subplot2(5,8,5,m),hist(log10(v(v(:,m)~=0,m)),128);
end