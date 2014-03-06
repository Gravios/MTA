MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');

Trial = MTATrial('jg05-20120317')

vel = MTADxyz([],[],Trial.vel,Trial.xyz.sampleRate,[],[],Trial.xyz.model);
vel.filter(gausswin(11)./sum(gausswin(11)));


edges = linspace(-2,2,64);
sbound = -40:40;
state = 'w';

v = vel(Trial.stc{state},1:8,1);

v = vel(:,1:8,1);




ixy = zeros(numel(sbound),size(v,2),size(v,2),1000);
v = vel(x:x+200,1:8,1);
vind = v(:,1)~=0;
nind = numel(vind);
if matlabpool('size')==0,matlabpool open 12,end

xc = 0;
for x = round(linspace(1,10*1000,1000)),
s = 1;
v = vel(x:x+200,1:8,1);
vind = v(:,1)~=0;
nind = numel(vind);
xc = xc+1;

for m = 1:size(v,2)
for o = 1:size(v,2)
for shift = sbound
[out,xb,yb,p]=hist2([log10(v(vind,m)),circshift(log10(v(vind,o)),shift)],edges,edges);
pxy = out./nind;
px = histc(log10(v(vind,m)),xb);
px = px(1:end-1)/nind;
py = histc(circshift(log10(v(vind,o)),shift),yb);
py = py(1:end-1)/nind;
pxyn = pxy.*log2(pxy./(px*py'));
ixy(s,m,o,xc) = nansum(pxyn(~isinf(pxyn)));
s = s+1;
end
s = 1;
end
end

end

[mixy,sixy] = max(ixy);
mixy = sq(mixy);
sixy = (sq(sixy)-ceil(numel(sbound)/2))./Trial.xyz.sampleRate*1000;


%figure
%subplot(1,2,1),imagesc(mixy(:,:))
%subplot(1,2,2),imagesc(sixy(:,:))


figure,plot(var(sq(reshape(sixy(:,:,:),[],1000))',[],2))
figure,plot(diff(Filter0(gausswin(7)./sum(gausswin(7)),reshape(mixy([1:3],[5:7],:),[],1000)')))
Lines(Trial.stc{'w',Trial.xyz.sampleRate/10}.data(1:10,1),[],'b');
Lines(Trial.stc{'w',Trial.xyz.sampleRate/10}.data(1:10,2),[],'b');
Lines(Trial.stc{'r',Trial.xyz.sampleRate/10}.data(1:10,1),[],'r');
Lines(Trial.stc{'r',Trial.xyz.sampleRate/10}.data(1:10,2),[],'r');

% $$$ figure,
% $$$ for f = 1:filtn,
% $$$ subplot2(filtn,2,f,1),imagesc(mixy(:,:,f))
% $$$ subplot2(filtn,2,f,2),imagesc(sixy(:,:,f))
% $$$ end
