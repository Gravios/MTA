
ang = Trial.ang.copy;
xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gausswin(61)./sum(gausswin(61)));
ang.create(Trial,xyz);

ind = xyz(:,1,1)~=0;

ad = circ_dist(ang(ind,1,2,1),ang(ind,1,3,1));
vl = [0;log10(sqrt(sum(diff(xyz(ind,1,:)).^2,3)))];
vld = [0;log10(sqrt(sum(diff(xyz(ind,7,:)-xyz(ind,1,:)).^2,3)))];
%vlm = [0;log10(sqrt(sum(diff(xyz(ind,1:8,:)).^2,3)))];

figure,
hist2([vl,ad],linspace(-2,2,64),linspace(-1.3,1.3,64));

vind = false([ang.size(1),1]);
vind(ind) = vl>0;

figure,hist(ang(vind,1,2,3),1000)
figure,hist(ang(vind,2,4,3),1000)
figure,hist(ang(vind,4,5,3),1000)

figure,hist2([ang(ind,1,2,3),ang(ind,4,5,3)],linspace(50,90,64),linspace(20,90,64));

hist2([ang(vind,1,2,3),ang(vind,2,4,3)],linspace(50,90,64),linspace(20,90,64));


xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gausswin(9)./sum(gausswin(9)));

c1 = cross(xyz(:,1,:)-xyz(:,2,:),xyz(:,3,:)-xyz(:,2,:),3);
c2 = cross(xyz(:,2,:)-xyz(:,3,:),xyz(:,4,:)-xyz(:,3,:),3);
a = acos(dot(c1,c2,3)./(sqrt(sum(c1.^2,3)).*sqrt(sum(c2.^2,3))));
figure,plot(a);Lines(Trial.stc{'w',xyz.sampleRate}.data(:),[],'k');Lines(Trial.stc{'r',xyz.sampleRate}.data(:),[],'r');


ind = a~=0&~isinf(a)&~isnan(a);
nind = sum(ind);
%a = a(ind);

shifts = -600:20:600;
ixy = nan(size(shifts));

for s = 1:numel(shifts)
[out,xb,yb,p]=hist2([a,circshift(a,shifts(s))],linspace(0,2*pi,32),linspace(0,2*pi,32));
pxy = out./nind;
px = histc(a,linspace(0,2*pi,32));
px = px(1:end-1)/nind;
py = histc(circshift(a,shifts(s)),linspace(0,2*pi,32));
py = py(1:end-1)/nind;
pxyn = pxy.*log2(pxy./(px*py'));
ixy(s) = nansum(pxyn(~isinf(pxyn)));
end

[~,arm] = WhitenSignal(a(ind));

[ya,fa,ta] = mtchglong(WhitenSignal(a,[],[],arm),2^8,xyz.sampleRate,2^7,2^7-32,[],[],[],[1,30]);






