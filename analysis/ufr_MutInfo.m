MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
Trial = MTATrial('jg05-20120317','all');
Trial.xyz.load(Trial);
Trial.ufr.create(Trial,Trial.xyz,'rear&walk');

edges = linspace(-2,2,64);
sbound = -10:10;

v = Trial.ufr.data;
v = log10(v);
v(isnan(v)) = -2;

%ixy = zeros(numel(sbound),size(v,2),size(v,2));
ixy = repmat({zeros(numel(sbound),size(v,2))},1,size(v,2));

nind = size(v,1);
s = 1;

if matlabpool('size')==0,matlabpool open 12,end
parfor m = 1:size(v,2)
s = 1
for o = 1:size(v,2)
for shift = sbound

[out,xb,yb,p]=hist2([v(:,m),circshift(v(:,o),shift)],edges,edges);
pxy = out./nind;
px = histc(v(:,m),xb);
px = px(1:end-1)/nind;
py = histc(circshift(v(:,o),shift),yb);
py = py(1:end-1)/nind;
pxyn = pxy.*log2(pxy./(px*py'));
ixy{m}(s,o) = nansum(pxyn(~isinf(pxyn)));
s = s+1;
end
s = 1;
end
end

sxy = reshape(cell2mat(ixy),[],size(v,2),size(v,2));
[mixy,sixy] = max(sxy);
%[mixy,sixy] = max(mint);
mixy = sq(mixy);
sixy = (sq(sixy)-ceil(numel(sbound)/2))./Trial.xyz.sampleRate*1000;


figure
subplot(1,2,1),imagesc(mixy(:,:))
subplot(1,2,2),imagesc(sixy(:,:))


pfr = MTAAknnpfs(Trial,[],'rear');
pfw = MTAAknnpfs(Trial,[],'walk');
pfh = MTAAknnpfs(Trial,[],'hwalk');
pfl = MTAAknnpfs(Trial,[],'lwalk');

us = [29,48];
us = [55,79];
us = [55,95];
us = [70,88];
us = [70,47];
us = [21,44];
us = [15,17];
us = [39,66];
us = [38,29];
us = [38,22];
us = [38,27];
us = [21,29];
us = [41,27];
us = [52,18];
us = [20,62];

us = [9,65];
us = [35,47];
us = [70,35];


figure,
subplot(221),pfw.plot(us(1));
subplot(222),pfw.plot(us(2));
subplot(223),pfr.plot(us(1));
subplot(224),pfr.plot(us(2));


figure,
subplot(221),pfh.plot(us(1));
subplot(222),pfh.plot(us(2));
subplot(223),pfl.plot(us(1));
subplot(224),pfl.plot(us(2));






