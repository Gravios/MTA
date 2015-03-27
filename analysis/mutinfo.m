

Trial = MTATrial('jg05-20120317');

marks = {'spine_lower','pelvis_root','spine_middle','spine_upper',...
         'head_back',  'head_left',  'head_front',  'head_right'};

vl = vel(Trial.load('xyz'),marks,[1,2]);
vl.data =  ButFilter(vl.data,3,4/(vl.sampleRate/2),'low');
v = log10(abs(vl(:,1:8)));


edges = linspace(-.5,2,64);
sbound = -30:30;
ixy = zeros([numel(sbound),size(v,2),size(v,2)]);

padding = [0,0];%[-.5,.5];
vind = logical(subsref(cast(resample(Trial.stc{'a'}+padding,xyz),'TimeSeries'),substruct('.',{'data'})));
nind = numel(vind);

s = 1;
for m = 1:size(v,2)
for o = 1:size(v,2)
for shift = sbound
[out,xb,yb,p]=hist2([v(vind,m),circshift(v(vind,o),shift)],edges,edges);
pxy = out./nind;
px = histc(v(vind,m),xb);
px = px(1:end-1)/nind;
py = histc(circshift(v(vind,o),shift),yb);
py = py(1:end-1)/nind;
ixy(s,m,o) = nansum(nansum(pxy.*log2(pxy./(px*py'))));
s = s+1;
end
s = 1;
end
end


[mixy,sixy] = max(ixy);
mixy = sq(mixy;)
sixy = sq(sixy)-ceil(numel(sbound)/2);


figure,
subplot2(1,2,1,1);
imagesc(mixy(:,:));
caxis([0,2]);
colorbar
title('mutual information between marker speeds');
set(gca,'YtickMode','manual');
set(gca,'Ytick',1:8);
set(gca,'YtickLabelMode','manual');
set(gca,'YtickLabel',vl.model.ml('short'));
subplot2(1,2,1,2);
imagesc(sixy(:,:)/vl.sampleRate*1000);
colorbar
title('time lag of maximum mutual information (ms)')
set(gcf,'position',[520, 443, 1046, 289]);


