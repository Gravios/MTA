

win = 12;
xs  = circshift(xyz.segs(1:size(xyz,1),win,nan),round(win/2),2);
vxs = bsxfun(@minus,xs,xs(round(win/2),:,:));
vxs = sqrt(sum(vxs.^2,3));


figure,plot(median(vxs)./log10(var(vxs)))
Lines(Trial.stc{'n'}(:),[],'g');Lines(Trial.stc{'w'}(:),[],'c');


lvx = xyz.copy;
lvx.data = clip(median(vxs)./log10(var(vxs)),-100,100);
figure,hist(lvx(Trial.stc{'a'}),100);
figure,hist(lvx(Trial.stc{'w'}),100);
figure,hist(lvx(Trial.stc{'p'}),100);

