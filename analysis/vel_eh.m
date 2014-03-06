MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
Trial = MTATrial('jg05-20120317','all');


edges = linspace(-2,2,32);
sbound = 60;
niter = 20000;

cv = clip(vel(:,7),0,100);
dcv = clip(diff(cv),-10,10);
ddcv = clip(diff(dcv),-5,5);
addcv = abs(ddcv(ddcv~=0));

vel = MTADxyz([],[],[Trial.vel([],[1,2]),Trial.vel([],[3])],Trial.xyz.sampleRate,[],[],Trial.xyz.model);
vel.filter(gausswin(61)./sum(gausswin(61)));


v = vel(:,:);
%vind = v(:,1)~=0;
v = log10(v(:,:));
nind = sbound;


%ixy = repmat({zeros(size(v,2),size(v,2))},1,niter);
ixy = repmat({zeros(size(v,2),size(v,2),sbound)},1,niter);

if matlabpool('size')==0,matlabpool open 12,end

parfor i  = 1:niter
for m = 1:size(v,2)
for o = 1:size(v,2)
% $$$ vseg = [v(i:(i+sbound),m),v((i:(i+sbound))+sbound,o)];
% $$$ if sum(isnan(vseg(:)))|sum(isinf(vseg(:))),continue,end
% $$$ [out,xb,yb,p]=hist2(vseg,edges,edges);
% $$$ pxy = out./nind;
% $$$ px = histc(vseg(:,1),xb);
% $$$ px = px(1:end-1)/nind;
% $$$ py = histc(vseg(:,2),yb);
% $$$ py = py(1:end-1)/nind;
% $$$ pxyn = pxy.*log2(pxy./(px*py'));
% $$$ ixy{i}(m,o) = nansum(pxyn(~isinf(pxyn)));

for s = 1:sbound
vseg = [v(i:(i+sbound),m),v((i:(i+sbound))+s,o)];
if sum(isnan(vseg(:)))|sum(isinf(vseg(:))),continue,end
[out,xb,yb,p]=hist2(vseg,edges,edges);
pxy = out./nind;
px = histc(vseg(:,1),xb);
px = px(1:end-1)/nind;
py = histc(vseg(:,2),yb);
py = py(1:end-1)/nind;
pxyn = pxy.*log2(pxy./(px*py'));
ixy{i}(m,o,s) = nansum(pxyn(~isinf(pxyn)));
end

end
end
end


mxy = permute(reshape(cell2mat(ixy),size(v,2),size(v,2),niter),[3,1,2]);
fxy = reshape(Filter0(gausswin(31)./sum(gausswin(31)),mxy),[],27,27);

[mixy,sixy] = max(ixy);

[mixy,sixy] = max(mint);
mixy = sq(mixy);
sixy = (sq(sixy)-ceil(numel(sbound)/2))./Trial.xyz.sampleRate*1000;


figure
subplot(1,2,1),imagesc(mixy(:,:))
subplot(1,2,2),imagesc(sixy(:,:))

figure,col = jet(27);hold on
for i = 18:27 ,plot(mxy(:,i,i),'Color',col(i,:)),end
Lines(Trial.stc{'r'}(1:10,1),[],'r');
Lines(Trial.stc{'r'}(1:10,2),[],'r');
Lines(Trial.stc{'w'}(1:20,1),[],'b');
Lines(Trial.stc{'w'}(1:20,2),[],'b');

figure,col = jet(27);hold on
s = 1;for i = 10:2:27 ,plot(mxy(:,i,i),'Color',col(s,:)),s=s+1;,end
Lines(Trial.stc{'r'}(1:15,1)+sbound,[],'r');
Lines(Trial.stc{'r'}(1:15,2)+sbound,[],'r');
Lines(Trial.stc{'w'}(1:30,1)+sbound,[],'b');
Lines(Trial.stc{'w'}(1:30,2)+sbound,[],'b');
hold on,plot((sbound+1):(niter+sbound),v(1:niter,:)-2)

figure,plot(mxy(1:20:end,1,1),mxy(1:20:end,26,21),'.')


figure,hist2([mxy(:,3,3)./mxy(:,6,6),mxy(:,3,3)./mxy(:,8,8)],200,200)

figure,hist2([mxy(1:20:end,1,1),mxy(1:20:end,7,7)],30,30)

figure,plot(mxy(1:20:end,10,10),mxy(1:20:end,25,25),'.');

