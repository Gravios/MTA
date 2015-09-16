





% $$$ nfet = MTADxyz('data',dot(repmat(circshift(xyz(:,1,[1,2]),-sh)-circshift(xyz(:,1,[1,2]),sh),[1,5,1]),fxyz(:,[2:5,7],[1,2])-repmat(fxyz(:,1,[1,2]),[1,5,1]),3),'sampleRate',Trial.xyz.sampleRate);
% $$$ nfet.filter('ButFilter',3,2.4,'low');
% $$$ nfet.data = nanmean(nfet.data,2);
% $$$ %nfet.data = log10(abs(nfet.data)+1).*sign(nfet.data);
nfet = MTADxyz('data',dot(circshift(xyz(:,1,[1,2]),-sh)-circshift(xyz(:,1,[1,2]),sh),fxyz(:,3,[1,2])-fxyz(:,1,[1,2]),3),'sampleRate',Trial.xyz.sampleRate);
nfet.filter('ButFilter',3,2.4,'low');
nfet.data = log10(abs(nfet.data)+1).*sign(nfet.data);



sh = 40;
fv = 2.5;
bb =bsxfun(@minus,[circshift(xyz(:,1,:),-sh),circshift(xyz(:,1,:),sh)],xyz(:,1,:));
mfet = Trial.xyz.copy;
mfet.data = sqrt(sum(diff(sq(cross(bb(:,1,:),bb(:,2,:),3))).^2,2));
mfet.filter('ButFilter',3,fv,'low');
mfet.data(mfet.data<0) = 1e-9;

mfet.data = [0;diff(mfet.data)];
mfet.data = sqrt(nansum(mfet.segs(1:mfet.size(1),60,nan).^2,1))';

Lines(Trial.stc{'w'}(:),[],'m');

man = Trial.load('fet','lsppc');
man.filter('ButFilter',3,2,'low');


man.resample(msr);


sh = 20
bb =bsxfun(@minus,[circshift(xyz(:,[1,4,7],:),-sh),circshift(xyz(:,[1,4,7],:),sh)],xyz(:,1,:));
bb = [bb,fxyz(:,4,[1,2]

mfet = Trial.xyz.copy;
mfet.data = sqrt(sum(diff(sq(cross(bb(:,1,:),bb(:,4,:),3))).^2,2));
