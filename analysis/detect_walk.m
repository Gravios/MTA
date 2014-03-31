function fet = detect_walk(Trial)


%Trial = MTATrial('jg05-20120317');
xyz = Trial.xyz.copy; 
xyz.load(Trial);
xyz.filter(gausswin(91)./sum(gausswin(91)));

% $$$ c='gbmr';
fet = [];
w = [60,120,180,240];
%figure
for m = 4:8,
for k = 1:4,
wins = round(linspace(10,w(k),30));
nw = numel(wins);
sxy = zeros([xyz.size([1,2]),2,nw]);

for i = 1:9,
for j = 1:nw
sxy(:,i,:,j) = circshift(xyz(:,i,[1,2])-circshift(xyz(:,i,[1,2]),-wins(j)),round(wins(j)/2));
end
end

ims = circshift(sq(sum(sxy(:,3,:,1:30).*sxy(:,m,:,1:30),3)),round(-.5*wins(end)))';
%figure,
%plot(mean(diff(diff(ims')'))',c(m-1)),
fet(:,k,m-3) = mean(diff(diff(ims')'))';
end
%hold on,plot(mean(fet,2),c(m-1));
end

fet = sq(mean(fet,2));
fet = sq(mean(fet,2));
fet = [0;fet];

% $$$ 
% $$$ Lines(Trial.stc{'w'}(:),[],'k');
% $$$ Lines([],0,'k');
% $$$ 
figure,plot(fetb,'b');
hold on,plot(feth,'r');
Lines(Trial.stc{'w'}(:),[],'k');
Lines([],0,'k');
% $$$ 
% $$$ ims(ims<1) = 1;
% $$$ %figure,imagesc(log10(ims)./repmat(max(log10(ims),[],2),[1,size(ims,2)])),axis xy;
% $$$ %figure,imagesc(diff(log10(ims)./repmat(max(log10(ims),[],2),[1,size(ims,2)]))),axis xy;
% $$$ 
% $$$ figure,plot(var(diff(log10(ims)./repmat(max(log10(ims),[],2),[1,size(ims,2)]))',[],2))
% $$$ hold on,plot(mean(diff(log10(ims)./repmat(max(log10(ims),[],2),[1,size(ims,2)]))',2),'r')
% $$$ Lines(Trial.stc{'w'}(:),[],'k');
% $$$ Lines(Trial.stc{'r'}(:),[],'r');
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ index = 26400:26700;
% $$$ pview = zeros([xyz.size(1),1]);
% $$$ pview(index) = 1;
% $$$ 
% $$$ 
% $$$ figure,plot3(xyz(index,1,1),xyz(index,1,2),xyz(index,1,3))
% $$$ hold on,plot3(xyz(index,7,1),xyz(index,7,2),xyz(index,7,3),'r')
% $$$ hold on,plot3(xyz.data(index(1),7,1),xyz.data(index(1),7,2),xyz.data(index(1),7,3),'.g')
% $$$ 
% $$$ 
% $$$ wper = Trial.stc{'w'}.copy;wper.cast('TimeSeries');
% $$$ wper.data = wper.data(1:end-1);
% $$$ windex = wper(:)&pview;
% $$$ 
% $$$ rper = Trial.stc{'r'}.copy;rper.cast('TimeSeries');
% $$$ rper.data = rper.data(1:end-1);
% $$$ rindex = rper(:)&pview;
% $$$ 
% $$$ hold on,plot3(xyz(windex,1,1),xyz(windex,1,2),xyz(windex,1,3),'.g')
% $$$ hold on,plot3(xyz(windex,7,1),xyz(windex,7,2),xyz(windex,7,3),'.g')
% $$$ 
% $$$ hold on,plot3(xyz(rindex,7,1),xyz(rindex,7,2),xyz(rindex,7,3),'.m')
% $$$ hold on,plot3(xyz(rindex,1,1),xyz(rindex,1,2),xyz(rindex,1,3),'.m')
% $$$ 
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ Trial = MTATrial('jg05-20120317');
% $$$ 
% $$$ xyzs = xyz.copy;
% $$$ xyzs.clear;
% $$$ xyzs.load(Trial);
% $$$ xyzs.filter(gausswin(7)./sum(gausswin(7)));
% $$$ 
% $$$ 
% $$$ xyzl = xyz.copy;
% $$$ xyzl.clear;
% $$$ xyzl.load(Trial);
% $$$ xyzl.filter(gausswin(91)./sum(gausswin(91)));
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ % $$$ figure,
% $$$ % $$$ hold on,
% $$$ % $$$ plot3(xyzl(index,1,1),xyzl(index,1,2),xyzl(index,1,3))
% $$$ % $$$ plot3(xyzl(index,7,1),xyzl(index,7,2),xyzl(index,7,3),'r')
% $$$ % $$$ plot3(xyzl.data(index(1),7,1),xyzl.data(index(1),7,2),xyzl.data(index(1),7,3),'.g')
% $$$ % $$$ 
% $$$ % $$$ plot3(xyzs(index,1,1),xyzs(index,1,2),xyzs(index,1,3))
% $$$ % $$$ plot3(xyzs(index,7,1),xyzs(index,7,2),xyzs(index,7,3),'r')
% $$$ % $$$ plot3(xyzs.data(index(1),7,1),xyzs.data(index(1),7,2),xyzs.data(index(1),7,3),'.g')
% $$$ 
% $$$ xyzc = xyzs.copy;
% $$$ xyzc.data = cat(2,xyzs(:,[1,3,4,7],:),xyzl(:,[1,3,4,7],:));
% $$$ 
% $$$ nang = Trial.ang.copy;
% $$$ nang.create(Trial,xyzc);
% $$$ 
% $$$ % $$$ figure,plot(diff(Filter0(gausswin(5)./sum(gausswin(5)),nang(:,1,3,3))))
% $$$ % $$$ %hold on,plot(diff(Filter0(gausswin(5)./sum(gausswin(5)),nang(:,2,4,3))),'r'),
% $$$ % $$$ Lines(Trial.stc{'w'}(:),[],'k');
% $$$ % $$$ Lines(Trial.stc{'r'}(:),[],'r');
% $$$ 
% $$$ 
% $$$ [ya,fa,ta] = mtchglong(WhitenSignal(diff(Filter0(gausswin(3)./sum(gausswin(3)),nang(:,1,5,3)))),2^8,nang.sampleRate,2^7,2^7-1,[],[],[],[.01,30]);
% $$$ 
% $$$ 
% $$$ vellf = sqrt(sum(diff(xyzl(:,:,:)).^2,3));
% $$$ velsf = sqrt(sum(diff(xyzs(:,:,:)).^2,3));
% $$$ 
% $$$ % $$$ sp = [];
% $$$ % $$$ figure,
% $$$ % $$$ sp(1) = subplot(211);imagesc(ta+1/(2^6),fa,log10(ya')),axis xy,
% $$$ % $$$ Lines(Trial.stc{'w',1}(:),[],'k');
% $$$ % $$$ sp(2) = subplot(212);plot(ta,vellf(round(ta.*nang.sampleRate)+1));
% $$$ % $$$ %plot((ta+1/(2^8),fa,mean(log10(ya(:,fa>')),axis xy,
% $$$ % $$$ Lines(Trial.stc{'w',1}(:),[],'k');
% $$$ % $$$ linkaxes(sp,'x') 
% $$$ 
% $$$ yap = cat(1,zeros(64,1),mean(ya(:,fa>3&fa<7),2),zeros(64,1));
% $$$ 
% $$$ % $$$ figure,plot(yap)
% $$$ % $$$ Lines(Trial.stc{'w'}(:)-64,[],'k');
% $$$ 
% $$$ 
% $$$ figure,
% $$$ vel = vellf(:,1);
% $$$ ves = velsf(:,1);
% $$$ ind = ~isnan(vel)&~isinf(vel)&~isnan(yap)&~isinf(yap);
% $$$ hist2([clip(log10(vel(ind)),-3,2),log10(yap(ind))],linspace(-2,1,60),linspace(-3,-.4,60))
