
Trial = MTATrial('jg05-20120310');
xyz = Trial.load('xyz');
ang = Trial.load('ang');

% $$$ sm = {'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_front'};
% $$$ a = [];
% $$$ for i = 2:5,
% $$$     xm = bsxfun(@minus,xyz(:,{sm{[i-1,i+1]}},:),xyz(:,sm{i},:));
% $$$     a(:,i-1) = abs(acos(dot(xm(:,1,:),xm(:,1,:),3)./prod(sqrt(sum(xm.^2,3)),2)));
% $$$ end
% $$$ 
% $$$ a = MTADxyz('data',a,'sampleRate',xyz.sampleRate);
% $$$ ind = nniz(a);
% $$$ figure,hist2([a(ind,1),a(ind,3)],0:.01:pi/2,0:.01:pi/2);caxis([0,500])
% $$$ ind = Trial.stc{'w'};
% $$$ figure,hist2([a(ind,1),a(ind,3)],0:.01:pi/2,0:.01:pi/2);caxis([0,500])
% $$$ ind = Trial.stc{'r'};
% $$$ figure,hist2([a(ind,1),a(ind,3)],0:.01:pi/2,0:.01:pi/2);caxis([0,500])


Trial = MTATrial('jg05-20120310');

xyz = Trial.load('xyz');
xyz.filter(gtwin(.05,xyz.sampleRate));
fet = xyz.acc([],3);
[ys,fs,ts] = fet_spec(Trial,fet,[],'wcsd','overwrite',true);
figure,imagesc(ts,fs,log10(ys.data(:,:,1,1))'),axis xy
hold on,Lines(Trial.stc{'w',ys.sampleRate}(:),[],'k');


chm = abs(ys(:,:,1,7))./sqrt(ys(:,:,1,1).*ys(:,:,7,7));
chs(:,1) = max(chm(:,fs<14&fs>6),[],2);
chm = abs(ys(:,:,2,4))./sqrt(ys(:,:,2,2).*ys(:,:,4,4));
chs(:,2) = max(chm(:,fs<14&fs>6),[],2);
chm = abs(ys(:,:,4,7))./sqrt(ys(:,:,7,7).*ys(:,:,4,4));
chs(:,3) = max(chm(:,fs<14&fs>6),[],2);

yhs = [];
yhm = log10(ys(:,:,1,1));
yhs(:,1) = median(yhm(:,fs<14&fs>6),2);
yhm = log10(ys(:,:,3,3));
yhs(:,2) = median(yhm(:,fs<14&fs>6),2);
yhm = log10(ys(:,:,5,5));
yhs(:,3) = median(yhm(:,fs<14&fs>6),2);
yhm = log10(ys(:,:,7,7));
yhs(:,4) = median(yhm(:,fs<14&fs>6),2);
yhs = MTADxyz('data',yhs,'sampleRate',ys.sampleRate);



dxy = sqrt(sum(diff(xyz(:,[1,7],[1,2]),1,2).^2,3));
dxy = MTADxyz('data',dxy,'sampleRate',xyz.sampleRate);


spk = Trial.spk.copy;
spk.create(Trial,dxy.sampleRate,'theta',[],'deburst');


units = spk.map(:,1)';
figure,
for u = units,
plot(ang(spk(u),4,5,3),dxy(spk(u)),'.'),xlim([20,100]),ylim([50,250])
waitforbuttonpress
end


[U,S,V] = svd(cov(log10(ys(Trial.stc{'a'},:))));
nfet = zeros([ys.size(1),1]);
nfet(nniz(ys)) = log10(ys(nniz(ys),:))*V(:,1);

figure,plot(nfet)
hold on,Lines(Trial.stc{'w',ys.sampleRate}(:),[],'k');


nfet = MTADxyz('data',nfet,'sampleRate',ys.sampleRate);


ufet = reshape([nfet.data';nanmean([nfet.data,circshift(nfet.data,-1)],2)'],1,[])';
ufet = MTADxyz('data',reshape([ufet';nanmean([ufet,circshift(ufet,-1)],2)'],1,[])',...
               'sampleRate',nfet.sampleRate*4);

dsx = xyz.copy;
dsx.resample(ufet);
dsa = ang.copy;

dsa.resample(ufet);

vel = dsx.vel([],[1,2]);
vel.data(vel.data<0.0001) = 0.0001;
vel.data = log10(vel.data);


wper = Trial.stc{'w'}.cast('TimeSeries');
wper.resample(ufet);
rper = Trial.stc{'r'}.cast('TimeSeries');
rper.resample(ufet);


lrfet = MTADxyz('data',[vel.data(:,'spine_lower'),...
                        vel.data(:,'head_front'),...
                        dsx(:,'spine_lower',3),...
                        dsa(:,'spine_lower','pelvis_root',3),...
                        dsa(:,'spine_lower','pelvis_root',2),...
                        dsa(:,'pelvis_root','spine_middle',4,3),...
                        dsa(:,'spine_lower','head_front',3)],...
                        ufet(:),...
                'sampleRate',ufet.sampleRate);
B = mnrfit(lrfet(Trial.stc{'a'},:),[wper(Trial.stc{'a'})+1+rper(Trial.stc{'a'}).*2],'model','nominal');

y = mnrval(B,lrfet.data);

figure,plot(y)
hold on,Lines(Trial.stc{'w',ufet.sampleRate}(:),[],'k');
hold on,Lines(Trial.stc{'r',ufet.sampleRate}(:),[],'r');

hold on,Lines(Trial.stc{'w'}(:),[],'k');

figure,plot(sum(a.data,2))


