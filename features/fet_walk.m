Trial = MTATrial('jg05-20120310');
Trial = MTATrial('Ed05-20140528');
Trial = MTATrial('jg05-20120311');

xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.1,xyz.sampleRate));
vel = xyz.vel([1:3,5,7],[1,2]);

ang = Trial.ang.copy;
ang.create(Trial,xyz);
vel = [circ_dist(ang(:,2,3,1),ang(:,1,2,1)),...
        circ_dist(ang(:,3,4,1),ang(:,2,3,1)),...
        circ_dist(ang(:,4,5,1),ang(:,3,4,1)),...
        circ_dist(ang(:,5,7,1),ang(:,4,5,1))];
vel = MTADxyz('data',vel,'sampleRate',ang.sampleRate);



N=10000;
shift=2000;
minD = 64;
X = GetSegs(vel(:,end),(1+shift):(N+shift),minD,0)'/sqrt(N);
[U,S,V] = svd(X);

figure,
subplot(121),imagesc(V')
subplot(122),plot(diag(S),'o')


Xa=[];
for i = 1:vel.size(2),
Xa(:,i) = GetSegs(vel(:,i),1:vel.size(1),minD,0)'/sqrt(N)*V(:,1);
end

mXa = abs(median(Xa,2));
vXa = var(Xa,[],2);
vXa(vXa<.000001) = .000001;


wper = Trial.stc{'w'}.cast('TimeSeries');
rper = Trial.stc{'r'}.cast('TimeSeries');

mXa = mXa(1:end-1);
vXa = vXa(1:end-1);

edm = linspace(-3,.5,100);
edv =linspace(-5,0,100);

% $$$ figure,hist2([log10(abs(Xa(nniz(Xa),1))),log10(abs(Xa(nniz(Xa),5)))],edm,edv);caxis([0,500])
% $$$ figure,hist2([log10(abs(Xa(nniz(Xa)&wper.data,1))),log10(abs(Xa(nniz(Xa)&wper.data,5)))],edm,edv);caxis([0,1000])
% $$$ figure,hist2([log10(abs(Xa(nniz(Xa)&~wper.data,1))),log10(abs(Xa(nniz(Xa)&~wper.data,5)))],edm,edv);caxis([0,1000])
% $$$ figure,hist2([log10(abs(Xa(nniz(Xa)&rper.data,1))),log10(abs(Xa(nniz(Xa)&rper.data,5)))],edm,edv);caxis([0,1000])


figure,hist2([log10(abs(mXa(nniz(mXa)))),log10(abs(vXa(nniz(mXa))))],edm,edv);caxis([0,500])
figure,hist2([log10(abs(mXa(nniz(mXa)&~wper.data))),log10(vXa(nniz(mXa)&~wper.data))],edm,edv);caxis([0,500])
figure,hist2([log10(abs(mXa(nniz(mXa)&wper.data))),log10(vXa(nniz(mXa)&wper.data))],edm,edv);
figure,hist2([log10(abs(mXa(nniz(mXa)&rper.data))),log10(vXa(nniz(mXa)&rper.data))],edm,edv);



%% ANG VEL
Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz').filter(.25,Trial.xyz.sampleRate);
vfet = xyz.vel([1,7],[1,2]);
vfet.data = diff(vfet.data,1,2);
ang = create(Trial.ang.copy,Trial,xyz);

afet = xyz.copy;
afet.data = circ_dist(circshift(ang(:,1,7,1),-120),ang(:,1,7,1));
bfet = xyz.copy;
bfet.data = circ_dist(circshift(ang(:,1,7,1),-60),ang(:,1,7,1));



[bas,fs,ts] = fet_spec(Trial,afet,'mtchglong',false);
[vas,fs,ts] = fet_spec(Trial,vfet,'mtchglong',false);

figure,imagesc(ts,fs,log10(circshift(bas.data,0)')),
axis xy,
colormap jet,
caxis([-9,-4])
Lines(Trial.stc{'w',bas.sampleRate}(:),[],'k');
Lines(Trial.stc{'n',bas.sampleRate}(:),[],'r');


ind = nniz(vas)&nniz(bas);
figure,hist2(log10([bas(ind,1),vas(ind,1)]),100,100)

ind = resample(Trial.stc{'w'}.cast('TimeSeries'),bas);
ind = resample(Trial.stc{'n'}.cast('TimeSeries'),bas);
ind = logical(ind.data)&nniz(vas)&nniz(bas);
ind = nniz(vas)&nniz(bas);
figure,hist2(log10([bas(ind,1),vas(ind,1)]),linspace(-10,-1,100),linspace(-7,-.5,100))


ind = Trial.stc{'a'};
ind = Trial.stc{'n'};
ind = Trial.stc{'w'};
figure,bar(aedgs,histc(log10(afet(ind)),aedgs),'histc')

aedgs = linspace(-4,1.7,100);
vedgs = linspace(-2,2,100);
figure,bar(aedgs,histc(log10(afet(ind)),aedgs),'histc')

edgs = {aedgs,vedgs};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
ns = 6;

sts = 'awnr';
stc = 'kbgr';
ind = Trial.stc{'a'};
hist2(log10([afet(ind),vfet(ind)]),aedgs,vedgs);
hold on,
for s = sts,
    i = Trial.stc{s};
    b = log10([afet(i),vfet(i)]);
    o = hist2(b,aedgs,vedgs);
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[20,20],'linewidth',1.5,'Color',stc(s==sts))
end


figure
for s = sts,
    subplot(1,numel(sts),find(s==sts));
    ind = Trial.stc{s};
    hist2(log10([afet(ind),vfet(ind)]),aedgs,vedgs);
end



figure
ind = Trial.stc{'w'};
hist2(log10([afet(ind),vfet(ind)]),aedgs,vedgs);colormap jet,caxis([0,100])
figure
ind = Trial.stc{'n'};
hist2(log10([afet(ind),vfet(ind)]),aedgs,vedgs);colormap jet,caxis([0,100])


vedgs = linspace(1.2,1.7,100);
figure
ind = Trial.stc{'w'};
hist2(log10([afet(ind),ang(ind,1,4,3)]),aedgs,vedgs);colormap jet,caxis([0,100])
figure
ind = Trial.stc{'n'};
hist2(log10([afet(ind),ang(ind,1,4,3)]),aedgs,vedgs);colormap jet,caxis([0,100])

vedgs = linspace(-1,0.5,100);
dat = MTADxyz('data',[afet(:),ang(:,1,2,2)+ang(:,2,3,2)],'sampleRate',xyz.sampleRate);
figure
ind = Trial.stc{'w'};
hist2(log10(dat(ind,:)),aedgs,vedgs);colormap jet,caxis([0,100])
figure
ind = Trial.stc{'n'};
hist2(log10(dat(ind,:)),aedgs,vedgs);colormap jet,caxis([0,100])


vedgs = linspace(.5,1.5,100);
figure
ind = Trial.stc{'w'};
hist2(log10([afet(ind),ang(ind,1,3,3)]),aedgs,vedgs);colormap jet,caxis([0,100])
figure
ind = Trial.stc{'n'};
hist2(log10([afet(ind),ang(ind,1,3,3)]),aedgs,vedgs);colormap jet,caxis([0,100])

d = circshift(circ_dist(circshift(ang(:,2,7,1),-30),ang(:,2,7,1)),15);
d = circshift(sum(GetSegs(d,1:size(d,1),30,nan))',15);

ds = circshift(circ_dist(circshift(ang(:,1,4,1),-30),ang(:,1,4,1)),15);
ds = circshift(sum(GetSegs(d,1:size(d,1),30,nan))',15);
d = abs(ds)+abs(d);

afet.data = abs(d);

figure,plot(abs(d)),Lines(Trial.stc{'n'}(:),[],'r');



ind = Trial.stc{'a'};
ind = Trial.stc{'n'};
ind = Trial.stc{'w'};
figure,bar(aedgs,histc(log10(afet(ind)),aedgs),'histc')
