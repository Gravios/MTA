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
