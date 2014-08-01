

Trial = MTATrial('jg05-20120310');
Trial = MTATrial('Ed05-20140528');

xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.1,xyz.sampleRate));
vel = xyz.vel([1:3,5,7],[1,2]);


vel = MTADxyz('data',[0;diff(circ_dist(s.ang(:,1,3,1),s.ang(:,5,7,1)))],'sampleRate',xyz.sampleRate);
vel = fet_rhm(Trial,[],'default');
ncp = fet_ncp(Trial,[],'default');

N=10000;
shift=2000;
%X = GetSegs(rhm,1:N,60,0)'/sqrt(N);
X = GetSegs(vel(:,end),(1+shift):(N+shift),60,0)'/sqrt(N);

%C = X'*X;

[U,S,V] = svd(X);

figure,
subplot(121),imagesc(V')
subplot(122),plot(diag(S),'o')
%figure,plot(X*V(:,1))

Xa = zeros(vel.size(1),60);
for i = 1;%:vel.size(2),
Xa(:,:,i) = GetSegs(vel(:,i),1:vel.size(1),60,0)'/sqrt(N)*V(:,:);
end
Xa(nniz(Xa),:) = unity(Xa(nniz(Xa),:));

Na = zeros(vel.size(1),60);
for i = 1:ncp.size(2),
Na(:,:,i) = GetSegs(ncp(:,i),1:ncp.size(1),60,0)'/sqrt(N)*V(:,:);
end
Na(nniz(Na),:) = unity(Na(nniz(Na),:));

figure
sp(1) = subplot(211);imagesc(Xa'),axis xy,caxis([1,2])
sp(2) = subplot(212);imagesc(Na'),axis xy,caxis([1.5,3])
linkaxes(sp,'xy')

Xa = GetSegs(vel(:,i),1:vel.size(1),60,0)'/sqrt(N)*V(:,5:14);
Na = GetSegs(ncp(:,i),1:ncp.size(1),60,0)'/sqrt(N)*V(:,5:14);

%[yv,fv,tv,phiv,fstv] = mtchglong(sum(Xa,2),2^8,vel.sampleRate,2^7,2^7*.875,3,'linear',[],[1,30]); 
% $$$ [yv,fv,tv,phiv,fstv] = mtchglong([unity(Na(:,5)),unity(Xa(:,5))],2^8,vel.sampleRate,2^7,2^7*.875,3,'linear',[],[1,30]); 
% $$$ figure,imagesc(tv,fv,(yv(:,:,1,1)')),
% $$$ imagesc(tv,fv,(yv(:,:,1,1)')),
% $$$ caxis([.5,1])


figure
plot(Xa)
Lines(Trial.stc{'w'}(:,1),[],'k');
Lines(Trial.stc{'w'}(:,2),[],'r');

mXa = abs(median(Xa,2));
vXa = var(Xa,[],2);
vXa(vXa<.0001) = .0001;



wper = Trial.stc{'w'}.cast('TimeSeries');
rper = Trial.stc{'r'}.cast('TimeSeries');

figure,hist(log10(abs(mXa(nniz(mXa)))),1000)

figure,hist2([log10(abs(Xa(nniz(Xa)&~wper.data,1))),log10(abs(Xa(nniz(Xa)&~wper.data,5)))],100,100);

figure,hist2([log10(abs(mXa(nniz(mXa)))),log10(vXa(nniz(mXa)))],100,100);caxis([0,400])
figure,hist2([log10(abs(mXa(nniz(mXa)&~wper.data))),log10(vXa(nniz(mXa)&~wper.data))],100,100);caxis([0,400])
figure,hist2([log10(abs(mXa(nniz(mXa)&wper.data))),log10(vXa(nniz(mXa)&wper.data))],100,100);
figure,hist2([log10(abs(mXa(nniz(mXa)&rper.data))),log10(vXa(nniz(mXa)&rper.data))],100,100);


figure,hist2([log10(abs(vel(nniz(mXa)&nniz(vel)&~wper.data,1))),log10(vXa(nniz(mXa)&nniz(vel)&~wper.data))],100,100);caxis([0,400])

figure,hist2([log10(abs(vel(nniz(mXa)&nniz(vel)&rper.data,1))),log10(vXa(nniz(mXa)&nniz(vel)&rper.data))],100,100);caxis([0,400])
figure,hist2([log10(abs(vel(nniz(mXa)&nniz(vel),1))),log10(vXa(nniz(mXa)&nniz(vel)))],100,100);caxis([0,400])

figure,hist2([log10(abs(Xa(nniz(mXa)&nniz(vel)&~wper.data,1))),log10(vXa(nniz(mXa)&nniz(vel)&~wper.data))],100,100);caxis([0,400])

figure,hist2([log10(abs(Xa(nniz(mXa)&nniz(vel),1))),log10(vXa(nniz(mXa)&nniz(vel)))],100,100);caxis([0,400])

figure,hist2([log10(abs(Xa(nniz(mXa)&nniz(vel),1))),log10(vXa(nniz(mXa)&nniz(vel)))],100,100);caxis([0,200])

figure,hist2([log10(abs(Xa(nniz(Xa)&nniz(vel),1))),log10(Xa(nniz(Xa)&nniz(vel)))],100,100);caxis([0,200])


