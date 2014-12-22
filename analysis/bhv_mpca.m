
Trial = MTATrial('jg05-20120310');
xyz = Trial.load('xyz');

win = 2^6;
hwin = 2^5;


smps = floor(xyz.size(1)/hwin)-1;
inds = 1;
pvecs = zeros([smps,xyz.size([2,3,3])]);
peigv = zeros([smps,xyz.size([2,3,3])]);
for s = 1:smps,
for i = 1:xyz.size(2);
    [~,peigv(s,i,:,:),pvecs(s,i,:,:)] = svd(cov(sq(xyz(inds:inds+win,i,:))));
end
inds = inds+hwin;
end

ssr = xyz.sampleRate/hwin;
pad = round([win/xyz.sampleRate,mod(xyz.size(1)-2^6,2^7)/xyz.sampleRate].*ssr)-[1,0];
szy = size(peigv);
peigv = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),peigv,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
pvecs = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),pvecs,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);

%figure,imagesc(log10(peigv(:,:,1,1))');

gind  = nniz(pvecs(:,1,1,[1,2]));
gind = Trial.stc{'w'};
figure,hist2([sum(pvecs(gind,1,[1,2],1).*pvecs(gind,7,[1,2],1),3),log10(mean(peigv(gind,1:4,1,1),2))],-1:.05:1,-2.8:.1:3.7)

gind  = nniz(pvecs(:,1,1,[1,2]));
gind = Trial.stc{'w'};
jnt = [];
for j = 1:xyz.size(2),
    for k = j:xyz.size(2),
     jnt(:,j,k) = sum(pvecs(gind,j,[1,2],1).*pvecs(gind,k,[1,2],1),3);
    end
end

txyz = multiprod(U,sq(xyz(:,1,:)),[1,2],2);
