s = MTASession.validate('RS0317-20170517.hcf.all');	
s = MTATrial.validate('jg05-20120317.cof.all');	
xyz = s.load('xyz');
flxyz = xyz.copy;
flxyz.filter('ButFilter',3,.1,'low');
flvxy = flxyz.vel([],[1,2]);
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,3,'low');
fvxy = fxyz.vel([],[1,2]);


stc = s.load('stc','hand_labeled_rev3_jg');
figure,
sp = []
sp(1) = subplot2(4,1,1:3,1);
plot(fvxy(:,5)),hold on,plot(flvxy(:,5))
sp(2) = subplot2(4,1,4,1);
plotSTC(stc,xyz.sampleRate,[],{'walk','rear','turn','pause','groom','sit'},'brgcmy');
linkaxes(sp,'x');

dflxyz = xyz.copy;
dflxyz.data = xyz.data-flxyz.data;
dflxyz.filter('ButFilter',3,8,'low');



dflxyz.resample(20);

sdx = dflxyz.copy;
sdx.data = dflxyz.segs([],100);
sdx.data = reshape(permute(sdx.data,[2,1,3,4]),size(sdx,2),100,18);

U = zeros([size(sdx,1),100,18]);
V = zeros([size(sdx,1),18,18]);
try
for i = 1:size(sdx,1),
    [U(i,:,:),~,V(i,:,:)] = svd(sq(sdx(i,:,:)),0);
end
end

figure,
plot(sq(dflxyz(:,2,:)))
hold on
plot(sq(dflxyz(:,5,:))+.6)
plot(sq(dflxyz(:,1,:))-.6)
plot(sq(dflxyz(:,3,:))-1.2)

figure,imagesc(V(1:i-1,:,1)')
figure,imagesc(sq(V(1:i-1,1,:))')

figure,plot(multiprod(reshape(dflxyz.data,[],18),sq(V(98000,:,1))',[1,2]))
hold on
plot(sq(dflxyz(:,5,:))+.6)

V(:,:,1)))
V(:,:,1)))

pv = V(:,:,1);

for i = 2:size(V,1)
    dpv = sqrt(sum([pv(i,:)-pv(i-1,:);pv(i,:)+pv(i-1,:)].^2,2));
    [~,dind] = min(dpv);
    if dind == 2,
        pv(i,:) = -pv(i,:);
    end
end

pvl = dflxyz.copy('empty');
pvl.data = zeros([1,size(dflxyz,1)]);
for i = 1:size(V,1)
    dfx = dflxyz(i,:,:);
    pvl.data(i) = dfx(:)'*sq(pv(i,:))';
end    

figure,
plot(pvl)


parspec = struct('nFFT',2^7,...
                 'Fs',  pvl.sampleRate,...
                 'WinLength',2^6,...
                 'nOverlap',2^6*.875,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[.1,10]);


[ys,fs,ts] = fet_spec(s,pvl,'mtchglong',false,[],parspec);


sp = [];
figure,
sp(1) = subplot2(4,1,1:3,1);
imagesc(ts,fs,log10(ys.data)');
axis xy;
colormap jet;
caxis([-10,0])
sp(2) = subplot2(4,1,4,1);
plot([1:size(pvl,1)]./pvl.sampleRate,pvl.data)
linkaxes(sp,'x');

figure,plot(multiprod(reshape(dflxyz.data,[],18),permute(pv,[1,3,2]),[2,3]))
hold on

ys.data(ys.data<=0) = eps;

sp = [];
figure,
sp(1) = subplot2(4,1,1:3,1);
imagesc(ts,fs,log10(ys.data)');
axis xy;
colormap jet;
caxis([-10,0])
sp(2) = subplot2(4,1,4,1);
plot([1:size(pvl,1)]./pvl.sampleRate,pvl.data)
linkaxes(sp,'x');


% Test with rat
Trial = MTASession.validate('RS0317-20170517.hcf.all');
Trial = MTATrial.validate('jg05-20120317.cof.all');
Trial = MTATrial.validate('Ed01-20140707.cof.all');
ncp = fet_ncp(Trial);
xyz = Trial.load('xyz');


[ys,fs,ts] = fet_srhmPCA(Trial,10,'mtchglong',false);
ys.data(ys.data<=0) = eps;

srhm = fet_srhmPCA(Trial,10,'mta')

sp = [];
figure,
sp(1) = subplot2(4,1,1:3,1);
imagesc(1:size(ys,1),fs,log10(ys.data)');
axis xy;
colormap jet;
caxis([-10,0])
sp(2) = subplot2(4,1,4,1);
plot([1:size(srhm,1)],srhm.data)
linkaxes(sp,'x');
