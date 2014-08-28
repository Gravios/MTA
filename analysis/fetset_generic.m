function fet = fetset_generic(Trial)

Trial = MTATrial('jg05-20120310');
Trial.load('ang');
Trial.load('xyz');

dwin= gausswin(61)./sum(gausswin(61));
fwin = gausswin(11)./sum(gausswin(11));

Trial.xyz.filter(fwin);

fpv = cat(1,[-2,-2,-2],log10(Filter0(dwin,Trial.vel({'spine_lower','spine_upper','head_front'},[1,2]))));
fpz = cat(1,[-2,-2],log10(Filter0(dwin,Trial.vel({'spine_lower','head_back'},[3]))));

fet =[];
fet = [fet,fpv,fpz];
fet = [fet,Trial.xyz(:,7,3)-Trial.xyz(:,1,3)];
fet = [fet,Trial.ang(:,3,4,2)];
fet = [fet,Trial.ang(:,5,7,2)];
fet = [fet,fet.^3];
fet = MTADxyz([],[],Filter0(fwin,fet),Trial.xyz.sampleRate);
fet.data(fet.data==0) = nan;



Trial = MTATrial('jg05-20120310');
Trial.load('xyz');
Trial.xyz.filter(gtwin(1.25,Trial.xyz.sampleRate));

fet = Trial.xyz.vel([],[1,2]);
fet.data = log10(fet.data);
wfet = zeros(fet.size);
wfet(nniz(fet.data),:) = WhitenSignal(fet(nniz(fet.data),:),[],1);

ys = [];
for i = 1:fet.size(2),
    [ys(:,:,i),fs,ts] = mtchglong(wfet(:,i),2^8,Trial.xyz.sampleRate,2^7,2^6,[],'linear',[],[1,30]);
end

ts = ts + diff(ts(1:2));
nsr = 1/diff(ts(1:2));
ys = cat(1,zeros([nsr*ts(1),size(ys,2),size(ys,3)]),ys,zeros([nsr*ts(1),size(ys,2),size(ys,3)]));

fets = MTADlfp('data',ys,'sampleRate',nsr);


sind = Trial.stc{'w'};
for f = 1:fets.size(2),
for m = 1:fets.size(3),
tfet = fets(nniz(fets(:,f,m)),f,m);
fedgs = prctile(tfet,[2,98]);
NA = histc(tfet,linspace(fedgs(1),fedgs(2),64));
NS = histc(fets(sind,f,m),linspace(fedgs(1),fedgs(2),64));
kld(f,m) = nansum(log((NS/sum(NS(:)))./(NA/sum(NA(:)))).*(NS/sum(NS(:))));
end
end


figure,imagesc(kld')



figure, imagesc(ts,fs,log10(ys')),axis xy
figure
hold on, plot(log10(ys(:,3)),'r')

    rhm = zeros(xyz.size(1),size(ys,2));
    rhm((2^6+1):(size(ys)+2^6),:) = ys;
    rhm = MTADlfp('data',rhm,'sampleRate',xyz.sampleRate);

