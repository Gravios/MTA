
Trial = MTATrial.validate('jg05-20120317.cof.all');

xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,4,'low');
fang = create(MTADang,Trial,fxyz);

chans = [65:1:96];
marker = 'spine_lower';

Trial = MTATrial.validate(sname);

xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);


Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',chans);

stc = Trial.load('stc','NN0317R');
states = {'loc','rear','pause','groom','sit'};

%wlfp = WhitenSignal(lfp.data,[],1);    
wlfp = WhitenSignal(lfp.data,round(180*lfp.sampleRate),true);    


%% compute low frequency stuff
try,delete(gcp('nocreate')),end
parp = parpool(10);

tl=[];    fl=[];    yl=[];
spectral.nfft = 2^11;
spectral.window = 2^10;
parfor i = 1:lfp.size(2),
    [yl(:,:,i),fl(:,i),tl(:,i)] = mtchglong(wlfp(:,i),...
                                            spectral.nfft,...
                                            lfp.sampleRate,...
                                            spectral.window,...
                                            spectral.window*0.875,...
                                            [],[],[],[1,40]);
end
fl = fl(:,1);
tl = tl(:,1);
tindShift = round(spectral.window/2/lfp.sampleRate*1/diff(tl(1:2,1))-1);
yld = MTADlfp('data',cat(1,zeros([tindShift,size(yl,2),size(yl,3)]),yl),...
              'sampleRate',1/diff(tl(1:2,1)));
yld.data(yld.data==0)=nan;
yld.data(yld.data==0)=1;
yld.data(isnan(yld.data))=1;
tpow = yld.copy;
tpow.clear;    
tpow.data = [nanmean(nanmean(yld(:,6<fl&fl<12,[6,7]),2)./nanmean(yld(:,fl<5,[6,7]),2),3),...
             nanmean(nanmean(yld(:,6<fl&fl<12,[16:20]),2)./nanmean(yld(:,fl<5,[16:20]),2),3)];    
yld.data = log10(yld.data);

nyld = yld.copy;
nyld.data = (clip(yld.data,0.01,100)-repmat(nanmedian(clip(yld.data,0.01,100)),[yld.size(1),1,1]))./...
    repmat(nanstd(clip(yld.data,0.01,100)),[yld.size(1),1,1]);    
tl = [[0:tindShift-1]'./yld.sampleRate;tl+spectral.window/2/lfp.sampleRate];    



figure,
imagesc(tl,fl,log10(yld(:,:,7)')),
axis xy
caxis(



figure,plot(unity(diff(circ_dist(fang(:,1,3,1),fang(:,1,4,1)))))
hold on,plot(unity(diff(fxyz(:,1,3))))
Lines(Trial.stc{'x'}(:),[],'r');


figure,
ind = Trial.stc{'a-x'};
plot(circ_dist(ang(ind,2,4,1),ang(ind,1,3,1)),xyz(ind,1,3),'.')
hold on
ind = Trial.stc{'x'};
dind = diff(ind.data,1,2);
ind.data = ind(dind>120,:);
plot(circ_dist(ang(ind,2,4,1),ang(ind,1,3,1)),xyz(ind,1,3),'.r')



Trial = MTATrial.validate('Ed05-20140529.ont.all');
ncp = fet_ncp(Trial);

xyz = Trial.load('xyz');
hxyz = xyz.copy;
hxyz.filter('ButFilter',3,[0.1,15],'bandpass');
ang = create(MTADang,Trial,xyz);


figure,plot((1:ncp.size(1))/ncp.sampleRate,nunity(ncp.data))
hold on, 
plot((1:ncp.size(1))/ncp.sampleRate,nunity(ang(:,4,5,3)).*50)
%plot((1:ncp.size(1))/ncp.sampleRate,[0;nunity(diff(hxyz(:,3,3)))].*25)
plot((1:ncp.size(1))/ncp.sampleRate,nunity(hxyz(:,5,3)).*100)
Lines(Trial.stc{'s',1}(:),[],'r');