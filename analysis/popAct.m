function popAct(Session)

Session = MTASession('jg05-20120310',{'ufr'});

% $$$ ufr = zeros(size(Session.ufr));
% $$$ 
% $$$ for unit = 1:size(ufr,2)
% $$$     ufr(Session.ufr(:,unit)>0.000003,unit) = unity(log(Session.ufr(Session.ufr(:,unit)>0.000003,unit)));
% $$$     ufr(Session.ufr(:,unit)<=0.000003,unit) = -5.7;
% $$$ end

%imagesc(ufr(1:20000,:)')




ufr = Session.ufr(1:500000,14)+1;

nbin = 1200;
trim = mod(size(ufr,1),nbin);

tpvec = zeros(nbin,round((size(ufr,1)-nbin-1)/nbin-1));
pvec = [];
for i = 0:nbin-1,
tpvec(:,:) = sq(reshape(ufr(1+i:end-trim-nbin+i,:),nbin,[],size(ufr,2)));
abin = tpvec(1:nbin/3,:);
bbin = tpvec(nbin/3+1:2*nbin/3,:);
cbin = tpvec(2*nbin/3+1:nbin,:);
pvec(1+i,:) = (sum(bbin)./sum(abin)+sum(bbin)./sum(cbin))./(sum(abin)+sum(bbin));
end

gvec = reshape(pvec,1,[]);


rind = round((Session.Bhv.getState('rear').state-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
wind = round((Session.Bhv.getState('walk').state-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
rind  = rind(rind(:,1)<500000,:);

clf
plot(gvec)
Lines(rind(:,1),[],'r');
Lines(rind(:,2),[],'g');



ufr = Session.ufr(1:200000,8:28);

nbin = 1000;
trim = mod(size(ufr,1),nbin);
downsample = 4;

round(((size(ufr,1)-nbin-1)/nbin))

tpvec = zeros(nbin/downsample,199,size(ufr,2));
ind = 1;
for i = 0:downsample:nbin-1,
tpvec(ind,:,:) = sq(sum(reshape(ufr(1+i:downsample:end-trim-nbin+i,:),nbin/downsample,[],size(ufr,2)),1));
ind=ind+1;
end

pvec = reshape(tpvec,[],size(ufr,2));



dvec = sqrt(sum((pvec(1:end-1,:)-pvec(2:end,:)).^2,2));



nbin = 125;
trim = mod(size(Session.ufr,1),nbin);


tpvec = zeros(nbin,round((size(Session.ufr,1)-nbin-1)/nbin),size(Session.ufr,2));
for i = 0:nbin-1,
tpvec(i+1,:,:) = sq(sum(reshape(Session.ufr(1+i:end-trim-nbin+i,:),nbin,[],size(Session.ufr,2)),1));
end
pvec = reshape(tpvec,[],size(Session.ufr,2));

dvec = sqrt(sum((pvec(1:end-1,1:100)-pvec(2:end,1:100)).^2,2));


rind = round((Session.Bhv.getState('rear').state-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
wind = round((Session.Bhv.getState('walk').state-1)./Session.xyzSampleRate.*Session.lfpSampleRate);

plot(dvec)

rind  = rind(rind(:,1)<200000,:);

imagesc(sdvec')
Lines(rind(:,1)/4,[],'r');
Lines(rind(:,2)/4,[],'g');


%plot(abs(diff(dvec)))
%Lines(rind(:,1),[],'r');
%Lines(rind(:,2),[],'g');


%ddvec = abs(diff(dvec));

limit = 10^7;

depth = 2000;
sdvec = zeros(size(pvec,1)-depth,depth);
for i = 31:depth,
tvec = sqrt(sum((pvec(1:end-i,:)-pvec(1+i:end,:)).^2,2));
sdvec(:,i) = tvec(1:end-depth+i);
end


rdur = diff(rind,1,2);
rind = rind(rdur>1875,:);
wdur = diff(wind,1,2);
wind = wind(wdur>1875,:);

shift = 0;
tshift = 625;
state = 2;
xlen = 1250;
figure, hold on
for i = 1:6,
subplotfit(i,6),
imagesc(1:xlen,1:500,sdvec(rind(i+shift,state)-tshift:rind(i+shift,state)+xlen-tshift,:)')
%imagesc(1:xlen,1:500,sdvec(wind(i+shift,state)-tshift:wind(i+shift,state)+xlen-tshift,:)')
end


wind = wind(2:end-1,:);



xlen = 2500;

ufet = zeros(size(wind,1),xlen,500);
rfet = zeros(size(rind,1),xlen,500);



state = 2;

for i = 1:size(wind,1),
ufet(i,:,:) = sdvec(wind(i,state)-tshift:wind(i,state)+xlen-tshift-1,:);
end


state = 1;
for i = 1:size(rind,1),
rfet(i,:,:) = sdvec(rind(i,state)-tshift:rind(i,state)+xlen-tshift-1,:);
end


figure,imagesc(sq(mean(ufet))')
figure,imagesc(sq(mean(rfet))')



nfft = 2^11;
win = 2^10;
fs = 1250;

ss = 400000;
[y,f,t] = mtchglong(WhitenSignal(sdvec(1+ss:300000+ss,75)),nfft,fs,win,win*.875,[],[],[],[1 60]);

[yl,fl,tl] = mtchglong(WhitenSignal(Session.lfp(1+ss:300000+ss)),nfft,fs,win,win*.875,[],[],[],[1 60]);


figure,imagesc(t,f,log10(y')),axis xy
figure,imagesc(tl,fl,log10(yl')),axis xy


upt = mean(y(:,f>13&f<24),2);
upl = mean(y(:,f>0&f<5),2);
lpt = mean(yl(:,fl>13&fl<24),2);
lpl = mean(yl(:,fl>0&fl<5),2);

figure
hold on,
plot(tl,lpt./lpl.*10),
plot(t,upt./upl,'r')

ss = 1
figure
hold on,
plot(lpt(ss:end)),
plot(upt.*10^10.7,'r')



msdv = mean(abs(diff(sdvec',1)))';

figure,plot(msdv)
Lines(rind(:,1),[],'r');
Lines(rind(:,2),[],'g');


rvec = pvec(rind(:,1),:);
rdvec = sqrt(sum((rvec(1:end-1,1:100)-rvec(2:end,1:100)).^2,2));

covufr = cov(ufr);

pc_ufr = princomp(covufr);

[max_ufr,ufr_order] = max(pc_ufr,[],2);

[~,uio] = sort(max_ufr);


rind = round((Session.Bhv.getState('rear').state-1)./Session.xyzSampleRate.*Session.lfpSampleRate);


unit = 1;
figure(1242883)
set(gcf,'CurrentCharacter','l');
while 1
clf
ufr_seg_on = sq(GetSegs(ufr(:,ufr_order),rind(unit,1)-3*Session.lfpSampleRate,6*Session.lfpSampleRate,[]));
ufr_seg_off = sq(GetSegs(ufr(:,ufr_order),rind(unit,2)-3*Session.lfpSampleRate,6*Session.lfpSampleRate,[]));
subplot(121);
imagesc(ufr_seg_on');
subplot(122);
imagesc(ufr_seg_off');
    %% Figure controls
    waitforbuttonpress
    whatkey = get(gcf,'CurrentCharacter');
    switch double(whatkey)
      case double('i')
        unit = input('Enter unit #: ');
      case double('n')
        unit = unit+1;
      case double('p')
        unit=unit-1;
      case double('q')
        return
    end
end





usn = GetSegs(ufr,rind(:,1)-0*Session.lfpSampleRate,.5*Session.lfpSampleRate,[]);
[u,s,v] = svd(sq(mean(usn)));

usf = GetSegs(ufr,rind(:,2)-0*Session.lfpSampleRate,.5*Session.lfpSampleRate,[]);



usn_cov = cov(usn);
pcuson = princomp(usn_cov);

