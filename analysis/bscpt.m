s = MTASession('jg05-20120315');
b1t = s.ang(:,4,5,3);
b2t = s.ang(:,5,7,2);
b1t(isnan(b1t))=1;
b2t(isnan(b2t))=1;
b1m = b1t/mean(b1t(5000:10000));
b1 = b1m-mean(b1m(b1m~=1));
b2 = b2t/mean(b2t(5000:10000));
figure,hold on,plot(b1*10),plot(b2,'r')
b = b1-b2;
wb = WhitenSignal(b);
nfft = 2^8;
window = 2^7;
[y,f,t] = mtchglong(wb,nfft,s.xyzSampleRate,window,(1-1/2^2)*window,[],[],[],[0.1 20]);
dt = 1/diff(t(1:2));
bpow = median(y(:,f<13&f>7),2);
bpow(bpow==0) = 1;
bp = zeros(size(bpow));
bp(log10(bpow)>-4.7) = 1;
bper = round(ThreshCross(bp,0.5,1)+dt/2);
bint = zeros(size(bper,1),1);
for i = 1:size(bper,1),
    bint(i) = sum(bpow(bper(i,1):bper(i,2)));
end
cbper = bper(bint<200,:);
cbint = bint(bint<200);
sniffingPeriods = cbper(log10(cbint)>-4,:);
sniffing = zeros(size(s.xyz,1),1);
maxSniff = LocalMinima(-bpow,5,-2e-5);
for i = 1:size(sniffingPeriods,1),
    sniffing(round(sniffingPeriods(i,1)/dt*s.xyzSampleRate):round(sniffingPeriods(i,2)/dt*s.xyzSampleRate))=1;
end
snPowPeaks = round((maxSniff+dt/2)/dt*s.xyzSampleRate);
bsel = b.*sniffing;
bPeaks = LocalMinima(bsel,6);

figure
Lines(st(bPeaks),[-10,10],'k');
hold on
plot(st,b)


figure
Lines(bPeaks,[-10,10],'k');
hold on
plot(b)
plot(vb,'g')

lfp = LoadBinary([s.path.nlx s.name '/' s.name '.lfp'],[65:96],96,[],[],[],[s.syncPeriods(1,1), s.syncPeriods(end,2)])';

wind = 1250
blfp=GetSegs(lfp,round(bPeaks/s.xyzSampleRate*1250-wind/2),wind,[]);
mbl =sq(max(blfp,[],1));
nbl =sq(min(blfp,[],1));

start=140600;
stop =141300;
bpind =find(bPeaks<stop&bPeaks>start);

chan =13;
figure,hist(mbl(:,chan),400)
hold on,Lines([1:5].*std(mbl(:,chan))+mean(mbl(:,chan)),[0,100],'r')

upthresh = 4*std(mbl(:,chan))+mean(mbl(:,chan));
downthresh = -4*std(nbl(:,chan))+mean(nbl(:,chan));
cblfp = blfp(:,mbl(:,chan)<upthresh&nbl(:,chan)>downthresh,:);
figure,imagesc(1:wind,1:length(bpind),blfp(:,bpind,7)');
colorbar
caxis([-2e4 2e4])

figure,plot(mean(cblfp(:,:,chan),2))


sl
btemp = s.ang(:,s.Model.gmi('spine_upper'),s.Model.gmi('head_back'),3);
btemp(isnan(btemp)) =mean(btemp(~isnan(btemp)));

kern = gausswin(11,0.1);
kern =kern/sum(kern);
bf = Filter0(kern,btemp);
vb = [0;diff(bf)];
vbf = Filter0(kern,vb);
ab = [0;diff(vbf)];
wab = WhitenSignal(ab);
wvb = WhitenSignal(vbf);
nfft = 2^8;
window = 2^7;
[yw,fw,tw] = mtchglong(wab,nfft,s.xyzSampleRate,window,(1-1/2^2)*window,[],[],[],[0.1 20]);
imagesc(tw,fw,log10(yw)'),axis xy

CheckEegStates('jg05-20120315','theta', {tw+s.syncPeriods(1,1)/s.lfpSampleRate,fw,log10(yw),'imagesc'},[],[],[], 'display', 0)

[yl,fl,tl] = mtchglong(wlfp,nfft,s.lfpSampleRate,window,(1-1/2^2)*window,[],[],[],[0.1 20]);

fmrlfp = ButFilter(mrlfp,2,[5 100]/625,'bandpass');

figure,hist(blfp(:,1,7),400)
figure,plot(cblfp(:,400:600,7))



s = MTASession('jg05-20120315');
t = MTATrial(s,'ctrl1')
lfp  = t.loadlfp(72);
[ys,fs,ts] = spect(t,t.ang(:,5,9,2));
[yl,fl,tl] = spect('lfp',lfp);


figure,
s1 = subplot(211);
imagesc(ts,fs,log10(ys)'),
axis xy,
s2 = subplot(212);
imagesc(tl,fl,log10(yl)'),
caxis([0.2,3.8])
axis xy,
linkaxes([s1,s2],'x')


kern = gausswin(11,0.1);
kern =kern/sum(kern);
vf = Filter0(kern,t.vel(7));
vb = [0;diff(vf)];
vbf = Filter0(kern,vb);
af = [0;diff(vbf);0];
aff = Filter0(kern,af);
afe =abs(hilbert(sqrt(aff.^2)));
afee =abs(hilbert(afe));
amind = LocalMinima(-afee,10);
ae = zeros(size(afee))';
ae(amind) =afee(amind);
aef = Filter0(kern,ae);

for i = 1:length(amind)-2,
 ae(amind(i):amind(i+1)) = (afee(amind(i))+afee(amind(i+1))+afee(amind(i+2)))/3;
end

figure,
s1 = subplot(311);
imagesc(ts,fs,log10(ys)'),
axis xy,
s2 = subplot(312);
imagesc(tl,fl,log10(yl)'),
caxis([0.2,3.8])
axis xy,
s3 = subplot(313);
rp = t.xyz(:,7,3).*t.ang(:,4,5,2).*t.ang(:,3,4,2);
rt = d2t(rp,t.xyzSampleRate,0);
plot(rt,ae)
linkaxes([s1,s3,s2],'x')

ysp = median(ys(:,fs<15),2);
ylp = median(yl(:,fl>4&fl<12),2);

ylp = interp1(1:size(ylp,1),ylp,1:size(ysp,1))';
lpa =  interp1(1:size(ylp),ylp,1:size(aef))';


figure,
p1 = subplot(211);
plot(rt,aef)
p2 = subplot(212);
plot(ts,ylp)
linkaxes([p1,p2],'x')


figure,
p1 = subplot(211);
plot(aef)
p2 = subplot(212);
plot(lpa)
linkaxes([p1,p2],'x')
