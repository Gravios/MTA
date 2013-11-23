function groom_state = groom(Session,varargin)


sl

head = s.transformOrigin('head_back','head_front',{'head_left','head_right'});
fet = head.roll;
fet(isnan(fet)) = 0;
fetw = WhitenSignal(fet);
tfet = d2t(fet,s.xyzSampleRate,0);

nfft = 2^9;
window = 2^8;
[yw,fw,tw] = mtchglong(fetw,nfft,s.xyzSampleRate,window,(1-1/2^2)*window,[],[],[],[0.1 30]);

fetpow = sum(yw(:,fw<2),2);


figure,
sp1 = subplot(311);
imagesc(tw+0.5*diff(tw(1:2)),fw,log10(yw)'),axis xy
sp2 = subplot(312);
plot(tw,fetpow);
sp3 = subplot(313);
plot(tfet,fet);
linkaxes([sp1,sp2,sp3],'x');

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


lfp = LoadBinary([s.spath.nlx s.name '.lfp'],78,96,[],[],[],[s.syncPeriods(1,1), s.syncPeriods(end,2)])';
wlfp = WhitenSignal(lfp);
nfft = 2^10;
window = 2^9;
[yl,fl,tl] = mtchglong(wlfp,nfft,s.lfpSampleRate,window,(1-1/2^2)*window,[],[],[],[0.1 20]);
%fmrlfp = ButFilter(mrlfp,2,[5 100]/625,'bandpass');


thpow = sum(yl(:,fl>5&fl<12),2);
fetpow = sum(yw(:,fw>2&fw<8),2);
plot(tl,thpow,tw,fetpow*10000000*5);

walkPeriods = s.Bhv.getState('walk').state;
for i = 1:size(walkPeriods,1),
