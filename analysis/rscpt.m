
s = MTASession('jg05-20120315');

ymas = {};
tmas = {};
fmas = {};
frange = [1 120; 1 40; 50 100; 100 150;150 200];
nffts = [2^11;2^11;2^9;2^8;2^7];
windows =[2^10;2^10;2^8;2^7;2^6];

rearings = s.xyz(:,s.Model.gmi('head_back'),3).*s.ang(:,s.Model.gmi('spine_middle'),s.Model.gmi('spine_upper'),2).*s.ang(:,s.Model.gmi('spine_upper'),s.Model.gmi('head_back'),2);
rear = zeros(size(rearings)); rear(rearings<44) = 0; rear(rearings>=44) = 1;
rper = ThreshCross(rear,0.5,60);
rper = rper(5:end-5,:);

lfp = LoadBinary([s.path.nlx s.name '/' s.name '.lfp'],[65:96],96,[],[],[],[s.syncPeriods(1,1), s.syncPeriods(end,2)])';
wlfp = WhitenSignal(lfp);


for b = 1:2,
    rlfp  = GetSegs(lfp ,round((rper(:,b)./s.xyzSampleRate-3).*1250),6*1250,[]);
    rwlfp = GetSegs(wlfp,round((rper(:,b)./s.xyzSampleRate-3).*1250),6*1250,[]);
    crlfp = rlfp(:,max(rlfp(:,:,7),[],1)<17800,:);
    crwlfp = rwlfp(:,max(rlfp(:,:,7),[],1)<17800,:);
    for k = 1:5,
        y = [];
        t = [];
        f = [];
        for j = 1:32;
            for i = 1:size(crwlfp,2),
                [y(:,:,i,j),f,t] = mtchglong(crwlfp(:,i,j),nffts(k),s.lfpSampleRate,windows(k),0.75*windows(k),[],[],[],frange(k,:));
            end
        end
        ymas{k,b}=y; 
        tmas{k,b}=t; 
        fmas{k,b}=f; 
    end
end


crper = rper(max(rlfp(:,:,7),[],1)<17800,:);
rears{1} = GetSegs(rearings,round(crper(:,1)-3*s.xyzSampleRate),round(6*s.xyzSampleRate),[]);
rears{2} = GetSegs(rearings,round(crper(:,2)-3*s.xyzSampleRate),round(6*s.xyzSampleRate),[]);

yd = sq(sum(y(t>-1&t<-.25,:,:,:),1)-sum(y(t>.25&t<1,:,:,:),1));
ydm = sq(median(yd,2));
figure,imagesc(f,1:32,ydm');


figure,hist(sq(mean(yd(f<12&f>6,:,20),1))')
ym = sq(mean(y,3));
figure,imagesc(t+diff(t(1:2))/2-3,f,log10(ym'));axis xy

mrlfp=sq(mean(crlfp,2));
fmrlfp = ButFilter(mrlfp,2,[5 100]/625,'bandpass');
PlotCSD(mrlfp,[1:size(mrlfp,1)]/1250,[],2)
imagesc(1:size(rlfp,1),1:size(rlfp,2),rlfp(:,rind,20)');axis xy



