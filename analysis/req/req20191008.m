;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.

MjgER2016_load_data();

Trial = Trials{20};
units = units{20}
sampleRate = 250;
nanColor = [1,1,1];
xyz = Trial.load('xyz','trb');
xyz.resample(250);
spk = create(copy(Trial.spk),Trial, sampleRate);
spikeWindow = 0.1;

pyrFr = Trial.load('ufr',xyz,spk,select_units(Trial,'pyr'),spikeWindow,true,'gauss');
intFr = Trial.load('ufr',xyz,spk,select_units(Trial,'int'),spikeWindow,true,'gauss');
pyrFr.data(pyrFr.data<0.2) = nan;
intFr.data(intFr.data<0.2) = nan;

stc = Trial.load('stc','msnn_ppsvd_raux');

dc2.ufr = cf(@(t,x,u,s)                                                         ...
              t.load('ufr',x,s,u,spikeWindow,true,'gauss'),                     ...
              Trials, dc2.xyz, units, dc2.spk);


channels = 66:80;
lfp = Trial.load('lfp',channels);
freqRange = [180,280];
lfpld = resample(copy(lfp),xyz);

channels = 57:64;
buz = Trial.load('lfp',channels);


ccr = Trial.load('lfp',35:8:64);
ccrlh = filter(copy(ccr),'ButFilter',3,[180,220],'bandpass');
ccrlb = Shilbert(ccrlh.data);
ccrlp = zeros([size(ccrlb,1),1]);
for t = 1:size(ccrlb,1),
    ccrlp(t) = PPC(angle(ccrlb(t,:)));
end


xyzlf = filter(copy(xyz),'ButFilter',3,[30],'low');
vxyz = vel(xyzlf,'head_back',[1,2,3]);
vxy = vel(xyzlf,'head_back',[1,2]);

flfp = filter(copy(lfp),'ButFilter',3,[180,280],'bandpass');
rpow = sq(sum(flfp.segs(1:7:size(flfp,1),51,0).^2));
lfpfh = filter(copy(lfp),'ButFilter',3,[300,350],'bandpass');
rpowh = sq(sum(lfpfh.segs(1:7:size(flfp,1),51,0).^2));

ts = ([1:7:size(lfp,1)]+26)'./lfp.sampleRate;

buzlfh = filter(copy(buz),'ButFilter',3,[300,350],'bandpass');
buzlf = filter(copy(buz),'ButFilter',3,[180,280],'bandpass');
bpow = sq(sum(buzlf.segs(1:7:size(buzlf,1),51,0).^2));
bpowh = sq(sum(buzlf.segs(1:7:size(buzlf,1),51,0).^2));
buzld = resample(copy(buz),xyz);

uinc = sum(pyrFr.data>0.2,2,'omitnan');

figure,plot
clrs = jet(11);
figure,hold('on');
for c = 1:11,
    plot(ts,rpow(:,c),'Color',clrs(c,:))
end

clrs = cool(15);
figure,
    sp = tight_subplot(6,1,0.01,0.1);
    s = 1;
    axes(sp(s));s=s+1;
    %imagesc(ts,[1:8,1:15],[nunity(bpow,[],[],[],[5,95]),nunity(rpow,[],[],[],[5,95])]');
        imagesc(ts,[1:8,1:15],[log10(bpow),log10(rpow)]');
        axis('ij');
        colormap('jet');
        caxis([-1,8]);
    axes(sp(s));s=s+1;
        imagesc(ts,[1:8,1:15],[log10(bpowh),log10(rpowh)]');
        %imagesc(ts,[1:8,1:15],[nunity(bpowh,[],[],[],[5,95]),nunity(rpowh,[],[],[],[5,95])]');
        axis('ij');
        colormap('jet');
        caxis([-1,8]);
    axes(sp(s));s=s+1;    
        hold('on');
        for c = 1:15;
            plot([1:size(lfpld,1)]./lfpld.sampleRate,lfpld(:,c)-c*2400,'Color',clrs(c,:));
        end
        for c = 1:8;
            plot([1:size(buzld,1)]./buzld.sampleRate,buzld(:,c)-c*2400+20000,'Color',clrs(c,:));
        end
        plot([1:size(vxyz,1)]./vxyz.sampleRate,vxyz.data.*1e3+20000,'r');
        plot([1:size(vxy,1)]./vxy.sampleRate,vxy.data.*1e3+20000,'k');        
        plot([1:size(pyrFr.data)]./pyrFr.sampleRate,uinc*1e3+20000,'c'),
        %plot([1:size(lfp.data)]./lfp.sampleRate,ccrlp*1e4-70000,'c'),
        
    axes(sp(s));s=s+1;
        imagescnan({[1:size(pyrFr.data)]./pyrFr.sampleRate,1:size(pyrFr,2),pyrFr.data'},...
                   [0,20],'nanRGB',nanColor);
    axes(sp(s));s=s+1;    
        imagescnan({[1:size(intFr.data)]./intFr.sampleRate,1:size(intFr,2),intFr.data'},...
                   [0,20],'nanRGB',nanColor);
    axes(sp(s));s=s+1;
        plot_stc(stc,1);
linkaxes(findobj(gcf,'Type','Axes'),'x');


cpow = clip(nunity(rpow),-1,40)+2;
cpow = (mean(cpow(:,1:4),2,'omitnan')+mean(cpow(:,9:12),2,'omitnan'))./(2*mean(cpow(:,5:8),2,'omitnan'));


% NORMALIZE ripple power
ripplePower = nunity(mean(mean(log10(ysg(:,fsg>freqRange(1)&fsg<freqRange(2),:)),2),3),...
                     @nan,...
                     ripMean,...
                     ripStd);
