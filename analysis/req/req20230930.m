Trial.spk.map% question: a sentence worded or expressed so as to elicit information.
%
% What information is questionable given the available information?
%    Given: timeseries of the local field potential as 4 sets of 64 equally spaced (20um) electrodes in distinct
%    brain regions (LEC, HIP, MEC, OFB) at a sampling rate of 1250Hz.
%    
%    if a network is pushing information there should be a lag
%    between information packets such that the reciever is slightly
%    delayed
%
%        - What is the physiologically expected timelag?
%
%        - does MECL3 gamma correspond to CA1LM gamma
%             - (if so) what is phase shift between activity?
%        
%

%= '/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.lfp'
parpath =  '/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.xml'

channelLEC = 121;
channelsLEC = [121,113,98,85,65]
channelMEC = 305;
channelHIP = 150;

channelsAll = [132,147,159,167,178,182,191,...
               128,119,114, 98, 85, 65,...
               320,311,299,284,271,257,];




lfpAll =  LoadBinary(filepath,channelsAll,par.nChannels,[],[], [],[])';

lfpLRCS = diff(LoadBinary(filepath,[128,114],par.nChannels,[],[], [],[])',1,2);

lfpMRC =  LoadBinary(filepath,[299,311],par.nChannels,[],[], [],[])';
lfpMRCD = diff(lfpMRC,1,2);



lfpHRCD = diff(LoadBinary(filepath,[158,167],par.nChannels,[],[], [],[])',1,2);
lfpHRCS = diff(LoadBinary(filepath,[140,153],par.nChannels,[],[], [],[])',1,2);
lfpHRCM = diff(LoadBinary(filepath,[144,158],par.nChannels,[],[], [],[])',1,2);



lfpRC =  LoadBinary(filepath,[159,167],par.nChannels,[],[], [],[])';
lfpRCD = diff(lfpRC,1,2);

par = LoadPar(parpath);    



figure();
hold('on');
plot(lpfTS,lfpAll(:,8));
plot(lpfTS,lfpAll(:,10));
plot(lpfTS,lfpAll(:,5)+1000);
plot(lpfTS,lfpLRCS(:,1)-1000)
Lines([],-1000,'k');


lfpVC =  LoadBinary(filepath,[129 ,178],par.nChannels,[],[], [],[])';


parspec = struct('nFFT',2^11,...
                 'Fs',  par.lfpSampleRate,...
                 'WinLength',2^10,...
                 'nOverlap',2^10*.5,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[1,30]);
wlfpLEC =  lfpVC;
%wlfpLEC =  WhitenSignal([lfpLEC(:,[1,2,5]),dlfpLEC(:,1)],[],false);
flagCrossSpec = true;
mode = 'mtcsdglong';
[vys,fs,ts] = spec(str2func(mode),wlfpLEC,parspec,flagCrossSpec);
ts = ts+(parspec.WinLength/2)/par.lfpSampleRate;

lfpAll = [lfpAll,lfpMRCD(:,1),lfpLRCD(:,1)];

ind = (log10(cys(:,11,2,2))-log10(cys(:,3,2,2)))>0.5;
frqRng = 5<fs&fs<10;
pmat = []
cmat = [];
for c1 = 1:size(lfpAll,2)-1
    for c2 = (c1+1):size(lfpAll,2)
        tic    
        parspec = struct('nFFT',2^11,...
                         'Fs',  par.lfpSampleRate,...
                         'WinLength',2^10,...
                         'nOverlap',2^10*.5,...
                         'NW',3,...
                         'Detrend',[],...
                         'nTapers',[],...
                         'FreqRange',[1,140]);
        %wlfpLEC =  [lfpHIP(:,[1,4]),lfpMECEND,lfpMECEND(:,1)-lfpMECEND(:,2),-lfpMRCD(:,1),lfpMRC(:,:)];
        %wlfpLEC =  WhitenSignal([lfpLEC(:,[1,2,5]),dlfpLEC(:,1)],[],false);
        wlfpLEC = lfpAll(:,[c1,c2]);
        flagCrossSpec = true;
        mode = 'mtcsdglong';
        [cys,fs,ts] = spec(str2func(mode),wlfpLEC,parspec,flagCrossSpec);
        ts = ts+(parspec.WinLength/2)/par.lfpSampleRate;
        cmat(c1,c2) = sum(mean(cys(ind,frqRng,1,2).*conj(cys(ind,frqRng,1,2))))./sum(mean(cys(ind,frqRng,1,1).*cys(ind,frqRng,2,2)));
        toc
    end
    pmat(c1,:) = log10(mean(cys(ind,:,1,1),1));    
end
pmat(c2,:) = log10(mean(cys(ind,:,2,2),1));    
    

figure,
imagesc([1:size(cmat,1)],[1:size(cmat,2)],cmat')
colormap('jet');
caxis([0,1]);
set(gca,'ytick',[1:size(cmat,2)]);
set(gca,'yticklabel',{'hpc:ca1-orid','hpc:ca1-oris','hpc:ca1-pyr','hpc:ca1-rad','hpc:ca1-lmd','hpc:ca1-lms','hpc:ca1-dg',...
        'lec-l1','lec-l2','lec-l3','lec-l4','lec-l5','lec-l6',...
        'mec-l1','mec-l2','mec-l3','mec-l4','mec-l5','mec-l6',...
        'mce-rc-1','mce-rc-2'});
set(gca,'xtick',[1:size(cmat,1)]);
set(gca,'xticklabel',{'hpc:ca1-orid','hpc:ca1-oris','hpc:ca1-pyr','hpc:ca1-rad','hpc:ca1-lmd','hpc:ca1-lms','hpc:ca1-dg',...
        'lec-l1','lec-l2','lec-l3','lec-l4','lec-l5','lec-l6',...
        'mec-l1','mec-l2','mec-l3','mec-l4','mec-l5','mec-l6',...
        'mce-rc-1','mce-rc-2'});
set(gca,'XTickLabelRotation',90);
cax = colorbar();
ylabel(cax,'Coherence')
title('Inter-Region/Layer Coherence for the Theta [6-10Hz] Frequency Band')



figure,
imagesc(fs,1:size(pmat,1),log10(bsxfun(@rdivide,10.^pmat,sum(10.^pmat))));
colormap('jet');
set(gca,'YTick',[1:size(pmat,1)]);
set(gca,'YTickLabel',...
        {'hpc:ca1-orid','hpc:ca1-oris','hpc:ca1-pyr','hpc:ca1-rad','hpc:ca1-lmd','hpc:ca1-lms','hpc:ca1-dg',...
        'lec-l1','lec-l2','lec-l3','lec-l4','lec-l5','lec-l6',...
        'mec-l1','mec-l2','mec-l3','mec-l4','mec-l5','mec-l6',...
         'mrc-1','mrc-2'});

title({'Mean Spectral Power','(Normalization: columns, Window: 2^{10}, fs:1250)'});
xlabel('Frequency (Hz)');
%set(gca,'XTickLabelRotation',90);


    
thetaVD = [log10(circshift(vys(:,11,1,1),-1))-log10(circshift(vys(:,11,1,1),1)),...
           log10(circshift(vys(:,11,2,2),-1))-log10(circshift(vys(:,11,2,2),1))];



thetaVD = [log10(circshift(cys(:,11,2,2),-1))-log10(circshift(cys(:,11,2,2),1)),...
           log10(circshift(cys(:,11,6,6),-1))-log10(circshift(cys(:,11,6,6),1))];
corr(thetaVD(ind,:))

figure,
subplot(211);polarhistogram(angle(cys(ind,11,1,6)),32);
subplot(212);polarhistogram(angle(cys(ind,11,2,6)),32);


circ_std(angle(cys(ind,11,2,6)))%  lm - DMEC
circ_std(angle(cys(ind,11,2,8)))%  lm - sup
circ_std(angle(cys(ind,11,1,6)))% pyr - DMEC
circ_std(angle(cys(ind,11,2,7)))%  lm - mid
circ_std(angle(cys(ind,11,7,8)))% mid - sup
circ_std(angle(cys(ind,11,1,8)))% pyr - sup

sum(cys(ind,11,2,6).*conj(cys(ind,11,2,6)))./sum(cys(ind,11,6,6).*cys(ind,11,2,2))

sum(cys(ind,11,2,8).*conj(cys(ind,11,2,8)))./sum(cys(ind,11,8,8).*cys(ind,11,2,2))
sum(cys(ind,11,2,7).*conj(cys(ind,11,2,7)))./sum(cys(ind,11,7,7).*cys(ind,11,2,2))
sum(cys(ind,11,1,2).*conj(cys(ind,11,1,2)))./sum(cys(ind,11,1,1).*cys(ind,11,2,2))
sum(cys(ind,11,1,4).*conj(cys(ind,11,1,4)))./sum(cys(ind,11,1,1).*cys(ind,11,4,4))
sum(cys(ind,11,1,6).*conj(cys(ind,11,1,6)))./sum(cys(ind,11,1,1).*cys(ind,11,6,6))


cmat = zeros([8,8]);
frqRng = 5<fs&fs<10;
for i = 1:8
    for j = i+1:8
        cmat(i,j) = sum(mean(cys(ind,frqRng,i,j).*conj(cys(ind,frqRng,i,j))))./sum(mean(cys(ind,frqRng,j,j).*cys(ind,frqRng,i,i)));
    end
end

figure,
imagesc(cmat')
colormap('jet');
caxis([0.5,1]);

figure,
subplot(211);polarhistogram(angle(cys(ind,11,7,8)),32);
subplot(212);polarhistogram(angle(cys(ind,11,7,6)),32);

circ_std(angle(vys(ind,11,1,2)))
circ_std(angle(cys(ind,11,2,6)))

n = 3
figure
chan = 1;    subplot(n,1,chan);    imagesc(ts,fs,log10(vys(:,:,chan,chan))');axis('xy');colormap(gca,'jet');title('Hippocampus PYR');caxis([0,4.5])
stspath = ['/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.sts.REM'];
ds = load(stspath);
Lines(ds(:)./par.lfpSampleRate,[],'m',[],3);
stspath = '/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.sts.RUN';
ds = load(stspath);
Lines(ds(:)./par.lfpSampleRate,[],'k',[],3); 
chan = 2;    subplot(n,1,chan);    imagesc(ts,fs,log10(cys(:,:,4,4))');axis('xy');colormap(gca,'jet');title('Hippocampus LM');caxis([0,4.5])
chan = 3;    subplot(n,1,chan);    imagesc(ts,fs,log10(cys(:,:,6,6))');axis('xy');colormap(gca,'jet');title('Hippocampus LM');caxis([0,4.5])
linkax('xy');


% LEC
lfpLEC = LoadBinary(filepath,channelLEC+[-4,4,-14,-24,-34,-44],par.nChannels,[],[], [],[])';
lfpLEC = LoadBinary(filepath,65:2:128,par.nChannels,[],[], [],[])';
lfpLEC = LoadBinary(filepath,channelsLEC,par.nChannels,[],[], [],[])';
dlfpLEC = diff(lfpLEC,1,2);

lpfTS = [1:size(dlfpLEC,1)]/1250;


ind = 4000001:4015001;
figure,
hold('on');
for chan = 1:size(lfpLEC,2)
    plot(lpfTS(ind),lfpLEC(ind,chan)+chan*300,'k')
end

lfpMECEND = LoadBinary(filepath,[257,320],par.nChannels,[],[], [],[])';

figure
hold('on');
plot(lpfTS,lfpLEC(:,1)+1400,'c');
plot(lpfTS,lfpLEC(:,end)+800,'m');
plot(lpfTS,lfpLEC(:,1)-lfpLEC(:,end),'k')
plot(lpfTS,lfpLEC(:,end-8)-lfpLEC(:,end)-400,'r')
plot(lpfTS,lfpHIP(:,4)-1400,'k')
plot(lpfTS,lfpMECEND(:,1)-2400,'c')
plot(lpfTS,lfpMECEND(:,2)-2800,'m')
plot(lpfTS,lfpMECEND(:,1)-lfpMECEND(:,2)-3200,'m')


figure
hold('on');
plot(lpfTS,-lfpMRCD(:,1)+1400,'c');
plot(lpfTS,lfpMRC(:,1),'k')
plot(lpfTS,lfpMRC(:,2)-200,'m')
plot(lpfTS,lfpHIP(:,4)-1400,'k')

ind = (log10(cys(:,11,2,2))-log10(cys(:,3,2,2)))>0.5;

thetaD = [log10(circshift(cys(:,11,2,2),-1))-log10(circshift(cys(:,11,2,2),1)),...
          log10(circshift(cys(:,11,6,6),-1))-log10(circshift(cys(:,11,6,6),1))];

figure,
ind = (log10(cys(:,11,2,2))-log10(cys(:,3,2,2)))>0.5;
subplot(211);
hist2(thetaD(ind,:),linspace(-1,1,50),linspace(-1,1,50))
subplot(212);
thetaDS = [circshift(thetaD(:,1),1),circshift(thetaD(:,2),0)];
hist2(thetaDS(ind,:),linspace(-1,1,50),linspace(-1,1,50));

corr(thetaD(ind,:))


figure,
subplot(211);
hist2(log10([cys(ind,11,2,2),cys(ind,11,6,6)]),linspace(3,4.2,50),linspace(2,4.2,50))
subplot(212);
stp = circshift(cys(:,11,2,2),-4);
hist2(log10([stp(ind),cys(ind,11,6,6)]),linspace(3,4.2,50),linspace(2,4.2,50))

n = 9
figure
chan = 1;    subplot(n,1,chan);    imagesc(ts,fs,log10(cys(:,:,chan,chan))');axis('xy');colormap(gca,'jet');title('Hippocampus PYR');caxis([0,4])
stspath = ['/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.sts.REM'];
ds = load(stspath);
Lines(ds(:)./par.lfpSampleRate,[],'m',[],3);
stspath = '/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.sts.RUN';
ds = load(stspath);
Lines(ds(:)./par.lfpSampleRate,[],'k',[],3); 
chan = 2;    subplot(n,1,chan);    imagesc(ts,fs,log10(cys(:,:,chan,chan))');axis('xy');colormap(gca,'jet');title('Hippocampus LM');caxis([0,4])
chan = 3;    subplot(n,1,chan);    imagesc(ts,fs,log10(cys(:,:,chan,chan))');axis('xy');colormap(gca,'jet');title('MEC Deep');caxis([0,4])
chan = 4;    subplot(n,1,chan);    imagesc(ts,fs,log10(cys(:,:,chan,chan))');axis('xy');colormap(gca,'jet');title('MEC Sup');caxis([0,4])
chan = 5;    subplot(n,1,chan);    imagesc(ts,fs,log10(cys(:,:,chan,chan))');axis('xy');colormap(gca,'jet');title('dMECSHANK');caxis([0,4])
chan = 6;    subplot(n,1,chan);    imagesc(ts,fs,log10(cys(:,:,chan,chan))');axis('xy');colormap(gca,'jet');title('dMEC');caxis([0,4])
chan = 7;    subplot(n,1,chan);    imagesc(ts,fs,log10(cys(:,:,chan,chan))');axis('xy');colormap(gca,'jet');title('MEC_MID');caxis([0,4])
chan = 8;    subplot(n,1,chan);    imagesc(ts,fs,((cys(:,:,1,5).*conj(cys(:,:,1,5)))./(cys(:,:,5,5).*cys(:,:,1,1)))');...
                                                                             axis('xy');colormap(gca,'jet');title('Hippocampus Pyr - dMEC');
chan = 9;    subplot(n,1,chan);    imagesc(ts,fs,((cys(:,:,2,5).*conj(cys(:,:,2,5)))./(cys(:,:,2,2).*cys(:,:,5,5)))');...
                                                                             axis('xy');colormap(gca,'jet');title('Hippocampus LM - dMEC');
linkax('xy');colormap(gca,'jet');

figure;
subplot(511);
imagesc(ts,fs,(log10(cys(:,:,6,6))-log10(cys(:,:,7,7)))');colormap(gca,'jet');axis('xy');caxis([-0.25,1]);
subplot(512);
imagesc(ts,fs,log10(cys(:,:,6,6))');colormap(gca,'jet');axis('xy');caxis([0,4])
subplot(513);
imagesc(ts,fs,log10(cys(:,:,7,7))');colormap(gca,'jet');axis('xy');caxis([0,4])
subplot(514);
imagesc(ts,fs,(log10(cys(:,:,7,7))-log10(cys(:,:,3,3)))');colormap(gca,'jet');axis('xy');caxis([-0.25,1]);
subplot(515);
imagesc(ts,fs,angle(cys(:,:,7,8))');colormap(gca,'hsv');axis('xy');
linkax('xy');
ForAllSubplots('axis(''xy'');colorbar();');


figure;
subplot(411);
imagesc(ts,fs,angle(cys(:,:,1,2))');colormap(gca,'hsv');axis('xy');
subplot(412);
imagesc(ts,fs,angle(cys(:,:,2,6))');colormap(gca,'hsv');axis('xy');
subplot(413);
imagesc(ts,fs,angle(cys(:,:,1,3))');colormap(gca,'hsv');axis('xy');
subplot(414);
imagesc(ts,fs,circ_dist(angle(cys(:,:,3,4)),angle(cys(:,:,2,5)))');colormap(gca,'hsv');axis('xy');
linkax('xy');
ForAllSubplots('axis(''xy'');colorbar();');



figure,
hold('on')
plot(lfpLEC(:,[1:2]),'b')
plot(-lfpLEC(:,3:end)+1000,'k')
plot(dlfpLEC(:,1),'r')

plot(dlfpLEC(:,1)+1000,'r')
plot(-lfpLEC(:,3)+1000,'c')



parspec = struct('nFFT',2^11,...
                 'Fs',  par.lfpSampleRate,...
                 'WinLength',2^10,...
                 'nOverlap',2^10*.5,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[1,140]);
wlfpLEC =  [lfpLEC(:,[1,2,5]),dlfpLEC(:,1)];
wlfpLEC =  WhitenSignal([lfpLEC(:,[1,2,5]),dlfpLEC(:,1)],[],false);
flagCrossSpec = true;
mode = 'mtcsdglong';
[lys,fs,ts] = spec(str2func(mode),wlfpLEC,parspec,flagCrossSpec);
ts = ts+(parspec.WinLength/2)/par.lfpSampleRate;


n = 6
figure
chan = 1;    subplot(n,1,chan);    imagesc(ts,fs,log10(lys(:,:,chan,chan))');axis('xy');colormap(gca,'jet');
stspath = ['/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.sts.REM'];
ds = load(stspath);
Lines(ds(:)./par.lfpSampleRate,[],'m',[],3);
stspath = '/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.sts.RUN';
ds = load(stspath);
Lines(ds(:)./par.lfpSampleRate,[],'k',[],3);

chan = 2;    subplot(n,1,chan);    imagesc(ts,fs,log10(lys(:,:,chan,chan))');axis('xy');colormap(gca,'jet');
chan = 3;    subplot(n,1,chan);    imagesc(ts,fs,log10(lys(:,:,chan,chan))');axis('xy');colormap(gca,'jet');
chan = 4;    subplot(n,1,chan);    imagesc(ts,fs,log10(lys(:,:,chan,chan))');axis('xy');colormap(gca,'jet');
chan = 5;    subplot(n,1,chan);    imagesc(ts,fs,((lys(:,:,3,4).*conj(lys(:,:,3,4)))./(lys(:,:,4,4).*lys(:,:,3,3)))');axis('xy');colormap(gca,'jet');
chan = 6;    subplot(n,1,chan);    imagesc(ts,fs,((lys(:,:,1,3).*conj(lys(:,:,1,3)))./(lys(:,:,1,1).*lys(:,:,3,3)))');axis('xy');colormap(gca,'jet');
linkax('xy');






n = 4
figure,
chan = 1;    subplot(n,1,1);    imagesc(ts,fs,log10(lys(:,:,chan, chan))');
chan = 3;    subplot(n,1,2);    imagesc(ts,fs,log10(lys(:,:,chan,chan))');axis('xy');
%chan = 1;    subplot(n,1,3);    imagesc(ts,fs,log10(lys(:,:,chan,chan)+10000)'./log10(lys(:,:,4,4)+10000)');
chan = 4;    subplot(n,1,3);    imagesc(ts,fs,log10(lys(:,:,chan,chan))');axis('xy');
ForAllSubplots('axis(''xy'');colormap(''jet'');');
stspath = ['/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.sts.REM'];
ds = load(stspath);
Lines(ds(:)./par.lfpSampleRate,[],'m',[],3);
stspath = '/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.sts.RUN';
ds = load(stspath);
Lines(ds(:)./par.lfpSampleRate,[],'k',[],3);
chan = 1;    subplot(n,1,4);    imagesc(ts,fs,angle(lys(:,:,1,2)'));axis('xy');colormap(gca,'hsv');
linkax('xy');
ForAllSubplots('colorbar()');










% MEC
lfpMEC = LoadBinary(filepath,channelMEC+[-4,4,-14,-24,-34,-44],par.nChannels,[],[], [],[])';
dlfpMEC = diff(lfpMEC,1,2);

figure,
hold('on')
plot(lfpMEC(:,[1:2]),'b')
plot(-lfpMEC(:,3:end)+1000,'k')
plot(dlfpMEC(:,1),'r')

plot(dlfpMEC(:,1)+1000,'r')
plot(-lfpMEC(:,3)+1000,'c')




parspec = struct('nFFT',2^11,...
                 'Fs',  par.lfpSampleRate,...
                 'WinLength',2^10,...
                 'nOverlap',2^10*.5,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[1,140]);
wlfpMEC =  [lfpMEC(:,[1,2,5]),dlfpMEC(:,1)];
%wlfpMEC =  WhitenSignal([lfpMEC(:,[1,2,5]),dlfpMEC(:,1)],[],true);
flagCrossSpec = true;
mode = 'mtcsdglong';
[mys,fs,ts] = spec(str2func(mode),wlfpMEC,parspec,flagCrossSpec);
ts = ts+(parspec.WinLength/2)/par.lfpSampleRate;



n = 6
figure
chan = 1;    subplot(n,1,chan);    imagesc(ts,fs,log10(mys(:,:,chan,chan))');
chan = 2;    subplot(n,1,chan);    imagesc(ts,fs,log10(mys(:,:,chan,chan))');axis('xy');
chan = 3;    subplot(n,1,chan);    imagesc(ts,fs,log10(mys(:,:,chan,chan))');axis('xy');
chan = 4;    subplot(n,1,chan);    imagesc(ts,fs,log10(mys(:,:,chan,chan))');axis('xy');
chan = 5;    subplot(n,1,chan);    imagesc(ts,fs,((mys(:,:,3,4).*conj(mys(:,:,3,4)))./(mys(:,:,4,4).*mys(:,:,3,3)))');axis('xy');
chan = 6;    subplot(n,1,chan);    imagesc(ts,fs,((mys(:,:,1,3).*conj(mys(:,:,1,3)))./(mys(:,:,1,1).*mys(:,:,3,3)))');axis('xy');
linkax('xy');
ForAllSubplots('axis(''xy'');colormap(''jet'');');





channelHIP+[-6,4,11,21,29,40]

channelsHIP = [147,151,163,178,182];
% LEC
lfpHIP = LoadBinary(filepath,channelsHIP,par.nChannels,[],[], [],[])';
lfpHIP = LoadBinary(filepath,channelHIP+[-4,4,11,21,29,40],par.nChannels,[],[], [],[])';
dlfpHIP = diff(lfpHIP,1,2);
% $$$ 
% $$$ figure,
% $$$ hold('on')
% $$$ plot(lfpHIP(:,[1:2]),'b')
% $$$ plot(-lfpHIP(:,3:end)+1000,'k')
% $$$ plot(dlfpHIP(:,1),'r')
% $$$ 
% $$$ plot(dlfpHIP(:,1)+1000,'r')
% $$$ plot(-lfpHIP(:,3)+1000,'c')



parspec = struct('nFFT',2^11,...
                 'Fs',  par.lfpSampleRate,...
                 'WinLength',2^10,...
                 'nOverlap',2^10*.5,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[1,140]);
wlfpHIP =  [lfpHIP(:,[1,2,4]),dlfpHIP(:,1)];
%wlfpHIP =  WhitenSignal([lfpHIP(:,[1,2,5]),dlfpHIP(:,1)],[],true);
flagCrossSpec = true;
mode = 'mtcsdglong';
[hys,fs,ts] = spec(str2func(mode),wlfpHIP,parspec,flagCrossSpec);
ts = ts+(parspec.WinLength/2)/par.lfpSampleRate;


wlfpHIP =  [lfpHIP(:,[1,2,4]),dlfpHIP(:,1)*2];
figure,
hold('on');
cclr = 'bbcr';
for chan = 1:4
plot([1:size(wlfpHIP,1)]./par.lfpSampleRate,wlfpHIP(:,chan)-chan*1000,['-',cclr(chan)]);
end


n = 6
figure
chan = 1;    subplot(n,1,chan);    imagesc(ts,fs,log10(hys(:,:,chan,chan))');
chan = 2;    subplot(n,1,chan);    imagesc(ts,fs,log10(hys(:,:,chan,chan))');axis('xy');
chan = 3;    subplot(n,1,chan);    imagesc(ts,fs,log10(hys(:,:,chan,chan))');axis('xy');
chan = 4;    subplot(n,1,chan);    imagesc(ts,fs,log10(hys(:,:,chan,chan))');axis('xy');
chan = 5;    subplot(n,1,chan);    imagesc(ts,fs,((hys(:,:,3,4).*conj(hys(:,:,3,4)))./(hys(:,:,4,4).*hys(:,:,3,3)))');axis('xy');
chan = 6;    subplot(n,1,chan);    imagesc(ts,fs,((hys(:,:,1,3).*conj(hys(:,:,1,3)))./(hys(:,:,1,1).*hys(:,:,3,3)))');axis('xy');
linkax('xy');



chan = 1;    subplot(n,1,1);    imagesc(ts,fs,log10(hys(:,:,chan,chan))');

n = 4
figure,
chan = 1;    subplot(n,1,1);    imagesc(ts,fs,log10(hys(:,:,chan, chan))');caxis([0.25,4.25]);
chan = 3;    subplot(n,1,2);    imagesc(ts,fs,log10(hys(:,:,chan,chan))');axis('xy');caxis([0.25,4.25]);
chan = 1;    subplot(n,1,3);    imagesc(ts,fs,log10(hys(:,:,chan,chan)+10000)'./log10(hys(:,:,4,4)+10000)');caxis([0.98,1.05]);
ForAllSubplots('axis(''xy'');colormap(''jet'');');

stspath = ['/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.sts.REM'];
ds = load(stspath);
Lines(ds(:)./par.lfpSampleRate,[],'m',[],3);
stspath = '/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.sts.RUN';
ds = load(stspath);
Lines(ds(:)./par.lfpSampleRate,[],'k',[],3);

stspath = '/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.sts.RUN';
ds = load(stspath);
Lines(ds(:)./par.lfpSampleRate,[],'k',[],3);

chan = 1;    subplot(n,1,4);    imagesc(ts,fs,angle(hys(:,:,1,2)'));axis('xy');colormap(gca,'hsv');
linkax('xy');

ForAllSubplots('colorbar()');






channelsGamma = [299,180];
lfpGamma = LoadBinary(filepath,channelsGamma,par.nChannels,[],[],[],[])';

wlfpGamma =  WhitenSignal(lfpGamma,[],true);
parspec = struct('nFFT',2^7,...
                 'Fs',  par.lfpSampleRate,...
                 'WinLength',2^8,...
                 'nOverlap',2^8*.5,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[20,140]);

[gys,gfs,gts] = spec(str2func(mode),wlfpLEC,parspec,flagCrossSpec);
gts = gts+(parspec.WinLength/2)/par.lfpSampleRate;


n = 5;
figure
chan = 1;    subplot(n,1,chan);    imagesc(gts,gfs,log10(gys(:,:,chan,chan))');axis('xy');colormap('jet');
chan = 2;    subplot(n,1,chan);    imagesc(gts,gfs,log10(gys(:,:,chan,chan))');axis('xy');colormap('jet');
chan = 3;    subplot(n,1,chan);    imagesc(gts,gfs,((gys(:,:,1,2).*conj(gys(:,:,1,2)))./(gys(:,:,1,1).*gys(:,:,2,2)))');axis('xy');colormap('jet');
chan = 4;    subplot(n,1,chan);    imagesc(gts,gfs,circ_dist(angle(gys(:,:,1,2))',pi));axis('xy');colormap(gca(),'hsv');
chan = 5;    subplot(n,1,chan);    imagesc(ts,fs,log10(hys(:,:,4, ...
                                                  4))');axis('xy');colormap(gca(),'jet');



figure,
imagesc(ts,fs,log10(hys(ind,:,4,4))');axis('xy');colormap(gca(),'jet');

hold('on');
plot(lpfTS,lfpGamma(:,1)/50+40,'k')
linkax('xy');


flfpGamma = ButFilter(lfpGamma,4,[6,11]/(par.lfpSampleRate/2),'bandpass');
phzGamma = phase(hilbert(flfpGamma(:,1)));




pts = round(gts.*par.lfpSampleRate);
gysThetaPhz = phzGamma(pts);

[~,gind] = NearestNeighbour(gts,ts(ind));


gtPhzInd = discretize(gysThetaPhz(gind),linspace(-pi,pi,32));
gtG = gys(gind,:,:,:);

my_cohere =  @(y) sum(y(:,1,2).*conj(y(:,1,2)))./sum(y(:,1,1).*y(:,2,2));
out = [];
for f = 1:numel(gfs)
    for p = 1:31
        out(p,f) = my_cohere(sq(gtG(gtPhzInd==p,f,:,:)));
    end
end


figure,
imagesc(linspace(-pi,pi,32),gfs,out')
axis('xy');
colormap('jet');
caxis([0.5,0.8])


% next get phz of lfp and bin coherence by theta phase