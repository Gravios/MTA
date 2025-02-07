% question: a sentence worded or expressed so as to elicit information.
%
% What hyptheseses are possible given the available information?
%    Given: timeseries of the local field potential as 4 sets of 64 equally spaced (20um) electrodes in distinct
%    brain regions (LEC, HIP, MEC, OFB) at a sampling rate of 1250Hz.
%    
%    if a network is pushing information there should be a lag
%    between information packets such that the reciever is slightly
%    delayed.
%        - What is the physiologically expected timelag?
%        - does MECL3 gamma correspond to CA1LM gamma
%             - (if so) what is phase shift between activity?
%        
%    Diffusion of charge.
%        - active? which channels?
%        - ion pumps?
%        - channel leakiness dependent upon membrane potential?
%
%
% *** Phase differences between local channels in the LM
%        - influenced by diffusion from radiatum and dentate???
%        - plot local lm theta phase differences as a function of
%          loacl radiatum phase differences.

path_lfp = '/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.lfp'
path_par = '/storage/share/IF/data/processed/ephys/IF14/IF14-20190717a/IF14-20190717a.xml'
path_xyz = '/storage/share/IF/data/processed/xyz/IF14/IF14-20190717a/IF14-20190717a.XYZ.mat';

par = LoadPar(path_par);    
sampleRate = par.lfpSampleRate;

channels_hpc = [132:2:191];
lfp_hpc =  LoadBinary(path_lfp, channels_hpc, par.nChannels,[],[],[],[])';
ts_lfp = [1:size(lfp_hpc,1)]'/sampleRate;
ts_xyz = ts_lfp;

xyz = load(path_xyz);

fvxy = ButFilter(sq(xyz.XYZ(:,1,[1,2])),4, 2.4./(sampleRate/2), 'low');
fvxy = log10(sqrt((sum(((circshift(fvxy,-1)-circshift(fvxy,1)).* (sampleRate/2)).^2,2))));
fvxy = interp1(ts_xyz, fvxy, lts);


channelLEC = 121;
channelsLEC = [121,113,98,85,65]
channelMEC = 305;
channelHIP = 150;

channelsAll = [132,147,159,167,178,182,191,...
               128,119,114, 98, 85, 65,...
               320,311,299,284,271,257,];

lfpHpcD =  diff(LoadBinary(filepath,channelsHpc,par.nChannels,[],[], [],[])',1,2);
lfpHpcDD =  diff(LoadBinary(filepath,152:2:168,par.nChannels,[],[], [],[])',1,2);

channelsLec = [65:2:128]
lfpLec =  LoadBinary(filepath,channelsLec,par.nChannels,[],[], [],[])';

flfpHpc = ButFilter(lfpHpc, 4, [6,11]./(1250/2),'bandpass');

dlHpc = flfpHpc


fwin = gausswin(2^10);
fwin = fwin./sum(fwin);
lpfet = ButFilter(lfp_hpc, 4, [20]./(1250/2),'low');

for chan = 1:30,
    lpfet(:,chan) = log10(conv(lpfet(:,chan).*lpfet(:,chan),fwin,'same'));
end

figure,imagesc(lpfet(1000001:2800001 ,:)')




figure();
hold('on');
plot(flfpHpc(lind,10),'b');
plot(diff(flfpHpc(lind,[8,11]),1,2)*3,'r');
plot(-diff(flfpHpc(lind,[27,29]),1,2),'m');

plot(diff(dlfp(:,[9,11]),1,2)-1000);
plot(diff(dlfp(:,[9,11]),1,2)*300);

lind = [ 1000001:2800001 ];
figure();hold('on');
plot(bsxfun(@minus,lfpHpc(lind,:),[1:30]*200)+8000,'b');
lind = ':';

plot(bsxfun(@minus,flfpHpc(lind,:),[1:30]*200)+8000,'r');
plot(flfpHpc(lind,10)+500,'b');
plot(diff(flfpHpc(lind,[8,11]),1,2)*3+500,'r');
plot(diff(lfpHpc(lind,[27,29]),1,2)+500,'m');
plot(diff(lfpHpc(lind,[18,21]),1,2)+500,'k');
plot(diff(lfpHpc(lind,[14,16]),1,2)*3-500,'k');

lfpHpcRC = [diff(lfp_hpc(:,[8,11]),1,2),...
            diff(lfpHpc(:,[14,16]),1,2),...
            diff(lfpHpc(:,[18,21]),1,2),...
            diff(lfpHpc(:,[27,29]),1,2)];

lfpHpcRC = [diff(lfp_hpc(:,[1,4]),1,2),...
            diff(lfpHpc(:,[8,11]),1,2),...
            diff(lfpHpc(:,[18,21]),1,2),...
            diff(lfpHpc(:,[27,29]),1,2)];

dlfp = diff(lfpHpc(lind,:));
dlfp = ButFilter(dlfp, 4, [6,11]./(1250/2),'bandpass');
dlfp = sq(sqrt(sum(GetSegs(dlfp.^2, 1:1:size(dlfp,1), 1250,0))));



figure();
imagesc(log10(bsxfun(@rdivide,log10(dlfp),sum(log10(dlfp)))'));
caxis([1,3]);
colormap('jet');

parspec = struct('nFFT',2^11,...
                 'Fs',  par.lfpSampleRate,...
                 'WinLength',2^10,...
                 'nOverlap',2^10*0.5,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[1,40]);
lind = [12000001:12800001];
lind = [ 1000001:2800001 ];
lind = ':';
%wlfpLEC = [lfpHpcDD(lind,:),lfpAll(lind,[2,5]),lfpMCD(lind),lfpHRCM(lind)];
wlfpLEC =  [lfpHpcRC(lind,:),lfpHpc(lind,[10]), dlfpHpcRC(lind,:)];
%wlfpLEC =  WhitenSignal(wlfpLEC,[],true);
flagCrossSpec = true;
mode = 'mtcsdglong';
[lys,lfs,lts] = spec(str2func(mode), wlfpLEC, parspec, flagCrossSpec);
lts = lts + (parspec.WinLength/2) / par.lfpSampleRate;


figure,
dlfpHpcRC = diff(lfpHpcRC(:,1:4),1,2);
dlfpHpcRC = diff(dlfpHpcRC(:,:),1,2);
figure,plot(dlfpHpcRC)

tdr =  (log10(lys(:,12,5,5))-log10(lys(:,2,5,5))) > 0.5;

figure(); hold('on');
for chan = 1:7;
subplot2(2,7,1,chan);
ind = fvxy > 0.1 & tdr;
binscatter(log10(lys(ind,3,chan,chan)), log10(lys(ind,12,chan,chan)),...
           'XLimits', [0,4.5], 'YLimits', [0,4.5]);
colormap(gca,'jet');grid(gca,'on');set(gca,'ColorScale','log');
line([0,4.5],[0,4.5],'color','k');
subplot2(2,7,2,chan);
ind = fvxy < 0.1 & tdr;
binscatter(log10(lys(ind,3,chan,chan)), log10(lys(ind,12,chan,chan)),...
           'XLimits', [0,4.5], 'YLimits', [0,4.5]);
colormap('jet');grid(gca,'on');set(gca,'ColorScale','log');
line([0,4.5],[0,4.5],'color','k');
linkxy();
end


figure,
subplot(211);
ind = fvxy > 0.1 & ~tdr;
hist2([log10(lys(ind,12,7,7)),angle(lys(ind,12,6,7))],...
      linspace(1.5,4.6,20),linspace(-pi,pi,20));
subplot(212);
ind = fvxy < 0.1 & ~tdr;
hist2([log10(lys(ind,12,7,7)),angle(lys(ind,12,6,7))],...
      linspace(1.5,4.6,20),linspace(-pi,pi,20));


figure();
subplot(221);
    ind = fvxy > 0.1 & tdr;
    hist2([log10(lys(ind,12,6,6)),log10(lys(ind,2,6,6))],...
          linspace(1.5,3.6,20),linspace(1.5,3.6,20));
subplot(222);
    ind = fvxy > 0.1 & tdr;
    hist2([log10(lys(ind,12,7,7)),log10(lys(ind,2,7,7))],...
          linspace(1.5,4.6,20),linspace(1.5,4.6,20));
subplot(223);
    ind = fvxy > 0.1 & ~tdr;
    hist2([log10(lys(ind,12,6,6)),log10(lys(ind,2,6,6))],...
          linspace(1.5,3.6,20),linspace(1.5,3.6,20));
subplot(224);
    ind = fvxy > 0.1 & ~tdr;
    hist2([log10(lys(ind,12,7,7)),log10(lys(ind,2,7,7))],...
          linspace(1.5,4.6,20),linspace(1.5,4.6,20));


figure();hold('on');
plot(log10(lys(:,12,7,7))-3,fvxy,'.');

chan = 7;
figure();hold('on');
plot(log10(lys(:,14,chan,chan))-log10(lys(:,1,chan,chan)), fvxy, '.');
xlim([-1,3]);
ylim([-3,3]);


chan = 3;
figure();hold('on');
%plot(log10(lys(:,14,chan,chan))-log10(lys(:,1,chan,chan)), fvxy, '.');
binscatter(log10(lys(:,14,chan,chan))-log10(lys(:,1,chan,chan)), fvxy);
set(gca,'ColorScale','log');
colormap(gca,'jet');
xlim(gca,[-1,3]);
ylim(gca,[-3,3]);


ind = ':';
ind = ~tdr;
figure();hold('on');
binscatter(log10(lys(ind,14,5,5))-log10(lys(ind,14,2,2)),...
     log10(lys(ind,14,3,3))-log10(lys(ind,1,3,3)));
xlim([-1,2.5]);
ylim([-1,3]);




figure();
subplot(231);
    ind = fvxy < 0.1 & tdr;
    binscatter(log10(lys(ind,14,5,5))./log10(lys(ind,14,2,2)),...
               log10(lys(ind,14,3,3))-log10(lys(ind,1,3,3)),[30,30]);
subplot(232);
    ind = 1 > fvxy & fvxy > 0.1 & tdr;
    binscatter(log10(lys(ind,14,5,5))./log10(lys(ind,14,2,2)),...
               log10(lys(ind,14,3,3))-log10(lys(ind,1,3,3)),[30,30]);
subplot(233);
    ind = fvxy > 1 & tdr;
    binscatter(log10(lys(ind,14,5,5))./log10(lys(ind,14,2,2)),...
               log10(lys(ind,14,3,3))-log10(lys(ind,1,3,3)),[30,30]);
subplot(234);
    ind = fvxy < 0.1 & ~tdr;
    binscatter(log10(lys(ind,14,5,5))./log10(lys(ind,14,2,2)),...
              log10(lys(ind,14,3,3))-log10(lys(ind,1,3,3)),[30,30]);
subplot(235);
    ind = 1 > fvxy & fvxy > 0.1 & -tdr;    
    binscatter(log10(lys(ind,14,5,5))./log10(lys(ind,14,2,2)),...
               log10(lys(ind,14,3,3))-log10(lys(ind,1,3,3)),[30,30]);
subplot(236);
    ind =  fvxy > 1 & ~tdr;
    binscatter(log10(lys(ind,14,5,5))./log10(lys(ind,14,2,2)),...
               log10(lys(ind,14,3,3))-log10(lys(ind,1,3,3)),[30,30]);
ForAllSubplots([...
    'xlim([0.9,2]);' ...
    'ylim([-1,3]);' ...
    'Lines([],1,''k'');',...
    'Lines(1.25,[],''k'');',...
    'set(gca,''ColorScale'',''log'');' ...
    'colormap(gca,''jet'');']);


% $$$ f = linspace(1,4,128);
% $$$ eul = zeros([size(lys,1),128]);
% $$$ figure,
% $$$ imagesc(log10(lys([10:20],:,5,5))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ eul(ind,fi) = bweuler(lys([ind-40:ind+40],:,5,5) < f(fi));


tdratio.sup = log10(lys(:,14,3,3))-log10(lys(:,1,3,3));
tdratio.dep = log10(lys(:,14,5,5))-log10(lys(:,1,5,5));


% Inter channel theta phase difference vs theta delta power of lfp_hpcRC (rad?)
chan1 = 1;
chan2 = 3;
cdiff = pi+1;
%cdiff = 0;
tdl = 'sup';
figure();
subplot(231);
    ind = fvxy < 0.1 & tdr;
    binscatter(circ_dist(angle(lys(ind,14,chan1,chan2)),cdiff),...
               tdratio.(tdl)(ind),...
               [30,30]);
subplot(232);
    ind = 1 > fvxy & fvxy > 0.1 & tdr;
    binscatter(circ_dist(angle(lys(ind,14,chan1,chan2)),cdiff),...
               tdratio.(tdl)(ind),...
               [30,30]);
subplot(233);
    ind = fvxy > 1 & tdr;
    binscatter(circ_dist(angle(lys(ind,14,chan1,chan2)),cdiff),...
               tdratio.(tdl)(ind),...
               [30,30]);
subplot(234);
    ind = fvxy < 0.1 & ~tdr;
    binscatter(circ_dist(angle(lys(ind,14,chan1,chan2)),cdiff),...
               tdratio.(tdl)(ind),...
               [30,30]);
subplot(235);
    ind = 1 > fvxy & fvxy > 0.1 & ~tdr;
    binscatter(circ_dist(angle(lys(ind,14,chan1,chan2)),cdiff),...
               tdratio.(tdl)(ind),...
               [30,30]);
subplot(236);
    ind = fvxy > 1 & ~tdr;
    binscatter(circ_dist(angle(lys(ind,14,chan1,chan2)),cdiff),...
               tdratio.(tdl)(ind),...
               [30,30]);
ForAllSubplots([...
    'xlim([-pi,pi]);' ...
    'ylim([-1,3]);' ...
    'Lines([],0.5,''k'');' ...
    'Lines(0,[],''k'');' ...    
    'set(gca,''ColorScale'',''log'');' ...
    'colormap(gca,''jet'');']);

figure();
ind = fvxy < 0.1 & ~tdr;
binscatter(circ_dist(angle(lys(ind,1,chan1,chan2)),cdiff),...
           log10(mean(lys(ind,lfs<18,1,1),2)), ...
               [30,30]);
ForAllSubplots([...
    'xlim([-pi,pi]);' ...
    'ylim([1,4]);' ...
    'set(gca,''ColorScale'',''log'');' ...
    'colormap(gca,''jet'');']);

figure();
ind = fvxy < 0.1 & ~tdr;
binscatter(circ_dist(angle(lys(ind,14,chan1,chan2)),cdiff),...
           circ_dist(angle(lys(ind,14,5,7)),pi/2), ...
               [30,30]);
ForAllSubplots([...
    'xlim([-pi,pi]);' ...
    'ylim([-pi,pi]);' ...
    'set(gca,''ColorScale'',''log'');' ...
    'colormap(gca,''jet'');']);


figure();
ind = tdr;
binscatter(log10(lys(ind,14,3,3))-log10(lys(ind,1,3,3)),...
           fvxy(ind),[20,20],'YLimits',[-3,3],'XLimits',[-1.25,3]);
ForAllSubplots([...
    'xlim([-1.25,3]);' ...
    'ylim([-3,3]);' ...
    'set(gca,''ColorScale'',''log'');' ...
    'colormap(gca,''jet'');']);

chan1 = 1;
chan2 = 3;
figure();
ind = tdr;
binscatter(circ_dist(angle(lys(ind,14,chan1,chan2)),pi),...
           fvxy(ind),'YLimits',[-3,3]);
ForAllSubplots([...
    'xlim([-pi,pi]);' ...
    'ylim([-3,3]);' ...
    'set(gca,''ColorScale'',''log'');' ...
    'colormap(gca,''jet'');']);





chan = 3;
figure();hold('on');
plot(log10(lys(:,14,chan,chan))-log10(lys(:,1,chan,chan)), mean(log10(lys(:,1:18,3,3)),2), '.');
xlim([-1,3]);
ylim([1,4]);

figure();hold('on');
plot(lts,log10(lys(:,12,7,7))-3)
plot(lts,log10(lys(:,12,5,5))-log10(lys(:,2,5,5)),'r')
Lines([],0.5,'k');

figure,
saxes = tight_subplot(14,1,[0.01],[0.1,0.1],[0.1,0.1]);
for k = 2:size(lys,3)
    axes(saxes(k-1));imagesc(lts, lfs, angle(lys(:,:,k-1,k))'); axis('xy'); colormap(gca(),'hsv');caxis([-pi,pi]);grid('on');colorbar();
end;
for k = 1:size(lys,3)
    axes(saxes(k+7));imagesc(lts, lfs, log10(abs(lys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');grid('on'); k=k+1;colorbar();
end; k = k+1;
linkx();

figure,
saxes = tight_subplot(6,1,[0.01],[0.1,0.1],[0.1,0.1]);
for k = 1:size(lys,3)
    axes(saxes(k));imagesc(lts, lfs, log10(abs(lys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');grid('on'); k=k+1;colorbar();caxis([-0.5,1.5]);ylim([20,130])
end; k = k+1;
hold on
plot([1:size(wlfpLEC,1)]./1250, lfpAll(lind,5)/10+80,'k')
plot([1:size(wlfpLEC,1)]./1250, lfpAll(lind,2)/10+80,'b')
linkxy();
ylim([20,130])
caxis([0,1.5]);

ind = 2.3e6:2.8e6;
ind = 1.3e7:1.3696e7;
figure();hold('on');
plot(bsxfun(@minus,lfpHpc(ind,:),[1:30]*200)+8000,'b');
plot(bsxfun(@minus,lfpLec(ind,:),[1:32]*200),'b');
plot(bsxfun(@minus,lfpHpc(ind,:),[1:30]*200)+8000,'b');
plot(bsxfun(@minus,lfpHpcD(ind,:),[1:29]*800)-400,'c');
plot(lfpHRCM(ind).*10+800,'r');
plot(bsxfun(@minus,lfpHpcDD(ind,:)*3,[1:8]*800)+8000,'m');

lfpHpcDDD =  diff(LoadBinary(filepath,181:182,par.nChannels,[],[], [],[])',1,2);


parspec = struct('nFFT',2^11,...
                 'Fs',  par.lfpSampleRate,...
                 'WinLength',2^10,...
                 'nOverlap',2^10*0.875,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[1,30]);
lind = [12000001:12800001];
lind = [ 1800001:2800001 ];
lind = ':';
%wlfpLEC = [lfpHpcDD(lind,:),lfpAll(lind,[2,5]),lfpMCD(lind),lfpHRCM(lind)];
wlfpLEC =  [lfp_hpc(lind,[4,8,11,14,18,22]),lfpHpcRC(lind,[1,3])];
%wlfpLEC =  WhitenSignal(wlfpLEC,[],true);
flagCrossSpec = true;
mode = 'mtcsdglong';
[hys,hfs,hts] = spec(str2func(mode),wlfpLEC,parspec,flagCrossSpec);
hts = hts+(parspec.WinLength/2)/par.lfpSampleRate;


parspec = struct('nFFT',2^12,...
                 'Fs',  par.lfpSampleRate,...
                 'WinLength',2^11,...
                 'nOverlap',2^11*0.875,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[1,30]);
lind = [12000001:12800001];
lind = [ 1800001:2800001 ];
lind = ':';
%wlfpLEC = [lfpHpcDD(lind,:),lfpAll(lind,[2,5]),lfpMCD(lind),lfpHRCM(lind)];
wlfpLEC =  [lfp_hpc(lind,[4,8,11,14,18,22]),lfpHpcRC(lind,[1,3])];
%wlfpLEC =  WhitenSignal(wlfpLEC,[],true);
flagCrossSpec = true;
mode = 'mtcsdglong';
[rys,rfs,rts] = spec(str2func(mode),wlfpLEC,parspec,flagCrossSpec);
rts = rts+(parspec.WinLength/2)/par.lfpSampleRate;



figure,
saxes = tight_subplot(15,1,[0.01],[0.1,0.1],[0.1,0.1]);
for k = 2:size(rys,3)
    axes(saxes(k-1));imagesc(rts, rfs, angle(rys(:,:,k-1,k))'); axis('xy'); colormap(gca(),'hsv');caxis([-pi,pi]);grid('on');colorbar();
end;
for k = 1:size(rys,3)
    axes(saxes(k+7));imagesc(rts, rfs, log10(abs(rys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');grid('on'); k=k+1;colorbar();
end; k = k+1;
linkx();



figure,
saxes = tight_subplot(15,1,[0.01],[0.1,0.1],[0.1,0.1]);
for k = 2:size(hys,3)
    %axes(saxes(k));imagesc(hts, hfs, log10(abs(hys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');grid('on'); k=k+1;colorbar();
    axes(saxes(k-1));imagesc(hts, hfs, angle(hys(:,:,k-1,k))'); axis('xy'); colormap(gca(),'hsv');caxis([-pi,pi]);grid('on');colorbar();
end;
for k = 1:size(hys,3)
    axes(saxes(k+7));imagesc(hts, hfs, log10(abs(hys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');grid('on'); k=k+1;colorbar();
end; k = k+1;
%axes(saxes(k));imagesc(hts, hfs, log10(abs(hys(:,:,9,9)))');axis('xy'); colormap(gca(),'jet');                grid('on'); k=k+1;colorbar();
%axes(saxes(k));imagesc(hts, hfs, (angle(hys(:,:,10,11)))');axis('xy'); colormap(gca(),'hsv');caxis([-pi,pi]);grid('on'); k=k+1;colorbar();
%axes(saxes(k));imagesc(hts, hfs, (angle(hys(:,:,4,5)))');  axis('xy'); colormap(gca(),'hsv');caxis([-pi,pi]);grid('on'); k=k+1;colorbar();
%axes(saxes(k));imagesc(hts,1:10,sq(angle(hys(:,22,end-3,end)))');            colormap(gca(),'hsv');caxis([-pi,pi]);grid('on');       colorbar();
linkx();


hfvxy = ButFilter(sq(xyz.XYZ(lind,1,[1,2])),4, 1.2./(sampleRate/2), 'low');
%hfvxy = sqrt((sum(((circshift(hfvxy,-1)-circshift(hfvxy,1)).* (sampleRate/2)).^2,2)));
hfvxy = log10(sqrt((sum(((circshift(hfvxy,-1)-circshift(hfvxy,1)).* (sampleRate/2)).^2,2))));
hfvxy = clip(interp1([1:size(hfvxy,1)]./sampleRate, hfvxy, hts),-3,3);;


% JPDF between the local phase differences within the radiatum and the lacosum moleculare.
figure,
binscatter(angle(hys(:,12,2,3)),...
           angle(hys(:,12,5,6)),...
           30,...
           'YLimits',[-pi,pi],...
           'XLimits',[-pi,pi]);
%daspect([0.5,1,1]);
% $$$ xlim(gca, [-1,0.5]);
% $$$ ylim(gca, [-0.5,2]);
set(gca,'ColorScale','log');
colormap(gca,'jet');



% JPDF between the local phase differences within the radiatum and the lacosum moleculare.
tdr = log10(hys(:,12,5,5))-log10(hys(:,1,5,5));
figure,
subplot(121);
ind = tdr<1;
binscatter(angle(hys(ind,12,2,3)),...
           angle(hys(ind,12,5,6)),...
           30,...
           'YLimits',[-pi,pi],...
           'XLimits',[-pi,pi]);
subplot(122);
ind = tdr>1;
binscatter(angle(hys(ind,12,2,3)),...
           angle(hys(ind,12,5,6)),...
           30,...
           'YLimits',[-pi,pi],...
           'XLimits',[-pi,pi]);
ForAllSubplots([...
    'xlim(gca, [-1,0.5]);'...
    'ylim(gca, [-0.5,2]);'...
    'set(gca,''ColorScale'',''log'');'...
    'colormap(gca,''jet'');']);
ForAllSubplots([...
    'xlim(gca, [-pi,pi]);'...
    'ylim(gca, [-pi,pi]);'...
    'set(gca,''ColorScale'',''log'');'...
    'colormap(gca,''jet'');']);



% JPDF between the local phase differences within the radiatum and the lacosum moleculare.
chanRC = 8;
figure,
binscatter(angle(hys(:,12,3,5)),...
           log10(hys(:,12,chanRC,chanRC))-log10(hys(:,1,chanRC,chanRC)),...
           30,...
           'YLimits',[-1,3],...
           'XLimits',[-pi,pi]);
%daspect([0.5,1,1]);
xlim(gca, [-pi,pi]);
ylim(gca, [-1,3]);
set(gca,'ColorScale','log');
colormap(gca,'jet');


% JPDF between the local phase differences within the radiatum and the lacosum moleculare.
figure,
subplot(121);
    ind = hfvxy<0;
    binscatter(angle(hys(ind,12,1,2)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-1,0.5],...
               'YLimits',[-1,3]);
subplot(122);
    ind = hfvxy>0;
    binscatter(angle(hys(ind,12,1,2)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-1,0.5],...
               'YLimits',[-1,3]);
ForAllSubplots([...
    'xlim(gca, [-1,0.5]);'...
    'ylim(gca, [-1,3]);'...
    'grid(gca,''on'');'...
    'set(gca,''ColorScale'',''log'');'...
    'colormap(gca,''jet'');']);


% JPDF between the local phase differences within the radiatum and the lacosum moleculare.
chanRC = 8;
figure,
subplot(141);
    ind = hfvxy<-0.5;
    binscatter(angle(hys(ind,12,4,6)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-pi,pi],...
               'YLimits',[-1,3]);
subplot(142);
    ind = hfvxy>-0.5 & hfvxy<0;
    binscatter(angle(hys(ind,12,4,6)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-pi,pi],...
               'YLimits',[-1,3]);
subplot(143);
    ind = hfvxy>0 & hfvxy<0.5;
    binscatter(angle(hys(ind,12,4,6)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-pi,pi],...
               'YLimits',[-1,3]);
subplot(144);
    ind = hfvxy>0.5;
    binscatter(angle(hys(ind,12,4,6)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-pi,pi],...
               'YLimits',[-1,3]);
    ForAllSubplots([...
    'xlim(gca, [-pi,pi]);'...
    'ylim(gca, [-1,3]);'...
    'grid(gca,''on'');'...    
    'set(gca,''ColorScale'',''log'');'...
    'colormap(gca,''jet'');']);


    
% JPDF between the local phase differences within the radiatum and the lacosum moleculare.
chanRC = 8;
figure,
subplot(141);
    ind = hfvxy<-0.5;
    binscatter(angle(hys(ind,12,3,5)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-pi,pi],...
               'YLimits',[-1,3]);
subplot(142);
    ind = hfvxy>-0.5 & hfvxy<0;
    binscatter(angle(hys(ind,12,3,5)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-pi,pi],...
               'YLimits',[-1,3]);
subplot(143);
    ind = hfvxy>0 & hfvxy<0.5;
    binscatter(angle(hys(ind,12,3,5)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-pi,pi],...
               'YLimits',[-1,3]);
subplot(144);
    ind = hfvxy>0.5;
    binscatter(angle(hys(ind,12,3,5)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-pi,pi],...
               'YLimits',[-1,3]);
    ForAllSubplots([...
    'xlim(gca, [-pi,pi]);'...
    'ylim(gca, [-1,3]);'...
    'grid(gca,''on'');'...    
    'set(gca,''ColorScale'',''log'');'...
    'colormap(gca,''jet'');']);


% JPDF between the local phase differences within the radiatum and the lacosum moleculare.
chanRC = 8;
figure,
subplot(141);
    ind = hfvxy<-0.5;
    binscatter(angle(hys(ind,12,5,6)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-1,0.5],...%               'XLimits',[-pi,pi],...
               'YLimits',[-1,3]);
subplot(142);
    ind = hfvxy>-0.5 & hfvxy<0;
    binscatter(angle(hys(ind,12,5,6)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-1,0.5],...%               'XLimits',[-pi,pi],...
               'YLimits',[-1,3]);
subplot(143);
    ind = hfvxy>0 & hfvxy<0.5;
    binscatter(angle(hys(ind,12,5,6)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-1,0.5],...%               'XLimits',[-pi,pi],...
               'YLimits',[-1,3]);
subplot(144);
    ind = hfvxy>0.5;
    binscatter(angle(hys(ind,12,5,6)),...
               log10(hys(ind,12,chanRC,chanRC))-log10(hys(ind,1,chanRC,chanRC)),...
               30,...
               'XLimits',[-1,0.5],...
               'YLimits',[-1,3]);
ForAllSubplots([...
    'xlim(gca, [-1,0.5]);'...
    'ylim(gca, [-1,3]);'...
    'grid(gca,''on'');'...    
    'set(gca,''ColorScale'',''log'');'...
    'colormap(gca,''jet'');']);


% JPDF between the local phase differences within the radiatum and the lacosum moleculare.
chanRC = 8;
figure,
subplot(121);
    ind = hfvxy<0;
    binscatter(log10(hys(ind,12,3,3))-log10(hys(ind,1,3,3)),...
               hys(ind,12,5,5)./hys(ind,12,6,6),...
               30,...
               'XLimits',[-pi,6],...
               'YLimits',[-1,3]);
subplot(122);
    ind = hfvxy>0;
    binscatter(log10(hys(ind,12,3,3))-log10(hys(ind,1,3,3)),...
               hys(ind,12,5,5)./hys(ind,12,6,6),...
               30,...
               'XLimits',[-pi,6],...
               'YLimits',[-1,3]);
ForAllSubplots([...
    'xlim(gca, [-pi,pi]);'...
    'ylim(gca, [-1,3]);'...
    'grid(gca,''on'');'...    
    'set(gca,''ColorScale'',''log'');'...
    'colormap(gca,''jet'');']);


% JPDF between the local phase differences within the radiatum and the lacosum moleculare.
figure,
subplot(121);
    ind = hfvxy<0;
    binscatter(log10(hys(ind,12,5,5))-log10(hys(ind,1,5,5)),...
               log10(hys(ind,12,6,6))-log10(hys(ind,1,6,6)),...
               30,...
               'XLimits',[-1,3],...
               'YLimits',[-1,3]);
subplot(122);
    ind = hfvxy>0;
    binscatter(log10(hys(ind,12,5,5))-log10(hys(ind,1,5,5)),...
               log10(hys(ind,12,6,6))-log10(hys(ind,1,6,6)),...
               30,...
               'XLimits',[-1,3],...
               'YLimits',[-1,3]);
ForAllSubplots([...
    'xlim(gca, [-1,3]);'...
    'ylim(gca, [-1,3]);'...
    'set(gca,''ColorScale'',''log'');'...
    'colormap(gca,''jet'');']);


chan1 = 2;
chan2 = 6;
figure,
ind = ':';
plot((log10(hys(ind,12,chan1,chan1))-log10(hys(ind,1,chan1,chan1)))+(log10(hys(ind,12,chan2,chan2))-log10(hys(ind,1,chan2,chan2))));
Lines([],0,'k');
hold on
plot((log10(hys(ind,12,chan1,chan1))-log10(hys(ind,1,chan1,chan1)))+4);
plot((log10(hys(ind,12,chan2,chan2))-log10(hys(ind,1,chan2,chan2)))+4);
Lines([],4,'k');
plot(hfvxy-5)
Lines([],-5,'k');

% JPDF between the local phase differences within the radiatum and the lacosum moleculare.
figure,
binscatter(angle(hys(:,12,1,2)),...
           angle(hys(:,12,2,5)),...
           30,...
           'YLimits',[-pi,pi],...
           'XLimits',[-pi,pi]);
daspect([0.5,1,1]);
xlim(gca, [-1,0.5]);
ylim(gca, [-pi,pi]);
set(gca,'ColorScale','log');
colormap(gca,'jet');


% JPDF between the local phase differences within the radiatum and the lacosum moleculare.
figure,
subplot(121);
binscatter(angle(hys(hfvxy<0.2,12,1,2)),...
           angle(hys(hfvxy<0.2,12,2,6)),...
           30,...
           'YLimits',[-pi,pi],...
           'XLimits',[-pi,pi]);
subplot(122);
binscatter(angle(hys(hfvxy>0.2,12,1,2)),...
           angle(hys(hfvxy>0.2,12,2,6)),...
           30,...
           'YLimits',[-pi,pi],...
           'XLimits',[-pi,pi]);
%daspect([0.5,1,1]);
ForAllSubplots([...
    'xlim(gca, [-1,0.5]);'...
    'ylim(gca, [-pi,pi]);'...
    'set(gca,''ColorScale'',''log'');'...
    'colormap(gca,''jet'');']);



rad_dphz_sub = discretize( angle(hys(:,14,1,2)), linspace([-1,0.2,30]));
lcm_dphz_sub = discretize( angle(hys(:,14,3,4)), linspace([-1,2,30]));
hfvxy_val = hfvxy(:);
ind = nniz([rad_dphz_sub, lcm_dphz_sub, hfvxy_val]);
A = accumarray( [rad_dphz_sub(ind), lcm_dphz_sub(ind)], hfvxy_val(ind), [29,29], @mean);
S = accumarray( [rad_dphz_sub(ind), lcm_dphz_sub(ind)], hfvxy_val(ind), [29,29], @std);

rad_dphz_sub = discretize( angle(hys(:,14,1,2)), linspace([-1,0.5,30]));
lcm_dphz_sub = discretize( angle(hys(:,14,2,6)), linspace([-pi,pi,30]));
hfvxy_val = hfvxy(:);
ind = nniz([rad_dphz_sub, lcm_dphz_sub, hfvxy_val]);
A = accumarray( [rad_dphz_sub(ind), lcm_dphz_sub(ind)], hfvxy_val(ind), [29,29], @mean);
S = accumarray( [rad_dphz_sub(ind), lcm_dphz_sub(ind)], hfvxy_val(ind), [29,29], @std);

figure
subplot(121);
    imagesc(A');
subplot(122);
    imagesc(S');
ForAllSubplots([...
'axis(gca(),''xy'');'...
'colormap(gca(),''jet'');'...
'colorbar(gca());']);




lfpAll =  LoadBinary(filepath,channelsAll,par.nChannels,[],[], [],[])';

lfpLC =  LoadBinary(filepath,[116,68],par.nChannels,[],[], [],[])';
lfpLCD =  diff(LoadBinary(filepath,[116,68],par.nChannels,[],[], [],[])',1,2);
lfpMCD =  diff(LoadBinary(filepath,[310,312],par.nChannels,[],[], [],[])',1,2);
lfpMCS =  diff(LoadBinary(filepath,[320,318],par.nChannels,[],[], [],[])',1,2);

lfpLRCS = diff(LoadBinary(filepath,[128,120],par.nChannels,[],[], [],[])',1,2);
lfpLRCS = diff(LoadBinary(filepath,[128,114],par.nChannels,[],[], [],[])',1,2);
lfpMRC =  LoadBinary(filepath,[299,311],par.nChannels,[],[], [],[])';
lfpMRCD = diff(lfpMRC,1,2);

lfpHRCD = diff(LoadBinary(filepath,[158,167],par.nChannels,[],[], [],[])',1,2);
lfpHRCS = diff(LoadBinary(filepath,[140,153],par.nChannels,[],[], [],[])',1,2);
%lfpHRCM = diff(LoadBinary(filepath,[140,143],par.nChannels,[],[], [],[])',1,2);
lfpHRCM = diff(LoadBinary(filepath,[150,156],par.nChannels,[],[], [],[])',1,2);

lfpRC =  LoadBinary(filepath,[159,167],par.nChannels,[],[], [],[])';
lfpRCD = diff(lfpRC,1,2);


lfpVC =  LoadBinary(filepath,[129 ,178],par.nChannels,[],[], [],[])';




figure();
hold('on');
plot(lfpTS,lfpAll(:,5)/3+1000,'b');
plot(lfpTS,lfpLCD-mean(lfpLCD)+lfpMCD-mean(lfpMCD)+500,'r');
plot(lfpTS,lfpLCD-mean(lfpLCD),'m');
plot(lfpTS,lfpMCD-mean(lfpMCD)+250,'g');

figure();
hold('on');
plot(lfpTS,lfpAll(:,8));
plot(lfpTS,lfpAll(:,10));
plot(lfpTS,lfpAll(:,5)+1000);
plot(lfpTS,lfpLRCS(:,1)*3-1000)
plot(lfpTS,lfpHRCM(:,1)-2000)
plot(lfpTS,(lfpLCD(:,1)*8)+4000,'m')
plot(lfpTS,(lfpMCD(:,1)*10)+1000,'g') 
Lines([],-1000,'k');


parspec = struct('nFFT',2^12,...
                 'Fs',  par.lfpSampleRate,...
                 'WinLength',2^10,...
                 'nOverlap',2^10*0.5,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[1,20]);
%lind = [2125001:3375001];
lind = ':';
wlfpLEC =  [lfpAll(lind, 5), lfpHRCM(lind), ...
            lfpAll(lind,10), lfpLCD(lind), ...
            lfpAll(lind,15), lfpMCD(lind)];
%wlfpLEC =  WhitenSignal([],[],false);
flagCrossSpec = true;
mode = 'mtcsdglong';
[vys,fs,ts] = spec(str2func(mode),wlfpLEC,parspec,flagCrossSpec);
ts = ts+(parspec.WinLength/2)/par.lfpSampleRate;


figure,plot(log10(mean(abs(vys(:,fs>6&fs<12,3,3)),2)),...
            log10(mean(abs(vys(:,fs>6&fs<12,7,7)),2)),'.');



tdratio = log10(mean(vys(:,18:28,1,1),2))-log10(mean(vys(:,2:6,1,1),2));
tdratio(~nniz(tdratio)) = 0;
ftdratio = ButFilter(tdratio,4,1/(9.765*0.5),'low');
figure,plot(ts,ftdratio)

chanlabels = {'CA1pyr','CA1lm','CA1pyrRC','CA1dgRC','MEC3RC'};
figure,hold on,
for chan1 = 1:size(hys,3),
    for chan2 = chan1+1:size(hys,3)
        subplot2(6,6,chan1,chan2);
coh = [];
k = 1;
for thresh = linspace(-1,1.5,15),
ind = ftdratio > thresh;
coh(k,:) = mean(abs(hys(ind,:,chan1,chan2)))./(mean(sqrt(hys(ind,:,chan1,chan1).*hys(ind,:,chan2,chan2))));
k=k+1;
end
imagesc(linspace(-1,1.5,15),fs,coh');
axis('tight');
axis('xy');
colorbar();
colormap('jet');
caxis([0.35,1]);
title([chanlabels{chan1} ' X ' chanlabels{chan2}]);
end
end

figure,
saxes = tight_subplot(11,1,[0.0],[0.1,0.1],[0.1,0.1]);
k = 1;
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,1,1)))'); axis('xy'); colormap(gca(),'jet');caxis([2,3.75]);grid('on'); k=k+1;colorbar();
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,2,2)))'); axis('xy'); colormap(gca(),'jet');caxis([2,3.25]);grid('on'); k=k+1;colorbar();
axes(saxes(k));imagesc(ts, fs, (angle(vys(:,:,1,2)))'); axis('xy'); colormap(gca(),'hsv');caxis([-pi/2,pi/2]);grid('on'); k=k+1;colorbar();
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,3,3)))'); axis('xy'); colormap(gca(),'jet');caxis([1.25,2.5]);grid('on'); k=k+1;colorbar();
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,4,4)))'); axis('xy'); colormap(gca(),'jet');caxis([0.25,1.75]);grid('on'); k=k+1;colorbar();
axes(saxes(k));imagesc(ts, fs, (angle(vys(:,:,3,4)))'); axis('xy'); colormap(gca(),'hsv');caxis([-pi,0]);grid('on'); k=k+1;colorbar();
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,5,5)))'); axis('xy'); colormap(gca(),'jet');caxis([1,3]);grid('on'); k=k+1;colorbar();
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,6,6)))'); axis('xy'); colormap(gca(),'jet');caxis([-0.5,2]);grid('on'); k=k+1;colorbar();
axes(saxes(k));imagesc(ts, fs, (angle(vys(:,:,5,6)))'); axis('xy'); colormap(gca(),'hsv');caxis([-pi/2,pi/2]);grid('on'); k=k+1;colorbar();
axes(saxes(k));imagesc(ts, fs, (angle(vys(:,:,2,6)))'); axis('xy'); colormap(gca(),'hsv');caxis([-pi/2,pi/2]);grid('on'); k=k+1;colorbar();
linkxy();
axes(saxes(k));plot(lfpTS,xyz.XYZ(:,:,3));grid('on'); ylim(saxes(k),[-50,300]);k=k+1;colorbar(); 
linkx();

figure,
binscatter(ftdratio,(angle(vys(:,22,1,2))),'YLimits',[-3,3],'XLimits',[-0.6,2]);set(gca,'ColorScale','log');colormap('jet')

figure,
binscatter(ftdratio,(angle(vys(:,22,2,6))),'YLimits',[-3,3],'XLimits',[-0.6,2]);set(gca,'ColorScale','log');colormap('jet')

figure,
binscatter(ftdratio,log10(abs(vys(:,22,2,2))),...
           'YLimits',[1.5,4.5],...
           'XLimits',[-0.6,2],...
           'ShowEmptyBins','OFF');
set(gca,'ColorScale','log');colormap('jet')


figure,
binscatter(ftdratio,log10(abs(vys(:,22,1,1))),...
           'YLimits',[1.5,4.5],...
           'XLimits',[-0.6,2],...
           'ShowEmptyBins','OFF');
set(gca,'ColorScale','log');colormap('jet')


figure,
binscatter(log10(abs(vys(:,3,1,1))),log10(abs(vys(:,22,1,1))),...
           'YLimits',[1.5,4.5],...
           'XLimits',[1.5,4.5],...
           'ShowEmptyBins','OFF');
set(gca,'ColorScale','log');colormap('jet')

fvxy = ButFilter(sq(xyz.XYZ(:,1,[1,2])),4, 2.4./(1250/2), 'low');
fvxy = log10(sqrt((sum(((circshift(fvxy,-1)-circshift(fvxy,1)).* (1250/2)).^2,2))));
fvxy = interp1(lfpTS,fvxy,lts);

figure,
ind = ftdratio>0.5 & fvxy<0.25;
%binscatter(log10(abs(vys(ind,22,1,1))),log10(abs(vys(ind,22,6,6))))
binscatter(log10(abs(vys(ind,22,6,6))),angle(vys(ind,22,1,6)))
set(gca,'ColorScale','log');colormap('jet')
xlim([0,2.5]);
%xlim([2.5,4]);
ylim([-pi,pi]);
%ylim([0,2.5]);


figure,
ind = ftdratio>0.5;
%binscatter(log10(abs(vys(ind,22,1,1))),log10(abs(vys(ind,22,6,6))))
binscatter(fvxy(ind),angle(vys(ind,22,1,6)))
set(gca,'ColorScale','log');colormap('jet')
xlim([-3,3]);
ylim([-pi,pi]);
%ylim([0,2.5]);


figure,
subplot(211);
ind = ftdratio>0.5 & fvxy<0.25;
binscatter(angle(vys(ind,22,1,2)),log10(abs(vys(ind,22,1,1)))-log10(abs(vys(ind,22,2,2))))
set(gca,'ColorScale','log');colormap('jet')
xlim([-pi,pi]);
ylim([-0.5,1.5]);
Lines([],0.75,'k');
Lines(0.5,[],'k');
subplot(212);
ind = ftdratio>0.5 & fvxy>0.25;
binscatter(angle(vys(ind,22,1,2)),log10(abs(vys(ind,22,1,1)))-log10(abs(vys(ind,22,2,2))))
set(gca,'ColorScale','log');colormap('jet')
xlim([-pi,pi]);
ylim([-0.5,1.5]);
Lines([],0.75,'k');
Lines(0.5,[],'k');


tdratio = log10(mean(vys(:,9:14,1,1),2))-log10(mean(vys(:,1:3,1,1),2));
tdratio(~nniz(tdratio)) = 0;
ftdratio = ButFilter(tdratio,4,1/(9.765*0.5),'low');

pyr lm hrc lmrc mrc

figure,
subplot(211); hold('on');
ind = ftdratio>0.5 & fvxy<0.5;
binscatter(angle(hys(ind,11,1,3)),log10(abs(hys(ind,11,5,5))))
set(gca,'ColorScale','log');colormap('jet')
xlim([-pi,pi]); ylim([0,4]); Lines([],2,'k'); Lines(0,[],'k');
grid('on');
subplot(212); hold('on');
ind = ftdratio>0.5 & fvxy>0.5;
binscatter(angle(hys(ind,11,1,3)),log10(abs(hys(ind,11,5,5))))
set(gca,'ColorScale','log');colormap('jet')
xlim([-pi,pi]); ylim([0,4]); Lines([],2,'k'); Lines(0,[],'k');
grid('on');


figure,
hold('on');
ind = ftdratio>0.5 & fvxy>0;
plot(circ_dist(mean(angle(hys(ind,10:15,1,3)),2),pi/2),log10(abs(hys(ind,11,5,5))),'.b');
ind = ftdratio>0.5 & fvxy<0;
plot(circ_dist(mean(angle(hys(ind,10:15,1,3)),2),pi/2),log10(abs(hys(ind,11,5,5))),'.r');
grid('on');

figure,
hold('on');
ind = ftdratio>0.5 & fvxy>0;
plot(circ_dist(circ_mean(angle(hys(ind,10:15,1,3)),[],2),pi/2),log10(mean(abs(hys(ind,:,1,1)),2)),'.b');
ind = ftdratio>0.5 & fvxy<0;
plot(circ_dist(circ_mean(angle(hys(ind,10:15,1,3)),[],2),pi/2),log10(mean(abs(hys(ind,:,1,1)),2)),'.r');
grid('on');


figure,
hold('on');
subplot(121)
ind = ftdratio>0.5 & fvxy>0;
hist2([circ_dist(circ_mean(angle(hys(ind,10:15,1,3)),[],2),pi),circ_dist(circ_mean(angle(hys(ind,10:15,1,5)),[],2),-pi)],linspace([-pi,pi,30]),linspace(-pi,pi,30))
Lines(0.5,[],'k');
Lines([],1,'k');
subplot(122);
ind = ftdratio>0.5 & fvxy<0;
hist2([circ_dist(circ_mean(angle(hys(ind,10:15,1,3)),[],2),pi),circ_dist(circ_mean(angle(hys(ind,10:15,1,5)),[],2),-pi)],linspace([-pi,pi,30]),linspace(-pi,pi,30))
Lines(0.5,[],'k');
Lines([],1,'k');



figure,
hold('on');
subplot(131)
ind = ftdratio>0.6 & fvxy>1;
hist2([log10(mean(abs(hys(ind,10:15,5,5)),2)), log10(mean(abs(hys(ind,10:15,1,1)),2))],linspace([0,3,30]),linspace(2,3.5,30))
Lines(1.75,[],'k');
Lines([],2.75,'k');
subplot(132)
ind = ftdratio>0.6 & fvxy>0&fvxy<1;
hist2([log10(mean(abs(hys(ind,10:15,5,5)),2)), log10(mean(abs(hys(ind,10:15,1,1)),2))],linspace([0,3,30]),linspace(2,3.5,30))
Lines(1.75,[],'k');
Lines([],2.75,'k');
subplot(133);
ind = ftdratio>0.6 & fvxy<0;
hist2([log10(mean(abs(hys(ind,10:15,5,5)),2)), log10(mean(abs(hys(ind,10:15,1,1)),2))],linspace([0,3,30]),linspace(2,3.5,30))
Lines(1.75,[],'k');
Lines([],2.75,'k');


figure,
saxes = tight_subplot(4,1,[0.0],[0.1,0.1],[0.1,0.1]);
k = 1;
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');caxis([1.5,3.5]);grid('on'); k=k+1;
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');caxis([2,3.75]);grid('on'); k=k+1;
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');caxis([0,2.0]);grid('on'); k=k+1
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');caxis([-0.5,1.5]);grid('on');k=k+1;
linkxy();


figure,
saxes = tight_subplot(7,1,[0.0],[0.1,0.1],[0.1,0.1]);
k = 1;
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');caxis([0,3.5]);grid('on'); k=k+1;
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');caxis([0,3.5]);grid('on'); k=k+1;
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');caxis([0,3.5]);grid('on'); k=k+1
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');caxis([0,3.5]);grid('on');k=k+1;
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');caxis([0,3.5]);grid('on');k=k+1;
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');caxis([0,3.5]);grid('on');k=k+1;
axes(saxes(k));imagesc(ts, fs, log10(abs(vys(:,:,k,k)))'); axis('xy'); colormap(gca(),'jet');caxis([0,3.5]);grid('on');
linkxy();
ForAllSubplots('ylim([1,35]);')


figure();
hold('on');
plot(lfpTS,lfpAll(:,8));
plot(lfpTS,lfpAll(:,10));
plot(lfpTS,lfpAll(:,5)+1000);
plot(lfpTS,lfpLRCS(:,1)-1000)
Lines([],-1000,'k');


parspec = struct('nFFT',2^11,...
                 'Fs',  par.lfpSampleRate,...
                 'WinLength',2^10,...
                 'nOverlap',2^10*.875,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[1,100]);
wlfpLEC =  lfpVC;
%wlfpLEC =  WhitenSignal([lfpLEC(:,[1,2,5]),dlfpLEC(:,1)],[],false);
flagCrossSpec = true;
mode = 'mtcsdglong';
[vys,fs,ts] = spec(str2func(mode),wlfpLEC,parspec,flagCrossSpec);
ts = ts+(parspec.WinLength/2)/par.lfpSampleRate;


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

lfpAll = [lfpAll,lfpMRCD(:,1),lfpLRCS(:,1)];
lfpAll = [lfpAll,lfpRCD(:,1)];

ind = (log10(vys(:,11,2,2))-log10(vys(:,3,2,2)))>0.5;
frqRng = 5<fs&fs<12;
pmat = []
cmat = [];
for c1 = 1:size(lfpAll,2)-1
    for c2 = (c1+1):size(lfpAll,2)
        tic    
        parspec = struct(     'nFFT' , 2^11,                    ...
                                'Fs' , par.lfpSampleRate,       ...
                         'WinLength' , 2^10,                    ...
                          'nOverlap' , 2^10 * 0.5,              ...
                                'NW' , 3,                       ...
                           'Detrend' , [],                      ...
                           'nTapers' , [],                      ...
                         'FreqRange' , [1, 20]                 ...
        );
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
    

figure();
subplot(122);
imagesc([1:size(cmat,1)],[1:size(cmat,2)],cmat')
colormap('jet');
caxis([0,1]);
set(gca,'ytick',[1:size(cmat,2)]);
set(gca,'yticklabel',{'hpc:ca1-orid','hpc:ca1-oris','hpc:ca1-pyr','hpc:ca1-rad','hpc:ca1-lmd','hpc:ca1-lms','hpc:ca1-dg',...
        'lec-l1','lec-l2','lec-l3','lec-l4','lec-l5','lec-l6',...
        'mec-l1','mec-l2','mec-l3','mec-l4','mec-l5','mec-l6',...
        'mec-rc-1','lec-rc-2','hpc-rc-1'});
set(gca,'xtick',[1:size(cmat,1)]);
set(gca,'xticklabel',{'hpc:ca1-orid','hpc:ca1-oris','hpc:ca1-pyr','hpc:ca1-rad','hpc:ca1-lmd','hpc:ca1-lms','hpc:ca1-dg',...
        'lec-l1','lec-l2','lec-l3','lec-l4','lec-l5','lec-l6',...
        'mec-l1','mec-l2','mec-l3','mec-l4','mec-l5','mec-l6',...
        'mec-rc-1','lec-rc-2','hpc-rc-1'});
set(gca,'XTickLabelRotation',90);
cax = colorbar();
ylabel(cax,'Coherence')
title('Inter-Region/Layer Coherence for the Theta [6-10Hz] Frequency Band')
subplot(121);
imagesc(fs,1:size(pmat,1),log10(bsxfun(@rdivide,10.^pmat,sum(10.^pmat))));
colormap('jet');
set(gca,'YTick',[1:size(pmat,1)]);
set(gca,'YTickLabel',...
        {'hpc:ca1-orid','hpc:ca1-oris','hpc:ca1-pyr','hpc:ca1-rad','hpc:ca1-lmd','hpc:ca1-lms','hpc:ca1-dg',...
        'lec-l1','lec-l2','lec-l3','lec-l4','lec-l5','lec-l6',...
        'mec-l1','mec-l2','mec-l3','mec-l4','mec-l5','mec-l6',...
         'mrc-1','lrc-2','hrc-1'});
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