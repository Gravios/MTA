%% Figure 14 - Spectral power for each behavior
%sname = 'jg05-20120317';
%sname = 'jg05-20120310';
%sname = 'jg05-20120309';
sname = 'er02-20110908.cof.all';
display = 0;
%chans = [1:2:8];
chans = [65:1:96];
marker = 'spine_lower';

Trial = MTATrial.validate(sname);

xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);


Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',chans);

%lfp.data = cat(2,lfp.data,fet_rhm(Trial,lfp,'raw'));

%stc = Trial.load('stc','nn0317_PP');
%stc = Trial.load('stc','nn0317');    

stc = Trial.load('stc','msnn_ppsvd_raux');
states = {'loc','rear','pause','groom','sit'};

%wlfp = WhitenSignal(lfp.data,[],1);    
wlfp = WhitenSignal(lfp.data,round(180*lfp.sampleRate),true);    
%wlfp = lfp.data;


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

%% Calculate clipped psd 1-40 Hz
dnanThreshl = [prctile(reshape(10.*log10(yld.data),[],size(yld,3)),75)'];
nanyld = yld.copy;
nanyld.data = 10*log10(nanyld.data)>repmat(permute(dnanThreshl(:,1),[3,2,1]),yld.size([1,2]));


%% compute high frequency stuff    
th=[];    fh=[];    yh=[];
spectral.nfft = 2^8;
spectral.window = 2^7;
spectral.freq = [30,130];
parfor i = 1:lfp.size(2),
    [yh(:,:,i),fh(:,i),th(:,i)] = mtchglong(wlfp(:,i),spectral.nfft,lfp.sampleRate,spectral.window,spectral.window*0.875,[],[],[],spectral.freq);
end
fh = fh(:,1);
th = th(:,1);
tindShift = round(spectral.window/2/lfp.sampleRate*1/diff(th(1:2,1))-1);        
yhd = MTADlfp('data',yh,'sampleRate',1/diff(th(1:2, 1)));
yhd.data(yhd.data==0)=nan;
tyhd = yhd.copy;    
yhd.data = clip(yhd.data,0.01,1e8);
nyhd = yhd.copy;    
nyhd.data = bsxfun(@rdivide,...
                   bsxfun(@minus,yhd.data,nanmedian(yhd.data)),...
                   nanstd(yhd.data));

th = [[0:tindShift-1]'./yhd.sampleRate;th+spectral.window/2/lfp.sampleRate];
try,delete(gcp('nocreate')),end

dnanThresh = [prctile(reshape(10.*log10(yhd.data(:,fh(:,1)<130,:,:)),[],size(yld,3)),75)'];
nanyhd = yhd.copy;
nanyhd.data = 10*log10(nanyhd.data)>repmat(permute(dnanThresh(:,1),[3,2,1]),yhd.size([1,2]));

fnyhd = nyhd.copy();
fnyhd.filter('RectFilter',3,3);
fnyhd.data = permute(fnyhd.data,[2,3,1]);
fnyhd.filter('RectFilter',3,3);
fnyhd.data = permute(fnyhd.data,[2,3,1]);
fnyhd.filter('RectFilter',3,3);
fnyhd.data = permute(fnyhd.data,[2,3,1]);

fdnanThresh = [prctile(reshape(10.*log10(fnyhd.data(:,fh(:,1)<130,:,:)),[],size(yld,3)),75)'];
fnanyhd = nyhd.copy();
fnanyhd.data = 10*log10(fnanyhd.data)>repmat(permute(fdnanThresh(:,1),[3,2,1]),yhd.size([1,2]));


sp = [];
figure,
sp(1) = subplot(211);  imagesc(th,fh,sq(fnyhd(:,:,5))'),  caxis([-1,3]),  colormap jet,  axis xy
sp(2) = subplot(212);  imagesc(th,fh,sq(fnyhd(:,:,12))'),caxis([ -1,3]),  colormap jet,  axis xy
linkaxes(sp,'xy');

[gmins,gvals] = LocalMinimaN(-fnyhd.data(:,:,:),0,[5,5,1]);




ind = 1:1000;
sp = [];
offset = 5
figure,
for c = 1:10
    sp(end+1) = subplot(10,1,c);  
    imagesc(th(ind),fh,sq(fnyhd(ind,:,c+offset))'),  caxis([-1,3]),  axis('xy');
    colormap('jet');  
    hold('on'); 
    plot(th(gmins(gmins(:,3)==c+offset,1)),fh(gmins(gmins(:,3)==c+offset,2)),'*m')
end
linkaxes(sp,'xy');



fxyz = xyz.copy;
fxyz.filter('RectFilter',3,3);

%% compute head speed    
vxy = xyz.copy;
vxy.filter('ButFilter',3,2.4);
vxy = vxy.vel([1,7],[1,2]);

vxyl = vxy.copy;    
vxyl.resample(yld);
vxyl.data(vxyl.data<1e-3)=1e-3;
vxyl.data = log10(vxyl.data);

vxyh = vxy.copy;        
vxyh.resample(yhd);
vxyh.data(vxyh.data<1e-3)=1e-3;
vxyh.data = log10(vxyh.data);

vxy.data(vxy.data<1e-3)=1e-3;
vxy.data = log10(vxy.data);




%% Computing theta phase for all channels 
tlfp = lfp.copy;
tlfp.data = lfp(:,1:numel(chans));
tlfp.resample(yhd);
thetaPhase =  tlfp.phase;       
tper = Trial.stc{'t'};
tper.cast('TimeSeries');
tper.resample(yhd);







sp = [];
figure,
sp(1) = subplot(211);  imagesc(th,fh,sq(nyhd(:,:,5))'),  caxis([-1,3]),  colormap jet,  axis xy
sp(2) = subplot(212);  imagesc(th,fh,sq(nanyhd(:,:,5))'),caxis([ 0,1]),  colormap jet,  axis xy
linkaxes(sp,'xy');



%% plot steady state spec count
hfig = figure(93993494);
hfig.Units = 'centimeters';
hfig.Position = [3,3,45,6];
hfig.PaperPositionMode = 'auto';
vstep = -2:0.5:1.5;
for i = 1:numel(vstep)    
    ind = stc{'r'}+[0.2,-0.2];
    ind.cast('TimeSeries',yld);
    ind = ind.data==1&vxyl(:,2)>0+vstep(i)&vxyl(:,2)<0.5+vstep(i);    
    subplot(1,numel(vstep),i);  
    imagesc(fl,1:32,sq(nansum(nanyld(ind,:,1:numel(chans))))'./sum(ind));
    caxis([0,1]);
    %colormap jet    
    colormap parula    
    xlabel('Frequency Hz')
end
suptitle(['steady state spec count: ',Trial.filebase,' (1-40Hz)'])



%% plot steady state spec count
hfig = figure(9399349);
hfig.Units = 'centimeters';
hfig.Position = [3,3,45,6];
hfig.PaperPositionMode = 'auto';
vstep = -2:0.5:1.5;
for i = 1:numel(vstep)    
    ind = stc{'r'}+[0.1,-0.1];
    ind.cast('TimeSeries',yhd);
    ind = ind.data==1&vxyh(:,2)>0+vstep(i)&vxyh(:,2)<0.5+vstep(i);    
    subplot(1,numel(vstep),i);  
    imagesc(fh,1:32,sq(nansum(nanyhd(ind,:,1:numel(chans))))'./sum(ind));
    caxis([0,.6]);
    %colormap jet    
    colormap parula    
    xlabel('Frequency Hz')
end
suptitle(['steady state spec count: ',Trial.filebase,' (1-40Hz)'])





f = fh<50;
%f = 50<fh&fh<120;
figure();
ch = [1:4:32];
for c = 1:numel(ch),    
    subplot(1,numel(ch),c);
    ind = sum(nanyhd(:,f,ch(c)),2)>0 & logical(tper.data);
    hist2([vxyh(ind,1),thetaPhase(ind,6)],linspace(-2,2,20),linspace(-pi,pi,20));
    colormap jet
end

bins = discretize(vxyh(:,2),linspace(-2,2,30));

gcount = [];
for f = 1:16,
    for b = 1:29,
        for c = 1:32,
        gcount(b,f,c) = sum(nanyhd(bins==b,f,c));
        end
    end
end

figure,
imagesc(linspace(-2,2,30),fh(1:16),gcount(:,:,10)');
        
    



figure,imagesc(th,fh,10*log10(tyhd(:,:,30))')



figure,
imagesc(1:sum(nniz(thetaPhase)),1:numel(chans),thetaPhase(nniz(thetaPhase),:)');
colormap hsv





% $$$ % nan thresholds
% $$$ dnanThresh = [prctile(reshape(10.*log10(yhd.data),[],size(yld,3)),90)'];
% $$$ %             prctile(reshape(10.*log10(yhd.data),[],size(yld,3)),90)'];
% $$$ nanyhd = yhd.copy;
% $$$ nanyhd.data = 10*log10(nanyhd.data)>repmat(permute(dnanThresh(:,1),[3,2,1]),yhd.size([1,2]));


%% plot diagnostic time resolved psd and patchSD
sp = [];
figure,
sp(1) = subplot(211);
imagesc(th,fh,nyhd(:,1:32,18)'),
caxis([-1,3]),
colormap jet, 
axis xy
sp(2) = subplot(212);
imagesc(th,fh,nanyhd(:,1:32,18)'),
caxis([0,1]),
colormap jet, 
axis xy
linkaxes(sp,'xy')



hfig = figure(93993498);clf
hfig.Units = 'centimeters';
hfig.Position = [3,10,45,4];
vstep = -2:0.5:1.5;
sts = 'lrp';
for s = 1:numel(sts),
for i = 1:numel(vstep)
    ind = stc{sts(s)}+[0.2,-0.2];
    ind.cast('TimeSeries',yhd);
    ind = ind.data==1&vxyh(:,2)>0+vstep(i)&vxyh(:,2)<0.5+vstep(i);    
    subplot2(3,numel(vstep),s,i);  
    imagesc(fh,1:32,sq(nansum(nanyhd(ind,:,1:numel(chans))))'./sum(ind))    
    imagesc(fh,1:32,sq(nansum(nanyhd(ind,:,1:numel(chans))))'./sum(ind))
    caxis([0,1]);
    colormap jet    
    %colormap parula
end
end









sts='r'
evt=2;
figure,cnt=1;
for i=linspace(-1,1,80).*yld.sampleRate
    imagesc(fl,1:numel(chans),...
            sq(nanmean(nyld(round(stc{sts,yld.sampleRate}.data(diff(stc{sts}.data,1,2)>150,evt)+i),:,1:numel(chans)),1))'),
    caxis([-2,2])
    text( 35,13,num2str((i)./yld.sampleRate),'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7])
    colorbar
    pause(.1)        
    frms(cnt) = getframe;
    cnt=cnt+1;
end

prm.fps = 10;
prm.loop = inf;
makeGif(frms,['/storage/gravio//ownCloud/MjgEdER2016/spect_rear_trans',num2str(evt),'_1:40.gif'],prm);


%% plot bhv triggered clipped psd 
OwnDir = '/storage/gravio/nextcloud/MjgER2016/';
FigDir = ['MeanPSD_bhv_',Trial.filebase];
for sts='wrnpms',
    for evt = 1:2,
        try,
            tlabel = {'onset','offset'};
            sper = stc{sts,yld.sampleRate};
            evttime = sort(sper.data(:,evt));
            evttime = unique(evttime(diff(stc{sts}.data,1,2)>150));
            ytrns = reshape(GetSegs(yld.data,round(evttime-1*yld.sampleRate),round(2*yld.sampleRate)),[],numel(evttime),numel(fl),numel(chans));

            set(0,'defaultAxesFontSize',8,...
                  'defaultTextFontSize',8)
            
            hfig = figure(82747472),clf

            set(hfig,'units','centimeters')
            set(hfig,'Position',[18,5,7,20])
            set(hfig,'PaperPositionMode','auto');
            
            imagesc(linspace(-1,1,10),1:numel(chans),sq(nanmean(nanmean(ytrns(:,:,6<fl&fl<13,:),2),3))')
            xlabel('Time(s)');
            ylabel('Channels');    
            title({['Mean theta power(6-13Hz) '],['@ ',sper.label,' ',tlabel{evt}]});
            cax = colorbar;
            caxis([-max(abs(caxis)),max(abs(caxis))])    
            ylabel(cax,'z-score');

            print(hfig,'-depsc2',fullfile(OwnDir,...
                                          [Trial.filebase,'-meanPSD_6-13Hz_',sper.label,'_',tlabel{evt},'.eps']));
            print(hfig,'-dpng',  fullfile(OwnDir,...
                                          [Trial.filebase,'-meanPSD_6-13Hz_',sper.label,'_',tlabel{evt},'.png']));
        end
    end
end



[Co,f] = Comodugram(Trial,lfp,{'walk','rear','turn','pause','groom','sit'},'mtchglong');
figure,imagesc(f,f,Co(:,:,12,33,2)')
axis xy

figure
plot([nanmean(nanmean(yld(:,6<fl&fl<12,[6,7]),3),2),nanmean(nanmean(yld(:,6<fl&fl<12,[16:20]),3),2)])

chanS = [7:8];
freqS = [6,12];    
freqS = [12,24];        
chanD = [16:19];
freqD = [6,12];
freqD = [12,24];

shift = [.30,-0.30];

figure,hold on
ind = Trial.stc{'w'}+shift;
ind.clean;
[~,sind] = sort(ind.data(:,1));
ind.data = ind.data(sind,:);
plot(nanmean(nanmean(nyld(ind,freqS(1)<fl&fl<freqS(2),chanS),3),2),...
     nanmean(nanmean(nyld(ind,freqD(1)<fl&fl<freqD(2),chanD),3),2),'.c');
[rho,pval] = corr([nanmean(nanmean(nyld(ind,freqS(1)<fl&fl<freqS(2),chanS),3),2),...
                   nanmean(nanmean(nyld(ind,freqD(1)<fl&fl<freqD(2),chanD),3),2)])
ind = Trial.stc{'r'}+shift;
ind.clean;    
plot(nanmean(nanmean(nyld(ind,freqS(1)<fl&fl<freqS(2),chanS),3),2),...
     nanmean(nanmean(nyld(ind,freqD(1)<fl&fl<freqD(2),chanD),3),2),'.r');
[rho,pval] = corr([nanmean(nanmean(nyld(ind,freqS(1)<fl&fl<freqS(2),chanS),3),2),...
                   nanmean(nanmean(nyld(ind,freqD(1)<fl&fl<freqD(2),chanD),3),2)])


figure,hold on
ind = Trial.stc{'w'}+shift;
ind.clean;    
[~,sind] = sort(ind.data(:,1));
ind.data = ind.data(sind,:);
plot((nanmean(nanmean(nyld(ind,freqS(1)<fl&fl<freqS(2),chanS),3),2)),...
     (nanmean(nanmean(nyld(ind,freqD(1)<fl&fl<freqD(2),chanD),3),2)),'.c');
ind = Trial.stc{'r'}+shift;
ind.clean;    
plot((nanmean(nanmean(nyld(ind,freqS(1)<fl&fl<freqS(2),chanS),3),2)),...
     (nanmean(nanmean(nyld(ind,freqD(1)<fl&fl<freqD(2),chanD),3),2)),'.r');



figure,hold on

ind = Trial.stc{'t'};
[~,sind] = sort(ind.data(:,1));
ind.data = ind.data(sind,:);    
plot(nanmean(nanmean(nyld(ind,6<fl&fl<12,[6,7]),3),2),...
     nanmean(nyld(ind,6<fl&fl<12,end),2),'.c');
[rho,pval] = corr([nanmean(nanmean(nyld(ind,6<fl&fl<12,[6,7]),3),2),...
                   nanmean(nyld(ind,6<fl&fl<12,end),2)])

ind = Trial.stc{'r'};
plot(nanmean(nanmean(nyld(ind,6<fl&fl<12,[16:20]),3),2),...
     nanmean(nyld(ind,6<fl&fl<12,end),2),'.r');
[rho,pval] = corr([nanmean(nanmean(nyld(ind,6<fl&fl<12,[6,7]),3),2),...
                   nanmean(nyld(ind,6<fl&fl<12,end),2)])





figure,hold on
ind = Trial.stc{'w'};
[~,sind] = sort(ind.data(:,1));
ind.data=  ind.data(sind,:);
%plot(tpow(ind,1),tpow(ind,2),'.c');        
plot(log10(tpow(ind,1)),log10(tpow(ind,2)),'.c');        
[rho,pval] = corr(tpow(ind,:))
ind = Trial.stc{'r'};
%plot(tpow(ind,1),tpow(ind,2),'.r');    
plot(log10(tpow(ind,1)),log10(tpow(ind,2)),'.r');            
[rho,pval] = corr(tpow(ind,:))





figure,
imagesc(tl,fl,yld(:,:,33)')
axis xy,
colormap jet,
caxis([-.5,.2])    

OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';
for sts='wrnpms',
    for evt = 1:2,
        try,
            tlabel = {'onset','offset'};
            sper = stc{sts,yld.sampleRate};
            evttime = sort(sper.data(:,evt));
            evttime = unique(evttime(diff(stc{sts}.data,1,2)>100));
            ytrns = yld(evttime,:,:,:);

            set(0,'defaultAxesFontSize',8,...
                  'defaultTextFontSize',8)
            
            hfig = figure(82747472),clf

            set(hfig,'units','centimeters')
            set(hfig,'Position',[ 18    13    18    12]);
            set(hfig,'PaperPositionMode','auto');
            
            imagesc(fl,1:numel(chans),sq(nanmean(ytrns))');
            xlabel('Frequency (Hz)');
            ylabel('Channels');    
            title({['Mean PSD '],['@ ',sper.label,' ',tlabel{evt}]});
            cax = colorbar;
            caxis([-max(abs(caxis)),max(abs(caxis))])    
            ylabel(cax,'z-score');

            print(hfig,'-depsc2',fullfile(OwnDir,...
                                          [Trial.filebase,'-meanPSD',sper.label,'_',tlabel{evt},'.eps']));
            print(hfig,'-dpng',  fullfile(OwnDir,...
                                          [Trial.filebase,'-meanPSD',sper.label,'_',tlabel{evt},'.png']));
        end
    end
end


OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';
for sts='wrnpms',
    for evt = 1:2,
        try,
            tlabel = {'onset','offset'};
            sper = stc{sts,yhd.sampleRate};
            evttime = sort(sper.data(:,evt));
            evttime = unique(evttime(diff(stc{sts}.data,1,2)>100));
            ytrns = log10(yhd(evttime,:,1:numel(chans)));
            ytrns = nyhd(evttime,:,1:numel(chans));
            %ytrns(yhd(evttime,:,:,:)==0,:,:)=[];
            

            set(0,'defaultAxesFontSize',8,...
                  'defaultTextFontSize',8)
            
            hfig = figure(82747472),clf

            set(hfig,'units','centimeters')
            set(hfig,'Position',[ 18    13    18    12]);
            set(hfig,'PaperPositionMode','auto');
            
            imagesc(fh,1:numel(chans),sq(nanmean(ytrns))');
            xlabel('Frequency (Hz)');
            ylabel('Channels');    
            title({['Mean PSD '],['@ ',sper.label,' ',tlabel{evt}]});
            cax = colorbar;
            %caxis([-max(abs(caxis)),max(abs(caxis))])    
            ylabel(cax,'z-score');

            print(hfig,'-depsc2',fullfile(OwnDir,...
                                          [Trial.filebase,'-meanPSD40-120',sper.label,'_',tlabel{evt},'.eps']));
            print(hfig,'-dpng',  fullfile(OwnDir,...
                                          [Trial.filebase,'-meanPSD40-120',sper.label,'_',tlabel{evt},'.png']));
        end
    end
end


sts='r'
evt=2;
figure,cnt=1;
for i=linspace(-2,2,200).*yhd.sampleRate
    imagesc(fh,1:yhd.size(3),sq(nanmean(yhd(round(Trial.stc{sts,yhd.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)+i-tshift),:,:),1))'),
    caxis([-.51,.51])
    text( 80,30,num2str(i./yhd.sampleRate),'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7])
    pause(.2)
    frms(cnt) = getframe;
    cnt=cnt+1;
end

prm.fps = 10;
prm.loop = inf;
makeGif(frms,'/gpfs01/sirota/bach/homes/gravio/figures/spect_rear_off_gamma.gif',prm);


figure,
s = 'r'
sp1 = subplot(121);hold on
sp2 = subplot(122);hold on
colors= jet(numel(chans))';

for c = colors
    plot(sp1,fl,-.5*find(ismember(colors',c','rows'))+nanmean(log10(yld(Trial.stc{s,yld.sampleRate}.data(2:end-1,:),:,ismember(colors',c','rows')))),'color',c)
    plot(sp2,fh,-.21*find(ismember(colors',c','rows'))+nanmean(log10(yhd(Trial.stc{s,yhd.sampleRate},:,ismember(colors',c','rows')))),'color',c)
end


figure,
%subplot(121),
chan = 15;
st1 = 'l'; st2 = 'g'; 
st1 = 'r'; st2 = 'w';
boundedline(fl,mean(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),2*std(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),'-r',fl,mean(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),2*std(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),'-b','alpha')

figure,
%subplot(121),
st1 = 'l'; st2 = 'g';
boundedline(fl,mean(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),2*std(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),'-r',fl,mean(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),2*std(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),'-b','alpha')


subplot(122),

boundedline(fh,mean(log10(yhd(Trial.stc{'r',yhd.sampleRate}.data(2:end-1,:),:,3))),std(log10(yhd(Trial.stc{'r',yhd.sampleRate}.data(2:end-1,:),:,3))),'-r',fh,mean(log10(yhd(Trial.stc{'w',yhd.sampleRate}.data(2:end-1,:),:,3))),std(log10(yhd(Trial.stc{'w',yhd.sampleRate}.data(2:end-1,:),:,3))),'-b','alpha')


evt = 1;
sts = 'r';
figure,imagesc(sq(nanmean(GetSegs(sq(nanmean(yld.data(:,fl<12&fl>6,:),2)),round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yld.sampleRate),round(4*yld.sampleRate),0),2))')

evt = 2;
sts = 'r';
figure,imagesc(sq(nanmean(GetSegs(sq(nanmean(yld.data(:,fl<25&fl>18,:),2)),round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yld.sampleRate),round(4*yld.sampleRate),0),2))')


evt = 1;
sts = 'r';
figure,imagesc(sq(nanmean(GetSegs(sq(nanmean(yhd.data(:,fh<120&fh>80,:),2)),round(Trial.stc{sts,yhd.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yhd.sampleRate),round(4*yhd.sampleRate),0),2))')

tdr = sq(mean(yl(:,fl>6&fl<12,:),2)./(mean(yl(:,fl<5,:),2)+mean(yl(:,fl<18&fl>13,:),2))/2);
tdr = MTADlfp('data',tdr,'sampleRate',1/diff(tl(1:2,1)));


tdr.data = (tdr.data-repmat(nanmedian(tdr.data),[tdr.size(1),1,1]))./repmat(nanstd(tdr.data),[tdr.size(1),1,1]);


evt = p1;
sts = 'r';
figure,imagesc(sq(nanmean(GetSegs(tdr.data,round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yld.sampleRate-tshift/yld.sampleRate),round(4*yld.sampleRate),0),2))')

sts = 'w';
[U,V,D] = svd(cov(yld(round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-tshift/yld.sampleRate),:,21)));
figure,imagesc(fl,fl,U),axis xy




ind = Trial.stc{'w'};
tpow = MTADxyz('data',sq(nanmean(yld(:,fl>6&fl<12,:),2)),'sampleRate',yld.sampleRate);
figure,
imagesc(sq(nanmean(tpow.segs(ind(:,1)-5,10),2))')
caxis([1,3])    
