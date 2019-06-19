

MjgER2016_load_data();

chans = [57:64,68:2:90];
freqRange = [50,1200];

% PREPROC xyz


Trial = MTATrial.validate('jg05-20120312.cof.all');

xyz = preproc_xyz(Trial,'trb');

ufr = Trial.load('ufr',xyz,[],[],0.02);


dat = Trial.lfp.copy();
dat.filename = [Trial.name,'.dat'];
dat.label = 'dat';
dat.ext = 'dat';
dat.sampleRate = Trial.sampleRate;
dat.load(Trial,chans,[],[4300,4600]);

sdat = dat.copy();
sdat.data = sdat.data(:,[2,4,5,7,9:18]);

specArgsGamma = struct('nFFT',2^11,...
                       'Fs',  dat.sampleRate,...
                       'WinLength',2^10,...
                       'nOverlap',2^10*0.875,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',freqRange);

[ysg,fsg,tsg] = fet_spec(Trial,sdat,'mtcsdglong',true,[],specArgsGamma,[],false);

specArgsTheta = struct('nFFT',2^16,...
                       'Fs',  dat.sampleRate,...
                       'WinLength',2^15,...
                       'nOverlap',2^15*0.875,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',[1,80]);

[yst,fst,tst] = fet_spec(Trial,sdat,'mtcsdglong',true,[],specArgsTheta,[],false);


specArgsM = struct('nFFT',2^14,...
                       'Fs',  dat.sampleRate,...
                       'WinLength',2^13,...
                       'nOverlap',2^13*0.875,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',[20,120]);

[ysm,fsm,tsm] = fet_spec(Trial,sdat,'mtcsdglong',true,[],specArgsM,[],false);


tdat = dat.copy();
tdat.data = sum(diff(bsxfun(@plus,dat.data(:,1:8),1000:2000:16000),1,2),2);
specArgsM = struct('nFFT',2^14,...
                       'Fs',  dat.sampleRate,...
                       'WinLength',2^13,...
                       'nOverlap',2^13*0.875,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',[20,120]);

[yslm,fslm,tslm] = fet_spec(Trial,tdat,'mtcsdglong',true,[],specArgsM,[],false);
specArgsM = struct('nFFT',2^16,...
                       'Fs',  dat.sampleRate,...
                       'WinLength',2^15,...
                       'nOverlap',2^15*0.875,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',[1,40]);

[yslt,fslt,tslt] = fet_spec(Trial,tdat,'mtcsdglong',true,[],specArgsM,[],false);

figure,
subplot(411);imagesc(tslm,fslm,log10(yslm(:,:)'));axis('xy');
subplot(412);imagesc(tslt,fslt,log10(yslt(:,:)'));axis('xy');caxis([-3,0]);
subplot(413);imagesc(tst,fst,log10(yst(:,:,1)'));axis('xy'); caxis([-2,0]); ylim([1,30]);
subplot(414);imagesc(tst,fst,log10(yst(:,:,4)'));axis('xy'); caxis([-2,0]); ylim([1,30]);
linkaxes(findobj(gcf,'Type','Axes'),'x');






% GOOD Feature for immobile theta detection
figure,
plot(log10(sum(yslt(:,fslt<20),2)),log10(yslt(:,15)./yslt(:,5)),'.'),
hold('on');
ind = 1:280;,plot(log10(sum(yslt(ind,fslt<20),2)),log10(yslt(ind,15)./yslt(ind,5)),'.r','Markersize',10)
ind = 812:820;plot(log10(sum(yslt(ind,fslt<20),2)),log10(yslt(ind,15)./yslt(ind,5)),'.g','Markersize',10)
ind = 1592:1632;,plot(log10(sum(yslt(ind,fslt<20),2)),log10(yslt(ind,15)./yslt(ind,5)),'.g','Markersize',10)
ind = 1960:2056;,plot(log10(sum(yslt(ind,fslt<20),2)),log10(yslt(ind,15)./yslt(ind,5)),'.g','MarkerSize',10)
ind = unique(round(tsg(log10(ysg(:,9))>0.5).*yslt.sampleRate));
plot(log10(sum(yslt(ind,fslt<20),2)),log10(yslt(ind,15)./yslt(ind,5)),'.c','MarkerSize',10)

figure,
plot(log10(sum(yst(:,fslt<20),2)),log10(yst(:,15)./yst(:,5)),'.'),
hold('on');
ind = 1:280;,plot(log10(sum(yst(ind,fslt<20),2)),log10(yst(ind,15)./yst(ind,5)),'.r','Markersize',10)
ind = 812:820;plot(log10(sum(yst(ind,fslt<20),2)),log10(yst(ind,15)./yst(ind,5)),'.g','Markersize',10)
ind = 1592:1632;,plot(log10(sum(yst(ind,fslt<20),2)),log10(yst(ind,15)./yst(ind,5)),'.g','Markersize',10)
ind = 1960:2056;,plot(log10(sum(yst(ind,fslt<20),2)),log10(yst(ind,15)./yst(ind,5)),'.g','MarkerSize',10)
ind = unique(round(tsg(log10(ysg(:,9))>0.5).*yst.sampleRate));
plot(log10(sum(yst(ind,fslt<20),2)),log10(yst(ind,15)./yst(ind,5)),'.c','MarkerSize',10)



gdat = dat.copy();
gdat.data = -sum(diff(bsxfun(@plus,dat(:,1:8),1000:2000:16000),1,2),2);
gdat.resample(ysg);
mdat = dat.copy();
mdat.data = mean(bsxfun(@plus,dat(:,1:8),1000:2000:16000),2);
mdat.resample(ysg);

ind = log10(ysg(:,9))>0.5;
ind = log10(ysg(:,9))>1;
ind = log10(ysg(:,9))>1.5;
ind = ':';
figure()
hold('on');
plot(gdat(ind),mdat(ind),'.g','MarkerSize',10)



ind = 1100:1200;
figure()
hold('on');
plot(gdat(ind),mdat(ind),'.b','MarkerSize',10)
figure,
plot(cos(gdat(ind)./mdat(ind)),'.b','MarkerSize',10)
figure,plot(cos(gdat(:)./mdat(:)))
figure,plot(nunity([gdat.data,mdat.data,sdat.data(:,1)]))
figure,plot([-nunity(-nunity(gdat.data)),nunity(sdat.data(:,1))])

figure,plot(gdat.data,log10(ysg(:,9)),'.')
hold on,plot(gdat.data,log10(ysg(:,19)),'.')
figure();
hold('on');
for j = 1:4,
    plot(log10(ysg(:,9,1)),log10(ysg(:,19,j)),'.')
end


figure,
sp = tight_subplot(size(sdat,2)*2,1,0,0.05);
for c = 1:size(sdat,2)
    axes(sp(c*2-1));
    imagesc(tsg,fsg,log10(ysg(:,:,c,c)'));
    axis('xy');
    caxis([-1,1.8]);
    colormap(sp(c*2-1),'jet');
    
    if c~=size(sdat,2),
        axes(sp(c*2));
        imagesc(tsg,fsg,angle(ysg(:,:,c,c+1)'));
        axis('xy');
        caxis([-pi,pi]);
        colormap(sp(c*2),'hsv');
    end
    
end
axes(sp(end));
hold('on');
plot((1:size(dat,1))./dat.sampleRate,bsxfun(@plus,dat.data(:,1:8),1000:2000:16000));
plot((1:size(dat,1))./dat.sampleRate,sum(diff(bsxfun(@plus,dat.data(:,1:8),1000:2000:16000),1,2),2)-2e4);
set(gca,'XTickMode','auto')
set(gca,'XTickLabelMode','auto')
set(gca,'YTickMode','auto')
set(gca,'YTickLabelMode','auto')
Lines([],-1.6e4,'k');
linkaxes(sp,'x');



unitSet = units{20};
spk = Trial.spk.copy();
spk.create(Trial,Trial.sampleRate,[],unitSet,'');
sWidth = 10./dat.sampleRate;

[spk.res,rind]  = SelectPeriods(spk.res,[4300,4600].*Trial.sampleRate,'d',1,1);
spk.clu = spk.clu(rind);




figure,     
spo = gobjects([0,1]);
spo(end+1) = subplot(511);
hold('on');
unitClr = repmat('b',[1,numel(unitSet)]);
for u = unitSet
    uind = find(u==unitSet);
    res = spk(u);
    cts = 1:max(res);
    sti = ismember(cts,res);
    patch(reshape(repmat(reshape(bsxfun(@plus,cts(sti)/dat.sampleRate,[-sWidth/2;sWidth/2]),[],1)',2,1),[],1),...
          repmat([0,1,1,0],[1,numel(res)])'+uind,...
          unitClr(uind),'EdgeColor',unitClr(uind));    
end
Lines([],1:numel(unitSet),'w');
% $$$ imagesc(tsg,fsg,log10(ysg(:,:,1)'));
% $$$ colormap(spo(end),'bone');
% $$$ caxis([-1,1.5]);
% $$$ Lines([],150,'m');
% $$$ Lines([],250,'m');
% $$$ ylim([100,600])
% $$$ axis('xy');
spo(end+1) = subplot(512);
imagesc(tsg,fsg,log10(ysg(:,:,3)'));
caxis([-1,1.5]);
colormap(spo(end),'bone');
Lines([],150,'m');
Lines([],250,'m');
axis('xy');
ylim([100,600]);
spo(end+1) = subplot(513);
imagesc(tsm,fsm,log10(ysm(:,:,1)'));
caxis([-1,0])
colormap(spo(end),'bone');
Lines([],30,'m');
axis('xy');
spo(end+1) = subplot(514);
imagesc(tst,fst,log10(yst(:,:,1)'));
caxis([-1.4,0]);
colormap(spo(end),'bone');
Lines([],30,'m');
ylim([0,30]);
axis('xy');
spo(end+1) = subplot(515);hold('on');
plot((1:size(dat,1))./dat.sampleRate,bsxfun(@plus,dat(:,1:8),1000:2000:16000));
plot((1:size(dat,1))./dat.sampleRate,-sum(diff(bsxfun(@plus,dat(:,1:8),1000:2000:16000),1,2),2)+2.2e4);
%plot((1:size(dat,1))./dat.sampleRate,sum(diff(bsxfun(@plus,dat.data(:,1:8),1000:2000:16000),1,2),2)-2.2e4);
plot((1:size(dat,1))./dat.sampleRate,bsxfun(@minus,dat(:,9:20)/2,1000:4000:48000)-2e4);
Lines([],-1.8e4,'k');
linkaxes(spo,'x');
ind = round(4300*ufr.sampleRate):round(4600*ufr.sampleRate);
plot(linspace(0,300,numel(ind)),sum(ufr(ind,units{20}),2)*5000+2e4);
plot(linspace(0,300,numel(ind)),sum(ufr(ind,unitsInt),2)*2000+2e4,'r');






xlim([199,203])
xlim([100.5,104.5])
xlim([113,118])
xlim([100,106])
xlim([14,20])
xlim([20,26])
xlim([23,29])
xlim([26,32])
xlim([28,34])
xlim([32,38])
xlim([34,40])
xlim([0,20]);
xlim([36,43])