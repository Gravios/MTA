

MjgER2016_load_data();

chans = [57:64];
freqRange = [50,1200];

% PREPROC xyz



Trial = MTATrial.validate('jg05-20120312.cof.all');

xyz = preproc_xyz(Trial,'trb');

ufr = Trial.load('ufr',xyz,[],units{20},0.02);



dat = Trial.lfp.copy();
dat.filename = [Trial.name,'.dat'];
dat.label = 'dat';
dat.ext = 'dat';
dat.sampleRate = Trial.sampleRate;
dat.load(Trial,chans,[],[4300,4500]);

sdat = dat.copy();
sdat.data = sdat.data(:,[2,4,5,7]);

specArgsGamma = struct('nFFT',2^11,...
                       'Fs',  dat.sampleRate,...
                       'WinLength',2^10,...
                       'nOverlap',2^10*0.875,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',freqRange);

[ysg,fsg,tsg] = fet_spec(Trial,sdat,'mtcsdglong',true,[],specArgsGamma,[],true);

specArgsTheta = struct('nFFT',2^16,...
                       'Fs',  dat.sampleRate,...
                       'WinLength',2^15,...
                       'nOverlap',2^15*0.875,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',[1,80]);

[yst,fst,tst] = fet_spec(Trial,sdat,'mtcsdglong',true,[],specArgsTheta,[],true);


specArgsM = struct('nFFT',2^14,...
                       'Fs',  dat.sampleRate,...
                       'WinLength',2^13,...
                       'nOverlap',2^13*0.875,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',[20,120]);

[ysm,fsm,tsm] = fet_spec(Trial,sdat,'mtcsdglong',true,[],specArgsM,[],false);



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
plot((1:size(dat,1))./dat.sampleRate,bsxfun(@plus,dat.data,1000:2000:16000));
plot((1:size(dat,1))./dat.sampleRate,sum(diff(bsxfun(@plus,dat.data,1000:2000:16000),1,2),2)-2e4);
set(gca,'XTickMode','auto')
set(gca,'XTickLabelMode','auto')
set(gca,'YTickMode','auto')
set(gca,'YTickLabelMode','auto')
Lines([],-1.6e4,'k');
linkaxes(sp,'x');


figure,     
spo = gobjects([0,1]);
spo(end+1) = subplot(511);
imagesc(tsg,fsg,angle(ysg(:,:,1,end)'));
colormap(spo(end),'hsv');
caxis([-pi,pi]);
Lines([],150,'m');
Lines([],250,'m');
ylim([100,600])
axis('xy');
spo(end+1) = subplot(512);
imagesc(tsg,fsg,log10(ysg(:,:,3,3)'));
caxis([-1,1.8]);
colormap(spo(end),'bone');
Lines([],150,'m');
Lines([],250,'m');
axis('xy');
ylim([100,600]);
spo(end+1) = subplot(513);
imagesc(tsm,fsm,log10(ysm(:,:,3)'));
caxis([-1,0])
colormap(spo(end),'bone');
Lines([],30,'m');
axis('xy');
spo(end+1) = subplot(514);
imagesc(tst,fst,log10(yst(:,:,3,3)'));
caxis([-1.4,0]);
colormap(spo(end),'bone');
Lines([],30,'m');
axis('xy');
spo(end+1) = subplot(515);hold('on');
plot((1:size(dat,1))./dat.sampleRate,sum(diff(bsxfun(@plus,dat.data,1000:2000:16000),1,2),2)-2.2e4);
plot((1:size(dat,1))./dat.sampleRate,bsxfun(@plus,dat.data,1000:2000:16000));
Lines([],-1.8e4,'k');
linkaxes(spo,'x');


unitSet = units{20};
spk = Trial.spk.copy();
spk.create(Trial,Trial.sampleRate,[],unitSet,'');
sWidth = 10./dat.sampleRate;

figure();
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


plot(linspace(0,200,round(200*ufr.sampleRate)+2),sum(ufr(round(4300*ufr.sampleRate): ...
                                                  round(4500*ufr.sampleRate),:),2)*5000-3e4);

