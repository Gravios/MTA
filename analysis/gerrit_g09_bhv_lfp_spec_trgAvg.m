
linkSession('g10-20130415',...
            '/storage/gerrit/data/project/Behavior/xyz',...
            '/storage/gerrit/data/prep/g10')
Session = MTASession('g10-20130415','cof',true,'0x0800','xyzSampleRate',119.881035);
QuickTrialSetup(Session,[25,-25],[4,5]);

Trial = MTATrial('g10-20130415');

Trial = MTATrial('g09-20120328');
Trial.load('stc','auto_wbhr');
Trial.load('stc','LGR-hand_labeled_rev1-wrsnkm');

xyz = Trial.load('xyz').filter(gtwin(.1,Trial.xyz.sampleRate));

[rhm,rfs,rts] = fet_rhm(Trial,[],'wcsd');
rhm.data = nanmean(log10(rhm(:,rfs<14&rfs>6)),2);
xyz.resample(rhm);
ang = create(Trial.ang.copy,Trial,xyz);


Data = Trial.xyz.copy;
Data.data = [ang(:,'head_back','head_front',2),rhm(:)];
Data.sampleRate = rhm.sampleRate;

wper = Trial.stc{'h'};
wper.cast('TimeSeries');
wper.resample(rhm);

nind = nniz(Data)&wper.data==1;
figure,hist2(Data(nind,:),linspace(-1.5,1.5,50),linspace(-7,-2,50))
cpts = ClusterH2(Data(nind,:),linspace(-1.5,1.5,50),linspace(-7,-2,50));


nind(nind) = cpts;



hwalk = zeros(size(wper));
lwalk = zeros(size(wper));
lwalk(wper(:)==1) = nind==1;
hwalk(wper(:)==1) = nind==0;

wper = Trial.stc{'w'};
Trial.stc.states(Trial.stc.gsi('h')) = [];
Trial.stc.addState(Trial.spath,...
             Trial.filebase,...
             ThreshCross(hwalk==1,.5,1),...                   
             rhm.sampleRate,...
             Trial.sync.copy,...
             Trial.sync.data(1),...
             'hswalk','h','TimePeriods');
Trial.stc.states{Trial.stc.gsi('h')} = Trial.stc{'h',wper.sampleRate}&wper;
Trial.stc.states{Trial.stc.gsi('h')} = Trial.stc{'h'}+[1/wper.sampleRate,-1/wper.sampleRate];


Trial.stc.states(Trial.stc.gsi('l')) = [];
Trial.stc.addState(Trial.spath,...
             Trial.filebase,...
             ThreshCross(lwalk==1,.5,1),...
             rhm.sampleRate,...
             Trial.sync.copy,...
             Trial.sync.data(1),...
             'lswalk','l','TimePeriods');
Trial.stc.states{Trial.stc.gsi('l')} = Trial.stc{'l',wper.sampleRate}&wper;


% Check segmentation of low and high walk
edges = 0:.5:150;
figure,
subplot(211)
N = histc(xyz(Trial.stc{'w'},7,3),edges);
bar(edges,N,'histc');

subplot(212),hold on,
N = histc(xyz(Trial.stc{'l'},7,3),edges);
bar(edges,N,'histc');
n
N = histc(xyz(Trial.stc{'h'},7,3),edges);
bax = bar(edges,N,'histc');
set(bax,'FaceColor',[1,0,0]);
set(bax,'EdgeColor',[1,0,0]);
set(bax,'FaceAlpha',.5);

% Okay looks like crap

% Triggered avg lfp HPC MEC LEC

%clabels = {'MEC3','CA1rad','DGgc','LEC3'};%g09-20120328
%chans = [10,33,40,85]; %g09-20120328

clabels = {'MEC3','CA1pyr','CA1rad','LEC3'};%g09-20120328
chans = [24,38,40,75]; %g10-20130415
nchan = numel(chans);
states = 'rhl';
nsts = numel(states);

lfp = load(Trial,'lfp',chans);

sparm = struct('nFFT',2^10,...
               'Fs',lfp.sampleRate,...
               'WinLength',2^9,...
               'nOverlap',2^9*.875,...
               'FreqRange',[1,50]);
[ys,fs,ts] = fet_spec(Trial,lfp,'mtchglong',true,lfp.sampleRate,sparm,true);




hfig = figure(99901);
for s = 1:nsts,
    sper = Trial.stc.filter(ys.sampleRate,{states(s),{'exclusion',{'hswalk','lswalk','rear'},1},{'duration',1}});
    for c = 1:nchan,
        ygs = GetSegs(ys(:,:,c,c),sper{1}-round(2*ys.sampleRate),round(4*ys.sampleRate)+1,nan);
        subplot2(nchan,nsts,c,s);
        imagesc(linspace(-2,2,round(4*ys.sampleRate)+1),fs,sq(nanmean(log10(ygs(:,sum(abs(nunity(log10(ygs(:,:,10))')),2)<50,:)),2))');axis xy
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title({'Behavior Onset Triggered Spec Avg',...
               [clabels{c},' : ' Trial.stc{states(s)}.label]});
    end
end



hfig = figure(99902);
for s = 1:nsts,
    sper = Trial.stc.filter(ys.sampleRate,{states(s),{'exclusion',{'hswalk','lswalk','rear'},1},{'duration',1}});
    for c = 1:nchan,
        ygs = GetSegs(ys(:,:,c,c),sper{2}-round(2*ys.sampleRate),round(4*ys.sampleRate)+1,nan);
        subplot2(nchan,nsts,c,s);
        imagesc(linspace(-2,2,round(4*ys.sampleRate)+1),fs,sq(nanmean(log10(ygs(:,sum(abs(nunity(log10(ygs(:,:,10))')),2)<50,:)),2))');axis xy
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title({'Behavior Offset Triggered Spec Avg',...
               [clabels{c},' : ' Trial.stc{states(s)}.label]});
    end
end






clabels = {'MEC3','CA1rad','DGgc','LEC3'};
chans = [10,33,40,85];
nchan = numel(chans);
states = {'rear','hswalk','lswalk'};
nsts = numel(states);

lfp = load(Trial,'lfp',chans);

sparm = struct('nFFT',2^8,...
               'Fs',lfp.sampleRate,...
               'WinLength',2^7,...
               'nOverlap',2^7*.875,...
               'FreqRange',[30,150]);
[gys,gfs,gts] = fet_spec(Trial,lfp,'mtchglong',true,lfp.sampleRate,sparm,true);




hfig = figure(99903);
for s = 1:nsts,
    sper = Trial.stc.filter(ys.sampleRate,{states{s},{'exclusion',states,1},{'duration',.5}});
    for c = 1:nchan,
        ygs = GetSegs(ys(:,:,c,c),sper{1}-round(2*ys.sampleRate),round(4*ys.sampleRate)+1,nan);
        subplot2(nchan,nsts,c,s);
        imagesc(linspace(-2,2,round(4*ys.sampleRate)+1),fs,sq(nanmean(log10(ygs(:,mean(abs(nunity(log10(ygs(:,:,10))')),2)<.71,:)),2))');axis xy
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title({'Behavior Onset Triggered Spec Avg',...
               [clabels{c},' : ' states{s}]});
    end
end
%ForAllSubplots('caxis([2.5980,3.3503])')



hfig = figure(99904);
for s = 1:nsts,
    sper = Trial.stc.filter(ys.sampleRate,{states{s},{'exclusion',states,1},{'duration',1}});
    for c = 1:nchan,
        ygs = GetSegs(ys(:,:,c,c),sper{2}-round(2*ys.sampleRate),round(4*ys.sampleRate)+1,nan);
        subplot2(nchan,nsts,c,s);
        gqual = mean(abs(nunity(log10(ygs(:,:,10))')),2);
        gthresh = mean(gqual)+std(gqual)/2;
        imagesc(linspace(-2,2,round(4*ys.sampleRate)+1),fs,sq(nanmean(log10(ygs(:,gqual<gthresh,:)),2))');axis xy
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title({'Behavior Offset Triggered Spec Avg',...
               [clabels{c},' : ' states{s}]});
        ca = caxis;
        caxis([2,ca(2)]);
    end
end
 




