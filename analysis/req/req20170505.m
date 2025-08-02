% req20170505
%
%   Status: active
%   Type: Analysis 
%   Final_Forms: NA
%   Description: 
% 
%  Bugs: NA
 



Trial = MTATrial.validate('jg05-20120317.cof.all');
stc = Trial.load('stc','hand_labeled_rev3_jg');

xyz = Trial.load('xyz');

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,0.5,'low');
figure,plot([xyz(:,5,3),fxyz(:,5,3)])
hold on,plot((diff(fxyz(:,5,3)))*100)


fang = create(MTADang,Trial,fxyz);

lfp = Trial.load('lfp',65:5:96);
lfp.resample(xyz);


parspec = struct('nFFT',2^9,...
                 'Fs',  lfp.sampleRate,...
                 'WinLength',2^8,...
                 'nOverlap',2^8-1,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[1,20]);


[ys,fs,ts] = fet_spec(Trial,lfp,'mtchglong',true,[],parspec);

figure,imagesc(ts,fs,log10(ys(:,:,1,1))');colormap jet,axis xy,caxis([3,5])
figure,imagesc(ts,fs,ys(:,:,1,end)');colormap jet,axis xy,

c = 1;
thetaBand = [5,13];
thetaBandFrequencies = fs(5<fs&fs<13);
thetaPeakFreq = xyz.copy('empty');
[~,thetaPeakFreq.data] = max(log10(ys(:,5<fs&fs<13,c,c)),[],2);
thetaPeakFreq.data = thetaBandFrequencies(thetaPeakFreq.data);
thetaPeakFreq.data = thetaPeakFreq.data(1:size(xyz,1),:); 

figure,plot(abs(diff(fang(:,5,7,2))))
figure,plot(abs(diff(fxyz(:,1,3))))
[headZshiftIndex,headZshift] = LocalMinima(-abs(diff(fxyz(:,5,3))),90,-0.1);
%[headZshiftIndex,headZshift] = LocalMinima(-diff(fang(:,1,3)),180,-0.02);
[headZshiftIndex,sind] = SelectPeriods(headZshiftIndex,stc{'t'}.data,'d',1,0);
headZshift = headZshift(sind);
[headZshiftSorted,headZshiftSortIndex] = sort(headZshift);

thetaPeakFreq.filter('ButFilter',3,0.5,'low');
thetaFreqSegs = thetaPeakFreq.segs(headZshiftIndex-120,240);
dtpf = thetaPeakFreq.copy;
dtpf.data = [0;diff(thetaPeakFreq.data).*thetaPeakFreq.sampleRate];
dtpfs = dtpf.segs(headZshiftIndex-120,240);


vfz = fxyz.copy;
vfz.data = [0;diff(fxyz(:,5,3))];


ind = stc{'w+n+p+r+t'};

ind = stc{'t&w'};
figure,
plot(vfz(ind),dtpf(ind),'.');
hold on
ind = stc{'t&r'};
plot(vfz(ind),dtpf(ind)+3,'.');

figure,
subplot(211); plot(nanmean(thetaFreqSegs(:,headZshiftSortIndex)')))

figure,
subplot(211); plot(nanmean(thetaFreqSegs(:,headZshiftSortIndex)')))
subplot(212); plot(nanstd(diff(thetaFreqSegs(:,headZshiftSortIndex)')))

figure,
subplot(211); plot(nanmean(dtpfs(:,headZshiftSortIndex)'))
subplot(212); plot(nanstd(dtpfs(:,headZshiftSortIndex)'))

figure, imagesc(thetaFreqSegs(:,headZshiftSortIndex)')









figure,
sf = {@uminus,@(x) x};

for s = 1:2,
    [headZshiftIndex,headZshift] = LocalMinima(sf{s}(diff(fxyz(:,5,3))),240,-0.1);
    headZshiftIndex(headZshift<-0.3) = [];
    [headZshiftIndex,sind] = SelectPeriods(headZshiftIndex,stc{'t'}.data,'d',1,0);
    yss = ys.segs(headZshiftIndex-120,240,nan);    
    for c = 1:7,
        subplot2(2,7,s,c);
        %imagesc(1:240,fs,log10(sq(nanstd(yss(:,:,:,c,c),[],2)))'),axis xy,caxis([2,5])
        imagesc(1:240,fs,log10(sq(nanmean(yss(:,:,:,c,c),2)))'),axis xy,caxis([3,5])
    end
end

figure,
subplot(121);imagesc(log10(sq(yss(60,:,:,5,5)))')
subplot(122);imagesc(log10(sq(yss(180,:,:,5,5)))')

figure,
imagesc(1:240,fs,(sq(nanmean(yss(:,:,:,2,3),2)))'),axis xy
colormap jet,
caxis([0.5,1])


