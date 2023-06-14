

global AP

configure_default_args();
MjgER2016_load_data();

Trial = Trials{7};
unitSubset = units{7};

rfq = fet_respiration_freq(Trial,[],65,'',false);
hzp = fet_hzp(Trial);
ncp = resample(Trial.load('lfp',65),250);
ncpFilt = filter(copy(ncp),'ButFilter',4,[2,16],'bandpass');

spk = Trial.load('spk',250,'walk+turn+pause',unitSubset);


figure();
plot(rfq.data);
hold('on');
plot(nunity(ncp.data));
plot(nunity(ncpFilt.data));

figure();
for u = unitSubset
clf();
hold('on');,
plot(rp(:,1))
plot(rp(:,2),'c');
plot(spk(u),rp(spk(u),2),'*m');
plot(nunity(ncpFilt.data)./4);
title(num2str(u));
waitforbuttonpress();
end


ind = [Trial.stc{'gper-sit-groom'}];
ind = [Trial.stc{'hloc&theta'}];
ind = [Trial.stc{'lloc&theta'}];



figure,
histogram(hzp(ind,2),linspace(0,280,100));

figure,
hist2([hzp(ind,:)],linspace(-2,0.8,100),linspace(20,250,100));
caxis([0,200]);

rp = fet_rfqXhp(Trial,'overwrite',true);
AP.fet_HB_pitchB.referenceTrial = 'Ed05-20140529.ont.all';
% compute_bhv_ratemaps -----------------------------------------------------------------------------
AP.compute_bhv_ratemaps =                                                                        ...
    struct('get_featureSet',            @fet_rfqXhp,                                             ...
           'sampleRate',                250,                                                     ...
           'pfsArgs',                   struct('states',           'walk+turn',                  ...
                                               'binDims',          [0.1,0.5],                    ...
                                               'SmoothingWeights', [1.8,1.8],                    ...
                                               'numIter',          1,                            ...
                                               'boundaryLimits',   [-2.5,0.8;0,14],              ...
                                               'halfsample',       false),                       ...
           'threshRate',                0.8,                                                     ...
           'threshDist',                250                                                      ...
           );
%---------------------------------------------------------------------------------------------------
bfs = compute_bhv_ratemaps(Trial,unitSubset);


figure();
for u = unitSubset,
    plot(bfs,u,1,'colorbar',[],false,'colorMap',@jet);
    title(num2str(u));
    waitforbuttonpress();
end

xyz = preproc_xyz(Trial);
ang = create(MTADang,Trial,filter(preproc_xyz(Trial),'ButFilter',4,20,'low'));

ang.filter('ButFilter',4,1,'high');


xts = [1:size(xyz)]./xyz.sampleRate;
lts = [1:size(rlfp)]./rlfp.sampleRate;



lfp = Trial.load('lfp',[57,64]);

rlfp = copy(lfp);

rlfp.data = diff(lfp.data,1,2);

[ys,fs,ts] = fet_spec(Trial,rlfp,'mtchglong',false,[],struct('nFFT',2^11,'Fs',rlfp.sampleRate,...
                            'WinLength',2^10,'nOverlap',2^10*.875,...
                            'FreqRange',[1,30]));

[mys,mfs,mts] = fet_spec(Trial,rlfp,'mtchglong',false,[],struct('nFFT',2^10,'Fs',rlfp.sampleRate,...
                            'WinLength',2^9,'nOverlap',2^9*.875,...
                            'FreqRange',[1,30]));

rhm = fet_rhm(Trial);

figure();
subplot(5,1,[1,2]);
    imagesc(mts,mfs,log10(mys.data)');
    axis('xy');
    colormap('jet');
    caxis([7,9])
subplot(5,1,[3,4]);
    hold('on');
    %plot(xts,MedianFilter(nunity(ang(:,3,4,2)),265)*3)
    plot(xts,MedianFilter(nunity(ang(:,3,4,3)),256)*3)
    plot(xts,MedianFilter(nunity(ang(:,2,4,3)),256)*3)    
    %plot(xts,MedianFilter(nunity(ang(:,4,5,3)),256)+0.05,'m')
    plot(xts,MedianFilter(nunity(ang(:,3,7,3)),256)-0.05,'g');
% $$$     plot(ts,nunity(log10(ys(:,5)))'.*0.01,'k');
    plot(ts,nunity(mean(log10(ys(:,fs<20)),2))'.*0.025,'r');  
    plot(ts,nunity((mean(log10(ys(:,fs>20)),2)./sum(log10(ys.data),2)))'.*0.01,'k');      
    %plot(xts,MedianFilter(nunity(ang(:,5,7,2)),256)+0.1,'m')    
    plot(xts,nunity(rhm(:,1))/20+0.1,'r')
    plot(lts,nunity(rlfp(:,1))./50+0.15,'c')
    plot(lts,nunity(MedianFilter(rlfp(:,1),64))./50+0.15,'k')
subplot(515);
    plotSTC(Trial.stc,1);
    ylim([0.5,7.5]);
linkx();


