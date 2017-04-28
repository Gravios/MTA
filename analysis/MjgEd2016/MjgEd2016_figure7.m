%% LFP state X theta analysis

Trial = MTATrial('jg05-20120310');
Stc = Trial.load('stc','NN0317');
labelTheta(Trial,Stc,72,true);

lfpChannels = 65:96;
embeddingWindow = round(0.5*xyz.sampleRate);

xyz = Trial.load('xyz');

Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',lfpChannels);
lfp.resample(xyz);

flfp = lfp.copy;
flfp.filter('ButFilter',3,[5,15],'bandpass');
flfp.unity;

rflfpOnset = flfp.segs(Stc{'r'}.data(:,1)-30,60,nan);

figure,
imagesc(rflfpOnset(:,:,20)')
thetaPeriods = [Stc{'theta',xyz.sampleRate}];

thetaPoints = [];
for t = thetaPeriods.data',
    thetaPoints = [thetaPoints; [t(1):embeddingWindow:t(2)]'];
end

    
thetaLFPembedding = flfp.segs(thetaPoints-30,60,0);
thetaLFPembedding = reshape(permute(thetaLFPembedding,[2,1,3]),[],embeddingWindow*numel(lfpChannels));

[U,S,V] = svd(thetaLFPembedding,'econ');

wts = (1:embeddingWindow)./xyz.sampleRate;
hfig = figure;
hfig.Units = 'centimeters';
hfig.Position(3:4) = [30,4];
hfig.PaperPositionMode = 'auto';
for i = 1:50,
    subplot(5,10,i);imagesc(wts,1:numel(lfpChannels),reshape(V(:,i),[],numel(lfpChannels))'),
    caxis([-0.08,0.08]);
    axis xy
end

flfp.data

% EMBEDDED lfp feature
tfs = flfp.segs([],embeddingWindow);
tfs = circshift(tfs,embeddingWindow/2,2);
tfs = MTADxyz('data',reshape(permute(tfs,[2,1,3]),size(tfs,2),[]),'sampleRate',xyz.sampleRate);
tfs.data(isnan(tfs.data(:)))=0;

% COMPUTE timeseries score for first 10 eigenvectors
fetT = MTADxyz('data',tfs.data*V(:,1),'sampleRate',xyz.sampleRate);
for i = 1:25,fetT.data(:,i) = tfs.data*V(:,i);end

fetTrms = fetT.segs([],embeddingWindow);
fetTrms = sq(rms(fetTrms,1));
fetTrms = circshift(fetTrms,embeddingWindow/2,1);
fetTrms = MTADxyz('data',fetTrms,'sampleRate',xyz.sampleRate);
fetTrms.data(isnan(fetTrms.data(:)))=0;
fetTrms.filter('ButFilter',3,2,'low');


states = {'walk','turn','pause','rear'};
sclr = 'bgcr';
ts = [1:size(xyz,1)]'./xyz.sampleRate;
sp = [];

figure,
sp(end+1) = subplot2(10,4,1:8,1:4); hold on
plot(ts,fetTrms(:,1:2:6))

sp(end+1)=subplot2(10,4,[9,10],[1:4]);
plotSTC(Stc,1,'text',states,sclr,[],false);
linkaxes(sp,'x');


figure, 
for v = 1:25,
    subplot(5,5,v);hold on;
edx = linspace(0,80,100);
ind = Stc{'w'};
h = bar(edx,histc(fetTrms(ind,v),edx),'histc');
h.FaceColor = 'r';h.EdgeColor = 'r';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;
hold on
ind = Stc{'r'};
h = bar(edx,histc(fetTrms(ind,v),edx),'histc');
h.FaceColor = 'c';h.EdgeColor = 'c';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;
end

figure
for v = 2:2:20
    for k = 2:2:20,
        if v>k
            ind = Stc{'w'};
        elseif k>v
            ind = Stc{'p'};
        else
            continue
        end
        subplot2(10,10,v/2,k/2);hold on;
        hist2([fetTrms(ind,v),fetTrms(ind,k)],edx/v*2,edx/k*2);axis tight
    end
end

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');
fvxy = fxyz.vel([],[1,2]);
fvxy.data(fvxy.data<1e-3) = 1e-3;
fvxy.data = log10(fvxy.data);

fetTrms.data(~nniz(fetTrms.data)) = 1;
fetTrms.data = log10(fetTrms.data);

edv = linspace(-2,2,50);
edt = linspace(.5,1.6,50);
figure
i = 1;
for s = 'arwnpms';
    subplot(3,3,i); i = i+1;
    ind = Stc{s};
    hist2([fetTrms(ind,3),fvxy(ind,1)],edt,edv);
end