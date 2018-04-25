
Trial = MTATrial.validate('jg05-20120317.cof.all');
Trial = MTATrial.validate('jg05-20120310.cof.all');

xyz = preproc_xyz(Trial,'trb');

spkw = Trial.spk.copy();
spkw.load_spk(Trial,1,'theta',[],'deburst');


lfp = Trial.load('lfp',78);
lfp.resample(xyz);
phz = lfp.phase([6,12]);

pft = pfs_2d_theta(Trial,[],false,true,1);
mrt = pft.maxRate();

drz = compute_drz(Trial,pft.data.clu,pft);

% PLOT mean substracted waveform binned by drz
unit = 68;
for unit = pft.data.clu(1:end),
numBins = 3;
figure(39203209);
clf();
subplot(numBins+1,numBins,1);
plot(pft,unit,1,false,[],true);
title(num2str(unit))


res = round(spkw(unit).*xyz.sampleRate)+1;
if numel(res) < 60 | pft.data.si(unit)<0.2, continue, end
wav = spkw.spk(spkw.clu==unit,1:8,:);

phzBin    = linspace(-pi,pi,numBins+1);
phzBinCen = (phzBin(1:end-1)+phzBin(2:end))/2;
phzBinInd = discretize(phz(res),phzBin);

drzBin    = linspace(-1,1,numBins+1);
drzBinCen = (drzBin(1:end-1)+drzBin(2:end))/2;
drzBinInd = discretize(drz(res,unit),drzBin);

meanWaveForm = sq(mean(wav));


subplot(numBins+1,numBins,2);
imagesc(meanWaveForm);
for binIndp = 1:numel(drzBinCen),
    for binIndd = 1:numel(drzBinCen),    
        try
        subplot2(numBins+1,numBins,binIndp+1,binIndd);
        %subplot(1,numBins+2,binInd+2);
        imagesc(sq(mean(wav(drzBinInd==binIndd&phzBinInd==binIndp,:,:)))-meanWaveForm);
        %imagesc(sq(mean(wav(phzBinInd==binInd,:,:)))-meanWaveForm); 
        caxis([-300,300]);
        end
    end
end

waitforbuttonpress();
end
 
