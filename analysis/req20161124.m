%% Spike sorting inquiry

%  from preclustered units 
%    - waveform realignment
%    - template matching 
%

Session = MTATrial.validate('Ed10-20140820.cof.all');
Par = LoadPar(fullfile(Session.spath, [Session.name,'.xml']));




%for i = 1:numel(Par.SpkGrps)
i = 1
Clu = LoadClu(fullfile(Session.spath, [Session.name,'.clu.',num2str(i)]));
cind = ~ismember(Clu,[0,1]);
Res = LoadRes(fullfile(Session.spath, [Session.name,'.clu.',num2str(i)]));
Spk = permute(LoadSpk(fullfile(Session.spath, [Session.name,'.spk.1']),numel(Par.SpkGrps(i).Channels),Par.SpkGrps(i).nSamples),...
              [3,1,2]);

Clu = Clu(cind);
Res = Res(cind);
Spk = Spk(cind,:,:);

clear('cind');


% plot the mean waveform of a cluster across the 8 channels of a buz32 shank
figure,plot(bsxfun(@plus,sq(mean(Spk(Clu==20,:,:)))',-[0:7].*2000))

% plot diff of 100 spike waveforms for 5th channel
figure,plot(diff(sq(Spk(find(Clu==20,100,'first'),5,:))'))


meanWaveForm = sq(mean(Spk(Clu==20,:,:)))';

figure,plot(xcorr(meanWaveForm(:,1),sq(Spk(find(Clu==20,1,'first'), 1,:))'))

unit = 20;
cind = find(Clu==unit);
chanResShift = nan([numel(cind),numel(Par.SpkGrps(i).Channels)]);
for c = 1:numel(cind),
    for chan = 1:numel(Par.SpkGrps(i).Channels),
        chanCorr = xcorr(meanWaveForm(:,chan),sq(Spk(cind(c),chan,:))');
        [~,chanResShift(c,chan)] = max(chanCorr);
    end
end
%Res(Clu==20)=Res(Clu==20)+chanResShift-52

newMeanWaveForm = 
for c = 1:numel(cind),
Spk(cind(c),:,52/2-)/numel(cind)

figure,hist(mean(chanResShift,2),100)
