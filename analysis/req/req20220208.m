
configure_default_args();
MjgER2016_load_data();

Co = {};
sampleRate = 250;
for tind = 1:numel(Trials)
Trial = Trials{tind};
rhm = fet_rhm(Trial,sampleRate);
Trial.lfp.filename = [Trial.name,'.lfp'];    
tlfp = Trial.load('lfp',sessionList(tind).subject.channelGroup.theta);
tlfp.resample(rhm);
x = copy(tlfp);
x.data = cat(2,x.data,rhm.data);
cstates = {'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta','lpause&theta'};
specArgsTheta = struct('nFFT',2^9,...
                       'Fs',  x.sampleRate,...
                       'WinLength',2^8,...
                       'nOverlap',2^8*0.5,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',[1,20]);
[Co{tind},Po{tind},fs] = Comodugram(Trial, x, cstates, 'mtcsglong', specArgsTheta);
end


tind = 6;
figure,
for sts = 1:6,
    subplot2(2,6,1,sts);
    imagesc(fs,fs,Co{tind}(:,:,1,2,sts));
    axis('xy');
    colormap('jet');
    colorbar();
    caxis([-1,1].*max(abs(caxis())));
    subplot2(2,6,2,sts);
    imagesc(fs,fs,Po{tind}(:,:,1,2,sts));
    axis('xy');
    colormap('jet');
    caxis([0,0.01]);    
end



