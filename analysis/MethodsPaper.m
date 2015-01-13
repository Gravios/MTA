
% figure 1


% subplot f


Trial = MTASession('Ed10-20140812');

xyz = Trial.load('xyz');




lfp = Trial.lfp.copy;
lfp.create(Trial,66);
lfp.resample(xyz);



specParms = struct('nFFT',2^11,...
                   'Fs',lfp.sampleRate,...
                   'WinLength',2^10,...
                   'nOverlap',2^10*.875,...
                   'FreqRange',[1,15]);


lfp.data = cat(2,lfp.data,xyz(:,1,3));
lfp.data = cat(2,lfp.data,xyz(:,3,3));

xyz.filter(gtwin(.5,xyz.sampleRate));
ang = Trial.ang.copy;
ang.create(Trial,xyz);


rhm = fet_rhm(Trial);
lfp.data = cat(2,lfp.data,rhm);

[ys,fs,ts,phi,fstat] = fet_spec(Trial,lfp,'mtchglong','overwrite',true);


figure,imagesc(ts,fs,ys(:,:,1,2)'),axis xy


