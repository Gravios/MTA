


Trial = MTATrial.validate('jg05-20120317.cof.all');
xyz = Trial.load('xyz');
ang = create(Trial.ang.copy,Trial,xyz);
dang = ang.copy;
dang.data = ang(Trial.stc{'a'},'head_back','head_front',3)-...
            median(ang(Trial.stc{'a'},'head_back','head_front',3));


pstr = struct('nFFT',2^11,'Fs',dang.sampleRate,...
              'WinLength',2^10,'nOverlap',2^10*.5,...
              'FreqRange',[0.1,60]);
[ys,fs,ts] = fet_spec(Trial,dang,'mtchglong',false,[],pstr);
pedges = linspace(0.0001,0.1,100);
N = histc(sqrt(ys(nniz(ys),:)),pedges,1);
figure,imagesc(fs,pedges,N),axis xy;
xlabel('Frequency (Hz)')
ylabel('mm/Hz')
title({'PSD of distance between two markers','of a ridgid body structure'})


