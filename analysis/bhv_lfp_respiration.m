Trial = MTATrial('jg05-20120317');
Trial.load('xyz');
periods = [268200,269700;294200,329900];

Trial.load('nq');
Trial.lfp.load(Trial,[71:3:84]);

tbp = ButFilter(Trial.lfp.data(:,:),3,[2,5]./(Trial.lfp.sampleRate/2),'bandpass');
tbp_hilbert = Shilbert(tbp);
tbp_phase = phase(tbp_hilbert);
tbp_phase = MTADlfp([],[],tbp_phase,Trial.lfp.sampleRate,Trial.lfp.sync.copy,Trial.lfp.origin);
tbp_phase.resample(Trial.xyz);


resp = Filter0(gausswin(61)./sum(gausswin(61)),Trial.xyz(periods,4,3));
rphs = tbp_phase(periods,:);

rpks = LocalMinima(-resp,40);

%plot(resp)
%Lines(rpks,[],'k')

figure,
histcirc(rphs(rpks(end-50:end)),60)
