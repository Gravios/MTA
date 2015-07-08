%Trial = MTATrial('jg05-20120310');
Trial = MTATrial('jg05-20120309');
Trial.stc.updateMode('auto_wbhr');
Trial.stc.load;
%'jg05-20120317'
%periods = [268200,269700;294200,329900];
freqRng = [6,12];
chans = [65:3:96];
Trial.load('nq');

xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.1,xyz.sampleRate));

ang = Trial.ang.copy;
ang.create(Trial,xyz);
ahp = ButFilter(ang(:,4,5,3),3,freqRng./(Trial.ang.sampleRate/2),'bandpass');
rhm = fet_rhm(Trial);
ahp = ButFilter(rhm,3,freqRng./(Trial.ang.sampleRate/2),'bandpass');
figure,plot(ahp)


lfp = Trial.lfp.copy;
lfp.load(Trial,chans);
lfp.resample(xyz);
tbp_phase = lfp.phase(freqRng);


sts = 'wgl';
sts = {'walk&theta','hwalk&theta','lwalk&theta'};
nsts = numel(sts);
nchan = numel(chans);

figure(2323923),clf
for s = 1:nsts
for c = 1:nchan
sper = Trial.stc{sts{s},xyz.sampleRate};
sper.cast('TimeSeries');
srpks = LocalMinima(-ahp(sper.data==1),10,-.15);

orpks = LocalMinima(-ahp(sper.data==1),10,-5);
srpks = srpks(~ismember(srpks,orpks));
%srpks = srpks(5:6:end);
%srpks = srpks(randi(numel(srpks),round(numel(srpks)/5),1000));
%srpks = sort(srpks);
%srpks = srpks([0;diff(srpks)]>xyz.sampleRate);
stbp_phase = tbp_phase(sper.data==1,:);

% $$$ for iter = 1:100
% $$$ [N(:,iter),b] =histcirc(stbp_phase(srpks(:,iter),c),20);
% $$$ end


subplot2(nchan,nsts,c,s);
%bar(b,mean(N,2),'histc'),axis tight
%bar(b,std(N,[],2),'histc'),axis tight
histcirc(stbp_phase(srpks(:,1),c),20);
% $$$ boundedline(b,mean(N,2),std(N,[],2)*2);
% $$$ xlim([0,710])
%bar(b,N,'histc'),axis tight

if c==nchan, xlabel(['theta(' num2str(freqRng(1)) '-' num2str(freqRng(2)) 'Hz) phase']),end
if c==1,     title (sper.label),end
if s==1, ylabel(['chan:' num2str(chans(c))]),end
end
end


% $$$ dsr = diff(srpks)*diff(srpks)';
% $$$ [ssrs,ind] = sort(dsr(:));
% $$$ ssr = diff(ssrs);
% $$$ x = zeros([prod(size(dsr)),1]);
% $$$ x(ind) = [false;ssr~=0];
% $$$ x = reshape(x,size(dsr));
% $$$ srpks(sum(x)==0) = [];

% $$$ rpks = LocalMinima(ahp,10,-.5);
% $$$ figure,
% $$$ histcirc(tbp_phase(rpks,2),15)
% $$$ %nni = nniz([tbp_phase.data,ahp]);
% $$$ figure,
% $$$ hist2([tbp_phase(rpks),log10(abs(ahp(rpks)))],30,30);


% $$$ sper = Trial.stc{'g'}.copy;
% $$$ sper.cast('TimeSeries');
% $$$ srpks = LocalMinima(ahp(sper.data==1),10,-.2);
% $$$ stbp_phase = tbp_phase(sper.data==1,:);
% $$$ figure,
% $$$ histcirc(stbp_phase(srpks,2),20)



%nni = nniz([tbp_phase.data,ahp]);
%figure,hist2([tbp_phase(rpks),log10(abs(ahp(rpks)))],30,30)
