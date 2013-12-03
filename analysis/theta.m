function theta(Session)

Session = MTASession('jg05-20120310');


[Res,Clu,Map] = LoadCluRes([Session.spath.nlx Session.name]);
Res = round(Res*Session.lfpSampleRate/Session.sampleRate);
[Res,Ind] = SelectPeriods(Res,[Session.syncPeriods(1),Session.syncPeriods(end)],'d',1);
Clu = Clu(Ind);

%% Calculate Lfp Spectrum Across Linear Shank
channels = 65:96;
lfp = Session.loadlfp(channels);
wlfp = WhitenSignal(lfp,[],1);
nfft = [2^11,2^9,2^7];
wind = [2^10,2^8,2^6];
fqrange = [1,30;30,60;60,120];
y = {};
f = {};
t = {};
for i = 1:3,
    for chan = 1:length(channels),
        [y{i}(:,:,chan),f{i},t{i}] = mtchglong(wlfp(:,chan),nfft(i),Session.lfpSampleRate,wind(i),wind(i)*.75,3,'linear',[],fqrange(i,:));
    end
end

channels = 65:96;
ytheta = [];
ftheta = [];
ttheta = [];
for chan = 1:length(channels),
    [ytheta(:,:,chan),ftheta,ttheta] = mtchglong(wlfp(:,chan),2^9,Session.lfpSampleRate,2^8,2^8*.75,3,'linear',[],[5 13]);
end

spec = load ('/data/homes/gravio/data/analysis/jg05-20120310/jg05-20120310.spectrums.mat');
ttheta = spec.ttheta+diff(spec.ttheta(1:2));


%% Calculate Instantaneous Firing Rate 

%% check to see how SelectPeriods really works
[Res,Clu,Map] = LoadCluRes([Session.spath.nlx Session.name]);
%[Res,Ind] = SelectPeriods(Res,[Session.syncPeriods(1),Session.syncPeriods(end)]./Session.lfpSampleRate.*Session.sampleRate,'d',1);
%Clu = Clu(Ind);

twin = 0.05; %ms
gwin =gausswin(round(twin*Session.sampleRate));
spks = zeros(max(Res),1);
ufr = [];
for unit = 1:size(Map,1),
    uRes = Res(Clu==unit);
    uClu = Clu(Clu==unit);
    spks(:)=0;
    spks(uRes) = 1;
    hufr = conv(spks,gwin);
    ufr(:,unit) =resample(hufr(round(Session.syncPeriods(1)/Session.lfpSampleRate*Session.sampleRate):end),Session.lfpSampleRate,Session.sampleRate);
end

for unit = 1:size(Map,1),
plot(ufr(ufr(:,unit)~=0,unit),frthpow(ufr(:,unit)~=0),'.');
pause(.9)
end

%% Calculate Instantaneous Firing Rate in 51.2 ms bins
binsz = 64;
unitfr = zeros(round(size(lfp,1)/binsz),size(Map,1));
for unit = 1:size(Map,1),
    subs = round(Res(Clu==unit)./binsz);
    vals = ones(size(subs,1),1);
    unitfr(1:max(subs),unit) = accumarray(subs,vals)./(binsz./Session.lfpSampleRate);
end
tu = (1:size(unitfr,1)).*(binsz./Session.lfpSampleRate);
tu = tu-tu(1)/2;
tu = tu./Session.lfpSampleRate;





for unit = 1:size(Map,1);
funitfr = Filter0(gausswin(3),unitfr(:,unit));
plot(log(funitfr(funitfr~=0)),frthpow(funitfr~=0),'.');
pause(0.9)
end



figure
sp1 = subplot2(3,1,1,1);
imagesc(t{3},f{3},log10(y{3}(:,:,8))'),axis xy
sp2 = subplot2(3,1,2,1);
imagesc(t{2},f{2},log10(y{2}(:,:,8))'),axis xy
sp3 = subplot2(3,1,3,1);
imagesc(t{1},f{1},log10(y{1}(:,:,8))'),axis xy
linkaxes([sp1,sp2,sp3],'x')