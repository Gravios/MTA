function OutArgs = rearss(Session,varargin)
[spkShift,numIterations,binSize,halfBins,normalization] = DefaultArgs(varargin,{0.25,10,100,50,'hz'});

%% Seed RandStream with system clock
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));


Session = MTASession(Session,{'CluRes','ang'});

%% Number of Units
numClu = size(Session.map,1);

%% Rearing onset/offset periods in xyzSampleRate
rper = Session.Bhv.getState('rear').state;
rdur = diff(rper,1,2);

fet = GetSegs(Session.ang(:,3,4,2),round(rper-3.*Session.xyzSampleRate),round(6*Session.xyzSampleRate),0);
fet = reshape(fet,round(6*Session.xyzSampleRate),size(rper,1),2);
gri = find(max(fet(1:300,:,1))<0&rdur'>1.5);
griu = find(max(fet(450:end,:,2))<0&rdur'>1.5);

spkShift = round(spkShift*Session.lfpSampleRate);

per_on  = rper(gri,1);
per_off = rper(griu,2);
poni  = round((per_on-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
poffi = round((per_off-1)./Session.xyzSampleRate.*Session.lfpSampleRate);

%% CCG of non-rearing periods and units

brccg = zeros(numIterations,halfBins*2+1,numClu,2,1);
for i=1:numIterations
    Session.res = Session.res+randi([-spkShift,spkShift],size(Session.res,1),1);
    Session.res(Session.res<1) = 1;
    for unit=1:numClu,
        [tccg tbin ] = Trains2CCG({poni,poffi,Session.res(Session.clu==unit)},{1,2,Session.clu(Session.clu==unit)},binSize,halfBins,Session.lfpSampleRate,normalization);
        brccg(i,:,unit,:,1) = sq(tccg(:,1:2,3:end));
    end
end

OutArgs = {brccg};
save([Session.spath.analysis Session.filebase '.rearss.' num2str(randi([100000,999999],1)) '.mat'],'OutArgs','-v7.3');
















