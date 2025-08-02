%Session = ...

% LOAD parameter file for the session
Par = LoadPar(fullfile(Session.spath, [Session.name '.xml']));

% SELECT the first channel of each Anatomical Group (AnatGrps)
channels = arrayfun(@(x) x.Channels(1),Par.AnatGrps)+1;

% LOAD lfp signal into an empty MTADlfp object
lfp = MTADlfp('data',LoadBinary(fullfile(Session.spath, [Session.name '.lfp']),...
                                channels,...
                                Par.nChannels)','sampleRate',Par.lfpSampleRate);
lfp.data(lfp.data==0)=1;

lfp.resample(30);

tpow=lfp.copy;
tpow.filter('ButFilter',3,[6,12],'bandpass');
dpow=lfp.copy;
dpow.filter('ButFilter',3,[1, 5],'bandpass');

tpow.data = sq(sum(tpow.segs(1:tpow.size(1),30).^2));
dpow.data = sq(sum(dpow.segs(1:dpow.size(1),30).^2));

tdr = lfp.copy;
tdr.data = log10(tpow.data./dpow.data);

% LOAD clusters and res
[Res, Clu, Map] = LoadCluRes(fullfile(Session.spath, Session.name));
Res = ceil(Res/Session.sampleRate*lfp.sampleRate);

% ALLOCATE output array
units = 1:Map(end,1);
Data.data = zeros([size(lfp,1),numel(units)]);

% CREATE convolution window 
twin = 1;
swin = round(twin*lfp.sampleRate);
gwin = ones([swin,1]);

% ACCUMULATE and convolve activity of each unit
for unit = units(:)'
    myRes = Res(Clu==unit);
    myRes(myRes==0) = 1;
    myRes(myRes>size(lfp,1)) = size(lfp,1);
    Data.data(:,unit==units) = conv(accumarray(myRes,1,[size(lfp,1),1])./twin,gwin,'same');
end

