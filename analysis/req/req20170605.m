% req20170605
%
%  Status: inactive
%  Type: Analysis 
%  Final_Forms: fet_rhmPCA
%  Description: comparison 
%  Bugs: NA
%

Trial = MTATrial('jg05-20120317');

Trial = MTATrial('Ed03-20140624');

Trial = MTATrial('Ed01-20140707');

% fet_rhmPCA
varargin = {};

% IF nasal cavity pressure sensor exists
ncp = fet_ncp(Trial);
rhm = fet_rhm(Trial);

% @ line 42
ncp.resample(xyz);

% @ line 56

sp = [];
figure,
sp(1) = subplot(311);
imagesc(V(:,:,1)');
sp(2) = subplot(312);
imagesc(pv(:,:,1)');
sp(3) = subplot(313);hold('on');
plot(fet.data);
plot(nunity(ncp.data))
linkaxes(sp,'x');

fet.resample(xyz);

fet.data = nunity([fet.data,ncp.data]);

parspec = struct('nFFT',2^8,...
                 'Fs',  xyz.sampleRate,...
                 'WinLength',2^7,...
                 'nOverlap',2^7*.875,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[.1,20]);

data = zeros(fet.size);
data(nniz(fet.data),:) = fet.data(nniz(fet.data),:);
mode = 'mtchglong';
[ys,fs,ts] = spec(str2func(mode),data,parspec);
ts = ts+(parspec.WinLength/2)/fet.sampleRate;
ssr = 1/diff(ts(1:2));
pad = round([ts(1),mod(fet.size(1)-round(parspec.WinLength/2),parspec.WinLength)/fet.sampleRate].*ssr)-[1,0];
szy = size(ys);
yhm = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);
ts = cat(2,([1:pad(1)]./ssr),ts',([pad(1)+size(ts,1)]+[1:pad(2)])./ssr)';

%fet.get_timestamps();
fts = [1:fet.size(1)]./fet.sampleRate;

sp = [];
figure,
sp(1) = subplot(411);
imagesc(ts,fs,log10(yhm(:,:,1,1))')
caxis([-3,-1]),colormap jet
axis xy
sp(2) = subplot(412);
imagesc(ts,fs,log10(yhm(:,:,2,2))')
caxis([-3,0]),colormap jet
axis xy
sp(3) = subplot(413);
imagesc(ts,fs,rhm(:,:,1,2)')
caxis([0.5,1]),colormap jet
axis xy
sp(4) = subplot(414);
plot(fts,nunity(fet.data));
linkaxes(sp,'x');


