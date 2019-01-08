% req20180629 ----------------------------------------------------
%  Status: retired
%  Type: Analysis
%  Author: Justin Graboski
%  Final_Forms: NA
%  Description: population phase precession as function of drz and behaviors
%  Bugs: NA

Trial = MTATrial.validate('jg05-20120312.cof.all');
channels = 64:96;
lfp = Trial.load('lfp',channels);


xyz = preproc_xyz(Trial,'trb');
lfp.resample(xyz);
phz = lfp.phase([6,12]);


defspec = struct('nFFT',2^7,'Fs',lfp.sampleRate,...
                 'WinLength',2^6,'nOverlap',2^6*.875,...
                 'FreqRange',[25,125]);

tlfp = lfp.copy();
tlfp.data = lfp.data(:,channels==76);
[ys,fs,ts] = fet_spec(Trial,tlfp,'mtchglong',true,[],defspec);


figure,
imagesc(ts,fs,log10(ys.data')),
axis('xy');
colormap(gca,'jet');
caxis([3,5])

ysus = ys.copy();
ysus.resample(xyz);

nysus = ysus.copy();
nysus.data = nysus.data>prctile(nysus.data(:),90);


figure,
imagesc(nysus.data');
axis('xy');
colormap(gca,'jet');
caxis([0,1])

xyz.filter('ButFilter',4,2.5,'low');
vxy = xyz.vel({'spine_lower','hcom'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);

binSizes = [25,25];

binsPhase = linspace(-pi,pi,binSizes(1));
binsVxy   = linspace(-1,2,binSizes(2));
%binsVxy   = linspace(0,40,binSizes(2));

tper = Trial.stc{'x&t',xyz.sampleRate};

indPhase = discretize(phz(tper,channels==72),binsPhase);
indVxy = discretize(vxy(tper,1),binsVxy);


gysus = nysus.copy();
gysus.data(:,:) = 0;
for f = 1:numel(fs),
gamats = ThreshCross(nysus(:,f),0.5,10);
gysus.data(gamats,f) = true;
end
 
gpow = gysus(tper,10);
%gpow = abs(log10(ysus(tper,20)));

nind = nniz(indPhase)&nniz(indVxy)&nniz(gpow);
%A = accumarray([indPhase(nind),indVxy(nind)],gpow(nind),[binSizes-1],@mean);
A = accumarray([indPhase(nind),indVxy(nind)],gpow(nind),[binSizes-1],@sum);


figure,imagesc(binsPhase,binsVxy,A'),axis('xy')
%figure,imagesc(binsPhase,binsVxy,bsxfun(@rdivide,A,sum(A))'),axis('xy')
colormap(gca,'jet')

