function [out,xedges,yedges,labels] = bhv_rhm_rrm_ncp(Trial,varargin)
%function [out,labels] = bhv_rhm_rrm(Trial,varargin)
%
%    varargin: (Name, Type, Default)
%
%        state:             (char)       'w+p+n'
%        stcMode:           (char)       'msnn_ppsvd'
%        freqRange:         (numeric)    [4,12],
%        rhmThreshold:      (numeric)    -0.2
%        referenceTrial:    (char)       'jg05-20120317.cof.all'
%

Trial = MTATrial.validate(Trial);

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('state',                  'w+p+n',                                              ...
                 'stcMode',                'msnn_ppsvd',                                         ...
                 'freqRange',              [4,12],                                               ...
                 'rhmThreshold',           -0.2,                                                 ...
                 'referenceTrial',         'jg05-20120317.cof.all'                               ...
);
[state,stcMode,freqRange,rhmThreshold,referenceTrial] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% LOAD behavioral state collection
stc = Trial.load('stc',stcMode);
% SET the state for subset selection
sper = stc{state};

% LOAD marker position data
xyz = preproc_xyz(Trial);

% CREATE Virtual marker, on axis orthogonal to the transverse plane of the rat's head
nz  = -cross(xyz(:,'head_back',:)-xyz(:,'hcom',:),xyz(:,'head_left',:)-xyz(:,'hcom',:));
nz  = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3)));
nm  = nz.*20+xyz(:,'hcom',:);
xyz.addMarker('htx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

% CREATE Virtual marker, a low pass filtered copy of the head's center of mass
xyz.addMarker('flhcom',[.7,1,.7],...
              {{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'hcom',:),3,freqRange(1)./(xyz.sampleRate/2),'low'));

% LOAD rhythmic head motion 
rhm = fet_rhm(Trial);

% COMPUTE intermarker angles
ang = create(MTADang,Trial,xyz);

% CREATE head pitch feature
[pitch,plabel] = fet_HB_pitch(Trial);
map_to_reference_session(pitch,Trial,referenceTrial);

% LOAD nasal cavity pressure signal
ncp = fet_ncp(Trial);

% COMPUTE longitudinal rhythmic rotational motion of the head
rtm = rhm.copy();
rtm.data = circ_dist(ang(:,'flhcom','htx',2),ang(:,'flhcom','head_front',2));
rtm.filter('ButFilter',3,freqRange,'bandpass');

% COMPUTE rhythmic pitch motion of the head
rpm = rhm.copy();
rpm.data = pitch(:,2);
rpm.filter('ButFilter',3,freqRange,'bandpass');

% COMPUTE phase of each timeseries
prhm = phase(rhm,freqRange);
prtm = phase(rtm,freqRange);
prpm = phase(rpm,freqRange);
pncp = phase(ncp,freqRange);

% DETECT peaks and troughs of respiration in the ncp signal
[rhmTroughs,rhmTroughsVals] = LocalMinima(rhm.data,5,rhmThreshold);
rhmTroughs(~WithinRanges(rhmTroughs,sper.data)) = [];
[rhmPeaks,rhmPeaksVals] = LocalMinima(-rhm.data,5,rhmThreshold);
rhmPeaks(~WithinRanges(rhmPeaks,sper.data)) = [];
rhmPeaksVals = -rhmPeaksVals;

% DETECT peaks and troughs of rhythmic head motion
ncpThresh = -2000;
[ncpTroughs,ncpTroughsVals] = LocalMinima(ncp.data,5,ncpThresh);
ncpTroughs(~WithinRanges(ncpTroughs,sper.data)) = [];
[ncpPeaks,ncpPeaksVals] = LocalMinima(-ncp.data,5,ncpThresh);
ncpPeaks(~WithinRanges(ncpPeaks,sper.data)) = [];
nncpPeaksVals = -ncpPeaksVals;

rhmPeaks(log10(yrhm(rhmPeaks,30))<-5) = [];
rhmTroughs(log10(yrhm(rhmTroughs,30))<-5) = [];
ncpPeaks(log10(yrhm(ncpPeaks,30))<-5) = [];
ncpTroughs(log10(yrhm(ncpTroughs,30))<-5) = [];

nbins = 30;
f = 2;
if (f == 1), dedge = linspace(0,0.5,nbins); elseif (f == 2), dedge = linspace(-1.2,0.8,nbins);end

pedge = linspace(-pi,pi,nbins);


out    = {};
xedges = {};
yedges = {};
labels = {};

% rhm peaks and troughs
% $$$ i = 1;
% $$$ [out{i},xedges{i},yedges{i}] = histcounts2(circ_dist(prtm(rhmTroughs),prhm(rhmTroughs)),pitch(rhmTroughs,f),pedge,dedge);
% $$$ labels{i} = [plabel{f},' versus rhm-rrm phase diff at rhm trough'];
% $$$ 
% $$$ i = i+1;
% $$$ [out{i},xedges{i},yedges{i}] = histcounts2(circ_dist(prtm(rhmPeaks),prhm(rhmPeaks)),pitch(rhmPeaks,f),pedge,dedge);
% $$$ labels{i} = 'head pitch versus rhm-rrm phase diff at rhm peak';
% $$$ 
% $$$ i = i+1;
% $$$ [out{i},xedges{i},yedges{i}] = histcounts2(circ_dist(prpm(rhmTroughs),prhm(rhmTroughs)),pitch(rhmTroughs,f),pedge,dedge);
% $$$ labels{i} = 'head pitch versus rhm-rpm phase diff at rhm trough';
% $$$ 
% $$$ i = i+1;
% $$$ [out{i},xedges{i},yedges{i}] = histcounts2(circ_dist(prpm(rhmPeaks),prhm(rhmPeaks)),pitch(rhmPeaks,f),pedge,dedge);
% $$$ labels{i} = 'head pitch versus rhm-rpm phase diff at rhm peak';
% $$$ 

% ncp peaks and troughs
% $$$ i = i+1;
% $$$ [out{i},xedges{i},yedges{i}] = histcounts2(circ_dist(prtm(ncpTroughs),pncp(ncpTroughs)),pitch(ncpTroughs,f),pedge,dedge);
% $$$ labels{i} = 'head pitch versus rhm-rrm phase diff at ncp trough ';
% $$$ 
% $$$ i = i+1;
% $$$ [out{i},xedges{i},yedges{i}] = histcounts2(circ_dist(prtm(ncpPeaks),pncp(ncpPeaks)),pitch(ncpPeaks,f),pedge,dedge);
% $$$ labels{i} = 'head pitch versus rhm-rrm phase diff at ncp peak';
% $$$ 
% $$$ i = i+1;
% $$$ [out{i},xedges{i},yedges{i}] = histcounts2(circ_dist(prpm(ncpTroughs),pncp(ncpTroughs)),pitch(ncpTroughs,f),pedge,dedge);
% $$$ labels{i} = 'head pitch versus rhm-rpm phase diff at ncp trough';
% $$$ 
% $$$ i = i+1;
% $$$ [out{i},xedges{i},yedges{i}] = histcounts2(circ_dist(prpm(ncpPeaks),pncp(ncpPeaks)),pitch(ncpPeaks,f),pedge,dedge);
% $$$ labels{i} = 'head pitch versus rhm-rpm phase diff at ncp peak';


i = 1;
[out{i},xedges{i},yedges{i}] = histcounts2(prtm(rhmTroughs),pitch(rhmTroughs,f),pedge,dedge);
labels{i} = 'rrm phase at rhm trough';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(prtm(rhmPeaks),pitch(rhmPeaks,f),pedge,dedge);
labels{i} = 'rrm phase at rhm peak';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(prhm(rhmTroughs),pitch(rhmTroughs,f),pedge,dedge);
labels{i} = 'rhm phase at rhm trough ';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(prhm(rhmPeaks),pitch(rhmPeaks,f),pedge,dedge);
labels{i} = 'rhm phase at rhm peak';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(prpm(rhmTroughs),pitch(rhmTroughs,f),pedge,dedge);
labels{i} = 'pitch phase at rhm trough';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(prpm(rhmPeaks),pitch(rhmPeaks,f),pedge,dedge);
labels{i} = 'pitch phase at rhm peak';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(pncp(rhmTroughs),pitch(rhmTroughs,f),pedge,dedge);
labels{i} = 'ncp phase at rhm trough';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(pncp(rhmPeaks),pitch(rhmPeaks,f),pedge,dedge);
labels{i} = 'ncp phase at rhm peak';




% ncp peaks and troughs
i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(prtm(ncpTroughs),pitch(ncpTroughs,f),pedge,dedge);
labels{i} = 'rrm phase at ncp trough ';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(prtm(ncpPeaks),pitch(ncpPeaks,f),pedge,dedge);
labels{i} = 'rrm phase at ncp peak';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(prhm(ncpTroughs),pitch(ncpTroughs,f),pedge,dedge);
labels{i} = 'rhm phase at ncp trough ';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(prhm(ncpPeaks),pitch(ncpPeaks,f),pedge,dedge);
labels{i} = 'rhm phase at ncp peak';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(prpm(ncpTroughs),pitch(ncpTroughs,f),pedge,dedge);
labels{i} = 'rpm phase at ncp trough';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(prpm(ncpPeaks),pitch(ncpPeaks,f),pedge,dedge);
labels{i} = 'rpm phase at ncp peak';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(pncp(ncpTroughs),pitch(ncpTroughs,f),pedge,dedge);
labels{i} = 'ncp phase at ncp trough';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(pncp(ncpPeaks),pitch(ncpPeaks,f),pedge,dedge);
labels{i} = 'ncp phase at ncp peak';




figure
for i = 1:16,
    subplot(4,4,i);
    imagesc(xedges{i},yedges{i},out{i}');
    axis('xy');
    title(labels{i});
    caxis([0,100])
end


[yrhm,fs,ts] = fet_rhm(Trial,[],'mtchglong');
yrhm.resample(ncp);
yrhm.data(yrhm.data<1e-10) = 1e-10;

[yncp,fs,ts] = fet_ncp(Trial,[],'mtchglong');
yncp.resample(ncp);
yncp.data(yncp.data<1e-10) = 1e-10;


np = yncp(ncpPeaks,30);
rp = yrhm(ncpPeaks,30);
rpgi = nniz(rp)&nniz(np);
rp = rp(rpgi);
np = np(rpgi);
pp = pitch(ncpPeaks,f);
pp = pp(rpgi);


figure,
subplot(221);
histogram2(log10(rp),pp,linspace(-8,-2,30),dedge,'DisplayStyle','tile');
subplot(222);
histogram2(log10(np),pp,linspace(2,6,30),dedge,'DisplayStyle','tile');
subplot(223);
histogram2(log10(yrhm(sper,30)),pitch(sper,2),linspace(-8,-2,30),dedge,'DisplayStyle','tile');
subplot(224);
histogram2(log10(yncp(sper,30)),pitch(sper,2),linspace(2,6,30),dedge,'DisplayStyle','tile');


sind{i} = sind{i}.^2/SmoothingWeights(i)^2/2;
Smoother = exp(sum(-cat(2,),ndims+1));
Smoother = Smoother./sum(Smoother(:));
SOcc = convn(Occupancy,Smoother,'same');
