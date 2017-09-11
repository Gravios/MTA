function [out,xedges,yedges,labels] = bhv_rhm_rrm(Trial,varargin)
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

out    = {};
xedges = {};
yedges = {};
labels = {};

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
pitch = fet_HB_pitch(Trial);
map_to_reference_session(pitch,Trial,referenceTrial);

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

% DETECT peaks and troughs of rhythmic head motion
[rhmTroughs,rhmTroughsVals] = LocalMinima(rhm.data,5,rhmThreshold);
rhmTroughs(~WithinRanges(rhmTroughs,sper.data)) = [];
[rhmPeaks,rhmPeaksVals] = LocalMinima(-rhm.data,5,rhmThreshold);
rhmPeaks(~WithinRanges(rhmPeaks,sper.data)) = [];
rhmPeaksVals = -rhmPeaksVals;

nbins = 30;
f = 2;
if (f == 1), dedge = linspace(0,0.5,nbins); elseif (f == 2), dedge = linspace(-.8,0.8,nbins);end


pedge = linspace(-pi,pi,nbins);

i = 1;
[out{i},xedges{i},yedges{i}] = histcounts2(circ_dist(prtm(rhmTroughs),prhm(rhmTroughs)),pitch(rhmTroughs,f),pedge,dedge);
labels{i} = 'head pitch versus rhm-rrm trough phase diff';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(circ_dist(prtm(rhmPeaks),prhm(rhmPeaks)),pitch(rhmPeaks,f),pedge,dedge);
labels{i} = 'head pitch versus rhm-rrm peak phase diff';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(circ_dist(prpm(rhmTroughs),prhm(rhmTroughs)),pitch(rhmTroughs,f),pedge,dedge);
labels{i} = 'head pitch versus rhm-rpm trough phase diff';

i = i+1;
[out{i},xedges{i},yedges{i}] = histcounts2(circ_dist(prpm(rhmPeaks),prhm(rhmPeaks)),pitch(rhmPeaks,f),pedge,dedge);
labels{i} = 'head pitch versus rhm-rpm peak phase diff';


% $$$ figure
% $$$ for i = 1:4,
% $$$     subplot(2,2,i);
% $$$     imagesc(xedges{i},yedges{i},out{i}');
% $$$     axis('xy');
% $$$     title(labels{i});
% $$$ end
