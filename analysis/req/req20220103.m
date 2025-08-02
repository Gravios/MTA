% req20220103
%     Tags: realtime state segmentation lm rad ratio thetarc
%     Status: Active
%     Type: Analysis
%     Author: Justin Graboski
%     Final_Forms: NA
%     Project: General
%     Description: theta return current feature used for realtime ephys state segmentation
%

% Notes:
%  quarter second time window does not capture the feature
%  of low <20Hz vs high >40 and <200 Hz as well
%  as the filtered means of the eigth of a second spectra.


%STARTHERE
% detect for other sessions
% there is a difference between high and low walk
% also do gamma detection across layers on the linear probe
% Use shift matching or some numerical depth to align clusters between sessions

%%%<<< Load Data ---------------------------------------------------------------
ThetaRC_load_data();
% - sampleRate
% - xyz
% - units
% ... etc

trialIndex = 18;

Trial   = Trials     { trialIndex };
units   = Units      { trialIndex };
meta    = sessionList( trialIndex );
pft = pfs_2d_theta(Trial,units);


txyz = preproc_xyz(Trial,'trb',sampleRate);
headYawCorrection = Trial.meta.correction.headYaw;
headCenterCorrection = Trial.meta.correction.headCenter;
% COMPUTE head basis
hvec = txyz(:,'nose',[1,2])-txyz(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(headYawCorrection),-sin(headYawCorrection);...
                  sin(headYawCorrection), cos(headYawCorrection)],...
                 [2,3],...
                 [1,2]);
hvfl = fet_href_HXY(Trial, sampleRate, [], 'trb', 2.4);
hafl = circshift(hvfl.data,1)-hvfl.data;
%bvfl = fet_bref_BXY(Trials{trlI}, sampleRate, [], 'trb', 5);
%bafl = circshift(bvfl(:,1),1)-bvfl(:,1);
hba = fet_hba(Trial,sampleRate); % Head to Body Angle
hav = fet_hbav(Trial,sampleRate);
pch = fet_HB_pitch(Trial,sampleRate);
pch.data = pch.data(:,3);



phz = load_theta_phase(Trial,Trial.lfp.sampleRate);
yphz = load_theta_phase(Trial,sampleRate);

%txyz = preproc_xyz(Trial,'trb',sampleRate);
fxyz = filter(txyz.copy(),'ButFilter',3,14,'low');
vxy = vel(filter(txyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
vxyz = vel(filter(txyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2,3]);
lvxy = copy(vxy);
lvxy.data(lvxy.data<=0.0001) = 0.0001;
lvxy.data = log10(lvxy.data);


lfpPyr = Trial.load('lfp', meta.subject.channelGroup.theta);

lfpPrc = Trial.load('lfp', [48,64]);
lfpPrc.data = diff(lfpPrc.data, 1, 2);
%lfp.resample(sampleRate);

cfp = Trial.load('lfp', 60);
%cfp.resample(sampleRate);

dfp = Trial.load('lfp', meta.subject.channelGroup.thetarc);
%rfp.resample(sampleRate);
rfp = copy(dfp)
rfp.data = diff(dfp.data, 1, 2);

ofp = Trial.load('lfp', [40,48,56]);

afp = Trial.load('lfp', [48,56,64]);

lfpRad = Trial.load('lfp', [73]);
lfpRrc = Trial.load('lfp', [72,74]);
lfpRrc.data = diff(lfpRrc.data, 1, 2);


lfpLm = Trial.load('lfp', [78]);
lfpLrc = Trial.load('lfp', [77,79]);
lfpLrc.data = diff(lfpLrc.data, 1, 2);

flfpRad = Trial.load('lfp', [76]);
flfpRad.filter('ButFilter',4,[6,12],'bandpass');
lfpRcsd = Trial.load('lfp', [75,77]);
lfpRcsd.filter('ButFilter',4,[6,12],'bandpass');
lfpRcsd.data = sum(lfpRcsd.data,2)-2*lfpRad.data;

lfpWhl = lfpRrc.copy();
lfpWhl.data = (lfpRrc.data+lfpLrc.data)/2;
phzWhl = lfpWhl.phase([5,12]);

tslfp = [1:size(lfpPyr,1)]./lfpPyr.sampleRate;
figure,hold('on');
plot(tslfp,lfpPyr.data+6000,'k');
plot(tslfp,lfpPrc.data+6000,'m');
plot(tslfp,lfpRad.data,'b');
plot(tslfp,lfpRrc.data,'r');
plot(tslfp,lfpLm.data-6000,'k');
plot(tslfp,lfpLrc.data-6000,'m');

plot(tslfp,(lfpRrc.data+lfpLrc.data)/2+10000,'g');

dlfpWhl = copy(lfpWhl);
dlfpWhl.resample(txyz);
dphzWhl = dlfpWhl.phase([5,12]);
dphzWhl.data = dphzWhl.data -pi/2-1;
dphzWhl.data(dphzWhl.data<0)  = dphzWhl.data(dphzWhl.data<0) +2*pi;

dlfpPyr = copy(lfpPyr);
dlfpPyr.resample(txyz);
dphzPyr = dlfpPyr.phase([5,12]);
dphzPyr.data(dphzPyr.data<0) = dphzPyr.data(dphzPyr.data<0)+2*pi;

dlfpPrc = copy(lfpPrc);
dlfpPrc.resample(txyz);
dphzPrc = dlfpPrc.phase([5,12]);
dphzPrc.data = dphzPrc.data -pi/2;
dphzPrc.data(dphzPrc.data<0) = dphzPrc.data(dphzPrc.data<0)+2*pi;


dlfpLm = copy(lfpLm);
dlfpLm.resample(txyz);
dphzLm = dlfpLm.phase([5,12]);
dphzLm.data = dphzLm.data -pi/2;
dphzLm.data(dphzLm.data<0) = dphzLm.data(dphzLm.data<0)+2*pi;

dlfpLrc = copy(lfpLrc);
dlfpLrc.resample(txyz);
dphzLrc = dlfpLrc.phase([5,12]);
dphzLrc.data = dphzLrc.data -pi/2-1;
dphzLrc.data(dphzLrc.data<0) = dphzLrc.data(dphzLrc.data<0)+2*pi;


dlfpRad = copy(lfpRad);
dlfpRad.resample(txyz);
dphzRad = dlfpRad.phase([5,12]);
dphzRad.data = dphzRad.data -pi/2;
dphzRad.data(dphzRad.data<0) = dphzRad.data(dphzRad.data<0)+2*pi;


dlfpRrc = copy(lfpRrc);
dlfpRrc.resample(txyz);
dphzRrc = dlfpRrc.phase([5,12]);
dphzRrc.data = dphzRrc.data -pi/2;
dphzRrc.data(dphzRrc.data<0) = dphzRrc.data(dphzRrc.data<0)+2*pi;

dxy = copy(txyz);
drz = compute_drz(Trial,units,pft,'feature',dxy);
ddz = compute_ddz(Trial,units,pft,'feature',dxy);
ghz = compute_ghz(Trial,units,pft,'feature',dxy);

pyr = Trial.load('spk', sampleRate, 'w+p&theta', units, 'deburst');

dphzWhl.data(dphzWhl.data>pi)= dphzWhl.data(dphzWhl.data>pi) -2*pi;
dphzPyr.data(dphzPyr.data>pi)= dphzPyr.data(dphzPyr.data>pi) -2*pi;


[PWhl,phzStatsWhl,RmaxWhl,rhoWhl] = MjgER2016_phasePrecession(Trial,ghz,ddz,dphzLrc,pyr,units,300,-pi:0.01:pi,'overwrite',true);
[PPyr,phzStatsPyr,RmaxPyr,rhoPyr] = MjgER2016_phasePrecession(Trial,ghz,ddz,dphzPyr,pyr,units,300,-pi:0.01:pi,'overwrite',true);

pyrS = copy(pyr);
pyrS.clu = pyrS.clu(abs(hba(pyrS.res))<0.3);
pyrS.res = pyrS.res(abs(hba(pyrS.res))<0.3);
pyrT = copy(pyr);
pyrT.clu = pyrT.clu(abs(hba(pyrT.res))>0.3);
pyrT.res = pyrT.res(abs(hba(pyrT.res))>0.3);

[PS,phzStatsS,RmaxS,rhoS] = MjgER2016_phasePrecession(Trial,ghz,ddz,dphzPyr,pyrS,units,300,-pi:0.01:pi,'overwrite',true);
[PT,phzStatsT,RmaxT,rhoT] = MjgER2016_phasePrecession(Trial,ghz,ddz,dphzPyr,pyrT,units,300,-pi:0.01:pi,'overwrite',true);

figure,hold('on');
plot(PS(:,1),PT(:,1),'.');
line([-5,5],[-5,5]);

figure,
hold('on');
plot(PS(:,2),PT(:,2),'.');
line([-pi,pi],[-pi,pi]);

figure,hold('on');
plot(rhoPyr(:,1,1), rhoWhl(:,1,1),'.')
line([-1,1],[-1,1]);

figure
histogram(rhoPyr(:,1,1)-rhoWhl(:,1,1),linspace(-1,1,50));


[h,p] = ttest(rhoPyr(:,1,1)-rhoWhl(:,1,1))

%    drz:        directional rate zones
%    ddz:        directional distance zones
%    phz:        Local field potential phase
%    spk:        MTASpk object which holds spike events
%    units:      list of units for computation


dphzWhl.data(dphzWhl.data<0)= dphzWhl.data(dphzWhl.data<0) + 2*pi;
dphzPyr.data(dphzPyr.data<0)= dphzPyr.data(dphzPyr.data<0) + 2*pi;


unit = 25;
res = pyr(unit,[Trial.stc{'w+p&t'}]);
figure,
dind = abs(ddz(res,units==unit))<300;
subplot(421);
plot(drz(res(dind),units==unit),dphzPrc(res(dind)),'.');
subplot(422);
plot(ghz(res(dind),units==unit),dphzPrc(res(dind)),'.');
subplot(423);
plot(drz(res(dind),units==unit),dphzPyr(res(dind)),'.');
subplot(424);
plot(ghz(res(dind),units==unit),dphzPyr(res(dind)),'.b');
hold('on');
plot(ghz(res(dind),units==unit),dphzPyr(res(dind))+2*pi,'.b');
plot(ghz(res(dind),units==unit),dphzPyr(res(dind))+4*pi,'.b');
line([-1,1],polyval(PPyr(units==unit,:,1)+[0,2*pi],[-1,1]),'Color','m');
subplot(425);
plot(drz(res(dind),units==unit),dphzRrc(res(dind)),'.');
subplot(426);
plot(ghz(res(dind),units==unit),dphzRrc(res(dind)),'.');
subplot(427);
plot(drz(res(dind),units==unit),dphzLrc(res(dind)),'.');
subplot(428);
plot(ghz(res(dind),units==unit),dphzLrc(res(dind)),'.');


dind = abs(ddz(res,units==unit))<300 ...
       & abs(hba(res))>0.4;
figure
scatter(ghz(res(dind),units==unit),dphzPyr(res(dind)),15,abs(hba(res(dind))),'filled');
colormap jet
caxis([0,1.2]);
colorbar();


dind = abs(ddz(res,units==unit))<300 ...
       & abs(hba(res))>0.4;



figure
hold('on');
dind = abs(ddz(res,units==unit))<300 & abs(hba(res))>0.3;
plot(ghz(res(dind),units==unit),dphzPyr(res(dind)),'.r')
xlim([-1,1]);
ylim([0,2*pi]);

figure,
dind = abs(ddz(res,units==unit))<300 & abs(hba(res))<0.3;
plot(ghz(res(dind),units==unit),dphzPyr(res(dind)),'.k')
xlim([-1,1]);
ylim([0,2*pi]);


figure
dind = abs(ddz(res,units==unit))<300 & abs(vxy(res))>5 & abs(hba(res))<0.2;
plot(ghz(res(dind),units==unit),dphzPyr(res(dind)),'.r')
xlim([-1,1]);
ylim([0,2*pi]);

figure
dind = abs(ddz(res,units==unit))<300 & abs(vxy(res))>5 & abs(hba(res))>0.2;
plot(ghz(res(dind),units==unit),dphzPyr(res(dind)),'.r')
xlim([-1,1]);
ylim([0,2*pi]);

figure,
dind = abs(ddz(res,units==unit))<300 & abs(vxy(res))<5 & abs(hba(res))<0.2;
plot(ghz(res(dind),units==unit),dphzPyr(res(dind)),'.k')
xlim([-1,1]);
ylim([0,2*pi]);

figure,
dind = abs(ddz(res,units==unit))<300 & abs(vxy(res))<5 & abs(hba(res))>0.2;
plot(ghz(res(dind),units==unit),dphzPyr(res(dind)),'.k')
xlim([-1,1]);
ylim([0,2*pi]);

figure,
plot(dphzWhl(res(dind)),dphzPyr(res(dind)),'.');

figure
scatter(dphzLm(res(dind)),dphzPyr(res(dind)),15,lvxy(res(dind)),'filled');
colormap jet
caxis([0,2]);

lfpRrc = Trial.load('lfp', [75,77]);
lfpRrc.data = diff(lfpR.data, 1, 2);





trialIndex = 20;
bhvState = 'walk+turn+pause&theta';
Trial   = Trials     { trialIndex };
units   = Units      { trialIndex };
meta    = sessionList( trialIndex );
Pft = pfs_2d_theta(Trial,units);
Pft     = placeFieldsNoRear{ trialIndex };
phz = load_theta_phase(Trial,sampleRate);
txyz = preproc_xyz(Trial,'trb',sampleRate);
headYawCorrection = Trial.meta.correction.headYaw;
headCenterCorrection = Trial.meta.correction.headCenter;
% COMPUTE head basis
hvec = txyz(:,'nose',[1,2])-txyz(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(headYawCorrection),-sin(headYawCorrection);...
                  sin(headYawCorrection), cos(headYawCorrection)],...
                 [2,3],...
                 [1,2]);
bvec = txyz(:,'spine_upper',[1,2])-txyz(:,'spine_lower',[1,2]);
bvec = sq(bsxfun(@rdivide,bvec,sqrt(sum(bvec.^2,3))));
bvec = cat(3,bvec,sq(bvec)*[0,-1;1,0]);
 
% $$$ hvfl = fet_href_HXY(Trial, sampleRate, [], 'trb', 2.4);
% $$$ hafl = circshift(hvfl.data,1)-hvfl.data;
%bvfl = fet_bref_BXY(Trials{trlI}, sampleRate, [], 'trb', 5);
%bafl = circshift(bvfl(:,1),1)-bvfl(:,1);
hba = fet_hba(Trial,sampleRate); % Head to Body Angle
hav = fet_hbav(Trial,sampleRate);
pch = fet_HB_pitch(Trial,sampleRate);
pch.data = pch.data(:,3);
pyr = Trial.load('spk', sampleRate, bhvState, units, 'deburst');
headAngle = sq(txyz(:,'nose',[1,2])-txyz(:,'hcom',[1,2]));
headAngle = atan2(headAngle(:,2),headAngle(:,1));
mazeAngle =  sq(txyz(:,'hcom',[1,2]));
mazeAngle = atan2(mazeAngle(:,2),mazeAngle(:,1));
hma = circ_dist(headAngle,mazeAngle);

% 11    29    33    42    49    52    54    60    75    78    80
%  6    22    24    25    37    43    44    61    63    68    77

%ED10-0817:   1    10    24    33    38    57    63    64    73   105   108
%ER06-0613:  28    35    54    61    76   107   119   175
%JG05-0316:  13    19    30    41    42    48    56    58    61    6

%bhvState = 'lbhv&theta';
bhvState = 'walk+turn+pause&theta';
unit = 25;
rmap = plot(Pft,unit);
[mxr,mxp] = Pft.maxRate(unit);
exy = [bsxfun(@plus,                                            ...
              multiprod(bsxfun(@minus,                          ...
                               mxp,                             ...
                               sq(txyz(:,'hcom',[1,2]))),  ...
                        hvec,2,[2,3]),                          ...
              0)];
% $$$ exy = [bsxfun(@plus,                                            ...
% $$$               multiprod(bsxfun(@minus,                          ...
% $$$                                mxp,                             ...
% $$$                                sq(txyz(:,'spine_upper',[1,2]))),  ...
% $$$                         bvec,2,[2,3]),                          ...
% $$$               0)];
sts = [Trial.stc{bhvState}];
sts.resample(txyz);
xts = copy(sts);
xts.cast('TimeSeries');
xts.resample(txyz);
xts.data = logical(xts.data);
 
figure
%for phzI = 1:bins.phz.count;
phzI = 3
pfsRadius = sqrt(sum(rmap(:)>2)*20*20/pi);
for hbaI = 1:bins.hba.count;
for hmaI = 1:4
ures = pyr(unit,sts);
ures = ures(within_ranges(dphzPyr(ures),bins.phz.edges([phzI,phzI+1])));
switch hmaI
  case 1
    spt = 'Center to Edge';
    ind = within_ranges(hma,[-pi/4, pi/4]) ...
          & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1]));
    ures = ures(ind(ures));%towards
    exyInd = ind & xts.data;
  case 2
    spt = 'Edge to Center';
    ind = within_ranges(hma,[-pi,-pi*3/4; pi*3/4,pi]) ...
          & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1]));
    ures = ures(ind(ures));%away
    exyInd = ind & xts.data;
  case 3
    spt = 'CW';
    ind =  within_ranges(hma,[-pi*3/4,-pi/4]) ...
         & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1]));
    ures = ures(ind(ures));%CW
    exyInd = ind & xts.data;
  case 4
    spt = 'CCW';
    ind = within_ranges(hma,[pi/4,pi*3/4]) ...
          & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1]));
    ures = ures(ind(ures)); %CCW
    exyInd = ind & xts.data;
end
subplot2(3,4,hbaI,hmaI); hold('on');
%plot(exy(exyInd,2),exy(exyInd,1),'.','Color',[0.8,0.8,0.8]);
scatter(exy(exyInd,2),exy(exyInd,1),2,hma(exyInd),'Filled');
scatter( exy(ures,2), exy(ures,1), 30, hba(ures), 'Filled','MarkerEdgeColor','k');
%plot( exy(ures,2), exy(ures,1),'.m','MarkerSize',10);
circle(0,0,pfsRadius,'-k');
Lines([],0,'k');
Lines(0,[],'k');
colorbar();colormap('hsv');caxis([-pi,pi]);xlim([-250,250]);ylim([-250,250]);daspect([1,1,1]);
title(['HBA: ',spt]);
end
end
% C2E - HBA_L(<>±Lat), HBA_R(<>±Lat), HBA_L||HBA_R (<|>±Lat)
%       rate(HBA_L(Lat<0)) vs rate(HBA_R(Lat<0)) Where Lat <= pfsRadius
%       rate(HBA_L(Lat<0)) vs rate(HBA_L(Lat>0))
%       rate(HBA_L(Lat<0)) vs rate(HBA_L(Lat>0))
%       rate(HBA_R(Lat<0)) vs rate(HBA_R(Lat>0))
% E2C - HBA_L(<>±Lat), HBA_R(<>±Lat), HBA_L||HBA_R (<|>±Lat)
% Estimated gaussian for forward field

% Head direction preference


ind_L = within_ranges(hba,bins.hba.edges([1,2])) & xts.data;
ind_R = within_ranges(hba,bins.hba.edges([3,4])) & xts.data;


usideleng = sqrt(sum(bsxfun(@minus, sq(txyz(:,'spine_upper',[1,2])),mxp).^2, 2));
tsideleng = sqrt(sum(bsxfun(@minus, sq(txyz(:,'spine_lower',[1,2])),mxp).^2, 2));
ssideleng = sqrt(sum(sq(txyz(:,'spine_upper',[1,2])-txyz(:,'spine_lower',[1,2])).^2, 2));

bma = acos( (ssideleng.^2 + tsideleng.^2 - usideleng.^2) ./ (2 .* ssideleng .* tsideleng));

headAngle = sq(txyz(:,'nose',[1,2])-txyz(:,'hcom',[1,2]));
headAngle = atan2(headAngle(:,2),headAngle(:,1));
fieldAngle =  bsxfun(@minus,sq(txyz(:,'hcom',[1,2])),mxp);
fieldAngle = atan2(fieldAngle(:,2),fieldAngle(:,1));
%hfa = circ_dist(headAngle,fieldAngle);
hfa = fieldAngle;

figure,scatter(txyz(:,'hcom',1),txyz(:,'hcom',2),30,fieldAngle,'filled');
colormap('hsv');
colorbar();





%bhvState = 'lbhv&theta';
bhvState = 'walk+turn+pause&theta';
unit = 25;
rmap = plot(Pft,unit);
[mxr,mxp] = Pft.maxRate(unit);
exy = [bsxfun(@plus,                                            ...
              multiprod(bsxfun(@minus,                          ...
                               mxp,                             ...
                               sq(txyz(:,'hcom',[1,2]))),  ...
                        hvec,2,[2,3]),                          ...
              0)];
sts = [Trial.stc{bhvState}];
sts.resample(txyz);
xts = copy(sts);
xts.cast('TimeSeries');
xts.resample(txyz);
xts.data = logical(xts.data);
headAngle = sq(txyz(:,'nose',[1,2])-txyz(:,'hcom',[1,2]));
headAngle = atan2(headAngle(:,2),headAngle(:,1));
fieldAngle =  bsxfun(@minus,sq(txyz(:,'hcom',[1,2])),mxp);
fieldAngle = atan2(fieldAngle(:,2),fieldAngle(:,1));
hfa = circ_dist(headAngle,fieldAngle);
pfa =  bsxfun(@minus,sq(txyz(:,'hcom',[1,2])),mxp);
pfa = atan2(pfa(:,2),pfa(:,1));
usideleng = sqrt(sum(bsxfun(@minus, sq(txyz(:,'hcom',[1,2])),mxp).^2, 2));
figure()
subplot(131)
%scatter(txyz(:,'hcom',1),txyz(:,'hcom',2),30,fieldAngle,'filled');
out = hist2([txyz(xts,'hcom',1),txyz(xts,'hcom',2)], ...
              bins.x.edges, ...
              bins.y.edges);
imagesc(bins.x.centers, bins.x.centers, (out./sampleRate)');
axis('xy');
colormap(gca(),'jet');
colorbar();
hold('on');
circle(mxp(1),mxp(2),pfsRadius,'r');
subplot(132)
ind =  xts.data & within_ranges(usideleng,radius);
out = hist2([hfa(ind), fieldAngle(ind)], bins.hfa.edges, bins.hfa.edges);
imagesc(bins.hfa.centers, bins.hfa.centers, (out./sampleRate)')
axis('xy');
colormap(gca(),'jet');
colorbar()
set(Lines(-pi/2,[],'k'),'LineWidth',2);
set(Lines(pi/2,[],'k'),'LineWidth',2);
subplot(133)
pfa = fieldAngle;
ind_L = within_ranges(hba,bins.hba.edges([1,2])) & ind;
ind_R = within_ranges(hba,bins.hba.edges([3,4])) & ind;
hfaInds = discretize(hfa(ind_L), bins.hfa.edges);
pfaInds = discretize(pfa(ind_L), bins.hfa.edges);
hbaLOcc = accumarray([hfaInds,pfaInds],hba(ind_L),[bins.hfa.count,bins.hfa.count],@numel);
hfaInds = discretize(hfa(ind_R), bins.hfa.edges);
pfaInds = discretize(pfa(ind_R), bins.hfa.edges);
hbaROcc = accumarray([hfaInds,pfaInds],hba(ind_R),[bins.hfa.count,bins.hfa.count],@numel);
imagescnan({bins.hfa.centers, bins.hfa.centers, (hbaROcc-hbaLOcc)'./sampleRate},'colorbarIsRequired',true,'colorMap',@jet)
axis('xy');
set(Lines(-pi/2,[],'k'),'LineWidth',2);
set( Lines(pi/2,[],'k'),'LineWidth',2);

% pfa


% $$$ radius = [0,100];
radius = [50,200];
hbaI = 1;
ures = pyr(unit,sts);
ures = ures(within_ranges(dphzPyr(ures),bins.phz.edges([phzI,phzI+2])));
% $$$ ind =  within_ranges(hma,[-pi*3/4,-pi/4]) ...
% $$$        & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
% $$$        & exy(:,2)>0 ...
% $$$        & exy(:,1)>-100; %CW
ind = within_ranges(hma,[pi/4,pi*3/4]) ... 
      & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
      & exy(:,2)>0 ...
      & exy(:,1)>-100; %CCW
ures = ures(ind(ures)); 
exyInd = ind & xts.data;
rexy = [exy(exyInd,1),exy(exyInd,2)];
sexy = exy(ures,:);
sum(within_ranges(sqrt(sum(sexy.^2,2)),radius))/sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))*sampleRate

%radius = [0,100];
radius = [25,pfsRadius];
hbaI = 1;
ures = pyr(unit,sts);
ures = ures(within_ranges(dphzPyr(ures),bins.phz.edges([phzI,phzI+2])));
% $$$ ind = within_ranges(hma,[-pi/4, pi/4]) ...
% $$$       & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
% $$$       & exy(:,2)<0; % C to E
ind = within_ranges(hma,[-pi,-pi*3/4; pi*3/4,pi]) ...
      & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
      & exy(:,2)>0; % E to C
ures = ures(ind(ures)); 
exyInd = ind & xts.data;
rexy = [exy(exyInd,1),exy(exyInd,2)];
sexy = exy(ures,:);
sum(within_ranges(sqrt(sum(sexy.^2,2)),radius))/sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))*sampleRate




%bsxfun(@minus, sq(txyz(:,'spine_upper',[1,2])),mxp)
bfa = atan2(txyz(:,'spine_upper',2)-mxp(2),txyz(:,'spine_upper',1)- mxp(1));
hfa = atan2(txyz(:,'spine_upper',2)-mxp(2),txyz(:,'spine_upper',1)-mxp(1));



figure,plot(bfa)



ind =   xts.data ...
        & within_ranges(usideleng, [0,400]) ;
%      & bma < pi/3 ...
ind_L = within_ranges(hba,bins.hba.edges([1,2])) & ind;
ind_R = within_ranges(hba,bins.hba.edges([3,4])) & ind;

figure();
subplot(121);
hold('on');
plot(txyz(ind_L,'hcom',1),txyz(ind_L,'hcom',2),'.b')
circle(mxp(1),mxp(2),pfsRadius,'-k');
subplot(122);
hold('on');
plot(txyz(ind_R,'hcom',1),txyz(ind_R,'hcom',2),'.r')
circle(mxp(1),mxp(2),pfsRadius,'-k');




ind =   xts.data;
ind_L = within_ranges(hba,bins.hba.edges([1,2])) & ind;
ind_C = within_ranges(hba,bins.hba.edges([2,3])) & ind;
ind_R = within_ranges(hba,bins.hba.edges([3,4])) & ind;

figure();
subplot(131);
hold('on');
plot(txyz(ind_L,'hcom',1),txyz(ind_L,'hcom',2),'.g')
circle(mxp(1),mxp(2),pfsRadius,'-k');
xlim([-500,500]),ylim([-500,500])
title(num2str(sum(ind_L)/sampleRate/60));
subplot(132);
hold('on');
plot(txyz(ind_C,'hcom',1),txyz(ind_C,'hcom',2),'.b')
circle(mxp(1),mxp(2),pfsRadius,'-k');
xlim([-500,500]),ylim([-500,500])
title(num2str(sum(ind_C)/sampleRate/60));
subplot(133);
hold('on');
plot(txyz(ind_R,'hcom',1),txyz(ind_R,'hcom',2),'.r')
circle(mxp(1),mxp(2),pfsRadius,'-k');
xlim([-500,500]),ylim([-500,500])
title(num2str(sum(ind_R)/sampleRate/60));


% bma as a function of hba sampled by the

% hma 

% place field direction





unit = 61;
figure,
for phzI = 1:bins.phz.count
    subplot(1,4,phzI); % ego field by theta phase
    plot(pfet{1}{phzI},unit,[],'colorbar',[2,15],false,[],true, ...
         'colorMap',@jet);
    xlim([-300,300]);
    ylim([-200,300]);
end
subplot(1,4,4); % allo field in theta state
plot(pft{20},unit,[],'colorbar','colorMap',@jet);


% PLOT  surface of rough ego field by theta phase
phzI = 3;
rmap =    plot(pfet{1}{phzI},unit,[],'colorbar',[2,15],false,[],true, ...
         'colorMap',@jet);
figure,
surface(linspace(-400,400,50),linspace(-400,400,50),rmap);
colormap jet;
    xlim([-300,300]);
    ylim([-200,300]);


% $$$ radius = [0,100];
radius = [50,200];
hbaI = 1;
ures = pyr(unit,sts);
ures = ures(within_ranges(dphzPyr(ures),bins.phz.edges([phzI,phzI+2])));
% $$$ ind =  within_ranges(hma,[-pi*3/4,-pi/4]) ...
% $$$        & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
% $$$        & exy(:,2)>0 ...
% $$$        & exy(:,1)>-100; %CW
ind = within_ranges(hma,[pi/4,pi*3/4]) ... 
      & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
      & exy(:,2)>0 ...
      & exy(:,1)>-100; %CCW
ures = ures(ind(ures)); 
exyInd = ind & xts.data;
rexy = [exy(exyInd,1),exy(exyInd,2)];
sexy = exy(ures,:);
sum(within_ranges(sqrt(sum(sexy.^2,2)),radius))/sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))*sampleRate

%radius = [0,100];
radius = [25,pfsRadius];
hbaI = 1;
ures = pyr(unit,sts);
ures = ures(within_ranges(dphzPyr(ures),bins.phz.edges([phzI,phzI+2])));
% $$$ ind = within_ranges(hma,[-pi/4, pi/4]) ...
% $$$       & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
% $$$       & exy(:,2)<0; % C to E
ind = within_ranges(hma,[-pi,-pi*3/4; pi*3/4,pi]) ...
      & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
      & exy(:,2)>0; % E to C
ures = ures(ind(ures)); 
exyInd = ind & xts.data;
rexy = [exy(exyInd,1),exy(exyInd,2)];
sexy = exy(ures,:);
sum(within_ranges(sqrt(sum(sexy.^2,2)),radius))/sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))*sampleRate


rmap = plot(Pft,11);
rmap(rmap<1) = 1;


figure,surface(log10(rmap));

xfp = Trial.load('lfp', [63,64]);
xfp.data = diff(xfp.data, 1, 2);


qfp= Trial.load('lfp', [59,60]);
zfp = Trial.load('lfp', [59,60]);
zfp.data = diff(zfp.data, 1, 2);
lts = [1:size(zfp,1)]./Trial.lfp.sampleRate;

dat = LoadBinary([Trial.spath,'/',Trial.name,'.dat'], [57:64], 96,[],[],[],[1,2^23]+round(Trial.sync.data(1)*Trial.sampleRate))';
dts = [1:size(dat,1)]./Trial.sampleRate;

qat = LoadBinary([Trial.spath,'/',Trial.name,'.dat'], [59,60], 96,[],[],[],[1,2^23]+round(Trial.sync.data(1)*Trial.sampleRate))';

% GET place field centers.
% FIND location with the lowest place cell denisity.
% GET Dat segments from that locations and determin if there is an lfp
%     signature for that location.



[res,clu,map] = LoadCluRes(fullfile(Trial.spath,Trial.name),8);
sres = res(clu==2);
sdat = LoadBinary([Trial.spath,'/',Trial.name,'.dat'], [57:64], 96,[],[],[],[-32,32]+sres)';
sdat = permute(reshape(sdat,65,[],8),[2,1,3]);


[U,S,V] = svd( reshape(sdat,[],65*8), 0);


rdat = reshape(sdat,[],65*8);

figure,
multiprod(rdat,V(:,1)

figure(),
plot(rdat*V(:,2), rdat*V(:,4),'.')


figure()
plot(rdat*V(:,2), rdat*V(:,4),'.')


figure()
hist(rdat*V(:,2), 100)


figure,plot(sq(mean(bsxfun(@minus, sdat, mean(sdat(:,[1:10],:),2)),1)));

figure();
plot(sq(std(bsxfun( @minus, sdat, mean(sdat(:,[1:10],:),2)),[],1)));


figure,plot(cumsum(diff(sq(sdat(100,:,1)))))
hold('on');,
%plot(ifft(mean(fft(cumsum(diff(sq(sdat(101,:,1))),2),[],2),1)))
plot(ifft(fft(mean(sdat(:,:,1)),[],2)))

sq(sdat(100,:,1))

figure(); hold('on');
plot(sdat(100,:,1))
plot(sdat(100,2:end,1)-[cumsum(diff(sq(sdat(100,:,1))))])


[yo,fo] = mtfft(cumsum(diff(sq(sdat(100,:,1)))), 2^7, Trial.sampleRate, 2^6);
[syo,sfo] = mtfft(sq(sdat(100,1:end-1,1)), 2^7, Trial.sampleRate, 2^6);


figure,
plot( sq( std( bsxfun( @minus, sdat, mean( sdat(:,[1:10],:), 2)), [], 1)));



figure,
plot( sq(bsxfun(@minus,sdat(10,:,:),))


figure,plot(mean(rdat(rdat*V(:,2)<0,:)))
figure,plot(mean(rdat(rdat*V(:,2)>0,:)))

mspk = sq(std(bsxfun( @minus, sdat, mean(sdat(:,[1:10],:),2)),[],1));

figure
hold('on');
for chan = 1:8
    [sys,sfs,sts] = mtcsdglong(mean(sdat(:,:,chan)),128,Trial.sampleRate,65);
    plot(sfs,log10(sys));
end
[sys,sfs,sts] = mtcsdglong(mspk(:,chan),128,Trial.sampleRate,65);
plot(sfs,log10(sys),'m');



[dys,dfs,dts] = mtcsdglong()


figure,
subplot(3,1,[1,2]);
%plot( dts, bsxfun( @plus, diff(dat,1,2)/5, 1000.*[1:7]), 'k')
hold('on');
plot( dts, bsxfun( @plus, dat, 9000 + 1000.*[1:8]));
plot( dts, ButFilter(bsxfun( @plus, diff(dat,1,2), 1000.*[1:7]), 4, 200/(0.5.*Trial.sampleRate),'low'), 'k')
plot(lts,rfp.data/10,'m')
plot( dts, ButFilter(bsxfun( @minus, diff(dat(:,[1,3]),1,2), 1000), 4, 200/(0.5.*Trial.sampleRate),'low'), 'r')
plot( dts, ButFilter(bsxfun( @minus, diff(dat(:,[2,4]),1,2), 2000), 4, 200/(0.5.*Trial.sampleRate),'low'), 'r')
plot( dts, ButFilter(bsxfun( @minus, dat(:,[4]), 3000), 4, 200/(0.5.*Trial.sampleRate),'low'), 'b')

subplot(3,1,3);
plotSTC( Trial.stc, ...
         1,...
         [],...
         {'lpause','lloc','hloc','hpause','rear','groom','sit'},...
         'cbkgrmy');
linkx()

plot( dts, ButFilter(bsxfun( @minus, diff(dat(:,[2,4]),1,2), 2000), 4, 200/(2.*Trial.sampleRate),'low'), 'r')

figure,
hold('on'),
plot( dts, dat/3+2000,    'r')
plot( dts, diff(qat,1,2),'k')
plot( dts, -diff(dat,1,2)-2000,'b')
plot( lts, rfp.data/9,    'm')
plot( lts, qfp.data/3+3000, 'g');


figure,
hold('on'),
plot(lts,xfp.data,'r')
plot(zfp.data,'k')
plot(lts,rfp.data/10,'m')
plot(qfp.data/10,'g');

%%%>>>--------------------------------------------------------------------------



%%%<<< Compute Spectrum of lfp and rfp -----------------------------------------

wfp = copy(rfp);
wfp.data = [lfp.data, rfp.data, cfp.data];

specArgsTheta = struct('nFFT', 2^9,                                          ...
                         'Fs', rfp.sampleRate,                               ...
                  'WinLength', 2^8,                                          ...
                   'nOverlap', 2^8*0.875,                                    ...
                         'NW', 3,                                            ...
                    'Detrend', [],                                           ...
                    'nTapers', [],                                           ...
                  'FreqRange', [1,32]);
overwriteARModel = false;
[mys,mfs,mts] = fet_spec(Trial, wfp, 'mtcsdglong', false, [], specArgsTheta, ...
                         overwriteARModel, true);
%%%>>>--------------------------------------------------------------------------


%%%<<< Compute Spectrum of lfp and rfp -----------------------------------------
wfp = copy(rfp);
%wfp.data = [afp];
wfp.data = wfp.data(1:600000,:);
specArgsTheta = struct('nFFT', 2^8,                                          ...
                         'Fs', rfp.sampleRate,                               ...
                  'WinLength', 2^7,                                          ...
                   'nOverlap', 2^7*0.875,                                    ...
                         'NW', 3,                                            ...
                    'Detrend', [],                                           ...
                    'nTapers', [],                                           ...
                  'FreqRange', [12,250]);
overwriteARModel = false;
[gys,gfs,gts] = fet_spec(Trial, wfp, 'mtcsdglong', true, [], specArgsTheta, ...
                         overwriteARModel, true);
%%%>>>--------------------------------------------------------------------------



%%%<<< Compute Spectrum of lfp and rfp -----------------------------------------
wfp = copy(rfp);
wfp.data = [ofp.data,afp.data];
wfp.data = wfp.data(1:2400000,:);
specArgsTheta = struct('nFFT', 2^8,                                          ...
                         'Fs', rfp.sampleRate,                               ...
                  'WinLength', 2^7,                                          ...
                   'nOverlap', 2^7*0.875,                                    ...
                         'NW', 3,                                            ...
                    'Detrend', [],                                           ...
                    'nTapers', [],                                           ...
                  'FreqRange', [12,250]);
overwriteARModel = false;
[agys,agfs,agts] = fet_spec(Trial, wfp, 'mtcsdglong', true, [], specArgsTheta, ...
                         overwriteARModel, false);
%%%>>>--------------------------------------------------------------------------


%%%<<< Compute Spectrum of lfp and rfp -----------------------------------------
wind = 1:size(rfp,1);
wfp = copy(rfp);
%wfp.data = wfp.data(wind,:);
specArgsTheta = struct('nFFT', 240,                                          ...
                         'Fs', rfp.sampleRate,                               ...
                  'WinLength', 80,                                           ...
                   'nOverlap', 80*0.875,                                     ...
                         'NW', 3,                                            ...
                    'Detrend', [],                                           ...
                    'nTapers', [],                                           ...
                  'FreqRange', [30,250]);
overwriteARModel = false;
[rhys,rhfs,rhts] = fet_spec(Trial, wfp, 'mtcsdglong', false, [], specArgsTheta, ...
                         overwriteARModel, false);





%%%<<< Compute Spectrum of lfp and rfp -----------------------------------------
wind = 1:size(rfp,1);
wfp = copy(rfp);
wfp.data = [ofp.data,afp.data];
wfp.data = wfp.data(wind,:);
specArgsTheta = struct('nFFT', 240,                                          ...
                         'Fs', rfp.sampleRate,                               ...
                  'WinLength', 80,                                           ...
                   'nOverlap', 80*0.875,                                     ...
                         'NW', 3,                                            ...
                    'Detrend', [],                                           ...
                    'nTapers', [],                                           ...
                  'FreqRange', [30,250]);
overwriteARModel = false;
[ahys,ahfs,ahts] = fet_spec(Trial, wfp, 'mtcsdglong', false, [], specArgsTheta, ...
                         overwriteARModel, false);


nhys = [];
for ciy = 1:6,
    nhys(:,:,ciy) = RectFilter(RectFilter(bsxfun(@rdivide,log10(ahys(:,:,ciy)),mean(log10(ahys(nniz(ahys(:,:,ciy)),:,ciy)))),3,3)',3,3);
end
nhys(isinf(nrhys)) = 1;

nrhys = RectFilter(RectFilter(bsxfun(@rdivide,log10(rhys.data),mean(log10(rhys(nniz(rhys.data),:,1)))),3,3)',3,3);
nrhys(isinf(nrhys)) = 1;


[mins,mval] = LocalMinimaN(-mean(nhys(2:end,:,:),3)',-1,5);

rmins = [];
rvals = [];
block.index = 1;
block.size = 2^16;
block.count = floor(rhys.size(1)/block.size);
for bindex = 1:block.count
    startIndex = (bindex - 1)* block.size + 1;
    stopIndex  = (bindex)    * block.size;
    timeIndex  = startIndex : stopIndex;
    [mins,vals] = LocalMinimaN( -nrhys(:, timeIndex, :)', -1, 5);
    rmins = cat( 1, rmins, mins + [startIndex,0]);
    rvals = cat( 1, rvals, vals);
end


ymins = [];
yvals = [];
block.index = 1;
block.size = 2^16;
block.count = floor(rhys.size(1)/block.size);
for bindex = 1:block.count
    startIndex = (bindex - 1)* block.size + 1;
    stopIndex  = (bindex)    * block.size;
    timeIndex  = startIndex : stopIndex;
    [mins,vals] = LocalMinimaN( -mean(nhys(:, timeIndex, :),3)', -1, 5);
    ymins = cat( 1, ymins, mins + [startIndex,0]);
    yvals = cat( 1, yvals, vals);
end


figure();
[nx, ny] = deal(1, 8);
for ciy = 1:6
    [ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);
    imagesc(ahts, ahfs, RectFilter(RectFilter(bsxfun(@rdivide,log10(ahys(:,:,ciy)),mean(log10(ahys(nniz(ahys(:,:,ciy)),:,ciy)))),3,3)',3,3));  Lines([],7,'k');
    axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[0.7,1.2]); grid('on');
if ciy ==1    
    hold(gca(), 'on'); plot([wind]./ffp.sampleRate, ffp.data(wind,2)./1e2+70,'-w');
    hold(gca(), 'on'); plot([wind]./rfp.sampleRate, rfp.data(wind)./1e2+90,'-m');
end
end

ciy = ciy + 1; 
[ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);  
imagesc(ahts, ahfs, mean(nhys(2:end,:,:),3));
axis(gca(), 'xy'); colormap(gca(), 'jet'); grid('on');
ciy = ciy + 1; 
[ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);  
plotSTC( Trial.stc, ...
         1,...
         [],...
         {'lpause','lloc','hloc','hpause','rear','groom','sit'},...
         'cbkgrmy');
hold(gca(), 'on');
%plot([1:size(lvxy,1)]./lvxy.sampleRate,lvxy.data);
plot([1:size(vxy,1)]./vxy.sampleRate,vxy.data/10);
ylim(gca(), [0,7]); grid('on');
title('Behavior States');
linkx()



figure();
[nx, ny] = deal(1, 4);
ciy = 0;

ciy = ciy + 1; 
[ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);
hold(gca(), 'on'); 
plot([wind]./ffp.sampleRate, dfp.data(wind,1)./1e2+70,'-r');
plot([wind]./ffp.sampleRate, cfp.data(wind,1)./1e2+70,'-k');
plot([wind]./ffp.sampleRate, dfp.data(wind,2)./1e2+70,'-b');
plot([wind]./rfp.sampleRate, rfp.data(wind)./1e2+90,'-m');
grid(gca(),'on');

ciy = ciy + 1; 
[ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);  
imagesc(rhts, rhfs, nrhys(2:end,:,:));
axis(gca(), 'xy'); colormap(gca(), 'jet'); grid('on');
hold('on');
plot(rhts(rmins(:,1)),rhfs(rmins(:,2)),'*m');
caxis([0.9,1.1])


ciy = ciy + 1; 
[ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);  
imagesc(ahts, ahfs, mean(nhys(2:end,:,:),3));
axis(gca(), 'xy'); colormap(gca(), 'jet'); grid('on');
hold('on');
plot(ahts(ymins(:,1)),ahfs(ymins(:,2)),'*m');
caxis([0.9,1.1])

ciy = ciy + 1; 
[ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);  
plotSTC( Trial.stc, ...
         1,...
         [],...
         {'lpause','lloc','hloc','hpause','rear','groom','sit'},...
         'cbkgrmy');
hold(gca(), 'on');
%plot([1:size(lvxy,1)]./lvxy.sampleRate,lvxy.data);
plot([1:size(vxy,1)]./vxy.sampleRate,vxy.data/10);
ylim(gca(), [0,7]); grid('on');
title('Behavior States');
linkx()


%STARTHERE
fid = 14;
figure,
rose(phz(round(ymins(ymins(:,2)==fid& logical(tper(ymins(:,1))),1)./rhys.sampleRate.*lfp.sampleRate)),32)

tper = Trial.stc{'t'};
tper = Trial.stc{'s&t'};
tper = [Trial.stc{'w+n+p&t'}];
tper = [Trial.stc{'a&t-s-m'}];
tper = [Trial.stc{'lloc&t'}];
%tper.cast('TimeSeries');

tper.resample(txyz);
rind =   WithinRanges(ymins(:,1), tper.data) ...
       & WithinRanges(lvxy(ymins(:,1),2), [0.5,2]);
       & WithinRanges(-yvals(:), [1.02,2]);
fid  = ymins(rind,2);
pind = ymins(rind,1);
figure,
hist2([yphz(pind),        ...
       ahfs(fid)],          ...
      linspace(0,2*pi,17),    ...
      ahfs,                 ...
      'xprob')
colormap('jet')



tper = Trial.stc{'t'};
tper = [Trial.stc{'a&t'}];
tper.cast('TimeSeries');


pvar = zeros([numel(rhfs),1]);
for fid = 1 : numel(rhfs)
    pind = round(rmins(rmins(:,2)==fid&rvals<-1.02,1)./rhys.sampleRate.*lfp.sampleRate);
    pind = pind(logical(tper(pind)));
    pvar(fid) = circ_mean(phz(pind));
end
pvar(pvar<0) = pvar(pvar<0) + 2*pi
figure,plot(rhfs,pvar);

%%%>>>--------------------------------------------------------------------------

%%%>>>--------------------------------------------------------------------------


figure,imagesc(ahts, ahfs, log10(ahys(:,:,1))');



%%%<<< Plot Spectrograph {lfp,rfp,cfp} -----------------------------------------

figure();
[nx, ny] = deal(1, 5);
 
[ix, iy] = deal( 1, 1);  subplot2(ny, nx, iy, ix);
imagesc(gts, gfs, RectFilter(RectFilter(log10(gys.data(:,:,1,1)),3,3)',3,3));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
title('CA1ori');
 
[ix, iy] = deal( 1, 2);  subplot2(ny, nx, iy, ix);
imagesc(gts, gfs, RectFilter(RectFilter(log10(gys.data(:,:,4,4)),3,3)',3,3));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
hold(gca(), 'on'); plot([1:2400000]./cfp.sampleRate, cfp.data(1:2400000)./1e2+70,'-w');
hold(gca(), 'on'); plot([1:2400000]./rfp.sampleRate, rfp.data(1:2400000)./1e2+90,'-m');
title('CA1pyr');
 
[ix, iy] = deal( 1, 3);  subplot2(ny, nx, iy, ix);
imagesc(gts, gfs, RectFilter(RectFilter(log10(gys.data(:,:,2,2)),3,3)',3,3));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
title('CA1rad');
 
[ix, iy] = deal( 1, 4);  subplot2(ny, nx, iy, ix);
imagesc(gts, gfs, angle(gys(:,:,1,4))'); Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'hsv'); caxis(gca(),[-pi,pi]); grid('on');
title('PhaseDiff(Pyr,Rec)');
 
[ix, iy] = deal( 1, 5);  subplot2(ny, nx, iy, ix);
plotSTC( Trial.stc, ...
         1,...
         [],...
         {'lpause','lloc','hloc','hpause','rear','groom','sit'},...
         'cbkgrmy');
hold(gca(), 'on');
%plot([1:size(lvxy,1)]./lvxy.sampleRate,lvxy.data);
plot([1:size(vxy,1)]./vxy.sampleRate,vxy.data/10);
ylim(gca(), [0,7]); grid('on');
title('Behavior States');
linkx()

%%%>>>--------------------------------------------------------------------------

%%%<<< Plot Spectrograph -------------------------------------------------------
figure();
[nx, ny] = deal(1, 7);

[ix, iy] = deal( 1, 1);  subplot2(ny, nx, iy, ix);
imagesc(gts, gfs, RectFilter(RectFilter(log10(gys.data(:,:,1)),3,3)',3,3));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
title('CA1ori');
 
[ix, iy] = deal( 1, 2);  subplot2(ny, nx, iy, ix);
imagesc(gts, gfs, RectFilter(RectFilter(log10(gys.data(:,:,2)),3,3)',3,3));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
hold(gca(), 'on'); plot([1:2400000]./cfp.sampleRate, dfp.data(1:2400000,2)./1e2+70,'-w');
hold(gca(), 'on'); plot([1:2400000]./rfp.sampleRate, rfp.data(1:2400000)./1e2+90,'-m');
title('CA1pyr');
 
[ix, iy] = deal( 1, 3);  subplot2(ny, nx, iy, ix);
imagesc(gts, gfs, RectFilter(RectFilter(log10(gys.data(:,:,3)),3,3)',3,3));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
title('CA1rad');
 
[ix, iy] = deal( 1, 4);  subplot2(ny, nx, iy, ix);
imagesc(gts, gfs, RectFilter(RectFilter(log10(gys.data(:,:,4)),3,3)',3,3));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
title('CA1ori');
 
[ix, iy] = deal( 1, 5);  subplot2(ny, nx, iy, ix);
imagesc(gts, gfs, RectFilter(RectFilter(log10(gys.data(:,:,5)),3,3)',3,3));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
hold(gca(), 'on'); plot([1:2400000]./ffp.sampleRate, ffp.data(1:2400000,2)./1e2+70,'-w');
hold(gca(), 'on'); plot([1:2400000]./rfp.sampleRate, rfp.data(1:2400000)./1e2+90,'-m');
title('CA1pyr');
 
[ix, iy] = deal( 1, 6);  subplot2(ny, nx, iy, ix);
imagesc(gts, gfs, RectFilter(RectFilter(log10(gys.data(:,:,6)),3,3)',3,3));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
title('CA1rad');
 
[ix, iy] = deal( 1, 7);  subplot2(ny, nx, iy, ix);
plotSTC( Trial.stc, ...
         1,...
         [],...
         {'lpause','lloc','hloc','hpause','rear','groom','sit'},...
         'cbkgrmy');
hold(gca(), 'on');
%plot([1:size(lvxy,1)]./lvxy.sampleRate,lvxy.data);
plot([1:size(vxy,1)]./vxy.sampleRate,vxy.data/10);
ylim(gca(), [0,7]); grid('on');
title('Behavior States');
linkx()
%%%>>>--------------------------------------------------------------------------


%%%<<< Plot Spectrograph -------------------------------------------------------
fgys = log10(cat(3,gys.data(:,:,1,1),gys.data(:,:,2,2),gys.data(:,:,3,3)));
fgys = RectFilter(permute(RectFilter(fgys,3,3),[2,1,3]),3,3);
fang = cat(3, ...
           (angle(gys.data(:,:,1,2))), ...
           (angle(gys.data(:,:,1,3))), ...
           (angle(gys.data(:,:,2,3))));
fang = abs(RectFilter(permute(RectFilter(fang,3,3),[2,1,3]),3,3));

fang = cat(3, ...
           (angle(gys.data(:,:,1,2))), ...
           (angle(gys.data(:,:,1,3))), ...
           (angle(gys.data(:,:,2,3))));
fang = abs(RectFilter(permute(RectFilter(fang,3,3),[2,1,3]),3,3));

fagys = log10(agys.data);
fagys = RectFilter(permute(RectFilter(fagys,3,3),[2,1,3]),3,3);


figure();
subplot(411);
imagesc(agts,agfs,mean(log10(agys(:,:,1:3)),3)');colormap(gca(),'jet');axis('xy');
subplot(412);
imagesc(agts,agfs,mean(fagys(:,:,1:6),3));colormap(gca(),'jet');axis('xy');
subplot(413);
hold(gca(), 'on'); plot([1:2400000]./ffp.sampleRate, ffp.data(1:2400000,2)./1e2+70,'-k');
hold(gca(), 'on'); plot([1:2400000]./rfp.sampleRate, rfp.data(1:2400000)./1e2+90,'-m');
subplot(414);
plotSTC( Trial.stc, ...
         1,...
         [],...
         {'lpause','lloc','hloc','hpause','rear','groom','sit'},...
         'cbkgrmy');
hold(gca(), 'on');
%plot([1:size(lvxy,1)]./lvxy.sampleRate,lvxy.data);
plot([1:size(vxy,1)]./vxy.sampleRate,vxy.data/10);
ylim(gca(), [0,7]); grid('on');
title('Behavior States');
linkx();

figure();
subplot(311);
imagesc(gts,gfs,mean(fgys,3));colormap(gca(),'jet');axis('xy');
subplot(312);
hold(gca(), 'on'); plot([1:600000]./ffp.sampleRate, ffp.data(1:600000,2)./1e2+70,'-k');
hold(gca(), 'on'); plot([1:600000]./rfp.sampleRate, rfp.data(1:600000)./1e2+90,'-m');
subplot(313);
imagesc(gts,gfs,mean(fgys,3)./(1+sum(fang,3)));colormap(gca(),'jet');axis('xy');
linkx();

figure();
[nx, ny] = deal(1, 6);
ciy = 1;
 
[ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);  ciy = ciy + 1;
imagesc(gts, gfs, mean(fgys,3));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
hold(gca(), 'on'); plot([1:600000]./ffp.sampleRate, ffp.data(1:600000,2)./1e2+70,'-w');
hold(gca(), 'on'); plot([1:600000]./rfp.sampleRate, rfp.data(1:600000)./1e2+90,'-m');
title('CA1radM');
 
[ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);  ciy = ciy + 1;
imagesc(gts, gfs, mean(fgys,3).*log10(mean(fgys,3)./std(fgys,[],3)));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); %caxis(gca(),[2.5,3.5]); grid('on');
title('CA1radM');
 
 
[ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);  ciy = ciy + 1;
imagesc(gts, gfs, fgys(:,:,1));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
title('CA1rad1');
 
[ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);  ciy = ciy + 1;
imagesc(gts, gfs, fgys(:,:,2));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
title('CA1rad2');
 
[ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);  ciy = ciy + 1;
imagesc(gts, gfs, fgys(:,:,3));  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet'); caxis(gca(),[2.5,3.5]); grid('on');
title('CA1rad3');
 
[ix, iy] = deal( 1, ciy);  subplot2(ny, nx, iy, ix);  ciy = ciy + 1;
plotSTC( Trial.stc, ...
         1,...
         [],...
         {'lpause','lloc','hloc','hpause','rear','groom','sit'},...
         'cbkgrmy');
hold(gca(), 'on');
%plot([1:size(lvxy,1)]./lvxy.sampleRate,lvxy.data);
plot([1:size(vxy,1)]./vxy.sampleRate,vxy.data/10);
ylim(gca(), [0,7]); grid('on');
title('Behavior States');
linkx()
%%%>>>--------------------------------------------------------------------------






figure,hist2([lvxy.data(:,2),rlpd],linspace(-2,2,50),linspace(-pi,pi,50))
figure,hist2([rtb,rlpd],linspace(4.6,6,30),linspace(-pi,pi,50))
figure,hist2([ltb,rlpd],linspace(4.2,6,30),linspace(-pi,pi,50))
figure,hist2([lvxy.data(:,2),rtb],linspace(-2,2,50),linspace(4.6,6,50))
figure,hist2([lvxy.data(:,2),ltb], linspace(-2,2,50), linspace(4.2,6,50))
figure,hist2([rtb,ltb],linspace(4.6,6,30),linspace(4.2,6,30))

rtfet = [llb,ltb,lhb,rlb,rtb,rhb,rlpd];




nind = [Trial.stc{'a'}];
nind.cast('TimeSeries');
nind.resample(mys);
nind = logical(nind.data & nniz(nrf.data));


xind = [1,4,5,6];
xind = [2,7];

[B,BINT,R,RINT,STATS] = regress(lvxy.data(nind,2),...
                                [ones([sum(nind),1]),nrf.data(nind,xind)]);

model.name = 'lfp2vel';
model.description = ['A linear multivariate model to estimate the ' ...
                    'speed of an animal based on normalized lfp spectral bands'];
model.b = B;
model.stats = STATS;
model.norm.labels = nrf.labels(xind);
model.norm.mean   = nrf.mean  (xind);

figure();
hold('on');
plot(10.^lvxy.data(nind,2));
plot(10.^lvxy.data(nind,2)+10.^R);

yhat = [ones([sum(nind),1]),nrf.data(nind,[1,4,5,6])]*B;
figure();
hold('on');
plot(10.^lvxy.data(nind,2));
plot((21.^yhat));

figure();
hold('on');
plot(lvxy.data(nind,2),yhat,'.');


corr([lvxy.data(nind,2),log10(21.^yhat)])
corrcoef([lvxy.data(nind,2),yhat])





figure,plot(lvxy.data(nind,2)-yhat);
hold('on');,plot(R);

lxyz = resample(copy(txyz),mys);
lvxy = resample(copy(vxy),mys);
lvxy.data(lvxy.data<=0.0001) = 0.0001;
lvxy.data = log10(lvxy.data);

vxy = vel(filter(resample(copy(txyz),mys),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);

fmys = copy(mys);
fmys.data = log10(mys.data(:,:,1,1));
fmys.data(:,:,2) = log10(mys.data(:,:,2,2));
fmys.data(:,:,3) = log10(mys.data(:,:,3,3));
fmys.filter('ButFilter',4,1,'low');

%%%<<< Plot Spectrograph -------------------------------------------------------
figure();
[nx, ny] = deal(1, 5);
[ix, iy] = deal( 1, 1);  subplot2(ny, nx, iy, ix);
imagesc(mts, mfs, log10(mys.data(:,:,1,1))');  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet')
[ix, iy] = deal( 1, 2);  subplot2(ny, nx, iy, ix);
imagesc(mts, mfs, log10(mys.data(:,:,2,2))');  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet')
[ix, iy] = deal( 1, 3);  subplot2(ny, nx, iy, ix);
imagesc(mts, mfs, log10(mys.data(:,:,3,3))');  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet')
[ix, iy] = deal( 1, 4);  subplot2(ny, nx, iy, ix);
imagesc(mts, mfs, circ_dist(angle(mys(:,:,1,2)),pi)'); Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'hsv')
[ix, iy] = deal( 1, 5);  subplot2(ny, nx, iy, ix);
plotSTC( Trial.stc, 1);
linkx()
%%%>>>--------------------------------------------------------------------------


figure();
subplot2(2,2,1,1);
    hold('on');
    ind = ':';
    plot(log10(mys(ind,14,2,2)),log10(mys(ind,2,2,2)),'.')
    ind = [Trial.stc{'w+p'}];
    plot(log10(mys(ind,14,2,2)),log10(mys(ind,2,2,2)),'.r')
    xlim([6, 9.5]);
    ylim([6, 9.5]);
subplot2(2,2,1,2);
hold('on');
    ind = ':';
    hist2([log10(mys(ind,14,2,2)),...
           log10(mys(ind,2,2,2))],...
          linspace(6, 9.5, 25),...
          linspace(6, 9.5, 25));
    axis('tight');
subplot2(2,2,2,1);
    hold('on');
    ind = ':';
    plot(log10(mys(ind,14,1,1)),log10(mys(ind,2,1,1)),'.')
    ind = [Trial.stc{'w+p'}];
    plot(log10(mys(ind,14,1,1)),log10(mys(ind,2,1,1)),'.r')
    xlim([6, 9.5]);
    ylim([6, 9.5]);
subplot2(2,2,2,2);
    hold('on');
    ind = ':';
    hist2([log10(mys(ind,14,1,1)),...
           log10(mys(ind, 2,1,1))],...
          linspace(6, 9.5, 25),...
          linspace(6, 9.5, 25));
    axis('tight');






figure();
subplot2(2,2,1,1);
    hold('on');
    ind = ':';
    plot(log10(mys(ind,14,2,2)),mean([log10(mys(ind,2,2,2)),log10(mys(ind,28,2,2))],2),'.')
    ind = [Trial.stc{'w+p'}];
    plot(log10(mys(ind,14,2,2)),mean([log10(mys(ind,2,2,2)),log10(mys(ind,28,2,2))],2),'.r')
    xlim([6, 9.5]);
    ylim([6, 9.5]);
subplot2(2,2,1,2);
hold('on');
    ind = ':';
    hist2([log10(mys(ind,14,2,2)),...
           mean([log10(mys(ind,2,2,2)),log10(mys(ind,28,2,2))],2)],...
          linspace(6, 9.5, 25),...
          linspace(6, 9.5, 25));
    axis('tight');
subplot2(2,2,2,1);
    hold('on');
    ind = ':';
    plot(log10(mys(ind,24,2,2)),log10(mys(ind,2,1,1)),'.')
    ind = [Trial.stc{'w+p'}];
    plot(log10(mys(ind,2,2,2)),log10(mys(ind,2,1,1)),'.r')
    xlim([6, 9.5]);
    ylim([6, 9.5]);
subplot2(2,2,2,2);
    hold('on');
    ind = ':';
    hist2([log10(mys(ind, 2,2,2)),...
           log10(mys(ind, 2,1,1))],...
          linspace(6, 9.5, 25),...
          linspace(6, 9.5, 25));
    axis('tight');



figure();
subplot2(2,2,1,1);
    hold('on');
    ind = ':';
    plot(log10(mys(ind,14,2,2)),mean([log10(mys(ind,2,2,2)),log10(mys(ind,28,2,2))],2),'.')
    ind = [Trial.stc{'w+p'}];
    plot(log10(mys(ind,14,2,2)),mean([log10(mys(ind,2,2,2)),log10(mys(ind,28,2,2))],2),'.r')
    xlim([6, 9.5]);
    ylim([6, 9.5]);
subplot2(2,2,1,2);
    hold('on');
    ind = ':';
    hist2([log10(mys(ind,14,2,2)),...
           mean([log10(mys(ind,2,2,2)),log10(mys(ind,28,2,2))],2)],...
          linspace(6, 9.5, 25),...
          linspace(6, 9.5, 25));
    axis('tight');
subplot2(2,2,2,1);
    hold('on');
    ind = ':';
    plot(log10(mys(ind,14,1,1)),mean([log10(mys(ind,2,1,1)),log10(mys(ind,28,1,1))],2),'.')
    ind = [Trial.stc{'w+p'}];
    plot(log10(mys(ind,14,1,1)),mean([log10(mys(ind,2,1,1)),log10(mys(ind,28,1,1))],2),'.r')
    xlim([6, 9.5]);
    ylim([6, 9.5]);
subplot2(2,2,2,2);
    hold('on');
    ind = ':';
    hist2([log10(mys(ind, 14, 1, 1)),...
           mean([log10(mys(ind,2,1,1)),log10(mys(ind,28,1,1))],2)],...
          linspace(6, 9.5, 25),...
          linspace(6, 9.5, 25));
    axis('tight');
    

%%%<<< LFP versus RFP : TD ratio -----------------------------------------------

figure();    
stss = {':','w+p+n+r','w','p&t','p-t','n','r','m','s&t','s-t'};
for ci = 1:2
    for si = 1:numel(stss)
subplot2(2,numel(stss),ci,si);        
    hold('on');
    if si == 1,
        ind = ':'
    else
        ind = [Trial.stc{stss{si}}];
    end
    hist2([log10(mys(ind,14,ci,ci)),...
           mean([log10(mys(ind,2,ci,ci)),log10(mys(ind,28,ci,ci))],2)],...
          linspace(6.5, 8.75, 25),...
          linspace(6.5, 8.75, 25));
    axis('tight');    
    Lines(8,[],'w');
    Lines([],7.5,'w');
    caxis([0,50]);    
    end
end

%%%>>>--------------------------------------------------------------------------

%%%<<< LFP versus RFP : PhaseDiff ----------------------------------------------
figure();    
stss = {':','w+p+n+r','w','p','n','r','m','s&t','s-t'};
for si = 1:numel(stss)
subplot2(1,numel(stss),1,si);        
    if si == 1,
        ind = ':'
    else
        ind = [Trial.stc{stss{si}}];
        ind.cast('TimeSeries',mys);
        ind = ind.data;
    end
    rose(mphi(ind,16,1,2),64);
end
%%%>>>--------------------------------------------------------------------------

%%%<<< LFP versus RFP : PhaseDiffChar ------------------------------------------

wscd = [];
for f = 1:numel(mfs)
    sind = [Trial.stc{'s-t'}];
    sind.cast('TimeSeries',mys);
    sind = sind.data;
    wind = [Trial.stc{'w+p+n&t'}];
    wind.cast('TimeSeries',mys);
    wind = wind.data;
    wscd(f) = circ_dist(            ...
        circ_mean(mphi(sind,f,1,2)),...
        circ_mean(mphi(wind,f,1,2)) ...
        );
    svar(f) = circ_var(mphi(sind,f,1,2));
    wvar(f) = circ_var(mphi(wind,f,1,2));
end

%%%>>>--------------------------------------------------------------------------

    
figure();
subplot(121);
    sind = [Trial.stc{'s-t'}];
    sind.cast('TimeSeries',mys);
    sind = sind.data;
    wind = [Trial.stc{'w+p+n&t'}];
    wind.cast('TimeSeries',mys);
    wind = wind.data;
    rind = [Trial.stc{'r&t'}];
    rind.cast('TimeSeries',mys);
    rind = rind.data;
    hold('on');
    plot(log10(abs(mys(sind,16,2,2))),circ_dist(angle(mys(sind,16,1,2)),pi),'.')
    plot(log10(abs(mys(wind,16,2,2))),circ_dist(angle(mys(wind,16,1,2)),pi),'.g')
    plot(log10(abs(mys(rind,16,2,2))),circ_dist(angle(mys(rind,16,1,2)),pi),'.r')
subplot(122);
    wind = [Trial.stc{'s&t'}];
    wind.cast('TimeSeries',mys);
    wind = wind.data;
    hold('on');
    plot(log10(abs(mys(sind,16,2,2))),circ_dist(angle(mys(sind,16,1,2)),pi),'.')
    plot(log10(abs(mys(wind,16,2,2))),circ_dist(angle(mys(wind,16,1,2)),pi),'.m')

    
figure();
subplot(221);
    sind = [Trial.stc{'s-t'}];
    sind.cast('TimeSeries');
    sind.resample(mys);
    sind = logical(sind.data);
    wind = [Trial.stc{'w+p+n',mys.sampleRate}];
    wind.cast('TimeSeries');
    wind.resample(mys);
    wind = logical(wind.data);
    rind = [Trial.stc{'r',mys.sampleRate}];
    rind.cast('TimeSeries');
    rind.resample(mys);    
    rind = logical(rind.data);
    hold('on');
    plot(log10(abs(mys(sind,16,2,2)))-mean([log10(mys(sind,2,2,2)),log10(mys(sind,24,2,2))],2),...
         circ_dist(angle(mys(sind,16,1,2)),pi),'.')
    plot(log10(abs(mys(wind,16,2,2)))-mean([log10(mys(wind,2,2,2)),log10(mys(wind,24,2,2))],2),...
         circ_dist(angle(mys(wind,16,1,2)),pi),'.g')
    plot(log10(abs(mys(rind,16,2,2)))-mean([log10(mys(rind,2,2,2)),log10(mys(rind,24,2,2))],2),...
         circ_dist(angle(mys(rind,16,1,2)),pi),'.r')
subplot(222);
    sind = [Trial.stc{'a-t'}];
    sind.cast('TimeSeries');
    sind.resample(mys);
    sind = logical(sind.data);
    wind = [Trial.stc{'a&t'}];
    wind.cast('TimeSeries');    
    wind.resample(mys);
    wind = logical(wind.data);
    hold('on');
    plot(log10(abs(mys(sind,16,2,2)))-mean([log10(mys(sind,2,2,2)),log10(mys(sind,24,2,2))],2),...
         circ_dist(angle(mys(sind,16,1,2)),pi),'.')
    plot(log10(abs(mys(wind,16,2,2)))-mean([log10(mys(wind,2,2,2)),log10(mys(wind,24,2,2))],2),...
         circ_dist(angle(mys(wind,16,1,2)),pi),'.m')
subplot(224);
    hold('on');
    plot(log10(abs(mys(sind,16,2,2)))-mean([log10(mys(sind,2,2,2)),log10(mys(sind,24,2,2))],2),...
         circ_dist(angle(mys(sind,16,1,2)),pi),'.')



figure();
subplot(221);
    sind = [Trial.stc{'a-t'}];
    sind.cast('TimeSeries');
    sind.resample(mys);
    sind = logical(sind.data);
    soutp = hist2([log10(abs(mys(sind,14,1,1)))-mean([log10(mys(sind,1:8,1,1)),log10(mys(sind,22:26,1,1))],2),...
           circ_dist(angle(mys(sind,14,1,2)),pi)],...
          linspace(-1.5,2,50),...
          linspace(-1.5,1.5,50));
subplot(223);
    soutr = hist2([log10(abs(mys(sind,14,2,2)))-mean([log10(mys(sind,1:8,2,2)),log10(mys(sind,22:26,2,2))],2),...
           circ_dist(angle(mys(sind,14,1,2)),pi)],...
          linspace(-1.5,2,50),...
          linspace(-1.5,1.5,50));
subplot(222);
    wind = [Trial.stc{'a&t'}];
    wind.cast('TimeSeries');    
    wind.resample(mys);
    wind = logical(wind.data);
    woutp = hist2([log10(abs(mys(wind,14,1,1)))-mean([log10(mys(wind,1:8,1,1)),log10(mys(wind,22:26,1,1))],2),...
           circ_dist(angle(mys(wind,14,1,2)),pi)],...
          linspace(-1.5,2,50),...
          linspace(-1.5,1.5,50));
subplot(224);
    woutr = hist2([log10(abs(mys(wind,14,2,2)))-mean([log10(mys(wind,1:8,2,2)),log10(mys(wind,22:26,2,2))],2),...
           circ_dist(angle(mys(wind,14,1,2)),pi)],...
          linspace(-1.5,2,50),...
          linspace(-1.5,1.5,50));
ForAllSubplots(...
    [ ...
        'Lines([],0.5,''w'');Lines(0.5,[],''w'');', ...
        'line([-0.5,2],[1.5,-1],''color'',''w'');' ...
        'caxis([0,250]);' ...
    ] ...
);

soutp = soutp./(sum(soutp(:)));
soutr = soutr./(sum(soutr(:)));
woutp = woutp./(sum(woutp(:)));
woutr = woutr./(sum(woutr(:)));

ind = nniz(woutr(:))&nniz(woutp(:));
-log(sum(sqrt(woutp(ind).*woutr(ind))))
sqrt(1-sum(sqrt(woutp(ind).*woutr(ind))))

ind = nniz(soutp(:))&nniz(woutp(:));
-log(sum(sqrt(woutp(ind).*soutp(ind))))
sqrt(1-sum(sqrt(woutp(ind).*soutp(ind))))

ind = nniz(soutp(:))&nniz(soutr(:));
-log(sum(sqrt(soutp(ind).*soutr(ind))))
sqrt(1-sum(sqrt(soutp(ind).*soutr(ind))))


ind = nniz(woutr(:))&nniz(soutr(:));
-log(sum(sqrt(woutr(ind).*soutr(ind))))
sqrt(1-sum(sqrt(woutr(ind).*soutr(ind))))


figure();
hold('on');
plot(log10(abs(mys(:,16,2,2)))-mean([log10(mys(:,1:8,2,2)),log10(mys(:,22:26,2,2))],2))
plot(log10(abs(mys(:,16,1,1)))-mean([log10(mys(:,1:8,1,1)),log10(mys(:,22:26,1,1))],2))
plot((log10(abs(mys(:,16,2,2)))-mean([log10(mys(:,1:8,2,2)),log10(mys(:,22:26,2,2))],2)) -(log10(abs(mys(:,16,1,1)))-mean([log10(mys(:,1:8,1,1)),log10(mys(:,22:26,1,1))],2)));



llb = mean(log10(mys(:,1:5,3,3)),2);
ltb = mean(log10(mys(:,10:18,3,3)),2);
lhb = mean(log10(mys(:,23:26,3,3)),2);
rlb = mean(log10(mys(:,1:5,2,2)),2);
rtb = mean(log10(mys(:,10:18,2,2)),2);
rhb = mean(log10(mys(:,23:26,2,2)),2);


llb = mean(fmys(:,1:5,3),2);
ltb = mean(fmys(:,10:18,3),2);
lhb = mean(fmys(:,23:26,3),2);
rlb = mean(fmys(:,1:5,2),2);
rtb = mean(fmys(:,10:18,2),2);
rhb = mean(fmys(:,23:26,2),2);

rlpd = circ_dist(circ_mean(angle(mys(:,14,3,2)),[],2),-pi/2);

figure,
subplot(211);
plot(mts,rlpd);
subplot(212);
plotSTC( Trial.stc, 1);
linkx();

figure,plot(lvxy(:,2),rlpd,'.')

sind = [Trial.stc{'a'}];
sind.cast('TimeSeries');
sind.resample(mys);
sind = logical(sind.data);
 
rtfet = [llb,ltb,lhb,rlb,rtb,rhb,rlpd];



%%%<<< EQUAL state Resampling and Normalization --------------------------------
sts = {'s-t','s&t','w','n','p','r','m'};
sampleInds = [];
% RESAMPLE data
for stssub = 1:numel(sts)
    sind = [Trial.stc{sts{stssub}}];
    sind.cast('TimeSeries',mys);
    sampleInds = cat(1, sampleInds, randsample( find(sind.data), 500));
end
sind = false(size(sind.data));
sind(sampleInds) = true;
srf.data = rtfet(sind,:);
srf.mean = mean(srf.data(2:end,:));
srf.std  = std(srf.data(2:end,:));
% NORMALIZE data
nsrf = bsxfun(@rdivide,                ...
              bsxfun(@minus,           ...
                     srf.data(2:end,:),...
                     srf.mean),        ...
              srf.std);
%%%>>>--------------------------------------------------------------------------

%%%<<< PCA-VARIMAX -------------------------------------------------------------
[LU,LR,FSr,VT] = erpPCA(nsrf',6);
figure()
subplot(121);
imagesc(FSr);
subplot(122);
plot(VT(:,4),'+-')
%%%>>>--------------------------------------------------------------------------

%%%<<< NORMALIZE full data set -------------------------------------------------

nrf.data = rtfet;
nrf.mean = srf.mean;
nrf.std  = srf.std;
nrf.labels = {'CA1pyr_DELTA',      ...
              'CA1pyr_THETA',      ...
              'CA1pyr_BETA' ,      ...
              'CA1prc_DELTA',      ...
              'CA1prc_THETA',      ...
              'CA1prc_BETA' ,      ...
              'CA1_PYRxPRC_PHASE'  ...
             };
% NORMALIZE data
nrf.data = bsxfun(@rdivide,                ...
              bsxfun(@minus,           ...
                     nrf.data,         ...
                     nrf.mean),        ...
              nrf.std);
nrf.data(~nniz(lxyz),:) = 0;

%%%>>>--------------------------------------------------------------------------

figure()
subplot(3,1,[1,2]);
hold('on');
plot(mts, nrf.data * FSr(:,1));
plot(mts, nrf.data * FSr(:,2),'r');
plot(mts, nrf.data * FSr(:,3),'g');
plot(mts, nrf.data * FSr(:,4),'k');
subplot(3,1,3);
plotSTC( Trial.stc, 1);
linkx();



xind = [1:size(nrf.data,2)];
clear('hmm');
updateOM = 1;
hmm.K = 5;
hmm = hmminit(nrf.data(nniz(nrf.data),xind),hmm,'full');
hmm.train.cyc = 100;
hmm.obsmodel='Gauss';
hmm.train.obsupdate = ones([ 1, hmm.K]) * updateOM;
hmm.train.init = 1;

sind = [Trial.stc{'a'}];
sind.cast('TimeSeries',mys);
sind = sind.data&nniz(nrf.data);

hmm = hmmtrain(nrf.data(sind,xind), sum(sind), hmm);


diag(hmm.P)

% COMPUTE hmm states
[decode] = hmmdecode(nrf.data(sind,xind), sum(sind), hmm);
decode.q_star = decode.q_star';
dstates = zeros([size(mys,1),1]);
dstates(sind) = decode.q_star;

% $$$ dstates = swap_state_vector_ids(dstates,3,1);
% $$$ dstates = swap_state_vector_ids(dstates,5,3);
% $$$ dstates = swap_state_vector_ids(dstates,5,2);
% $$$ dstates = swap_state_vector_ids(dstates,4,1);
% $$$ dstates(dstates==5) = 4;

figure()
subplot(5,1,1);
imagesc(mts, mfs, log10(mys.data(:,:,1,1))');  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet')
subplot(5,1,2);
imagesc(mts, mfs, log10(mys.data(:,:,2,2))');  Lines([],7,'k');
axis(gca(), 'xy'); colormap(gca(), 'jet')
subplot(5,1,3);
hold('on');
plot(mts, nrf.data * FSr(:,1));
plot(mts, nrf.data * FSr(:,2),'r');
plot(mts, nrf.data * FSr(:,3),'g');
plot(mts, nrf.data * FSr(:,4),'k');
subplot(5,1,4);
plot(mts, dstates);
ylim([0,hmm.K+1]);
subplot(5,1,5);
plotSTC( Trial.stc, 1);
ylim([1,7]);
linkx();




%    1,  2,  3,  4,  5,  6,   7
% [llb,ltb,lhb,rlb,rtb,rhb,rlpd];

figure,
for sid = 1:hmm.K
    sax = subplot(1,hmm.K,sid);
    hist2( nrf.data(dstates==sid,[1,3]),...
           linspace(-5,5,25),...
           linspace(-5,5,25));
    grid(sax,'on');
    sax.GridAlphaMode = 'manual';  sax.GridAlpha = 0.5;
    sax.GridColorMode = 'manual';  sax.GridColor = [1,1,1];
    sax.MinorGridAlphaMode = 'manual';  sax.MinorGridAlpha = 0.5;
    sax.MinorGridColorMode = 'manual';  sax.MinorGridColor = [1,1,1];
end

rdr = rtb - mean([ rlb, rhb],2);ldr = ltb - mean([ llb, lhb],2);
lrd = rdr - ldr;



bins.rdr.name = 'rdr';
bins.rdr.description = 'rec theta delta ratio';
bins.rdr.edges = linspace(-1, 1.75, 25);
bins.rdr.centers = (bins.rdr.edges(1:end-1)+bins.rdr.edges(2:end))./2;
bins.rdr.count = numel(bins.rdr.centers);

bins.rdr.name = 'rdr';
bins.rdr.description = 'rec theta delta ratio';
bins.rdr.edges = linspace(-1, 1.75, 25);
bins.rdr.centers = (bins.rdr.edges(1:end-1)+bins.rdr.edges(2:end))./2;
bins.rdr.count = numel(bins.rdr.centers);

bins.ldr.name = 'ldr';
bins.ldr.description = 'lfp theta delta ratio';
bins.ldr.edges = linspace(-1, 1.75, 25);
bins.ldr.centers = (bins.ldr.edges(1:end-1)+bins.ldr.edges(2:end))./2;
bins.ldr.count = numel(bins.ldr.centers);

bname = 'rtb';
bins.(bname).name = bname;
bins.(bname).description = 'rec theta power';
bins.(bname).edges = linspace(4, 6.75, 25);
bins.(bname).centers = (bins.(bname).edges(1:end-1)+bins.(bname).edges(2:end))./2;
bins.(bname).count = numel(bins.(bname).centers);

bins.lrd.name = 'lrd';
bins.lrd.description = 'rdr ldr difference';
bins.lrd.edges = linspace(-1, 1.75, 25);
bins.lrd.centers = (bins.lrd.edges(1:end-1)+bins.lrd.edges(2:end))./2;
bins.lrd.count = numel(bins.lrd.centers);


xvar = rtb; xlab = 'rtb';
yvar = lrd; ylab = 'lrd';
sts = {'s-t','s&t','w','n','p','r','m'};
nx = ceil(sqrt(numel(sts)));
ny = nx;
figure();
stssub = 1;
for xsub = 1 : nx
    for ysub = 1: ny
        subplot2( ny, nx, xsub, ysub);
        sind = [Trial.stc{sts{stssub}}];
        sind.cast('TimeSeries',mys);
        sind = sind.data;
        hist2([ xvar(sind),    ...
                yvar(sind) ],  ...
              bins.(xlab).edges,   ...
              bins.(ylab).edges    ...
              );
        Lines(bins.(xlab).centers(round(bins.(xlab).count.*[0.25, 0.5, 0.75])),[],'w');
        Lines([],bins.(ylab).centers(round(bins.(ylab).count.*[0.25, 0.5, 0.75])),'w');
        line(bins.(xlab).edges([1,end]), bins.(ylab).edges([1,end]),'color','w');
        title(sts{stssub});
        if stssub == numel(sts),break;end
        stssub = stssub + 1;
    end
end

figure,
subplot(211);
plot(mts,log10(mean(mys(:,[1:5],2,2),2)));
hold('on');
plot(mts,(mean(fmys(:,[1:5],2),2)),'r');
plot(mts,(mean(fmys(:,[11:20],1),2)),'k');
plot(mts,(mean(fmys(:,[11:20],2),2)),'g');
subplot(212);
plotSTC( Trial.stc, 1);
linkx()

figure();
subplot(221);
    sind = [Trial.stc{'s&t'}];
    sind.cast('TimeSeries');
    sind.resample(mys);
    sind = logical(sind.data);
    hist2([log10(abs(mys(sind,14,2,2)))-mean(log10(mys(sind,1:40,1,1)),2),...
           circ_dist(angle(mys(sind,14,1,2)),-pi/2)],...
          linspace(-1.5,2,50),...
          linspace(-1.5,1.5,50))
subplot(223);
    hist2([log10(abs(mys(sind,14,2,2)))-mean(log10(mys(sind,1:40,1,1)),2),...
           circ_dist(angle(mys(sind,14,1,2)),-pi/2)],...
          linspace(-1.5,2,50),...
          linspace(-1.5,1.5,50))
subplot(222);
    wind = [Trial.stc{'w+n+p&t'}];
    wind.cast('TimeSeries');    
    wind.resample(mys);
    wind = logical(wind.data);
    hist2([log10(abs(mys(wind,14,1,1)))-mean(log10(mys(wind,1:40,1,1)),2),...
           circ_dist(angle(mys(wind,14,1,2)),-pi/2)],...
          linspace(-1.5,2,50),...
          linspace(-1.5,1.5,50))
subplot(224);
    hist2([log10(abs(mys(wind,14,2,2)))-mean(log10(mys(wind,1:40,1,1)),2),...
           circ_dist(angle(mys(wind,14,1,2)),-pi/2)],...
          linspace(-1.5,2,50),...
          linspace(-1.5,1.5,50))
ForAllSubplots(...
    [ ...
        'Lines([],0.5,''w'');Lines(0.5,[],''w'');', ...
        'line([-0.5,2],[1.5,-1],''color'',''w'');' ...
        'caxis([0,250]);' ...
    ] ...
);



% If I see a clear shift in the
% ind = [Data.parent.stc{statePeriods}];


figure,
subplot(131);
sind = [Trial.stc{'s-t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([log10(mys(sind,16,1,1))  ,...
       log10(mys(sind,16,2,2)) ],...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(132);
sind = [Trial.stc{'s&t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([log10(mys(sind,16,1,1))  ,...
       log10(mys(sind,16,2,2)) ],...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(133);
sind = [Trial.stc{'w+p+n+r&t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([log10(mys(sind,16,1,1))  ,...
       log10(mys(sind,16,2,2)) ],...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
ForAllSubplots(...
    [ ...
        'Lines([],8,''w'');Lines(7.5,[],''w'');', ...
        'line([6.5,8.5],[6.75,8.75],''color'',''w'');' ...
         ... % 'caxis([0,150]);' ...
    ] ...
);


figure,
subplot(131);
sind = [Trial.stc{'s-t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([log10(mys(sind,16,1,1))  ,...
       log10(mys(sind,50,1,1)) ],...
      linspace(6.5, 8.75, 25),...
      linspace(4.5, 7, 25)...
);
subplot(132);
sind = [Trial.stc{'s&t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([log10(mys(sind,16,1,1))  ,...
       log10(mys(sind,50,1,1)) ],...
      linspace(6.5, 8.75, 25),...
      linspace(4.5, 7, 25)...
);
subplot(133);
sind = [Trial.stc{'w+p+n+r&t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([log10(mys(sind,16,1,1))  ,...
       log10(mys(sind,50,1,1)) ],...
      linspace(6.5, 8.75, 25),...
      linspace(4.5, 7, 25)...
);
ForAllSubplots(...
    [ ...
        'Lines([],8,''w'');Lines(7.5,[],''w'');', ...
        'line([6.5,8.5],[6.75,8.75],''color'',''w'');' ...
         ... % 'caxis([0,150]);' ...
    ] ...
);



figure,
subplot(131);
sind = [Trial.stc{'s-t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([ log10(fmys(sind,16,1))  ,...
             log10(fmys(sind,16,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(132);
sind = [Trial.stc{'s&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([ log10(fmys(sind,16,1))  ,...
             log10(fmys(sind,16,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(133);
sind = [Trial.stc{'w+p+n+r&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([ log10(fmys(sind,16,1))  ,...
             log10(fmys(sind,16,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
ForAllSubplots(...
    [ ...
        'Lines([],8,''w'');Lines(7.5,[],''w'');', ...
        'line([6.5,8.5],[6.75,8.75],''color'',''w'');' ...
         ... % 'caxis([0,150]);' ...
    ] ...
);




figure,
subplot(131);
sind = [Trial.stc{'s-t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([ log10(fmys(sind,24,1))  ,...
             log10(fmys(sind,24,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(132);
sind = [Trial.stc{'s&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([ log10(fmys(sind,24,1))  ,...
             log10(fmys(sind,24,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(133);
sind = [Trial.stc{'w+p+n+r&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([ log10(fmys(sind,24,1))  ,...
             log10(fmys(sind,24,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
ForAllSubplots(...
    [ ...
        'Lines([],8,''w'');Lines(7.5,[],''w'');', ...
        'line([6.5,8.5],[6.75,8.75],''color'',''w'');' ...
         ... % 'caxis([0,150]);' ...
    ] ...
);

figure,
subplot(131);
sind = [Trial.stc{'s-t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([ log10(fmys(sind,2,1))  ,...
             log10(fmys(sind,2,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(132);
sind = [Trial.stc{'s&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([ log10(fmys(sind,2,1))  ,...
             log10(fmys(sind,2,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(133);
sind = [Trial.stc{'w+p+n+r&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([ log10(fmys(sind,2,1))  ,...
             log10(fmys(sind,2,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
ForAllSubplots(...
    [ ...
        'Lines([],8,''w'');Lines(7.5,[],''w'');', ...
        'line([6.5,8.5],[6.75,8.75],''color'',''w'');' ...
         ... % 'caxis([0,150]);' ...
    ] ...
);

figure,
subplot(221);
sind = [Trial.stc{'s-t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([mean([log10(mys(sind,2,2,2)),log10(mys(sind,28,2,2))],2),...
       log10(mys(sind,16,2,2)) ],...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(222);
sind = [Trial.stc{'w&t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([mean([log10(mys(sind,2,2,2)),log10(mys(sind,28,2,2))],2),...
       log10(mys(sind,16,2,2)) ],...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(223);
sind = [Trial.stc{'p-t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([mean([log10(mys(sind,2,2,2)),log10(mys(sind,28,2,2))],2),...
       log10(mys(sind,16,2,2)) ],...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(224);
sind = [Trial.stc{'p&t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([mean([log10(mys(sind,2,2,2)),log10(mys(sind,28,2,2))],2),...
       log10(mys(sind,16,2,2)) ],...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);

fmys = copy(mys);
fmys.data = mys.data(:,:,1,1);
fmys.data(:,:,2) = mys.data(:,:,2,2);
fmys.filter('ButFilter',4,3,'low');

figure,
subplot(231);
sind = [Trial.stc{'s-t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([real(mean([log10(fmys(sind,2,2)),log10(fmys(sind,28,2))],2)),...
       real(log10(fmys(sind,16,2))) ],...
       linspace(6.5, 8.75, 25),...
       linspace(6.5, 8.75, 25)...
);
subplot(232);
sind = [Trial.stc{'w&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([mean([log10(fmys(sind,2,2)),log10(fmys(sind,28,2))],2),...
       log10(fmys(sind,16,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(233);
sind = [Trial.stc{'p-t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2([real(mean([log10(fmys(sind,2,2)),log10(fmys(sind,28,2))],2)),...
       real(log10(fmys(sind,16,2))) ],...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
      );
subplot(234);
sind = [Trial.stc{'s&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([mean([log10(fmys(sind,2,2)),log10(fmys(sind,28,2))],2),...
       log10(fmys(sind,16,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(235);
sind = [Trial.stc{'r&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([mean([log10(fmys(sind,2,2)),log10(fmys(sind,28,2))],2),...
       log10(fmys(sind,16,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(236);
sind = [Trial.stc{'p&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([mean([log10(fmys(sind,2,2)),log10(fmys(sind,28,2))],2),...
       log10(fmys(sind,16,2)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
ForAllSubplots('Lines([],8,''w'');Lines(7.5,[],''w'');');




figure,
subplot(231);
sind = [Trial.stc{'s-t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([real(mean([log10(fmys(sind,2,1)),log10(fmys(sind,28,1))],2)),...
       real(log10(fmys(sind,16,1))) ],...
       linspace(6.5, 8.75, 25),...
       linspace(6.5, 8.75, 25)...
);
subplot(232);
sind = [Trial.stc{'w&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([mean([log10(fmys(sind,2,1)),log10(fmys(sind,28,1))],2),...
       log10(fmys(sind,16,1)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(233);
sind = [Trial.stc{'p-t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2([real(mean([log10(fmys(sind,2,1)),log10(fmys(sind,28,1))],2)),...
       real(log10(fmys(sind,16,2))) ],...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
      );
subplot(234);
sind = [Trial.stc{'s&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([mean([log10(fmys(sind,2,1)),log10(fmys(sind,28,1))],2),...
       log10(fmys(sind,16,1)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(235);
sind = [Trial.stc{'r&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([mean([log10(fmys(sind,2,1)),log10(fmys(sind,28,1))],2),...
       log10(fmys(sind,16,1)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(236);
sind = [Trial.stc{'p&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([mean([log10(fmys(sind,2,1)),log10(fmys(sind,28,1))],2),...
       log10(fmys(sind,16,1)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
ForAllSubplots('Lines([],8,''w'');Lines(7.5,[],''w'');');


figure,
subplot(231);
sind = [Trial.stc{'s-t'}];
sind.cast('TimeSeries',mys);
sind = sind.data;
hist2([real(mean([log10(fmys(sind,2,2)),log10(fmys(sind,28,2))],2)),...
       real(log10(fmys(sind,16,1))) ],...
       linspace(6.5, 8.75, 25),...
       linspace(6.5, 8.75, 25)...
);
subplot(232);
sind = [Trial.stc{'w&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([mean([log10(fmys(sind,2,2)),log10(fmys(sind,28,2))],2),...
       log10(fmys(sind,16,1)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(233);
sind = [Trial.stc{'p-t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2([real(mean([log10(fmys(sind,2,2)),log10(fmys(sind,28,2))],2)),...
       real(log10(fmys(sind,16,2))) ],...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
      );
subplot(234);
sind = [Trial.stc{'s&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([mean([log10(fmys(sind,2,2)),log10(fmys(sind,28,2))],2),...
       log10(fmys(sind,16,1)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(235);
sind = [Trial.stc{'r&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([mean([log10(fmys(sind,2,2)),log10(fmys(sind,28,2))],2),...
       log10(fmys(sind,16,1)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
subplot(236);
sind = [Trial.stc{'p&t'}];
sind.cast('TimeSeries',fmys);
sind = sind.data;
hist2(real([mean([log10(fmys(sind,2,2)),log10(fmys(sind,28,2))],2),...
       log10(fmys(sind,16,1)) ]),...
      linspace(6.5, 8.75, 25),...
      linspace(6.5, 8.75, 25)...
);
ForAllSubplots('Lines([],8,''w'');Lines(7.5,[],''w'');');


figure,plot(log10(mys(:,16,1,1))),hold('on');plot(log10(fmys(:,16,1,1)));


xyz = preproc_xyz(Trial,'trb',sampleRate);
fxyz = filter(xyz.copy(),'ButFilter',3,14,'low');
vxy = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
vxyz = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2,3]);
lvxy = copy(vxy);
lvxy.data(lvxy.data<=0.0001) = 0.0001;
lvxy.data = log10(lvxy.data);

figure();
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxyz)]./vxy.sampleRate,vxyz(:,2)/5+1,'r');
plot([1:size(vxyz)]./vxy.sampleRate,vxyz(:,1)/5+1,'g');
xlabel('Time (s)');


Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',meta.subject.channelGroup.theta);


flfp = diff(get(Trial.load('lfp',[meta.subject.channelGroup.thetarc]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2);
%flfp = diff(get(Trial.load('lfp',[1,8]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2);

% $$$ flfp = [diff(get(Trial.load('lfp',[33,40]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2),...
% $$$         diff(get(Trial.load('lfp',[41,48]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2),...
% $$$         diff(get(Trial.load('lfp',[49,56]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2),...
% $$$         diff(get(Trial.load('lfp',[57,64]).filter('ButFilter',4,[1,200],'bandpass'),'data'),1,2)];

rlfp = Trial.load('lfp',[61]);

%rlfp & flfp
% $$$ 
% $$$ 
% $$$ figure,
% $$$ subplot(211);
% $$$ imagesc([1:size(fslfp,1)]./lfp.sampleRate,1:8,nunity(fslfp(:,:))')
% $$$ hold('on');
% $$$ plot([1:size(fslfp,1)]./lfp.sampleRate,nunity(diff(fslfp(:,[1,end]),1,2))+4.5,'m','LineWidth',2);
% $$$ plot([1:size(fslfp,1)]./lfp.sampleRate,sfslfp(:,1)+4.5,'r','LineWidth',2);
% $$$ colormap('jet');
% $$$ axis('xy');
% $$$ caxis([-3,3]);
% $$$ Lines([],4.5,'k');
% $$$ subplot(212);
% $$$ hold('on');
% $$$ plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
% $$$ ylim([1,9]);
% $$$ plot([1:size(vxyz)]./vxy.sampleRate,vxyz(:,2)/5+1,'r');
% $$$ plot([1:size(vxyz)]./vxy.sampleRate,vxyz(:,1)/5+1,'g');
% $$$ xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ 
% $$$ figure,plot(slfp(1:1e5,:))
% $$$ hold('on');
% $$$ plot(flfp(1:1e5,4),'m','LineWidth',2)








unitsInt = Trial.spk.get_unit_set(Trial,'interneurons');

int = Trial.load('spk', Trial.lfp.sampleRate, '', unitsInt, '');

% $$$ % LOAD theta phase
% $$$ pchan = [1,8,9,16,17,24,25,32,33,40,41,48,49,56,57,64];
% $$$ for c = 1:numel(pchan)
% $$$ phz = load_theta_phase(Trial,...
% $$$                        Trial.lfp.sampleRate,...
% $$$                        pchan(c),...
% $$$                        meta.subject.correction.thetaPhase);
% $$$ for ii = 1:numel(unitsInt),
% $$$ phzMean(c,ii) = circ_mean(phz(int(unitsInt(ii))));
% $$$ phzR(c,ii) = circ_r(phz(int(unitsInt(ii))));
% $$$ end
% $$$ end


%sort_interneurons_by_tpp(Trial,meta)

phz = load_theta_phase(Trial,...
                       Trial.lfp.sampleRate,...
                       meta.subject.channelGroup.theta,...
                       meta.subject.correction.thetaPhase);

% ORDER interneurons based on theta phase preference
intPhzPref = zeros([numel(unitsInt),1]);
for ii = 1:numel(unitsInt),
    intPhzPref(ii) = circ_mean(phz(int(unitsInt(ii))));
    intPhzR(ii) = circ_r(phz(int(unitsInt(ii))));
end
intPhzPref(intPhzPref<0) = intPhzPref(intPhzPref<0)+2*pi;
[mpv,mpi] = sort(intPhzPref,'descend');
mrv = intPhzR(mpi);
unitsInt = unitsInt(mpi);


%     5.22099827436091       9                     4.84224310697765    11
%     4.43219199687354      10                     4.79422023722339     7
%     3.93546757990167       7                     4.64561500480617    13
%     3.83702698028012      14                     3.65005269713263     1
%     3.78686563321746      13                     2.93272011688267    14
%     3.54861455639414       6                     2.69081170479908    10
%     3.44208801573508       8                     2.54160077488562    12
%     3.36913474655993      12                      2.3985019436955     6
%     3.02889248200655       5                     2.38542610007183     4
%     2.85983023151636       3                     2.29636991494833     9
%     2.68403989509002      11                     2.23340627953952     3
%     2.59567869954665       2                     2.19082410252469     8
%     2.59076386975001       4                     1.98523738734222     2
%     2.19652065581479       1                     1.54319131641419     5

ufr = Trial.load('ufr', lvxy, [], unitsInt, 0.12, 'boxcar');
% $$$ fufr = Trial.load('ufr', lvxy, [], unitsInt, 0.5, 'gauss');
fufr = Trial.load('ufr', lvxy, [], unitsInt, 1.5, 'gauss');
ufrp = Trial.load('ufr', lvxy, [], unitsSubset, 0.12, 'boxcar');
%ufrl = Trial.load('ufr', lvxy, [], unitsInt, 0.24, 'boxcar');
% $$$ 
% $$$ phzThresh = 3;
% $$$ figure
% $$$ subplot(211);hold('on');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,RectFilter(sum(ufr.data,2)./numel(unitsInt),9,3));
% $$$ % $$$ plot([1:size(vxy)]./vxy.sampleRate,RectFilter(sum(ufr(:,mpv>phzThresh),2)./sum(mpv>phzThresh),9,3));
% $$$ % $$$ plot([1:size(vxy)]./vxy.sampleRate,RectFilter(sum(ufr(:,mpv<3),2)./sum(mpv<3),9,3));
% $$$ plot([1:size(vxy)]./vxy.sampleRate,log10(RectFilter(sum(ufr(:,mpv>phzThresh),2)./sum(mpv>phzThresh),9,3)./RectFilter(sum(ufr(:,mpv<2.5),2)./sum(mpv<2.5),9,3)));
% $$$ plot([1:size(vxy)]./vxy.sampleRate,RectFilter(sum(ufr.data,2)./numel(unitsInt),9,3)- ...
% $$$      log10(RectFilter(sum(ufr(:,mpv>phzThresh),2)./sum(mpv>phzThresh),9,3)./RectFilter(sum(ufr(:,mpv<3),2)./sum(mpv<3),9,3)));
% $$$ Lines([],0,'k');
% $$$ Lines([],0.5,'k');
% $$$ subplot(212);
% $$$ hold('on');
% $$$ plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
% $$$ ylim([1,9]);
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
% $$$ xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ 
% $$$ 
% $$$ 

figure,
udcor = [];
for u = 1:numel(unitsInt),
    subplot(2,7,u);
% $$$     sufr = fufr(dc.ind(((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>1)),u);
% $$$     dcom = sqrt(sum(dc.ecom((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>1,:).^2,2));
    sufr = ufr(dc.ind((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1)),u);
    dcom = sqrt(sum(dc.ecom((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1),:).^2,2));
    nind = nniz(sufr) & nniz(dcom) & sufr>0.0001;
    plot(sufr(nind),dcom(nind),'.');
    udcor(u) = corr(sufr(nind),dcom(nind));
end


figure,
uvcor = [];
for u = 1:numel(unitsInt),
    subplot(2,7,u);
    sufr = fufr(dc.ind(((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2)),u);
    dcom = sqrt(sum(dc.ecom((dc.stcm(:,3)==3|dc.stcm(:,5)==5)&dc.stcm(:,1)==1&dc.ucnt>2,:).^2,2));
% $$$     sufr = ufr(dc.ind((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1)),u);
% $$$     dcom = lvxy(dc.ind((dc.stcm(:,8)==8&dc.stcm(:,1)~=1&dc.ucnt>1)),2);
    nind = nniz(sufr) & nniz(dcom) & sufr>0.0001;
    plot(sufr(nind),dcom(nind),'.');
    uvcor(u) = corr(sufr(nind),dcom(nind));
end




stcm = stc2mat(Trial.stc,xyz,{'theta','rear','loc','pause','sit','groom'});
gufr = Trial.load('ufr',xyz,int,unitsInt,0.125,'gauss');
%ind = (stcm(:,3)==3|stcm(:,4)==4);
%ind = (stcm(:,5)==5);
ind = (stcm(:,4)==4) & stcm(:,1)~=1;
suI = gufr(ind,:);
svI = lvxy(ind,:);
nClust = 2;
[idx,C,SUMD,D] = kmeans(suI,nClust,'Distance','correlation');

kCov = zeros(size(gufr,2),size(gufr,2),nClust);
for vind = 1:nClust,
    kCov(:,:,vind) = (bsxfun(@minus,suI(idx==vind,:),median(suI(idx==vind,:)))' ...
                       *bsxfun(@minus,suI(idx==vind,:),median(suI(idx==vind,:)))) ...
                      ./sum(idx==vind);
end

velBins = linspace(-2.5,1.8,20);
velInds = discretize(lvxy(:,2),velBins);
figure,
for vind = 1:nClust,
    subplot2(nClust,2,vind,1);
        imagesc(kCov(:,:,vind));
    subplot2(nClust,2,vind,2);
    histogram(svI(idx==vind,2),velBins);
end

out = [];
for vind =  1:nClust
    out(:,vind) = -.5*log(det(kCov(:,:,vind)))-0.5*(multiprod(bsxfun(@minus,gufr(:,:),median(suI(idx==vind,:))),...
                  multiprod(inv(kCov(:,:,vind)),bsxfun(@minus,gufr(:,:),median(suI(idx==vind,:))),[1,2],[2]),2,2));
end
figure,plot(out)


% $$$ figure();
% $$$ ind = (stcm(:,5)==5);
% $$$ plot(sum(gufr(ind,:),2)./14,lvxy(ind,2),'.')
% $$$ hold('on');
% $$$ ind = (stcm(:,3)==3|stcm(:,4)==4);
% $$$ plot(sum(gufr(ind,:),2)./14,lvxy(ind,2),'.r')
% $$$ 

udThresh = 0.1;
udThresh = 0.0;
irat = copy(vxy);
irat.data = log10(  sum(fufr(:,udcor>0&abs(udcor)>udThresh),2)./sum(udcor>0&abs(udcor)>udThresh) ...
                  ./sum(fufr(:,udcor<0&abs(udcor)>udThresh),2)./sum(udcor<0&abs(udcor)>udThresh));
% $$$ 
% $$$ frat = copy(vxy);
% $$$ frat.data = log10(  RectFilter(sum(ufr(:,udcor>0&abs(udcor)>udThresh),2)./sum(udcor>0&abs(udcor)>udThresh),21,3) ...
% $$$                   ./RectFilter(sum(ufr(:,udcor<0&abs(udcor)>udThresh),2)./sum(udcor<0&abs(udcor)>udThresh),21,3));
% $$$ 
% $$$ uvThresh = 0.0;
% $$$ vrat = copy(vxy);
% $$$ vrat.data = log10(  RectFilter(sum(ufr(:,udcor>0&abs(udcor)>uvThresh),2)./sum(udcor>0&abs(udcor)>uvThresh),9,3) ...
% $$$                   ./RectFilter(sum(ufr(:,udcor<0&abs(udcor)>uvThresh),2)./sum(udcor<0&abs(udcor)>uvThresh),9,3));
% $$$ 
% $$$ rat = copy(vxy);
% $$$ nrat.data = log(RectFilter(p1,21,3)./RectFilter(p2,21,3));
% $$$ 
% $$$ prat = copy(vxy);
% $$$ prat.data = log(RectFilter(p1,31,3)+RectFilter(p2,31,3));
% $$$ 
% $$$ 
% $$$ rdThresh = randperm(14);
% $$$ irat.data = log10(  RectFilter(sum(ufr(:,rdThresh(1:7)),2)./7,9,3) ...
% $$$                   ./RectFilter(sum(ufr(:,rdThresh(8:14)),2)./7,9,3));
% $$$ 
phzThresh = 2.5;
%phzThresh = 3;
trat = copy(vxy);
trat.data = log10(  (sum(fufr(:,mpv>phzThresh),2)./sum(mpv>phzThresh))  ...
                  ./(sum(fufr(:,mpv<phzThresh),2)./sum(mpv<phzThresh)));

figure
plot(unity(hfet.data),unity(sum(fufr.data,2)),'.');

figure
subplot(211);hold('on');
%plot([1:size(vxy)]./vxy.sampleRate,unity(sum(fufr.data,2)));
plot([1:size(vxy)]./vxy.sampleRate,trat.data);
Lines([],0,'k');
Lines([],0.5,'k');
subplot(212);
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

% DECODING vars -----------------------------------------------------------------------------

dc = accumulate_decoding_vars(Trial,                               ...
                              unitSubset,                          ...
                              meta.subject.channelGroup.theta,     ...
                              meta.subject.correction.thetaPhase,  ...
                              meta.subject.correction.headYaw,     ...
                              meta.subject.correction.headBody);


% $$$ figure
% $$$ plot(glfp(1:1e6,:))

% $$$ figure
% $$$ plot(diff(glfp(1:1e6,[1,2]),1,2))
% $$$ hold('on')
% $$$ plot(glfp(1:1e6,1))

% $$$ elfp = copy(lfp);
% $$$ elfp.data = diff(glfp(1:1e6,[1,3]),1,2);
% $$$ elfp = copy(lfp);
% $$$ elfp.data = RectFilter(glfp,3,3);
% $$$ specArgsTheta = struct('nFFT',2^8,...
% $$$                   'Fs',  elfp.sampleRate,...
% $$$                   'WinLength',2^7,...
% $$$                   'nOverlap',2^7*0.875,...
% $$$                   'NW',3,...
% $$$                   'Detrend',[],...
% $$$                   'nTapers',[],...
% $$$                   'FreqRange',[1,250]);
% $$$ [mys,mfs,mts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);
% $$$ 
% $$$ mxyz = preproc_xyz(Trial,'trb');
% $$$ mxyz.resample(mys);
% $$$ mfxyz = filter(mxyz.copy(),'ButFilter',3,14,'low');
% $$$ mvxy = vel(filter(mxyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
% $$$ mlvxy = copy(mvxy);
% $$$ mlvxy.data(mlvxy.data<=0.0001) = 0.0001;
% $$$ mlvxy.data = log10(mlvxy.data);
% $$$ 
% $$$ 
% $$$ mNBins = 31;
% $$$ mlvxyBinEdgs = linspace(-3,2,mNBins);
% $$$ mlvxyBinCntr = mean([mlvxyBinEdgs(1:end-1);mlvxyBinEdgs(2:end)]);
% $$$ mlvxyBinInds = discretize(mlvxy(:,2),mlvxyBinEdgs);
% $$$ 
% $$$ mspc = [];
% $$$ for v = 1:mNBins-1
% $$$     tspc = log10(mys(v == mlvxyBinInds,:,4));
% $$$     mspc(:,v) = mean(tspc(nniz(tspc),:));
% $$$ end
% $$$ 
% $$$ figure;imagesc(mfs,mlvxyBinCntr,bsxfun(@rdivide,mspc,sum(mspc,2))');axis('xy');colormap('jet');colorbar();
% $$$ figure;imagesc(mfs,mlvxyBinCntr,mspc');axis('xy');colormap('jet');colorbar();



fwin = gausswin(2^10);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = rlfp(:,1);
elfp.filter('ButFilter', 4, [20], 'low');
lpfet = copy(elfp);
lpfet.data = log10(conv(lpfet.data.^2,fwin,'same'));
lpfet.resample(xyz);

% $$$ fwin = gausswin(128);
% $$$ elfp = copy(lfp);
% $$$ elfp.data = flfp(:,1);
% $$$ elfp.filter('ButFilter', 4, [1,20], 'bandpass');
% $$$ lfet = copy(elfp);
% $$$ lfet.data = log10(conv(lfet.data.^2,fwin,'same'));
% $$$ lfet.resample(xyz);
% $$$ lfet.filter('ButFilter', 4, [1], 'low');

fwin = gausswin(64);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = flfp(:,1);
%elfp.filter('ButFilter', 4, [50,200], 'bandpass');
%elfp.filter('ButFilter', 4, [50,100], 'bandpass');
%elfp.filter('ButFilter', 4, [125,200], 'bandpass');
elfp.filter('ButFilter', 4, [75,125], 'bandpass');
hfet = copy(elfp);
hfet.data = log10(conv(hfet.data.^2,fwin,'same'));
hfet.filter('ButFilter', 4, [1], 'low');
hfet.resample(xyz);

% ifet 
ifet = copy(fufr);
ifet.data = log10(sum(fufr.data,2)./size(fufr,2));

% rfet
fwin = gausswin(2^8);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = rlfp(:,1);
elfp.filter('ButFilter', 4, [140,220], 'bandpass');
rfet = copy(elfp);
rfet.data = log10(conv(rfet.data.^2,fwin,'same'));
rfet.resample(xyz);

%dfet
fwin = gausswin(2^12);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = rlfp(:,1);
elfp.filter('ButFilter', 4, [4], 'low');
dfet = copy(elfp);
dfet.data = log10(conv(dfet.data.^2,fwin,'same'));
ttdfet = copy(dfet);
dfet.resample(xyz);
%tfet
elfp = copy(lfp);
elfp.data = rlfp(:,1);
elfp.filter('ButFilter', 4, [5,10], 'bandpass');
tfet = copy(elfp);
tfet.data = log10(conv(tfet.data.^2,fwin,'same'));
ttfet = copy(tfet);
tfet.resample(xyz);

tdfet = copy(ttfet);
tdfet.data = log10(ttfet.data./ttdfet.data);
tdfet.resample(xyz);


% drfet
fwin = gausswin(2^12);
fwin = fwin./sum(fwin);
elfp = copy(lfp);
elfp.data = flfp(:,1);
elfp.filter('ButFilter', 4, [4], 'low');
drfet = copy(elfp);
drfet.data = log10(conv(drfet.data.^2,fwin,'same'));
ttdrfet = copy(drfet);
drfet.resample(xyz);
% trfet
elfp = copy(lfp);
elfp.data = flfp(:,1);
elfp.filter('ButFilter', 4, [5,10], 'bandpass');
trfet = copy(elfp);
trfet.data = log10(conv(trfet.data.^2,fwin,'same'));
ttrfet = copy(trfet);
trfet.resample(xyz);

tdrfet = copy(ttdrfet);
tdrfet.data = log10(ttrfet.data./ttdrfet.data);
tdrfet.resample(xyz);

xts = [1:size(trfet,1)]./sampleRate;
figure,
subplot(211);hold('on');
plot(xts,tfet.data)
plot(xts,trfet.data)
Lines([],-0.05,'k')
subplot(212);
    hold('on');
    plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


figure,
ind = stcm(:,1)==1&(stcm(:,3)==3|stcm(:,4)==4|stcm(:,2)==2);
%ind = [Trial.stc{'w+p+r+n+s'}];
subplot(221);
hist2([trfet(ind),tfet(ind)],linspace(5,8,50),linspace(5,7,50));
grid('on');set(gca,'GridColor','w');
subplot(222);
hist2([dfet(ind),tfet(ind)],linspace(4,7,50),linspace(4.5,7,50));
grid('on');set(gca,'GridColor','w');
ind = stcm(:,1)==1 & stcm(:,5)==5;
subplot(223);
hist2([trfet(ind),tfet(ind)],linspace(5,8,50),linspace(5,7,50));
grid('on');set(gca,'GridColor','w');
subplot(224);
hist2([dfet(ind),tfet(ind)],linspace(4,7,50),linspace(4.5,7,50));
grid('on');set(gca,'GridColor','w');

% $$$ drfet = rfet.copy();
% $$$ drfet.data = diff([0;rfet.data]);
% $$$ ind = [stc{'s'}];
% $$$ ind = [stc{'w+n+p+r'}];
% $$$ %plot(rfet(ind),lpfet(ind),'.');
% $$$ sum(rfet(ind)<5.5&rfet(ind)>4.25 &abs(drfet(ind))<0.2)
% $$$ sum((rfet(ind)>5.5|rfet(ind)<4.25)& abs(drfet(ind))<0.2)
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ ind = [stc{'s'}];
% $$$ %plot(rfet(ind),lpfet(ind),'.');
% $$$ plot(rfet(ind),diff([0;rfet(ind)]),'.');
% $$$ ind = [stc{'w+n+p+r'}];
% $$$ %plot(rfet(ind),lpfet(ind),'.');
% $$$ plot(rfet(ind),diff([0;rfet(ind)]),'.');
% $$$ 
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ ind = [stc{'s'}];
% $$$ %plot(lpfet(ind),lpfet(ind),'.');
% $$$ plot(lpfet(ind),diff([0;lpfet(ind)]),'.');
% $$$ ind = [stc{'w+n+p+r'}];
% $$$ %plot(lpfet(ind),lpfet(ind),'.');
% $$$ plot(lpfet(ind),diff([0;lpfet(ind)]),'.');
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ ind = [stc{'s'}];
% $$$ plot(rfet(ind),lpfet(ind),'.');
% $$$ %plot(rfet(ind),diff([0;rfet(ind)]),'.');
% $$$ ind = [stc{'w+n+r'}];
% $$$ plot(rfet(ind),lpfet(ind),'.');



% $$$ fwin = gausswin(2^10);
% $$$ fwin = fwin./sum(fwin);
% $$$ elfp = copy(lfp);
% $$$ elfp.data = flfp(:,1);
% $$$ elfp.filter('ButFilter', 4, [5,12], 'bandpass');
% $$$ dfet = copy(elfp);
% $$$ dfet.data = log10(conv(dfet.data.^2,fwin,'same'));
% $$$ dfet.resample(xyz);

% $$$ sfet = copy(sfslfp);
% $$$ sfet.resample(xyz);
% $$$ dsfet = sfet.copy;
% $$$ dsfet.data = nunity(circshift(dsfet.data,1)-dsfet.data);


% $$$ 
% $$$ figure();
% $$$ subplot(211);
% $$$     hold('on');
% $$$     plot([1:size(hfet)]./lfet.sampleRate,hfet.data)
% $$$     plot([1:size(lfet)]./hfet.sampleRate,lfet.data)
% $$$ subplot(212);
% $$$     hold('on');
% $$$     plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
% $$$     ylim([1,9]);
% $$$     plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
% $$$     plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
% $$$     xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ 
% $$$ ind = [Trial.stc{'s+w+n+r+p'}];
% $$$ ind = [Trial.stc{'w+n+p'}];
% $$$ [Uf,Sf,Vf] = svd([lfet(ind)-mean(lfet(ind)),hfet(ind)-mean(hfet(ind))],0);
% $$$ Vf
% $$$ % $$$ Vf = [0.771753449181619,        -0.635921861297655;
% $$$ % $$$      0.635921861297655,         0.771753449181619];
% $$$ Vf = [-0.83534705515978,         0.549722927879022; ...
% $$$         -0.549722927879022,         -0.83534705515978];
% $$$ Vf = [0.916079244726519,        -0.400997278521052;...
% $$$          0.400997278521052,         0.916079244726519];
% $$$ 
% $$$ aind = [Trial.stc{'w+p'}];
% $$$ aind.cast('TimeSeries');
% $$$ aind.resample(xyz);
% $$$ Vv = [];
% $$$ for v = 1:19,
% $$$     ind = v==velInds & logical(aind.data);
% $$$     ml(v) = mean(lfet(ind));
% $$$     hl(v) = mean(hfet(ind));
% $$$     [Uv,Sv,Vv(:,:,v)] = svd([lfet(ind)-ml(v),hfet(ind)-hl(v)],0);
% $$$ end    
% $$$ figure,plot(abs(sq(Vv(:,1,:)))')
% $$$ figure,plot(ml,hl,'-+');
% $$$ 
% $$$ mlfet = lfet.copy();  mlfet.data = mlfet.data - mean(mlfet(ind));
% $$$ mhfet = hfet.copy();  mhfet.data = mhfet.data - mean(mhfet(ind));
% $$$ 
% $$$ mlhfet = lfet.copy();
% $$$ mlhfet.data = multiprod([mlfet.data,mhfet.data],Vf,2,[1,2]);
% $$$ 
% $$$ inds = {[Trial.stc{'s+w+r+n+p'}],[Trial.stc{'s'}],[Trial.stc{'w'}],[Trial.stc{'p'}]};
% $$$ figure,
% $$$ normf = '';      clim = [0,300];
% $$$ for s = 1:numel(inds),
% $$$     subplot(2,2,s);
% $$$     ind = inds{s};
% $$$     hist2(multiprod([mlfet(ind),mhfet(ind)],Vf,2,[1,2]),linspace(-2,2,50),linspace(-1,1,50),normf);
% $$$     Lines([],0,'m');
% $$$     Lines(0,[],'m');
% $$$     caxis(clim);
% $$$ end




elfp = copy(lfp);
elfp.data = flfp(:,1);
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  elfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[0.5,25]);
[lys,lfs,lts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);


elfp = copy(lfp);
elfp.data = rlfp(:,1);
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  elfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[0.5,25]);
[rys,rfs,rts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);

figure
subplot(311);
imagesc(lts,lfs,log10(lys.data)');
axis('xy');
colormap('jet');
caxis([3,5.5])
subplot(312);
imagesc(rts,rfs,log10(rys.data)');
axis('xy');
colormap('jet');
caxis([4,6.5])
subplot(313);
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

% $$$ % $$$ elfp = copy(lfp);
% $$$ % $$$ elfp.data = diff(glfp(1:1e6,[1,3]),1,2);
% $$$ elfp = copy(lfp);
% $$$ specArgsTheta = struct('nFFT',2^11,...
% $$$                   'Fs',  elfp.sampleRate,...
% $$$                   'WinLength',2^1s0,...
% $$$                   'nOverlap',2^10*0.875,...
% $$$                   'NW',3,...
% $$$                   'Detrend',[],...
% $$$                   'nTapers',[],...
% $$$                   'FreqRange',[1,30]);
% $$$ [rys,rfs,rts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);


% $$$ elfp = copy(lfp);
% $$$ specArgsTheta = struct('nFFT',2^8,...
% $$$                   'Fs',  elfp.sampleRate,...
% $$$                   'WinLength',2^7,...
% $$$                   'nOverlap',2^7*0.875,...
% $$$                   'NW',3,...
% $$$                   'Detrend',[],...
% $$$                   'nTapers',[],...
% $$$                   'FreqRange',[50,280]);
% $$$ [tys,tfs,tts] = fet_spec(Trial,elfp,[],true,[],specArgsTheta);
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure
% $$$ subplot(411);
% $$$ imagesc(mts,mfs,log10(mean(mys.data,3))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(412);
% $$$ imagesc(mts,mfs,log10(std(mys.data,[],3))');
% $$$ %imagesc(tts,tfs,log10(tys.data)');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(413);
% $$$ imagesc(mts,mfs,log10(mean(mys.data,3))'.^2./log10(std(mys.data,[],3))');
% $$$ %imagesc(tts,tfs,log10(tys.data)');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(414);
% $$$ hold('on');
% $$$ plotSTC(Trial.stc,1,'text',states,'krggbb');
% $$$ ylim([1,9]);
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
% $$$ xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ 
% $$$ figure
% $$$ subplot(511);
% $$$ imagesc(mts,mfs,nunity(log10(mys(:,:,1)))');
% $$$ hold('on');
% $$$ plot((1:size(elfp,1))./elfp.sampleRate,nunity(elfp.data(:,3))*10+30,'w');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(512);
% $$$ imagesc(mts,mfs,nunity(log10(mys(:,:,2)))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(513);
% $$$ imagesc(mts,mfs,nunity(log10(mys(:,:,3)))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(514);
% $$$ imagesc(mts,mfs,nunity(log10(mys(:,:,4)))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ subplot(515);
% $$$ hold('on');
% $$$ plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
% $$$ ylim([1,9]);
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
% $$$ xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ 
% $$$ figure
% $$$ subplot(511);
% $$$ imagesc(mts,mfs,(log10(mys(:,:,2)))');
% $$$ hold('on');
% $$$ plot((1:size(elfp,1))./elfp.sampleRate,nunity(elfp.data(:,3))*10+30,'w');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ grid('on')
% $$$ caxis([3.75,5.25])
% $$$ subplot(512);
% $$$ imagesc(mts,mfs,(log10(mys(:,:,3)))');
% $$$ %imagesc(mts,mfs,(log10(mys(:,:,2)))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ grid('on')
% $$$ caxis([3.75,5.25])
% $$$ subplot(513);
% $$$ imagesc(mts,mfs,(log10(mys(:,:,4)))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ grid('on')
% $$$ caxis([4,5.25])
% $$$ subplot(514);
% $$$ hold('on');
% $$$ imagesc(mts,mfs,bsxfun(@rdivide,RectFilter(log10(mys(:,:,4)),3,1),sum(RectFilter(log10(mys(:,:,4)),3,1),2))');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ grid('on');
% $$$ caxis([0.019,0.022]);
% $$$ plot((1:size(elfp,1))./elfp.sampleRate,nunity(elfp.data(:,3))*10+30,'w');
% $$$ subplot(515);
% $$$ hold('on');
% $$$ plotSTC(Trial.stc,1,'text',states,'krggbb');
% $$$ ylim([1,9]);
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,(tdRatio(:,1)+1)*5,'r');
% $$$ xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ figure
% $$$ subplot(211);
% $$$ %plot(mts,mean((log10(mys(:,:,3))),2)');
% $$$ hold('on');
% $$$ plot(mts,mean((log10(mys(:,mfs<200&mfs>40,3))),2));
% $$$ plot(mts,mean(RectFilter(log10(mys(:,mfs<200&mfs>40,3)),11,3),2));
% $$$ %plot(mts,mean(RectFilter(log10(mys(:,mfs<100,4)),7,3),2));
% $$$ plot(mts,mean(log10(mys(:,mfs<20,3)),2),'r');
% $$$ plot(mts,mean(RectFilter(log10(mys(:,mfs<20,3)),11,3),2),'k');
% $$$ plot(mts,mean(RectFilter(log10(mys(:,mfs<20,3)),21,3),2),'c');
% $$$ Lines([],4,'k');
% $$$ Lines([],4.5,'k');
% $$$ Lines([],5,'k');
% $$$ subplot(212);
% $$$ hold('on');
% $$$ plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
% $$$ ylim([1,9]);
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
% $$$ plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
% $$$ xlabel('Time (s)');
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
% $$$ 
% $$$ 
% $$$ fmys = copy(mys);
% $$$ fmys.data =  RectFilter(log10(mys(:,:,:)),21,3);
% $$$ 
% $$$ rcInd = 1;
% $$$ tpm = copy(mys);
% $$$ tpm.data = mean(fmys(:,mfs<20,rcInd),2);
% $$$ tpm.resample(vxy);
% $$$ 
% $$$ mpm = copy(mys);
% $$$ mpm.data = mean(fmys(:,mfs>140 & mfs<200,rcInd),2);
% $$$ mpm.resample(vxy);
% $$$ 
% $$$ ipm = copy(mys);
% $$$ ipm.data = mean(fmys(:,mfs>40 & mfs<200,rcInd),2);
% $$$ ipm.resample(vxy);
% $$$ 
% $$$ gpm = copy(mys);
% $$$ gpm.data = mean(fmys(:,mfs>50 & mfs<100,rcInd),2);
% $$$ gpm.resample(vxy);
% $$$ 
% $$$ hpm = copy(mys);
% $$$ hpm.data = mean(fmys(:,mfs>200,rcInd),2);
% $$$ hpm.resample(vxy);
% $$$ 
% $$$ tdpow = copy(rys);
% $$$ tdpow.resample(vxy);
% $$$ tdpow.data = log(abs(log(mean(tdpow(:,rfs>5 & rfs<12),2))-log(mean(tdpow(:,rfs<5 | (rfs>12 & rfs<15)),2))));
% $$$ 
% $$$ tpow = copy(rys);
% $$$ tpow.resample(vxy);
% $$$ tpow.data = log(mean(tpow(:,rfs>5 & rfs<12),2));
% $$$ 
% $$$ 
% $$$ % ERROR in stc s+t is not correct
% $$$ figure
% $$$ normf = '';      clim = [0,500];
% $$$ % $$$ normf = '';      clim = [0,100];
% $$$ % $$$ normf = 'xprob'; clim = [0,0.15];
% $$$ ind = [Trial.stc{'s+w+p+r+n'}];
% $$$ %ind = [Trial.stc{'w+p+r+n'}];
% $$$ %ind = [Trial.stc{'r'}];
% $$$ %ind = [Trial.stc{'s'}];
% $$$ subplot(251);
% $$$ hist2([tpm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(252);
% $$$ hist2([gpm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(253);
% $$$ hist2([mpm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(254);
% $$$ hist2([ipm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(255);
% $$$ hist2([hpm(ind),lvxy(ind,2)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(256);
% $$$ hist2([tpm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(257);
% $$$ hist2([gpm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(258);
% $$$ hist2([mpm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(259);
% $$$ hist2([ipm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ subplot(2,5,10);
% $$$ hist2([hpm(ind),lvxy(ind,1)],linspace(3.25,5.5,50),linspace(-3,2,50),normf);
% $$$ caxis(clim)
% $$$ 
% $$$ 
% $$$ 
% $$$ ind = [Trial.stc{'s+w+r+n+p'}];
% $$$ tpmm = copy(tpm); tpmm.data = tpmm.data-mean(tpmm(ind));
% $$$ ipmm = copy(ipm); ipmm.data = ipmm.data-mean(ipmm(ind));
% $$$ %mpmm = copy(mpm); mpmm.data = mpmm.data-mean(mpmm(ind));
% $$$ 
% $$$ 
% $$$ %ind = [Trial.stc{'s'}];
% $$$ figure
% $$$ hist2([tpmm(ind),ipmm(ind)],linspace(-1.25,1.25,50),linspace(-1.25,1.25,50),normf);
% $$$ caxis(clim)
% $$$ 
% $$$ ind = [Trial.stc{'s+w+r+n+p'}];
% $$$ [U,S,V] = svd([tpm(ind)-mean(tpm(ind)),ipm(ind)-mean(ipm(ind))],0);
% $$$ V = [0.771753449181619,        -0.635921861297655;
% $$$      0.635921861297655,         0.771753449181619];
% $$$ 
% $$$ 
% $$$ %[U,S,V] = svd([tpm(ind)-mean(tpm(ind)),mpm(ind)-mean(mpm(ind))],0);
% $$$ % $$$ figure,
% $$$ % $$$ plot(U(:,1),U(:,2),'.')
% $$$ % $$$ 
% $$$ % $$$ figure
% $$$ % $$$ hist2([U(:,1),U(:,2)],linspace(-3.2e-3,-2.2e-3,50),linspace(-8e-3,8e-3,50))
% $$$ 
% $$$ 
% $$$ figure,
% $$$ normf = '';
% $$$ subplot(221);
% $$$ ind = [Trial.stc{'w+r+n+p+s'}];
% $$$ hist2(multiprod([tpmm(ind),ipmm(ind)],V,2,[1,2]),linspace(-1.25,1.25,50),linspace(-.5,.5,50),normf);
% $$$ %line([-6.85,-5.5],[-0.4,0.4],'Color','r');
% $$$ line([-.7,0.8],[0.5,-0.5],'Color','r');
% $$$ subplot(222);
% $$$ ind = [Trial.stc{'p'}];
% $$$ hist2(multiprod([tpmm(ind),ipmm(ind)],V,2,[1,2]),linspace(-1.25,1.25,50),linspace(-.5,.5,50),normf);
% $$$ %line([-6.85,-5.5],[-0.4,0.4],'Color','r');
% $$$ line([-.7,0.8],[0.5,-0.5],'Color','r');
% $$$ subplot(223);
% $$$ ind = [Trial.stc{'w+r+n'}];
% $$$ hist2(multiprod([tpmm(ind),ipmm(ind)],V,2,[1,2]),linspace(-1.25,1.25,50),linspace(-.5,.5,50),normf);
% $$$ %line([-6.85,-5.5],[-0.4,0.4],'Color','r');
% $$$ line([-.7,0.8],[0.5,-0.5],'Color','r');
% $$$ subplot(224);
% $$$ ind = [Trial.stc{'s'}];
% $$$ hist2(multiprod([tpmm(ind),ipmm(ind)],V,2,[1,2]),linspace(-1.25,1.25,50),linspace(-.5,.5,50),normf);
% $$$ %line([-6.85,-5.5],[-0.4,0.4],'Color','r');
% $$$ line([-.7,0.8],[0.5,-0.5],'Color','r');
% $$$ 
% $$$ 
% $$$ dc = accumulate_decoding_vars(Trial,                               ...
% $$$                               units{trialId},                      ...
% $$$                               sessionList(trialId).thetaRefGeneral,...
% $$$                               phzCorrection(trialId),              ...
% $$$                               headRotation{trialId},               ...
% $$$                               hbangCorrection{trialId});
% $$$ 
% $$$ figure
% $$$ subplot(141);hold('on');
% $$$     ind = [Trial.stc{'s'}]; 
% $$$     set(histogram(lpfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'); 
% $$$     ind = [Trial.stc{'w'}]; 
% $$$     set(set(histogram(lpfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     xlim([5.5,7.5]);
% $$$ subplot(142);hold('on');
% $$$     ind = [Trial.stc{'s'}]; 
% $$$     set(histogram(tdfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'); 
% $$$     ind = [Trial.stc{'w'}]; 
% $$$     set(set(histogram(tdfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     xlim([-0.1,0.2]);
% $$$ subplot(143);hold('on');
% $$$     ind = [Trial.stc{'s'}]; 
% $$$     set(histogram(rfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'); 
% $$$     ind = [Trial.stc{'w'}]; 
% $$$     set(set(histogram(rfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     xlim([4,6]);
% $$$ subplot(144);hold('on');
% $$$     ind = [Trial.stc{'s'}]; 
% $$$     set(histogram(dfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'); 
% $$$     ind = [Trial.stc{'w'}]; 
% $$$     set(set(histogram(dfet(ind,1),100,'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     xlim([5,7.5]);

    
% $$$ 
% $$$ mfufr = copy(fufr);
% $$$ mfufr.data = sum(fufr(ind,:),2)./size(fufr,2);
% $$$ figure
% $$$ fufrSteps = 0:4:16;
% $$$ for f = 1:size(afet,2)
% $$$ subplot(1,size(afet,2),f);hold('on');
% $$$     for u = 1:numel(fufrSteps)-1
% $$$         ind = mfufr(:,1)>fufrSteps(u) & mfufr(:,1)<fufrSteps(u+1);
% $$$         set(set(histogram(afet(ind,f)-afetMeanW(f),linspace([afetLims(f,:),100]),'Normalization', ...
% $$$                           'probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     end
% $$$     ind = [Trial.stc{'w'}]; 
% $$$     set(set(histogram(afet(ind,f)-afetMeanW(f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     xlim(afetLims(f,:));
% $$$ end    


% $$$ figure,
% $$$ ind = ':';
% $$$ hist2([afet(ind,3),sum(fufr(ind,:),2)./size(ufr,2)],25,25);


afet = copy(lpfet);    
%afet.data = [lpfet(:,1),tdfet(:,1),rfet(:,1)];
%afet.data = [lpfet(:,1),tdfet(:,1),rfet(:,1),dfet(:,1),mlhfet(:,1),mlhfet(:,2)];
%afet.data = [lpfet(:,1),tdfet(:,1),rfet(:,1),dfet(:,1),lfet(:,1),hfet(:,1)];
%afet.data = [lpfet(:,1),tdfet(:,1),rfet(:,1),lfet(:,1),hfet(:,1)];
%afet.data = [lpfet(:,1),tdfet(:,1),tdrfet(:,1),rfet(:,1),lfet(:,1),hfet(:,1)];
%afet.data = [lpfet(:,1),tdfet(:,1),tdrfet(:,1),rfet(:,1),tfet(:,1),trfet(:,1),lfet(:,1),hfet(:,1)];
%afet.data = [dfet(:,1), lpfet(:,1), rfet(:,1), trfet(:,1), hfet(:,1), trat(:,1), ifet(:,1)];
afet.data = [dfet(:,1), lpfet(:,1), rfet(:,1), trfet(:,1), hfet(:,1)];
afetMeanW = median(afet(Trial.stc{'w'},:));
afet.data = bsxfun(@minus,afet.data,afetMeanW);
afetLims = [-1,1.5; -1,2;  -0.2,0.2; -1.5,1.25; -1.5,1.5; -1.5,1.5;-1.5,0.5; -1,0.8; -1.5,0.25  ];




figure();
for a = 1:size(afet,2)
subplot(2,4,a);
hold('on');
ind = stcm(:,1)==1&(stcm(:,3)==3|stcm(:,4)==4|stcm(:,2)==2);
set(histogram(afet(ind,a),linspace([afetLims(a,:),50]),'Normalization','probability'),'EdgeColor','none'); 
ind = stcm(:,1)==1 & stcm(:,5)==5;
set(histogram(afet(ind,a),linspace([afetLims(a,:),50]),'Normalization','probability'),'EdgeColor','none'); 
ind = stcm(:,1)~=1 & stcm(:,5)==5;
set(histogram(afet(ind,a),linspace([afetLims(a,:),50]),'Normalization','probability'),'EdgeColor','none'); 
end

%afetLims = [-1,1.5;  -0.175,0.1; -1,2; -1.5,1.25; -1.5,1.5; -1.5,1.5,];
% $$$ 
% $$$ 
% $$$ figure
% $$$ for f = 1:size(afet,2)
% $$$ subplot(1,size(afet,2),f);hold('on');
% $$$     ind = [Trial.stc{'s'}]; 
% $$$     set(histogram(afet(ind,f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'); 
% $$$     ind = [Trial.stc{'w'}]; 
% $$$     set(set(histogram(afet(ind,f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     ind = [Trial.stc{'n+p+r'}]; 
% $$$     set(set(histogram(afet(ind,f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     %ind = abs(nunity(afet(:,3)))<2 & sum(ufr(:,:),2)./size(ufr,2)>0.8;
% $$$     %ind = sum(fufr(ind,:),2)./size(fufr,2)>10;
% $$$     set(set(histogram(afet(ind,f),linspace([afetLims(f,:),100]),'Normalization','probability'),'EdgeColor','none'),'FaceAlpha',0.5); 
% $$$     xlim(afetLims(f,:));
% $$$ end    
% $$$ 
% $$$ nfet = copy(lpfet);
% $$$ nfet.data = [lpfet(:,1),tdrfet(:,1),rfet(:,1),lfet(:,1),hfet(:,1)];
% $$$ figure
% $$$ nind = nniz(nfet);
% $$$ for a = 1:size(nfet,2)
% $$$     for b = a+1:size(nfet,2)
% $$$         subplot2(size(nfet,2),size(nfet,2),a,b);
% $$$         hist2([nfet(nind,a),nfet(nind,b)],50,50);
% $$$         caxis([0,1000]);
% $$$     end
% $$$ end
% $$$ 
% $$$ 

figure
%nind = nniz(afet);
%nind = stcm(:,1)==1 & stcm(:,5)==5;
nind = stcm(:,1)~=1 & stcm(:,5)==5;
nind = stcm(:,1)==1 & (stcm(:,3)==3|stcm(:,4)==4);
for a = 1:size(afet,2)
    for b = a+1:size(afet,2)
        subplot2(size(afet,2),size(afet,2),a,b);
        hist2([afet(nind,a),afet(nind,b)],linspace([afetLims(a,:),30]),linspace([afetLims(b,:),30]));
        caxis([0,100]);        
    end
end

nind = stcm(:,1)==1 & (stcm(:,3)==3|stcm(:,4)==4);
[U,S,V] = svd(afet(nind,:),0);

[LU, LR, FSr, VT] = erpPCA( afet(nind,:), 5);

figure,
subplot(121);
imagesc(LR')
subplot(122);
plot(VT(:,4))


clear('hmm');
updateOM = 1;
hmm.K = 6;
hmm = hmminit(afet.data,hmm,'full');
hmm.train.cyc = 100;
hmm.obsmodel='Gauss';
hmm.train.obsupdate=ones([1,hmm.K])*updateOM;
hmm.train.init = 1;
hmm = hmmtrain(afet.data,size(afet,1),hmm);

% trcpow hrcpow tpow

save(fullfile(Trial.path.project,'analysis',['test_hmm_model.mat']),'hmm');
load(fullfile(Trial.path.project,'analysis',['test_hmm_model.mat']),'hmm');
 
diag(hmm.P)

% COMPUTE hmm states
[decode] = hmmdecode(afet.data,size(afet,1),hmm);
decode.q_star = decode.q_star';


figure
subplot(311);
ind = stcm(:,1)~=1&stcm(:,5)==5;
bar(accumarray(decode.q_star(ind),ones([sum(ind),1]),[hmm.K,1],@sum))
subplot(312);
ind = stcm(:,1)==1&stcm(:,5)==5;
bar(accumarray(decode.q_star(ind),ones([sum(ind),1]),[hmm.K,1],@sum))
subplot(313);
ind = stcm(:,1)==1&(stcm(:,2)==2|stcm(:,3)==3|stcm(:,4)==4);
bar(accumarray(decode.q_star(ind),ones([sum(ind),1]),[hmm.K,1],@sum))


afetGmean = [];
for grp = 1:hmm.K
    for aft = 1:size(afet,2)
        afetGmean(grp,aft) = mean(afet(decode.q_star==grp,aft));
    end
end
figure,
imagesc(afetGmean');

grpSymbol = '*o^s+<d';
figure();hold('on');
    ind = [Trial.stc{'w+r+n+p+s'}]; cast(ind,'TimeSeries'); resample(ind,vxy); ind = logical(ind.data);
    hist2([xcomp.data(ind),ycomp.data(ind)],xcomp.edgs,ycomp.edgs);
    grid('on'); Lines([],0,'w'); Lines(0,[],'w'); colormap('parula'); colorbar();
    for grp = 1:hmm.K
        plot(mean(xcomp.data(decode.q_star==grp)),mean(ycomp.data(decode.q_star==grp)),['m',grpSymbol(grp)]);
    end
    xlim(xcomp.edgs([1,end]));
    ylim(ycomp.edgs([1,end]));    

% ASSIGN new group order (automate later)
grpNewOrder = [5,1,4,2,6,3];


grpLabels = {'activeTheta','ImmobileTheta','LowTheta','REMTheta','ImmobileStuff','SWR'};

% REORDER grps
nhmmState = decode.q_star+hmm.K;
for grp = 1:hmm.K
    hmmState((grpNewOrder(grp)+hmm.K)==hmmState) = grp;
end
thmmState = hmmState;

thmmState(thmmState==2) = 1;
thmmState(thmmState==3) = 1;
thmmState(thmmState==4) = 1;

thmmState(~xor(thmmState==1,thmmState==4)) = 3;
%figure,plot(thmmState)
thmmState(thmmState==4) = 2;

thmmStateSegs = circshift(GetSegs(thmmState,1:size(thmmState,1),301,0),151,2);
thmmStateSegs = circshift(GetSegs(thmmState,1:size(thmmState,1),301,0),151,2);

mhmmState = sq(mode(thmmStateSegs))';

fhmmState = hmmState;
fhmmState(thmmState==2&mhmmState==1) = 3;
fhmmState(thmmState==1&mhmmState==2) = 4;
fhmmState(hmmState==5&mhmmState==3) = 5;
fhmmState(hmmState==5&mhmmState==1) = 5;
fhmmState(hmmState==6&(mhmmState==1|mhmmState==2)) = 6;


figure,
subplot(311);
    hold('on');
    plot([1:size(hmmState,1)]./sampleRate,hmmState)
    plot([1:size(hmmState,1)]./sampleRate,thmmState)
    %plot([1:size(hmmState,1)]./sampleRate,sq(mode(thmmStateSegs)))
subplot(312);
    hold('on');
    plot([1:size(hmmState,1)]./sampleRate,hmmState)
subplot(313);
    hold('on');
    plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


derr = nan([size(xyz,1),1]);
derr(dc.ind(dc.ucnt>2)) = sqrt(sum(dc.esax(dc.ucnt>2,:).^2,2));
derf = nan([size(xyz,1),1]);
derf(dc.ind(dc.ucnt>2)) = dc.ecom(dc.ucnt>2,1);
figure();
subplot(511);
    imagesc(rts,rfs,log10(rys.data)');
    axis('xy');
    colormap('jet');
    caxis([4,6.5])
subplot(512);
    imagesc(lts,lfs,log10(lys.data)');
    axis('xy');
    colormap('jet');
    caxis([2,6])
subplot(513);
    hold('on');
    plot([1:size(derr,1)]./sampleRate,derf)
    plot([1:size(derr,1)]./sampleRate,derr)
    Lines([],0,'k')
    Lines([],100,'k')
subplot(514);
    hold('on');
    plot([1:size(hmmState,1)]./sampleRate,hmmState)
    plot([1:size(hmmState,1)]./sampleRate,thmmState)    
    plot([1:size(hmmState,1)]./sampleRate,rfet.data)
    %plot([1:size(hmmState,1)]./sampleRate,fhmmState)
    ylim([0.5,hmm.K+0.5]);
    grid('on');
subplot(515);
    hold('on');
    plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

nfet = copy(afet);
nfet.data = [lpfet(:,1),hfet(:,1),trat(:,1) ifet(:,1)];
nfetMeanW = median(nfet(Trial.stc{'w'},:));
nfet.data = bsxfun(@minus,nfet.data,nfetMeanW);
nfetLims = [-1,1.5; -1.5,1.5; -1.5,1.25; -1.5,1.5; -1.5,1.5];


figure
%nind = thmmState==1;
%nind = thmmState==1 & stcm(:,3)==3;
nind = thmmState==1 & stcm(:,5)==5;
for a = 1:size(nfet,2)
    for b = a+1:size(nfet,2)
        subplot2(size(nfet,2),size(nfet,2),a,b);
        hist2([nfet(nind,a),nfet(nind,b)],linspace([nfetLims(a,:),30]),linspace([nfetLims(b,:),30]));
        caxis([0,2000]);        
    end
end

clear('hmm2');
updateOM = 1;
hmm2.K = 3;
hmm2 = hmminit(nfet(thmmState==1,:),hmm2,'full');
hmm2.train.cyc = 100;
hmm2.obsmodel='Gauss';
hmm2.train.obsupdate=ones([1,hmm2.K])*updateOM;
hmm2.train.init = 1;
hmm2 = hmmtrain(nfet.data,size(nfet,1),hmm2);

diag(hmm2.P)

% COMPUTE hmm states
[decode2] = hmmdecode(nfet.data,size(nfet,1),hmm2);
decode2.q_star = decode2.q_star';

%hmmState2 = decode2.q_star+hmm2.K;
hmmState2 = decode2.q_star;
hmmState2(thmmState~=1) = 0;




figure
subplot(311);
ind = hmmState==1 & stcm(:,5)==5;
bar(accumarray(decode2.q_star(ind),ones([sum(ind),1]),[hmm.K,1],@sum))
subplot(312);
ind = hmmState==1 &(stcm(:,2)==2|stcm(:,3)==3|stcm(:,4)==4);
bar(accumarray(decode2.q_star(ind),ones([sum(ind),1]),[hmm.K,1],@sum))
subplot(313);
ind = hmmState==1 &(stcm(:,6)==6);
bar(accumarray(decode2.q_star(ind),ones([sum(ind),1]),[hmm.K,1],@sum))


derr = nan([size(xyz,1),1]);
derr(dc.ind(dc.ucnt>2)) = sqrt(sum(dc.esax(dc.ucnt>2,:).^2,2));
derf = nan([size(xyz,1),1]);
derf(dc.ind(dc.ucnt>2)) = dc.ecom(dc.ucnt>2,1);
figure();
subplot(511);
    imagesc(rts,rfs,log10(rys.data)');
    axis('xy');
    colormap('jet');
    caxis([4.25,6.25])
subplot(512);
    imagesc(lts,lfs,log10(lys.data)');
    axis('xy');
    colormap('jet');
    caxis([3,5])
subplot(513);
    hold('on');
    plot([1:size(derr,1)]./sampleRate,derf)
    plot([1:size(derr,1)]./sampleRate,derr)
    Lines([],0,'k')
    Lines([],100,'k')
subplot(514);
    hold('on');
    plot([1:size(hmmState,1)]./sampleRate,hmmState)
    plot([1:size(hmmState,1)]./sampleRate,hmmState2)    
    plot([1:size(hmmState,1)]./sampleRate,rfet.data)
    %plot([1:size(hmmState,1)]./sampleRate,fhmmState)
    ylim([0.5,hmm.K+0.5]);
    grid('on');
subplot(515);
    hold('on');
    plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

derr = zeros([size(xyz,1),1]);
derr(dc.ind(dc.ucnt>2)) = sqrt(sum(dc.esax(dc.ucnt>2,:).^2,2));

figure();
subplot(411)
    ind = thmmState==1&hmmState2==3 & (stcm(:,3)==3|stcm(:,4)==4);
    histogram(derr(ind),linspace(0,600,100),'Normalization','probability');
subplot(412)
    ind = thmmState==1&hmmState2==2 & (stcm(:,3)==3|stcm(:,4)==4);
    histogram(derr(ind),linspace(0,600,100),'Normalization','probability');
subplot(413)
    ind = thmmState==1&hmmState2==1 & (stcm(:,3)==3|stcm(:,4)==4);
    histogram(derr(ind),linspace(0,600,100),'Normalization','probability');
subplot(414)
    ind = thmmState==1&hmmState2==2 & stcm(:,5)==5;
    histogram(derr(ind),linspace(0,600,100),'Normalization','probability');
ForAllSubplots('ylim([0,0.4])');


figure,
subplot(141);
ind = stcm(:,1)==1&stcm(:,3)==3;
hist2([tdfet(ind),tdrfet(ind)],linspace(-0.05,0.2,40),linspace(-0.05,0.2,40));
subplot(142);
ind = stcm(:,1)==1&stcm(:,4)==4;
hist2([tdfet(ind),tdrfet(ind)],linspace(-0.05,0.2,40),linspace(-0.05,0.2,40));
subplot(143);
ind = stcm(:,1)==1&stcm(:,5)==5;
hist2([tdfet(ind),tdrfet(ind)],linspace(-0.05,0.2,40),linspace(-0.05,0.2,40));
subplot(144);
ind = stcm(:,1)~=1&stcm(:,5)==5&hmmState==3;
hist2([tdfet(ind),tdrfet(ind)],linspace(-0.05,0.2,40),linspace(-0.05,0.2,40));


figure,
subplot(141);
ind = stcm(:,1)==1&stcm(:,3)==3;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);
subplot(142);
ind = stcm(:,1)==1&stcm(:,4)==4;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);
subplot(143);
ind = stcm(:,1)==1&stcm(:,5)==5;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);
subplot(144);
ind = stcm(:,1)~=1&stcm(:,5)==5&hmmState==3;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);

figure,
subplot(141);
ind = stcm(:,1)==1&stcm(:,3)==3&hmmState==4;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);
subplot(142);
ind = stcm(:,1)==1&stcm(:,4)==4&hmmState==4;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);
subplot(143);
ind = stcm(:,1)==1&stcm(:,5)==5&hmmState==4;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);
subplot(144);
ind = stcm(:,1)~=1&stcm(:,5)==5&hmmState==4;
hist2([afet(ind,[5,6])],linspace(-1,1,40),linspace(-1,1,40));
line([-1.5,1.5],[-1.5,1.5],'Color','w');
xlim([-1,1]);ylim([-1,1]);






clear('xcomp','ycomp');
% $$$ xcomp.data = rfet(:,1);
% $$$ xcomp.edgs = linspace(3.5,6.5,30);
% $$$ xcomp.data = dfet(:,1);
% $$$ xcomp.edgs = linspace(5,8,40);
xcomp.data = lpfet(:,1);
xcomp.edgs = linspace(5.5,7.5,40);
ycomp.data = tdfet(:,1);
ycomp.edgs = linspace(-0.1,0.2,40);
% $$$ ycomp.data = dfet(:,1);
% $$$ ycomp.edgs = linspace(5,8,40);

figure,
for grp = 1:hmm.K
% $$$     mind = dc.ucnt>=1 ...
% $$$     & (  dc.stcm(:,3)==3 ...
% $$$        | dc.stcm(:,5)==5 ...
% $$$        | dc.stcm(:,2)==2 ...
% $$$        | dc.stcm(:,4)==4 ...
% $$$        | dc.stcm(:,6)==6) ...
% $$$     & decode.q_star(dc.ind)==grp;
% $$$     mind = dc.ucnt>1 ...
% $$$     & (  dc.stcm(:,3)==3 ...
% $$$        | dc.stcm(:,5)==5) ...
% $$$     & decode.q_star(dc.ind)==grp;
% $$$     mind = dc.ucnt>=1 ...
% $$$            & dc.stcm(:,1)~=1 ...
% $$$     & (  dc.stcm(:,4)==4 ...
% $$$        | dc.stcm(:,6)==6) ...
% $$$     & decode.q_star(dc.ind)==grp;
% $$$ mind = dc.ucnt>=1 & (dc.stcm(:,1)~=1&dc.stcm(:,8)==8 ) & decode.q_star(dc.ind)==grp;
% $$$ mind = dc.ucnt>1 & (dc.stcm(:,1)==1&dc.stcm(:,8)==8 ) & decode.q_star(dc.ind)==grp;


mind = dc.ucnt>1 & (dc.stcm(:,1)~=1&dc.stcm(:,7)==7 ) & decode.q_star(dc.ind)==grp;
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,derr);
climM = [0,600];climS = [0,250];
zmean = zcomp.mean;
zmean(zcomp.count<20) = nan;
zstd  = zcomp.std;
zstd(zcomp.count<20) = nan;
zcount = zcomp.count;
zcount(zcomp.count<20) = nan;
subplot2(4,hmm.K,1,grp);
    imagesc(xcomp.ctrs,ycomp.ctrs,zmean');
    Lines([],0,'w');Lines(0,[],'w');
    axis('xy'); colormap('jet'); colorbar();caxis([climM]);
subplot2(4,hmm.K,2,grp);
    imagesc(xcomp.ctrs,ycomp.ctrs,zstd');
    Lines([],0,'w'); Lines(0,[],'w');
    axis('xy'); colormap('jet'); colorbar();caxis([climS]);
subplot2(4,hmm.K,3,grp);
    imagesc(xcomp.ctrs,ycomp.ctrs,zcount');
    Lines([],0,'w'); Lines(0,[],'w');
    axis('xy'); colormap('jet'); colorbar();
subplot2(4,hmm.K,4,1);
    ind = [Trial.stc{'w+r+n+p+s'}]; cast(ind,'TimeSeries'); resample(ind,vxy); ind = logical(ind.data);
    hist2([xcomp.data(ind),ycomp.data(ind)],xcomp.edgs,ycomp.edgs);
    grid('on'); Lines([],0,'w'); Lines(0,[],'w'); colormap('jet'); colorbar();
subplot2(4,hmm.K,4,2);
    ind = [Trial.stc{'w+r+n'}]; cast(ind,'TimeSeries'); resample(ind,vxy); ind = logical(ind.data);
    hist2([xcomp.data(ind),ycomp.data(ind)],xcomp.edgs,ycomp.edgs);
    grid('on'); Lines([],0,'w'); Lines(0,[],'w'); colormap('jet'); colorbar();
subplot2(4,hmm.K,4,3);
    ind = [Trial.stc{'s'}]; cast(ind,'TimeSeries'); resample(ind,vxy); ind = logical(ind.data);
    hist2([xcomp.data(ind),ycomp.data(ind)],xcomp.edgs,ycomp.edgs);
    grid('on'); Lines([],0,'w'); Lines(0,[],'w'); colormap('jet'); colorbar();
end



stsInds = {};
for grp = 1:hmm.K
stsInds{grp} = { dc.ucnt>1 & ( dc.stcm(:,2)==2 ) & decode.q_star(dc.ind)==grp,...
            ...
            dc.ucnt>1 & ( dc.stcm(:,3)==3  ...
                        | dc.stcm(:,5)==5) ...
           & decode.q_star(dc.ind)==grp,...
            ...
            dc.ucnt>1 & ( dc.stcm(:,4)==4  ...
                        | dc.stcm(:,6)==6) ...
           & decode.q_star(dc.ind)==grp,...
           ...
           dc.ucnt>1 & ( dc.stcm(:,1)==1 & dc.stcm(:,7)==7 ) & decode.q_star(dc.ind)==grp,...
           dc.ucnt>1 & ( dc.stcm(:,1)~=1 & dc.stcm(:,7)==7 ) & decode.q_star(dc.ind)==grp,...            
           dc.ucnt>1 & ( dc.stcm(:,1)==1 & dc.stcm(:,8)==8 ) & decode.q_star(dc.ind)==grp,...
           dc.ucnt>1 & ( dc.stcm(:,1)~=1 & dc.stcm(:,8)==8 ) & decode.q_star(dc.ind)==grp ...
};
end
stsLabels = {'Rear','Loc','Pause','GroomT','GroomN','SitT','SitN'};

figure,
for sts = 1:numel(stsInds{1})
for grp = 1:hmm.K
    mind = stsInds{grp}{sts};
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
    subplot2(numel(stsInds{1}),hmm.K,sts,grp);
    set(histogram(derr(nniz(derr)),linspace(0,600,50)),'EdgeColor','none');
    if grp==1,ylabel(stsLabels{sts});end
    if sts==1,title(['HMM Grp: ',num2str(grp)]);end
end    
end

ss = 6
sind = stcm(:,ss)==ss;
%sind = true([size(stcm,1),1]);
cmat = [sum(stcm(:,1)==1 & sind & decode.q_star==1),...
        sum(stcm(:,1)==1 & sind & decode.q_star==2),...
        sum(stcm(:,1)==1 & sind & decode.q_star==3),...
        sum(stcm(:,1)==1 & sind & decode.q_star==4),...
        sum(stcm(:,1)==1 & sind & decode.q_star==5);...
        ...
        sum(stcm(:,1)~=1 & sind & decode.q_star==1),...
        sum(stcm(:,1)~=1 & sind & decode.q_star==2),...
        sum(stcm(:,1)~=1 & sind & decode.q_star==3),...
        sum(stcm(:,1)~=1 & sind & decode.q_star==4),...
        sum(stcm(:,1)~=1 & sind & decode.q_star==5)]
bsxfun(@rdivide,cmat,sum(cmat))

sum(cmat)./sum(sum(cmat))


clear('xcomp','ycomp','zcomp','wcomp');
nbins = 15;
xcomp.data = tipmmv(:,1);
xcomp.edgs = linspace(-1,1,nbins);
ycomp.data = tipmmv(:,2);
ycomp.edgs = linspace(-.4,.4,nbins);
% $$$ zcomp.data = tdpow(:,1);
% $$$ zcomp.edgs = linspace(-5,1,nbins);
%zcomp.data = lvxy(:,2);
%zcomp.edgs = linspace(-3,2,nbins);
zcomp.data = irat(:,1);
zcomp.edgs = linspace(-1,2,nbins);

cmatY = [];
cmatN = [];
for s = 1:size(stcm,2)-1,
    for g = 1:hmm.K
        cmatY(g,s) = sum(stcm(:,1)==1 & stcm(:,s+1)==s+1 & decode.q_star==g);
        cmatN(g,s) = sum(stcm(:,1)~=1 & stcm(:,s+1)==s+1 & decode.q_star==g);
    end
end





bsxfun(@rdivide,cmat,sum(cmat))





figure
hold('on');
mind = dc.ucnt>2 & dc.stcm(:,1)==1;
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
pos = cell([1,3]);
[pos{:}] = ndgrid(xcomp.ctrs,ycomp.ctrs,zcomp.ctrs);
phl = patch(isosurface(pos{:},wcomp.count,10));
phl.FaceColor = [1,0,0];
%phl.FaceAlpha = 0.3;

mind = dc.ucnt>2 & dc.stcm(:,1)~=1;
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,10));
phl.FaceColor = [0,0,1];
phl.FaceAlpha = 0.3;


figure,imagesc(wcomp.mean(:,:,15)');axis('xy');colorbar();colormap('jet');

figure
hold('on');
% $$$ mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5 | dc.stcm(:,4)==4 |dc.stcm(:,6)==6);
% $$$ mind = dc.ucnt>3 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5 );
mind = dc.ucnt>3 & (dc.stcm(:,4)==4 |dc.stcm(:,6)==6);
derr = zeros(size(xcomp.data));
%derr(dc.ind(mind)) = dc.ucnt(mind);
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
derr(dc.ind(mind)) = dc.ecom(mind,1);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
pos = cell([1,3]);
[pos{:}] = ndgrid(xcomp.ctrs,ycomp.ctrs,zcomp.ctrs);
phl = patch(isosurface(pos{:},wcomp.count,30));
phl.FaceColor = [1,0,0];

mind = dc.ucnt>3 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5 );
derr = zeros(size(xcomp.data));
%derr(dc.ind(mind)) = dc.ucnt(mind);
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
derr(dc.ind(mind)) = dc.ecom(mind,1);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
pos = cell([1,3]);
[pos{:}] = ndgrid(xcomp.ctrs,ycomp.ctrs,zcomp.ctrs);
phl = patch(isosurface(pos{:},wcomp.count,30));
phl.FaceColor = [1,0,1];

mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)==1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomvp.count,30));
phl.FaceColor = [0,0,1];
phl.FaceAlpha = 0.3;

mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,30));
phl.FaceColor = [0,1,0];
phl.FaceAlpha = 0.3;
% $$$ phl = patch(isosurface(pos{:},wcomp.count,50));
% $$$ phl.FaceColor = [0,1,0];
% $$$ phl.FaceAlpha = 0.3;



clear('xcomp','ycomp','zcomp','wcomp');
nbins = 15;
xcomp.data = lpfet(:,1);
xcomp.edgs = linspace(3.5,7.5,nbins);
% $$$ xcomp.data = mlhfet(:,1);
% $$$ xcomp.edgs = linspace(-2,2,nbins);
% $$$ xcomp.data = mlhfet(:,2);
% $$$ xcomp.edgs = linspace(-1.5,1.5,nbins);
% $$$ ycomp.data = mlhfet(:,2);
% $$$ ycomp.edgs = linspace(-1.5,1.5,nbins);
ycomp.data = rfet(:,1);
ycomp.edgs = linspace(2.5,7,nbins);
% $$$ zcomp.data = lvxy(:,2);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
% $$$ zcomp.data = log10(sqrt(sfet(:,1).^2)).*cos(pi/6)+log10(sqrt(dsfet(:,1).^2)).*sin(pi/6);
% $$$ zcomp.edgs = linspace(-5,2,nbins);
% $$$ zcomp.data = rfet(:,1);
% $$$ zcomp.edgs = linspace(2.5,6,nbins);
zcomp.data = tdfet(:,1);
zcomp.edgs = linspace(-0.1,0.2,nbins);
% $$$ zcomp.data = irat(:,1);
% $$$ zcomp.edgs = linspace(-3,3,nbins);
% $$$ zcomp.data = tfet(:,1);
% $$$ zcomp.edgs = linspace(-5.5,8,nbins);
% $$$ zcomp.data = frat(:,1);
% $$$ zcomp.edgs = linspace(-3,3,nbins);
% $$$ zcomp.data = vrat(:,1);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
% $$$ zcomp.data = trat(:,1);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
% $$$ zcomp.data = nrat(:,1);
% $$$ zcomp.edgs = linspace(-2,2,nbins);
% $$$ zcomp.data = prat(:,1);
% $$$ zcomp.edgs = linspace(0.5,6,nbins);

figure();
hold('on');
ind = stc{'s'};
plot3(sdfet(ind),rfet(ind,1),tdfet(ind,1),'.')
ind = stc{'w+n+r'};
plot3(sdfet(ind),rfet(ind,1),tdfet(ind,1),'.')

figure();
hold('on');
ind = stc{'s'};
plot(rfet(ind),mlhfet(ind,2),'.')
ind = stc{'w+n+r'};
plot(rfet(ind),mlhfet(ind,2),'.')


isothresh = [500];
isothresh = [300];
isothresh = [250];
isothresh = [200];
isothresh = [100];
isothresh = [50];
isothresh = [20];

figure();
hold('on');
mvec = [];

mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
[pos{:}] = ndgrid(xcomp.ctrs,ycomp.ctrs,zcomp.ctrs);
%phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.std));
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
phl.FaceAlpha = 0.5;
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
caxis([0,400])
%caxis([0,16])
mvec(end+1,:) = [mean(xcomp.data(nniz(derr))),...
                 mean(ycomp.data(nniz(derr))),...
                 mean(zcomp.data(nniz(derr)))];
scatter3(mvec(end,1),...
         mvec(end,2),...
         mvec(end,3),...
         100,...
         'k',...
         'Filled');



mind = dc.ucnt>1 ...
       & dc.stcm(:,1)==1 ...
       & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2 );
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
phl.FaceAlpha = 0.5;
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
mvec(end+1,:) = [mean(xcomp.data(nniz(derr))),...
                 mean(ycomp.data(nniz(derr))),...
                 mean(zcomp.data(nniz(derr)))];
scatter3(mvec(end,1),...
         mvec(end,2),...
         mvec(end,3),...
         100,...
         'k',...
         'Filled');

bvec = diff(mvec(1:2,:));
bvec = bvec./sqrt(sum(bvec.^2,2));

tfet = sum(bsxfun(@times,bsxfun(@minus,[xcomp.data,ycomp.data,zcomp.data],mvec(2,:)),bvec),2);
t2fet = tfet;

figure,
subplot(211);
    mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2 );
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
    hist2([tfet(nniz(derr)),t2fet(nniz(derr))],linspace(-2,3,30),linspace(-1,0.6,30));
    Lines(-0.5,[],'k');Lines([],-0.25,'k');
    caxis([0,1000])
subplot(212);
    mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
    hist2([tfet(nniz(derr)),t2fet(nniz(derr))],linspace(-2,3,30),linspace(-1,0.6,30));
    Lines(-0.5,[],'k');Lines([],-0.25,'k');
    caxis([0,1000])

mind = dc.ucnt>1 ...
       & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2|dc.stcm(:,8)==8);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));

cmat = [sum(stcm(:,1)==1 & (tfet(:)>-0.3&t2fet(:)>-0.75)),...
        sum(stcm(:,1)~=1 & (tfet(:)>-0.3&t2fet(:)>-0.75));...
        sum(stcm(:,1)==1 & (tfet(:)<-0.3|t2fet(:)<-0.75)),...
        sum(stcm(:,1)~=1 & (tfet(:)<-0.3|t2fet(:)<-0.75))]
bsxfun(@rdivide,cmat,sum(cmat))
    
figure();
    mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2|dc.stcm(:,8)==8);
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
    hist2([tfet(nniz(derr)),t2fet(nniz(derr))],linspace(-2,3,30),linspace(-1,0.6,30));

    

figure();
mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2 );
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
hhdl = histogram(tfet(nniz(derr)),linspace(-2.5,4,100));
hhdl.EdgeColor = 'none';
hold('on');
mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
hhdl = histogram(tfet(nniz(derr)),linspace(-2.5,4,100));
hhdl.EdgeColor = 'none';
hhdl.FaceAlpha = 0.5;

figure();
subplot(211);
plot([1:size(tfet,1)]./sampleRate,tfet);
subplot(212);
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


clear('xcomp','ycomp','zcomp','wcomp');
nbins = 15;

sdfet = sfet.copy();
sdfet.data = log10(sqrt(sfet(:,1).^2)).*cos(pi/6)+log10(sqrt(dsfet(:,1).^2)).*sin(pi/6);

afet = [tdfet(:,1),lpfet(:,1)];
%afet = [mlhfet(:,1),mlhfet(:,2),tdfet(:,1),irat(:,1)];
%afet = [tdfet(:,1),sdfet(:,1),irat(:,1)];


mvec = [];
mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
mvec(end+1,:)= mean(afet(nniz(derr),:));
mind = dc.ucnt>1 ...
       & dc.stcm(:,1)==1 ...
       & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2 );
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
mvec(end+1,:)= mean(afet(nniz(derr),:));


wfet = bsxfun(@minus,afet,mvec(2,:));
ind = stcm(:,1)==1 & (stcm(:,2)==2 | stcm(:,3)==3 | stcm(:,4)==4);
wsigma = cov(wfet(ind,:));
wscore = copy(xyz);
wscore.data = exp(-1/2*(multiprod(wfet,multiprod((wsigma.*40).^-1,wfet,[1,2],2),2,2)))...
          ./sqrt((2*pi).^size(wfet,2).*det(wsigma));

ssfet = bsxfun(@minus,afet,mvec(1,:));
ind = stcm(:,1)~=1 & (stcm(:,5)==5);
ssigma = cov(wfet(ind,:));
sscore = copy(xyz);
sscore.data = exp(-1/2*(multiprod(ssfet,multiprod((ssigma.*40).^-1,ssfet,[1,2],2),2,2)))...
          ./sqrt((2*pi).^size(ssfet,2).*det(ssigma));


figure,
hold('on');
ind = stcm(:,1)~=1 & (stcm(:,5)==5);
plot(log10(sscore(ind)),log10(wscore(ind)),'.')
ind = stcm(:,1)==1 & (stcm(:,2)==2 | stcm(:,3)==3 | stcm(:,4)==4);
plot(log10(sscore(ind)),log10(wscore(ind)),'.')

figure,
subplot(211);
    ind = stcm(:,1)~=1 & (stcm(:,5)==5);
    hist2([log10(sscore(ind)),log10(wscore(ind))],linspace(-20,10,30),linspace(-100,150,30));
    caxis([0,1000]);
    Lines([],25,'k');Lines(0,[],'k');    
subplot(212);
    ind = stcm(:,1)==1 & (stcm(:,2)==2 | stcm(:,3)==3 | stcm(:,4)==4);
    hist2([log10(sscore(ind)),log10(wscore(ind))],linspace(-20,10,30),linspace(-100,150,30));
    caxis([0,1000]);
    Lines([],25,'k');Lines(0,[],'k');

arot = 0.05;
nsscore = log10(sscore(:,1)).*cos(arot)+log10(wscore(:,1)).*sin(arot);   
nwscore = log10(sscore(:,1)).*-sin(arot)+log10(wscore(:,1)).*cos(arot);   


figure,
hold('on');
ind = stcm(:,1)~=1 & (stcm(:,5)==5);
plot(nsscore(ind),nwscore(ind),'.')
ind = stcm(:,1)==1 & (stcm(:,5)==5);
plot(nsscore(ind),nwscore(ind),'.g')
ind = stcm(:,1)==1 & (stcm(:,2)==2 | stcm(:,3)==3 | stcm(:,4)==4);
plot(nsscore(ind),nwscore(ind),'.r')
ind = stcm(:,1)~=1 & (stcm(:,2)==2 | stcm(:,3)==3 | stcm(:,4)==4);
plot(nsscore(ind),nwscore(ind),'.m')

wthresh = -50;
sthresh = 1;
nind = nniz(xyz);
cmat = [sum(stcm(nind,1)==1&nsscore(nind)<sthresh&nwscore(nind)>wthresh),...
        sum(stcm(nind,1)~=1&nsscore(nind)<sthresh&nwscore(nind)>wthresh);...
        sum(stcm(nind,1)==1&(nsscore(nind)>sthresh|nwscore(nind)<wthresh)),...
        sum(stcm(nind,1)~=1&(nsscore(nind)>sthresh|nwscore(nind)<wthresh))]
bsxfun(@rdivide,cmat,sum(cmat))



    
wthresh = 1e4;
nind = nniz(xyz);
cmat = [sum(stcm(nind,1)==1&wscore(nind)<wthresh),sum(stcm(nind,1)~=1&wscore(nind)<wthresh);...
        sum(stcm(nind,1)==1&wscore(nind)>wthresh),sum(stcm(nind,1)~=1&wscore(nind)>wthresh)]
bsxfun(@rdivide,cmat,sum(cmat))


derr = zeros(size(xcomp.data));
derr(dc.ind) = sqrt(sum(dc.ecom(:,:).^2,2));

figure();
subplot(211);
hold('on');
plot([1:size(wscore,1)]./sampleRate,double(nsscore(:)<sthresh&nwscore(:)>wthresh).*600);
plot([1:size(wscore,1)]./sampleRate,derr,'r');
%plot([1:size(wscore,1)]./sampleRate,log10(wscore.data));
subplot(212);
hold('on');
plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
ylim([1,9]);
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


nsscore(:)<sthresh&nwscore(:)>wthresh

figure();
mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2 );
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
hhdl = histogram(log10(wscore(nniz(derr))),linspace(-200,200,100));
hhdl.EdgeColor = 'none';
hold('on');
mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
hhdl = histogram(log10(wscore(nniz(derr))),linspace(-200,200,100));
hhdl.EdgeColor = 'none';
hhdl.FaceAlpha = 0.5;



mind = dc.ucnt>1 & (dc.stcm(:,8)==8&dc.stcm(:,1)==1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
phl.FaceAlpha = 0.5;

mind = dc.ucnt>1 & (dc.stcm(:,6)==6);
derr = zeros(size(xcomp.data));
%derr(dc.ind(mind)) = dc.ecom(mind,1);
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
%phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.std));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
%phl.FaceAlpha = 0.5;


mind = dc.ucnt>2 & (dc.stcm(:,4)==4);
derr = zeros(size(xcomp.data));
%derr(dc.ind(mind)) = dc.ecom(mind,1);
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
%phl.FaceAlpha = 0.5;

mind = dc.ucnt>2 & (dc.stcm(:,2)==2);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
%phl.FaceAlpha = 0.5;

%mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5 );
mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5|dc.stcm(:,4)==4 |dc.stcm(:,6)==6|dc.stcm(:,2)==2 );
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ucnt(mind);
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
phl.FaceColor = 'interp';
phl.EdgeColor = 'none';
colormap('jet');
sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
%phl.FaceAlpha = 0.5;



%norm = 'yprob';
norm = '';
figure
subplot(311);
ind = stcm(:,1)==1|stcm(:,5)==5;
hist2([log10(sfet(ind,1).^2),log10(dsfet(ind,1).^2)],linspace(-5,1.5,30),linspace(-5,1.5,30),norm);
caxis([0,200])
%caxis([0,50])
Lines([],0,'k');
Lines([],-1,'k');
Lines(0,[],'k');
subplot(312);
ind = stcm(:,1)~=1&stcm(:,5)==5;
%ind(ind) = rand([sum(ind),1])>0.86;
hist2([log10(sfet(ind,1).^2),log10(dsfet(ind,1).^2)],linspace(-5,1.5,30),linspace(-5,1.5,30),norm);
caxis([0,200])
%caxis([0,50])
Lines([],0,'k');
Lines([],-1,'k');
Lines(0,[],'k');
subplot(313);
ind = stcm(:,1)==1&stcm(:,2)==2;
hist2([log10(sfet(ind,1).^2),log10(dsfet(ind,1).^2)],linspace(-5,1.5,30),linspace(-5,1.5,30),norm);
caxis([0,200])
Lines([],0,'k');
Lines([],-1,'k');
Lines(0,[],'k');


figure
mind = dc.ucnt>1 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5 );
derr = zeros(size(sfet.data));
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
%derr(dc.ind(mind)) = dc.ecom(mind,1);
%hist2([log10(sfet(nniz(derr),1).^2),derr(nniz(derr))],linspace(-5,1.5,30),linspace(-150,150,30),'xprob');
hist2([log10(sfet(nniz(derr),1).^2),derr(nniz(derr))],linspace(-5,1.5,30),linspace(0,200,30),'xprob');
corr(log10(sfet(nniz(derr),1).^2),derr(nniz(derr)))


xlim(xcomp.edgs([1,end]));ylim(ycomp.edgs([1,end]));zlim(zcomp.edgs([1,end]));

xlim([-1.4,1.4]);ylim([-0.5,0.5]);zlim([-1,2]);

mind = dc.ucnt>2 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
mind = dc.ucnt>2 & (dc.stcm(:,8)==8&dc.stcm(:,1)==1);
mind = dc.ucnt>2 & (dc.stcm(:,3)==3&dc.stcm(:,1)==1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.esax(mind,:).^2,2));
figure,plot3(trat(dc.ind(mind)),irat(dc.ind(mind)),derr(dc.ind(mind)),'.r')
xlim([-2,1]);ylim([-1,1]);zlim([0,800]);


mind = dc.ucnt>2 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
mind = dc.ucnt>2 & (dc.stcm(:,8)==8&dc.stcm(:,1)==1);
mind = dc.ucnt>2 & (dc.stcm(:,3)==3&dc.stcm(:,1)==1);
derr = zeros(size(xcomp.data));
derr(dc.ind(mind)) = sqrt(sum(dc.esax(mind,:).^2,2));
figure,plot3(log10(abs(out(dc.ind(mind),1))),log10(abs(out(dc.ind(mind),2))),derr(dc.ind(mind)),'.r')
xlim([-2,1]);ylim([-1,1]);zlim([0,800]);

view(-169.499998031745,       2.40000023759164);

 
 
 
 



% $$$ phl = patch(isosurface(pos{:},wcomp.count,50));
% $$$ phl.FaceColor = [0,1,0];
% $$$ phl.FaceAlpha = 0.3;




%mind = dc.ucnt>2 & (dc.stcm(:,4)==4 |dc.stcm(:,6)==6);
%mind = dc.ucnt>3 & (dc.stcm(:,3)==3 |dc.stcm(:,5)==5);
mind = dc.ucnt>3 & (dc.stcm(:,1)~=1 & dc.stcm(:,8)==8);
derr = zeros(size(xcomp.data));
%derr(dc.ind(mind)) = dc.ecom(mind,1);
derr(dc.ind(mind)) = sqrt(sum(dc.ecom(mind,:).^2,2));
[xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);

figure,
wmean = wcomp.mean;
wmean(wcomp.count(:)<20) = nan;
for s = 1:nbins-1
    subplot(2,7,s);
    set(pcolor(xcomp.ctrs,ycomp.ctrs,wmean(:,:,s)'),'EdgeColor','none');
    axis('xy')
    colormap('jet');
    caxis([0,400]);
    %caxis([-40,80]);
end



swrper = [Trial.stc{'R'}];
swrper.cast('TimeSeries');
swrper.resample(xyz);

wper = [Trial.stc{'w'}];
wper.cast('TimeSeries');
wper.resample(xyz);


figure();
hold('on');
plot(mpv,mean(ufr(logical(swrper.data),:))')
plot(mpv,mean(ufr(logical(wper.data),:))')
xlim([0,2*pi]);



clear('xcomp','ycomp','zcomp','wcomp');
nbins = 15;
xcomp.data = mlhfet(:,1);
xcomp.edgs = linspace(-2,2,nbins);
ycomp.data = mlhfet(:,2);
ycomp.edgs = linspace(-1.5,1.5,nbins);
% $$$ zcomp.data = lvxy(:,2);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
zcomp.data = irat(:,1);
zcomp.edgs = linspace(-3,3,nbins);
% $$$ zcomp.data = frat(:,1);
% $$$ zcomp.edgs = linspace(-3,3,nbins);
% $$$ zcomp.data = vrat(:,1);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
% $$$ zcomp.data = trat(:,1);
% $$$ zcomp.edgs = linspace(-3,2,nbins);
%fet = dc.ecom(:,1);
fet = sqrt(sum(dc.ecom(:,:).^2,2)); climM = [0,400]; climS = [0,200];
%fet = dc.ucnt; climM = [0,16]; climS = [0,8];


isothresh = [100];
isothresh = [50];
figure();
hold('on');
sp = gobjects([0,0]);

sp(end+1) = subplot(221);
    mind = dc.ucnt>2 & (dc.stcm(:,4)==4|dc.stcm(:,6)==6)&dc.stcm(:,1)==1;
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = fet(mind,1);
    [xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
    phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
    phl.FaceColor = 'interp';
    phl.EdgeColor = 'none';
    colormap('jet');
    sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
    caxis(climM)
    colorbar();    
    xlim(xcomp.edgs([1,end]));ylim(ycomp.edgs([1,end]));zlim(zcomp.edgs([1,end]));    
sp(end+1) = subplot(223);
    mind = dc.ucnt>2 & (dc.stcm(:,4)==4|dc.stcm(:,6)==6)&dc.stcm(:,1)==1;
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = fet(mind,1);    
    [xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
    phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.std));
    phl.FaceColor = 'interp';
    phl.EdgeColor = 'none';
    colormap('jet');
    sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
    caxis(climS)
    colorbar();    
sp(end+1) = subplot(222);
    mind = dc.ucnt>2 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = fet(mind,1);        
    [xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
    phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.mean));
    phl.FaceColor = 'interp';
    phl.EdgeColor = 'none';
    colormap('jet');
    sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
    caxis(climM)
    colorbar();
sp(end+1) = subplot(224);
    mind = dc.ucnt>2 & (dc.stcm(:,8)==8&dc.stcm(:,1)~=1);
    derr = zeros(size(xcomp.data));
    derr(dc.ind(mind)) = fet(mind,1);    
    [xcomp,ycomp,zcomp,wcomp] = compute_3d_discrete_stats(xcomp,ycomp,zcomp,derr);
    phl = patch(isosurface(pos{:},wcomp.count,isothresh,wcomp.std));
    phl.FaceColor = 'interp';
    phl.EdgeColor = 'none';
    colormap('jet');
    sum(wcomp.count(wcomp.count(:)>isothresh))./sum(wcomp.count(:))
    caxis(climS)
    colorbar();    
    %linkprop(sp,{'XLim','YLim','ZLim','CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle'});
linkprop(sp,{'XLim','YLim','ZLim','View'});
    %linkprop(sp,{'View'});    
