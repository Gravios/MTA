function Trial = bhv_auto(Trial,Stc)
% Basic Set of heuristics for the most basic behavioral states
% 
% Note: Requires an update to include MoCap sampling rates other than 120Hz

%% Testing Vars

Trial = MTATrial('jg05-20120317');
marks = Trial.xyz.model.gmi({'spine_lower','pelvis_root','spine_middle','head_back','head_front'});
%%%%%%%%%%%%%%%%


% Load marker positions and calculate inter marker angles
% Down-sample to 120 Hz if necessary
xyz = Trial.load('xyz');
if xyz.sampleRate>120,
    xyz.resample(120);
end
ang = create(Trial.ang.copy,Trial,xyz);

vfet = xyz.vel(marks,[1,2]);
vfet.data = ButFilter(vfet.data,3,8/(xyz.sampleRate/2));
ind = nniz(vfet);
figure,hist2(log10([mean(vfet(ind,1),2),mean(vfet(ind,4),2)]),linspace(-.5,2,100),linspace(-.5,2,100))
ind = Trial.stc{'w'};
figure,hist2(log10([mean(vfet(ind,1),2),mean(vfet(ind,4),2)]),linspace(-.5,2,100),linspace(-.5,2,100))
ind = Trial.stc{'n'};
figure,hist2(log10([mean(vfet(ind,1),2),mean(vfet(ind,4),2)]),linspace(-.5,2,100),linspace(-.5,2,100))

% Movement of markers in general body direction feature
% Get segments for .5 second periods
winds = round(.3*xyz.sampleRate);
hwin = round(winds/2);
dtraj = sq(circshift(xyz.data,-hwin)-circshift(xyz.data,hwin));
otraj = ButFilter(xyz(:,4,:),3,4/xyz.sampleRate*.5,'low')-ButFilter(xyz(:,1,:),3,4/xyz.sampleRate*.5,'low');
ndvtm = dot(dtraj,repmat(otraj,[1,size(dtraj,2),1]),3);
ndat = [ndvtm(:,1),ndvtm(:,4)];
sd = sign(ndat);
ndat(ndat<1&ndat>-1)=1;
ndat(~nniz(ndat(:))) = 1;
ndat = log10(abs(ndat)).*sd;
tfet = xyz.copy;
tfet.data = ndat*[.5;.5];

% figure,plot(tfet.data*[.5;.5])
% figure,bar(linspace(-5,5,1000),hist(tfet.data*[.5;.5],linspace(-5,5,1000)),'histc'),ylim([0,1000])
% figure,bar(linspace(-5,5,1000),hist(tfet(Trial.stc{'w'},:)*[.5;.5],linspace(-5,5,1000)),'histc'),ylim([0,1000])
% figure,bar(linspace(-5,5,1000),hist(tfet(Trial.stc{'r'},:)*[.5;.5],linspace(-5,5,1000)),'histc'),ylim([0,1000])
% figure,bar(linspace(-5,5,1000),hist(tfet(Trial.stc{'n'},:)*[.5;.5],linspace(-5,5,1000)),'histc'),ylim([0,1000])
% figure,bar(linspace(-5,5,1000),hist(tfet(Trial.stc{'m'},:)*[.5;.5],linspace(-5,5,1000)),'histc'),ylim([0,1000])
% figure,bar(linspace(-5,5,1000),hist(tfet(Trial.stc{'s'},:)*[.5;.5],linspace(-5,5,1000)),'histc'),ylim([0,1000])
% nind = Trial.stc{'w'};
% figure,hist2(tfet(nind,:),linspace(-5,5,100),linspace(-5,5,100))
% caxis([0,200])

rfet = rear(Trial,'MTA');

% Total spine curvature
sfet = ang.copy;
sfet.data = sum(abs([circ_dist(ang(:,2,3,1),ang(:,1,2,1)),...
                     circ_dist(ang(:,3,4,1),ang(:,2,3,1)),...
                     circ_dist(ang(:,4,5,1),ang(:,3,4,1)),...
                     circ_dist(ang(:,5,7,1),ang(:,4,5,1))]),2);


%% BASE_FEATURES
vmv = vtrajMeanD.*vtrajVarD;
wf = mean(log10(vmv(:,1:2)),2);
af =  Filter0(gausswin(21)./sum(gausswin(21)),circ_mean(atrajMean,[],2).*mean(atrajVarD,2));

sf =  Filter0(gausswin(21)./sum(gausswin(21)),circ_mean(strajMean,[],2).*mean(strajVarD,2));
bf =  Filter0(gausswin(21)./sum(gausswin(21)),circ_mean(btrajMean,[],2).*mean(btrajVarD,2));
hf =  Filter0(gausswin(21)./sum(gausswin(21)),circ_mean(htrajMean,[],2).*mean(htrajVarD,2));
ind = ~isnan(wf)&~isnan(af)&wf>1.5&(af>.001|af<-.001);
 


%% @(BHV-BASE-FILTER,MOVEMENT)
wft = 1.6;
wfl = wf>wft;
wfp = ThreshCross(wfl,0.5,5);

%% @(BHV-BASE-FILTER,ANGULAR_TRAVEL)
aft = 0.003;
afl = af>aft|af<-aft;
afp = ThreshCross(afl,0.5,5);

hft = 0.003;
hfl = hf>hft|hf<-hft;
hfp = ThreshCross(hfl,0.5,5);


%% @(BHV-FILTER,DIST_TRAV_DUR_WPER)
cur_wper = round((wfp+0.25*trajSampleRate).*Trial.xyz.sampleRate./trajSampleRate);
cur_wper(end) =size(xyz,1);
fxyz = Filter0(gausswin(65)./sum(gausswin(65)),sq(xyz(:,1,[1,2])));
wdist_score = zeros(size(cur_wper,1),1);
for i=1:size(cur_wper,1)
    wdist_score(i) = sum(sqrt(sum(sq(diff(fxyz(cur_wper(i,1):cur_wper(i,2),:))).^2,2)));
end
wdist_score = log10(wdist_score);
wdist_thresh = 1.85;
wper = cur_wper(wdist_score>wdist_thresh,:);
bdwper = cur_wper(wdist_score<wdist_thresh,:);


%% @(BHV-FILTER,SPINE_CURV_THRESH)
cur_wper = wper;
wang_score = zeros(size(cur_wper,1),1);
for i=1:size(cur_wper,1)
    wang_score(i) = sum(sum(sfet(cur_wper(i,1):cur_wper(i,2),:),2))/(cur_wper(i,2)-cur_wper(i,1));
end
wang_score = abs(wang_score);
wang_thresh = 0.9;
wper = cur_wper(wang_score<wang_thresh,:);
bawper = cur_wper(wang_score>wang_thresh,:);


%% @(BHV-FILTER,ANG_SUM_THRESH)
cur_btper = round(wper./Trial.xyz.sampleRate.*trajSampleRate);
btpa = zeros(size(cur_btper,1),1);
for i=1:size(cur_btper,1)
    btpa(i) = max(abs(bf(cur_btper(i,1):cur_btper(i,2))));
end
btdt = 0.02;
cur_btper = cur_btper(btpa>btdt,:);
baf = zeros(size(bf));
for i=1:size(cur_btper,1)
    baf(cur_btper(i,1):cur_btper(i,2)) = abs(bf(cur_btper(i,1):cur_btper(i,2)));
end
bwnbtper = ThreshCross(baf>.02,0.5,1);
wper = SubstractRanges(wper,round((bwnbtper+0.25*trajSampleRate)./trajSampleRate.*Trial.xyz.sampleRate+repmat([0,-0.25*Trial.xyz.sampleRate],size(bwnbtper,1),1)));


%% @(BHV-FILTER,REMOVE_SHORT_PERIODS,MIN_32_FRAMES)
wper_dur_thresh = 32;
wper = wper(diff(wper,1,2)>wper_dur_thresh,:);


%% [wper,bdwper] = @({BHV-FILTER,DIST_TRAV_DUR_WPER},wper,xyz)
cur_wper = wper;
fxyz = Filter0(gausswin(65)./sum(gausswin(65)),sq(xyz(:,1,[1,2])));
wdist_score = zeros(size(cur_wper,1),1);
for i=1:size(cur_wper,1)
    wdist_score(i) = sum(sqrt(sum(sq(diff(fxyz(cur_wper(i,1):cur_wper(i,2),:))).^2,2)));
end
wdist_score = log10(wdist_score);
wdist_thresh = 1.6;
wper = cur_wper(wdist_score>wdist_thresh,:);
bdwper = cur_wper(wdist_score<wdist_thresh,:);


%% @(BHV-FILL-GAPS,MAX_16_FRAMES)
cur_wper = wper;
dnpi = find(cur_wper(2:end,1)-cur_wper(1:end-1,2)<16);
for i = 1:length(dnpi),
newper = [cur_wper(dnpi(i),1),cur_wper(dnpi(i)+1,2)];
cur_wper(dnpi(i),:)=newper;
cur_wper(dnpi(i)+1,:)=[];
dnpi=dnpi-1;
end
wper = cur_wper;


%% @(BHV-FILTER,REMOVE_SHORT_PERIODS,MAX_32_FRAMES)
wper_dur_thresh = 32;
wper = wper(diff(wper,1,2)>wper_dur_thresh,:);


%% @(BHV-FILTER,NON_FORWARD_WALK_PERIODS)
cur_wper = round(wper./Trial.xyz.sampleRate.*trajSampleRate);
twpa=zeros(size(cur_wper,1),1);
for i=1:size(cur_wper,1)
    twpa(i) = median(sum(dvtm(cur_wper(i,1):cur_wper(i,2),[2,3]),2));
end
twdt = 0;
wper = round(cur_wper(twpa>twdt,:)./trajSampleRate.*Trial.xyz.sampleRate);


%% @(BHV-FILTER,trajectory path vector(1,end) projection onto the final body direction
%% of each trajectory for 0.5 second trajectory segments)
cur_wper = round(wper./Trial.xyz.sampleRate.*trajSampleRate);
ntwpa=zeros(size(cur_wper,1),1);
for i=1:size(cur_wper,1)
ntwpa(i) = median(sum(ndvtm(cur_wper(i,1):cur_wper(i,2),[2,3]),2))/(cur_wper(i,2)-cur_wper(i,1));
end
ntwdt = 10;
wper = round(cur_wper(ntwpa>ntwdt,:)./trajSampleRate.*Trial.xyz.sampleRate);


%% @(BHV-FILTER)
cur_wper = round(wper./Trial.xyz.sampleRate.*trajSampleRate);
vdplts = zeros(size(dtraj,1),1);
for i = 1:size(cur_wper,1),
    vdplts(cur_wper(i,1):cur_wper(i,2)) = 1;
end
vdplts(mean(ndvtm(:,[2:size(ndvtm,2)]),2)<0) = 0;
wper = round(ThreshCross(vdplts,0.5,2)./trajSampleRate.*Trial.xyz.sampleRate);


%% @(BHV-FILTER,DIST_TRAV_DUR_WPER)
cur_wper = wper;
wdist_score = [];
fxyz = Filter0(gausswin(65)./sum(gausswin(65)),sq(xyz(:,1,[1,2])));
for i=1:size(cur_wper,1)
    wdist_score(i) = sqrt(sum(sq(fxyz(cur_wper(i,2),:)-fxyz(cur_wper(i,1),:)).^2,2));
end
wdist_score = log10(wdist_score);
wdist_thresh = 1.8;
cur_wper = cur_wper(wdist_score>wdist_thresh,:);
wper = cur_wper;


%% @(BHV-FILTER,SPINE_CURV_THRESH)
cur_wper = wper;
wang_score = zeros(size(cur_wper,1),1);
for i=1:size(cur_wper,1)
    wang_score(i) = sum(sum(sfet(cur_wper(i,1):cur_wper(i,2),:),2))/(cur_wper(i,2)-cur_wper(i,1));
end
wang_score = abs(wang_score);
wang_thresh = 0.9;
wper = cur_wper(wang_score<wang_thresh,:);
bawper = cur_wper(wang_score>wang_thresh,:);


nwinlen = 16;
nnOverlap = 2;
ntrlen = xyzlen/nwinlen*nnOverlap;
ntrajSampleRate = (Trial.xyz.sampleRate/nwinlen)*nnOverlap;

ftraj =[];
for i = 1:nnOverlap,
tftraj = reshape(circshift(xyz,-(i-1).*nwinlen/nnOverlap),[],xyzlen/nwinlen,size(xyz,2),size(xyz,3));
tftraj = reshape(tftraj,size(tftraj,1),size(tftraj,2),[]);
ftraj(:,i:nnOverlap:ntrlen,:) = tftraj-repmat(tftraj(1,:,:),nwinlen,1);
end
ftrajCov  = zeros(size(ftraj,2),size(ftraj,3),size(ftraj,3));
for i=1:size(ftraj,2),
ftrajCov(i,:,:) = cov(sq(ftraj(:,i,:)));
end
ftrajVar =  zeros(size(ftraj,2),size(ftraj,3));
 for i=1:size(ftrajCov,1),
ftrajVar(i,:) = diag(sq(ftrajCov(i,:,:)));
end
ftrajVarD = sqrt(sum(reshape(ftrajVar,size(ftrajVar,1),size(xyz,2),2).^2,3));
ftrajMean  = zeros(size(ftraj,2),size(ftraj,3));
for i=1:size(ftraj,2),
   ftrajMean(i,:) = mean(sq(ftraj(:,i,:)));
end
ftrajMeanD = sqrt(sum(reshape(ftrajMean,size(ftrajMean,1),size(xyz,2),2).^2,3));

fmv = ftrajMeanD.*ftrajVarD;
fwf = Filter0(gausswin(11)./sum(gausswin(11)),mean(log10(fmv(:,1:2)),2));


%% @(BHV-FILTER,MOVEMENT)
cur_wper = round(wper./Trial.xyz.sampleRate.*ntrajSampleRate);
wf_mov = zeros(size(fwf,1),1);
for i = 1:size(cur_wper,1),
    wf_mov(cur_wper(i,1):cur_wper(i,2)) = 1;
end
wf_mov_thresh=1.49;
wf_mov(fwf<wf_mov_thresh) = 0;
wper = round(ThreshCross(wf_mov,0.5,5)./ntrajSampleRate.*Trial.xyz.sampleRate);



%% @(BHV-FILTER,MAX_REAR_HEIGHT)
cur_rper = rear(Trial,'com',45); %50
mrh = zeros(size(cur_rper,1),1);
for i = 1:size(cur_rper,1),
    mrh(i) = max(Trial.xyz(cur_rper(i,1):cur_rper(i,2),Trial.model.gmi('head_front'),3));
end
rper_height_thresh = 170;
rper = cur_rper(mrh>rper_height_thresh,:);

%% @(BHV-FILTER,REMOVE_REARING_PERIODS)
wper = SubstractRanges(wper,rper+repmat([-64,24],size(rper,1),1));


Stc.addState(Trial.spath,...
                   Trial.filebase,...
                   wper,...
                   Trial.xyz.sampleRate,...
                   Trial.xyz.sync.copy,...
                   Trial.xyz.origin,...
                   'walk','w');

Stc.addState(Trial.spath,...
                   Trial.filebase,...
                   round((afp+0.25*trajSampleRate).*Trial.xyz.sampleRate./trajSampleRate),...
                   Trial.xyz.sampleRate,...
                   Trial.xyz.sync.copy,...
                   Trial.xyz.origin,...
                   'bturn','b');

Stc.addState(Trial.spath,...
                   Trial.filebase,...
                   round((hfp+0.25*trajSampleRate).*Trial.xyz.sampleRate./trajSampleRate),...
                   Trial.xyz.sampleRate,...
                   Trial.xyz.sync.copy,...
                   Trial.xyz.origin,...
                   'hturn','h');

Stc.addState(Trial.spath,...
                   Trial.filebase,...
                   rper,...
                   Trial.xyz.sampleRate,...
                   Trial.xyz.sync.copy,...
                   Trial.xyz.origin,...
                   'rear','r');
               
Stc.save(1);
Trial.stc = Stc;
Trial.save;
               

