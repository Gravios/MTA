function Trial = label_behavior_with_heuristics(Trial,varargin)
% function Trial = labelBhv(Trial,Stc,varargin)
% Basic Set of heuristics for the most basic behavioral states
% 
% Note: Requires an update to include MoCap sampling rates other than 120Hz
[Stc,winlen,nOverlap,nwinlen,nnOverlap] = DefaultArgs(varargin,{[],64,8,16,2});

if isempty(Stc),
    Stc = Trial.stc.copy;
    Stc.updateMode('auto_wbhr');
    Stc.states = {};
end


xyz = Trial.xyz.copy;
xyzSampleRate = Trial.xyz.sampleRate;

Stc.states = {};
txyz = Trial.xyz.copy;
if txyz.isempty,txyz.load(Trial); end    
if txyz.sampleRate>120,txyz.resample(120);end
tang = create(MTADang,Trial,txyz);

txyz.filter('ButFilter',3,50,'low');
xyzlen = size(txyz,1);

trajSampleRate = (txyz.sampleRate/winlen)*nOverlap;

zpad = mod(xyzlen,winlen);
if zpad~=0,
xyz = txyz(1:end-zpad,{'spine_lower','pelvis_root','spine_middle','head_back'},[1,2]);
ang = tang(1:end-zpad,:,:,:);
else 
xyz = txyz(:,{'spine_lower','pelvis_root','spine_middle','head_back'},[1,2]);
ang = tang.data;
end 


xyzlen = size(xyz,1);
trlen = xyzlen/winlen*nOverlap;


sfet = [circ_dist(ang(:,2,3,1),ang(:,1,2,1)),...
        circ_dist(ang(:,3,4,1),ang(:,2,3,1)),...
        circ_dist(ang(:,4,5,1),ang(:,3,4,1)),...
        circ_dist(ang(:,5,7,1),ang(:,4,5,1))];
bfet = [ang(:,1,2,1),ang(:,2,3,1),ang(:,3,4,1),ang(:,4,5,1),ang(:,5,7,1)];
afet = [ang(:,1,2,1),ang(:,2,3,1),ang(:,3,4,1)];
hfet = [ang(:,3,4,1),ang(:,5,7,1)];
rfet = rear(Trial,'fet');
rfet = rfet(1:end-zpad);


vtraj =[];
for i = 1:nOverlap,
tvtraj = reshape(circshift(xyz,-(i-1).*winlen/nOverlap),[],xyzlen/winlen,size(xyz,2),size(xyz,3));
tvtraj = reshape(tvtraj,size(tvtraj,1),size(tvtraj,2),[]);
vtraj(:,i:nOverlap:trlen,:) = tvtraj-repmat(tvtraj(1,:,:),winlen,1);
end
vtrajCov  = zeros(size(vtraj,2),size(vtraj,3),size(vtraj,3));
for i=1:size(vtraj,2),
vtrajCov(i,:,:) = cov(sq(vtraj(:,i,:)));
end
vtrajVar =  zeros(size(vtraj,2),size(vtraj,3));
 for i=1:size(vtrajCov,1),
vtrajVar(i,:) = diag(sq(vtrajCov(i,:,:)));
end
vtrajVarD = sqrt(sum(reshape(vtrajVar,size(vtrajVar,1),size(xyz,2),2).^2,3));
vtrajMean  = zeros(size(vtraj,2),size(vtraj,3));
for i=1:size(vtraj,2),
   vtrajMean(i,:) = mean(sq(vtraj(:,i,:)));
end
vtrajMeanD = sqrt(sum(reshape(vtrajMean,size(vtrajMean,1),size(xyz,2),2).^2,3));


straj =[];
for i = 1:nOverlap,
tstraj = reshape(circshift(sfet,-(i-1).*winlen/nOverlap),[],xyzlen/winlen,size(sfet,2));
tstraj = reshape(tstraj,size(tstraj,1),size(tstraj,2),[]);
straj(:,i:nOverlap:trlen,:) = tstraj-repmat(tstraj(1,:,:),winlen,1);
end
strajVar  = zeros(size(straj,2),size(straj,3));
for i=1:size(straj,2),
    for j=1:size(straj,3),
        strajVar(i,j) = circ_var(sq(straj(:,i,j)));
    end
end
strajVarD = sqrt(sum(strajVar.^2,2));
strajMean  = zeros(size(straj,2),size(straj,3));
for i=1:size(straj,2),
    strajMean(i,:) = circ_mean(sq(straj(:,i,:)));
end
strajMeanD = sqrt(sum(strajMean.^2,2));


atraj =[];
for i = 1:nOverlap,
tatraj = reshape(circshift(afet,-(i-1).*winlen/nOverlap),[],xyzlen/winlen,size(afet,2));
tatraj = reshape(tatraj,size(tatraj,1),size(tatraj,2),[]);
atraj(:,i:nOverlap:trlen,:) = tatraj-repmat(tatraj(1,:,:),winlen,1);
end
atrajVar  = zeros([size(atraj,2),size(atraj,3)]);
for i=1:size(atraj,2),
    for j=1:size(atraj,3),
        atrajVar(i,j) = circ_var(sq(atraj(:,i,j)));
    end
end
atrajVarD = sqrt(sum(atrajVar.^2,2));
atrajMean  = zeros(size(atraj,2),size(atraj,3));
for i=1:size(atraj,2),
    atrajMean(i,:) = circ_mean(sq(atraj(:,i,:)));
end
atrajMeanD = sqrt(sum(atrajMean.^2,2));


htraj =[];
for i = 1:nOverlap,
thtraj = reshape(circshift(hfet,-(i-1).*winlen/nOverlap),[],xyzlen/winlen,size(hfet,2));
thtraj = reshape(thtraj,size(thtraj,1),size(thtraj,2),[]);
htraj(:,i:nOverlap:trlen,:) = thtraj-repmat(thtraj(1,:,:),winlen,1);
end
htrajVar  = zeros(size(htraj,2),size(htraj,3));
for i=1:size(htraj,2),
    for j=1:size(htraj,3),
        htrajVar(i,j) = circ_var(sq(htraj(:,i,j)));
    end
end
htrajVarD = sqrt(sum(htrajVar.^2,2));
htrajMean  = zeros(size(htraj,2),size(htraj,3));
for i=1:size(htraj,2),
    htrajMean(i,:) = circ_mean(sq(htraj(:,i,:)));
end
htrajMeanD = sqrt(sum(htrajMean.^2,2));


btraj =[];
for i = 1:nOverlap,
tbtraj = reshape(circshift(hfet,-(i-1).*winlen/nOverlap),[],xyzlen/winlen,size(hfet,2));
tbtraj = reshape(tbtraj,size(tbtraj,1),size(tbtraj,2),[]);
btraj(:,i:nOverlap:trlen,:) = tbtraj-repmat(tbtraj(1,:,:),winlen,1);
end
btrajVar  = zeros(size(btraj,2),size(btraj,3));
for i=1:size(btraj,2),
    for j=1:size(btraj,3),
        btrajVar(i,j) = circ_var(sq(btraj(:,i,j)));
    end
end
btrajVarD = sqrt(sum(btrajVar.^2,2));
btrajMean  = zeros(size(btraj,2),size(btraj,3));
for i=1:size(btraj,2),
    btrajMean(i,:) = circ_mean(sq(btraj(:,i,:)));
end
btrajMeanD = sqrt(sum(btrajMean.^2,2));


%% BASE_FEATURES
vmv = vtrajMeanD.*vtrajVarD;
wf = mean(log10(vmv(:,1:2)),2);

af =  log10(abs(get(MTADxyz('data',circ_mean(atrajMean,[],2).*atrajVarD,'sampleRate',trajSampleRate).filter('ButFilter',3,3,'low'),'data')));

sf =  Filter0(gausswin(21)./sum(gausswin(21)),circ_mean(strajMean,[],2).*mean(strajVarD,2));
bf =  Filter0(gausswin(21)./sum(gausswin(21)),circ_mean(btrajMean,[],2).*mean(btrajVarD,2));
hf =  Filter0(gausswin(21)./sum(gausswin(21)),circ_mean(htrajMean,[],2).*mean(htrajVarD,2));
ind = ~isnan(wf)&~isnan(af)&wf>1.5&(af>.001|af<-.001);
 

% vtraj projected onto a com vector 
dtraj= xyz(round(linspace(winlen,size(xyz,1),size(vtraj,2))),:,:)-repmat(xyz(round(linspace(winlen,size(xyz,1),size(vtraj,2))),1,:),1,size(xyz,2));
mdtraj = permute(reshape(repmat(sum(dtraj.*dtraj,3),2,1),size(dtraj,1),2,size(xyz,2)),[1,3,2]);
ndtraj = dtraj./mdtraj;
nvtrajMean = reshape(vtrajMean,[],size(xyz,2),2);
mvtrajm = permute(reshape(repmat(sum(nvtrajMean.*nvtrajMean,3),2,1),size(nvtrajMean,1),2,size(xyz,2)),[1,3,2]);
nvtrajm = nvtrajMean./mvtrajm;
dvtm = dot(nvtrajm,ndtraj,3);
ndvtm = dot(dtraj,nvtrajMean,3);


%% @(BHV-BASE-FILTER,MOVEMENT)
wft = 1.6;
wfl = wf>wft;
wfp = ThreshCross(wfl,0.5,5);

%% @(BHV-BASE-FILTER,ANGULAR_TRAVEL)
aft = -1.6;
afl = af>aft;
afp = ThreshCross(afl,0.5,5);

hft = 0.003;
hfl = hf>hft|hf<-hft;
hfp = ThreshCross(hfl,0.5,5);


%% @(BHV-FILTER,DIST_TRAV_DUR_WPER)
cur_wper = round((wfp+0.25*trajSampleRate).*txyz.sampleRate./trajSampleRate);
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
cur_btper = round(wper./txyz.sampleRate.*trajSampleRate);
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
wper = SubstractRanges(wper,round((bwnbtper+0.25*trajSampleRate)./trajSampleRate.*txyz.sampleRate+repmat([0,-0.25*txyz.sampleRate],size(bwnbtper,1),1)));


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
cur_wper = round(wper./txyz.sampleRate.*trajSampleRate);
twpa=zeros(size(cur_wper,1),1);
for i=1:size(cur_wper,1)
    twpa(i) = median(sum(dvtm(cur_wper(i,1):cur_wper(i,2),[2,3]),2));
end
twdt = 0;
wper = round(cur_wper(twpa>twdt,:)./trajSampleRate.*txyz.sampleRate);


%% @(BHV-FILTER,trajectory path vector(1,end) projection onto the final body direction
%% of each trajectory for 0.5 second trajectory segments)
cur_wper = round(wper./txyz.sampleRate.*trajSampleRate);
ntwpa=zeros(size(cur_wper,1),1);
for i=1:size(cur_wper,1)
ntwpa(i) = median(sum(ndvtm(cur_wper(i,1):cur_wper(i,2),[2,3]),2))/(cur_wper(i,2)-cur_wper(i,1));
end
ntwdt = 10;
wper = round(cur_wper(ntwpa>ntwdt,:)./trajSampleRate.*txyz.sampleRate);


%% @(BHV-FILTER)
cur_wper = round(wper./txyz.sampleRate.*trajSampleRate);
vdplts = zeros(size(dtraj,1),1);
for i = 1:size(cur_wper,1),
    vdplts(cur_wper(i,1):cur_wper(i,2)) = 1;
end
vdplts(mean(ndvtm(:,[2:size(ndvtm,2)]),2)<0) = 0;
wper = round(ThreshCross(vdplts,0.5,2)./trajSampleRate.*txyz.sampleRate);


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


%nwinlen = 16;
%nnOverlap = 2;
%keyboard
ntrlen = xyzlen/nwinlen*nnOverlap;
ntrajSampleRate = (txyz.sampleRate/nwinlen)*nnOverlap;

ftraj =[];
for i = 1:nnOverlap,
tftraj = reshape(circshift(xyz,round(-(i-1).*nwinlen/nnOverlap)),[],xyzlen/nwinlen,size(xyz,2),size(xyz,3));
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
cur_wper = round(wper./txyz.sampleRate.*ntrajSampleRate);
wf_mov = zeros(size(fwf,1),1);
for i = 1:size(cur_wper,1),
    wf_mov(cur_wper(i,1):cur_wper(i,2)) = 1;
end
wf_mov_thresh=1.49;
wf_mov(fwf<wf_mov_thresh) = 0;
wper = round(ThreshCross(wf_mov,0.5,5)./ntrajSampleRate.*txyz.sampleRate);



%% @(BHV-FILTER,MAX_REAR_HEIGHT)
%rear_feature = abs(txyz(:,'head_front',3)-txyz(:,'spine_lower',3)).*tang(:,'spine_middle','spine_upper',2);
rear_feature = abs(txyz(:,'head_front',3)-txyz(:,'spine_lower',3)).*tang(:,'pelvis_root','spine_middle',2);
rear_feature(isnan(rear_feature))=0;
rearThresh = 45;
minimum_interval = 64;
rear_feature = MTADxyz('data',rear_feature,'sampleRate',txyz.sampleRate);
rearPeriods = ThreshCross(rear_feature.filter('ButFilter',3,1)>rearThresh,...
                          0.5,minimum_interval);

cur_rper = rearPeriods;%rear(Trial,'com',45); %50
mrh = zeros(size(cur_rper,1),1);
for i = 1:size(cur_rper,1),
    mrh(i) = max(txyz(cur_rper(i,1):cur_rper(i,2),Trial.model.gmi('head_front'),3));
end
rper_height_thresh = 170;
rper = cur_rper(mrh>rper_height_thresh,:);

%% @(BHV-FILTER,REMOVE_REARING_PERIODS)
if ~isempty(rper),
    wper = SubstractRanges(wper,rper+repmat([-64,24],size(rper,1),1));
end

aper = resample(Trial.sync.copy,txyz.sampleRate);
aper.data = aper.data-aper.data(1)+1;
aper = aper+[5,-5];

Stc.addState(Trial.spath,...
                   Trial.filebase,...
                   aper.data,...
                   txyz.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'gper','a');

Stc{'a'}.save(1);
Stc.addState(Trial.spath,...
                   Trial.filebase,...
                   wper,...
                   txyz.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'walk','w');

Stc.addState(Trial.spath,...
                   Trial.filebase,...
                   rper,...
                   txyz.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'rear','r');

for i= 1:numel(Stc.states),
   Stc.states{i}.data(Stc.states{i}.data(:,1)>txyz.size(1),:)=[];
   Stc.states{i}.data(end,Stc.states{i}.data(end,:)>txyz.size(1))=txyz.size(1);
   Stc.states{i}.resample(xyzSampleRate);
end






% $$$ wind = Stc{'w'}.copy;
% $$$ sind = Stc{'w'}.copy;
% $$$ ang = tang.copy;
% $$$ ang.load(Trial);
% $$$ sind.data = ThreshCross(ang(:,5,7,2)<-angThresh,.5,20);
% $$$ sind.label = 'lang';
% $$$ sind.key = 'p';
% $$$ sind.filename = [Trial.filebase '.sst.lang.p.mat'];
% $$$ sind = sind&wind.data;
% $$$ Stc.states{end+1} = sind;
% $$$ 
% $$$ 
% $$$ sind = Stc{'w'}.copy;
% $$$ sind.data = ThreshCross(ang(:,5,7,2)>-angThresh,.5,20);
% $$$ sind.label = 'hang';
% $$$ sind.key = 'c';
% $$$ sind.filename = [Trial.filebase '.sst.hang.c.mat'];
% $$$ sind = sind&wind.data;
% $$$ Stc.states{end+1} = sind;



Stc.save(1);
Trial.stc = Stc;
Trial.save;
               

