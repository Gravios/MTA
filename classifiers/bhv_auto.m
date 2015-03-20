function Trial = bhv_auto(Trial,Stc)
% Basic Set of heuristics for the most basic behavioral states
% 
% Note: Requires an update to include MoCap sampling rates other than 120Hz

%% Testing Vars

Trial = MTATrial('jg05-20120317');
marks = Trial.xyz.model.gmi({'spine_lower','pelvis_root','head_back','head_front'});
%%%%%%%%%%%%%%%%



% Load Marker Positions and down sample to 30 Hz

vfet = vel(resample(Trial.load('xyz'),8),marks,[1,2]);

xyz = resample(Trial.load('xyz'),36);
% Get segments for .5 second periods
xyzs = xyz.segs(1:xyz.size(1),9,0);
xyzs = xyzs(:,:,:,[1,2]);
% Calculate distances
dxyzs = sqrt(sum(bsxfun(@minus,xyzs,xyzs(5,:,:,:)).^2,4));

% Movement Feature 
%a bit better than the sum of the marker vels
mfet = xyz.copy;
mfet.data = circshift(log10(abs(sq(sum(prod(bsxfun(@minus,dxyzs(:,:,marks),median(dxyzs(:,:,marks))),3)))').*sq(mean(dxyzs(:,:,marks(1))))'),5);

rfet = rear(Trial,'MTA')

vfet.resample(mfet);

vvf = xyz.copy;
vvf.data = var(log10(vfet.data),[],2);
% $$$ nind = Trial.stc{'a'};
% $$$ nind = Trial.stc{'m'};
% $$$ nind = Trial.stc{'n'};
% $$$ nind = Trial.stc{'w'};
% $$$ nind = Trial.stc{'a-r'};
% $$$ nind = Trial.stc{'a-r-w'};
% $$$ figure,bar(linspace(-8,12,1000),histc(mfet(nind)-log10(vvf(nind)),linspace(-8,12,1000)),'histc')
% $$$ % $$$ nind = nniz(vfet)&nniz(vvf);
% $$$ % $$$ nind = Trial.stc{'w'};
% $$$ % $$$ nind = Trial.stc{'a-r-w'};
% $$$ figure,hist2([sum(log10(vfet(nind,:)),2),log10(vvf(nind))],linspace(-8,8,100),linspace(-5,.5,100))
% $$$ 
% $$$ figure,hist2([sum(log10(vfet(nind,:)),2),log10(xyz(nind,1,3))],linspace(-8,8,100),linspace(0,2,100))
% $$$ figure,hist2([mfet(nind,:),log10(vvf(nind))],linspace(-8,8,100),linspace(-5,.5,100))


ang = create(Trial.ang.copy,Trial,xyz);

sfet = xyz.copy;
sfet.data = [circ_dist(ang(:,2,3,1),ang(:,1,2,1)),...
        circ_dist(ang(:,3,4,1),ang(:,2,3,1)),...
        circ_dist(ang(:,4,5,1),ang(:,3,4,1)),...
        circ_dist(ang(:,5,7,1),ang(:,4,5,1))];


figure,plot(diff(sum( sfet,2)))
Lines(Trial.stc{'w',ang.sampleRate}(:),[],'k');
Lines(Trial.stc{'n',ang.sampleRate}(:),[],'g');
Lines(Trial.stc{'m',ang.sampleRate}(:),[],'m');

san = xyz.copy;
san.data = [ang(:,1,4,1),ang(:,1,3,1),ang(:,2,4,1),ang(:,4,7,1)];
san = san.segs(1:san.size(1),9,0);
san = bsxfun(@circ_dist,san,san(5,:,:,:));

afet = xyz.copy;
afet.data = circshift(log10(abs(sum(prod(bsxfun(@circ_dist,san,circ_mean(san)),3))').*abs(sq(circ_mean(san(:,:,1)))')),5);





nind = Trial.stc{'a'};
nind = Trial.stc{'m'};
nind = Trial.stc{'n'};
nind = Trial.stc{'w'};
nind = Trial.stc{'a-r'};
nind = Trial.stc{'a-r-w'};
figure,hist2([mfet(nind),afet(nind)],linspace(-8,8,100),linspace(-20,2,100))

figure,hist2(log10([vfet(nind,1),vfet(nind,4)]),linspace(-2,2,100),linspace(-2,2,100))


xyzs = xyz.segs(1:xyz.size(1),18,0);
xyzs = xyzs(:,:,:,[1,2]);

dtraj = sq(xyzs(end,:,:,:)-xyzs(1,:,:,:));
mdtraj = permute(reshape(repmat(sum(dtraj.*dtraj,3),2,1),size(dtraj,1),2,size(xyz,2)),[1,3,2]);
ndtraj = dtraj./mdtraj;
ndvtm = dot(dtraj,sq(mean(xyzs)),3);

%% START HERE
figure,plot(ndvtm(:,1:8))
Lines(Trial.stc{'w',ang.sampleRate}(:),[],'k');
Lines(Trial.stc{'n',ang.sampleRate}(:),[],'g');
Lines(Trial.stc{'m',ang.sampleRate}(:),[],'m');
Lines(Trial.stc{'r',ang.sampleRate}(:),[],'r');




%% BASE_FEATURES
vmv = vtrajMeanD.*vtrajVarD;
wf = mean(log10(vmv(:,1:2)),2);
af =  Filter0(gausswin(21)./sum(gausswin(21)),circ_mean(atrajMean,[],2).*mean(atrajVarD,2));

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
               

