function revcor(Trial)


%Trial = MTATrial('jg05-20120317',{'ang'},'all');
%Trial = MTATrial('jg05-20120311',{'ang'},'all');
%Trial = MTATrial('jg05-20120315',{'ang'},'all');
Trial = MTATrial('jg05-20120310',{'ang'},'all');
%Trial = MTATrial('jg05-20120309',{'ang'},'all');

fwin = 9;
Trial.xyz = reshape(Filter0(gausswin(fwin)./sum(gausswin(fwin)),Trial.xyz),size(Trial.xyz,1),size(Trial.xyz,2),size(Trial.xyz,3));
xyzlen = size(Trial.xyz,1);

if isempty(Trial.ang),
    Trial = Trial.load_ang(0);
end

winlen = 64;
nOverlap = 8;
trajSampleRate = (Trial.xyzSampleRate/winlen)*nOverlap;

zpad = mod(xyzlen,winlen);
if zpad~=0,
xyz = Trial.xyz(1:end-zpad,Trial.Model.gmi({'spine_lower','pelvis_root','spine_middle','head_back'}),[1,2]);
ang = Trial.ang(1:end-zpad,:,:,:);
else 
xyz = Trial.xyz(:,Trial.Model.gmi({'spine_lower','pelvis_root','spine_middle','head_back'}),[1,2]);
ang = Trial.ang;
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
atrajVar  = zeros(size(atraj,2),size(atraj,3));
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




Trial = Trial.load_ufr(trajSampleRate,length(wf));
Trial.ufr = (Trial.ufr+min(Trial.ufr(:))).*20000;


figure,plot(Trial.ufr(:,11))
Lines(Trial.Bhv.getState('rear').state(:,1),[],'k');
Lines(Trial.Bhv.getState('rear').state(:,2),[],'r');


unit = 1;
while unit~=-1,
idx = LocalMinima(-Trial.ufr(:,unit),1,1);
ufet = [];
ufet(:,1) = wf(idx);
ufet(:,2) = af(idx)';
ufet(:,3) = sf(idx)';
ufet(:,4) = bf(idx)';
ufet(:,5) = hf(idx)';
clf
fcom = size(ufet,2)-1;
r = 1;
o = 1:fcom+1;
o(r)=[];
for i = 1:fcom,
subplot2(fcom,1,i,1);
plot(ufet(:,r),ufet(:,o(i)),'.')
xlim([-6,7])
ylim([-.3,.3])
end
title(num2str(unit));
unit = figure_controls(gcf,unit);
end


ucov = [];
for unit = 1:size(Trial.ufr,2),
idx = LocalMinima(-Trial.ufr(:,unit),1,1);
for shift = [-fliplr(1:12),0 1:12,]
idx(find(idx<=shift|idx-shift>length(wf)))=[];
ufet = [];
ufet(:,1) = wf(idx-shift);
ufet(:,2) = af(idx-shift)';
ufet(:,3) = sf(idx-shift)';
ufet(:,4) = bf(idx-shift)';
ufet(:,5) = hf(idx-shift)';
ufet = ufet(~isnan(ufet(:,2)),:);
ucov(:,:,unit,shift+13) = cov(ufet);
end
end

unit = 1;
d = 1;
while unit~=-1,
clf
plot(sq(ucov(d,2:end,unit,:))')
title(num2str(unit));
unit = figure_controls(gcf,unit);
end
