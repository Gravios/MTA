% strajcat.m
% Some method of examining the embeded covariance matricies of many
% features to discriminate between behaviors.

Trial = MTATrial('jg05-20120317');
Trial.xyz.load(Trial);
Trial.ang.load(Trial);

Trial.filter('xyz');

comb = zeros(size(Trial.xyz,1),1,3);
comb(:,1,:) = Trial.com(Trial.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'}));
comh = zeros(size(Trial.xyz,1),1,3);
comh(:,1,:) = Trial.com(Trial.model.rb({'head_back','head_left','head_front','head_right'}));

Trial = Trial.addMarker(Trial.xyz,'body_com','[0,0,255]',...
                {{'spine_lower','body_com',[0,0,255]},...
                 {'pelvis_root','body_com',[0,0,255]},...
                 {'spine_middle','body_com',[0,0,255]},...
                 {'spine_upper','body_com',[0,0,255]}},comb);
Trial = Trial.addMarker(Trial.xyz,'head_com','[0,0,0]',...
                {{'head_back','head_com',[0,0,255]},...
                 {'head_left','head_com',[0,0,255]},...
                 {'head_front','head_com',[0,0,255]},...
                 {'head_right','head_com',[0,0,255]}},comh);
%Trial = Trial.ang.create(Trial);

sfet = [];

% $$$ sfet = [circ_dist(Trial.ang(:,2,3,1),Trial.ang(:,1,2,1)),...
% $$$ %        circ_dist(Trial.ang(:,2,3,1),Trial.ang(:,1,10,1)),...
% $$$         circ_dist(Trial.ang(:,3,4,1),Trial.ang(:,2,3,1)),...
% $$$ %        circ_dist(Trial.ang(:,3,4,1),Trial.ang(:,1,10,1)),...
% $$$         circ_dist(Trial.ang(:,4,5,1),Trial.ang(:,3,4,1)),...
% $$$ %        circ_dist(Trial.ang(:,4,5,1),Trial.ang(:,1,10,1)),...
% $$$         circ_dist(Trial.ang(:,5,7,1),Trial.ang(:,4,5,1))];,
% $$$ %        circ_dist(Trial.ang(:,5,7,1),Trial.ang(:,1,11,1))];

pfet = [Trial.ang(:,1,2,2),...
        Trial.ang(:,2,3,2),...
        Trial.ang(:,3,4,2),...
        Trial.ang(:,4,5,2),...
        Trial.ang(:,5,6,2),...
        Trial.ang(:,5,7,2),...
        Trial.ang(:,5,8,2)];

%xyz = cat(2,sfet,pfet);
markers = [1,3,4,5,7,10,11];
xyz = Filter0(gausswin(31)./gausswin(31),[zeros(1,numel(markers));Trial.vel(markers)]);
xyz = cat(2,xyz,pfet);

xyzlen = size(xyz,1);
winInterval = 10;
zpad = mod(xyzlen,winInterval);


fetseg = GetSegs(xyz,1:winInterval:size(xyz,1)-zpad,32,0);
fetSR = Trial.xyz.sampleRate/winInterval;

Trial.stc.load('auto_wbhr');
wper = Trial.stc{'w'}.copy;;
wper.resample(fetSR);
wper.cast('TimeSeries');

wfet = fetseg(:,wper.data,:);
wm = sq(mean(sq(mean(wfet))));


trajCov  = zeros(size(fetseg,2),size(fetseg,3),size(fetseg,3));

for i=1:size(fetseg,2),
trajCov(i,:,:) = (sq(fetseg(:,i,:))-repmat(wm,32,1))'*(sq(fetseg(:,i,:))-repmat(wm,32,1));
end


wtrajCov =  trajCov(wper.data,:,:);

wcoffset = 2000;
wscore = zeros(4000,1000);%zeros(size(trajCov,1),1);

for i=1:size(wscore,1)
    if sum(sum(sq(trajCov(i,:,:))))==0
        wscore(i,:) = -1;
        continue,
    end
    tcis = sq(trajCov(i,:,:));
    for j = 1:size(wscore,2);
    wtc = sq(wtrajCov(j+wcoffset,:,:));
    if sum(sum(wtc))==0|sum(isinf(wtc(:)))>=1|sum(isnan(wtc(:)))>=1,
        wscore(i,:) = -1;
        continue,
    end
        wscore(i,j) = log(det((tcis+wtc)./2))-0.5*log(det(tcis*wtc));
        %wscore(i,j) = norm(tcis*wtc*tcis,'fro');
        %dcov = tcis*wtc*tcis;
        
        %deig = real(eig(dcov));
        %wscore(i,j) =log(max(deig))-log(min(deig));
    end
end

wscore(wscore==-1)=nan;

figure,plot(1./min(real(wscore),[],2).*sq(mean(fetseg(:,1:4000,1)))')
hold on,plot(1./Filter0(gausswin(11)./sum(gausswin(11)),min(real(wscore),[],2)).*sq(mean(fetseg(:,1:4000,1)))')
%hold on,plot(max(real(wscore),[],2),'r')
%hold on,plot(mean(real(wscore),2),'g')
Lines(round(Trial.stc{'w'}(:,1)./Trial.xyz.sampleRate.*fetSR),[],'k');
Lines(round(Trial.stc{'w'}(:,2)./Trial.xyz.sampleRate.*fetSR),[],'r');
hold on,plot(sq(mean(fetseg(:,1:4000,1))),'r')

figure,plot(log10(min(real(wscore),[],2).*var(real(wscore),[],2)))
Lines(round(Trial.stc{'w'}(:,1)./Trial.xyz.sampleRate.*fetSR),[],'k');
Lines(round(Trial.stc{'w'}(:,2)./Trial.xyz.sampleRate.*fetSR),[],'r');


figure,plot(log10(var(real(wscore),[],2)./median(real(wscore),2)))
hold on,plot(Filter0(gausswin(11)./sum(gausswin(11)),var(real(wscore),[],2)./log10(median(real(wscore),2))))
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r');

figure,plot(var(real(wscore),[],2))
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r');
Lines([],10^2.03,'k');

figure,plot(std(real(wscore),[],2))
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r');


slws = sort(log10(real(wscore)));
figure,plot(slws(:,1)
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r');





twtrajCov = permute(wtrajCov,[2,3,1]);
wCovScore = zeros(size(twtrajCov,3),size(twtrajCov,3));
for i=1:size(wCovScore,1)
    tcis = sq(twtrajCov(:,:,i));
    for j = 1:size(wCovScore,2);
    wtc = sq(twtrajCov(:,:,j));
         wCovScore(i,j) = norm(tcis*wtc*tcis,'fro');
    end
end

figure,imagesc(log(wCovScore)')


rtrajCov = SelectPeriods(trajCov,round(Trial.Bhv.getState('turnR').state./Trial.xyzSampleRate.*fetSR),'c',1,1);
rtrajCov = reshape(rtrajCov,size(rtrajCov,1),size(trajCov,2),size(trajCov,2));
%rCov = cov(reshape(rtrajCov,size(rtrajCov,1),[]));
rmCov =sq(mean(rtrajCov));
rsCov =sq(std(rtrajCov));
figure,imagesc(rmCov'),colorbar
figure,imagesc(rsCov'),colorbar

iwCovSet = zeros(size(wtrajCov,1),size(wtrajCov,3),size(wtrajCov,3));
for i=1:size(wtrajCov,1)
iwCovSet(i,:,:) = inv(sq(wtrajCov(i,:,:)));
end

for i=1:size(trajCov),
wscore(i) = rtrajCov(i,:)*iwCov*rtrajCov(i,:)';
end



wCov = cov(reshape(wtrajCov,size(wtrajCov,1),[]));
iwCov = inv(wCov);
rtrajCov = reshape(trajCov,size(trajCov,1),[]);

wscore = zeros(size(trajCov,1),1);
for i=1:size(trajCov),
wscore(i) = rtrajCov(i,:)*iwCov*rtrajCov(i,:)';
end


figure,imagesc(log(abs(wscore)))

figure,plot(log(abs(wscore)))
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k')
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r')

U = zeros(size(wtrajCov,1),size(wtrajCov,2),size(wtrajCov,3));
S = zeros(size(wtrajCov,1),size(wtrajCov,2),size(wtrajCov,3));
V = zeros(size(wtrajCov,1),size(wtrajCov,2),size(wtrajCov,3));

for i=12536:size(trajCov,1),
    if sum(sum(isnan(sq(trajCov(i,:,:))))),
        continue
    end
    [U(i,:,:),S(i,:,:),V(i,:,:)] = svd(sq(trajCov(i,:,:)));
end


wU = zeros(size(wtrajCov,1),size(wtrajCov,2),size(wtrajCov,3));
wS = zeros(size(wtrajCov,1),size(wtrajCov,2),size(wtrajCov,3));
wV = zeros(size(wtrajCov,1),size(wtrajCov,2),size(wtrajCov,3));

for i=1:size(wtrajCov,1),[wU(i,:,:),wS(i,:,:),wV(i,:,:)] = svd(sq(wtrajCov(i,:,:)));,end


nwcoffset = 1001;
nwscore = zeros(4000,1000);
for i=1:size(nwscore,1)
    tcis = sq(V(i,1:3,:));
    for j = 1:size(nwscore,2);
        wtc = sq(wV(j+nwcoffset,1:3,:));
        nwscore(i,j) = sum(dot(wtc,tcis));
    end
end

figure,plot(log(abs(min(nwscore))))
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k')
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r')



%% 20130417 - varscore 

fwin = 9;
Trial.xyz = reshape(Filter0(gausswin(fwin)./sum(gausswin(fwin)),Trial.xyz),size(Trial.xyz,1),size(Trial.xyz,2),size(Trial.xyz,3));

sfet = [circ_dist(Trial.ang(:,2,3,1),Trial.ang(:,1,2,1)),...
        circ_dist(Trial.ang(:,3,4,1),Trial.ang(:,2,3,1)),...
        circ_dist(Trial.ang(:,4,5,1),Trial.ang(:,3,4,1)),...
        circ_dist(Trial.ang(:,5,7,1),Trial.ang(:,4,5,1))];

afet = [Trial.ang(:,1,2,1),Trial.ang(:,2,3,1),Trial.ang(:,3,4,1)];


rbb = Trial.Model.rb({'spine_lower','pelvis_root','spine_middle'});
rbh = Trial.Model.rb({'head_back','head_left','head_front','head_right'});
comb = Filter0(gausswin(fwin)./sum(gausswin(fwin)),Trial.com(rbb));
comh = Filter0(gausswin(fwin)./sum(gausswin(fwin)),Trial.com(rbh));
comb = sqrt(sum(diff(comb(:,[1,2])).^2,2));
comh = sqrt(sum(diff(comh(:,[1,2])).^2,2));
comv = [comb,comh];
for i =1:2,
comv(:,i) = ButFilter(comv(:,i),11,2./(Trial.xyzSampleRate/2), ...
                      'low');
comv(:,i) = clip(comv(:,i),0.003,30);
end
comv = cat(1,comv(1,:),comv);


combvf = Filter0(gausswin(121)./sum(gausswin(121)),comv(:,1))*12;

markers = [1,2,3];
xyz = Trial.xyz(:,markers,[1,2]);
lwin = 121;
xyz = reshape(Filter0(gausswin(lwin)./sum(gausswin(lwin)),xyz),size(xyz,1),length(markers),size(xyz,3));


xyzlen = size(xyz,1);
winlen = 64;
zpad = mod(xyzlen,winlen);
if zpad~=0,
%% add script to pad xyz with zeros so future matrix
%% reshaping functions will work.
xyzlen = size(xyz,1);
end 
nOverlap = 8;
fetSR = (Trial.xyzSampleRate/winlen)*nOverlap;
trlen = xyzlen/winlen*nOverlap;


vtraj =[];
for i = 1:nOverlap,
tvtraj = reshape(circshift(xyz,-(i-1).*winlen/nOverlap),[],xyzlen/winlen,size(xyz,2),size(xyz,3));
tvtraj = reshape(tvtraj,size(tvtraj,1),size(tvtraj,2),[]);
vtraj(:,i:nOverlap:trlen,:) = tvtraj-repmat(tvtraj(1,:,:),winlen,1);
end


atraj =[];
for i = 1:nOverlap,
tatraj = reshape(circshift(afet,-(i-1).*winlen/nOverlap),[],xyzlen/winlen,size(afet,2));
tatraj = reshape(tatraj,size(tatraj,1),size(tatraj,2),[]);
atraj(:,i:nOverlap:trlen,:) = tatraj-repmat(tatraj(1,:,:),winlen,1);
end


straj =[];
for i = 1:nOverlap,
tstraj = reshape(circshift(sfet,-(i-1).*winlen/nOverlap),[],xyzlen/winlen,size(sfet,2));
tstraj = reshape(tstraj,size(tstraj,1),size(tstraj,2),[]);
straj(:,i:nOverlap:trlen,:) = tstraj-repmat(tstraj(1,:,:),winlen,1);
end


%nvtraj = vtraj./repmat(max(vtraj),winlen,1);

vtrajCov  = zeros(size(vtraj,2),size(vtraj,3),size(vtraj,3));
for i=1:size(vtraj,2),
vtrajCov(i,:,:) = cov(sq(traj(:,i,:)));
end


atrajVar  = zeros(size(atraj,2),size(atraj,3));
for i=1:size(atraj,2),
    for j=1:size(atraj,3),
        atrajCov(i,j) = circ_var(sq(atraj(:,i,j)));
    end
end


strajVar  = zeros(size(straj,2),size(straj,3));
for i=1:size(straj,2),
    for j=1:size(straj,3),
        strajCov(i,j) = circ_var(sq(straj(:,i,j)));
    end
end


vtrajMean  = zeros(size(vtraj,2),size(vtraj,3));
for i=1:size(vtraj,2),
vtrajMean(i,:) = mean(sq(vtraj(:,i,:)));
end

vtrajMeanD = sqrt(sum(reshape(vtrajMean,size(vtrajMean,1),3,2).^2,3));


strajVarD = sqrt(sum(strajVar.^2,3));


strajMean  = zeros(size(straj,2),size(straj,3));
for i=1:size(straj,2),
strajMean(i,:) = circ_mean(sq(straj(:,i,:)));
end

strajMeanD = sqrt(sum(reshape(strajMean,size(strajMean,1),3,2).^2,3));


atrajMean  = zeros(size(atraj,2),size(atraj,3));
for i=1:size(atraj,2),
atrajMean(i,:) = circ_mean(sq(atraj(:,i,:)));
end

atrajMeanD = sqrt(sum((atrajMean.^2,2));


atrajVarD = sqrt(sum(atrajVar.^2,2));
strajVarD = sqrt(sum(strajVar.^2,2));



varscore = zeros(size(vtrajCov,1),1);
for i = 1:size(vtrajCov,1),
varscore(i) = trace(sq(vtrajCov(i,:,:)));
end

vmv = vtrajMeanD.*vtrajVarD;

figure,plot(vtrajMeanD.*vtrajVarD)
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k')
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r')


dvs = diff(varscore);

wondvsd = mean(GetSegs(dvs,ceil(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR)-7,8));
wofdvsd = mean(GetSegs(dvs,ceil(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR)-7,8));


figure,
plot(log10(clip(Filter0(gausswin(3)./sum(gausswin(3)),varscore),0,100000))/4)
Lines([],.4,'k');

Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r');
Lines(round(Trial.Bhv.getState('turnR').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g');
Lines(round(Trial.Bhv.getState('turnR').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m');
Lines(round(Trial.Bhv.getState('turnL').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g');
Lines(round(Trial.Bhv.getState('turnL').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m');
Lines([],.2,'k');
Lines([],-.2,'k');

figure,plot(Filter0(gausswin(21)./sum(gausswin(21)),circ_mean(atrajMean,[],2)))
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r');
Lines(round(Trial.Bhv.getState('turnR').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g');
Lines(round(Trial.Bhv.getState('turnR').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m');
Lines(round(Trial.Bhv.getState('turnL').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g');
Lines(round(Trial.Bhv.getState('turnL').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m');
Lines([],.2,'k');
Lines([],-.2,'k');

figure,plot(Filter0(gausswin(31)./sum(gausswin(31)),circ_mean(atrajMean,[],2).*mean(atrajVarD,2)))





figure,plot(Filter0(gausswin(21)./sum(gausswin(21)),circ_mean(atrajMean,[],2).*mean(atrajVarD,2)))
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r');
Lines(round(Trial.Bhv.getState('turnR').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g');
Lines(round(Trial.Bhv.getState('turnR').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m');
Lines(round(Trial.Bhv.getState('turnL').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g');
Lines(round(Trial.Bhv.getState('turnL').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m');
Lines([],.01,'k');
Lines([],-.01,'k');

detscore = zeros(size(vtrajCov,1),1);
for i = 1:size(vtrajCov,1),
detscore(i) = det(sq(vtrajCov(i,:,:)));
end

figure,
plot(log10(clip(detscore,0,100000)))
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k')
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r')

%% dsvtraj
dsvtraj = abs(sq(diff(sum(vtraj))));

figure,plot(Filter0(gausswin(11)./sum(gausswin(11)),log10(sum(dsvtraj,2))))
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k')
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r')
Lines(round(Trial.Bhv.getState('rear').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g')
Lines(round(Trial.Bhv.getState('rear').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m')
Lines([],2.5,'k');


dvs = diff(Filter0(gausswin(11)./sum(gausswin(11)),log10(sum(dsvtraj,2))));

wondvsd = mean(GetSegs(dvs,ceil(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR)-7,8));
wofdvsd = mean(GetSegs(dvs,ceil(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR)-7,8));

figure
subplot(211);
hist(wondvsd,100);
sp(1) = gca;
xl(1,:) = xlim;
subplot(212);
hist(wofdvsd,100)
sp(2) = gca;
xl(1,:) = xlim;
for s =1:length(sp),
    xlim(sp(s),[min(xl(:)),max(xl(:))])
end


tdv = sum(sq(diff(sum(vtraj))),2);
plot(abs(tdv))
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k')
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r')
Lines(round(Trial.Bhv.getState('rear').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g')
Lines(round(Trial.Bhv.getState('rear').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m')

intdv = sum(GetSegs(abs(tdv),6:size(tdv,1)-5,10));
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k')
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r')
Lines(round(Trial.Bhv.getState('rear').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g')
Lines(round(Trial.Bhv.getState('rear').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m')

plot(log10(intdv))
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k')
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r')
Lines(round(Trial.Bhv.getState('rear').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g')
Lines(round(Trial.Bhv.getState('rear').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m')


something = log10(clip(Filter0(gausswin(11)./sum(gausswin(11)),abs(sum(sq(diff(sum(vtraj))),2))),0.0001,5000));

combvf = log10(abs(tdv));
figure
bhvl = 
yl = zeros(length(bhvl),2);
for i = 1:length(bhvl),
bhvdist = SelectPeriods(combvf(:,1),round(Trial.Bhv.getState(bhvl{i}).state./Trial.xyzSampleRate.*fetSR),'c',1,1);
notbhvper = SubstractRanges([1,size(combvf,1)],round(Trial.Bhv.getState(bhvl{i}).state./Trial.xyzSampleRate.*fetSR));
notbhvdist = SelectPeriods(combvf(:,1),notbhvper,'c',1,1);
edges = [-1.4:7/1000:2];
nnwd = histc(log10(abs(notbhvdist)),edges);
nwd  = histc(log10(abs(bhvdist)),edges);
subplot(length(bhvl)+1,1,i)
bar(edges,nwd,'histc')
title(bhvl{i});
ylabel('count')
xlabel(['com body log10 speed cm/s']);
yl(i,:) = ylim;
end
subplot(length(bhvl)+1,1,i+1)
bar(edges,nnwd,'histc')
ylabel('count')
xlabel(['com body log10 speed cm/s']);
ylm = max(yl(:));
ylim([0,ylm])



tdvs = sq(sum(vtraj));


wper = round(Trial.Bhv.getState('walk').state./Trial.xyzSampleRate.*fetSR);
wper = wper(find(diff(wper,1,2)>8),:);
wper(1)=1;
dtrajCov = SelectPeriods(vmv,wper+repmat([1,-0.5*nOverlap],size(wper,1),1),'c',1,1);
vmvWalkCov = cov(dtrajCov);


tws = size(vmv,1);
fpt = zeros(size(vmv,1),1);

for i = 1:size(wper,1),
fpt(wper(i,1):wper(i,2)) = 1;
end

ymfsm = log10(vmv(1:tws,:))-repmat(mean(vmv(fpt(1:size(vmv,1))==1,:)),tws,1);
%lwpri = log(sum(fpt(:,s)==1)/length(fpt(:)));
isbcm = inv(vmvWalkCov);

twi = zeros(tws,1);
for i = 1:tws,
twi(i) = -.5*log(det(vmvWalkCov))-0.5*ymfsm(i,:)*isbcm*ymfsm(i,:)';
end

figure,plot(real(twi))
figure,plot(Filter0(gausswin(9)./sum(gausswin(9)),real(twi)))
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r');
Lines(round(Trial.Bhv.getState('rear').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g');
Lines(round(Trial.Bhv.getState('rear').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m');

wonvmn = vmv(wper(:,1),:);
wofvmn = vmv(wper(:,2),:);
figure,hist(log10(wonvmn(:)),100)
figure,hist(log10(wofvmn(:)),100)


figure,plot(mean(log10(vmv(:,1:2)),2)),Lines([],2,'k');
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r');
Lines(round(Trial.Bhv.getState('turnR').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g');
Lines(round(Trial.Bhv.getState('turnR').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m');
Lines(round(Trial.Bhv.getState('turnL').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g');
Lines(round(Trial.Bhv.getState('turnL').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m');



wf = mean(log10(vmv(:,1:2)),2);
af =  Filter0(gausswin(21)./sum(gausswin(21)),circ_mean(atrajMean,[],2).*mean(atrajVarD,2));
ind = ~isnan(wf)&~isnan(af)&wf>1.5&(af>.001|af<-.001;
figure,hist2([wf(ind),clip(af(ind),-.1,.1)],50,100)
figure,plot(wf(ind),clip(af(ind),-.1,.1),'.')

wft = 1.6;
aft = 0.003;

wfl = wf>wft;
afl = af>aft|af<-aft;

wfp = ThreshCross(wfl,0.5,5);
afp = ThreshCross(afl,0.5,5);

Trial.Bhv = MTABhv(Trial,'altest20130424',1)
Trial.Bhv.States{1} = MTAState('w','walk',round((wfp+0.25*fetSR).*Trial.xyzSampleRate./fetSR))
Trial.Bhv.States{2} = MTAState('t','turn',round((afp+0.25*fetSR).*Trial.xyzSampleRate./fetSR))
Trial.Bhv.save(Trial,1);
Trial.save

wper = Trial.Bhv.getState('walk').state;
wper(end) =size(Trial.xyz,1);
xyz = Filter0(gausswin(65)./sum(gausswin(65)),sq(Trial.xyz(:,1,[1,2])));
for i=1:size(wper,1)
wpd(i) = sum(sqrt(sum(sq(diff(xyz(wper(i,1):wper(i,2),:))).^2,2)));
end
wpd = log10(wpd);

wdt = 1.85
nwper = wper(wpd>wdt,:);



for i=1:size(wper,1)
wpdm(i) = mean(sqrt(sum(sq(diff(xyz(wper(i,1):wper(i,2),:))).^2,2)));
end
wpdm = log10(wpdm);


for i=1:size(wper,1)
wpdm(i) = mean(sqrt(sum(sq(diff(xyz(wper(i,1):wper(i,2),:))).^2,2)));
end
wpdm = log10(wpdm);


Trial.Bhv = MTABhv(Trial,'altest20130424',1)
Trial.Bhv.States{1} = MTAState('w','walk',nwper);
Trial.Bhv.States{2} = MTAState('n','old_walk',round((wfp+0.25*fetSR).*Trial.xyzSampleRate./fetSR));
Trial.Bhv.States{3} = MTAState('t','turn',round((afp+0.25*fetSR).*Trial.xyzSampleRate./fetSR));
Trial.Bhv.save(Trial,1);
Trial.save


figure,plot(af)
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r');
Lines(round(Trial.Bhv.getState('turnR').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g');
Lines(round(Trial.Bhv.getState('turnR').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m');
Lines(round(Trial.Bhv.getState('turnL').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g');
Lines(round(Trial.Bhv.getState('turnL').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m');



vind = [round(1:size(Trial.xyz,1)/size(vmv,1):size(Trial.xyz,1))]';

rear_feature = abs(Trial.xyz(vind,Trial.Model.gmi('head_front'),3)-Trial.xyz(vind,Trial.Model.gmi('spine_lower'),3)).*Trial.ang(vind,Trial.Model.gmi('spine_middle'),Trial.Model.gmi('spine_upper'),2);

rifet =-1./rear_feature;
rifet(rifet>0)=1;
wfet = vmv(:,1).*rifet;
wfet(wfet<=0) = 0.00001;
wfet = log10(wfet);
wfet(wfet<=0) = 0.00001;
figure,plot(wfet)
Lines(round(Trial.Bhv.getState('walk').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'r');
Lines(round(Trial.Bhv.getState('rear').state(:,1)./Trial.xyzSampleRate.*fetSR),[],'g');
Lines(round(Trial.Bhv.getState('rear').state(:,2)./Trial.xyzSampleRate.*fetSR),[],'m');




hang = Trial.transformOrigin('head_back','head_front',{'head_left','head_right'});

figure,plot(hang.roll)
Lines(round(Trial.Bhv.getState('walk').state(:,1)),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)),[],'r');
Lines(round(Trial.Bhv.getState('rear').state(:,1)),[],'g');
Lines(round(Trial.Bhv.getState('rear').state(:,2)),[],'m');

figure,plot(Filter0(gausswin(121)./sum(gausswin(121)),hang.roll))
Lines(round(Trial.Bhv.getState('walk').state(:,1)),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)),[],'r');
Lines(round(Trial.Bhv.getState('rear').state(:,1)),[],'g');
Lines(round(Trial.Bhv.getState('rear').state(:,2)),[],'m');



figure,plot(diff(Filter0(gausswin(121)./sum(gausswin(121)),hang.roll)))
Lines(round(Trial.Bhv.getState('walk').state(:,1)),[],'k');
Lines(round(Trial.Bhv.getState('walk').state(:,2)),[],'r');
Lines(round(Trial.Bhv.getState('rear').state(:,1)),[],'g');
Lines(round(Trial.Bhv.getState('rear').state(:,2)),[],'m');
