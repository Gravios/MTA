% asFIG-20130507.m
% investigate xy distance variance for 0.5 sec trajectories
% for walking and rearing


slist= {'jg04-20120130','jg04-20120128','jg05-20120309','jg05-20120310','jg05-20120311','jg05-20120315','jg05-20120316','jg05-20120317'}; 
wmVMVper = {[]};
rmVMVper = {[]};      

for j=1:length(slist)
s = MTASession(slist{j},{},'cof');
Trial = MTATrial(s,{},'all');

%% investigate xy distance variance for 0.5 sec trajectories
%% for walking and rearing

fwin = 9;
% add filter function Trial = Trial.filter('xyz',9,'gausswin');
Trial.xyz = reshape(Filter0(gausswin(fwin)./sum(gausswin(fwin)),Trial.xyz),size(Trial.xyz,1),size(Trial.xyz,2),size(Trial.xyz,3));

xyzlen = size(Trial.xyz,1);
winlen = 64;
nOverlap = 8;
trajSampleRate = (Trial.xyzSampleRate/winlen)*nOverlap;

zpad = mod(xyzlen,winlen);
if zpad~=0,
xyz = Trial.xyz(1:end-zpad,Trial.Model.gmi({'head_back','head_left','head_front','head_right'}),[1,2]);
else
xyz = Trial.xyz(1:end,Trial.Model.gmi({'head_back','head_left','head_front','head_right'}),[1,2]);
end

xyzlen = size(xyz,1);
trlen = xyzlen/winlen*nOverlap;

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


vmv = vtrajMeanD.*vtrajVarD;

%% FIG vmv distr
min_vmv = min(vmv(vmv(:)~=0));
vmv(vmv==0) = min_vmv;
figure,
 hist(log10(vmv),1000)
title({'Log 10 Distributions of the vmv';'jg05-20120309 head markers'})
%% FIG - END

figure
sp=[];
wper = round(Trial.Bhv.getState('walk').state./Trial.xyzSampleRate.*trajSampleRate);
wmVMVper{j} = SelectPeriods(mean(vmv,2),wper,'c',1,1);
sp(1)=subplot2(2,1,1,1);
hist(log10(wmVMVper{j}),500)

rper = round(Trial.Bhv.getState('rear').state./Trial.xyzSampleRate.*trajSampleRate);
rmVMVper{j} = SelectPeriods(mean(vmv,2),rper,'c',1,1);
sp(2)=subplot2(2,1,2,1);
hist(log10(rmVMVper{j}),500)

linkaxes(sp,'x')

end


lmw = cellfun(@log10,wmVMVper,'UniformOutput',false);
lmw = cellfun(@mean,lmw);
lmr = cellfun(@log10,rmVMVper,'UniformOutput',false);
lmr = cellfun(@mean,lmr);




%% smoothed hists
for j = 1:length(rmVMVper),
figure
rmpl = length(rmVMVper{j});
[ncount,vals] = hist(log10(rmVMVper{j}),round(rmpl/log10(rmpl)));           
gwo = round(10^(log10(rmpl)/2));
gwo = gwo+(0==mod(gwo,2));
sp(1)=subplot2(2,1,1,1);
bar(vals,Filter0(gausswin(gwo)./sum(gausswin(gwo)),ncount))
wmpl = length(wmVMVper{j});
[ncount,vals] = hist(log10(wmVMVper{j}),round(wmpl/log10(wmpl)));           
gwo = round(10^(log10(wmpl)/2));
gwo = gwo+(0==mod(gwo,2));
sp(2)=subplot2(2,1,2,1);
bar(vals,Filter0(gausswin(gwo)./sum(gausswin(gwo)),ncount))
linkaxes(sp,'x')
end