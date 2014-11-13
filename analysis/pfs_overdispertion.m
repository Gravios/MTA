
Session = 'jg05-20120309';
Session = 'ER06-20130612';
trialName = 'all';
mazeName = 'cof';

downSampleRate = 20;



%% Load Session if Session is not already a MTASession
if ~isa(Session,'MTASession'),
    Trial = MTATrial(Session,trialName,mazeName);
end

if ~strcmp(Trial.stc.mode,'auto_wbhr'),Trial.stc.updateMode('auto_wbhr');Trial.stc.load;end


%% Get units which are of good enough quality
units = select_units(Trial,18,'pyr');
numClu = numel(units);

newSampleRate = downSampleRate;

xyz = Trial.load('xyz');
myxyz = xyz.copy;
myxyz.filter(gtwin(.1,xyz.sampleRate));
myxyz.resample(newSampleRate);
myxyz.data = sq(myxyz(:,Trial.trackingMarker,[1,2]));



myufr = Trial.ufr.copy;
myufr.create(Trial,xyz,[],units);
myufr.resample(myxyz);

smyxyz = myxyz.copy;
smyufr = myufr.copy;


states = {'theta','rear','walk','hswalk','lswalk'};
expr = {};
ufrwd ={};
pfc = {};
ethresh =5;
for s = 1:numel(states),
smyxyz.data = myxyz(Trial.stc{states{s}},:);
smyufr.data = myufr(Trial.stc{states{s}},:);


%% Get the expected ufr for each xy 
%% Substract the expected ufr from the observed
%pfc = MTAAknnpfs(Trial,units,states{s});
pfc{s} = MTAApfs(Trial,units,states{s},'binDims',[20,20],'SmoothingWeights',[1.8,1.8]);
twpmr = ones(smyxyz.size(1),numel(units));
[~,indx] = min(abs(repmat(pfc{s}.adata.bins{1}',smyxyz.size(1),1)-repmat(smyxyz(:,1),1,numel(pfc{s}.adata.bins{1}))),[],2);
[~,indy] = min(abs(repmat(pfc{s}.adata.bins{2}',smyxyz.size(1),1)-repmat(smyxyz(:,2),1,numel(pfc{s}.adata.bins{2}))),[],2);
indrm = sub2ind(pfc{s}.adata.binSizes',indx,indy);

% $$$ figure
% $$$ subplot(121)
% $$$ pfc{s}.plot(units(1),[],1);
% $$$ subplot(122)
% $$$ scatter(indx(nniz(smyufr(:,1)),1),indy(nniz(smyufr(:,1)),1),smyufr(nniz(smyufr(:,1)),1))

for unit = units,
    rateMap = pfc{s}.plot(unit);
    twpmr(:,unit==units) = rateMap(indrm)/smyxyz.sampleRate;
end

chunksz = myxyz.sampleRate*ethresh;
trim = mod(size( twpmr,1),chunksz);
wpmr = reshape(twpmr(1:end-trim,:),chunksz,[],numel(units));
wufr = reshape(smyufr.data(1:end-trim,:),chunksz,[],numel(units));


expr{s} = sq(sum(wufr))/smyxyz.sampleRate;

ufrwd{s} = (expr{s}-sq(sum(wpmr)))./sqrt(expr{s});

end

%twpmr = wpmr;
%wpmr = pfc{s}.data.rateMap(repmat(indrm,[1,numel(units)])+size(pfc{s}.data.rateMap,1)*(ones(size(indrm,1),1)*[0:numel(units)-1]));
%[wpmr(1:10,1),twpmr(1:10,1)]

%ufrwd = (smyufr.data/smyxyz.sampleRate-wpmr)./sqrt(wpmr);


unit =22;
u =find(units==unit);
for s = 1:numel(states),
var(ufrwd{s}(nniz(ufrwd{s}(:,u))&expr{s}(:,u)>ethresh,u))
end


hfig =figure,
set(hfig,'position',[360,377,1070,324])
unit = units(1);
while unit~=-1,
u =find(units==unit);
for s = 1:numel(states);
subplot2(2,numel(states),1,s),
pfc{s}.plot(unit,[],1);
title([states{s},': ',num2str(unit)])
subplot2(2,numel(states),2,s),
hist(ufrwd{s}(expr{s}(:,u)>2,u),100)
title(['var: ',num2str(var(ufrwd{s}(nniz(ufrwd{s}(:,u))&expr{s}(:,u)>ethresh,u)))])
end
unit = figure_controls(hfig,unit,units);
end


for s =1:5,var(ufrwd{s}(nniz(ufrwd{s}(:))&expr{s}(:)>ethresh)),end

