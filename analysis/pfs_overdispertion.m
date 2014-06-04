
%Session = 'jg05-20120317';
pfcbhv = 'hwalk&theta';
downSampleRate = 20;


%% Load Session if Session is not already a MTASession
% $$$ if ~isa(Session,'MTASession'),
% $$$     Trial = MTATrial(Session);
% $$$     %{'Pfs',   {'walk','rear'}}},...
% $$$ end
if isempty(Trial.nq),Trial.load('nq');end
if ~strcmp(Trial.stc.mode,'auto_wbhr'),Trial.stc.updateMode('auto_wbhr');Trial.stc.load;end


%% Get units which are of good enough quality
units = select_units(Trial,25);
numClu = numel(units);

newSampleRate = downSampleRate;

myxyz = Trial.xyz.copy;
myxyz.load(Trial);
myxyz.filter(gausswin(9)./sum(gausswin(9)));
myxyz.resample(newSampleRate);
myxyz.data = sq(myxyz(:,7,[1,2]));



myufr = Trial.ufr.copy;
myufr.create(Trial,myxyz,[],units);


myxyz.data = myxyz(Trial.stc{pfcbhv},:);
myufr.data = myufr(Trial.stc{pfcbhv},:);


%% Get the expected ufr for each xy 
%% Substract the expected ufr from the observed
pfc = MTAAknnpfs(Trial,units,pfcbhv);
wpmr = ones(myxyz.size(1),numel(units));
[~,indx] = min(abs(repmat(pfc.adata.bins{1}',myxyz.size(1),1)-repmat(myxyz(:,1),1,numel(pfc.adata.bins{1}))),[],2);
[~,indy] = min(abs(repmat(pfc.adata.bins{2}',myxyz.size(1),1)-repmat(myxyz(:,2),1,numel(pfc.adata.bins{2}))),[],2);
indrm = sub2ind(pfc.adata.binSizes',indx,indy);
for unit = units,
    rateMap = pfc.plot(unit);
    wpmr(:,unit==units) = rateMap(indrm);
end
ufrwd = (myufr.data-wpmr)./sqrt(wpmr);


var(ufrwd(~isinf(ufrwd(:,3))&~isnan(ufrwd(:,3)),3))