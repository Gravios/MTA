function Anal_Struct = newcufrccg(Session,varargin)
[ bhvs,           pfcbhv, surbhv, downSampleRate, thresh_rad, test_sample_size, min_sample_size, niter] = DefaultArgs(varargin,...
{{'rear&theta','walk&theta'}, 'walk&theta', 'walk&theta',             20,       1000,               20,              10, 1000});

Session = 'jg05-20120317';
thresh_rad = 150;
test_sample_size = 7;
min_sample_size = 4;



%% Load Session if Session is not already a MTASession
if ~isa(Session,'MTASession'),
    Trial = MTATrial(Session);
    Trial.load('nq');
    %{'Pfs',   {'walk','rear'}}},...
end



%% Get units which are of good enough quality
units = find(Trial.nq.SpkWidthR>0.3&Trial.nq.eDist>19)';
numClu = numel(units);

newSampleRate = downSampleRate;

myxyz = Trial.xyz.copy;
myxyz.load(Trial);
myxyz.filter(gausswin(9)./sum(gausswin(9)));
myxyz.resample(newSampleRate);

myufr = Trial.ufr.copy;
myufr.create(Trial,myxyz,[],units);


%% Get the expected ufr for each xy 
%% Substract the expected ufr from the observed
pfc = MTAAknnpfs(Trial,units,pfcbhv,0,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',70);
wpmr = ones(myxyz.size(1),numel(units));
[~,indx] = min(abs(repmat(pfc.adata.bins{1}',myxyz.size(1),1)-repmat(myxyz(:,7,1),1,numel(pfc.adata.bins{1}))),[],2);
[~,indy] = min(abs(repmat(pfc.adata.bins{2}',myxyz.size(1),1)-repmat(myxyz(:,7,2),1,numel(pfc.adata.bins{2}))),[],2);
indrm = sub2ind(pfc.adata.binSizes',indx,indy);
for unit = units,
    rateMap = pfc.plot(unit);
    wpmr(:,unit==units) = rateMap(indrm);
end
ufrwd = myufr.data-wpmr;

%% Reduce Bhv events to a subset 
StcFilters = {{'rear',{'exclusion','rear',2},{'select_boarder_states','walk',3},{'duration',1.5}},...
              {'walk',{'exclusion','walk',2},{'duration',0.75},{'complete'}}};

sts_evts = {};
evt_name = {};
for f = 1:numel(StcFilters),
    [f_evts,filterName] = Trial.stc.filter(newSampleRate,StcFilters{f});
    sts_evts =  cat(2,sts_evts,f_evts);
    evt_name =  cat(2,evt_name,filterName);
end

sufrs = {};

sufrs = cat(2,sufrs,GetSegs(ufrwd,round(sts_evts{i}-3*newSampleRate),round(6*newSampleRate),0);
