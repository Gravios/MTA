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
Trial.stc.updateMode('auto_wbhr');
Trial.stc.load;


%% Get units which are of good enough quality
units = find(Trial.nq.SpkWidthR>0.5&Trial.nq.eDist>30)';
numClu = numel(units);

newSampleRate = downSampleRate;

myxyz = Trial.xyz.copy;
myxyz.load(Trial);
myxyz.filter(gausswin(9)./sum(gausswin(9)));
myxyz.resample(newSampleRate);
myxyz.data = sq(myxyz(:,7,[1,2]));

myufr = Trial.ufr.copy;
myufr.create(Trial,myxyz,[],units);


%% Get the expected ufr for each xy 
%% Substract the expected ufr from the observed
pfc = MTAAknnpfs(Trial,units,pfcbhv,0,'numIter',0,'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',70);
wpmr = ones(myxyz.size(1),numel(units));
[~,indx] = min(abs(repmat(pfc.adata.bins{1}',myxyz.size(1),1)-repmat(myxyz(:,1),1,numel(pfc.adata.bins{1}))),[],2);
[~,indy] = min(abs(repmat(pfc.adata.bins{2}',myxyz.size(1),1)-repmat(myxyz(:,2),1,numel(pfc.adata.bins{2}))),[],2);
indrm = sub2ind(pfc.adata.binSizes',indx,indy);
for unit = units,
    rateMap = pfc.plot(unit);
    wpmr(:,unit==units) = rateMap(indrm);
end
ufrwd = myufr.data-wpmr;

%% Reduce Bhv events to a subset 
StcFilters = {{'rear',{'exclusion','rear',2},{'select_boarder_states','walk',3},{'duration',1.5}},...
              {'walk',{'exclusion','walk',2},{'duration',0.75},{'complete'}}};

          
%for f = numel(StcFilters),
[rear_evts,filtName] = Trial.stc.filter(newSampleRate,StcFilters{1});
%end
[walk_evts,filtName] = Trial.stc.filter(newSampleRate,StcFilters{2});


%% Get the ufr for each behavioral res +- 5 s
rufrs = GetSegs(ufrwd,round(rear_evts{1}-3*newSampleRate),round(6*newSampleRate),[]);
dufrs = GetSegs(ufrwd,round(rear_evts{2}-3*newSampleRate),round(6*newSampleRate),[]);
wufrs = GetSegs(ufrwd,round(walk_evts{1}-3*newSampleRate),round(6*newSampleRate),[]);
sufrs = GetSegs(ufrwd,round(walk_evts{2}-3*newSampleRate),round(6*newSampleRate),[]);

%pfMaxPosRear = zeros(size(pfr.cluMap,1),2);
[urate,upos] = pfc.maxRate(units);

rposon =  myxyz(rear_evts{1},:);
rposoff = myxyz(rear_evts{2},:);
wposon = myxyz(walk_evts{1},:);
wposoff= myxyz(walk_evts{2},:);

%% time bins for display pourposes
nbins = size(rufrs,1);
tbins = linspace(-3000,3000,nbins);

%% START - rear walk premutation stats
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

ntrans = 5;
perm_stat = zeros(nbins,numClu,ntrans,ntrans,niter);
diff_stat = zeros(nbins,numClu,ntrans,ntrans,niter);

%% Generate surrogate walk segments from place field max rate pos

wsur = [];
for i = 1:numel(walk_evts{1})
    wsur = cat(1,wsur,[walk_evts{1}(i):walk_evts{2}(i)]');
end


%% Generate random samples from the walk surrogates
wsss = size(wposon,1);
gss = zeros(numClu,ntrans);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for unit = 1:size(rufrs,3),
    fprintf('Starting unit: %i \n',unit)
    iter = 1;
    %% Rear onset pfc distance
    rwdist = sqrt(sum((rposon-repmat(upos(unit,:),size(rposon,1),1)).^2,2));
    %% Rear offset pfc distance
    rwdisto = sqrt(sum((rposoff-repmat(upos(unit,:),size(rposoff,1),1)).^2,2));
    %% Walk onset pfc distance
    wwdist = sqrt(sum((wposon-repmat(upos(unit,:),size(wposon,1),1)).^2,2));
    %% Walk offset pfc distance
    wwdisto = sqrt(sum((wposoff-repmat(upos(unit,:),size(wposoff,1),1)).^2,2));
    wsur_rand_ind = randi([1,length(wsur)],wsss,niter);
    wsursamp = wsur(wsur_rand_ind(:,iter));
    wsurpos = myxyz(wsursamp,:);
    wsurdist = sqrt(sum((wsurpos-repmat(upos(unit,:),wsss,1)).^2,2));
    %% resample if too few are found
    timeout_countdown = 10;
    while sum(wsurdist<thresh_rad)<test_sample_size && timeout_countdown ~= 0 ,
        %% Draw new random indicies for surrogate selection
        wsur_rand_ind(:,iter) = randi([1,length(wsur)],wsss,1);
        %% Pull random sample from the walk surrogate tarjectories
        wsursamp = wsur(wsur_rand_ind(:,iter));
        %% Get the position of the random samples from the walk surrogate tarjectories
        wsurpos = myxyz(wsursamp,:);
        %% Calc the distance of each sample to the center of the current units placefield
        wsurdist = sqrt(sum((wsurpos-repmat(upos(unit,:),wsss,1)).^2,2));                        
        timeout_countdown = timeout_countdown - 1;
    end
    %% give up on resampling if too few are found too many times in a row
    if timeout_countdown==0,
        continue
    end
    %% Get the unit firing rates for the surrogate walking segments
    wsurufrs = GetSegs(ufrwd(:,unit),...
                       round(wsursamp-3*newSampleRate),...
                       round(6*newSampleRate));
    %% Sample Sets
    cufrs{1} = rufrs(:,rwdist<thresh_rad,unit);
    cufrs{2} = dufrs(:,rwdisto<thresh_rad,unit);
    cufrs{3} = wufrs(:,wwdist<thresh_rad,unit);
    cufrs{4} = sufrs(:,wwdisto<thresh_rad,unit);
    cufrs{5} = wsurufrs(:,wsurdist<thresh_rad);

    %% max number of trajectories for each behavioral transition
    sample_size = cellfun(@size,cufrs,repmat({2},1,ntrans));
    gss(unit,:) = sample_size;
    rand_samp_ind = zeros(test_sample_size,ntrans,2,niter);
    for i = 1:ntrans,
        for j = 1:2,
            if sample_size(i)<min_sample_size,continue,end
            rand_samp_ind(:,i,j,:) = randi([1,sample_size(i)],test_sample_size,1,1,niter);
        end
    end
    tic
    for p1 = 1:ntrans ,
        for p2 = p1:ntrans ,             

            if p1~=ntrans & p2~=ntrans ,
                %% IF too few samples, don't bother with permutations
                if sample_size(p1)<min_sample_size|sample_size(p2)<min_sample_size,
                    continue
                end
                for iter = 1:niter,
                    perm_stat(:,unit,p1,p2,iter) = nanmean(cat(2,cufrs{p1}(:,rand_samp_ind(:,p1,1,iter)),...
                                                               cufrs{p2}(:,rand_samp_ind(:,p2,2,iter))),2);
                end
            elseif p1==1 & p2==ntrans ,
                for iter = 1:niter ,
                    for p3 = 1:ntrans ,
                        %% IF too few samples, don't bother with permutations
                        if sample_size(ntrans)<min_sample_size||sample_size(p3)<min_sample_size,
                            continue
                        end
                        %% Pull random sample from the walk surrogate tarjectories
                        wsursamp = wsur(wsur_rand_ind(:,iter));
                        %% Get the position of the random samples from the walk surrogate tarjectories
                        wsurpos = myxyz(round((wsursamp+0.5*Trial.xyz.sampleRate/newSampleRate)./Trial.xyz.sampleRate.*newSampleRate),:);
                        %% Calc the distance of each sample to the center of the current units placefield
                        wsurdist = sqrt(sum((wsurpos-repmat(upos(unit,:),wsss,1)).^2,2));
                        %% resample if too few are found
                        timeout_countdown = 10;
                        while sum(wsurdist<thresh_rad)<test_sample_size && timeout_countdown ~= 0 ,
                            %% Draw new random indicies for surrogate selection
                            wsur_rand_ind(:,iter) = randi([1,length(wsur)],wsss,1);
                            %% Pull random sample from the walk surrogate tarjectories
                            wsursamp = wsur(wsur_rand_ind(:,iter));
                            %% Get the position of the random samples from the walk surrogate tarjectories
                            wsurpos = myxyz(wsursamp,:);
                            %% Calc the distance of each sample to the center of the current units placefield
                            wsurdist = sqrt(sum((wsurpos-repmat(upos(unit,:),wsss,1)).^2,2));                        
                            timeout_countdown = timeout_countdown - 1;
                        end
                        %% give up on resampling and just fill the
                        %% iterations row with nan's which can give us
                        %% an idea how 'cramped/undersampled' an area is
                        if timeout_countdown==0,
                            perm_stat(:,unit,p3,ntrans,iter) = nan(numel(tbins),1);
                            diff_stat(:,unit,p3,ntrans,iter) = nan(numel(tbins),1);
                            continue
                        end

                        %% Get the unit firing rate time series around each sample +-5 seconds
                        wsurufrs = GetSegs(ufrwd(:,unit),...
                                           round(wsursamp-3*newSampleRate),...
                                           round(6*newSampleRate));
                        %% Select sample rates base on sample distance from placefield center
                        cufrs{5} = wsurufrs(:,wsurdist<thresh_rad);
                        %% Adjust the stored value for walk surrogate sample size, which varies with each sample
                        sample_size(5) = size(cufrs{5},2);
                        %% Adjust the random selection indicies for sample permutation sub sampling
                        rand_samp_ind(:,5,1,:) = randi([1,size(cufrs{ntrans},2)],test_sample_size,1,1,niter);
                        rand_samp_ind(:,5,2,:) = randi([1,size(cufrs{ntrans},2)],test_sample_size,1,1,niter);

                        perm_stat(:,unit,p3,ntrans,iter) = nanmean(cat(2,cufrs{p3}(:,rand_samp_ind(:,p3,1,iter)),...
                                                                       cufrs{ntrans}(:,rand_samp_ind(:,ntrans,2,iter))),2);
                        diff_stat(:,unit,ntrans,p3,iter) =  nanmean(cufrs{p3}(:,rand_samp_ind(:,p3,1,iter)),2)...
                                                           -nanmean(cufrs{ntrans}(:,rand_samp_ind(:,ntrans,2,iter)),2);
                    end
                end
            end

        end
    end
    toc
end




count = 1;
if verbose, fprintf('Calculating random permuted differences for obvStates VS surState\n'),end
tic
for p1 = 1:ntrans,
    for p2 = p1:ntrans,
        count = 1;
        if p1==p2|p2==ntrans,
            for i =  randi(1000,[1,110]),
                for j = randi(1000,[1,110]),
                    if i==j,continue,end
                    diff_stat(:,:,p1,p2,count) = perm_stat(:,:,p1,p2,i)-perm_stat(:,:,p1,p2,j);
                    if count == niter , break , end
                    count = count+1;
                end
                if count == niter , break , end
            end
        end
    end
end
toc




psig = .05;
wrsps = sort(diff_stat(:,:,:,:,1:niter),5);
nnsamp = zeros(nbins,numClu,ntrans,ntrans);
tic
for unit = 1:numClu,
    for p1 = 1:ntrans,
        for p2 = 1:ntrans,
            for t = 1:nbins,
                wnind = sq(~isnan(wrsps(t,unit,p1,p2,1:niter)));
                nwnan = sum(~wnind);
                nnsamp(t,unit,p1,p2) = niter-nwnan; 
                if nwnan == niter, continue,end
                wrucl(t,unit,p1,p2) = wrsps(t,unit,p1,p2,round((niter-nwnan)*(1-psig)));
                wrlcl(t,unit,p1,p2) = wrsps(t,unit,p1,p2,round((niter-nwnan)*psig));
            end
        end
    end
end
toc



%% Calculating random permuted differences for obvStates VS surState
prankstat = zeros(numel(tbins),numClu,ntrans-1);
hrankstat = zeros(numel(tbins),numClu,ntrans-1);
tic
for u = 1:numClu;
    for p1 = 1:ntrans,
        for t = 1:nbins,
            [prankstat(t,u,p1),hrankstat(t,u,p1)] = ranksum(sq(diff_stat(t,u,p1,p1,~isnan(diff_stat(t,u,p1,p1,:)))),sq(diff_stat(t,u,p1,ntrans,~isnan(diff_stat(t,u,p1,ntrans,:)))));
        end
    end
end
toc



auxdata.numUnits = numClu;
auxdata.nbins = nbins;
auxdata.tbins = tbins;
auxdata.tbound = find(auxdata.tbins>-2000&auxdata.tbins<2000);
auxdata.bname = [];
auxdata.ntrans = ntrans;
auxdata.newSampleRate = newSampleRate;

%% metadata
data.filebase = repmat({Trial.filebase},1,numClu);
data.clu = units;
data.el = Trial.spk.map(units,2)';


%% data
data.idd = zeros(nbins,numClu,ntrans,niter);
data.idd_p95 = zeros(nbins,numClu,ntrans);
data.idd_p05 = zeros(nbins,numClu,ntrans);
data.pdd = zeros(nbins,numClu,ntrans,niter);
data.pdd_p95 = zeros(nbins,numClu,ntrans);
data.pdd_p05 = zeros(nbins,numClu,ntrans);
data.edd = zeros(nbins,numClu,ntrans,niter);
data.edd_p95 = zeros(nbins,numClu,ntrans);
data.edd_p05 = zeros(nbins,numClu,ntrans);

data.midd = zeros(nbins,numClu,ntrans);
data.mpdd = zeros(nbins,numClu,ntrans);
data.medd = zeros(nbins,numClu);

data.idd_pdd_prs = prankstat;
data.idd_pdd_hrs = hrankstat;
data.bhv_sample_size = gss';
          
for p1 = 1:ntrans,data.idd(:,:,p1,:) = sq(diff_stat(:,:,p1,p1,:));end
for p1 = 1:ntrans,data.idd_p95(:,:,p1) = wrucl(:,:,p1,p1);,end
for p1 = 1:ntrans,data.idd_p05(:,:,p1) = wrlcl(:,:,p1,p1);,end

for p1 = 1:ntrans,data.pdd(:,:,p1,:) = sq(diff_stat(:,:,p1,ntrans,:));end
for p1 = 1:ntrans,data.pdd_p95(:,:,p1) = wrucl(:,:,p1,ntrans);,end
for p1 = 1:ntrans,data.pdd_p05(:,:,p1) = wrlcl(:,:,p1,ntrans);,end

for p1 = 1:ntrans,data.edd(:,:,p1,:) = sq(diff_stat(:,:,ntrans,p1,:));,end
for p1 = 1:ntrans,data.edd_p95(:,:,p1) = wrucl(:,:,ntrans,p1);,end
for p1 = 1:ntrans,data.edd_p05(:,:,p1) = wrlcl(:,:,ntrans,p1);,end

for p1 = 1:ntrans,data.midd(:,:,p1) = sq(nanmean(diff_stat(:,:,p1,p1,:),5));end
for p1 = 1:ntrans,data.mpdd(:,:,p1) = sq(nanmean(diff_stat(:,:,p1,ntrans,:),5));end
for p1 = 1:ntrans,data.medd(:,:,p1) = sq(nanmean(diff_stat(:,:,ntrans,p1,:),5));end




%% Calculate the SNR between the intra-distribution-difference and
%% the permuted-distribution-difference
data.pdd_idd_snr = zeros(nbins,numClu,ntrans);
tic
for u = 1:numClu,
    for p1 = 1:ntrans,
        data.pdd_idd_snr(:,u,p1) = (nanmean(diff_stat(:,u,p1,p1,1:niter),5)-nanmean(diff_stat(:,u,ntrans,p1,1:niter),5))./(nanvar(diff_stat(:,u,p1,p1,1:niter),[],5)+nanvar(diff_stat(:,u,ntrans,p1,1:niter),[],5));
    end
end
toc
%% Calculate the SNR between the intra-distribution-difference and
%% the permuted-distribution-difference
data.pdd_edd_snr = zeros(nbins,numClu,ntrans);
tic
for u = 1:numClu,
    for p1 = 1:ntrans,
        data.pdd_edd_snr(:,u,p1) = (nanmean(diff_stat(:,u,ntrans,ntrans,1:niter),5)-nanmean(diff_stat(:,u,ntrans,p1,1:niter),5))./(nanvar(diff_stat(:,u,ntrans,ntrans,1:niter),[],5)+nanvar(diff_stat(:,u,ntrans,p1,1:niter),[],5));
    end
end
toc


%% PVAL of mean of edd vs dist of surrogate pdd
data.medd_pval = zeros(nbins,numClu,ntrans);
medd_pval = sum((repmat(diff_stat(:,:,ntrans,ntrans,1:niter),[1,1,ntrans,ntrans,1])<repmat(nanmean(diff_stat(:,:,ntrans,:,1:niter),5),[1,1,ntrans,1,niter])),5)./nnsamp;

medd_pval(medd_pval==1) = medd_pval(medd_pval==1) - 1./nnsamp(medd_pval==1);
medd_pval(medd_pval==0) = medd_pval(medd_pval==0) + 1./nnsamp(medd_pval==0);

medd_pval = 0.5-abs(medd_pval-0.5);

tic
for u = 1:numClu,
    for p1 = 1:ntrans,
        data.medd_pval(:,u,p1) = medd_pval(:,u,p1,p1);
    end
end
toc

data.medd_pval(data.medd_pval>0.5) =data.medd_pval(data.medd_pval>0.5)-1;
data.nnsamp = nnsamp;


data = StructArray(data,2);

save(fullfile(Trial.spath, [Trial.filebase '.cufrccg.' ...
      num2str(thresh_rad) '.pfc_' pfcbhv '.sur_' surbhv ...
      '.mat']), 'data','auxdata','-v7.3')
          
          % 

% sts_evts = {};
% evt_name = {};
% for f = 1:numel(StcFilters),
%     [f_evts,filterName] = Trial.stc.filter(newSampleRate,StcFilters{f});
%     sts_evts =  cat(2,sts_evts,f_evts);
%     evt_name =  cat(2,evt_name,[filterName1 '_on'],[filterName1 '_off']);
% end
% 
% sufrs = {};
% 
% for i = 1:numel(sts_evts)
%     sufrs = cat(2,sufrs,GetSegs(ufrwd,round(sts_evts{i}-3*newSampleRate),round(6*newSampleRate),0));
% end
% 
% 
% [urate,upos] = pfc.maxRate(units);
% 
% for i = 1:numel(sts_evts)
%     spos{i} = myxyz(rear_evts{i},:);
% end
% 
% nbins = size(sufrs{1},1);
% tbins = linspace(-3000,3000,nbins);
% 
% %% START - rear walk premutation stats
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
% 
% ntrans = 5;
% perm_stat = zeros(nbins,numClu,ntrans,ntrans,niter);
% diff_stat = zeros(nbins,numClu,ntrans,ntrans,niter);
% 
% 
% wsur = [];
% for i = 1:numel(sts_evts{3})
%     wsur = cat(1,wsur,[sts_evts{3}(i):sts_evts{4}(i)]');
% end
% 
% wsss = size(wposon,1);



