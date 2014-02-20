function Anal_Struct = Corrected_UFR_CCG(Session,varargin)
[ bhvs,           pfcbhv, surbhv, downSampleRate, thresh_rad, test_sample_size, min_sample_size, niter] = DefaultArgs(varargin,...
{{'rear','walk'}, 'walk', 'walk',             20,       1000,               20,              10, 1000});

Session = 'jg05-20120317';
thresh_rad = 150;
test_sample_size = 7;
min_sample_size = 4;



%% Load Session if Session is not already a MTASession
if ~isa(Session,'MTASession'),
    Trial = MTATrial(Session);
    Trial.spk.create(Trial.xyz.sampleRate);
    %{'Pfs',   {'walk','rear'}}},...
end





%% Number of units
numClu = size(Trial.map,1);

newSampleRate = downSampleRate;

myxyz = Trial.xyz.copy;
myxyz.load(Trial);
myxyz.filter(gausswin(9)./sum(gausswin(9)));
myxyz.resample(newSampleRate);

myufr = Trial.ufr.copy;
myufr.create(Trial,myxyz);


%% Get the expected ufr for each xy 
%% Substract the expected ufr from the observed
pfc = MTAAknnpfs(Trial,units,pfcbhv,0,'numIter',1000,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
wpmr = ones(myxyz.size(1),numel(units));
[~,indx] = min(abs(repmat(pfc.adata.bins{1}',myxyz.size(1),1)-repmat(myxyz(:,7,1),1,numel(pfc.adata.bins{1}))),[],2);
[~,indy] = min(abs(repmat(pfc.adata.bins{2}',myxyz.size(1),1)-repmat(Trial.xyz(:,7,2),1,numel(pfc.adata.bins{2}))),[],2);
indrm = sub2ind(pfc.adata.binSizes',indx,indy);
for unit = units,
    rateMap = pfc.plot(unit);
    wpmr(:,unit==units) = rateMap(indrm);
end
ufrwd = myufr-wpmr;



%% Reduce Bhv events to a subset 
StcFilters = {{'rear',{'exclusion','rear',2},{'select_boarder_states','walk',3},{'duration',1.5}},...
              {'walk',{'exclusion','walk',2},{'duration',0.75},{'complete'}}};

%for f = numel(StcFilters),
[rear_evts,filtName] = Trial.stc.filter(StcFilters{1});
%end
[walk_evts,filtName] = Trial.stc.filter(StcFilters{2});


%% Get the ufr for each behavioral res +- 5 s
rufrs = GetSegs(ufrwd,round(((rear_evts{1}+0.5*Trial.xyz.sampleRate/newSampleRate)./Trial.xyz.sampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));
dufrs = GetSegs(ufrwd,round(((rear_evts{2}+0.5*Trial.xyz.sampleRate/newSampleRate)./Trial.xyz.sampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));
wufrs = GetSegs(ufrwd,round(((walk_evts{1}+0.5*Trial.xyz.sampleRate/newSampleRate)./Trial.xyz.sampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));
sufrs = GetSegs(ufrwd,round(((walk_evts{2}+0.5*Trial.xyz.sampleRate/newSampleRate)./Trial.xyz.sampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));

%pfMaxPosRear = zeros(size(pfr.cluMap,1),2);
pfMaxPos = zeros(size(pfc.cluMap,1),2);

for unit = 1:size(pfc.cluMap,1),
    try, pfMaxPos(unit,:) = pfc.maxRatePos{unit}(pfc.maxRateMax{unit},:);,end
    %try, pfMaxPosRear(unit,:) = pfr.maxRatePos{unit}(pfr.maxRateMax{unit},:);,end
    upos(unit,:) = pfMaxPos(unit,:);
    %upos(unit,:) = pfMaxPosRear(unit,:);
end

rposon = myxyz(round((rear_evts{1}+0.5*Trial.xyz.sampleRate/newSampleRate)./Trial.xyz.sampleRate.*newSampleRate),:);
rposoff = myxyz(round((rear_evts{2}+0.5*Trial.xyz.sampleRate/newSampleRate)./Trial.xyz.sampleRate.*newSampleRate),:);
wposon = myxyz(round((walk_evts{1}+0.5*Trial.xyz.sampleRate/newSampleRate)./Trial.xyz.sampleRate.*newSampleRate),:);
wposoff= myxyz(round((walk_evts{2}+0.5*Trial.xyz.sampleRate/newSampleRate)./Trial.xyz.sampleRate.*newSampleRate),:);

%% time bins for display pourposes
nbins = size(rufrs,1);
tbins = linspace(-5000,5000,nbins);

%% START - rear walk premutation stats
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

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
    wsurpos = myxyz(round((wsursamp+0.5*Trial.xyz.sampleRate/newSampleRate)./Trial.xyz.sampleRate.*newSampleRate),:);
    wsurdist = sqrt(sum((wsurpos-repmat(upos(unit,:),wsss,1)).^2,2));
    %% resample if too few are found
    timeout_countdown = 10;
    while sum(wsurdist<thresh_rad)<test_sample_size & timeout_countdown ~= 0 ,
        %% Draw new random indicies for surrogate selection
        wsur_rand_ind(:,iter) = randi([1,length(wsur)],wsss,1);
        %% Pull random sample from the walk surrogate tarjectories
        wsursamp = wsur(wsur_rand_ind(:,iter));
        %% Get the position of the random samples from the walk surrogate tarjectories
        wsurpos = myxyz(round((wsursamp+0.5*Trial.xyz.sampleRate/newSampleRate)./Trial.xyz.sampleRate.*newSampleRate),:);
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
                       round(((wsursamp+0.5*Trial.xyz.sampleRate/newSampleRate)...
                              ./Trial.xyz.sampleRate.*newSampleRate)...
                             -5*newSampleRate),...
                       round(10*newSampleRate));
    %% Sample Sets
    cufrs{1} = rufrs(:,rwdist<thresh_rad,unit);
    cufrs{2} = dufrs(:,rwdisto<thresh_rad,unit);
    cufrs{3} = wufrs(:,wwdist<thresh_rad,unit);
    cufrs{4} = sufrs(:,wwdisto<thresh_rad,unit);
    cufrs{5} = wsurufrs(:,wsurdist<thresh_rad);

    %% max number of trajectories for each behavioral transition
    sample_size = cellfun(@size,cufrs,repmat({2},1,ntrans));
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
                        if sample_size(ntrans)<min_sample_size|sample_size(p3)<min_sample_size,
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
                        while sum(wsurdist<thresh_rad)<test_sample_size & timeout_countdown ~= 0 ,
                            %% Draw new random indicies for surrogate selection
                            wsur_rand_ind(:,iter) = randi([1,length(wsur)],wsss,1);
                            %% Pull random sample from the walk surrogate tarjectories
                            wsursamp = wsur(wsur_rand_ind(:,iter));
                            %% Get the position of the random samples from the walk surrogate tarjectories
                            wsurpos = myxyz(round((wsursamp+0.5*Trial.xyz.sampleRate/newSampleRate)./Trial.xyz.sampleRate.*newSampleRate),:);
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
                                           round(((wsursamp+0.5*Trial.xyz.sampleRate/newSampleRate)...
                                                  ./Trial.xyz.sampleRate.*newSampleRate)...
                                                 -5*newSampleRate),...
                                           round(10*newSampleRate));
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
            for i = 1:110,
                for j = i:110,
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

figure,
unit = 14;
while unit~=-1;
clf
for k = 1:ntrans-1
subplot(ntrans-1,1,k)
hold on,
plot(nanmean(diff_stat(:,unit,5,5,:),5),'b')
plot(wrucl(:,unit,5,5),'b')
plot(wrlcl(:,unit,5,5),'b')
plot(nanmean(diff_stat(:,unit,k,k,:),5),'r')
plot(wrucl(:,unit,k,k),'r')
plot(wrlcl(:,unit,k,k),'r')
plot(nanmean(diff_stat(:,unit,5,k,:),5),'g')
plot(wrucl(:,unit,5,k),'g')
plot(wrlcl(:,unit,5,k),'g')
end
title(num2str(unit))
unit = figure_controls(gcf,unit)
end

for unit = 1:numClu;
    pron = zeros(numel(tbins),numClu);
    hron = zeros(numel(tbins),numClu);
    hwr = zeros(numel(tbins),numClu);
    pwr = zeros(numel(tbins),numClu);

        [pron(t,unit),hron(t,unit)] = ranksum(sq(diff_stat(t,unit,1,1,~isnan(diff_stat(t,unit,1,1,:)))),sq(diff_stat(t,unit,1,5,~isnan(diff_stat(t,unit,1,5,:)))));
        [pwsu(t,unit),hwsu(t,unit)] = ranksum(sq(diff_stat(t,unit,5,5,~isnan(diff_stat(t,unit,5,5,:)))),sq(diff_stat(t,unit,1,5,~isnan(diff_stat(t,unit,1,5,:)))));
        [pwr(t,unit)  ,hwr(t,unit)] = ranksum(sq(diff_stat(t,unit,5,5,~isnan(diff_stat(t,unit,5,5,:)))),sq(diff_stat(t,unit,1,1,~isnan(diff_stat(t,unit,1,1,:)))));
    end
end

if verbose, fprintf('Calculating random permuted differences for obvStates VS surState\n'),end
prankstat = zeros(numel(tbins),numClu,ntrans-1);
hrankstat = zeros(numel(tbins),numClu,ntrans-1);
tic
for unit = 1:numClu;
for p1 = 1:ntrans,
    for t = 1:numel(tbins),
        [prankstat(t,unit,p1),hrankstat(t,unit,p1)] = ranksum(sq(diff_stat(t,unit,p1,p1,~isnan(diff_stat(t,unit,p1,p1,:)))),sq(diff_stat(t,unit,p1,ntrans,~isnan(diff_stat(t,unit,p1,ntrans,:)))));
    end
end
end
toc



cufrccg.idd = zeros(nbins,numClu,ntrans,niter);
cufrccg.pdd = zeros(nbins,numClu,ntrans,niter);
cufrccg.edd = zeros(nbins,numClu,     1,niter);

cufrccg.midd = zeros(nbins,numClu,ntrans);
cufrccg.mpdd = zeros(nbins,numClu,ntrans);
cufrccg.medd = zeros(nbins,numClu);

cufrccg.idd_edd_prs = prankstat;
cufrccg.idd_edd_hrs = hrankstat;


for p1 = 1:ntrans,cufrccg.idd(:,:,p1,:) = sq(diff_stat(:,:,p1,p1,:));end
for p1 = 1:ntrans,cufrccg.pdd(:,:,p1,:) = sq(diff_stat(:,:,p1,ntrans,:));end
for p1 = 1:ntrans,cufrccg.midd(:,:,p1) = sq(nanmean(diff_stat(:,:,p1,p1,:),5));end
for p1 = 1:ntrans,cufrccg.mpdd(:,:,p1) = sq(nanmean(diff_stat(:,:,p1,ntrans,:),5));end
cufrccg.edd(:,:,1,:) = sq(diff_stat(:,:,ntrans,1,:));
cufrccg.medd(:,:) = sq(nanmean(diff_stat(:,:,ntrans,1,:),5));



%% Calculate the SNR between the intra-distribution-difference and
%% the permuted-distribution-difference
cufrccg.pdd_idd_snr = zeros(nbins,numClu,ntrans);
tic
for unit = 1:numClu,
for p1 = 1:ntrans,
cufrccg.pdd_idd_snr(:,unit,p1) = (nanmean(diff_stat(:,unit,p1,p1,1:niter),5)-nanmean(diff_stat(:,unit,ntrans,p1,1:niter),5))./(nanvar(diff_stat(:,unit,p1,p1,1:niter),[],5)+nanvar(diff_stat(:,unit,ntrans,p1,1:niter),[],5));
end
end
toc
%% Calculate the SNR between the intra-distribution-difference and
%% the permuted-distribution-difference
cufrccg.pdd_edd_snr = zeros(nbins,numClu,ntrans);
tic
for unit = 1:numClu,
for p1 = 1:ntrans,
cufrccg.pdd_edd_snr(:,unit,p1) = (nanmean(diff_stat(:,unit,ntrans,ntrans,1:niter),5)-nanmean(diff_stat(:,unit,ntrans,p1,1:niter),5))./(nanvar(diff_stat(:,unit,ntrans,ntrans,1:niter),[],5)+nanvar(diff_stat(:,unit,ntrans,p1,1:niter),[],5));
end
end
toc


%% PVAL of mpdd <
cufrccg.mpdd_pval = zeros(nbins,numClu,ntrans);
mpdd_pval = sum((repmat(diff_stat(:,:,ntrans,ntrans,1:niter),[1,1,ntrans,ntrans,1])<repmat(nanmean(diff_stat(:,:,ntrans,:,1:niter),5),[1,1,ntrans,1,niter])),5)./nnsamp;
mpdd_pval(mpdd_pval==1) = mpdd_pval(mpdd_pval==1) - (1/nnsamp(mpdd_pval==1))';
mpdd_pval(mpdd_pval==0) = mpdd_pval(mpdd_pval==0) + (1/nnsamp(mpdd_pval==0))';
tic
for unit = 1:numClu,
for p1 = 1:ntrans,
cufrccg.mpdd_pval(:,unit,p1) = mpdd_pval(:,unit,p1,p1);
end
end
toc

cufrccg.mpdd_pval(cufrccg.mpdd_pval>0.5) =cufrccg.mpdd_pval(cufrccg.mpdd_pval>0.5)-1;


anal_struct.data = cufrccg;
anal_struct.vars = 
