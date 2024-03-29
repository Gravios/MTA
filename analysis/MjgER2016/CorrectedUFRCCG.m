function CorrectedUFRCCG(Session,mode,varargin)
%function [data,parameters] = cufrccg(Session,varargin)
[ bhvf,           pfcbhv, surbhv, klen, overlap, thresh_rad, test_sample_size, min_sample_size, niter,states] = DefaultArgs(varargin,...
{{{'rear',{'exclusion','rear',2},{'select_boarder_states','walk',3},{'duration',1.5}},...
  {'walk',{'exclusion','walk',2},{'duration',0.75},{'complete'}}},...
                  'walk', 'walk',   64,      16,       200,                7,               4, 5000,'theta'});


parameters.bhvf = bhvf;
parameters.pfcbhv = pfcbhv;
parameters.surbhv = surbhv;
parameters.klen = klen;
parameters.overlap = overlap;
parameters.thresh_rad = thresh_rad;
parameters.test_sample_size = test_sample_size;
parameters.min_sample_size = min_sample_size;
parameters.niter = niter;
parameters.states = states;


switch mode

case 'compute'


%% TESTING VARS
% Session = 'jg05-20120309';
% $$$ 
% $$$ bhvf = {{'rear',{'exclusion','rear',2},{'select_boarder_states','walk',3},{'duration',1.5}},...
% $$$   {'walk',{'exclusion','walk',2},{'duration',0.75},{'complete'}}};
% $$$  
% $$$ thresh_rad = 200;
% $$$ test_sample_size = 7;
% $$$ min_sample_size = 4;
%% END TESTING VARS


%% Load Session if Session is not already a MTASession
if ~isa(Session,'MTASession'),
    s = MTASession(Session);
    Trial = MTATrial(s,...
                    {{'CluRes',s.xyzSampleRate,[],states}},...
                    'all');
    clear('s');
end

%%C JUST IN CASE
Trial.Bhv.sampleRate = Trial.xyzSampleRate;


%%M SELECT units were selected for computation based on their firing
%%M rate, and separation from other clusters.
query = {{MTAPlaceField(Trial,[],'walk'),'maxRate',@gt,3},...
          @or,...
         {MTAPlaceField(Trial,[],'rear'),'maxRate',@gt,3}};
units = Trial.selectUnits(query);


numUnits = numel(units);


%%M The rate correction a.... FINISH THIS DESCRIPTION LATER

%%C Get the placeField whos ratemap centers will be used for
%%C spatially selective cufrccg analysis
pfc = Trial.getPfs(pfcbhv);
pfs = Trial.getPfs(surbhv);



%%M The xyz coordinates of the tracking marker's trajectory were
%%M smoothed and downsampled by averaging over a moving window.

%%C Filter xyz 
Trial = Trial.filter();



%%G Generate myxyz & newSampleRate
%%C Kernal with klen number of bins
kern = ones(klen,1);
%%C Get the mean position of each bin in the xy coordinates
%%C Resample xyz to new bin intervals
myxyz = sq(Trial.xyz(:,7,[1,2]));
t =         permute(reshape(          myxyz(1:size(myxyz,1)-mod(size(myxyz,1),klen),:),klen,[],2),[4,1,2,3]);
for shift = 1:klen/overlap-1,
t = cat(1,t,permute(reshape(circshift(myxyz(1:size(myxyz,1)-mod(size(myxyz,1),klen),:),-overlap*shift),klen,[],2),[4,1,2,3]));
end
myxyz = t;
myxyz = reshape(sq(sum(repmat(permute(repmat(permute(repmat(kern./sum(kern),1,size(myxyz,3)),[5,4,1,2,3]),size(myxyz,4),1),[2,3,4,1]),klen/overlap,1).*myxyz,2)),[],2);
%% Calculate the new sampling rate for indexing purposes 
newSampleRate = 1/((size(Trial.xyz,1)-mod(size(Trial.xyz,1),klen))/Trial.xyzSampleRate/length(myxyz));



%% Get the unit firing rate (ufr) for each unit
%% Resample ufr to new bin intervals
myufr = {};
for unit = units'
fprintf('Calculating ufr for unit: %i \n',unit)
myres = round(Trial.res(Trial.clu==unit));
myufr{end+1} = zeros(size(Trial.xyz,1),1);
if ~isempty(myres),
   myufr{end}(1:myres(end)) = accumarray(myres,1);
end
t =         permute(reshape(myufr{end}(1:size(myufr{end},1)-mod(size(myufr{end},1),klen)),klen,[]),[3,1,2]);
for shift = 1:klen/overlap-1
t = cat(1,t,permute(reshape(circshift(myufr{end}(1:size(myufr{end},1)-mod(size(myufr{end},1),klen)),-overlap*shift),klen,[]),[3,1,2]));
end
myufr{end} = t; 
myufr{end} = reshape(sq(sum(repmat(permute(repmat(kern,1,size(myufr{end},3)),[3,1,2]),size(myufr{end},1),1).*myufr{end},2))/(klen/Trial.xyzSampleRate),[],1);
end
myufr = cell2mat(myufr);



%% Get the expected ufr for each xy 
wpmr = zeros(size(myxyz,1),numUnits);
for u =1:numUnits
fprintf(['Correcting ufr for unit: %i \nUsing expected rates from the ' ...
'%s placefield\n'],units(u),pfs.stateLabel)
    [~,indx] = min(abs(repmat(pfs.xbin,size(myxyz,1),1)-repmat(myxyz(:,1),1,numel(pfs.xbin))),[],2);
    [~,indy] = min(abs(repmat(pfs.ybin,size(myxyz,1),1)-repmat(myxyz(:,1),1,numel(pfs.ybin))),[],2);
    indrm = sub2ind([numel(pfs.xbin),numel(pfs.ybin)],indx,indy);
    wpmr(:,u) = pfs.rateMap{units(u)}(indrm);
end
%% Substract the expected ufr from the observed
ufrwd = myufr-wpmr;



%% Reduce Bhv events to a subset 
bhv_evts = {};
filtName = {};
for i = 1:numel(bhvf),
    [bhv_evts{end+1},filtName{end+1}] = Trial.Bhv.filter(bhvf{i});
end


%% Get ufr segments for each behavioral res +- 5 s
bufrs = {};
bst = {'on','off'};
bname = {};
for i = 1:numel(bhv_evts),
    for j = 1:numel(bhv_evts{i}),
        bufrs{end+1} = GetSegs(ufrwd,round(((bhv_evts{i}{j}+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));
        bname{end+1} = [filtName{i} '_' bst{j}];
    end
end



%% Find the position of the maximum expected rate for each unit
pfMaxPos = zeros(numUnits,2);
for unit = 1:numUnits
    try, pfMaxPos(u,:) = pfc.maxRatePos{units(u)}(pfc.maxRateMax{units(u)},:);,end
end



%% Find the position of each bhv event
bpos = {};
for i = 1:numel(bhv_evts),
    for j = 1:numel(bhv_evts{i}),
        bpos{end+1} = myxyz(round((bhv_evts{i}{j}+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);
    end
end


%% Create surrogate 
% TODO add some other checks for more complicated surrogates
surind = []
for i = 1:numel(bhvf)
    if ischar(bhvf{i}{1}),
        if strcmp(bhvf{i}{1},surbhv),
            surind = i;
        end
    end
end
bsname = [filtName{surind} '_sur'];
bssur = [];
for i = 1:numel(bhv_evts{surind}{1})
    bssur = cat(1,bssur,[bhv_evts{surind}{1}(i):bhv_evts{surind}{2}(i)]');
end


%% time bins for display pourposes
nbins = size(bufrs{1},1);
tbins = linspace(-5000,5000,nbins);


%% START - rear walk premutation stats
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
bname{end+1} = bsname;

ntrans = numel(bpos)+1;
perm_stat = zeros(nbins,numUnits,ntrans,ntrans,niter);
diff_stat = zeros(nbins,numUnits,ntrans,ntrans,niter);

%% Generate surrogate walk segments from place field max rate pos


%% Generate random samples from the walk surrogates
%% max samplesize
wsss = round(mean([size(bhv_evts{surind}{1},1),size(bhv_evts{surind}{2},1)]));
gss = zeros(numUnits,ntrans);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for u = 1:numUnits,
    fprintf('Starting unit: %i \n',units(u))
    iter = 1;
    
    %% bhv onset/offset pfc distance
    bdist = {};
    for i = 1:numel(bpos),
        bdist{end+1} =  sqrt(sum((bpos{i}-repmat(pfMaxPos(u,:),size(bpos{i},1),1)).^2,2));
    end
 
    bssur_rand_ind = randi([1,length(bssur)],wsss,niter);
    bssursamp = bssur(bssur_rand_ind(:,iter));
    bssurpos = myxyz(round((bssursamp+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);
    bssurdist = sqrt(sum((bssurpos-repmat(pfMaxPos(u,:),wsss,1)).^2,2));
    %% resample if too few are found
    timeout_countdown = 10;
    while sum(bssurdist<thresh_rad)<test_sample_size & timeout_countdown ~= 0 ,
        %% Draw new random indicies for surrogate selection
        bssur_rand_ind(:,iter) = randi([1,length(bssur)],wsss,1);
        %% Pull random sample from the walk surrogate tarjectories
        bssursamp = bssur(bssur_rand_ind(:,iter));
        %% Get the position of the random samples from the walk surrogate tarjectories
        bssurpos = myxyz(round((bssursamp+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);
        %% Calc the distance of each sample to the center of the current units placefield
        bssurdist = sqrt(sum((bssurpos-repmat(pfMaxPos(u,:),wsss,1)).^2,2));                        
        timeout_countdown = timeout_countdown - 1;
    end
    %% give up on resampling if too few are found too many times in a row
    if timeout_countdown==0,
        continue
    end
    %% Get the unit firing rates for the surrogate walking segments
    bssurufrs = GetSegs(ufrwd(:,u),...
                       round(((bssursamp+0.5*Trial.xyzSampleRate/newSampleRate)...
                              ./Trial.xyzSampleRate.*newSampleRate)...
                             -5*newSampleRate),...
                       round(10*newSampleRate));
    %% Sample Sets
    cufrs = {};
    for i = 1:numel(bpos),
        cufrs{i} = bufrs{i}(:,bdist{i}<thresh_rad,u);
    end
    cufrs{ntrans} = bssurufrs(:,bssurdist<thresh_rad);

    %% max number of trajectories for each behavioral transition
    sample_size = cellfun(@size,cufrs,repmat({2},1,ntrans));
    gss(u,:) = sample_size;
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
                    perm_stat(:,u,p1,p2,iter) = nanmean(cat(2,cufrs{p1}(:,rand_samp_ind(:,p1,1,iter)),...
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
                        bssursamp = bssur(bssur_rand_ind(:,iter));
                        %% Get the position of the random samples from the walk surrogate tarjectories
                        bssurpos = myxyz(round((bssursamp+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);
                        %% Calc the distance of each sample to the center of the current units placefield
                        bssurdist = sqrt(sum((bssurpos-repmat(pfMaxPos(u,:),wsss,1)).^2,2));
                        %% resample if too few are found
                        timeout_countdown = 10;
                        while sum(bssurdist<thresh_rad)<test_sample_size & timeout_countdown ~= 0 ,
                            %% Draw new random indicies for surrogate selection
                            bssur_rand_ind(:,iter) = randi([1,length(bssur)],wsss,1);
                            %% Pull random sample from the walk surrogate tarjectories
                            bssursamp = bssur(bssur_rand_ind(:,iter));
                            %% Get the position of the random samples from the walk surrogate tarjectories
                            bssurpos = myxyz(round((bssursamp+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);
                            %% Calc the distance of each sample to the center of the current units placefield
                            bssurdist = sqrt(sum((bssurpos-repmat(pfMaxPos(u,:),wsss,1)).^2,2));                        
                            timeout_countdown = timeout_countdown - 1;
                        end
                        %% give up on resampling and just fill the
                        %% iterations row with nan's which can give us
                        %% an idea how 'cramped/undersampled' an area is
                        if timeout_countdown==0,
                            perm_stat(:,u,p3,ntrans,iter) = nan(numel(tbins),1);
                            diff_stat(:,u,p3,ntrans,iter) = nan(numel(tbins),1);
                            continue
                        end

                        %% Get the unit firing rate time series around each sample +-5 seconds
                        bssurufrs = GetSegs(ufrwd(:,u),...
                                           round(((bssursamp+0.5*Trial.xyzSampleRate/newSampleRate)...
                                                  ./Trial.xyzSampleRate.*newSampleRate)...
                                                 -5*newSampleRate),...
                                           round(10*newSampleRate));
                        %% Select sample rates base on sample distance from placefield center
                        cufrs{ntrans} = bssurufrs(:,bssurdist<thresh_rad);
                        %% Adjust the stored value for walk surrogate sample size, which varies with each sample
                        sample_size(ntrans) = size(cufrs{ntrans},2);
                        %% Adjust the random selection indicies for sample permutation sub sampling
                        rand_samp_ind(:,ntrans,1,:) = randi([1,size(cufrs{ntrans},2)],test_sample_size,1,1,niter);
                        rand_samp_ind(:,ntrans,2,:) = randi([1,size(cufrs{ntrans},2)],test_sample_size,1,1,niter);

                        perm_stat(:,u,p3,ntrans,iter) = nanmean(cat(2,cufrs{p3}(:,rand_samp_ind(:,p3,1,iter)),...
                                                                       cufrs{ntrans}(:,rand_samp_ind(:,ntrans,2,iter))),2);
                        diff_stat(:,u,ntrans,p3,iter) =  nanmean(cufrs{p3}(:,rand_samp_ind(:,p3,1,iter)),2)...
                                                           -nanmean(cufrs{ntrans}(:,rand_samp_ind(:,ntrans,2,iter)),2);
                    end
                end
            end

        end
    end
    toc
end




%% Calculating random permuted differences for obvStates VS surState
count = 1;
tic
for p1 = 1:ntrans,
    for p2 = p1:ntrans,
        count = 1;
        if p1==p2|p2==ntrans,
            for i = 1:round(sqrt(niter*2)),
                for j = i:round(sqrt(niter*2)),
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
wrucl =  zeros(nbins,numUnits,ntrans,ntrans);
wrlcl =  zeros(nbins,numUnits,ntrans,ntrans);
nnsamp = zeros(nbins,numUnits,ntrans,ntrans);
tic
for u = 1:numUnits,
    for p1 = 1:ntrans,
        for p2 = 1:ntrans,
            for t = 1:nbins,
                wnind = sq(~isnan(wrsps(t,u,p1,p2,1:niter)));
                nwnan = sum(~wnind);
                nnsamp(t,u,p1,p2) = niter-nwnan; 
                if nwnan == niter, continue,end
                wrucl(t,u,p1,p2) = wrsps(t,u,p1,p2,ceil((niter-nwnan)*(1-psig)));
                wrlcl(t,u,p1,p2) = wrsps(t,u,p1,p2,ceil((niter-nwnan)*psig));
            end
        end
    end
end
toc


%% Calculating random permuted differences for obvStates VS surState
prankstat = zeros(numel(tbins),numUnits,ntrans-1);
hrankstat = zeros(numel(tbins),numUnits,ntrans-1);
tic
for u = 1:numUnits;
    for p1 = 1:ntrans,
        for t = 1:nbins,
            [prankstat(t,u,p1),hrankstat(t,u,p1)] = ranksum(sq(diff_stat(t,u,p1,p1,~isnan(diff_stat(t,u,p1,p1,:)))),sq(diff_stat(t,u,p1,ntrans,~isnan(diff_stat(t,u,p1,ntrans,:)))));
        end
    end
end
toc

auxdata.numUnits = numUnits;
auxdata.nbins = nbins;
auxdata.tbins = tbins;
auxdata.tbound = find(auxdata.tbins>-2000&auxdata.tbins<2000);
auxdata.bname = bname;
auxdata.ntrans = ntrans;
auxdata.newSampleRate = newSampleRate

%% metadata
data.filebase = repmat({Trial.filebase},1,numUnits);
data.clu = units';
data.el = Trial.map(units,2)';


%% data
data.idd = zeros(nbins,numUnits,ntrans,niter);
data.idd_p95 = zeros(nbins,numUnits,ntrans);
data.idd_p05 = zeros(nbins,numUnits,ntrans);
data.pdd = zeros(nbins,numUnits,ntrans,niter);
data.pdd_p95 = zeros(nbins,numUnits,ntrans);
data.pdd_p05 = zeros(nbins,numUnits,ntrans);
data.edd = zeros(nbins,numUnits,ntrans,niter);
data.edd_p95 = zeros(nbins,numUnits,ntrans);
data.edd_p05 = zeros(nbins,numUnits,ntrans);

data.midd = zeros(nbins,numUnits,ntrans);
data.mpdd = zeros(nbins,numUnits,ntrans);
data.medd = zeros(nbins,numUnits);

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
data.pdd_idd_snr = zeros(nbins,numUnits,ntrans);
tic
for u = 1:numUnits,
    for p1 = 1:ntrans,
        data.pdd_idd_snr(:,u,p1) = (nanmean(diff_stat(:,u,p1,p1,1:niter),5)-nanmean(diff_stat(:,u,ntrans,p1,1:niter),5))./(nanvar(diff_stat(:,u,p1,p1,1:niter),[],5)+nanvar(diff_stat(:,u,ntrans,p1,1:niter),[],5));
    end
end
toc
%% Calculate the SNR between the intra-distribution-difference and
%% the permuted-distribution-difference
data.pdd_edd_snr = zeros(nbins,numUnits,ntrans);
tic
for u = 1:numUnits,
    for p1 = 1:ntrans,
        data.pdd_edd_snr(:,u,p1) = (nanmean(diff_stat(:,u,ntrans,ntrans,1:niter),5)-nanmean(diff_stat(:,u,ntrans,p1,1:niter),5))./(nanvar(diff_stat(:,u,ntrans,ntrans,1:niter),[],5)+nanvar(diff_stat(:,u,ntrans,p1,1:niter),[],5));
    end
end
toc


%% PVAL of mean of edd vs dist of surrogate pdd
data.medd_pval = zeros(nbins,numUnits,ntrans);
medd_pval = sum((repmat(diff_stat(:,:,ntrans,ntrans,1:niter),[1,1,ntrans,ntrans,1])<repmat(nanmean(diff_stat(:,:,ntrans,:,1:niter),5),[1,1,ntrans,1,niter])),5)./nnsamp;

medd_pval(medd_pval==1) = medd_pval(medd_pval==1) - 1./nnsamp(medd_pval==1);
medd_pval(medd_pval==0) = medd_pval(medd_pval==0) + 1./nnsamp(medd_pval==0);

medd_pval = 0.5-abs(medd_pval-0.5);

tic
for u = 1:numUnits,
    for p1 = 1:ntrans,
        data.medd_pval(:,u,p1) = medd_pval(:,u,p1,p1);
    end
end
toc

data.medd_pval(data.medd_pval>0.5) =data.medd_pval(data.medd_pval>0.5)-1;
data.nnsamp = nnsamp;


data = StructArray(data,2);

save([Trial.spath.analysis Trial.filebase '.cufrccg.' ...
      num2str(parameters.thresh_rad) '_' num2str(parameters.klen) '.pfc_' pfcbhv '.sur_' surbhv ...
      '.mat'], 'data','auxdata','parameters','-v7.3')

case 'display'

%% Load Session if Session is not already a MTASession
if ~isa(Session,'MTASession'),
    s = MTASession(Session);
    Trial = MTATrial(s,{},'all');
    clear('s');
end


load([Trial.spath.analysis Trial.filebase '.cufrccg.' ...
      num2str(parameters.thresh_rad) '_' num2str(parameters.klen) '.pfc_' parameters.pfcbhv '.sur_' parameters.surbhv ...
      '.mat']);

cmap = ['bgr'];
figure,
u=1;
report = 0;
while u~=-1
clf
for i = 1:auxdata.ntrans-1,
subplot(auxdata.ntrans-1,1,i);
boundedline(auxdata.tbins,data(u).midd(:,i),abs([data(u).idd_p05(:,i),data(u).idd_p95(:,i)]-repmat(data(u).midd(:,i),1,2)),cmap(1))
boundedline(auxdata.tbins,data(u).medd(:,i),abs([data(u).edd_p05(:,i),data(u).edd_p95(:,i)]-repmat(data(u).medd(:,i),1,2)),cmap(3),'alpha')
%boundedline(auxdata.tbins,data(u).mpdd(:,i),abs([data(u).pdd_p05(:,i),data(u).pdd_p95(:,i)]-repmat(data(u).mpdd(:,i),1,2)),cmap(2),'alpha')
boundedline(auxdata.tbins,data(u).mpdd(:,auxdata.ntrans),abs([data(u).pdd_p05(:,auxdata.ntrans),data(u).pdd_p95(:,auxdata.ntrans)]-repmat(data(u).mpdd(:,auxdata.ntrans),1,2)),cmap(2),'alpha')

%plot(auxdata.tbins,wrucl(:,u,i,i),'.')
%plot(auxdata.tbins,wrlcl(:,u,i,i),'.')

yl = ylim;
oym = yl(1);
imbins = round(abs(oym)./6);
yl = yl.*[2,1];
yr = [oym:oym+imbins]-imbins;

isignif = repmat(log10(1./[data(u).idd_pdd_prs(:,i).*data(u).idd_pdd_hrs(:,i)]'),numel(yr),1);
isignif(isinf(isignif)) = 0;
isignif(isignif<log10(20)) = 0;
imagesc(auxdata.tbins,yr,isignif)

yr = yr(1)-numel(yr)-1:yr(1);
isignif = repmat(log10(abs(1./[data(u).medd_pval(:,i).*(data(u).medd_pval(:,i)<0.05)]')),numel(yr),1);
isignif(isinf(isignif)) = 0;
isignif(isignif<log10(20)) = 0;
imagesc(auxdata.tbins,yr,isignif)

Lines(auxdata.tbins(auxdata.tbound(1)),[],'k');
Lines(auxdata.tbins(auxdata.tbound(end)),[],'k');
ylim([yr(1),yl(2)])
xlim([auxdata.tbins(1),auxdata.tbins(end)])
caxis([0,3])
colorbar
%Lines(auxdata.tbins(data(u).idd_ed_hrs
%imagecsc(
end
title(num2str(data(u).clu));


if report, 
    Trial.printFig([],[],data(u).clu,['cufrccg_' 'testing']);
    %                        num2str(parameters.thresh_rad) '_' num2str(parameters.klen)]);
    u = u+1;
    if u>auxdata.numUnits,
        break,
    end
else
    u = figure_controls(gcf,u);
end
end


case 'stats'



%% Load Session if Session is not already a MTASession
if ~isa(Session,'MTASession'),
    s = MTASession(Session);
    Trial = MTATrial(s,{},'all');
    clear('s');
end

load([Trial.spath.analysis Trial.filebase '.cufrccg.' ...
      num2str(parameters.thresh_rad) '.pfc_' parameters.pfcbhv '.sur_' parameters.surbhv ...
      '.mat']);

  
%% perm distr diff pval clusters
pddsc = cell(auxdata.numUnits,auxdata.ntrans-1);
mpddsc = cell(auxdata.numUnits,auxdata.ntrans-1);
mmpddsc = nan(auxdata.numUnits,auxdata.ntrans-1);
for u = 1:auxdata.numUnits,
    for i = 1:auxdata.ntrans-1,
        pddsc{u,i} = IntersectRanges(auxdata.tbound([1,end]),ThreshCross(data(u).idd_pdd_hrs(:,i),0.5,0));
        if ~isempty(pddsc{u,i}),
            for j = 1:size(pddsc{u,i},1),
                mpddsc{u,i} = [];
                mpddsc{u,i}(j) = mean(data(u).idd_pdd_prs(pddsc{u,i}(j,1):pddsc{u,i}(j,1),i));
                mmpddsc(u,i) = min(mpddsc{u,i});
            end
        end
    end
end



%% perm distr diff pval clusters
mepsc = cell(auxdata.numUnits,auxdata.ntrans-1);
dmepsc = cell(auxdata.numUnits,auxdata.ntrans-1);
mepdsc = cell(auxdata.numUnits,auxdata.ntrans-1);
mpddsc = cell(auxdata.numUnits,auxdata.ntrans-1);
mmpddsc = nan(auxdata.numUnits,auxdata.ntrans-1);
smmepsc = nan(auxdata.numUnits,auxdata.ntrans-1);
gdmepsc = nan(auxdata.numUnits,auxdata.ntrans-1);
igmepsc = nan(auxdata.numUnits,auxdata.ntrans-1);
immepsc = nan(auxdata.numUnits,auxdata.ntrans-1);
tigmepsc = nan(auxdata.numUnits,auxdata.ntrans-1);
timmepsc = nan(auxdata.numUnits,auxdata.ntrans-1);
sgmepsc = nan(auxdata.numUnits,auxdata.ntrans-1);
for u = 1:auxdata.numUnits,
    for i = 1:auxdata.ntrans-1,
        mepd_hs = abs(data(u).medd_pval(:,i))<0.05;
        mepd_ss = sign(data(u).medd_pval(:,i));
        mepsc{u,i} = IntersectRanges(auxdata.tbound([1,end]),ThreshCross(mepd_hs,0.5,0));
        if ~isempty(mepsc{u,i}),
            dmepsc{u,i} = diff(mepsc{u,i},1,2);
            for j = 1:size(mepsc{u,i},1),
                mmepsc{u,i} = [];
                mmepsc{u,i}(j) = mean(abs(data(u).medd_pval(mepsc{u,i}(j,1):mepsc{u,i}(j,2),i)));
            end

            [gdmepsc(u,i),igmepsc(u,i)] = max(dmepsc{u,i});
            [mmmepsc(u,i),immepsc(u,i)] = min(mmepsc{u,i});

            %% cluster index of smallest pvalue and greatest cluster size
            timmepsc(u,i) = round(mean([mepsc{u,i}(immepsc(u,i),1),mepsc{u,i}(immepsc(u,i),2)]));
            tigmepsc(u,i) = round(mean([mepsc{u,i}(igmepsc(u,i),1),mepsc{u,i}(igmepsc(u,i),2)]));

            %% sign of smallest pvalue and greatest cluster size
            smmepsc(u,i) = mepd_ss(timmepsc(u,i));
            sgmepsc(u,i) = mepd_ss(tigmepsc(u,i));

            %% mean pvalue of smallest pvalue and greatest cluster size
            mgmpv(u,i) =  mmepsc{u,i}(igmepsc(u,i));
            mmmpv(u,i) =  mmepsc{u,i}(immepsc(u,i));
        end
    end
end


end
