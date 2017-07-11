%rear_ufr.m
%Precursor to CorrectedUFRCCG.m


s = MTASession('jg05-20120317');

Trial = MTATrial(s,{{'CluRes',s.xyzSampleRate}},'all');

Trial = Trial.filter();

klen = 64;
kern = ones(klen,1);
overlap = 16;
nnn = 81;
nxbins = 50;
nybins = 50;
xbins = linspace(Trial.Maze.boundaries(1,1),Trial.Maze.boundaries(1,2),nxbins)';
ybins = linspace(Trial.Maze.boundaries(2,2),Trial.Maze.boundaries(2,1),nxbins)';


myxyz = sq(Trial.xyz(:,7,[1,2]));
t =         permute(reshape(          myxyz(1:size(myxyz,1)-mod(size(myxyz,1),klen),:),klen,[],2),[4,1,2,3]);
for shift = 1:klen/overlap-1,
t = cat(1,t,permute(reshape(circshift(myxyz(1:size(myxyz,1)-mod(size(myxyz,1),klen),:),-overlap*shift),klen,[],2),[4,1,2,3]));
end
myxyz = t;
myxyz = reshape(sq(sum(repmat(permute(repmat(permute(repmat(kern./sum(kern),1,size(myxyz,3)),[5,4,1,2,3]),size(myxyz,4),1),[2,3,4,1]),klen/overlap,1).*myxyz,2)),[],2);

newSampleRate = 1/((size(Trial.xyz,1)-mod(size(Trial.xyz,1),klen))/Trial.xyzSampleRate/length(myxyz));

myrx = SelectPeriods(myxyz,round((Trial.Bhv.getState('rear').state+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),'c',1,1);
mywx = SelectPeriods(myxyz,round((Trial.Bhv.getState('walk').state+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),'c',1,1);




% $$$ ufrs = {};
% $$$ %%iu = [12,14,24];
% $$$ while unit ~=-1
myufr = {};
for unit = 1:76

myres = Trial.res(Trial.clu==unit);
myufr{unit} = zeros(size(Trial.xyz,1),1);
myufr{unit}(1:myres(end)) = accumarray(myres,ones(size(myres),1));
t =         permute(reshape(myufr{unit}(1:size(myufr{unit},1)-mod(size(myufr{unit},1),klen)),klen,[]),[3,1,2]);
for shift = 1:klen/overlap-1
t = cat(1,t,permute(reshape(circshift(myufr{unit}(1:size(myufr{unit},1)-mod(size(myufr{unit},1),klen)),-overlap*shift),klen,[]),[3,1,2]));
end
myufr{unit} = t; 
myufr{unit} = reshape(sq(sum(repmat(permute(repmat(kern,1,size(myufr{unit},3)),[3,1,2]),size(myufr{unit},1),1).*myufr{unit},2))/(klen/Trial.xyzSampleRate),[],1);
end

myufr = cell2mat(myufr);

pfw = MTAPlaceField(Trial,[],'walk');

wpmr = zeros(size(myxyz,1),76);

for unit =1:76,
for i=1:size(myxyz,1),
wpmr(i,unit) = pfw.rateMap{unit}(find(pfw.xbin<round(myxyz(i,1)),1,'last'),find(pfw.ybin<round(myxyz(i,2)),1,'last'));
end
end

ufrwd = myufr-wpmr;






non_immobile_bhvs = {'walk','bturn','hturn','rear'};
iper = [1,size(Trial.xyz,1)];
for i = 1:length(non_immobile_bhvs),
iper = SubstractRanges(iper,Trial.Bhv.getState(non_immobile_bhvs(i)).state);
end
Trial.Bhv = Trial.Bhv.addState('i','immobile',iper);
iperdur = diff(Trial.Bhv.getState('immobile').state,1,2);



%%Filtered behaviors 
rper = Trial.Bhv.getState('rear').state;
interrdur = [rper(2:end,1)-rper(1:end-1,2)];
grperind = find(interrdur>3*Trial.xyzSampleRate);
%rearing offset: filtered out any with rearing offsets preceding by
%3 seconds
no_fol_rear_rearoff = rper(grperind,2); 
%rearing offset: filtered out any with rearing offsets following by
%3 seconds
no_pre_rear_rearon =  rper([1;grperind+1],1); 




wper = Trial.Bhv.getState('walk').state;
interwdur = [wper(2:end,1)-wper(1:end-1,2)];
gwperind = find(interwdur>3*Trial.xyzSampleRate);
%walking offset: filtered out any with walking offsets preceding by
%3 seconds
no_fol_walk_walkoff = wper(gwperind,2); 
%walking offset: filtered out any with walking offsets following by
%3 seconds
no_pre_walk_walkon =  wper([1;gwperind+1],1); 


[~,nfroffi,npwoni] = NearestNeighbour(no_fol_rear_rearoff,no_pre_walk_walkon);


interwrdur = abs(no_pre_walk_walkon(npwoni)-no_fol_rear_rearoff(nfroffi));
bwrperind = find(interwrdur<3*Trial.xyzSampleRate);
badr = unique(nfroffi(bwrperind));
no_fol_rearwalk_rearoff = no_fol_rear_rearoff;
no_fol_rearwalk_rearoff(badr) = [];
ys_fol_rearwalk_rearoff =no_fol_rear_rearoff(badr);


[~,nfroni,npwoffi] = NearestNeighbour(no_pre_rear_rearon,no_fol_walk_walkoff);

interwrdur = abs(no_fol_walk_walkoff(npwoffi)-no_pre_rear_rearon(nfroni));
bwrperind = find(interwrdur<3*Trial.xyzSampleRate);
badr = unique(nfroni(bwrperind));
no_pre_rearwalk_rearon = no_pre_rear_rearon;
no_pre_rearwalk_rearon(badr) = [];
ys_pre_rearwalk_rearon =no_pre_rear_rearon(badr);


btper = Trial.Bhv.getState('bturn').state;

no_pre_rear_bturn =[];
for i=1:size(btper,1),
if isempty(find(abs(repmat(btper(i,1),size(rper,1),1)-rper(:,1))<2*Trial.xyzSampleRate,1,'first'))
no_pre_rear_bturn(end+1,:) = btper(i,:);
end
end

no_pre_rear_bturnon = no_pre_rear_bturn(:,1);
no_fol_rear_bturnoff = no_pre_rear_bturn(:,2);


no_pre_rearturn_bturn =[];
for i=1:size(no_pre_rear_bturn,1),
if isempty(find(abs(repmat(no_pre_rear_bturn(i,1),size(no_pre_rear_bturn,1),1)-no_pre_rear_bturn(:,1))<2*Trial.xyzSampleRate&abs(repmat(no_pre_rear_bturn(i,1),size(no_pre_rear_bturn,1),1)-no_pre_rear_bturn(:,1))>0,1,'first'))
no_pre_rearturn_bturn(end+1,:) = no_pre_rear_bturn(i,:);
end
end

no_pre_rearturn_bturnon = no_pre_rearturn_bturn(:,1);
no_fol_rearturn_bturnoff = no_pre_rearturn_bturn(:,2);



group_resdescrip{1} = 'no_pre_rear_rearon';
group_restrains{1} =   no_pre_rear_rearon;

group_resdescrip{2} = 'no_fol_rear_rearoff';
group_restrains{2} =   no_fol_rear_rearoff;

group_resdescrip{3} = 'no_pre_rearwalk_rearon';
group_restrains{3} =   no_pre_rearwalk_rearon;

group_resdescrip{4} = 'no_fol_rearwalk_rearoff';
group_restrains{4} =   no_fol_rearwalk_rearoff;

group_resdescrip{5} = 'ys_pre_rearwalk_rearon';
group_restrains{5} =   ys_pre_rearwalk_rearon;

group_resdescrip{6} = 'ys_fol_rearwalk_rearoff';
group_restrains{6} =   ys_fol_rearwalk_rearoff;

group_resdescrip{7} = 'no_pre_walk_walkon';
group_restrains{7} =   no_pre_walk_walkon;

group_resdescrip{8} = 'no_fol_walk_walkoff';
group_restrains{8} =   no_fol_walk_walkoff;

group_resdescrip{9} = 'no_pre_rear_bturnon';
group_restrains{9} =   no_pre_rear_bturnon;

group_resdescrip{10} = 'no_fol_rear_bturnoff';
group_restrains{10} =   no_fol_rear_bturnoff;

group_resdescrip{11} = 'no_pre_rearturn_bturnon';
group_restrains{11} =   no_pre_rearturn_bturnon;

group_resdescrip{12} = 'no_fol_rearturn_bturnoff';
group_restrains{12} =   no_fol_rearturn_bturnoff;

ronpw = group_restrains{1};
roffw = group_restrains{2};
wonprm = group_restrains{7};
woffprm = group_restrains{8};

rufrs = GetSegs(ufrwd,round(((ronpw+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));
dufrs = GetSegs(ufrwd,round(((roffw+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));
wufrs = GetSegs(ufrwd,round(((wonprm+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));
sufrs = GetSegs(ufrwd,round(((woffprm+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));

%ufrs = GetSegs(myufr,round(((ronpw+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));
%ufrf = GetSegs(myufr,round(((roffw+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));

pfw = MTAPlaceField(Trial,[],'walk');
pfr = MTAPlaceField(Trial,[],'rear');

pfMaxPosRear = zeros(size(pfr.cluMap,1),2);
pfMaxPosWalk = zeros(size(pfr.cluMap,1),2);

for unit = 1:size(pfw.cluMap,1),
    try, pfMaxPosWalk(unit,:) = pfw.maxRatePos{unit}(pfw.maxRateMax{unit},:);,end
    try, pfMaxPosRear(unit,:) = pfr.maxRatePos{unit}(pfr.maxRateMax{unit},:);,end
    upos(unit,:) = pfMaxPosWalk(unit,:);
    %upos(unit,:) = pfMaxPosRear(unit,:);
end

rxyzs = GetSegs(myxyz,round(((ronpw+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));
wxyzs = GetSegs(myxyz,round(((Trial.Bhv.getState('walk').state(:,1)+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate)-5*newSampleRate),round(10*newSampleRate));
rposon = myxyz(round((ronpw+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);
rposoff = myxyz(round((roffw+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);
wposon = myxyz(round((wonprm+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);
wposoff= myxyz(round((woffprm+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);

tbins = linspace(-5000,5000,size(rufrs,1));

% $$$ flen = 5
% $$$ unit = 29;
% $$$ while unit~=-1,
% $$$ rwdist = sqrt(sum((rposon-repmat(upos(unit,:),size(rposon,1),1)).^2,2));
% $$$ [srwdist,srwdind] = sort(rwdist);
% $$$ rwdisto = sqrt(sum((rposoff-repmat(upos(unit,:),size(rposoff,1),1)).^2,2));
% $$$ [srwdisto,srwdindo] = sort(rwdisto);
% $$$ 
% $$$ wwdist = sqrt(sum((wposon-repmat(upos(unit,:),size(wposon,1),1)).^2,2));
% $$$ [swwdist,swwdind] = sort(wwdist);
% $$$ wwdisto = sqrt(sum((wposoff-repmat(upos(unit,:),size(wposoff,1),1)).^2,2));
% $$$ [swwdisto,swwdindo] = sort(wwdisto);
% $$$ 
% $$$ 
% $$$ clf
% $$$ subplot2(4,5,[1,2],1);
% $$$ pfr.plot(unit,1);
% $$$ title(num2str(unit))
% $$$ subplot2(4,5,[3,4],1);
% $$$ pfw.plot(unit,1);
% $$$ 
% $$$ 
% $$$ bs(1) =subplot2(4,5,[1,2],[2,3]);
% $$$ hold on
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(rufrs(:,rwdist<200,unit),2)))
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(rufrs(:,rwdist>200&rwdist<400,unit),2)),'r')
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(rufrs(:,rwdist>400&rwdist<600,unit),2)),'g')
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(rufrs(:,rwdist>600&rwdist<800,unit),2)),'c')
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(rufrs(:,:,unit),2)),'k')
% $$$ Lines([],0,'k');
% $$$ bs(2) =subplot2(4,5,[3,4],[2,3]);
% $$$ hold on
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(dufrs(:,rwdisto<200,unit),2)))
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(dufrs(:,rwdisto>200&rwdisto<400,unit),2)),'r')
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(dufrs(:,rwdisto>400&rwdisto<600,unit),2)),'g')
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(dufrs(:,rwdisto>600&rwdisto<800,unit),2)),'c')
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(dufrs(:,:,unit),2)),'k')
% $$$ Lines([],0,'k');
% $$$ 
% $$$ bs(3) =subplot2(4,5,[1,2],[4,5]);
% $$$ hold on
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(wufrs(:,wwdist<200,unit),2)))
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(wufrs(:,wwdist>200&wwdist<400,unit),2)),'r')
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(wufrs(:,wwdist>400&wwdist<600,unit),2)),'g')
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(wufrs(:,wwdist>600&wwdist<800,unit),2)),'c')
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(wufrs(:,:,unit),2)),'k')
% $$$ Lines([],0,'k');
% $$$ bs(4) =subplot2(4,5,[3,4],[4,5]);
% $$$ hold on
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(sufrs(:,wwdisto<200,unit),2)))
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(sufrs(:,wwdisto>200&wwdisto<400,unit),2)),'r')
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(sufrs(:,wwdisto>400&wwdisto<600,unit),2)),'g')
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(sufrs(:,wwdisto>600&wwdisto<800,unit),2)),'c')
% $$$ plot(tbins,Filter0(gausswin(flen)./sum(gausswin(flen)),nanmean(sufrs(:,:,unit),2)),'k')
% $$$ Lines([],0,'k');
% $$$ linkaxes(bs,'xy');
% $$$ xlim([tbins(1),tbins(end)])
% $$$ title(num2str(unit))
% $$$ %subplot(313),
% $$$ %plot(rufrs(:,rwdist>500&rwdist<600,unit))
% $$$ %plot(rufrs(:,rwdist<200,unit))
% $$$ unit = figure_controls(gcf,unit)
% $$$ end








%% START - rear walk premutation stats
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
thresh_rad = 1000;
test_sample_size  = 20;
min_sample_size = 10;
niter = 1000;
ntrans = 5;
perm_stat = zeros(size(rufrs,1),size(rufrs,3),ntrans,ntrans,niter);
diff_stat = zeros(size(rufrs,1),size(rufrs,3),ntrans,ntrans,niter);


%% Generate surrogate walk segments from place field max rate pos

wsur = [];
for i = 1:size(wper,1)
wsur = cat(1,wsur,[wper(i,1):wper(i,2)]');
end





%% Generate random samples from the walk surrogates
wsss = size(wposon,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for unit = 1:size(rufrs,3),
fprintf('Starting unit: %i \n',unit)
%% Rear onset pfc distance
rwdist = sqrt(sum((rposon-repmat(upos(unit,:),size(rposon,1),1)).^2,2));
%% Rear offset pfc distance
rwdisto = sqrt(sum((rposoff-repmat(upos(unit,:),size(rposoff,1),1)).^2,2));
%% Walk onset pfc distance
wwdist = sqrt(sum((wposon-repmat(upos(unit,:),size(wposon,1),1)).^2,2));
%% Walk offset pfc distance
wwdisto = sqrt(sum((wposoff-repmat(upos(unit,:),size(wposoff,1),1)).^2,2));

iter = 1;
wsur_rand_ind = randi([1,length(wsur)],wsss,niter);
wsursamp = wsur(wsur_rand_ind(:,iter));
wsurpos = myxyz(round((wsursamp+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);
wsurdist = sqrt(sum((wsurpos-repmat(upos(unit,:),wsss,1)).^2,2));

%% resample if too few are found
timeout_countdown = 10;
while sum(wsurdist<thresh_rad)<test_sample_size & timeout_countdown ~= 0 ,
    %% Draw new random indicies for surrogate selection
    wsur_rand_ind(:,iter) = randi([1,length(wsur)],wsss,1);
    %% Pull random sample from the walk surrogate tarjectories
    wsursamp = wsur(wsur_rand_ind(:,iter));
    %% Get the position of the random samples from the walk surrogate tarjectories
    wsurpos = myxyz(round((wsursamp+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);
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
                   round(((wsursamp+0.5*Trial.xyzSampleRate/newSampleRate)...
                                     ./Trial.xyzSampleRate.*newSampleRate)...
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
                diff_stat(:,unit,p1,p2,iter) =  nanmean(cufrs{p1}(:,rand_samp_ind(:,p1,1,iter)),2)...
                    -nanmean(cufrs{p2}(:,rand_samp_ind(:,p2,2,iter)),2);
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
                    wsurpos = myxyz(round((wsursamp+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);
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
                        wsurpos = myxyz(round((wsursamp+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),:);
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
                                       round(((wsursamp+0.5*Trial.xyzSampleRate/newSampleRate)...
                                              ./Trial.xyzSampleRate.*newSampleRate)...
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
                    diff_stat(:,unit,p3,ntrans,iter) =  nanmean(cufrs{p3}(:,rand_samp_ind(:,p3,1,iter)),2)...
                                                   -nanmean(cufrs{ntrans}(:,rand_samp_ind(:,ntrans,2,iter)),2);
                end
            end
        end

    end
end
toc

end



% $$$ psig = .05
% $$$ sps = sort(perm_stat(:,:,5,5,:),5);
% $$$ ucl = sps(:,:,:,:,round(niter*(1-psig)));
% $$$ lcl = sps(:,:,:,:,round(niter*psig));



% $$$ 
% $$$ psig = .05
% $$$ sps = sort(perm_stat(:,:,5,5,:),5);
% $$$ ucl = max(sps(35:40,:,:,:,round(niter*(1-psig))));
% $$$ lcl = min(sps(35:40,:,:,:,round(niter*psig)));

% $$$ 
% $$$ psig = .05;
% $$$ sps = sort(diff_stat(:,:,5,5,1:niter),5);
% $$$ ucl = max(sps(35:40,:,:,:,round(niter*(1-psig))));
% $$$ lcl = min(sps(35:40,:,:,:,round(niter*psig)));


% $$$ 
% $$$ psig = .05;
% $$$ sps = sort(diff_stat(:,:,5,5,1:niter),5);
% $$$ ucl = sps(:,:,:,:,round(niter*(1-psig)));
% $$$ lcl = sps(:,:,:,:,round(niter*psig));


psig = .05;
wrsps = sort(diff_stat(:,:,:,:,1:niter),5);
wrucl = wrsps(:,:,:,:,round(niter*(1-psig)));
wrlcl = wrsps(:,:,:,:,round(niter*psig));

% $$$ 
% $$$ psig = .05;
% $$$ sps = sort(diff_stat(:,36,5,5,1:niter),5);
% $$$ ucl = max(sps(35:40,:,:,:,round(niter*(1-psig))));
% $$$ lcl = min(sps(35:40,:,:,:,round(niter*psig)));


%figure

% $$$ if ~exist('pfw','var')
% $$$ pfw = MTAPlaceField(Trial,[],'walk');
% $$$ pfr = MTAPlaceField(Trial,[],'rear');
% $$$ [accg,tbin] = autoccg(Trial);

unit=32;
while unit~=-1
rwdist = sqrt(sum((rposon-repmat(upos(unit,:),size(rposon,1),1)).^2,2));
%% Rear offset pfc distance
rwdisto = sqrt(sum((rposoff-repmat(upos(unit,:),size(rposoff,1),1)).^2,2));
%% Walk onset pfc distance
wwdist = sqrt(sum((wposon-repmat(upos(unit,:),size(wposon,1),1)).^2,2));
%% Walk offset pfc distance
wwdisto = sqrt(sum((wposoff-repmat(upos(unit,:),size(wposoff,1),1)).^2,2));
cufrs{1} = rufrs(:,rwdist<thresh_rad,unit);
cufrs{2} = dufrs(:,rwdisto<thresh_rad,unit);
cufrs{3} = wufrs(:,wwdist<thresh_rad,unit);
cufrs{4} = sufrs(:,wwdisto<thresh_rad,unit);
cufrs{5} = wsurufrs(:,wsurdist<thresh_rad);
clf
for x = 1:ntrans,
for y = x:ntrans,
sph(x,y) = subplotfit(x,y,[ntrans,ntrans]);
%plot(tbins,mean(perm_stat(:,unit,x,y,:),5))
plot(tbins,nanmean(diff_stat(:,unit,x,y,:),5))
hold on

if x==y,
plot(tbins,nanmean(diff_stat(:,unit,x,ntrans,:),5),'-r')
plot(tbins,wrucl(:,unit,x,5),'-r')
plot(tbins,wrlcl(:,unit,x,5),'-r')
plot(tbins,wrucl(:,unit,x,y),'-b')
plot(tbins,wrlcl(:,unit,x,y),'-b')
end
%plot(tbins,ucl(:,unit),'-r')
%plot(tbins,lcl(:,unit),'-r')
%Lines([],ucl(unit),'r');
%Lines([],lcl(unit),'r');
%plot(sq(perm_stat(:,unit,x,y,:)))
%imagesc(tbins,1:niter,sq(perm_stat(:,unit,x,y,:))'),axis xy
xlim([tbins(1),tbins(end)])
end
end
title(num2str(unit))
subplot2(5,5,[2,3],1);
pfr.plot(unit,1);
subplot2(5,5,[4,5],1);
pfw.plot(unit,1);
subplot2(5,5,[4,5],2);
bar(tbin,accg(:,unit)),axis tight


% $$$ pron = zeros(numel(tbins),1);
% $$$ hron = zeros(numel(tbins),1);
% $$$ for t = 1:numel(tbins),
% $$$ [pron(t),hron(t)] = ranksum(sq(diff_stat(t,unit,1,1,~isnan(diff_stat(t,unit,1,1,:)))),sq(diff_stat(t,unit,1,5,~isnan(diff_stat(t,unit,1,5,:)))));
% $$$ end
% $$$ plot(pron(:))


%xlim([1,numel(tbins)])
unit = figure_controls(gcf,unit)
end




%% calculate random permuted differences for obvStates VS surState
wsdd = zeros(numel(tbins),size(rufrs,3),ntrans,niter);
count = 1;
for k = 1:ntrans,
count = 1;
for i = 1:110,
for j = i:110,
if i==j,continue,end
wsdd(:,:,k,count) = perm_stat(:,:,k,5,i)-perm_stat(:,:,k,5,j);
if count == niter , break , end
count = count+1;
end
end
end

unit = 26;
pron = zeros(numel(tbins),1);
hron = zeros(numel(tbins),1);
hwr = zeros(numel(tbins),1);
pwr = zeros(numel(tbins),1);
for t = 1:numel(tbins),
[pron(t),hron(t)] = ranksum(sq(diff_stat(t,unit,1,1,~isnan(diff_stat(t,unit,1,1,:)))),sq(wsdd(t,unit,1,~isnan(wsdd(t,unit,1,:)))));
[pwsu(t),hwsu(t)] = ranksum(sq(diff_stat(t,unit,5,5,~isnan(diff_stat(t,unit,5,5,:)))),sq(wsdd(t,unit,1,~isnan(wsdd(t,unit,1,:)))));
[pwr(t),hwr(t)] = ranksum(sq(diff_stat(t,unit,5,5,~isnan(diff_stat(t,unit,5,5,:)))),sq(diff_stat(t,unit,1,1,~isnan(diff_stat(t,unit,1,1,:)))));
end


%mdsup = 1./sum(repmat(nanmedian(diff_stat(:,:,1,1,:),5),[1,1,1,1,niter])>diff_stat(:,:,1,5,:),5);
%uds = (diff_stat(:,:,1,5,:)-repmat(nanmean(diff_stat(:,:,1,5,:),5),[1,1,1,1,niter]))./repmat(nanstd(diff_stat(:,:,1,5,:),[],5),[1,1,1,1,niter]);
%figure,plot(sq(uds(:,32,1,1,:)))

figure,
plot(pwr(:))

figure,
plot(pron(:))
hold on
plot(pwsu,'g')

t = 36;
[countw,bin] = hist(sq(diff_stat(t,unit,5,5,1:niter)),100);
[countr,bin1] = hist(sq(wsdd(t,unit,1,1:niter)),100);
[countp,bin2] = hist(sq(diff_stat(t,unit,1,1,1:niter)),100);
figure, 
plot(bin,countw,bin1,countr,bin2,countp)


