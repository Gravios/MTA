% script
% name: MjgER2016F3V1_ccg_distXbhv.m

% ANALYZE Pairwise unit ccg 
% 1. find placefields' centers
% 2. compute distance matrix of placefield centers
% 3. find midpoint between placefields which have centers less than 20 cm apart
% 4. construct basis from midpoint to one field center
% 5. construct rotation matrix to transform trajectories from room to field pair frame of reference
% 6. compute derivative of primary axis from midpoint to field
% 7. convert transformed trajectories to polar coordinates
% 8. compute ccg asymmetry for angle-distance pairs or just time-distance 
%    separeted by travel direction


pfindex = 1;

% LOAD variables:
%  Trials
%  units
%  cluSessionMap
%  pitchReferenceTrial
%  FigDir
%  sessionListName
%  sessionList
%  states
%  numStates
%  interpParPfsp
%  interpParDfs
% LOAD helper functions:
%  reshape_eigen_vector()
MjgER2016_load_data();
configure_default_args();

% LOAD theta state placefields
pft = cf(@(T,u)  pfs_2d_theta(T,u,'overwrite',false),  Trials,units);
% LOAD drzbhv placefields
bfs   = cf(@(t,u)   compute_bhv_ratemaps(t,u),  Trials, units);

% LABEL rippels (this goes somewhere else)
af(@(s) label_ripples(s,[],s.rippleDetectionChannels), sessionList);
label_ripples(sessionList(26),[],sessionList(26).rippleDetectionChannels);

% $$$ [pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = ...
% $$$     req20180123_ver5(Trials,[],[],false,false);

% LOAD bhv erpPCA scores
% REQUIRED Vars
% pfindex     - Numeric[1x1]                   index of the drzbhv field
% pfd         - cellarray[T x p]{MTAApfs}      drz restricted bhv field
% units       - cellarray[1 x T]{NumericArray} list session of units
% unitSubsets - cellarray[1 x p]{NumericArray} concatenated list of units
% eigVec      - cellarray[1 x p]{NumericArray} eigenvectors
% validDims   - cellarray[1 x p]{NumericArray} valid eigenvector indicies
% LOAD variables - see script:  MjgER2016_load_bhv_erpPCA_scores.m
%     FSrC     - matrix[U x V](Numeric); zscored fscores
%     fsrcz    - matrix[U x V](Numeric); zscore of FSrC
%     rmaps    - matrix[D x S](numeric); rate maps corresponding to the valid eigenvector dims

% $$$ MjgER2016_load_bhv_erpPCA_scores();
% $$$ clear('clu','tlu','rind','D','LR','FSCFr','rsMean','rsStd','pfdShuffled','FSrM','FSrS',...
% $$$       'fsrsMean','fsrsStd','rmapsShuffledMean','rmapsShuffled');

[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfs, units);



% GET pft rate maps
pmaps = cf(@(p,u)  mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'),pft,units);
clu   = cf(@(p,u)  p.data.clu(:,ismember(p.data.clu,u),:), pft,units);
tlu   = cf(@(i,u)  repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
pmaps = cat(2, pmaps{:});
clu   = cat(2, clu{:});
tlu   = cat(2, tlu{:});
clu   = [tlu',clu'];
[~,pind] = sortrows(clu);
pmaps = pmaps(:,pind);
pmaps = pmaps(nniz(pmaps),unitSubset);

bmaps = cf(@(p,u)  mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'),bfs,units);
clu   = cf(@(p,u)  p.data.clu(:,ismember(p.data.clu,u),:), pft,units);
tlu   = cf(@(i,u)  repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
bmaps = cat(2, bmaps{:});
clu   = cat(2, clu{:});
tlu   = cat(2, tlu{:});
clu   = [tlu',clu'];
[~,pind] = sortrows(clu);
bmaps = bmaps(:,pind);
bmaps = bmaps(validDims,unitSubset);
bmaps(isnan(bmaps)) = 0;


% INITIALIZE variables for collecting neigboring placefields
cluSessionMapSubset = cluSessionMap(unitSubset,:);
pairUids   = [];
pairClus   = [];
pairSDist  = [];
pairBDist  = [];
pairMmAmp  = [];
pairPkAmp  = [];
pairPkTime = [];
pairEigs   = [];
pairFSrC   = [];
pairMDSs   = [];

% INITIALIZE parameters for collecting neigboring placefields
pairCnt = 1;
unitCumCnt = 0;
epochs = [];
binSize = 1;
halfBins = 35;
normalization = 'hz';
eds = linspace(-pi,pi,2);
filt_fun = @(x) RectFilter(x,5,3);
mccg = zeros([2*halfBins+1,10000,2]);
distThreshold = 400;
sampleRate = 250;
state = 'theta-groom-sit';
state = 'R-theta';

for tind = 17:24
    %for tind = 1:numel(Trials),    
    disp(['MjgER2016F3V1_ccg_distXbhv_fine:tind:',num2str(tind)]);
    tic
    xyz = preproc_xyz(Trials{tind},'trb');
    xyz.resample(sampleRate);
    %xyz.filter('RectFilter');
    xyz.filter('ButFilter',5,1,'low');
    vxy = xyz.vel('spine_lower');
    feature = xyz.copy();
    feature.data = sq(feature(:,'hcom',[1,2]));
    %feature.data = sq(diff(feature(:,{'hcom','nose'},[1,2]),1,2));
    spk = Trials{tind}.spk.copy();
    spk = create(spk, Trials{tind}, xyz.sampleRate, state, units{tind}, ''); 
    


    unitsOfTrial = cluSessionMapSubset(cluSessionMapSubset(:,1)==tind,2);
    ddz = compute_ddz(Trials{tind},unitsOfTrial',pft{tind},'sampleRate',sampleRate);
    
    [mrate,mpos] = pft{tind}.maxRate(unitsOfTrial,true,'mean',0.9);

    nunits = numel(unitsOfTrial);

    for i=1:nunits-1,
        for j = i+1:nunits,
            
            d = pdist2(mpos(i,:),mpos(j,:),'euclidean');
            if d<=distThreshold,

                pairUids  = cat(1,pairUids,[i,j]+unitCumCnt);
                pairClus = cat(1,pairClus,[tind,unitsOfTrial([i,j])']);
                pairSDist = cat(1,pairSDist,d);
                pairEigs = cat(1,pairEigs,permute(eigScrs(pairUids(end,:),1:3),[3,1,2]));
                %pairMDSs = cat(1,pairMDSs,permute(mapa(pairUids(end,:),:),[3,1,2]));
                pairBDist = cat(1,pairBDist,sqrt(sum(diff(sq(pairEigs(end,:,:))).^2)));
                midpoint = sum(mpos([i,j],:))./2;
% COMPUTE the bais which represents the midpoint pair of placefields with major axis in the direction of second field
                pfsPairBasis = mpos(j,:)-midpoint;
                pfsPairBasis = pfsPairBasis./sqrt(sum(pfsPairBasis.^2));
                pfsPairBasis = [pfsPairBasis',pfsPairBasis([2,1])'];
% ROTATE traj coordinates
                pfhxy = multiprod(feature.data,pfsPairBasis,[2],[1,2]);
                pfhxy = cat(2,permute(pfhxy,[1,3,2]),circshift(permute(pfhxy,[1,3,2]),round(feature.sampleRate/5)));
% COMPUTE derivative of trajectory in pfs reference frame
                dpfhxy = sq(diff(pfhxy,1,2));
                pcor = cell([1,2]);  [pcor{:}] = cart2pol(dpfhxy(:,1),dpfhxy(:,2));
                %pcor = cell([1,2]);  [pcor{:}] = cart2pol(pfhxy(:,1),pfhxy(:,2));                
                th = circ_dist(pcor{1},-pi/2);
% $$$                 pfhxyH = multiprod(xyz(:,{'hcom','head_front'},:),pfsPairBasis,[2],[1,2]);
% COMPUTE derivative of trajectory in pfs reference frame
% $$$                 pcorH = cell([1,2]);
% $$$                 [pcorH{:}] = cart2pol(pfhxyH(:,1),pfhxyH(:,2));
% $$$                 thH = pcorH{1};
% SELECT spikes of units
                ii = unitsOfTrial(i)==unitsOfTrial;
                jj = unitsOfTrial(j)==unitsOfTrial;
                iRes = spk(unitsOfTrial(i));
                jRes = spk(unitsOfTrial(j));
% RESTRICT to spikes within distThreshold
                iRes = iRes(abs(ddz(iRes,i))<=distThreshold);
                jRes = jRes(abs(ddz(jRes,j))<=distThreshold);
% SPLIT the trajectories into paralel or perpendicular travel relative to place field pair
                for b = 1:numel(eds)-1,
                    grind = eds(b) <= th(iRes,1) & th(iRes,1) <= eds(b+1);
                    grjnd = eds(b) <= th(jRes,1) & th(jRes,1) <= eds(b+1);
% COMPUTE CCG of units for trajectories moving in specified direction
                    if sum(grind)&sum(grjnd),
                        [tccg,tbin] = CCG([iRes(grind);jRes(grjnd)],...
                                          [ones([sum(grind),1]);2*ones([sum(grjnd),1])],...
                                          binSize,halfBins,spk.sampleRate,[1,2],normalization);
                    else
                        tccg = zeros([halfBins*2+1,2,2]);
                    end
                    mccg(:,pairCnt,b) = tccg(:,1,2);
                end%for b

% ITERATE pair counter
                pairCnt = pairCnt + 1;
            end%if dist
        end%for j
    end%for i
     
% ADD nunits to cumulative unit counter
    unitCumCnt = unitCumCnt + nunits;

    toc();%Trials
end%for Trials



mccg(:,[size(pairUids,1)+1]:end,:) = [];

[pairPkAmp,pairPkInd] = max(cat(3,filt_fun(mccg(:,:,1)),...
                                  filt_fun(mccg(:,:,2))));

pairPkAmp  = sq(pairPkAmp);  
pairPkInd  = sq(pairPkInd);
pairMmAmp  = max(mean(cat(3,mccg(1:halfBins,:,1),mccg(halfBins+2:end,:,1),mccg(1:halfBins,:,2),mccg(halfBins+2:end,:,2))),[],3)';
pairPkTime = tbin(pairPkInd);

pairDCC = [];
for u = 1:size(pairUids,1),
    pairDCC(u,:,:) = corrcoef(pmaps(:,pairUids(u,:)));
end
pairDCC = pairDCC(:,1,2);

pairBCC = [];
for u = 1:size(pairUids,1),
    pairBCC(u,:,:) = corrcoef(bmaps(:,pairUids(u,:)));
end
pairBCC = pairBCC(:,1,2);


% CHECKPOINT --------------------------------------------
switch state
  case 'theta-groom-sit',
    datapath = fullfile(Trials{1}.path.project,'analysis','ccg_distXbhv_fine_theta.mat');
  case 'R-theta',
    datapath = fullfile(Trials{1}.path.project,'analysis','ccg_distXbhv_fine_SPW.mat');
end

save(datapath,...    
'binSize',...
'halfBins',...
'normalization',...
'eds',...
'filt_fun',...
'distThreshold',...
'sampleRate',...
'pairUids',...
'pairClus',...
'pairSDist',...
'pairBDist',...
'pairMmAmp',...
'pairPkAmp',...
'pairPkTime',...
'pairEigs',...
'pairFSrC',...
'pairMDSs',...
'mccg',...
'tbin',...
'pairPkAmp',...
'pairPkInd',...
'pairMmAmp',...
'pairPkTime',...
'pairDCC',...
'pairBCC');
% -------------------------------------------------------



load(datapath);



edy = linspace(-1,1,8);
edx = linspace(-1,1,8);

indx = discretize(pairDCC,edy);
indy = discretize(pairBCC,edy);

nxbins = numel(edx)-1;
nybins = numel(edy)-1;

% COMPUTE Sum of squares for 
% requires:
%   mean ccg magnitude of each bin
condExpCcgsFSdistXBdist = nan([nxbins,nybins]);
condExpCcgtFSdistXBdist = nan([nxbins,nybins]);
condExpCcgsFSdistXBstd = nan([nxbins,nybins]);
condExpCcgtFSdistXBstd = nan([nxbins,nybins]);
for x = 1:nxbins,
    for y = 1:nybins,
        xind = indx==x&indy==y;
        if sum(xind)>1,
            [mpk,mpki] = max(pairPkAmp(xind,:),[],2);
            nind = nniz(mpk)&mpk>1;
            tshift = sum(abs([tbin(pairPkInd(xind,1))',tbin(pairPkInd(xind,2))']).*~[mpki-1,mpki-2],2);
            condExpCcgtFSdistXBdist(x,y) = mean(tshift(nind));
            condExpCcgtFSdistXBstd(x,y)  = std(tshift(nind));            
            condExpCcgsFSdistXBdist(x,y) = mean(mpk);            
            condExpCcgsFSdistXBstd(x,y)  = std(mpk);                        
        end
    end
    xind = indx==x;
end


%figure();

clf();
subplot(221);
    hax = gca();
    imagesc(edx,edy,condExpCcgsFSdistXBdist');
    ca = colorbar();
    axis('xy');
    title({'Expected peak CCG Conditioned on','Spatial and Behaivoral Correlations'});
    xlabel('Spatial Correlation');
    ylabel('Behavioral Correlation');
    ylabel(ca,'mean peak ccg (Hz)');
    hax.Units = 'centimeters';
    hax.Position = [hax.Position(1:2),4,4];
    caxis([0,max(caxis)]);
subplot(222);
    hax = gca();
    imagesc(edx,edy,(condExpCcgtFSdistXBdist'));
    ca = colorbar();
    axis('xy');
    title({'Expected Time Lag CCG Conditioned on','Spatial and Behaivoral Correlations'});
    xlabel('Spatial Correlation');
    ylabel('Behavioral Correlation');
    ylabel(ca,'mean peak lag (ms)');
    hax.Units = 'centimeters';
    hax.Position = [hax.Position(1:2),4,4];
    caxis([0,max(caxis)]);
subplot(223);
    hax = gca();
    imagesc(edx,edy,condExpCcgsFSdistXBstd');
    ca = colorbar();
    axis('xy');
    title({'Expected peak CCG Conditioned on','Spatial and Behaivoral Correlations'});
    xlabel('Spatial Correlation');
    ylabel('Behavioral Correlation');
    ylabel(ca,'mean peak ccg (Hz)');
    hax.Units = 'centimeters';
    hax.Position = [hax.Position(1:2),4,4];
    caxis([0,max(caxis)]);
subplot(224);
    hax = gca();
    imagesc(edx,edy,(condExpCcgtFSdistXBstd'));
    ca = colorbar();
    axis('xy');
    title({'Expected Time Lag CCG Conditioned on','Spatial and Behaivoral Correlations'});
    xlabel('Spatial Correlation');
    ylabel('Behavioral Correlation');
    ylabel(ca,'mean peak lag (ms)');
    hax.Units = 'centimeters';
    hax.Position = [hax.Position(1:2),4,4];
    caxis([0,max(caxis)]);
colormap('jet');

print(gcf,'-dpng',fullfile('/storage/share/Projects/BehaviorPlaceCode/ccg',...
                             'ccg_bhv_short_timescale.png'));
print(gcf,'-depsc2',fullfile('/storage/share/Projects/BehaviorPlaceCode/ccg',...
                             'ccg_bhv_short_timescale.eps'));





% COMPUTE multi variate linear regression

ind = pairMmAmp~=0;
ind = max(pairPkAmp(ind,:),[],2) > 0;
ind = ':';

ind = max(pairPkAmp(ind,:),[],2) > 0;

xlabel('DistCC');ylabel('BhvCC');zlabel('max(CCG)');

[beta,Sigma,E,CovB,logL] = mvregress([pairDCC(ind),pairCC(ind)],log10(median(pairPkAmp(ind,:),2)+1));

[beta,Sigma,E,CovB,logL] = mvregress([pairSDist(ind)/400,pairCC(ind)],median(pairPkAmp(ind,:),2));


figure,plot3(pairSDist(ind)./400,pairCC(ind),log10(median(pairPkAmp(ind,:),2)+1),'.');
figure,plot3(pairSDist(ind)./400,pairCC(ind),median(pairPkAmp(ind,:),2),'.');

[beta,Sigma,E,CovB,logL] = mvregress([pairDCC(ind),pairCC(ind)],log10(median(pairPkAmp(ind,:),2)+1));
[beta,Sigma,E,CovB,logL] = mvregress([pairDCC(ind),pairCC(ind)],mean(pairPkAmp(ind,:),2));


beta
CovB
for j = 1:2,  [beta(j)/sqrt(CovB(j,j)), tcdf(beta(j)/sqrt(CovB(j,j)),numel(E)-2,'upper')],  end


