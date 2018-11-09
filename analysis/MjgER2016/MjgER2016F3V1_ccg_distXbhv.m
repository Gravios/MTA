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

% LOAD theta state placefields
pft = cf(@(T,u)  pfs_2d_theta(T,u,'overwrite',false),  Trials,units);
% LOAD drzbhv placefields
[pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = ...
    req20180123_ver5(Trials,[],[],false,false);

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
MjgER2016_load_bhv_erpPCA_scores();
clear('clu','tlu','rind','D','LR','FSCFr','rsMean','rsStd','pfdShuffled','FSrM','FSrS',...
      'fsrsMean','fsrsStd','rmapsShuffledMean','rmapsShuffled');


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
pmaps = pmaps(nniz(pmaps),unitSubsets{pfindex});



% INITIALIZE variables for collecting neigboring placefields
cluSessionMapSubset = cluSessionMap(unitSubsets{1},:);
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
binSize = 10;
halfBins = 50;
normalization = 'hz';
eds = linspace(-pi,pi,3);
filt_fun = @(x) RectFilter(x,11,3);
mccg = zeros([2*halfBins+1,10000,2]);
distThreshold = 400;


for tind = 17:numel(Trials),
    disp(['MjgER2016F3V1_ccg_distXbhv:tind:',num2str(tind)]);
    tic
    xyz = preproc_xyz(Trials{tind},'trb');
    xyz.filter('RectFilter');
    xyz.filter('ButFilter',5,1,'low');
    vxy = xyz.vel('spine_lower');
    feature = xyz.copy();
    feature.data = sq(feature(:,'nose',[1,2]));
    spk = Trials{tind}.spk.copy();
    spk = create(spk,Trials{tind},xyz.sampleRate,'theta-groom-sit',units{tind},'deburst'); 
    

    unitSubset = cluSessionMapSubset(cluSessionMapSubset(:,1)==tind,2);
    ddz = compute_ddz(Trials{tind},unitSubset',pft{tind});
    
    [mrate,mpos] = pft{tind}.maxRate(unitSubset,true,'mean',0.9);

    nunits = numel(unitSubset);

    for i=1:nunits-1,
       for j = i+1:nunits,
           d = pdist2(mpos(i,:),mpos(j,:),'euclidean');
           if d<=distThreshold,
               keyboard();
               pairUids  = cat(1,pairUids,[i,j]+unitCumCnt);
               pairClus = cat(1,pairClus,[tind,unitSubset([i,j])']);
               pairSDist = cat(1,pairSDist,d);
               pairEigs = cat(1,pairEigs,permute(eigScore{1}(pairUids(end,:),1:3),[3,1,2]));
               pairFSrC = cat(1,pairFSrC,permute(FSrC(pairUids(end,:),1:3),[3,1,2]));
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
               th = pcor{1};
               pfhxyH = multiprod(xyz(:,{'hcom','head_front'},:),pfsPairBasis,[2],[1,2]);

% COMPUTE derivative of trajectory in pfs reference frame
               pcorH = cell([1,2]);
               [pcorH{:}] = cart2pol(pfhxyH(:,1),pfhxyH(:,2));
               thH = pcorH{1};

% SELECT spikes of units
               ii = unitSubset(i)==unitSubset;
               jj = unitSubset(j)==unitSubset;
               iRes = spk(unitSubset(i));
               jRes = spk(unitSubset(j));
% RESTRICT to spikes within distThreshold
               iRes = iRes(abs(ddz(iRes,i))<=distThreshold);
               jRes = jRes(abs(ddz(jRes,j))<=distThreshold);

% SPLIT the trajectories into paralel or perpendicular travel relative to place field pair
               for b = 1:numel(eds)-1,
                   grind = eds(b) <= th(iRes,1) & th(iRes,1) <= eds(b+1);
                   grjnd = eds(b) <= th(jRes,1) & th(jRes,1) <= eds(b+1);

% COMPUTE CCG of units for trajectories moving in specified direction
                    if sum(grind)&sum(grjnd),
                        [tccg,tbin,pairs] = CCG([iRes(grind);jRes(grjnd)],...
                                          [ones([sum(grind),1]);2*ones([sum(grjnd),1])],...
                                          binSize,halfBins,spk.sampleRate,[1,2],normalization);
                    else
                        tccg = zeros([halfBins*2+1,2,2]);
                    end
                    mccg(:,pairCnt,b) = filt_fun(tccg(:,1,2));
                end
% ITERATE pair counter
                pairCnt = pairCnt + 1;
            end
        end
    end
     
% ADD nunits to cumulative unit counter
    unitCumCnt = unitCumCnt + nunits;

    toc();%Trials
end%for Trials



mccg(:,[size(pairUids,1)+1]:end,:) = [];

[pairPkAmp,pairPkInd] = max(cat(3,mccg(1:halfBins,:,1),...
                                  mccg(halfBins+2:end,:,1),...
                                  mccg(1:halfBins,:,2),...
                                  mccg(halfBins+2:end,:,2)));
pairPkAmp = sq(pairPkAmp);  pairPkInd = sq(pairPkInd);

pairMmAmp = max(mean(cat(3,mccg(1:halfBins,:,1),mccg(halfBins+2:end,:,1),mccg(1:halfBins,:,2),mccg(halfBins+2:end,:,2))),[],3)';
pairPkTime = tbin(pairPkInd);

pairCC = [];
for u = 1:size(pairUids,1),
    pairCC(u,:,:) = corrcoef(rmaps(:,pairUids(u,:)));
end
pairCC = pairCC(:,1,2);

pairDCC = [];
for u = 1:size(pairUids,1),
    pairDCC(u,:,:) = corrcoef(pmaps(:,pairUids(u,:)));
end
pairDCC = pairDCC(:,1,2);







edy = linspace(-1,1,7);
edx = linspace(-1,1,7);
%edx = linspace(0,40,7);
indx = discretize(pairDCC,edy);
%indx = discretize(pairSDist/10,edx);
indy = discretize(pairCC,edy);

nxbins = numel(edy)-1;
%nxbins = numel(edx)-1;
nybins = numel(edy)-1;

% COMPUTE Sum of squares for 
% requires:
%   mean ccg magnitude of each bin

condExpCcgsFSdistXBdist = nan([nxbins,nybins]);
%rssModel = nan([nxbins,nybins]);
%rssModelCol = nan([nxbins,1]);
for x = 1:nxbins,
    for y = 1:nybins,
        xind = indx==x&indy==y;

        if sum(xind)>1,
            %condExpCcgsFSdistXBdist(x,y) = mean(pairMmAmp(xind));            
            condExpCcgsFSdistXBdist(x,y) = mean(max(pairPkAmp(xind,:),[],2));
            %condExpCcgsFSdistXBdist(x,y) = mean(median(pairPkAmp(xind,:),2));
            %rssModel(x,y) = sum((pairMmAmp(xind)-condExpCcgsFSdistXBdist(x,y)).^2);
        end
    end
    xind = indx==x;
    %rssModelCol(x) = sum((pairMmAmp(xind)-condExpCcgsFSdistXBdist(sub2ind(size(condExpCcgsFSdistXBdist),repmat(x,[sum(xind),1]),indy(xind)))).^2,'omitnan');
end


figure();
hax = gca();
imagesc(edx,edy,condExpCcgsFSdistXBdist');
ca = colorbar();
axis('xy');
title({'Expected peak CCG Conditioned on','Spatial and Behaivoral Correlations'});
xlabel('Spatial Correlation');
ylabel('Behavioral Correlation');
ylabel(ca,'mean peak ccg (Hz)');
hax.Units = 'centimeters';
hax.Position = [hax.Position(1:2),6,6];


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


% RASTER STUFF
Trial = Trials{20};    % jg05-20120312.cof.all
unts  = units{20};     % good units
xyz = preproc_xyz(Trial,'trb');


spk = Trial.spk.copy();
spk = create(spk,Trial,Trial.lfp.sampleRate,'theta-groom-sit',unts,'deburst');


 
unitSubset = cluSessionMapSubset(cluSessionMapSubset(:,1)==20,2);
unitBhvScores  = fsrcz(cluSessionMapSubset(:,1)==20,1:3);


figure();
subplot(3,1,[1,2]);
hold('on');
unitSet = unts;
unitClr = repmat('b',[1,numel(unts)]);
for u = unitSet
    uind = find(u==unitSet);
    res = spk(u);
    ts = 1:max(res);
    sWidth = 0.4;
    sti = ismember(ts,res);
    patch(reshape(repmat(reshape(bsxfun(@plus,ts(sti),[-sWidth/2;sWidth/2]),[],1)',2,1),[],1),...
          repmat([0,1,1,0],[1,numel(res)])'+uind,...
          unitClr(uind),'EdgeColor',unitClr(uind));
end
subplot(3,1,3);
plotSTC(Trial.stc,Trial.lfp.sampleRate,[],fliplr({'rear','hloc','hpause','lloc','lpause'}),fliplr('rbcgk'));
linkaxes(get(gcf,'Children'),'x');