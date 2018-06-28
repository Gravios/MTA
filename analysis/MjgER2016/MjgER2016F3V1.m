

% F3A1 unit ccg
% F3A2 

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



% F3A1 Pairwise unit ccg 
% 1. find placefields' centers
% 2. compute distance matrix of placefield centers
% 3. find midpoint between placefields which have centers less than 20 cm apart
% 4. construct basis from midpoint to one field center
% 5. construct rotation matrix to transform trajectories from room to field pair frame of reference
% 6. compute derivative of primary axis from midpoint to field
% 7. convert transformed trajectories to polar coordinates
% 8. compute ccg asymmetry for angle-distance pairs or just time-distance 
%    separeted by travel direction



% LOAD theta state placefields
pft = cf(@(T,u)  pfs_2d_theta(T,u,'overwrite',false),  Trials,units);
% LOAD drzbhv placefields
[pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = ...
    req20180123_ver5(Trials,[],[],false,false);


tind = 20;


% GET Behavior fields
xyz = preproc_xyz(Trials{tind},'trb');
xyz.filter('RectFilter');
xyz.filter('ButFilter',5,2.5,'low');
vxy = xyz.vel('spine_lower');
feature = xyz.copy();
feature.data = sq(feature(:,'hcom',[1,2]));
[mrate,mpos] = pft{tind}.maxRate(units{tind},true,'mean',0.9);
spk = Trials{tind}.spk.copy();
spk = create(spk,Trials{tind},xyz.sampleRate,'theta-groom-sit',units{tind},'deburst'); 

epochs = [];

binSize = 3;
halfBins = 40;

% $$$ binSize = 1;
% $$$ halfBins = 12;
% $$$ normalization = 'hz';

figure
eds = linspace(0,2*pi,2);
filt_fun = @(x) RectFilter(x,5,3);
%filt_fun = @(x) x;
mccg = [];sccg = [];occg = [];
nunits = numel(units{tind});
for i=1:nunits-1,
    for j = i+1:nunits,
        
        d = pdist2(mpos(i,:),mpos(j,:),'euclidean');
        if d<=200,
            
            midpoint = sum(mpos([i,j],:))./2;
            
            clf();
            subplot(431);
            plot(pft{tind},units{tind}(i),1,true,[],true);
            hold('on');
            plot(mpos(i,1),mpos(i,2),'*m');
            plot(midpoint(1),midpoint(2),'*g');            
            subplot(433);
            plot(pft{tind},units{tind}(j),1,true,[],true);
            hold('on');            
            plot(mpos(j,1),mpos(j,2),'*m');
            plot(midpoint(1),midpoint(2),'*g');



            subplot(434);
            plot(pfd{tind,1},units{tind}(i),1,true,[],false);            
            subplot(436);
            plot(pfd{tind,1},units{tind}(j),1,true,[],false);


            pfsPairBasis = mpos(j,:)-midpoint;
            pfsPairBasis = pfsPairBasis./sqrt(sum(pfsPairBasis.^2));
            pfsPairBasis = [pfsPairBasis',pfsPairBasis([2,1])'];
% ROTATE traj coordinates
            pfhxy = multiprod(feature.data,pfsPairBasis,[2],[1,2]);
            pfhxy = cat(2,permute(pfhxy,[1,3,2]),circshift(permute(pfhxy,[1,3,2]),round(feature.sampleRate/5)));
            % COMPUTE derivative of trajectory in pfs reference frame
            dpfhxy = sq(diff(pfhxy,1,2));
            pcor = cell([1,2]);
            [pcor{:}] = cart2pol(dpfhxy(:,1),dpfhxy(:,2));
            th = pcor{1}+pi+eds(2)/2;

            pfhxyH = multiprod(xyz(:,{'hcom','head_front'},:),pfsPairBasis,[2],[1,2]);
            % COMPUTE derivative of trajectory in pfs reference frame
            pcorH = cell([1,2]);
            [pcorH{:}] = cart2pol(pfhxyH(:,1),pfhxyH(:,2));
            thH = pcorH{1}+pi+eds(2)/2;
            
            ii = units{tind}(i)==units{tind};
            jj = units{tind}(j)==units{tind};

            iRes = spk(units{tind}(i));
%            iRes(vxy(iRes)<2) = [];
            jRes = spk(units{tind}(j));
%            jRes(vxy(jRes)<2) = [];
            % GET indicies within certain range
            % abs(x)+100
            
% $$$             grind = sign(dpfhxy(iRes,1))==1;
% $$$             grjnd = sign(dpfhxy(jRes,1))==1;
            

            for b = 1:numel(eds)-1,
% $$$                 for h = 1:numel(eds)-1,
% $$$                     grind = eds(b) <= mod(th(iRes,1),2*pi) & mod(th(iRes,1),2*pi) <= eds(b+1) & ...
% $$$                             eds(h) <= mod(thH(iRes,1),2*pi) & mod(thH(iRes,1),2*pi) <= eds(h+1);
% $$$                     grjnd = eds(b) <= mod(th(jRes,1),2*pi) & mod(th(jRes,1),2*pi) <= eds(b+1) & ...
% $$$                             eds(h) <= mod(thH(jRes,1),2*pi) & mod(thH(jRes,1),2*pi) <= eds(h+1);

                grind = eds(b) <= th(iRes,1) & th(iRes,1) <= eds(b+1);
                grjnd = eds(b) <= th(jRes,1) & th(jRes,1) <= eds(b+1);
% $$$             grind = ( abs(pfhxy(iRes,1))+100-pfhxy(iRes,2))>0&...
% $$$                     (-abs(pfhxy(iRes,1))-100-pfhxy(iRes,2))<0;
% $$$             grjnd = ( abs(pfhxy(jRes,1))+100-pfhxy(jRes,2))>0&...
% $$$                     (-abs(pfhxy(jRes,1))-100-pfhxy(jRes,2))<0;
                            
                if sum(grind)&sum(grjnd),
                    [tccg,tbin] = CCG([iRes(grind);jRes(grjnd)],...
                                      [ones([sum(grind),1]);2*ones([sum(grjnd),1])],...
                                      binSize,halfBins,spk.sampleRate,[1,2],normalization);
                else
                    tccg = zeros([halfBins*2+1,2,2]);
                end
% $$$                 mccg(:,b,h) = filt_fun(tccg(:,1,2));
% $$$                 sccg(:,b,h) = filt_fun(tccg(:,1,1));
% $$$                 occg(:,b,h) = filt_fun(tccg(:,2,2));

                mccg(:,b) = filt_fun(tccg(:,1,2));
                sccg(:,b) = filt_fun(tccg(:,1,1));
                occg(:,b) = filt_fun(tccg(:,2,2));
% $$$                 end
            end
% $$$             for h = 1:numel(eds)-1,
% $$$                 subplot2(8,3,h+4,1);
% $$$                 imagesc(tbin,eds-pi,sccg(:,:,h)');axis('tight');
% $$$                 subplot2(8,3,h+4,2);
% $$$                 %imagesc(tbin,eds-pi,bsxfun(@rdivide,mccg(:,:,h),max(mccg(:,:,h)))');axis('tight');
% $$$                 imagesc(tbin,eds-pi,mccg(:,:,h)');axis('tight');
% $$$                 %imagesc(tbin,eds,mccg');axis('tight');
% $$$                 Lines(0,[],'r');                
% $$$                 subplot2(8,3,h+4,3);
% $$$                 imagesc(tbin,eds-pi,occg(:,:,h)');axis('tight');
% $$$             end                
% $$$             
            subplot2(8,3,[5:8],1);
            imagesc(tbin,eds-pi-eds(2)/2,sccg');axis('tight');axis('xy');
            subplot2(8,3,[5:8],2);
            %imagesc(tbin,eds-pi-eds(2)/2,bsxfun(@rdivide,mccg,max(mccg))');axis('tight');axis('xy')
            bar(tbin,mccg(:,1,1));axis('tight');
            %imagesc(tbin,eds-pi-eds(2)/2,mccg');axis('tight');axis('xy');
            Lines(0,[],'m');                
            subplot2(8,3,[5:8],3);
            imagesc(tbin,eds-pi-eds(2)/2,occg');axis('tight');axis('xy');
            
            
% $$$             
% $$$             
% $$$ % $$$             grind = sign(dpfhxy(iRes,1))==-1;
% $$$ % $$$             grjnd = sign(dpfhxy(jRes,1))==-1;
% $$$             % negative direction            
% $$$             grind = eds(1) >= pcor{1}(iRes,1) | pcor{1}(iRes,1) >= eds(4);
% $$$             grjnd = eds(1) >= pcor{1}(jRes,1) | pcor{1}(jRes,1) >= eds(4);
% $$$             if sum(grind)&sum(grjnd),            
% $$$                 [tccg,tbin] = CCG([iRes(grind);jRes(grjnd)],...
% $$$                                   [ones([sum(grind),1]);2*ones([sum(grjnd),1])],...
% $$$                                   binSize,halfBins,spk.sampleRate,[1,2],normalization,...
% $$$                                   epochs);
% $$$ 
% $$$                 subplot2(8,3,6,1);
% $$$                 bar(tbin,tccg(:,1,1));axis('tight');
% $$$                 subplot2(8,3,6,2);                
% $$$                 bar(tbin,filt_fun(tccg(:,1,2)));axis('tight');
% $$$                 Lines(0,[],'r');
% $$$                 subplot2(8,3,6,3);
% $$$                 bar(tbin,tccg(:,2,2));axis('tight');
% $$$             end            
% $$$ 
% $$$             
% $$$             % orto positive
% $$$             grind = eds(3) <= pcor{1}(iRes,1) & pcor{1}(iRes,1) <= eds(4);
% $$$             grjnd = eds(3) <= pcor{1}(jRes,1) & pcor{1}(jRes,1) <= eds(4);
% $$$                         if sum(grind)&sum(grjnd),
% $$$                 [tccg,tbin] = CCG([iRes(grind);jRes(grjnd)],...
% $$$                                   [ones([sum(grind),1]);2*ones([sum(grjnd),1])],...
% $$$                                   binSize,halfBins,spk.sampleRate,[1,2],normalization);
% $$$                 subplot2(8,3,7,1);
% $$$                 bar(tbin,tccg(:,1,1));axis('tight');
% $$$                 subplot2(8,3,7,2);
% $$$                 bar(tbin,filt_fun(tccg(:,1,2)));axis('tight');
% $$$                 Lines(0,[],'r');                
% $$$                 subplot2(8,3,7,3);
% $$$                 bar(tbin,tccg(:,2,2));axis('tight');
% $$$             end
% $$$             
% $$$ % $$$             grind = sign(dpfhxy(iRes,1))==-1;
% $$$ % $$$             grjnd = sign(dpfhxy(jRes,1))==-1;
% $$$             
% $$$             % ortho negative
% $$$             grind = eds(1) <= pcor{1}(iRes,1) & pcor{1}(iRes,1) <= eds(2);
% $$$             grjnd = eds(1) <= pcor{1}(jRes,1) & pcor{1}(jRes,1) <= eds(2);
% $$$             if sum(grind)&sum(grjnd),            
% $$$                 [tccg,tbin] = CCG([iRes(grind);jRes(grjnd)],...
% $$$                                   [ones([sum(grind),1]);2*ones([sum(grjnd),1])],...
% $$$                                   binSize,halfBins,spk.sampleRate,[1,2],normalization,...
% $$$                                   epochs);
% $$$                 subplot2(8,3,8,1);
% $$$                 bar(tbin,tccg(:,1,1));axis('tight');
% $$$                 subplot2(8,3,8,2);
% $$$                 bar(tbin,filt_fun(tccg(:,1,2)));axis('tight');
% $$$                 Lines(0,[],'r');                
% $$$                 subplot2(8,3,8,3);
% $$$                 bar(tbin,tccg(:,2,2));axis('tight');
% $$$             end            

            colormap('jet');
            drawnow();
            waitforbuttonpress();
        
        end
    end
end




%% SEE MjgER2016F3V1_ccgXstate.m
% EXA State segmented ccgs



% LONG time scale

% GET Behavior scores ------------------------------------------------------------------------------
MjgER2016_load_bhv_erpPCA_scores();
% END GET Behavior scores ------------------------------------------------------------------------------


%
% GET bhv rate maps
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
%pmaps(isnan(pmaps)) = 0;



% $$$ sigUnits = any(abs(fsrcz(:,1:3))>=1.96,2);
% $$$ cc = FSrC(:,[2,1,3])+0.75;
% $$$ cc(~sigUnits,:) = repmat([0.75,0.75,0.75],[sum(~sigUnits),1]);
% $$$ 
% $$$ % GET mds map of behavior space
% $$$ D = pdist(FSrC(:,1:3));
% $$$ mapa = mdscale(D,2);
% $$$ 
% $$$ % $$$ scatter(-mapa(:,2),mapa(:,1),10,cc,'filled')
% $$$ mapa = [-mapa(:,2),mapa(:,1)];

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



pairCnt = 1;
unitCumCnt = 0;
epochs = [];
binSize = 10;
halfBins = 50;
normalization = 'hz';
eds = linspace(-pi,pi,3);
filt_fun = @(x) RectFilter(x,11,3);
mccg = zeros([2*halfBins+1,10000,2]);

for tind = 1:numel(Trials),
    disp(['MjgER2016F4V1:tind:',num2str(tind)]);
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
    [mrate,mpos] = pft{tind}.maxRate(unitSubset,true,'mean',0.9);

    nunits = numel(unitSubset);

     for i=1:nunits-1,
        for j = i+1:nunits,
            
            d = pdist2(mpos(i,:),mpos(j,:),'euclidean');
            if d<=400,
                pairUids  = cat(1,pairUids,[i,j]+unitCumCnt);
                pairClus = cat(1,pairClus,[tind,unitSubset([i,j])']);
                pairSDist = cat(1,pairSDist,d);
                pairEigs = cat(1,pairEigs,permute(eigScore{1}(pairUids(end,:),1:3),[3,1,2]));
                pairFSrC = cat(1,pairFSrC,permute(FSrC(pairUids(end,:),1:3),[3,1,2]));
                pairMDSs = cat(1,pairMDSs,permute(mapa(pairUids(end,:),:),[3,1,2]));
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


% F3A1 --------------------------------------------------------------------
% The analysis of unit mean xcorr as a function of physical and behavioral
% distance.
%
% Physical distance is measured in centimeters
% Behavioral distance is measured as the angle computed by unit behavior 
% F-scores computed in MTA:analysis:MjgER2016:MjgER2016F2V2:F2D1


% HELPER FUNCITONS --------------------------------------------------------
shuffle = @(x) x(randperm(numel(x)));
%--------------------------------------------------------------------------

% $$$ figure();
% $$$ plot3(reshape(pairFSrC(:,:,1),[],1),...
% $$$       reshape(pairFSrC(:,:,2),[],1),...
% $$$       reshape(pairFSrC(:,:,3),[],1),...
% $$$       '.');
% $$$ 
% $$$ shiftx = 1;
% $$$ shifty = 1;
% $$$ mdist = circ_dist(atan2(pairFSrC(:,1,1)+shiftx,pairFSrC(:,1,3)+shifty),...
% $$$                   atan2(pairFSrC(:,2,1)+shiftx,pairFSrC(:,2,3)+shifty));
% $$$ mdistLabel = 'FSr';
% $$$ shiftx = 0;
% $$$ shifty = 0;
% $$$ mdist = pairCC;
% $$$ mdistLabel = 'cc';
% $$$ 
% $$$ shiftx = 0;
% $$$ shifty = 0;
% $$$ mdist = sqrt(sum(([pairFSrC(:,1,1)+shiftx,pairFSrC(:,1,3)+shifty] ...
% $$$                  -[pairFSrC(:,2,1)+shiftx,pairFSrC(:,2,3)+shifty]).^2,2));
% $$$ mdistLabel = 'FSrd';
% $$$ 
% $$$ shiftx = 0;
% $$$ shifty = 0;
% $$$ mdist = circ_dist(atan2(pairMDSs(:,1,1)+shiftx,pairMDSs(:,1,2)+shifty),...
% $$$                   atan2(pairMDSs(:,2,1)+shiftx,pairMDSs(:,2,2)+shifty));
% $$$ mdistLabel = 'MDS';
% $$$ ind =  log10(pairMmAmp)>-3;
%ind = ~any(pairEigs(:,:,2)>1,2)&log10(pairMmAmp)>-3;
mdist = pairCC;

figure();
subplot(321);
hold('on');
plot(pairFSrC(ind,1,1)+shiftx,pairFSrC(ind,1,3)+shifty,'.');
plot(pairFSrC(ind,2,1)+shiftx,pairFSrC(ind,2,3)+shifty,'.');
% $$$ plot(pairMDSs(ind,1,1),pairMDSs(ind,1,2),'.');
% $$$ plot(pairMDSs(ind,2,1),pairMDSs(ind,2,2),'.');
% $$$ plot(pairEigs(ind,1,1)+shiftx,pairEigs(ind,1,3)+shifty,'.');
% $$$ plot(pairEigs(ind,2,1)+shiftx,pairEigs(ind,2,3)+shifty,'.');
grid('on');
subplot(322);
hist(mdist,100)
drawnow();


edx = [0,15,25,32.5,37,40];
nx = numel(edx);
ny = 11;
%edy = linspace(-pi,pi,ny);
edy = linspace(-1,1,ny);

indx = discretize(pairSDist/10,edx);
indy = discretize(mdist,edy);
ind = ~isnan(indx)&~isnan(indy)&ind;


subs = [indx(ind),indy(ind)];
vals = pairMmAmp(ind);
sz   = [numel(edx)-1,numel(edy)-1];

ecratem = accumarray(subs,vals,sz,@mean);
ecrates = accumarray(subs,vals,sz,@std);
ecratec = accumarray(subs,ones([sum(ind),1]),sz,@sum);

ecratemr = [];
ecratesr = [];
for i = 1:1000,
    vals = pairMmAmp(ind);
    for j = 1:nx-1,  vals(subs(:,1)==j) = shuffle(vals(subs(:,1)==j));  end
    ecratemr = cat(3,ecratemr,accumarray(subs,vals,sz,@mean));
    ecratesr = cat(3,ecratesr,accumarray(subs,vals,sz,@std));
end

ecratemdr = [];
ecratesdr = [];
for i = 1:1000,
    vals = pairMmAmp(ind);
    for j = 1:ny-1,  vals(subs(:,2)==j) = shuffle(vals(subs(:,2)==j));  end
    ecratemdr = cat(3,ecratemdr,accumarray(subs,vals,sz,@mean));
    ecratesdr = cat(3,ecratesdr,accumarray(subs,vals,sz,@std));
end

xbins = edx(1:end-1)+abs(diff(edx)/2);
ybins = edy(1:end-1)+abs(diff(edy)/2);


subplot(3,2,3);
imagescnan(1:numel(xbins),1:numel(ybins),ecratem');
colorbar();
caxis([0,1.5]);
set(gca,'XTick',[0:numel(xbins)]+0.5)
set(gca,'XTickLabels',edx);
xlabel('Distance (cm)');
set(gca,'YTick',[1:numel(ybins)])
set(gca,'YTickLabels',round(ybins,3));
ylabel([mdistLabel,' inter-placefield corrcoef'])
axis('xy');
title('Mean XCORR magnitude')

subplot(3,2,4);
imagesc(1:numel(xbins),1:numel(ybins),ecrates');
colorbar();
caxis([0,1.5]);
set(gca,'XTick',[0:numel(xbins)]+0.5)
set(gca,'XTickLabels',edx);
xlabel('Distance (cm)');
set(gca,'YTick',[1:numel(ybins)])
set(gca,'YTickLabels',round(ybins,3));
ylabel([mdistLabel,' inter-placefield corrcoef'])
axis('xy');
title('Std XCORR magnitude')

subplot(3,2,5);
%plot(1:numel(ybins),((ecratem-mean(ecratemr,3))./std(ecratemr,[],3))');grid('on');
imagesc(1:numel(xbins),1:numel(ybins),((ecratem-mean(ecratemr,3))./std(ecratemr,[],3))');
colorbar();
caxis([-4,4]);
set(gca,'XTick',[0:numel(xbins)]+0.5)
set(gca,'XTickLabels',edx);
xlabel('Distance (cm)');
set(gca,'YTick',[1:numel(ybins)])
set(gca,'YTickLabels',round(ybins,3));
ylabel([mdistLabel,' inter-placefield corrcoef'])
axis('xy');
title('zscore of shuffled Mean XCORR magnitude')

subplot(3,2,6);
%plot(1:numel(ybins),((ecrates-mean(ecratesr,3))./std(ecratesr,[],3))');grid('on')
imagesc(1:numel(xbins),1:numel(ybins),((ecrates-mean(ecratesr,3))./std(ecratesr,[],3))');
colorbar();
caxis([-4,4]);
set(gca,'XTick',[0:numel(xbins)]+0.5)
set(gca,'XTickLabels',edx);
set(gca,'YTick',[1:numel(ybins)])
set(gca,'YTickLabels',round(ybins,3));
axis('xy');
title('zscore of shuffled Std XCORR magnitude')


 



% COMPUTE R stats for circlar distances

indx = discretize(pairSDist/10,edx);
ind = ~isnan(indx);

r = [];
for x = 1:numel(edx),
    xind = indx==x;
    r(x) = sum(pairMmAmp(xind).*exp(i.*mdist(xind)),'omitnan')./sum(pairMmAmp(xind),'omitnan');
end

rs = [];
sMdist = nan([2000,numel(edx),1000]);
sMmAmp = nan([2000,numel(edx),1000]);
for x = 1:numel(edx),
    xind = indx==x;   
    for j = 1:1000,
        sMdist(1:sum(xind),x,j) = shuffle(mdist(xind));
        sMmAmp(1:sum(xind),x,j) = pairMmAmp(xind);
        rs(x,j) = sum(pairMmAmp(xind).*exp(i.*sMdist(nniz(sMdist(:,x,j)),x,j)),'omitnan')./...
                  sum(pairMmAmp(xind),'omitnan');
    end
end
sMdist = reshape(sMdist,[],1000);
sMmAmp = reshape(sMmAmp,[],1000);
sMdist(isnan(sMdist(:,1)),:) = [];
sMmAmp(isnan(sMmAmp(:,1)),:) = [];
rps = sum(sMmAmp.*exp(i.*sMdist),'omitnan')./sum(sMmAmp,'omitnan');
rp = sum(pairMmAmp.*exp(i.*mdist),'omitnan')./sum(pairMmAmp,'omitnan');
rpz = (abs(rp)-mean(abs(rps)))./std(abs(rps));

% RESULTS F3A1

%  zscore               rpz:15.2
%  zscore   edx: 0  -  15    rz: 7.09
%               15  -  25        8.74
%               25  -  32.5      7.11
%               32.5-  37        5.29
%               37  -  40        5.27



figure,hist(abs(rs(1,:)),100)
rz = (abs(r')-mean(abs(rs),2))./std(abs(rs),[],2);

figure,plot(edx,rz)

ecrater = accumarray(subs,vals,sz,);

% $$$ ecrates = accumarray(subs,vals,sz,@std);
% $$$ ecratec = accumarray(subs,ones([sum(ind),1]),sz,@sum);

ecratemr = [];
ecratesr = [];
for i = 1:1000,
    %subs = [shuffle(subs(:,1)),shuffle(subs(:,2))];
    vals = pairMmAmp(ind);
    for j = 1:nx-1
        vals(subs(:,1)==j) = shuffle(vals(subs(:,1)==j));
    end
    ecratemr = cat(3,ecratemr,accumarray(subs,vals,sz,@mean));
    ecratesr = cat(3,ecratesr,accumarray(subs,vals,sz,@std));
end


ecratemdr = [];
ecratesdr = [];
for i = 1:1000,
    %subs = [shuffle(subs(:,1)),shuffle(subs(:,2))];
    vals = pairMmAmp(ind);
    for j = 1:ny-1,
        vals(subs(:,2)==j) = shuffle(vals(subs(:,2)==j));
    end
    ecratemdr = cat(3,ecratemdr,accumarray(subs,vals,sz,@mean));
    ecratesdr = cat(3,ecratesdr,accumarray(subs,vals,sz,@std));
end

%pecrate = accumarray([indx,indy],pairPkAmp(ind),[10,10],@mean);
xbins = edx(1:end-1)+abs(diff(edx)/2);
ybins = edy(1:end-1)+abs(diff(edy)/2);



% COMPUTE R stats for circlar distances -------------------------------------------------------

edy = linspace(-1,1,7);
edx = linspace(0,40,7);

%indx = discretize(pairDCC,edy);
indx = discretize(pairSDist/10,edx);
indy = discretize(mdist,edy);

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
            %rssModel(x,y) = sum((pairMmAmp(xind)-condExpCcgsFSdistXBdist(x,y)).^2);
        end
    end
    xind = indx==x;
    %rssModelCol(x) = sum((pairMmAmp(xind)-condExpCcgsFSdistXBdist(sub2ind(size(condExpCcgsFSdistXBdist),repmat(x,[sum(xind),1]),indy(xind)))).^2,'omitnan');
end

% $$$ [ti,tj] = ind2sub(size(condExpCcgsFSdistXBdist),(indx(1)-1)*nxbins+indy(1))
% $$$ [ti]    = sub2ind(size(condExpCcgsFSdistXBdist),repmat(x,[sum(xind),1]),indy(xind));
% $$$ [ti]    = (indx(1)-1)*nxbins+indy(1);

figure();imagesc(edx,edy,condExpCcgsFSdistXBdist');axis('xy');



sMmAmp = nan([2000,nxbins,1000]);
sMdist = nan([2000,nxbins,1000]);
sSdist = nan([2000,nxbins,1000]);
for x = 1:nxbins,
    xind = indx==x;
    nXind = sum(xind);
    for j = 1:1000,
        sMmAmp(1:nXind,x,j) = shuffle(pairMmAmp(xind));
        sMdist(1:nXind,x,j) =             mdist(xind);
        sSdist(1:nXind,x,j) =         pairSDist(xind)/10;
    end
end


sMmAmp = reshape(sMmAmp,[],1000);
sMdist = reshape(sMdist,[],1000);
sSdist = reshape(sSdist,[],1000);

sMmAmp(isnan(sMmAmp(:,1)),:) = [];
sMdist(isnan(sMdist(:,1)),:) = [];
sSdist(isnan(sSdist(:,1)),:) = [];

% COMPUTE residual sum of squares for each shuffle along the rows of condExpCcgsFSdistXBdist 
% RSS = ∑(x-x̄)² 


rssShuffled = nan([nxbins,nybins,1000]);
rssShuffledCol = nan([nxbins,1000]);
for j = 1:1000,
    indx = discretize(sSdist(:,j),edx);
    indy = discretize(sMdist(:,j),edy);
    for x = 1:nxbins,
        for y = 1:nybins,
            xyind = indx==x&indy==y;
            rssShuffled(x,y,j) = sum((sMmAmp(xyind,j)-condExpCcgsFSdistXBdist(x,y)).^2);
        end
        xind = indx==x;
        rssShuffledCol(x,j) = sum((sMmAmp(xind,j)-condExpCcgsFSdistXBdist(sub2ind(size(condExpCcgsFSdistXBdist),repmat(x,[sum(xind),1]),indy(xind)))).^2,'omitnan');
    end
end

figure();
subplot(221);  imagesc(mean(rssModel,3)');       axis('xy');  colorbar();
subplot(223);  imagesc(mean(rssShuffled,3)');    axis('xy');  colorbar();
subplot(224);  imagesc(std(rssShuffled,[],3)');  axis('xy');  colorbar();


zsModelCol = (rssModelCol-mean(rssShuffledCol,2))./std(rssShuffledCol,[],2);

figure();
plot(mean(condExpCcgsFSdistXBdist,2));



X = [pairSDist/100,mdist];
bhat = inv(X'*X)*X'*pairMmAmp;
yhat = X*bhat;

rerr = pairMmAmp - yhat;

ind = pairMmAmp~=0;

[beta,Sigma,E,CovB,logL] = mvregress([pairSDist(ind)/10,mdist(ind)],pairMmAmp(ind));
[beta,Sigma,E,CovB,logL] = mvregress([pairSDist(ind)/10,mdist(ind)],log10(pairMmAmp(ind)+1));
[beta,Sigma,E,CovB,logL] = mvregress([pairSDist(ind)/10,mdist(ind)],log10(max(pairPkAmp(ind,:),[],2)+1));
[beta,Sigma,E,CovB,logL] = mvregress([pairDCC(ind),mdist(ind)],log10(pairMmAmp(ind)+1));

[beta,Sigma,E,CovB,logL] = mvregress([pairSDist/10,mdist],pairMmAmp);

%mse = mean(sqrt(E.^2));

figure,plot3(pairDCC(ind),mdist(ind),log10(pairMmAmp(ind)+1),'.');
figure,plot3(pairDCC(ind),mdist(ind),log10(max(pairPkAmp(ind,:),[],2)+1),'.');
figure,plot3(pairSDist(ind)./10,mdist(ind),log10(max(pairPkAmp(ind,:),[],2)+1),'.');
xlabel('DistCC');ylabel('BhvCC');zlabel('log10(CCG+1)');
figure,plot3(pairSDist/10,mdist,pairMmAmp,'.');
figure,plot3(pairSDist/10,mdist,pairMmAmp,'.');

for j = 1:1000,
    [betaShuffled(:,j),~,~,~,lugLShuffled(j)] = mvregress([sSdist(:,j)/10,sMdist(:,j)],sMmAmp(:,j));
end

(beta-mean(betaShuffled,2))./std(betaShuffled,[],2)

rps = sum(sMmAmp.*exp(i.*sMdist),'omitnan')./sum(sMmAmp,'omitnan');
rp = sum(pairMmAmp.*exp(i.*mdist),'omitnan')./sum(pairMmAmp,'omitnan');
rpz = (abs(rp)-mean(abs(rps)))./std(abs(rps));


% RESULTS F3A1


figure,hist(abs(rs(1,:)),100)
rz = (abs(r')-mean(abs(rs),2))./std(abs(rs),[],2);

figure,plot(edx,rz)

ecrater = accumarray(subs,vals,sz,);

% $$$ ecrates = accumarray(subs,vals,sz,@std);
% $$$ ecratec = accumarray(subs,ones([sum(ind),1]),sz,@sum);

ecratemr = [];
ecratesr = [];
for i = 1:1000,
    %subs = [shuffle(subs(:,1)),shuffle(subs(:,2))];
    vals = pairMmAmp(ind);
    for j = 1:nx-1
        vals(subs(:,1)==j) = shuffle(vals(subs(:,1)==j));
    end
    ecratemr = cat(3,ecratemr,accumarray(subs,vals,sz,@mean));
    ecratesr = cat(3,ecratesr,accumarray(subs,vals,sz,@std));
end


ecratemdr = [];
ecratesdr = [];
for i = 1:1000,
    %subs = [shuffle(subs(:,1)),shuffle(subs(:,2))];
    vals = pairMmAmp(ind);
    for j = 1:ny-1,
        vals(subs(:,2)==j) = shuffle(vals(subs(:,2)==j));
    end
    ecratemdr = cat(3,ecratemdr,accumarray(subs,vals,sz,@mean));
    ecratesdr = cat(3,ecratesdr,accumarray(subs,vals,sz,@std));
end

%pecrate = accumarray([indx,indy],pairPkAmp(ind),[10,10],@mean);
xbins = edx(1:end-1)+abs(diff(edx)/2);
ybins = edy(1:end-1)+abs(diff(edy)/2);
