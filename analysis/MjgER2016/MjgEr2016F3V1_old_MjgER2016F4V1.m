


MjgER2016_load_data



% Pairwise unit ccg 
% 1. find placefields' centers
% 2. compute distance matrix of placefield centers
% 3. find midpoint between placefields which have centers less than 20 cm apart
% 4. construct basis from midpoint to one field center
% 5. construct rotation matrix to transform trajectories from room to field pair frame of reference
% 6. compute derivative of primary axis from midpoint to field
% 7. convert transformed trajectories to polar coordinates
% 8. compute ccg asymmetry for angle-distance pairs or just time-distance 
%    separeted by travel direction
%

pft = cf(@(T,u)  pfs_2d_theta(T,u,'overwrite',false),  Trials,units);
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






binSize = 1;
halfBins = 12;
normalization = 'hz';

pfs = cf(@(t,u)  pfs_2d_states(t,u),  Trials,units);


%filt_fun = @(x) x;
mccg = zeros([2*halfBins+1,0,1,1]);
sccg = zeros([2*halfBins+1,0,1,1]);
occg = zeros([2*halfBins+1,0,1,1]);
unitPairs = [];

for tind = 1:numel(Trials),
nunits = numel(units{tind});

for i=1:nunits-1,
    for j = i+1:nunits,
        d = pdist2(mpos(i,:),mpos(j,:),'euclidean');
        if d<=200,
            unitPairs = cat(1,unitPairs,[tind,units{tind}([i,j])]);
            midpoint = sum(mpos([i,j],:))./2;
            
% $$$             clf();
% $$$ % PLOT Place fields
% $$$             subplot(431); hold('on');
% $$$             plot(pft{tind},units{tind}(i),1,true,[],true);
% $$$             plot(mpos(i,1),mpos(i,2),'*m');
% $$$             plot(midpoint(1),midpoint(2),'*g');
% $$$             subplot(433); hold('on');
% $$$             plot(pft{tind},units{tind}(j),1,true,[],true);
% $$$             plot(mpos(j,1),mpos(j,2),'*m'); 
% $$$             plot(midpoint(1),midpoint(2),'*g');

% PLOT Behavior fields
% $$$             subplot(434); plot(pfd{1},units{tind}(i),1,true,[],false);            
% $$$             subplot(436); plot(pfd{1},units{tind}(j),1,true,[],false);
% CONSTRUCT basis based on the midpoint and second placefield center
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

            mccg(:,end+1,b,s) = zeros([2*halfBins+1,1,1,1]);
            sccg(:,end+1,b,s) = zeros([2*halfBins+1,1,1,1]);
            occg(:,end+1,b,s) = zeros([2*halfBins+1,1,1,1]);

            for s = 1:numel(states),
                iRes = spk(units{tind}(i));
                jRes = spk(units{tind}(j));
                iRes = SelectPeriods(iRes,Trials{tind}.stc{states{s},xyz.sampleRate},'d',1,0);
                jRes = SelectPeriods(jRes,Trials{tind}.stc{states{s},xyz.sampleRate},'d',1,0);
                for b = 1:numel(eds)-1,
                    grind = eds(b) <= th(iRes,1) & th(iRes,1) <= eds(b+1);
                    grjnd = eds(b) <= th(jRes,1) & th(jRes,1) <= eds(b+1);
                    if sum(grind)&sum(grjnd),
                        [tccg,tbin] = CCG([iRes(grind);jRes(grjnd)],...
                                          [ones([sum(grind),1]);2*ones([sum(grjnd),1])],...
                                          binSize,halfBins,spk.sampleRate,[1,2],normalization);
                    else
                        tccg = zeros([halfBins*2+1,2,2]);
                    end
                    mccg(:,end,b,s) = filt_fun(tccg(:,1,2));
                    sccg(:,end,b,s) = filt_fun(tccg(:,1,1));
                    occg(:,end,b,s) = filt_fun(tccg(:,2,2));
                end
            end

% $$$             subplot2(8,3,[5:8],1);
% $$$             imagesc(tbin,eds-pi-eds(2)/2,sccg');axis('tight');axis('xy');
% $$$             subplot2(8,3,[5:8],2);
% $$$             imagesc(tbin,eds-pi-eds(2)/2,bsxfun(@rdivide,mccg,max(mccg))');axis('tight');axis('xy')
% $$$ 
% $$$             Lines(0,[],'m');                
% $$$             subplot2(8,3,[5:8],3);
% $$$             imagesc(tbin,eds-pi-eds(2)/2,occg');axis('tight');axis('xy');
% $$$             
% $$$             colormap('jet');
% $$$             drawnow();
% $$$             waitforbuttonpress();
        end
    end
end
end



figure,
for u = 1:size(unitPairs,1),
    for s = 1:numel(states),
        subplot(numel(states),1,s);
        imagesc(tbin,eds-pi-eds(2)/2,sq(mccg(:,u,:,s))');axis('tight');axis('xy')
        title([states{s},': ',num2str(unitPairs(u,1)),' vs ',num2str(unitPairs(u,2))]);
    end
    waitforbuttonpress();
end




% LONG time scale


% GET Behavior scores

pfindex = 1;

pfdShuffled =  cf(@(t,g) MTAApfs(t,'tag',g), Trials, repmat({'HBPITCHxBPITCH_shuffled'},size(Trials)));

% GET bhv rate maps
rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfd(:,pfindex),units');
clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfd(:,pfindex),units');    
tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
rmaps = cat(2, rmaps{:});   
clu = cat(2, clu{:});    
tlu = cat(2, tlu{:});    
clu = [tlu',clu'];
[~,rind] = sortrows(clu);
rmaps = rmaps(:,rind);
rmaps = rmaps(validDims{pfindex},unitSubsets{pfindex});
rmaps(isnan(rmaps)) = 0;

rmapsShuffledMean = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfdShuffled',units');
clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfdShuffled',units');    
tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
rmapsShuffledMean = cat(2, rmapsShuffledMean{:});   
clu = cat(2, clu{:});    
tlu = cat(2, tlu{:});    
clu = [tlu',clu'];
[~,rind] = sortrows(clu);
rmapsShuffledMean = rmapsShuffledMean(:,rind);
rmapsShuffledMean = rmapsShuffledMean(validDims{pfindex},unitSubsets{pfindex});
rmapsShuffledMean(isnan(rmapsShuffledMean)) = 0;

rmapsShuffled = cf(@(p,u) p.data.rateMap(:,ismember(p.data.clu,u),:), pfdShuffled',units');
clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfdShuffled',units');    
tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
rmapsShuffled = cat(2, rmapsShuffled{:});   
clu = cat(2, clu{:});    
tlu = cat(2, tlu{:});    
clu = [tlu',clu'];
[~,rind] = sortrows(clu);
rmapsShuffled = rmapsShuffled(:,rind,:);
rmapsShuffled = rmapsShuffled(validDims{pfindex},unitSubsets{pfindex},:);
rmapsShuffled(isnan(rmapsShuffled)) = 0;

D = cov(rmapsShuffledMean');
LR = eigVec{pfindex};
% compute rotated FS coefficients        
FSCFr = LR * inv(LR' * LR);          % this is pseudo-inverse of LR
% rescale rotated FS coefficients by the corresponding SDs 
rk = size(FSCFr,2);%rank(D,1e-4);       % why not on D? would save on corrcoef(X) computation
FSCFr = FSCFr .* repmat(sqrt(diag(D)),1,rk);
% compute rotated factor scores from the normalized raw data and  the
% corresponding rescaled factor score coefficients
rsMean = mean(rmapsShuffledMean');
rsStd  = std( rmapsShuffledMean');

% MEAN shuffled score
FSrM =  ((rmapsShuffledMean'-rsMean)./rsStd) * FSCFr;

% MEAN normal score
FSrC =  ((rmaps'-rsMean)./rsStd) * FSCFr;

FSrS = [];
for i = 1:pfd{1}.parameters.numIter
    FSrS(:,:,end+1) =  ((rmapsShuffled(:,:,i)'-rsMean)./rsStd) * FSCFr;
end
fsrsMean = mean(FSrS,3);
fsrsStd = std(FSrS,[],3);

fsrcz = (FSrC-fsrsMean)./fsrsStd;


sigUnits = any(abs(fsrcz(:,1:3))>=1.96,2);
cc = FSrC(:,[2,1,3])+0.75;
cc(~sigUnits,:) = repmat([0.75,0.75,0.75],[sum(~sigUnits),1]);


% GET mds map of behavior space
D = pdist(FSrC(:,1:3));
mapa = mdscale(D,2);

% $$$ figure();
% $$$ scatter(-mapa(:,2),mapa(:,1),10,cc,'filled')
mapa = [-mapa(:,2),mapa(:,1)];



cluSessionMapSubset = cluSessionMap(unitSubsets{1},:);
pairClus   = [];
pairSDist  = [];
pairBDist  = [];
pairMmAmp  = [];
pairPkAmp  = [];
pairPkTime = [];
pairEigs   = [];
pairFSrC   = [];
pairMDSs   = [];


for tind = 1:23,
    disp(['MjgER2016F4V1:tind:',num2str(tind)]);
xyz = preproc_xyz(Trials{tind},'trb');
xyz.filter('RectFilter');
xyz.filter('ButFilter',5,2.5,'low');
vxy = xyz.vel('spine_lower');
feature = xyz.copy();
feature.data = sq(feature(:,'nose',[1,2]));
spk = Trials{tind}.spk.copy();
spk = create(spk,Trials{tind},xyz.sampleRate,'theta-groom-sit',units{tind},'deburst'); 

unitSubset = cluSessionMapSubset(cluSessionMapSubset(:,1)==tind,2);
[mrate,mpos] = pft{tind}.maxRate(unitSubset,true,'mean',0.9);

epochs = [];
binSize = 5;
halfBins = 100;
normalization = 'hz';

% figure;
eds = linspace(0,2*pi,3);
filt_fun = @(x) RectFilter(x,11,3);
mccg = [];sccg = [];occg = [];

nunits = numel(unitSubset);
for i=1:nunits-1,
    for j = i+1:nunits,
        
        d = pdist2(mpos(i,:),mpos(j,:),'euclidean');
        if d<=400,
            pairClus = cat(1,pairClus,[tind,unitSubset([i,j])']);
            pairSDist = cat(1,pairSDist,d);
            
            pairEigs = cat(1,pairEigs,permute(eigScore{1}(ismember(cluSessionMapSubset, ...
                                                              [[tind;tind],unitSubset([i,j])],'rows'),1:3),[3,1,2]));
            pairFSrC = cat(1,pairFSrC,permute(FSrC(ismember(cluSessionMapSubset,[[tind;tind],unitSubset([i,j])],'rows'),1:3),[3,1,2]));
            pairMDSs = cat(1,pairMDSs,permute(mapa(ismember(cluSessionMapSubset,[[tind;tind],unitSubset([i,j])],'rows'),:),[3,1,2]));
            pairBDist = cat(1,pairBDist,...
                            sqrt(sum(diff(sq(pairEigs(end,:,:))...
                                          ).^2)));
            midpoint = sum(mpos([i,j],:))./2;
            
% $$$             clf();
% $$$             subplot(431);
% $$$             plot(pft{tind},unitSubset(i),1,true,[],true);
% $$$             hold('on');
% $$$             plot(mpos(i,1),mpos(i,2),'*m');
% $$$             plot(midpoint(1),midpoint(2),'*g');            
% $$$             subplot(433);
% $$$             plot(pft{tind},unitSubset(j),1,true,[],true);
% $$$             hold('on');            
% $$$             plot(mpos(j,1),mpos(j,2),'*m');
% $$$             plot(midpoint(1),midpoint(2),'*g');
% $$$ 
% $$$             subplot(434);
% $$$             plot(pfd{tind,1},unitSubset(i),1,true,[],false);            
% $$$             subplot(436);
% $$$             plot(pfd{tind,1},unitSubset(j),1,true,[],false);

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
            
            ii = unitSubset(i)==unitSubset;
            jj = unitSubset(j)==unitSubset;

            iRes = spk(unitSubset(i));
            jRes = spk(unitSubset(j));

            for b = 1:numel(eds)-1,
                grind = eds(b) <= th(iRes,1) & th(iRes,1) <= eds(b+1);
                grjnd = eds(b) <= th(jRes,1) & th(jRes,1) <= eds(b+1);
                            
                if sum(grind)&sum(grjnd),
                    [tccg,tbin] = CCG([iRes(grind);jRes(grjnd)],...
                                      [ones([sum(grind),1]);2*ones([sum(grjnd),1])],...
                                      binSize,halfBins,spk.sampleRate,[1,2],normalization);
                else
                    tccg = zeros([halfBins*2+1,2,2]);
                end
                mccg(:,b) = filt_fun(tccg(:,1,2));
                sccg(:,b) = filt_fun(tccg(:,1,1));
                occg(:,b) = filt_fun(tccg(:,2,2));
            end
            
            [pkAmp,pkInd] = max(mean([mccg(:,1),flipud(mccg(:,2))],2));
            mmAmp = max(mean([mccg(:,1),flipud(mccg(:,2))]));
            pairPkAmp = cat(1,pairPkAmp,pkAmp);
            pairMmAmp = cat(1,pairMmAmp,mmAmp);
            pairPkTime = cat(1,pairPkTime,tbin(pkInd));
% $$$ 
% $$$             subplot2(8,3,[5:8],1);
% $$$             imagesc(tbin,eds-pi-eds(2)/2,sccg');axis('tight');axis('xy');
% $$$             subplot2(8,3,[5:8],2);
% $$$             hold('on');
% $$$             %imagesc(tbin,eds-pi-eds(2)/2,bsxfun(@rdivide,mccg,max(mccg))');axis('tight');axis('xy')
% $$$             %bar(tbin,mccg(:,1,1));axis('tight');            
% $$$             stairs(tbin,mccg(:,1));axis('tight');
% $$$             stairs(tbin,mccg(:,2));axis('tight');             
% $$$             %imagesc(tbin,eds-pi-eds(2)/2,mccg');axis('tight');axis('xy');
% $$$             Lines(0,[],'m');                
% $$$             subplot2(8,3,[5:8],3);
% $$$             imagesc(tbin,eds-pi-eds(2)/2,occg');axis('tight');axis('xy');
% $$$ 
% $$$             colormap('jet');
% $$$             drawnow();
% $$$             waitforbuttonpress();

        end
    end
end

end
  






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



shiftx = 1;
shifty = 1;
mdist = circ_dist(atan2(pairFSrC(:,1,1)+shiftx,pairFSrC(:,1,3)+shifty),...
                  atan2(pairFSrC(:,2,1)+shiftx,pairFSrC(:,2,3)+shifty));
mdistLabel = 'FSr';

shiftx = 0;
shifty = 0;
mdist = circ_dist(atan2(pairMDSs(:,1,1)+shiftx,pairMDSs(:,1,2)+shifty),...
                  atan2(pairMDSs(:,2,1)+shiftx,pairMDSs(:,2,2)+shifty));
mdistLabel = 'MDS';

% $$$ shiftx = 1.1;
% $$$ shifty = 0.6;
% $$$ mdist = circ_dist(atan2(pairEigs(:,1,1)+shiftx,pairEigs(:,1,3)+shifty),...
% $$$                   atan2(pairEigs(:,2,1)+shiftx,pairEigs(:,2,3)+shifty));


ind =  log10(pairMmAmp)>-3;
% $$$ ind = ~any(pairEigs(:,:,2)>0&pairEigs(:,:,1)<-0.5,2);
% $$$ ind = all([ismember(pairClus(:,[1,2]) ,cluSessionMapSubset(sigUnits,:)),...
% $$$            ismember(pairClus(:,[1,3]) ,cluSessionMapSubset(sigUnits,:))],2);
% $$$ ind = any(pairFSrC(:,:,2)<0,2);

figure();
subplot(321);
hold('on');
% $$$ plot(pairFSrC(ind,1,1)+shiftx,pairFSrC(ind,1,2)+shifty,'.');
% $$$ plot(pairFSrC(ind,2,1)+shiftx,pairFSrC(ind,2,2)+shifty,'.');
plot(pairMDSs(ind,1,1),pairMDSs(ind,1,2),'.');
plot(pairMDSs(ind,2,1),pairMDSs(ind,2,2),'.');
% $$$ plot(pairEigs(ind,1,1)+shiftx,pairEigs(ind,1,3)+shifty,'.');
% $$$ plot(pairEigs(ind,2,1)+shiftx,pairEigs(ind,2,3)+shifty,'.');
grid('on');
subplot(322);
rose(mdist(ind)),
drawnow();


%ind = pairPkAmp>0.5;
%ind = log10(pairMmAmp)>-3&all(pairEigs(:,:,2)<1,2);
%ind = all(pairEigs(:,:,2)<0.5,2);
%ind = all(pairEigs(:,:,1)>0,2)&all(pairEigs(:,:,3)>0,2);
%edx = linspace(0,40,5);

edx = [0,15,25,32.5,37,40];
%edx = [0,20,27.5,32.5,35,38,40];
nx = numel(edx);
% $$$ nx = 2;
% $$$ edx = [0,20];

%edy = linspace(-pi,pi,21);

ny = 16
edy = linspace(-pi,pi,ny);
% $$$ XMin = -pi;
% $$$ XMax =  pi;
% $$$ mycdf = @(x) cdf('norm',x,0,1);
% $$$ probspace = @(CDF, XMin, XMax, N) arrayfun( @(p) max(fsolve( @(x) CDF(x)-p,...
% $$$                                                   [-0.01, 0.01],...
% $$$                                                   optimoptions('fsolve','Algorithm','trust-region'))), linspace(0,1,N) );
% $$$ edy = probspace(mycdf,0.02,0.98,ny);


indx = discretize(pairSDist/10,edx);
%indy = discretize(pairBDist(ind),edy);
indy = discretize(mdist,edy);

ind = ~isnan(indx)&~isnan(indy)&ind;

%A = ACCUMARRAY(SUBS,VAL,SZ,FUN)

%vals = pairMmAmp(ind);
subs = [indx(ind),indy(ind)];
vals = pairMmAmp(ind);
%vals = pairPkAmp(ind);
sz   = [numel(edx)-1,numel(edy)-1];
ecratem = accumarray(subs,vals,sz,@mean);
ecrates = accumarray(subs,vals,sz,@std);
ecratec = accumarray(subs,ones([sum(ind),1]),sz,@sum);

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

 
subplot(3,2,3);
imagesc(1:numel(xbins),1:numel(ybins),ecratem');
colorbar();
caxis([0,1]);
set(gca,'XTick',[0:numel(xbins)]+0.5)
set(gca,'XTickLabels',edx);
xlabel('Distance (cm)');
set(gca,'YTick',[1:numel(ybins)])
set(gca,'YTickLabels',round(ybins,3));
ylabel([mdistLabel,' inter-placefield angle (rad)'])
axis('xy');
title('Mean XCORR magnitude')

subplot(3,2,4);
imagesc(1:numel(xbins),1:numel(ybins),ecrates');
colorbar();
caxis([0,1]);
set(gca,'XTick',[0:numel(xbins)]+0.5)
set(gca,'XTickLabels',edx);
xlabel('Distance (cm)');
set(gca,'YTick',[1:numel(ybins)])
set(gca,'YTickLabels',round(ybins,3));
ylabel([mdistLabel,' inter-placefield angle (rad)'])
axis('xy');


subplot(3,2,5);
%plot(1:numel(ybins),((ecratem-mean(ecratemr,3))./std(ecratemr,[],3))');grid('on');
imagesc(1:numel(xbins),1:numel(ybins),((ecratem-mean(ecratemr,3))./std(ecratemr,[],3))');
colorbar();
caxis([-4,4]);
set(gca,'XTick',[0:numel(xbins)]+0.5)
set(gca,'XTickLabels',edx);
set(gca,'YTick',[1:numel(ybins)])
set(gca,'YTickLabels',round(ybins,3));
axis('xy');


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


subplot(3,2,5);
imagesc(1:numel(xbins),1:numel(ybins),((ecratem-mean(ecratemr,3))./std(ecratemr,[],3))');
colorbar();
caxis([-4,4]);
set(gca,'XTick',[0:numel(xbins)]+0.5)
set(gca,'XTickLabels',edx);
set(gca,'YTick',[1:numel(ybins)])
set(gca,'YTickLabels',round(ybins,3));
axis('xy');

subplot(3,2,6);
imagesc(1:numel(xbins),1:numel(ybins),((ecrates-mean(ecratesr,3))./std(ecratesr,[],3))');
colorbar();
caxis([-4,4]);
set(gca,'XTick',[0:numel(xbins)]+0.5)
set(gca,'XTickLabels',edx);
set(gca,'YTick',[1:numel(ybins)])
set(gca,'YTickLabels',round(ybins,3));
axis('xy');






figure();
hold('on');
plot(ybins,ecratem(1,:));
plot(ybins,ecratemr(1,:));
grid('on');
daspect([1,1,1])





figure,hist2([shuffle(pairPkAmp(ind)),shuffle(mdist(ind))],50,50)



figure,hold('on');
plot(ecratem(:,:,1),'b')
plot(sq(ecratemr(:,:,1:20)),'r')

