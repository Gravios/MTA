% LOAD Data
MjgER2016_load_data
%  Variables:
%      Trials
%      units
%      cluSessionMap
%      pitchReferenceTrial
%      FigDir
%      sessionListName
%      sessionList
%      states
%      numStatesg
%      interpParPfsp
%      interpParDfs
%
%  Functions:
%      reshape_eigen_vector



% EXA State segmented ccgs

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



figure();
eds = linspace(0,2*pi,2);
filt_fun = @(x) RectFilter(x,3,3);

binSize = 1;
halfBins = 12;
normalization = 'hz';


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

            mccg(:,end+1,1,1) = zeros([2*halfBins+1,1,1,1]);
            sccg(:,end+1,1,1) = zeros([2*halfBins+1,1,1,1]);
            occg(:,end+1,1,1) = zeros([2*halfBins+1,1,1,1]);

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
