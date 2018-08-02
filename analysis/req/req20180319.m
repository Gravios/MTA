% req20180319 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Author: Justin Graboski
%  Final_Forms: NA
%  Description: Compute z-score of randomly shuffled spike time chunks 
%               projected on each eigen vector, as computed by 
%               req20180123_ver5.m
%  Bugs: NA



sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
Trials = af(@(s) MTATrial.validate(s), sessionList);
units = cf(@(T)  select_placefields(T),  Trials); 
units = req20180123_remove_bad_units(units);

cluMap = [];
for u = 1:numel(units)
    cluMap = cat(1,cluMap,[u*ones([numel(units{u}),1]),units{u}(:)]);
end




pft = cf(@(T,u)  pfs_2d_theta(T,u,'overwrite',false),  Trials,units);

[pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(Trials,[],[],false,true);

pfdShuffled = {};


%% COMPUTE/LOAD pfd shuffled %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tind = 19;
pfindex = 1;
tags={}; fetSets={}; fetInds={}; pfdParam={}; states={}; ranges={};
tags{end+1}        = ['HBPITCHxBPITCH_shuffled'];
fetSets{end+1}     = 'fet_HB_pitchB';
fetInds{end+1}     = [1,2];
boundaryLimits     = [-2,2;-2,2];
binDims            = [0.1,0.1];
SmoothingWeights   = [1.5,1.5];
states{end+1}      = 'theta-groom-sit';
ranges{end+1}      = [1100,1450];

for tind = 1:numel(Trials),
    Trial = Trials{tind};
    
    xyz = preproc_xyz(Trial,'trb');
    drz = compute_drz(Trial,units{tind},pft{tind});
    tper =[Trial.stc{states{pfindex}}];        
    tper.resample(xyz);        
    fet = feval(fetSets{pfindex},Trial);
    fet.data = fet(:,fetInds{pfindex});

    pargs = get_default_args('MjgER2016','MTAApfs','struct');        
    pargs.tag              = tags{1};
    pargs.units            = units{tind};
    pargs.numIter          = 1001;
    pargs.posShuffle       = true;
    pargs.halfsample       = false;
    pargs.overwrite        = true;
    pargs.boundaryLimits   = boundaryLimits;
    pargs.binDims          = binDims;
    pargs.SmoothingWeights = SmoothingWeights;
    pargs.xyzp             = fet.copy();
    pargs.bootstrap        = false;
    pargs.autoSaveFlag     = false;

    drzState = {};
    u = 1;        
    dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                     fet.sampleRate,fet.sync.copy(),fet.origin,'TimePeriods','sts',[],'tdrz','d');
    drzState{u} = dper&tper;
    pargs.units  = units{tind}(u);
    pargs.states = drzState{u};
    pfsArgs = struct2varargin(pargs);
    pfTemp = MTAApfs(Trial,pfsArgs{:});
    pfTemp.save();
    for u = 2:numel(units{tind});
        dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                         fet.sampleRate,fet.sync.copy(),fet.origin,'TimePeriods','sts',[],'tdrz','d');
        drzState{u} = dper&tper;
        pargs.units  = units{tind}(u);
        pargs.states = drzState{u};
        pfsArgs = struct2varargin(pargs);
        try
            pfTemp = MTAApfs(pfTemp,pfsArgs{:});            
        catch err,
            disp(err);
            pfTemp = MTAApfs(pfTemp,pfsArgs{:});
        end
    end
    pfTemp.save();        
end

pfdShuffled =  cf(@(t,g) MTAApfs(t,'tag',g), Trials, repmat(tags(pfindex),size(Trials)));
pfdShuffled =  cf(@(t,g) MTAApfs(t,'tag',g), Trials, repmat({'HBPITCHxBPITCH_shuffled'},size(Trials)));


%% COMPUTE ZSCORE OF erpPCA FSTAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SET helper function
reshape_eigen_vector = @(V,pfd) reshape(V(:),pfd{1}.adata.binSizes')';

numComp = 5;
pfindex = 1;

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

figure,plot(fsrcz(:,2),fsrcz(:,1),'.');hold('on');circle(0,0,2.5,'r');
figure,plot(fsrcz(:,2),fsrcz(:,3),'.');hold('on');circle(0,0,2.5,'r');

mapz = tsne([fsrcz(:,1:3),si(unitSubsets{pfindex})'],[],2,4,30);
figure,plot(mapz(:,1),mapz(:,2),'.');




% RECONSTRUCT eigenvectors into bhv space
numComp = 5;
pfindex = 1;
fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan(size(validDims{pfindex}));
    fpc{i}(validDims{pfindex}) = eigVec{pfindex}(:,i)./sqrt(nansum(eigVec{pfindex}(:,1).^2));
end
% PLOT eigenvectors
figure
for i = 1:numComp,
    subplot(1,numComp,i); imagesc(reshape(fpc{i}(:),pfd{1,pfindex}.adata.binSizes')');axis('xy');
end








% $$$ validDims{1} = vDims;
% $$$ eigScore{1} = FSr;
% $$$ eigVec{1} = LR;
% $$$ unitSubsets{1} = unitSubset;
% $$$ pfindex = 1;
% $$$ 
% $$$ validDims{2} = vDims;
% $$$ eigScore{2} = FSr;
% $$$ eigVec{2} = LR;
% $$$ unitSubsets{2} = unitSubset;
% $$$ pfindex = 2;

tind = 20;
numComp = 5;

fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan(size(validDims{pfindex}));
    fpc{i}(validDims{pfindex}) = eigVec{pfindex}(:,i)./sqrt(nansum(eigVec{pfindex}(:,1).^2));
    %fpc{i}(validDims{pfindex}) = fpc{i}(validDims{pfindex})-mean(fpc{i}(validDims{pfindex}));
end
figure
for i = 1:numComp,
    subplot(1,numComp,i); imagesc(reshape(fpc{i}(:),pfd{1,pfindex}.adata.binSizes')');axis('xy');
end

% GET auxiliary features
clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfd(:,pfindex),units');    
si  =  cf(@(p,u) p.data.si(:,ismember(p.data.clu,u),:), pfd(:,pfindex),units');
tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
si  = cat(2, si{:});   
clu = cat(2, clu{:});    
tlu = cat(2, tlu{:});    
clu = [tlu',clu'];


[~,rind] = sortrows(clu);
si  = si(:,rind);

%cf(@(T) T.load('nq'), Trials);
%anq = cellfun(@(t,u), t.nq.

% EIGSCORE 

mazeMask = fpc{i};
mazeMask(~isnan(mazeMask)) = 1;
mazeMask = reshape(mazeMask(:),pfd{1}.adata.binSizes');

map = tsne(eigScore{pfindex}(:,1:5),[],2,5,18);

mapa = tsne([FSrC(:,1:3),si(unitSubsets{pfindex})'],[],2,3,25);
%mapa = tsne([FSrC(:,1:3)],[],2,3,25);
figure();scatter(mapa(:,1)/10, mapa(:,2)/10,10,cc,'o','filled'); grid('on');
figure();scatter(mapa(:,1)/10+randn([size(mapa,1),1])/10, mapa(:,2)/10+randn([size(mapa,1),1])/10,10,cc,'o','filled'); grid('on');

mapa = tsne([FSrC(:,1:3)],[],2,3,30);
mapa = tsne([fsrcz(:,1:3)],[],2,4,18);
mapa = tsne([eigScore{pfindex}(:,1:3),si(unitSubsets{pfindex})'],[],2,3,25);
figure,plot(mapa(:,1),mapa(:,2),'.');

mapz = tsne([fsrcz(:,1:3),si(unitSubsets{pfindex})',eigScore{pfindex}(:,1:3)],[],2,7,25);
figure,plot(mapz(:,1),mapz(:,2),'.');

%pfi = cf(@(us,ui) ismember(us,ui), unitSubsets,repmat({unitIntersection},size(unitSubsets)));
%mapa = tsne([eigScore{1}(pfi{1},1:5),eigScore{2}(pfi{2},1:2)],[],2,7,20);

hfig = figure(3203920);
clf();
sp = gobjects([1,0]);
% EIGENVECTORS 
subplot(351);       imagesc(pfd{1}.adata.bins{:},reshape(fpc{1}(:),pfd{1,pfindex}.adata.binSizes')');axis('xy'); title('rear');
subplot(352);       imagesc(pfd{1}.adata.bins{:},reshape(fpc{2}(:),pfd{1,pfindex}.adata.binSizes')');axis('xy'); title('low prone');
subplot(353);       imagesc(pfd{1}.adata.bins{:},reshape(fpc{3}(:),pfd{1,pfindex}.adata.binSizes')');axis('xy'); title('high prone');
% FSCORE 
sp(end+1)=subplot(356); hold('on');plot(eigScore{pfindex}(:,2),eigScore{pfindex}(:,1),'.'); title('fscore low prone VS rear');      grid('on'); 
sp(end+1)=subplot(357); hold('on');plot(eigScore{pfindex}(:,2),eigScore{pfindex}(:,3),'.'); title('fscore low prone VS high prone');grid('on');
sp(end+1)=subplot(358); hold('on');plot(eigScore{pfindex}(:,3),eigScore{pfindex}(:,1),'.'); title('fscore high prone VS rear');     grid('on');
% ZSCORED - shuffled fscore
sp(end+1)=subplot(3,5,11); hold('on');plot(fsrcz(:,2),fsrcz(:,1),'.'); title('zscore low prone VS rear');      grid('on'); 
sp(end+1)=subplot(3,5,12); hold('on');plot(fsrcz(:,2),fsrcz(:,3),'.'); title('zscore low prone VS high prone');grid('on');
sp(end+1)=subplot(3,5,13); hold('on');plot(fsrcz(:,3),fsrcz(:,1),'.'); title('zscore high prone VS rear');     grid('on');


%cc = tiedrank(eigScore{pfindex}(:,1:3))./numel(unitSubsets{pfindex});
cc = eigScore{pfindex}(:,[2,1,3])./3+0.5;
%cc = tiedrank(FSrC(:,1:3))./numel(unitSubsets{pfindex});
%cc = tiedrank(fsrcz(:,1:3))./numel(unitSubsets{pfindex});
cc(all(abs(fsrcz(:,1:3))<2,2),:) = repmat([0.75,0.75,0.75],[sum(all(abs(fsrcz(:,1:3))<2,2)),1]);


sp(end+1)=subplot(354); hold('on');scatter(mapa(:,1)/10, mapa(:,2)/10,10,cc,'o','filled'); grid('on');
title('tsne on fscores')
xlim([-6,6])
sp(end+1)=subplot(359); hold('on');scatter(mapz(:,1)/10,mapz(:,2)/10,10,cc,'o','filled'); grid('on');
title('tsne on zscores')
%sp(5)=subplot(259); hold('on');scatter(mapa(:,1)/10,mapa(:,2)/10,20,cc,'filled'); grid('on');
%sp(4)=subplot(259); hold('on');plot(map(:,1)/10,map(:,2)/10,'.');grid('on');
%sp(5)=subplot(254); hold('on');plot(mapl(:,1),unity(mapl(:,2)),'.');grid('on');

drawnow();
c = gobjects([1,5]);;
xy = [0,0];
while isempty(hfig.CurrentCharacter)||hfig.CurrentCharacter~='q',
    waitforbuttonpress();
% REMOVE old unit marker
    delete(c);     
% GET current subplot index    
    axind = find(arrayfun(@isequal,sp,repmat([gca],size(sp))));
% GET xy position of currsor on left mouse button down within current subplot    
    xy = sp(axind).CurrentPoint(1,1:2);
% GET axes Data
    axData = [sp(axind).Children.XData',sp(axind).Children.YData'];
% FIND closest point to currsor on left mouse button down
    [~,mind]=min(sqrt(sum(bsxfun(@minus,xy,axData).^2,2)));
% PLOT pft and pfd
    subplot(355);plot(pfd{cluMap(unitSubsets{1}(mind),1),pfindex},cluMap(unitSubsets{1}(mind),2),'mean',true,[],false,'mazeMask',mazeMask);
                 %imagesc(pfd{1}.adata.bins{:},reshape(rmaps(:,unitSubsets{1}(mind)),pfd{tind,pfindex}.adata.binSizes')');axis('xy'); 
    subplot(3,5,10);plot(pft{cluMap(unitSubsets{1}(mind),1)},cluMap(unitSubsets{1}(mind),2),'mean',true,[],true);
    axes(sp(axind));
% HIGHLIGHT selected unit on current subplot
    c(axind) = circle(axData(mind,1),axData(mind,2),0.25,'g');
    for a = find(~ismember(1:numel(sp),axind))
        axes(sp(a));        
        c(a) = circle(sp(a).Children.XData(mind),sp(a).Children.YData(mind),0.25,'k');
    end
end


figure();
hold('on');
for i = 1:3,
    [F,X] = ecdf(abs(fsrcz(:,i)));
    cdfplot(X);
end
legend({'rear','low prone','high prone'},'Location','southeast');
Lines(1.96,[],'k');
Lines(1.96,[],'k');



figure();
hold('on');
for i = 1:3,
    [F,X] = ecdf(fsrcz(:,i));
    cdfplot(X);
end
legend({'rear','low prone','high prone'},'Location','southeast');
Lines(-1.96,[],'k');
Lines(1.96,[],'k');

Lines(-3.1,[],'k');
Lines(3.1,[],'k');


figure,plot(eigScore{pfindex}(eigScore{pfindex}(:,1)<0.5,2),eigScore{pfindex}(eigScore{pfindex}(:,1)<0.5,3),'.'); ...
    title('low prone VS high prone');grid('on');


figure,
for i = 1:5,
subplot(1,5,i); scatter(map(:,1)/10,map(:,2)/10,20,eigScore{1}(:,i),'Filled');
colormap jet
end


%% 
tind = 20;
pfindex = 1;
interpParDfs = struct('bins',{{linspace(-2,2,100)',linspace(-2,2,100)'}},...
                   'nanMaskThreshold', 0,...
                   'methodNanMap',     'linear',...
                   'methodRateMap',    'linear');
feature = fet_HB_pitchB(Trials{tind});

drzPfs = compute_drz(Trials{tind},units{tind},pft{tind});
drzBhv = compute_drz(Trials{tind},units{tind},pfd{tind,pfindex},[],[],[],interpParDfs,feature);

xyz = preproc_xyz(Trials{tind});
xyz.filter('ButFilter',3,2.5,'low');
vxy = xyz.vel({'spine_lower','hcom'});
vxy.data(vxy.data<1e-5) = 1e-5;
vxy.data = log10(vxy.data);
vxy.data = vxy.data.^(1/5);
tper = [Trials{tind}.stc{'theta-groom-sit'}];
tper.resample(xyz);

spk = Trials{tind}.spk;
spk.create(Trials{tind},xyz.sampleRate,tper,units{tind},'deburst');

u = 1;
;

boundaryLimits = [-2.5,2.5];
binDims = [0.1];
SmoothingWeights = [1.4];

fet = xyz.copy('empty');
fet.data = vxy(:,2);

% Try and use MTAApfs not sure if it works for 1D
pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.tag              = 'veltestbody';
pargs.tag              = 'veltesthead';
pargs.units            = units{tind};
pargs.numIter          = 1;
pargs.numIter          = 1001;
pargs.posShuffle       = false;
pargs.halfsample       = false;
pargs.overwrite        = true;
pargs.boundaryLimits   = boundaryLimits;
pargs.binDims          = binDims;
pargs.SmoothingWeights = SmoothingWeights;
pargs.xyzp             = fet.copy();
pargs.bootstrap        = false;
pargs.autoSaveFlag     = false;

drzState = {};
u = 1;        
% SELECT drz periods
dper = MTADepoch([],[],ThreshCross(-0.5<drzPfs(:,u)&drzPfs(:,u)<0.5&-0.5<drzBhv(:,u)&drzBhv(:,u)<0.5,0.5,0),...
                 fet.sampleRate,fet.sync.copy(),fet.origin,'TimePeriods','sts',[],'tdrz','d');
drzState{u} = dper&tper;
pargs.units  = units{tind}(u);
pargs.states = drzState{u};
pfsArgs = struct2varargin(pargs);
pfTemp = MTAApfs(Trials{tind},pfsArgs{:});
pfTemp.save();
for u = 2:numel(units{tind});
    % SELECT drz periods
    dper = MTADepoch([],[],ThreshCross(-0.5<drzPfs(:,u)&drzPfs(:,u)<0.5&-0.5<drzBhv(:,u)&drzBhv(:,u)<0.5,0.5,1),...
                     fet.sampleRate,fet.sync.copy(),fet.origin,'TimePeriods','sts',[],'tdrz','d');
    drzState{u} = dper&tper;
    pargs.units  = units{tind}(u);
    pargs.states = drzState{u};
    pfsArgs = struct2varargin(pargs);
    try
        pfTemp = MTAApfs(pfTemp,pfsArgs{:});            
    catch err,
        disp(err);
        pfTemp = MTAApfs(pfTemp,pfsArgs{:});
    end
end
pfTemp.save();        


pfdbVelB =  MTAApfs(Trials{tind},'tag','veltestbody');
pfdbVel =  MTAApfs(Trials{tind},'tag','veltesthead');
pfdVel =  MTAApfs(Trials{tind},'tag','BSPEEDxHSPEED');

figure();
for u = 1:numel(units{tind}),
subplot(151);
plot(pft{tind},units{tind}(u),1,true,[],true);
subplot(152);
plot(pfd{tind,pfindex},units{tind}(u),1,true,[],false);
subplot(153);
plot(pfd{tind,3},units{tind}(u),1,true,[],false);
subplot(154);
plot(pfdbVelB.adata.bins{:},pfdbVelB.data.rateMap(:,u));
xlim([-2,2]);
subplot(155);
plot(pfdbVel.adata.bins{:},pfdbVel.data.rateMap(:,u));
xlim([-2,2]);
waitforbuttonpress();
end


nb= 25
hl = 2.5
figure,
sts = {'a-m-s','p','w','r','n','s','m'};
for s = 1:numel(sts),
    subplot(3,3,s);
    ind = Trials{tind}.stc{sts{s}};
    %hist2(log10([vxy(ind,:)]),linspace(-2,2,50),linspace(-2,2,50));caxis([0,1000])
    %hist2([vxy(ind,:).^(1/3)],linspace(0,4,50),linspace(0,4,50));caxis([0,1000])    
    hist2([vxy(ind,:).^(1/5)],linspace(0,hl,nb),linspace(0,hl,nb));caxis([0,1000])        
end


pfindex = 1;
interpParDfs = struct('bins',{{linspace(-2,2,100)',linspace(-2,2,100)'}},...
                   'nanMaskThreshold', 0,...
                   'methodNanMap',     'linear',...
                   'methodRateMap',    'linear');

for tind = 1:numel(Trials),

    feature = fet_HB_pitchB(Trials{tind});
    drzPfs = compute_drz(Trials{tind},units{tind},pft{tind});
    drzBhv = compute_drz(Trials{tind},units{tind},pfd{tind,pfindex},[],[],[],interpParDfs,feature);


    % 2d with head and body
    xyz = preproc_xyz(Trials{tind});
    xyz.filter('ButFilter',3,2.5,'low');
    vxy = xyz.vel({'spine_lower','hcom'});
    vxy.data(vxy.data<1e-5) = 1e-5;
    vxy.data = log10(vxy.data);
    ufr = Trials{tind}.ufr.copy();
    ufr.create(Trials{tind},vxy,'theta-groom-sit',units{tind},0.5);
    %vxy.data = vxy.data.^(1/5);
    tper = [Trials{tind}.stc{'theta-groom-sit'}];
    tper.resample(xyz);

    spk = Trials{tind}.spk;
    spk.create(Trials{tind},xyz.sampleRate,tper,units{tind},'deburst');

    u = 1;

    boundaryLimits = [-2,2;-2,2];
    binDims = [0.1,0.1];
    SmoothingWeights = [1.1,1.1];

    fet = xyz.copy('empty');
    fet.data = vxy.data;


    % Try and use MTAApfs not sure if it works for 1D
    pargs = get_default_args('MjgER2016','MTAApfs','struct');        
    pargs.tag              = 'DRZBHV_BSxHS';
    pargs.numIter          = 1001;
    pargs.posShuffle       = false;
    pargs.halfsample       = true;
    pargs.overwrite        = true;
    pargs.boundaryLimits   = boundaryLimits;
    pargs.binDims          = binDims;
    pargs.SmoothingWeights = SmoothingWeights;
    pargs.xyzp             = fet.copy();
    pargs.bootstrap        = false;
    pargs.autoSaveFlag     = false;

    drzState = {};
    u = 1;        
    % SELECT drz periods
    dper = MTADepoch([],[],ThreshCross(-0.5<drzPfs(:,u)&drzPfs(:,u)<0.5&-0.5<drzBhv(:,u)&drzBhv(:,u)<0.5,0.5,0),...
                     fet.sampleRate,fet.sync.copy(),fet.origin,'TimePeriods','sts',[],'tdrz','d');
    drzState{u} = dper&tper;
    pargs.units  = units{tind}(u);
    pargs.states = drzState{u};
    pfsArgs = struct2varargin(pargs);
    pfTemp = MTAApfs(Trials{tind},pfsArgs{:});
    pfTemp.save();
    for u = 2:numel(units{tind});
        % SELECT drz periods
        dper = MTADepoch([],[],ThreshCross(-0.5<drzPfs(:,u)&drzPfs(:,u)<0.5&-0.5<drzBhv(:,u)&drzBhv(:,u)<0.5,0.5,1),...
                         fet.sampleRate,fet.sync.copy(),fet.origin,'TimePeriods','sts',[],'tdrz','d');
        drzState{u} = dper&tper;
        pargs.units  = units{tind}(u);
        pargs.states = drzState{u};
        pfsArgs = struct2varargin(pargs);
        try
            pfTemp = MTAApfs(pfTemp,pfsArgs{:});            
        catch err,
            disp(err);
            pfTemp = MTAApfs(pfTemp,pfsArgs{:});
        end
    end
    pfTemp.save();        
    
% $$$ pfdVel =  MTAApfs(Trials{tind},'tag','BSXHS_test');
% $$$ 
% $$$ figure();
% $$$ for u = 1:numel(units{tind}),
% $$$ subplot(151);
% $$$ plot(pft{tind},units{tind}(u),1,true,[],true);
% $$$ subplot(152);
% $$$ plot(pfd{tind,pfindex},units{tind}(u),1,true,[],false);
% $$$ subplot(153);
% $$$ plot(pfdVel,units{tind}(u),'mean',true,[],false,0.9);
% $$$ subplot(154);
% $$$ plot(pfdbVelB.adata.bins{:},pfdbVelB.data.rateMap(:,u));
% $$$ xlim([-2,2]);
% $$$ subplot(155);
% $$$ plot(pfdbVel.adata.bins{:},pfdbVel.data.rateMap(:,u));
% $$$ xlim([-2,2]);
% $$$ waitforbuttonpress();
% $$$ end



    %for tind = 1:numel(Trials),
    % Try and use MTAApfs not sure if it works for 1D
    pargs = get_default_args('MjgER2016','MTAApfs','struct');        
    pargs.tag              = 'DRZBHV_BSxHS_shuffled';
    pargs.numIter          = 1001;
    pargs.posShuffle       = true;
    pargs.halfsample       = false;
    pargs.overwrite        = true;
    pargs.boundaryLimits   = boundaryLimits;
    pargs.binDims          = binDims;
    pargs.SmoothingWeights = SmoothingWeights;
    pargs.xyzp             = fet.copy();
    pargs.bootstrap        = false;
    pargs.autoSaveFlag     = false;

    drzState = {};
    u = 1;        
    % SELECT drz periods
    dper = MTADepoch([],[],ThreshCross(-0.5<drzPfs(:,u)&drzPfs(:,u)<0.5&-0.5<drzBhv(:,u)&drzBhv(:,u)<0.5,0.5,0),...
                     fet.sampleRate,fet.sync.copy(),fet.origin,'TimePeriods','sts',[],'tdrz','d');
    drzState{u} = dper&tper;
    pargs.units  = units{tind}(u);
    pargs.states = drzState{u};
    pfsArgs = struct2varargin(pargs);
    pfTemp = MTAApfs(Trials{tind},pfsArgs{:});
    pfTemp.save();
    for u = 2:numel(units{tind});
        % SELECT drz periods

        dper = MTADepoch([],[],ThreshCross(-0.5<drzPfs(:,u)&drzPfs(:,u)<0.5&-0.5<drzBhv(:,u)&drzBhv(:,u)<0.5,0.5,1),...
                         fet.sampleRate,fet.sync.copy(),fet.origin,'TimePeriods','sts',[],'tdrz','d');
        drzState{u} = dper&tper;

        pargs.units  = units{tind}(u);
        pargs.states = drzState{u};
        pfsArgs = struct2varargin(pargs);
        try
            pfTemp = MTAApfs(pfTemp,pfsArgs{:});            
        catch err,
            disp(err);
            pfTemp = MTAApfs(pfTemp,pfsArgs{:});
        end
    end
    pfTemp.save();        
end


tind = 20;
pfdPchC = MTAApfs(Trials{tind},'tag','HBPITCHxBPITCH_v7');
pfdPchS = MTAApfs(Trials{tind},'tag','HBPITCHxBPITCH_shuffled');   



pfdVelC = MTAApfs(Trials{tind},'tag','DRZBHV_BSxHS');   
pfdVelS =  MTAApfs(Trials{tind},'tag','DRZBHV_BSxHS_shuffled');


figure();
for u = 1:numel(units{tind}),
subplot(251);
plot(pft{tind},units{tind}(u),1,true,[],true);
subplot(252);
plot(pfdVelC,units{tind}(u),'mean',true,[],false,0.9);
subplot(253);
zsc = (plot(pfdVelC,units{tind}(u),'mean',true,[],false,0.9)-plot(pfdVelS,units{tind}(u),'mean',true,[],false,0.9))./plot(pfdVelS,units{tind}(u),'std',true,[],false,0.9);
imagesc(pfdVelS.adata.bins{:},zsc');
axis('xy');
caxis([2,3]);
colorbar();

subplot(254);
plot(pfdbVelB.adata.bins{:},pfdbVelB.data.rateMap(:,u));
xlim([-2,2]);
subplot(255);
plot(pfdbVel.adata.bins{:},pfdbVel.data.rateMap(:,u));
xlim([-2,2]);

subplot(259);
plot(pfdPchC,units{tind}(u),'mean',true,[],false,0.9);
subplot(2,5,10);
zsc = (plot(pfdPchC,units{tind}(u),'mean',true,[],false,0.9)-plot(pfdPchS,units{tind}(u),'mean',true,[],false,0.9))./plot(pfdPchS,units{tind}(u),'std',true,[],false,0.9);
imagesc(pfdPchS.adata.bins{:},zsc');
axis('xy');
caxis([2,3]);
colorbar();

waitforbuttonpress();
end





% Unit firing rate vs speed
MjgER2016_load_data();

[pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(Trials);

% INITIALIZE output vars
cumMedLocPas = [];
cumPvalLocPas = [];
sampleRate = 8;

% LOOP over all trilas
for tind = 1:numel(Trials),

% LOAD position data
    xyz = preproc_xyz(Trials{tind});

% LOAD instananeous unit firing rates
    ufr = Trials{tind}.ufr.copy();
    ufr.create(Trials{tind},xyz,'theta-groom-sit',units{tind},0.5);
        
% FILTER xyz data lowpass 2.5Hz cutoff
% RESAMPLE xyz data
% SELECT target marker
    xyz.filter('ButFilter',3,2.5,'low');
    xyz.resample(sampleRate);
    xyz.data = sq(xyz(:,'nose',[1,2]));

% RESAMPLE ufr to match xyz
    ufr.resample(xyz);
    
% COMPUTE drz restriction on theta field
    drzPfs = compute_drz(Trials{tind},units{tind},pft{tind},[],[],[],[],xyz);
% COMPUTE drz restriction on bhv field
    feature = fet_HB_pitchB(Trials{tind},sampleRate);
    drzBhv = compute_drz(Trials{tind},units{tind},pfd{tind,pfindex},[],[],[],interpParDfs,feature);

    
% LOAD behavioral states
    tper = [Trials{tind}.stc{'theta-groom-sit'}];
    tper.cast('TimeSeries');
    tper.resample(xyz);            
    lper = resample(cast([Trials{tind}.stc{'loc&gper'}],'TimeSeries'),xyz.sampleRate);
    pper = resample(cast([Trials{tind}.stc{'pause&gper'}] ,'TimeSeries'),xyz.sampleRate);

% INITIALIZE temporary vars
    pvalLocPas = nan([numel(units{tind}),1]);    
    medLocPas  = nan([numel(units{tind}),2]);
    for u = 1:numel(units{tind});
    % SELECT drz periods
        dper = MTADepoch([],[],-0.5<drzPfs(:,u)&drzPfs(:,u)<0.5&-0.5<drzBhv(:,u)&drzBhv(:,u)<0.5,...
                         xyz.sampleRate,...
                         xyz.sync.copy(),...
                         xyz.origin,...
                         'TimeSeries','sts',[],'tdrz','d');


        indLoc = dper&tper.data&lper.data;
        ufrLoc =  ufr(indLoc,u);
        indPas = dper&tper.data&pper.data;    
        ufrPas =  ufr(logical(indPas.data),u);    
    
        try,
        pvalLocPas(u) = ranksum(ufrPas,ufrLoc);
        medLocPas(u,:) = [median(ufrPas),median(ufrLoc)];
        
        end
    end
    cumPvalLocPas = cat(1,cumPvalLocPas,pvalLocPas);
    cumMedLocPas = cat(1,cumMedLocPas,medLocPas);    
end


figure,plot(diff(cumMedLocPas,1,2)./sum(cumMedLocPas,2),log10(cumPvalLocPas),'.')


figure();
cumPvalLocPasSubset =  cumPvalLocPas(unitSubsets{1});
cumMedLocPasSubset  =  cumMedLocPas(unitSubsets{1},:);
ind = cumPvalLocPasSubset<=(0.05/numel(unitSubsets{1}));
%plot(cumMedLocPasSubset(~ind,1),cumMedLocPasSubset(~ind,2),'.');
plot(log10(cumMedLocPasSubset(~ind,1)),log10(cumMedLocPasSubset(~ind,2)),'.');
hold('on');
plot(log10(cumMedLocPasSubset(ind,1)),log10(cumMedLocPasSubset(ind,2)),'.');
%plot(cumMedLocPasSubset(ind,1),cumMedLocPasSubset(ind,2),'.r');
daspect([1,1,1]);
xlim([-1,2]);
ylim([-1,2]);
line([-1,2],[-1,2]);
title({'Median Firing Rate:','Pause Vs Locomotion'})
xlabel('Pause log10(spk/s)');
ylabel('Locomation log10(spk/s)');


eds = linspace(-2,2,50);
figure();
hold('on');
hax = bar(eds,histc(log10(ufrLoc),eds),'histc');
hax.FaceColor = 'c';
hax.EdgeColor = 'c';
hax.FaceAlpha = 0.4;
hax.EdgeAlpha = 0.4;
hax = bar(eds,histc(log10(ufrPas),eds),'histc');
hax.FaceColor = 'r';
hax.EdgeColor = 'r';
hax.FaceAlpha = 0.4;
hax.EdgeAlpha = 0.4;
