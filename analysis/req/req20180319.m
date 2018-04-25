


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

[pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection] = req20180123_ver5(Trials,'overwriteErpPCAFlag',true);

pfdShuffled = {};



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



figure,
u = 36;
subplot(131);
plot(pfd{tind},        u,'mean',true,[],false,0.5,false,[],@jet);
subplot(132);
plot(pfdShuffled{tind},u,'mean',true,[],false,0.5,false,[],@jet);
subplot(133);
plot(pfdShuffled{tind},u,'std', true,[],false,0.5,false,[],@jet);




reshape_eigen_vector = @(V,pfd) reshape(V(:),pfd{1}.adata.binSizes')';




numComp = 5;
pfindex = 1;
fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan(size(validDims{pfindex}));
    fpc{i}(validDims{pfindex}) = eigVec{pfindex}(:,i)./nansum(eigVec{pfindex}(:,i));
    fpc{i}(validDims{pfindex}) = fpc{i}(validDims{pfindex})-mean(fpc{i}(validDims{pfindex}));
end


figure
subplot(121); imagesc(reshape(fpc{eigInd}(:),pfd{tind,pfindex}.adata.binSizes')');axis('xy');
subplot(122); imagesc(reshape(pfd{tind,pfindex}.data.rateMap(:,pfd{tind,pfindex}.data.clu==u,1),pfd{tind,pfindex}.adata.binSizes')');axis('xy');


figure,imagesc(reshape_eigen_vector(fpc{1},pfd(1,pfindex)));axis('xy');

reshape(,pfTemp.adata.binSizes')'.*reshape_eigen_vector(V,pfd(1,pfindex));



% Compute z-score of randomly shuffled spike time chunks 
%               projected on each eigen vector, as computed by 
%               req20180123_ver5.m% Score from


score = [];
for u = pfdShuffled{tind}.data.clu,
    cpi = ~isnan(pfd{tind,pfindex}.data.rateMap(:,pfd{tind,pfindex}.data.clu==u,1));    
    for eigInd = 1:3,
        pfdistrib = sq(nansum(bsxfun(@times,fpc{eigInd}(cpi),pfdShuffled{tind}.data.rateMap(cpi,pfdShuffled{tind}.data.clu==u,:))))';

        es = sq(nansum(fpc{eigInd}(cpi).*pfd{tind,pfindex}.data.rateMap(cpi,pfd{tind,pfindex}.data.clu==u,1)))';

        score(u==units{tind},eigInd) = (es-mean(pfdistrib))./std(pfdistrib);
    end
end

figure,hist(pfdistrib,100);

for eigInd = 1:3,
    eigScore{pfindex}(find(unitSubsets{pfindex}==find(ismember(cluMap,[tind,u],'rows'))),eigInd)
end

(eigScore{pfindex}(find(unitSubsets{pfindex}==find(ismember(cluMap,[tind,u],'rows'))),eigInd)-mean(pfdistrib))./std(pfdistrib)


bsxfun(@times,pfTemp.data.rateMap(:,1,:),eigVec{1}(:,1)');