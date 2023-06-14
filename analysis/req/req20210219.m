% req20210219
%     Tags: rhm bfs
%     Status: Active
%     Type: Analysis
%     Author: Justin Graboski
%     Final_Forms: NA
%     Project: General
%     Description:  This script was made to examine the posibility of using only the theta-trough spikes for theta place field computation.




MjgER2016_load_data();
configure_default_args();


overwrite = false;
pft = cf(@(t,u)  pfs_2d_theta(t,u,'overwrite',overwrite),          Trials, units);
pfs = cf(@(t,u)  pfs_2d_states(t,u,'overwrite',overwrite),         Trials, units);
bfs = cf(@(t,u)  compute_bhv_ratemaps(t,u,'overwrite',overwrite),  Trials, units);
pftt = cf(@(t,u,tr,pc)                                             ...
          pfs_2d_thetaT(t,u,tr,pc,'overwrite',overwrite),          ...
          Trials, units, num2cell([sessionList.thetaRefGeneral]), num2cell(phzCorrection));


t = 18;
[mxr,mxp] = maxRate(pft{t},units{t});
[mxrt,mxpt] = maxRate(pftt{t},units{t});
figure;
for u = units{t},
    subplot(121);    
    plot(pft{t},u,1,'text');
    Lines(mxp(u==units{t},1),[],'k');
    Lines([],mxp(u==units{t},2),'k');    
    
    Lines(mxpt(u==units{t},1),[],'m');
    Lines([],mxpt(u==units{t},2),'m');    
    
    subplot(122);
    plot(pftt{t},u,1,'text');    
    Lines(mxpt(u==units{t},1),[],'k');
    Lines([],mxpt(u==units{t},2),'k');    
    title(num2str(u));
    waitforbuttonpress();
end


% This script was made to examine the posibility of using only the theta-trough spikes for theta place
% field computation.

% This may mean new place field centers (PFCs) for each unit which may result in better decoding 
% BUT Decoding may not be determined by rate probability.
% ALTERNATIVE : use mode PFC of the active placecell set using the mean-shift procedure by Fukunaga and hostetler (1975)
% MODIFICATION : use aformentioned mode with a weighted position history which has "momentum" -> refine definition

% EXAMINATION OF THEORY - bayesian decoding 
% IF the phase code and rate code are independent then bayesian decoding in an improper estimator of position. 
% At the network level, what is the implication for phase precession if it depends on firing rate?



%figure,plot(mxrt-mxr,vecnorm(mxp-mxpt,2,2),'.');XB

% 20210331
% NOTE - the decoding computations need to be recomptuted for all sessions

% start with the theta state

% 20210402
% 
%
% p(s|x)
% p(s|b)
% p(s|x,b)


clearvars('-GLOBAL','AP');
activeState = 'theta-groom-sit';

stc = copy(Trial.stc);


xyz = preproc_xyz(Trial,'trb');
xyz.resample(2);

myInts = unitsInts{20};
nUnits = numel(myInts);

% MASK POS
width = pft{1}.adata.binSizes(1);
height = pft{1}.adata.binSizes(2);
radius = round(width/2)-find(pft{1}.adata.bins{1}<=-450,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
maskPos = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);


% MASK BHV
bfs = cf(@(t,u)  compute_bhv_ratemaps(t,u,'overwrite',false),            Trials, units);
[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfs, units, [], [], false);
maskBhv = reshape(validDims,bfs{1}.adata.binSizes');

% SUBSET of xyz within active stae
sxyz = sq(xyz([stc{activeState}],'hcom',[1,2]));

ufr = Trial.load('ufr',xyz,[],unitsInts{20},0.5);

% p(s)
pfsArgs = struct('states',           activeState,                  ...
                 'binDims',          [1000,1000],                  ...
                 'SmoothingWeights', [],                           ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-500,500;-500,500],          ...
                 'halfsample',       false);
pfi_ps = compute_ratemaps(Trial,unitsInts{20},[],[],pfsArgs);

rss = {};
ufrstate = ufr([stc{activeState}],:);
nSamples = size(ufrstate,1);
rss{1} = sqrt(sum(bsxfun(@minus,ufrstate,pfi_ps.data.rateMap).^2)./size(ufrstate,1));


%% p(s|x) --------------------------------------------------------------------------------
pfsArgs = struct('states',           activeState,                  ...
                 'binDims',          [20,20],                      ...
                 'SmoothingWeights', [3.5,3.5],                    ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-500,500;-500,500],          ...
                 'halfsample',       false);
pfi_psgx = compute_ratemaps(Trial,myInts,@fet_xy,[],pfsArgs,'overwrite',true);

fet = fet_xy(Trial);
fet.resample(2);

exRate = get_expected_rate(pfi_psgx,myInts,fet,[stc{activeState}],maskPos);

rss{2} = sqrt(bsxfun(@rdivide,sum((ufr([stc{activeState}],:)-exRate).^2,'omitnan'),sum(~isnan(exRate))));


%% p(s|b) --------------------------------------------------------------------------------
pfsArgs = struct('states',           activeState,                  ...
                 'binDims',          [0.1,0.1],                    ...
                 'SmoothingWeights', [1.8,1.8],                    ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-2,0.8;-0.8,2],              ...
                 'halfsample',       false,                        ...
                 'compute_pfs',      @PlotPF_Int);
pfi_psgb = compute_bhv_ratemaps(Trial,myInts,[],[],[],pfsArgs,'threshRate',1,'threshDist',inf);

% SUBSET of xyz within active stae
fet = fet_HB_pitchB(Trial);
fet.resample(2);

exRate = get_expected_rate(pfi_psgb,myInts,fet,[stc{activeState}],maskBhv);

rss{3} = sqrt(bsxfun(@rdivide,sum((ufr([stc{activeState}],:)-exRate).^2,'omitnan'),sum(~isnan(exRate))));


%% p(s|xb) --------------------------------------------------------------------------------
pfsArgs = struct('states',           activeState,                       ...
                 'binDims',          [20,20,0.1,0.1],                   ...
                 'SmoothingWeights', [2.5,2.5,1.8,1.8],                 ...
                 'numIter',          1,                            ...
                 'boundaryLimits',   [-500,500;-500,500;-2,0.8;-0.8,2],...
                 'halfsample',       false,                        ...
                 'compute_pfs',      @PlotPF_Int);
pfi_psgxb = compute_ratemaps(Trial,myInts,@fet_xyhb,[],pfsArgs);

% SUBSET of xyz within active stae
fet = fet_xyhb(Trial);
fet.resample(2);

exRate = get_expected_rate(pfi_psgxb,myInts,fet,[stc{activeState}]);

rss{4} = sqrt(bsxfun(@rdivide,sum((ufr([stc{activeState}],:)-exRate).^2,'omitnan'),sum(~isnan(exRate))));


hfig = figure();










% $$$ edx = [pfi_psgx.adata.bins{1}-pfi_psgx.parameters.binDims(1)/2;pfi_psgx.adata.bins{1}(end)+pfi_psgx.parameters.binDims(1)/2];
% $$$ edy = [pfi_psgx.adata.bins{2}-pfi_psgx.parameters.binDims(2)/2;pfi_psgx.adata.bins{2}(end)+pfi_psgx.parameters.binDims(2)/2];

