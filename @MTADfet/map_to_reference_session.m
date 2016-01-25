function features = map_to_reference_session(features,Trial,RefTrial,varargin)
% function features = normalize_features_to_reference_trial(Trial,features,referenceTrial)
% This Function is only meant for use with fet_tsne at the moment.

[normEachSyncEpoch] = DefaultArgs(varargin,{false},1);

% If Trial is not a MTASession try loading it.
if ischar(Trial),
    Trial = MTATrial(Trial);
elseif iscell(Trial),
    Trial = MTATrial(Trial{:});
end

tempFet = features.copy;
features.resample(Trial.xyz.sampleRate);

% If Trial is not a MTASession try loading it.
if ischar(RefTrial),
    RefTrial = MTATrial(RefTrial);
elseif iscell(RefTrial),
    RefTrial = MTATrial(RefTrial{:});
end
rfet = feval(features.label,RefTrial,features.sampleRate);

switch features.label
  case 'fet_tsne_rev10',
    fetInds =      [   1                    , 7:9                     , 10:14 ];
    stdThresh = cat(2, repmat({10},    1,1) , repmat({.1},        1,3), repmat({10},    1,5));
    diffFun =   cat(2, repmat({@minus},1,1) , repmat({@circ_dist},1,3), repmat({@minus},1,5));
  case 'fet_tsne_rev9',
    fetInds =      [   1                    , 7:8                     , 9:13 ];
    stdThresh = cat(2, repmat({10},    1,1) , repmat({.1},        1,2), repmat({10},    1,5));
    diffFun =   cat(2, repmat({@minus},1,1) , repmat({@circ_dist},1,2), repmat({@minus},1,5));

  case 'fet_tsne_rev8',
    fetInds =      [  1                   ,7:8                  ,9:11 ];
    stdThresh = cat(2,repmat({10},1,1)    ,repmat({.1},1,2)     ,repmat({10},1,3));
    diffFun =   cat(2,repmat({@minus},1,1),repmat({@circ_dist},1,2),repmat({@minus},1,3));

  case 'fet_tsne_rev7',
    fetInds =      [  1                   ,7:8                    ];
    stdThresh = cat(2,repmat({10},1,1)    ,repmat({.1},1,2)        );
    diffFun =   cat(2,repmat({@minus},1,1),repmat({@circ_dist},1,2));

  case 'fet_tsne_rev6',
    fetInds =      [  1                   ,7,8                   ];
    stdThresh = cat(2,repmat({10},1,2)    ,repmat({.1},1,2)        );
    diffFun =   cat(2,repmat({@minus},1,2),repmat({@circ_dist},1,2));

  case 'fet_tsne_rev5',
    fetInds =      [  1                   ,7:8                    ];
    stdThresh = cat(2,repmat({10},1,1)    ,repmat({.1},1,2)        );
    diffFun =   cat(2,repmat({@minus},1,1),repmat({@circ_dist},1,2));

  case 'fet_tsne_rev4',
    fetInds =      [  1:5                 ,13:15                   ,20];
    stdThresh = cat(2,repmat({10},1,5)    ,repmat({.2},1,3)        ,{10});
    diffFun =   cat(2,repmat({@minus},1,5),repmat({@circ_dist},1,3),{@minus});
    %diffFun =   cat(2,repmat({@minus},1,5),repmat({@circ_dist},1,3),{@minus});
  case 'fet_tsne_rev3',
    fetInds =      [  1:5                 ,13:15                  ];
    stdThresh = cat(2,repmat({10},1,5)    ,repmat({.2},1,3)        );
    diffFun =   cat(2,repmat({@minus},1,5),repmat({@circ_dist},1,3));
    %diffFun =   cat(2,repmat({@minus},1,5),repmat({@circ_dist},1,3),{@minus});
  otherwise ,
    fets = [];
end  

[tarMean,tarStd] = mean_embeded_feature_vbvh(features,   Trial,fetInds);
[refMean,refStd] = mean_embeded_feature_vbvh(rfet,    RefTrial,fetInds);

if normEachSyncEpoch,
    inSync = features.sync&features.sync.sync;
    inSync.resample(features);
    inSync = inSync.data'-inSync.data(1)+1;
    inSync(end) = features.size(1);
else
    inSync = [1;features.size(1)];
end


for ind = inSync,
    for f = fetInds;
        nnz = nniz([tarMean{f}(:),refMean{f}(:)])...
              & tarStd{f}(:)<stdThresh{f==fetInds}...
              & tarStd{f}(:)<stdThresh{f==fetInds};
        mzd = diffFun{f==fetInds}(tarMean{f}(:),refMean{f}(:));
        
        mshift = nanmean(mzd(nnz));
        features.data(ind(1):ind(2),f) = bsxfun(diffFun{f==fetInds},...
                                               features.data(ind(1):ind(2),f),...
                                               mshift);
    end
end
features.resample(tempFet);