function features = map_to_reference_session(features,Trial,RefTrial,varargin)
% function features = normalize_features_to_reference_trial(Trial,features,referenceTrial)
% This Function is only meant for use with fet_tsne at the moment.

[normEachSyncEpoch] = DefaultArgs(varargin,{false},1);
% If Trial is not a MTASession try loading it.
Trial = MTATrial.validate(Trial);

tempFet = features.copy;
features.resample(Trial.xyz.sampleRate);

RefTrial = MTATrial.validate(RefTrial);
if ~strcmp(features.label,'fet_all'), % Sorry this exists ...
    rfet = feval(features.label,RefTrial,features.sampleRate);
end

switch features.label
  case 'fet_all'

    xyz = preproc_xyz(RefTrial,'spline_spine');

    rfet = xyz.copy;
    rfet.data = rfet(:,:,3);

    fetInds = 1:rfet.size(2);
    stdThresh = repmat({10},1,rfet.size(2));
    diffFun =   repmat({@minus},1,rfet.size(2));
    
  case 'fet_tsne_rev15',
    fetInds =      [   1:4                    , 7:9                  ];
    stdThresh = cat(2, repmat({10},    1,4)   , repmat({.1},      1,3));
    diffFun =   cat(2, repmat({@minus},1,4) , repmat({@circ_dist},1,3));

  case 'fet_tsne_rev14',
    fetInds =    4:6;
    stdThresh =  repmat({.1},      1,3);
    diffFun =    repmat({@circ_dist},1,3);
  
  case 'fet_tsne_rev13',
    fetInds =      [   1:4                    , 8:10                  ];
    stdThresh = cat(2, repmat({10},    1,4)   , repmat({.1},      1,3));
    diffFun =   cat(2, repmat({@minus},1,4) , repmat({@circ_dist},1,3));
% $$$   case 'fet_tsne_rev13',
% $$$     fetInds =      [   1:4                    , 8:10                  ,[1:4]+18                , [8:10]+18               ,[1:4]+18*2                , [8:10]+18*2                  ];
% $$$     stdThresh = cat(2, repmat({10},    1,4)   , repmat({.1},      1,3),repmat({10},    1,4)   , repmat({.1},      1,3),repmat({10},    1,4)   , repmat({.1},      1,3));
% $$$     diffFun =   cat(2, repmat({@minus},1,4) , repmat({@circ_dist},1,3),repmat({@minus},1,4) , repmat({@circ_dist},1,3),repmat({@minus},1,4) , repmat({@circ_dist},1,3));

  case 'fet_tsne_rev12',
    fetInds =      [   1:3                    , 9:11                     , 12:16 ];
    stdThresh = cat(2, repmat({10},    1,3) , repmat({.1},        1,3), repmat({10},    1,5));
    diffFun =   cat(2, repmat({@minus},1,3) , repmat({@circ_dist},1,3), repmat({@minus},1,5));

  case 'fet_tsne_rev11',
    fetInds =      [   1                    , 7:9                     , 10:14 ];
    stdThresh = cat(2, repmat({10},    1,1) , repmat({.1},        1,3), repmat({10},    1,5));
    diffFun =   cat(2, repmat({@minus},1,1) , repmat({@circ_dist},1,3), repmat({@minus},1,5));

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