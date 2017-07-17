function features = map_to_reference_session(features,Trial,RefTrial,varargin)
%function features = map_to_reference_session(features,Trial,RefTrial,varargin)
%[normEachSyncEpoch,verbose] = DefaultArgs(varargin,{false,false},1);

[normEachSyncEpoch,verbose] = DefaultArgs(varargin,{false,false},1);
% If Trial is not a MTASession try loading it.
Trial = MTATrial.validate(Trial);

RefTrial = MTATrial.validate(RefTrial);
if strcmp(Trial.filebase,RefTrial.filebase),
    return
end


tempFet = features.copy;
features.resample(Trial.xyz.sampleRate);
minimumOccupancy = 1; % seconds

rfet = feval(features.label,RefTrial,features.sampleRate);

switch features.label
  case 'fet_bref_emb'
    fetInds = [];
    for i = 0:(size(features,2)/30)-1,
        fetInds = [fetInds,[1:15]+30*i];
    end
% $$$     stdThresh = repmat({30},1,10);
% $$$     kurThresh = repmat({20},1,10);
    diffFun   = repmat({@minus},1,numel(fetInds));
  
  case 'fet_bref_rev4'
    fetInds   = [1:9];
    diffFun   = repmat({@minus},1,9);

  case 'fet_bref_rev3'
    fetInds   = [1:12];
    diffFun   = repmat({@minus},1,12);
  
  case 'fet_bref_rev2'
    fetInds   = [1:12];
% $$$     stdThresh = repmat({30},1,10);
% $$$     kurThresh = repmat({20},1,10);
    diffFun   = repmat({@minus},1,12);

  case 'fet_bref'
    fetInds   = [1:15];
% $$$     stdThresh = repmat({30},1,10);
% $$$     kurThresh = repmat({20},1,10);
    diffFun   = repmat({@minus},1,15);
  case 'fet_bref_exp'
    fetInds   = [1:15];
% $$$     stdThresh = repmat({30},1,10);
% $$$     kurThresh = repmat({20},1,10);
    diffFun   = repmat({@minus},1,15);

  case 'fet_svd_mis'
    fetInds =      [   1:5                      , 7:10                ];
    stdThresh = cat(2, repmat({.2},    1,5)     , repmat({15},    1,4));
    kurThresh = cat(2, repmat({5},     1,5)     , repmat({15},    1,4));
    diffFun =   cat(2, repmat({@circ_dist},1,5) , repmat({@minus},1,4));
    
  case 'fet_rear_desc'
    fetInds =      [   1:4                      , 5:11                ];
    stdThresh = cat(2, repmat({.2},    1,4)     , repmat({15},    1,7));
    kurThresh = cat(2, repmat({5},     1,4)     , repmat({15},    1,7));
    diffFun =   cat(2, repmat({@circ_dist},1,4) , repmat({@minus},1,7));

  case 'fet_all2'
    fetInds =      [1:12,[40:44]];
    stdThresh = cat(2,repmat({35},1,12),repmat({.2},1,12),repmat({15},1,5));
    diffFun = cat(2, ...
              repmat({@minus},1,12),...      % 1:12
              repmat({@circ_dist},1,12),...  %13:24
              repmat({@minus},1,5)...       %52:56
    );    
    
  case 'fet_raw'
    fetInds =      [1:12,13:24,52:56];
    stdThresh = cat(2,repmat({35},1,12),repmat({.2},1,12),repmat({15},1,5));
    diffFun = cat(2, ...
              repmat({@minus},1,12),...      % 1:12
              repmat({@circ_dist},1,12),...  %13:24
              repmat({@minus},1,5)...       %52:56
    );
    
  case 'fet_head_pitch'
    fetInds   = [1];
    stdThresh = {0.1};
    diffFun   = {@circ_dist};
    
  case 'fet_mis'
    fetInds =      [   1:5                      , 7:10                ];
    stdThresh = cat(2, repmat({.2},    1,5)     , repmat({15},    1,4));
    kurThresh = cat(2, repmat({5},     1,5)     , repmat({15},    1,4));
    diffFun =   cat(2, repmat({@circ_dist},1,5) , repmat({@minus},1,4));

  case 'fet_all'

    xyz = preproc_xyz(RefTrial,'spline_spine');

    rfet = xyz.copy;
    rfet.data = xyz(:,{'spine_lower','pelvis_root','spine_middle','spine_upper',...
                       'head_back',  'head_left'  ,'head_front',  'head_right',...
                       'bcom',       'hcom',       'acom'},3);

    fetInds = 1:rfet.size(2);
    stdThresh = repmat({15},1,rfet.size(2));
    diffFun =   repmat({@minus},1,rfet.size(2));
    
  case 'fet_tsne_rev15',
    fetInds =      [   1:4                    , 7:9                  ];
    stdThresh = cat(2, repmat({15},    1,4)   , repmat({.1},      1,3));
    diffFun =   cat(2, repmat({@minus},1,4) , repmat({@circ_dist},1,3));

  case 'fet_tsne_rev14',
    fetInds =    4:6;
    stdThresh =  repmat({.1},      1,3);
    diffFun =    repmat({@circ_dist},1,3);
  
  case 'fet_tsne_rev13',
    fetInds =      [   1:4                    , 8:10                  ];
    stdThresh = cat(2, repmat({15},    1,4)   , repmat({.1},      1,3));
    diffFun =   cat(2, repmat({@minus},1,4) , repmat({@circ_dist},1,3));
% $$$   case 'fet_tsne_rev13',
% $$$     fetInds =      [   1:4                    , 8:10                  ,[1:4]+18                , [8:10]+18               ,[1:4]+18*2                , [8:10]+18*2                  ];
% $$$     stdThresh = cat(2, repmat({15},    1,4)   , repmat({.1},      1,3),repmat({15},    1,4)   , repmat({.1},      1,3),repmat({15},    1,4)   , repmat({.1},      1,3));
% $$$     diffFun =   cat(2, repmat({@minus},1,4) , repmat({@circ_dist},1,3),repmat({@minus},1,4) , repmat({@circ_dist},1,3),repmat({@minus},1,4) , repmat({@circ_dist},1,3));

  case 'fet_tsne_rev12',
    fetInds =      [   1:3                    , 9:11                     , 12:16 ];
    stdThresh = cat(2, repmat({15},    1,3) , repmat({.1},        1,3), repmat({15},    1,5));
    diffFun =   cat(2, repmat({@minus},1,3) , repmat({@circ_dist},1,3), repmat({@minus},1,5));

  case 'fet_tsne_rev11',
    fetInds =      [   1                    , 7:9                     , 10:14 ];
    stdThresh = cat(2, repmat({15},    1,1) , repmat({.1},        1,3), repmat({15},    1,5));
    diffFun =   cat(2, repmat({@minus},1,1) , repmat({@circ_dist},1,3), repmat({@minus},1,5));

  case 'fet_tsne_rev10',
    fetInds =      [   1                    , 7:9                     , 10:14 ];
    stdThresh = cat(2, repmat({15},    1,1) , repmat({.1},        1,3), repmat({15},    1,5));
    diffFun =   cat(2, repmat({@minus},1,1) , repmat({@circ_dist},1,3), repmat({@minus},1,5));
  case 'fet_tsne_rev9',
    fetInds =      [   1                    , 7:8                     , 9:13 ];
    stdThresh = cat(2, repmat({15},    1,1) , repmat({.1},        1,2), repmat({15},    1,5));
    diffFun =   cat(2, repmat({@minus},1,1) , repmat({@circ_dist},1,2), repmat({@minus},1,5));

  case 'fet_tsne_rev8',
    fetInds =      [  1                   ,7:8                  ,9:11 ];
    stdThresh = cat(2,repmat({15},1,1)    ,repmat({.1},1,2)     ,repmat({15},1,3));
    diffFun =   cat(2,repmat({@minus},1,1),repmat({@circ_dist},1,2),repmat({@minus},1,3));

  case 'fet_tsne_rev7',
    fetInds =      [  1                   ,7:8                    ];
    stdThresh = cat(2,repmat({15},1,1)    ,repmat({.1},1,2)        );
    diffFun =   cat(2,repmat({@minus},1,1),repmat({@circ_dist},1,2));

  case 'fet_tsne_rev6',
    fetInds =      [  1                   ,7,8                   ];
    stdThresh = cat(2,repmat({15},1,2)    ,repmat({.1},1,2)        );
    diffFun =   cat(2,repmat({@minus},1,2),repmat({@circ_dist},1,2));

  case 'fet_tsne_rev5',
    fetInds =      [  1                   ,7:8                    ];
    stdThresh = cat(2,repmat({15},1,1)    ,repmat({.1},1,2)        );
    diffFun =   cat(2,repmat({@minus},1,1),repmat({@circ_dist},1,2));

  case 'fet_tsne_rev4',
    fetInds =      [  1:5                 ,13:15                   ,20];
    stdThresh = cat(2,repmat({15},1,5)    ,repmat({.2},1,3)        ,{15});
    diffFun =   cat(2,repmat({@minus},1,5),repmat({@circ_dist},1,3),{@minus});
    %diffFun =   cat(2,repmat({@minus},1,5),repmat({@circ_dist},1,3),{@minus});
  case 'fet_tsne_rev3',
    fetInds =      [  1:5                 ,13:15                  ];
    stdThresh = cat(2,repmat({15},1,5)    ,repmat({.2},1,3)        );
    diffFun =   cat(2,repmat({@minus},1,5),repmat({@circ_dist},1,3));
    %diffFun =   cat(2,repmat({@minus},1,5),repmat({@circ_dist},1,3),{@minus});
  otherwise ,
    fets = [];
end  

% This section was a waste of time
% $$$ targetFeatureDomainBoundaries = prctile(features.data(nniz(features.data),fetInds'),[5,95]);
% $$$ referencFeatureDomainBoundaries = prctile(rfet.data(nniz(rfet.data),fetInds'),[5,95]);
% $$$ featureDomainBoundaries = [min([targetFeatureDomainBoundaries(1,:);referencFeatureDomainBoundaries(1,:)]);...
% $$$                            max([targetFeatureDomainBoundaries(2,:);referencFeatureDomainBoundaries(2,:)])]';
% $$$ [tarMean,tarStd,tarCnt,tarDom] = mean_embeded_feature_vbvh(features,   Trial,fetInds,featureDomainBoundaries);
% $$$ [refMean,refStd,refCnt,refDom] = mean_embeded_feature_vbvh(rfet,    RefTrial,fetInds,featureDomainBoundaries);
% -----
%[tarMean,tarStd,tarKur] = mean_embeded_feature_vbvh(features,   Trial,fetInds,[],verbose);
%[refMean,refStd,refKur] = mean_embeded_feature_vbvh(rfet,    RefTrial,fetInds,[],verbose);
%[tarMean,tarStd,tarKur,tarCnt] = mean_embeded_feature_vbvhzbzh(features,   Trial,fetInds,[],verbose);
%[refMean,refStd,refKur,refCnt] = mean_embeded_feature_vbvhzbzh(rfet,    RefTrial,fetInds,[],verbose);
[tarMean,tarStd,tarKur,tarCnt] = mean_embeded_feature_vbvhza(features,   Trial,fetInds,[],verbose);
[refMean,refStd,refKur,refCnt] = mean_embeded_feature_vbvhza(rfet,    RefTrial,fetInds,[],verbose);



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
        nnz = nniz([tarMean{f}(:),refMean{f}(:)])&...
              tarCnt{f}(:)>minimumOccupancy*features.sampleRate&...
              refCnt{f}(:)>minimumOccupancy*features.sampleRate;%,...
% $$$               & tarStd{f}(:)<stdThresh{f==fetInds}...
% $$$               & refStd{f}(:)<stdThresh{f==fetInds}...
% $$$               & tarKur{f}(:)<kurThresh{f==fetInds}...
% $$$               & refKur{f}(:)<kurThresh{f==fetInds};

% $$$         nnz = nniz([tarMean{f}(:),refMean{f}(:)])...
% $$$               & tarStd{f}(:)<stdThresh{f==fetInds}...
% $$$               & refStd{f}(:)<stdThresh{f==fetInds};
        mzd = diffFun{f==fetInds}(tarMean{f}(:),refMean{f}(:));
        
        mshift = nanmedian(mzd(nnz));

        sind = ind(1):ind(2);
        zind = features.data(sind,f)==0;
        features.data(sind,f) = bsxfun(diffFun{f==fetInds},...
                                       features.data(ind(1):ind(2),f),...
                                       mshift);

        features.data(sind(zind),f) =0; 
        
    end
end
features.resample(tempFet);


