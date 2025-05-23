function MjgER2016_figure6_args_alt1(section)
% MjgER2016_figure6_args



clearvars('-GLOBAL','AP');

global AP

% GENERAL ------------------------------------------------------------------------------------------
sessionListName = 'MjgER2016';
activeState     = 'theta-groom-sit';

AP.fet_HB_pitchB.referenceTrial = 'Ed05-20140529.ont.all';
%---------------------------------------------------------------------------------------------------

switch section
  case 'section 1',

% compute_bhv_ratemaps -----------------------------------------------------------------------------
% $$$ AP.compute_bhv_ratemaps =                                                                        ...
% $$$     struct('get_featureSet',            @fet_HB_pitchB,                                          ...
% $$$            'sampleRate',                16,                                                      ...
% $$$            'pfsArgs',                   struct('states',           activeState,                  ...
% $$$                                                'binDims',          [0.2,0.2],                    ...
% $$$                                                'SmoothingWeights', [1.2,1.2],                    ...
% $$$                                                'numIter',          1,                            ...
% $$$                                                'boundaryLimits',   [-2,0.8;-0.8,2],              ...
% $$$                                                'halfsample',       false),                       ...
% $$$            'threshRate',                0.8,                                                     ...
% $$$            'threshDist',                250                                                      ...
% $$$            );
AP.compute_bhv_ratemaps =                                                                        ...
    struct('get_featureSet',            @fet_HB_pitchB,                                          ...
           'sampleRate',                16,                                                      ...
           'pfsArgs',                   struct('states',           activeState,                  ...
                                               'binDims',          [0.1,0.1],                    ...
                                               'SmoothingWeights', [1.8,1.8],                    ...
                                               'numIter',          1,                            ...
                                               'boundaryLimits',   [-2,0.8;-0.8,2],              ...
                                               'halfsample',       false),                       ...
           'threshRate',                0.8,                                                     ...
           'threshDist',                250                                                      ...
           );
%---------------------------------------------------------------------------------------------------



% compute_bhv_ratemaps_shuffled --------------------------------------------------------------------
% $$$ AP.compute_bhv_ratemaps_shuffled =                                                               ...
% $$$     struct('get_featureSet',            @fet_HB_pitchB,                                          ...
% $$$            'sampleRate',                16,                                                      ...
% $$$            'pfsArgs',                   struct('states',           activeState,                  ...
% $$$                                                'binDims',          [0.2,0.2],                    ...
% $$$                                                'SmoothingWeights', [1.2,1.2],                    ...
% $$$                                                'numIter',          1001,                         ...
% $$$                                                'posShuffle',       true,                         ...
% $$$                                                'boundaryLimits',   [-2,0.8;-0.8,2],              ...
% $$$                                                'halfsample',       false),                       ...
% $$$            'threshRate',                0.8,                                                     ...
% $$$            'threshDist',                250                                                      ...
% $$$            );
AP.compute_bhv_ratemaps_shuffled =                                                               ...
    struct('get_featureSet',            @fet_HB_pitchB,                                          ...
           'sampleRate',                16,                                                      ...
           'pfsArgs',                   struct('states',           activeState,                  ...
                                               'binDims',          [0.1,0.1],                    ...
                                               'SmoothingWeights', [1.8,1.8],                    ...
                                               'numIter',          1001,                         ...
                                               'posShuffle',       true,                         ...
                                               'boundaryLimits',   [-2,0.8;-0.8,2],              ...
                                               'halfsample',       false),                       ...
           'threshRate',                0.8,                                                     ...
           'threshDist',                250                                                      ...
           );
%---------------------------------------------------------------------------------------------------



% compute_bhv_ratemap_erpPCA -----------------------------------------------------------------------
% $$$ AP.comput_bhv_ratemaps_erpPCA =                                                                  ...
% $$$     struct('numComp',                   5,                                                       ...
% $$$            'range',                     [120,356]                                                ...
% $$$           );
AP.comput_bhv_ratemaps_erpPCA =                                                                  ...
    struct('numComp',                   5,                                                       ...
           'range',                     [120,356]                                                ...
          );
%---------------------------------------------------------------------------------------------------



% compute_bhv_ratemap_nnmf -------------------------------------------------------------------------
% $$$ AP.comput_bhv_ratemaps_nnmf =                                                                    ...
% $$$     struct('numComp',                   5,                                                       ...
% $$$            'range',                     [120,356]                                                ...
% $$$           );
AP.comput_bhv_ratemaps_erpPCA =                                                                  ...
    struct('numComp',                   5,                                                       ...
           'range',                     [120,356]                                                ...
          );
%---------------------------------------------------------------------------------------------------





  case 'section 2', %*******************************************************************************

% compute_bhv_ratemaps -----------------------------------------------------------------------------
% $$$ AP.compute_bhv_ratemaps =                                                                        ...
% $$$     struct('get_featureSet',            @fet_xyhb,                                               ...
% $$$            'sampleRate',                16,                                                      ...
% $$$            'pfsArgs',                   struct('states',           activeState,                  ...
% $$$                                                'binDims',          [50,50,0.2,0.2],              ...
% $$$                                                'SmoothingWeights', [1.2,1.2,1.2,1.2],            ...
% $$$                                                'numIter',          1,                            ...
% $$$                                                'boundaryLimits',   [-500,500;-500,500;           ...
% $$$                                                                     -2,0.8;-0.8,2],              ...
% $$$                                                'halfsample',       false),                       ...
% $$$            'threshRate',                0.8,                                                     ...
% $$$            'threshDist',                250                                                      ...
% $$$            );
% $$$ AP.compute_bhv_ratemaps =                                                                        ...
% $$$     struct('get_featureSet',            @fet_xyhb,                                               ...
% $$$            'sampleRate',                16,                                                      ...
% $$$            'pfsArgs',                   struct('states',           activeState,                  ...
% $$$                                                'binDims',          [100,100,0.2,0.2],            ...
% $$$                                                'SmoothingWeights', [1.2,1.2,1.2,1.2],            ...
% $$$                                                'numIter',          1,                            ...
% $$$                                                'boundaryLimits',   [-500,500;-500,500;           ...
% $$$                                                                     -2,0.8;-0.8,2],              ...
% $$$                                                'halfsample',       false),                       ...
% $$$            'threshRate',                0.8,                                                     ...
% $$$            'threshDist',                250                                                      ...
% $$$            );
AP.compute_bhv_ratemaps =                                                                        ...
    struct('get_featureSet',            @fet_xyhb,                                               ...
           'sampleRate',                16,                                                      ...
           'pfsArgs',                   struct('states',           activeState,                  ...
                                               'binDims',          [20,20,0.1,0.1],              ...
                                               'SmoothingWeights', [4.8,4.8,1.8,1.8],            ...
                                               'numIter',          1,                            ...
                                               'boundaryLimits',   [-500,500;-500,500;           ...
                                                                    -2,0.8;-0.8,2],              ...
                                               'halfsample',       false),                       ...
           'threshRate',                0.8,                                                     ...
           'threshDist',                250                                                      ...
           );
%---------------------------------------------------------------------------------------------------

  case 'section 3', %*******************************************************************************

% compute_bhv_ratemaps -----------------------------------------------------------------------------
AP.compute_bhv_ratemaps =                                                                        ...
    struct('get_featureSet',            @fet_xyhb,                                               ...
           'sampleRate',                16,                                                      ...
           'pfsArgs',                   struct('states',           activeState,                  ...
                                               'binDims',          [50,50,0.2,0.2],              ...
                                               'SmoothingWeights', [1.2,1.2,1.2,1.2],            ...
                                               'numIter',          1,                            ...
                                               'boundaryLimits',   [-500,500;-500,500;           ...
                                                                    -2,0.8;-0.8,2],              ...
                                               'halfsample',       false),                       ...
           'threshRate',                0.8,                                                     ...
           'threshDist',                250                                                      ...
           );
%---------------------------------------------------------------------------------------------------

end