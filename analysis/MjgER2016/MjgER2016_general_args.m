function MjgER2016_general_args(section)
% MjgER2016_figure6_args



clearvars('-GLOBAL','AP');

global AP

% GENERAL ------------------------------------------------------------------------------------------
sessionListName = 'MjgER2016';
activeState     = 'theta-groom-sit';

AP.fet_HB_pitchB.referenceTrial = 'Ed05-20140529.ont.all';
%---------------------------------------------------------------------------------------------------



switch section
    
% SECTION 1 ------------------------------------------------------------------------------------------    
  case 'section 1',
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
AP.compute_bhv_ratemaps_erpPCA =                                                                 ...
    struct('numComp',                   5,                                                       ...
           'range',                     [250,600]                                                ...
          );
%---------------------------------------------------------------------------------------------------


% compute_bhv_ratemap_nnmf -------------------------------------------------------------------------
AP.compute_bhv_ratemaps_nnmf =                                                                   ...
    struct('numComp',                   5,                                                       ...
           'range',                     [250,600]                                                ...
          );
%---------------------------------------------------------------------------------------------------


end