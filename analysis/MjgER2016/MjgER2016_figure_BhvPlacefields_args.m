% MjgER2016_figure_BhvPlacefields_args.m
% SET global overide pramaters for the listed functions.
% - fet_HB_pitchB
% - compute_bhv_ratemaps
% - compute_bhv_ratemaps_shuffled
% - compute_ratemaps
% - compute_bhv_ratemap_erpPCA

clearvars('-GLOBAL','AP');

global AP

% GENERAL ------------------------------------------------------------------------------------------
sessionListName = 'MjgER2016';
activeState     = 'theta-groom-sit';
AP.fet_HB_pitchB.referenceTrial = 'Ed05-20140529.ont.all';
%---------------------------------------------------------------------------------------------------

% >>> compute_bhv_ratemaps >>> -------------------------------------------------
AP.compute_bhv_ratemaps =                                                    ...
    struct('get_featureSet',    @fet_HB_pitchB,                              ...
           'sampleRate',        16,                                          ...
           'pfsArgs',           struct('states',           activeState,      ...
                                       'binDims',          [0.1,0.1],        ...
                                       'SmoothingWeights', [1.8,1.8],        ...
                                       'numIter',          1,                ...
                                       'boundaryLimits',   [-2,0.8;-0.8,2],  ...
                                       'halfsample',       false),           ...
           'threshRate',        0.8,                                         ...
           'threshDist',        250                                          ...
           );
% <<< compute_bhv_ratemaps <<< -------------------------------------------------

% >>> compute_bhv_ratemaps_shuffled >>> ----------------------------------------
AP.compute_bhv_ratemaps_shuffled =                                           ...
    struct('get_featureSet',    @fet_HB_pitchB,                              ...
           'sampleRate',        16,                                          ...
           'pfsArgs',           struct('states',           activeState,      ...
                                       'binDims',          [0.1,0.1],        ...
                                       'SmoothingWeights', [1.8,1.8],        ...
                                       'numIter',          1001,             ...
                                       'posShuffle',       true,             ...
                                       'boundaryLimits',   [-2,0.8;-0.8,2],  ...
                                       'halfsample',       false),           ...
           'threshRate',        0.8,                                         ...
           'threshDist',        250                                          ...
           );
% <<< compute_bhv_ratemaps_shuffled <<< ----------------------------------------

% >>> compute_ratemaps >>> -----------------------------------------------------
AP.compute_ratemaps =                                                        ...
    struct('get_featureSet',  @fet_xy,                                       ...
           'sampleRate',      16,                                            ...
           'pfsArgs',         struct('states',           activeState,        ...
                                     'binDims',          [50,50],            ...
                                     'SmoothingWeights', [2.4,2.4],          ...
                                     'numIter',          1,                  ...
                                     'boundaryLimits',   [-500,500;-500,500],...
                                     'halfsample',       false)              ...
);
% <<< compute_ratemaps <<< -----------------------------------------------------

% >>> compute_bhv_ratemap_erpPCA >>> ------------------------------------------
AP.compute_bhv_ratemaps_erpPCA =                                            ...
    struct('numComp',                   5,                                  ...
           'range',                     [250,600]                           ...
);
% <<< compute_bhv_ratemap_erpPCA <<< ------------------------------------------

