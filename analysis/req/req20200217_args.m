
% req20200217_args


clearvars('-GLOBAL','AP');

global AP

% GENERAL ------------------------------------------------------------------------------------------
sessionListName = 'MjgER2016';
activeState     = 'theta-groom-sit-turn-rear';
%activeState     = 'theta-groom-sit';

AP.fet_HB_pitchB.referenceTrial = 'Ed05-20140529.ont.all';
AP.fet_hbp_hba.referenceTrial = 'Ed05-20140529.ont.all';
AP.fet_hbp_hba.offset = [0,-0.2];
%---------------------------------------------------------------------------------------------------



% compute_bhv_ratemaps -----------------------------------------------------------------------------
AP.compute_bhv_ratemaps =                                                                        ...
    struct('get_featureSet',            @fet_hbp_hbd,                                            ...
           'sampleRate',                16,                                                      ...
           'pfsArgs',                   struct('states',           activeState,                  ...
                                               'binDims',          [0.1,5],                      ...
                                               'SmoothingWeights', [1.8,1.8],                    ...
                                               'numIter',          1,                            ...
                                               'boundaryLimits',   [-2,0.8;110,160],             ...
                                               'halfsample',       false),                       ...
           'threshRate',                0.8,                                                     ...
           'threshDist',                250                                                      ...
           );
% $$$ AP.compute_bhv_ratemaps =                                                                        ...
% $$$     struct('get_featureSet',            @fet_hbp_hba,                                            ...
% $$$            'sampleRate',                16,                                                      ...
% $$$            'pfsArgs',                   struct('states',           activeState,                  ...
% $$$                                                'binDims',          [0.1,0.1],                    ...
% $$$                                                'SmoothingWeights', [1.8,1.8],                    ...
% $$$                                                'numIter',          1,                            ...
% $$$                                                'boundaryLimits',   [-2,0.8;-1.2,1.2],            ...
% $$$                                                'halfsample',       false),                       ...
% $$$            'threshRate',                0.8,                                                     ...
% $$$            'threshDist',                250                                                      ...
% $$$            );
%---------------------------------------------------------------------------------------------------
