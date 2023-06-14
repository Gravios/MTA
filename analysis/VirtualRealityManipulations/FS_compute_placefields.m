% ROOM frame of reference  (X,Y)
pfsArgs = struct('states',           activeState,                  ... Computational Periods 
                 'tag',              'theta',                      ...
                 'binDims',          [20,20],                      ... Physical size of bins in milimeters
                 'SmoothingWeights', [3.5,3.5],                    ... Gaussian smoother prameters, std deviation in bins 
                 'numIter',          1,                            ... number of bootstraps
                 'boundaryLimits',   [-500,500;-1000,1000],        ... Computational domain
                 'halfsample',       false);                      %... throw out half of the data each iteration ... or don't

pfs = compute_ratemaps(Trial,                                      ... MTATrial
                       [],                                         ... Unit list (e.g. [1, 2, ... , N])
                       @fet_rbo_AC,                                ... Function handle, ratemap space
                       [],                                         ... sampleRate, default 16Hz to speed it up
                       pfsArgs,                                    ... Arguments to the MTAApfs 
                       'overwrite',true);


% ROOM frame of reference  (X,Y)
pfsArgs = struct('states',           'theta&gper&rear',                  ... Computational Periods 
                 'tag',              'rear',                       ...
                 'binDims',          [20,20],                      ... Physical size of bins in milimeters
                 'SmoothingWeights', [3.5,3.5],                    ... Gaussian smoother prameters, std deviation in bins 
                 'numIter',          1,                            ... number of bootstraps
                 'boundaryLimits',   [-500,500;-1000,1000],        ... Computational domain
                 'halfsample',       false);                      %... throw out half of the data each iteration ... or don't

pfr = compute_ratemaps(Trial,                                      ... MTATrial
                       [],                                         ... Unit list (e.g. [1, 2, ... , N])
                       @fet_rbo_AC,                                ... Function handle, ratemap space
                       [],                                         ... sampleRate, default 16Hz to speed it up
                       pfsArgs,                                    ... Arguments to the MTAApfs 
                       'overwrite',true);

                       
% ROOM frame of reference  (X,Y)
pfsArgs = struct('states',           'theta&gper&approach-rear',   ... Computational Periods
                 'tag',              'approach',                   ...
                 'binDims',          [20,20],                      ... Physical size of bins in milimeters
                 'SmoothingWeights', [3.5,3.5],                    ... Gaussian smoother prameters, std deviation in bins 
                 'numIter',          1,                            ... number of bootstraps
                 'boundaryLimits',   [-500,500;-1000,1000],        ... Computational domain
                 'halfsample',       false);                      %... throw out half of the data each iteration ... or don't

pfa = compute_ratemaps(Trial,                                      ... MTATrial
                       [],                                         ... Unit list (e.g. [1, 2, ... , N])
                       @fet_rbo_AC,                                ... Function handle, ratemap space
                       [],                                         ... sampleRate, default 16Hz to speed it up
                       pfsArgs,                                    ... Arguments to the MTAApfs 
                       'overwrite',true);

                       
% ROOM frame of reference  (X,Y)
pfsArgs = struct('states',           'theta&gper&depart-rear',     ... Computational Periods 
                 'tag',              'depart',                     ...
                 'binDims',          [20,20],                      ... Physical size of bins in milimeters
                 'SmoothingWeights', [3.5,3.5],                    ... Gaussian smoother prameters, std deviation in bins 
                 'numIter',          1,                            ... number of bootstraps
                 'boundaryLimits',   [-500,500;-1000,1000],        ... Computational domain
                 'halfsample',       false);                      %... throw out half of the data each iteration ... or don't

pfd = compute_ratemaps(Trial,                                      ... MTATrial
                       [],                                         ... Unit list (e.g. [1, 2, ... , N])
                       @fet_rbo_AC,                                ... Function handle, ratemap space
                       [],                                         ... sampleRate, default 16Hz to speed it up
                       pfsArgs,                                    ... Arguments to the MTAApfs 
                       'overwrite',true);
