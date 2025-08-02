function [pfs] = compute_xyhb_ratemaps(Trial,varargin);
% function [pfs,auxData] = compute_bhv_ratemaps(Trial,varargin);
% 
% Compute spatially restricted behavior ratemaps
%
% CODE STRUCTURE
% 
%  1. ATTEMPT to load existing MTAApfs object from filesytem
%  2. IF all units exist return object
%  2. ELSE Compute
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',                         [],                                            ...
                 'get_featureSet',                @fet_xyhb,                                     ...
                 'sampleRate',                    16,                                            ...
                 'pfsArgs',                       struct('states',           'theta-groom-sit',  ...
                                                         'binDims',          [0.2,0.2],          ...
                                                         'SmoothingWeights', [2,2],              ...
                                                         'numIter',          1,                  ...
                                                         'boundaryLimits',   [-2,0.8;-0.8,2],    ...
                                                         'halfsample',       false),             ...
                 'overwrite',                     false                                          ...
);
[units,get_featureSet,sampleRate,pfsArgs,overwrite] =                  ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% TAGS ---------------------------------------------------------------------------------------------
pfsTag = DataHash({get_featureSet,sampleRate,...
                   pfsArgs.states,pfsArgs.binDims,pfsArgs.numIter,pfsArgs.boundaryLimits});
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------


% ATTEMPT to load existing MTAApfs object from filesytem
pfs = [];
if ~overwrite && MTAApfs.exist(Trial,pfsTag)
% LOAD existing MTAApfs object from filesytem    
    pfs = MTAApfs(Trial,'tag',pfsTag);
    if all(ismember(units,pfs.data.clu)),
% RETURN MTAApfs object        
        return;
    end
end

% LOAD Feature space
fet  = get_featureSet(Trial,sampleRate);

% ASSIGN name and units to MTAApfs computational arguments
pfsArgs.tag   = pfsTag;
pfsArgs.units = units;

% MAYBE overwrite all values
pfsArgs.overwrite = overwrite;
pfsArgs.units  = units;
pfsArgs.autoSaveFlag = true;
pfsArgs.xyzp = fet;
pfsArgsArray = struct2varargin(pfsArgs);
pfs = MTAApfs(Trial,pfsArgsArray{:});

% $$$ if nargout==1,
% $$$     return;
% $$$ else
% $$$     auxData.periods     = {drzState};
% $$$     auxData.sampleRate  = sampleRate;
% $$$     auxData.get_featureSet = fun2str(get_featureSet);
% $$$     auxData.threshRate  = threshRate;
% $$$     auxData.threshDist  = threshDist;    
% $$$ end


% END MAIN -----------------------------------------------------------------------------------------