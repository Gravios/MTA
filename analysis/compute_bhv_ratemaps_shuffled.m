function [pfs] = compute_bhv_ratemaps_shuffled(Trial,varargin);
% function [pfs,auxData] = compute_bhv_ratemaps_shuffled(Trial,varargin);
% 
% Compute spatially restricted behavior ratemaps with shuffled spike positions
%
% CODE STRUCTURE
% 
%  1. ATTEMPT to load existing MTAApfs object from filesytem
%  2. IF all units exist return object
%  2. ELSE Compute missing 
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',                         [],                                            ...
                 'get_featureSet',                @fet_HB_pitchB,                                ...
                 'sampleRate',                    16,                                            ...
                 'pft',                           [],                                            ...
                 'pfsArgs',                       struct('states',           'theta-groom-sit',  ...
                                                         'binDims',          [0.2,0.2],          ...
                                                         'SmoothingWeights', [2,2],              ...
                                                         'numIter',          1002,               ...
                                                         'posShuffle',       true,               ...
                                                         'boundaryLimits',   [-2,0.8;-0.8,2],    ...
                                                         'halfsample',       false),             ...
                 'threshRate',                    [],                                            ...
                 'threshDist',                    [],                                            ...
                 'overwrite',                     false,                                         ...
                 'purge',                         false                                          ...
);
[units, get_featureSet, sampleRate, pft, pfsArgs, threshRate, threshDist, overwrite, purge] =    ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

if isempty(units)
    pfs = [];
    return;
end


% TAGS ---------------------------------------------------------------------------------------------
pfsTag = DataHash({func2str(get_featureSet),sampleRate,threshRate,threshDist,...
                   pfsArgs.states,pfsArgs.binDims,pfsArgs.numIter,pfsArgs.boundaryLimits,...
                   pfsArgs.posShuffle,pfsArgs.halfsample});
%---------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

if purge
    MTAApfs(Trial,'tag',pfsTag).purge_savefile();
end

    
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


states = pfsArgs.states;
if isempty(pft),
% LOAD spatial placefields
    pft = pfs_2d_theta(Trial,units);
end

% LOAD Feature space
fet  = get_featureSet(Trial,sampleRate);
% LOAD distance metrics
xyz  = resample(preproc_xyz(Trial,'trb'),sampleRate);
drz  = compute_drz(Trial,units,pft,'feature',xyz);
ddz  = compute_ddz(Trial,units,pft,'feature',xyz);
% LOAD theta periods
tper = resample(cast([Trial.stc{states}],'TimeSeries'),xyz);
% LOAD Spike clusters
spk  = create(copy(Trial.spk),Trial,sampleRate,states,units,'deburst');


% ASSIGN name and units to MTAApfs computational arguments
pfsArgs.tag   = pfsTag;
pfsArgs.units = units;
pfsArgs.spk   = spk;

% MAYBE overwrite all values
pfsArgs.overwrite = overwrite;        
    
drzState = {};
u = 1;        
% RESTRICT periods by drz and ddz
dper = MTADepoch([],[],abs(drz(:,u))<threshRate & ddz(:,u)<threshDist,                           ...
                 fet.sampleRate,fet.sync.copy(),fet.origin,'TimeSeries','sts',[],'tdrz','d');
drzState{u} = cast(dper&tper,'TimePeriods');

pfsArgs.units  = units(u);
pfsArgs.states = drzState{u};
pfsArgs.autoSaveFlag = false;
pfsArgs.xyzp = fet;

pfsArgsArray = struct2varargin(pfsArgs);
pfs = Trial;
pfs = MTAApfs(pfs,pfsArgsArray{:});
pfs.save();
for u = 1:numel(units);
% RESTRICT periods by drz and ddz
    dper = MTADepoch([],[],abs(drz(:,u))<threshRate & ddz(:,u)<threshDist,...
                     fet.sampleRate,fet.sync.copy(),fet.origin,'TimeSeries','sts',[],'tdrz','d');
    drzState{u} = cast(dper&tper,'TimePeriods');
    pfsArgs.units  = units(u);
    pfsArgs.states = drzState{u};
    pfsArgsArray = struct2varargin(pfsArgs);
    pfs = MTAApfs(pfs,pfsArgsArray{:});
end
pfs.save();

% $$$ if nargout==1,
% $$$     return;
% $$$ else
% $$$     auxData.periods     = {drzState};
% $$$     auxData.sampleRate  = sampleRate;
% $$$     auxData.get_featureSet = fun2str(get_featureSet);
% $$$     auxData.featureInds = featureInds;
% $$$     auxData.threshRate  = threshRate;
% $$$     auxData.threshDist  = threshDist;    
% $$$ end


% END MAIN -----------------------------------------------------------------------------------------