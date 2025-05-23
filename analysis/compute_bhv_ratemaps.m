function [pfs] = compute_bhv_ratemaps(Trial,varargin);
% function [pfs,auxData] = compute_bhv_ratemaps(Trial,varargin);
% 
% Compute spatially restricted behavior ratemaps
% VARARGIN :
%    units          - Array(Numeric): {[]            }
%    get_featureSet - FuncHandle:     {@fet_HB_pitchB}
%    sampleRate     - Numeric:        {16            }
%    pft            - MTAApfs:        {[]            }
%    pfsArgs        - Struct:         {struct('states',           'theta-groom-sit', 
%                                             'binDims',          [0.2,0.2],         
%                                             'SmoothingWeights', [2,2],             
%                                             'numIter',          1,                 
%                                             'boundaryLimits',   [-2,0.8;-0.8,2],   
%                                             'halfsample',       false))}
%    threshRate     - Numeric:        {[]            }
%    threshDist     - Numeric:        {[]            }
%    overwrite      - Logicial:       {false         }



% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',                         [],                                            ...
                 'get_featureSet',                @fet_HB_pitchB,                                ...
                 'sampleRate',                    16,                                            ...
                 'pft',                           [],                                            ...
                 'pfsArgs',                       struct('states',           'theta-groom-sit',  ...
                                                         'binDims',          [0.2,0.2],          ...
                                                         'SmoothingWeights', [2,2],              ...
                                                         'numIter',          1,                  ...
                                                         'boundaryLimits',   [-2,0.8;-0.8,2],    ...
                                                         'halfsample',       false),             ...
                 'threshRate',                    [],                                            ...
                 'threshDist',                    [],                                            ...
                 'overwrite',                     false,                                         ...
                 'purge',                         false                                          ...
);
[units,get_featureSet,sampleRate,pft,pfsArgs,threshRate,threshDist,overwrite,purge] =            ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% TAGS ---------------------------------------------------------------------------------------------
pfsTag = DataHash({func2str(get_featureSet),sampleRate,threshRate,threshDist,...
                   pfsArgs.states,pfsArgs.binDims,pfsArgs.numIter,pfsArgs.boundaryLimits});
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------


% ATTEMPT to load existing MTAApfs object from filesytem
pfs = [];
if purge
% LOAD existing MTAApfs object from filesytem    
    pfs = MTAApfs(Trial,'tag',pfsTag);
    pfs.purge_savefile();
end


if ~overwrite && MTAApfs.exist(Trial,pfsTag)
% LOAD existing MTAApfs object from filesytem    
    pfs = MTAApfs(Trial,'tag',pfsTag);
    if all(ismember(units,pfs.data.clu)),
% RETURN MTAApfs object        
        return;
    end
end


states = pfsArgs.states;
if ischar(pft) && strcmp(pft,'none'),
    pft = [];
elseif isempty(pft),
% LOAD spatial placefields
    pft = pfs_2d_theta(Trial,units);
end

% LOAD Feature space
fet  = feval(get_featureSet,Trial,sampleRate);
% LOAD distance metrics
xyz  = resample(preproc_xyz(Trial,'trb'),sampleRate);
% TODO add logical statement to ...
%      IF threshRate >= 1
%       ... ( skip load ) assign zeros to drz 
%      IF threshDist >= 
%       ... ( skip load ) assign zeros to ddz
if threshRate < 1,
    drz = compute_drz(Trial,units,pft,'sampleRate',sampleRate);
else
    drz = zeros([size(xyz,1),numel(units)]);
end
if threshDist < sqrt(sum(diff(Trial.maze.boundaries(1:2,:),1,2).^2)),
    ddz = compute_ddz(Trial,units,pft,'sampleRate',sampleRate);
else
    ddz = zeros([size(xyz,1),numel(units)]);
end

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
if sum(double(dper.data))==0,
    dper.data(:) = true;
end

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
    if sum(double(dper.data))==0,
        dper.data(:) = true;
    end
    
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
% $$$     auxData.threshRate  = threshRate;
% $$$     auxData.threshDist  = threshDist;    
% $$$ end


% END MAIN -----------------------------------------------------------------------------------------