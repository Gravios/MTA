function [rmapA] = compute_patch_ratemap(Trial,varargin);
% function [pfs,auxData] = compute_bhv_ratemaps(Trial,varargin);
% 
% Compute spatially restricted behavior ratemaps
% VARARGIN :
%    units          - Array(Numeric): {[]            }
%    get_featureSet - FuncHandle:     {@fet_HB_pitchB}
%    sampleRate     - Numeric:        {16            }
%    patchCenter    - Numeric:        {[x,y,...      }
%    marker         - String:         'hcom'
%    pfsArgs        - Struct:         {struct('states',           'theta-groom-sit', 
%                                             'binDims',          [0.2,0.2],         
%                                             'SmoothingWeights', [2,2],             
%                                             'numIter',          1,                 
%                                             'boundaryLimits',   [-2,0.8;-0.8,2],   
%                                             'halfsample',       false))}
%    threshRate     - Numeric:        {[]            }
%    threshDist     - Numeric:        {[]            }
%    xyz            - MTADxyz:        {[]            }
%    fet            - MTADxyz:        {[]            }
%    tag            - String:         {''            }
%    overwrite      - Logicial:       {false         }



% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',                         [],                                            ...
                 'get_featureSet',                @fet_HB_pitchB,                                ...
                 'sampleRate',                    16,                                            ...
                 'patchCenter',                   [],                                            ...
                 'marker',                        'hcom',                                        ...
                 'pfsArgs',                       struct('states',           'theta-groom-sit',  ...
                                                         'binDims',          [0.2,0.2],          ...
                                                         'SmoothingWeights', [2,2],              ...
                                                         'numIter',          1,                  ...
                                                         'boundaryLimits',   [-2,0.8;-0.8,2],    ...
                                                         'mask',             [],                 ...
                                                         'halfsample',       false),             ...
                 'threshRate',                    [],                                            ...
                 'threshDist',                    [],                                            ...
                 'xyz',                           [],                                            ...
                 'fet',                           [],                                            ...
                 'spk',                           [],                                            ...
                 'tag',                           '',                                            ...
                 'overwrite',                     false                                          ...
);
[units,get_featureSet,sampleRate,patchCenter,marker,pfsArgs,threshRate,threshDist,xyz,fet,...
 spk,tag,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% TAGS ---------------------------------------------------------------------------------------------
pfsTag = DataHash({func2str(get_featureSet),sampleRate,threshRate,threshDist,...
                   pfsArgs.states,pfsArgs.binDims,pfsArgs.numIter,pfsArgs.boundaryLimits});
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------


% ATTEMPT to load existing MTAApfs object from filesytem
% $$$ pfs = [];
% $$$ if ~overwrite && MTAApfs.exist(Trial,pfsTag)
% $$$ % LOAD existing MTAApfs object from filesytem    
% $$$     pfs = MTAApfs(Trial,'tag',pfsTag);
% $$$     if all(ismember(units,pfs.data.clu)),
% $$$ % RETURN MTAApfs object        
% $$$         return;
% $$$     end
% $$$ end

    
% LOAD Feature space
if isempty(fet)
    fet = feval(get_featureSet,Trial,sampleRate);
end
% LOAD distance metrics
if isempty(xyz)
    xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
end

if threshDist < sqrt(sum(diff(Trial.maze.boundaries(1:2,:),1,2).^2)),
    ddz = sqrt(sum((repmat(xyz(:,marker,[1,2]),[1,2,1])-repmat(patchCenter,[size(xyz,1),1,1])).^2,3));
else
    ddz = zeros([size(xyz,1),numel(units)]);
end


% LOAD theta periods
tper = resample(cast([Trial.stc{pfsArgs.states}],'TimeSeries'),xyz);
% LOAD Spike clusters
if isempty(spk),
    spk  = create(copy(Trial.spk),Trial,sampleRate,pfsArgs.states,units,'');
end


% RESTRICT periods by drz and ddz

states{1} = MTADepoch([],[], ddz(:,1)<threshDist & tper.data,sampleRate,fet.sync.copy(),fet.origin,'TimeSeries','sts',[],'tdrz','d');

stateDurA = sum(states{1}.data);

pooledState = double(states{1}.data);
pooledState = find(pooledState);

res = spk(units);

sstposA = fet(logical(states{1}.data),:);

spkindA = ismember(res,find(states{1}.data==1));
if sum(spkindA)~=0, 
    spkposA = fet(res(spkindA),:);
else
    spkposA = [];
end

[rmapA, bins] = PlotPF(Trial,spkposA,sstposA,pfsArgs.binDims,pfsArgs.SmoothingWeights,'hb',pfsArgs.boundaryLimits,sampleRate);





% END MAIN -----------------------------------------------------------------------------------------