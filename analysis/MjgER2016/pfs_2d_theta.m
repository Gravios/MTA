function pfs = pfs_2d_theta(Trial,varargin)
%function pft = pfs_2d_theta(Trial,varargin)
%
% Load/compute specified unit placefields of a Trial (MTATrial) object.
% placefields are computed with the horizontal (xy) plane using periods during the hippocampal 
% theta state.
%
%  varargin:
%    units:      numeric, subset of unit cluster ids.
% 
%    states:     string or MTADepoch
%
%    overwrite:  logical, flag 1: recompute place fields
%                              0: load place fields from file
%
%    pfsArgsOverride: struct, list of vars and arguments to be overridden when computing placefield
%
%        units      -  Array,        list of units
%        states     -  String,       computation periods
%        overwrite  -  Logical,      overwrite rate maps
%        tag        -  String,       optional string for filename 
%        binDims    -  Array,        dimensions of ratemap bins (mm)
%        SmoothingWeights - Array,   std of gaussian kernel (bins)
%        type       -  String,       dimension discriptive characters (e.g. xy, xyz)
%        spkShuffle -  Logical,      FLAG, shuffle spiketimes  randomly
%        posShuffle -  Logical,      FLAG, shuflle position (blocks) randomly
%        numIter    -  Integer,      number of bootstrap iterations
%        xyzp       -  MTAData,      Alternative data object for ratemap space
%        boundaryLimits - Matrix,    Computational boundaries (default->Trial.maze.boundaries)
%        bootstrap  -  Logical,      FLAG, Compute bootstrap ratemaps
%        halfsample -  Logical,      FLAG, Compute each bootstrap with half the data
%        shuffleBlockSize - Double,  Data block size in seconds
%        trackingMarker   - String,  Marker name within MTADxyz object 
%                                    (ignored if xyzp not empty)
%        autoSaveFlag - Logical,     FLAG save between unit computations
%        spkMode      - Logical,     FLAG remove burst (ignored if spk is provided)
%        spk          - MTASpk,      contains unit res and clu
%        compute_pfs  - function,    function handle to computational core
%        loadMethod   - RESERVED FOR FUTURE USE
%

% DEFARGS ---------------------------------------------------------------------------------------
defargs = struct('units',          Trial.spk.map(:,1)',                                       ...
                 'states',         'theta-sit-groom',                                         ...
                 'overwrite',      0,                                                         ...
                 'pfsArgsOverride',[],                                                        ...
                 'purge',          false                                                      ...
);
[units,states,overwrite,pfsArgsOverride,purge] = DefaultArgs(varargin,defargs,'--struct');
% END DEFARGS -----------------------------------------------------------------------------------



% DEFVARS ----------------------------------------------------------------------
Trial= MTATrial.validate(Trial);    


% load theta periods
if isempty(Trial.stc.gsi('t')),Trial = labelTheta(Trial);end

% END DEFVARS ----------------------------------------------------------------------



% MAIN -------------------------------------------------------------------------    
%% compute 2d place fields for the theta state

try, 
    Trial.stc{states};
catch err
    disp(err)
    Trial.stc{'walk&gper&theta'};
    states = 'walk&gper&theta';
end

pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units = units;
pargs.states = states;
pargs.overwrite = overwrite;
pargs.binDims = [20,20];
pargs.SmoothingWeights = [3,3];
pargs.numIter = 1;
pargs.halfsample = false;
if ~isempty(pfsArgsOverride) && isstruct(pfsArgsOverride),
    for f = fieldnames(pfsArgsOverride)'
        pargs.(f{1}) = pfsArgsOverride.(f{1});
    end
end
pfsArgs = struct2varargin(pargs);
disp(['[INFO] MTA:analysis:placefields:pfs_2d_theta:processing:' Trial.filebase]);
pfs = MTAApfs(Trial,'purge',purge,pfsArgs{:});


% END MAIN -------------------------------------------------------------------------    