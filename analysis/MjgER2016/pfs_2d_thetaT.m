function pfs = pfs_2d_thetaT(Trial,varargin)
%function pft = pfs_2d_theta(Trial,varargin)
%
% Load/compute specified unit placefields of a Trial (MTATrial) object.
% placefields are computed with the horizontal (xy) plane using periods during the hippocampal 
% theta state.
%
%  varargin:
%    units:      numeric, subset of unit cluster ids.
% 
%    reportFig:  logical, flag 1: save figures as pngs of each place field in a predefined location
%                              0: return placefield object without creating figures
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
defargs = struct('units',              Trial.spk.map(:,1)',                                   ...
                 'thetaRef',           1,                                                     ...
                 'phzCorrection',      0,                                                     ...
                 'reportFig',          0,                                                     ...
                 'overwrite',          0,                                                     ...
                 'pfsArgsOverride',    []                                                     ...
);
[units,thetaRef,phzCorrection,reportFig,overwrite,pfsArgsOverride] = DefaultArgs(varargin,defargs,'--struct');
% END DEFARGS -----------------------------------------------------------------------------------



% DEFVARS ----------------------------------------------------------------------
Trial= MTATrial.validate(Trial);    
sampleRate = 250;%Hz

% load theta periods
if isempty(Trial.stc.gsi('t')),Trial = labelTheta(Trial);end

%% Setup figure paths
if reportFig
    OwnDir = '/storage/gravio/nextcloud/MjgER2016/';
    FigDir = ['pfs_2d_theta_',Trial.filebase];
    mkdir(fullfile(OwnDir,FigDir));
end
% END DEFVARS ----------------------------------------------------------------------


% MAIN -------------------------------------------------------------------------    
%% compute 2d place fields for the theta state

try, 
    state = 'theta-sit-groom-rear';
    Trial.stc{state};
catch
    state = 'theta';
end

tper = Trial.stc{state};


if ~overwrite,
    pargs = get_default_args('MjgER2016','MTAApfs','struct');
    pargs.units = units;
    pargs.tag = 'thetaT';
    if ~isempty(pfsArgsOverride) && isstruct(pfsArgsOverride),
        for f = fieldnames(pfsArgsOverride)'
            pargs.(f{1}) = pfsArgsOverride.(f{1});
        end
    end
    pfsArgs = struct2varargin(pargs);    
    disp(['MTA:analysis:placefields:pfs_2d_theta:processing:' Trial.filebase]);
    pfs = MTAApfs(Trial,pfsArgs{:});
else
% COMPUTE theta phase
Trial.lfp.filename = [Trial.name,'.lfp'];% Seriously just fix that one session
lfp = Trial.load('lfp',thetaRef);
phz = lfp.phase([5,13]);    
phz.data = unwrap(phz.data);
phz.resample(sampleRate);    
%phz.data = mod(phz.data+pi,2*pi)-pi;
phz.data = mod(phz.data+pi,2*pi)-pi + phzCorrection; 
phz.data(phz.data>pi) = phz.data(phz.data>pi)-2*pi;


tper.cast('TimeSeries');
tper.resample(phz);
tpT = phz.data>pi*7/8 | phz.data<-pi*7/8;

tperT = tper&tpT;
tperT.cast('TimePeriods');
tperT.label = 'thetaT';

spk = Trial.spk.create(Trial,sampleRate,tperT,units,'');

pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units = units;
pargs.tag = 'thetaT';
pargs.binDims = [20,20];
pargs.SmoothingWeights = [3,3];
pargs.states = tperT;
pargs.overwrite = overwrite;
pargs.xyzp = fet_xy(Trial,sampleRate);
pargs.spk = spk;
pargs.numIter = 1;
pargs.halfsample = false;
if ~isempty(pfsArgsOverride) && isstruct(pfsArgsOverride),
    for f = fieldnames(pfsArgsOverride)'
        pargs.(f{1}) = pfsArgsOverride.(f{1});
    end
end
pfsArgs = struct2varargin(pargs);
disp(['MTA:analysis:placefields:pfs_2d_theta:processing:' Trial.filebase]);
pfs = MTAApfs(Trial,pfsArgs{:});
end

% END MAIN -------------------------------------------------------------------------    