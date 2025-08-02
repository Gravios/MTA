function pfs = pfd_2d_pp(Trial,varargin)

overwrite = true;
pitchReferenceTrial = 'Ed05-20140529.ont.all';

units = select_placefields(Trial);
xyz = preproc_xyz(Trial,'trb');
pft = pfs_2d_theta(Trial,'overwrite',false);
pch = fet_HB_pitchB(Trial,[],[],[],pitchReferenceTrial);
drz = compute_drz(Trial,units,pft);%,pfstats);

tper = [Trial.stc{'theta-groom-sit'}];
tper.resample(xyz);

drzState = {};
for u = 1:numel(units);
    dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...
                     xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
    drzState{u} = dper&tper;    
end


% CE DRZxPITCHxHEIGHT
pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units      = units;
pargs.numIter = 1001;
pargs.halfsample = true;
pargs.tag            = 'DRZxHBPITCHxBPITCH';
pargs.boundaryLimits = [-pi/2,pi/2;-pi/2,pi/2];
pargs.binDims        = [0.05,0.05];
pargs.SmoothingWeights = [3,3];
if overwrite,
    pargs.overwrite = true;    
    pargs.xyzp = MTADxyz('data',pch.data,'sampleRate',xyz.sampleRate);
    for u = 1:numel(units);
        pargs.units  = units(u);        
        pargs.states = drzState{u};
        pfsArgs = struct2varargin(pargs);
        MTAApfs(Trial,pfsArgs{:});
    end
end
pargs.states    = 'tdrz';
pargs.units     = units;
pargs.overwrite = false;
pfsArgs = struct2varargin(pargs);
pfs = MTAApfs(Trial,pfsArgs{:});
