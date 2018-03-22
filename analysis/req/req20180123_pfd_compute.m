function [drzState] = req20180123_pfd_compute(Trial,xyzp,tag,units,drz,tper,boundaryLimits,binDims,SmoothingWeights);


pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.tag              = tag;
pargs.units            = units;
pargs.numIter          = 1001;
pargs.halfsample       = true;
pargs.overwrite        = true;
pargs.boundaryLimits   = boundaryLimits;
pargs.binDims          = binDims;
pargs.SmoothingWeights = SmoothingWeights;
pargs.xyzp             = xyzp;
pargs.autoSaveFlag     = false;

drzState = {};
u = 1;        
dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                 xyzp.sampleRate,xyzp.sync.copy(),xyzp.origin,'TimePeriods','sts',[],'tdrz','d');
drzState{u} = dper&tper;
pargs.units  = units(u);
pargs.states = drzState{u};
pfsArgs = struct2varargin(pargs);
pfTemp = MTAApfs(Trial,pfsArgs{:});
pfTemp.save();
for u = 1:numel(units);
    dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                     xyzp.sampleRate,xyzp.sync.copy(),xyzp.origin,'TimePeriods','sts',[],'tdrz','d');
    drzState{u} = dper&tper;
    pargs.units  = units(u);
    pargs.states = drzState{u};
    pfsArgs = struct2varargin(pargs);
    pfTemp = MTAApfs(pfTemp,pfsArgs{:});            
end
pfTemp.save();        
