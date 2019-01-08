function [drzState] = req20180123_pfd_shuffled(Trial,xyzp,tag,units,drz,ddz,tper,boundaryLimits,binDims,SmoothingWeights);


Trial = MTATrial.validate(Trial);

pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.tag              = [tag,'-shuffled'];
pargs.units            = units;
pargs.numIter          = 1001;
pargs.posShuffle       = true;
pargs.halfsample       = false;
pargs.overwrite        = true;
pargs.boundaryLimits   = boundaryLimits;
pargs.binDims          = binDims;
pargs.SmoothingWeights = SmoothingWeights;
pargs.xyzp             = xyzp.copy();
pargs.bootstrap        = false;
pargs.autoSaveFlag     = false;

drzState = {};
u = 1;        
dper = MTADepoch([],[],ThreshCross(abs(drz(:,u))<0.8 & ddz(:,u)<250,0.5,1),...% SELECT periods where drz
                 xyzp.sampleRate,xyzp.sync.copy(),xyzp.origin,'TimePeriods','sts',[],'tdrz','d');
drzState{u} = dper&tper;
pargs.units  = units(u);
pargs.states = drzState{u};
pfsArgs = struct2varargin(pargs);
pfTemp = MTAApfs(Trial,pfsArgs{:});
pfTemp.save();
for u = 2:numel(units);
    dper = MTADepoch([],[],ThreshCross(abs(drz(:,u))<0.8 & ddz(:,u)<250,0.5,1),...% SELECT periods where drz
                     xyzp.sampleRate,xyzp.sync.copy(),xyzp.origin,'TimePeriods','sts',[],'tdrz','d');
    drzState{u} = dper&tper;
    pargs.units  = units(u);
    pargs.states = drzState{u};
    pfsArgs = struct2varargin(pargs);
    try
        pfTemp = MTAApfs(pfTemp,pfsArgs{:});            
    catch err,
        disp(err);
        pfTemp = MTAApfs(pfTemp,pfsArgs{:});
    end
end
pfTemp.save();        

