function [pfs] = compute_ratemap_ghz_rhm(Trial,units,xyz,spk,pft,rot,hbaCorrection,headCenterCorrection,overwrite)

sampleRate = xyz.sampleRate;
binPhzs = linspace(0,2*pi,6);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
pfs = cell([1,numel(binPhzc)]);
verbose = true;

if verbose,
    disp(['[status]        compute_ego_ratemap: processing trial: ',Trial.filebase]);
end

if isempty(units),
    return;
end;% if

% COMPUTE head basis
hvec = xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(rot),-sin(rot);sin(rot),cos(rot)],...
                 [2,3],...
                 [1,2]);

% GET theta state behaviors, minus rear
thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);

pfTemp = Trial;

pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.units        = units;
pargs.tag          = 'egofield';
pargs.binDims      = [20, 20]; % X Y 
pargs.SmoothingWeights = [3, 3];
pargs.halfsample   = false;
pargs.numIter      = 1;   
pargs.boundaryLimits = [-410,410;-410,410];
pargs.states       = '';
pargs.overwrite    = true;
pargs.autoSaveFlag = false;    
electrode = 0;

% CHECK existence of pfs object
if pargs.SmoothingWeights(1)~=2,
    pargs.tag = ['egofield_SW',num2str(pargs.SmoothingWeights(1))];
else
    pargs.tag = 'egofield';
end
filepath = fullfile(Trial.spath,...
                    [Trial.filebase,'.Pfs.',pargs.tag,'.mat']);
    
if exist(filepath,'file'),
    pfs = load(filepath);
    pfs = pfs.Pfs;
    if overwrite & isa(pfs,'MTAApfs'),
        pfs.purge_savefile();
    else,
        return;
    end;% if overwrite
end;% if exist


for unit = 1:numel(units),
    if unit==1 | electrode ~= spk.map(spk.map(:,1)==units(unit),2), % update phase state
        pargs.spk = copy( spk );
        electrode = 1;
        pargs.states = copy( thetaState );
        pargs.states.label = ['theta'];
        cast( pargs.states, 'TimePeriods' );
        resInd = WithinRanges( pargs.spk.res, pargs.states.data );
        pargs.spk.res = pargs.spk.res( resInd );
        pargs.spk.clu = pargs.spk.clu( resInd );
    end;% if
    
    [mxr,mxp] = pft.maxRate(units(unit));
    pfsCenterHR = MTADfet.encapsulate(Trial,                                            ...
                                      [bsxfun(@plus,                                    ...
                                              multiprod(bsxfun(@minus,                  ...
                                                               mxp,                     ...
                                                               sq(xyz(:,'hcom',[1,2]))),...
                                                         hvec,2,[2,3]),                 ...
                                              headCenterCorrection)],                   ...
                                       sampleRate,                                      ...
                                      'egocentric_placefield',                          ...
                                      'egopfs',                                         ...
                                      'p'                                               ...
                                      );
    pargs.xyzp = pfsCenterHR;
    pargs.units  = units(unit);
    pfsArgs = struct2varargin(pargs);
    pfTemp = MTAApfs(pfTemp,pfsArgs{:});
    if unit==1,
        try,
            pfTemp.purge_savefile();
        end;
        pfTemp.save();        
    end;% if 
end;% for unit
pfTemp.save();
pfs = pfTemp;    
pfTemp = Trial;

