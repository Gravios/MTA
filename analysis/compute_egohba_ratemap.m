function [pfs] = compute_egohba_ratemap(Trial,units,varargin)
% function [pfs] = compute_egohb_ratemap(Trial,units,varargin)
% 
% Compute the mean firing rate of a neuron given spike position as the place field center 
% in the 2d coordinate system of the head, theta phase and head body angle.
%


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct(             'xyz', [],                                                         ...
                              'spk', [],                                                         ...
                              'pft', [],                                                         ...
                'headYawCorrection', Trial.meta.correction.headYaw,                              ...
             'headCenterCorrection', [0,0],                                                      ...
                       'overwrite' , false                                                       ...
);
[xyz, spk, pft, headYawCorrection, headCenterCorrection, overwrite] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

sampleRate = xyz.sampleRate;
binPhzs = linspace(0.5,2*pi-0.5,4);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
binHbas = [-1.2,-0.2,0.2,1.2];
binHbac = (binHbas(1:end-1)+binHbas(2:end))./2;
pfs = cell([numel(binPhzc),numel(binHbac)]);
verbose = true;

if verbose,
    disp(['[status] compute_egohba_ratemap: processing trial: ',Trial.filebase]);
end

if isempty(units), return; end;

if overwrite,
% COMPUTE anglular distance between the head and body
    headBodyAng = [xyz(:,'spine_upper',[1,2])-xyz(:,'bcom',[1,2]),...
                   xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2])];
    headBodyAng = sq(bsxfun(@rdivide,headBodyAng,sqrt(sum(headBodyAng.^2,3))));
    headBodyAng = cart2pol(headBodyAng(:,:,1),headBodyAng(:,:,2));
    headBodyAng = circ_dist(headBodyAng(:,2),headBodyAng(:,1));
    headBodyAng = MTADfet.encapsulate(Trial,...
                                      -(headBodyAng+Trial.meta.correction.headBody),...
                                      sampleRate,...
                                      'hba','hba','h');
% TRANSFORM Local Field Potential -> theta phase
    phz = load_theta_phase(Trial,xyz);
% COMPUTE head basis
    hvec = xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    hvec = multiprod(hvec,...
                     [cos(headYawCorrection),-sin(headYawCorrection);...
                      sin(headYawCorrection),cos(headYawCorrection)],...
                     [2,3],...
                     [1,2]);
% GET theta state behaviors, minus rear
    thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);
end

pfTemp = Trial;

pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.units        = units;
pargs.tag          = 'egofield';
pargs.binDims      = [20, 20];                           % X Y HBA
pargs.SmoothingWeights = [2, 2];                     % X Y HBA
pargs.halfsample   = false;
pargs.numIter      = 1;   
pargs.boundaryLimits = [-410,410;-410,410];
pargs.states       = '';
pargs.overwrite    = true;
pargs.autoSaveFlag = false;    


% Don't judge me
if pargs.SmoothingWeights(1)~=2,
    stag = ['_SW',num2str(pargs.SmoothingWeights(1))];
else
    stag = '';
end

for phzI = 1:numel(binPhzc)
    for hbaI = 1:numel(binHbac)
% CHECK existence of pfs object
    pargs.tag = ['egofield_theta_phase_hba_',num2str(hbaI),'_',num2str(phzI),stag];
    filepath = fullfile(Trial.spath, [Trial.filebase,'.Pfs.',pargs.tag,'.mat']);
    
    if exist(filepath,'file'),
        pfs{phzI,hbaI} = load(filepath).Pfs;
        if overwrite
            pfs{phzI,hbaI}.purge_savefile();
        else
            continue;
        end;% if
    end;% if exist
        
    for unit = 1:numel(units)
        if unit==1
            pargs.spk = copy(spk);
            pargs.states = copy(thetaState);
            pargs.states.label = ['thetaPhz_',num2str(phzI)];
            pargs.states.data(   phz(:,1) < binPhzs(phzI)             ...
                               | phz(:,1) >= binPhzs(phzI+1)          ... 
                               | headBodyAng.data < binHbas(hbaI)       ...
                               | headBodyAng.data >= binHbas(hbaI+1)) = 0;
            cast(pargs.states,'TimePeriods');
            resInd = WithinRanges( pargs.spk.res, pargs.states.data);
            pargs.spk.res = pargs.spk.res(resInd);
            pargs.spk.clu = pargs.spk.clu(resInd);
        end% if
        
        [mxr,mxp] = pft.maxRate(units(unit));
        pfsCenterHR = MTADfet.encapsulate(Trial,                                           ...
                                          bsxfun(                                         ...
                                              @plus,                                       ...
                                              multiprod(bsxfun(@minus,                     ...
                                                          mxp,                             ...
                                                          sq(xyz(:,'hcom',[1,2]))),        ...
                                                     hvec,2,[2,3]),                        ...
                                              headCenterCorrection),                       ...
                                          sampleRate,                                      ...
                                          'egocentric_placefield',                         ...
                                          'egopfs',                                        ...
                                          'p'                                              ...
                                          );
        pargs.xyzp = pfsCenterHR;
        pargs.units  = units(unit);
        pfsArgs = struct2varargin(pargs);
        pfTemp = MTAApfs(pfTemp,pfsArgs{:});
        if unit==1
            try
                pfTemp.purge_savefile();
            end
            pfTemp.save();        
        end% if 
    end% for unit
    pfTemp.save();
    pfs{phzI,hbaI} = pfTemp;    
    pfTemp = Trial;
    end% for hbaI
end% for phzI
