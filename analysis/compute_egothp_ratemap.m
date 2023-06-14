function [pfs] = compute_egothp_ratemap(Trial,units,varargin)
% function [pfs] = compute_egothp_ratemap(Trial,units,varargin)
% 
% Compute the mean firing rate of a neuron given spike position as the place field center 
% in the 2d coordinate system of the head and theta phase.
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
% $$$ binPhzs = linspace(0,2*pi,6);
% $$$ binPhzc = mean([binPhzs(1:end-1);binPhzs(2:end)]);
pfs = cell([1,numel(binPhzc)]);
verbose = true;

if verbose,
    disp(['[status] compute_egothp_ratemap: processing trial: ',Trial.filebase]);
end

if isempty(units), pfs = []; return; end;

if overwrite
% TRANSFORM Local Field Potential -> theta phase
    phz = load_theta_phase(Trial,xyz);
% COMPUTE head basis 
    hvec = xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    hvec = multiprod(hvec,...
                     [cos(headYawCorrection),-sin(headYawCorrection); ...
                      sin(headYawCorrection), cos(headYawCorrection)],...
                     [2,3],...
                     [1,2]);
% GET theta state behaviors, minus rear
    thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);
end

pfTemp = Trial;

pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.units        = units;
pargs.tag          = 'egofield';
%pargs.tag          = 'egofieldR';
%pargs.tag          = 'egofieldM';
pargs.binDims      = [20, 20];                           % X Y
pargs.SmoothingWeights = [2, 2];                         % X Y
%pargs.SmoothingWeights = [0.8, 0.8];                         % X Y egofieldR
%pargs.SmoothingWeights = [1.2, 1.2];                         % X Y egofieldM
pargs.halfsample   = false;
pargs.numIter      = 1;   
pargs.boundaryLimits = [-500,500;-500,500];
pargs.states       = '';
pargs.overwrite    = true;
pargs.autoSaveFlag = false;    
electrode = 0;

% Don't judge me
if pargs.SmoothingWeights(1)~=2,
    stag = ['_SW',num2str(pargs.SmoothingWeights(1))];
else
    stag = '';
end

for phase = 1:numel(binPhzc)
% CHECK existence of pfs object
    pargs.tag = ['egofield_theta_phase_',num2str(phase),stag];

    filepath = fullfile(Trial.spath, [Trial.filebase,'.Pfs.',pargs.tag,'.mat']);
    
    if exist(filepath,'file'),
        pfs{phase} = load(filepath).Pfs;
        if overwrite,
            pfs{phase}.purge_savefile();
        else,
            continue;
        end;% if
    end;% if exist
        
    for unit = 1:numel(units),
        if unit==1 | electrode~=spk.map(spk.map(:,1)==units(unit),2), % update phase state            
            pargs.spk = copy(spk);
            pargs.states = copy(thetaState);
            pargs.states.label = ['thetaPhz_',num2str(phase)];
            pargs.states.data((phz(:,1) < binPhzs(phase) )    ...
                              | (phz(:,1) >= binPhzs(phase+1)) ) = 0;
            cast(pargs.states,'TimePeriods');
            resInd = WithinRanges(pargs.spk.res,pargs.states.data);
            pargs.spk.res = pargs.spk.res(resInd);
            pargs.spk.clu = pargs.spk.clu(resInd);
        end;% if
        
        [mxr,mxp] = pft.maxRate(units(unit));
        pfsCenterHR = MTADfet.encapsulate(Trial,                                           ...
                                          [bsxfun(                                         ...
                                              @plus,                                       ...
                                              multiprod(bsxfun(@minus,                     ...
                                                               mxp,                        ...
                                                               sq(xyz(:,'hcom',[1,2]))),   ...
                                                         hvec,2,[2,3]),                    ...
                                              headCenterCorrection)],                      ...
                                          sampleRate,                                      ...
                                          'egocentric_placefield',                         ...
                                          'egopfs',                                        ...
                                          'p'                                              ...
                                          );
        pargs.xyzp = pfsCenterHR;
        pargs.units  = units(unit);
        pfsArgs = struct2varargin(pargs);
        pfTemp = MTAApfs(pfTemp,pfsArgs{:});
        if unit==1,
            try
                pfTemp.purge_savefile();
            end
            pfTemp.save();        
        end% if 
    end;% for unit
    pfTemp.save();
    pfs{phase} = pfTemp;    
    pfTemp = Trial;
end;% for phase
