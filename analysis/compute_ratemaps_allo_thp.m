function [pfs] = compute_ratemaps_allo_thp(Trial,varargin);
% function [pfs,auxData] = compute_bhv_ratemaps(Trial,varargin);
% 
% Compute spatially restricted behavior ratemaps
%
% CODE STRUCTURE
% 
%  1. ATTEMPT to load existing MTAApfs object from filesytem
%  2. IF all units exist return object
%  2. ELSE Compute
%


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct(           'units', [],                                                         ...
                              'xyz', [],                                                         ...
                              'spk', [],                                                         ...
                'headYawCorrection', Trial.meta.correction.headYaw,                              ...
             'headCenterCorrection', [0,0],                                                      ...
                       'overwrite' , false                                                       ...
);
[units, xyz, spk, headYawCorrection, headCenterCorrection, overwrite] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

sampleRate = xyz.sampleRate;
binPhzs = linspace(0.5,2*pi-0.5,4);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
pfs = cell([1,numel(binPhzc)]);


verbose = true;

if verbose,
    disp(['[status] compute_ratemaps_allothp: processing trial: ',Trial.filebase]);
end

if isempty(units), pfs = []; return; end;



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


if overwrite
% TRANSFORM Local Field Potential -> theta phase
    phz = load_theta_phase(Trial,xyz);
% GET theta state behaviors, minus rear
    thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);
end



pfTemp = Trial;

pargs.units        = units;
pargs.tag          = '';
%pargs.tag          = 'allofieldR';
%$pargs.tag          = 'allofieldM';
pargs.binDims      = [20, 20];                           % X Y
%pargs.SmoothingWeights = [0.8, 0.8];                         % X Y allofieldR
%pargs.SmoothingWeights = [1.2, 1.2];                         % X Y allofieldM
pargs.SmoothingWeights = [2, 2];                         % X Y 
pargs.halfsample   = false;
pargs.numIter      = 1;   
pargs.boundaryLimits = [-500,500;-500,500];
pargs.states       = '';
pargs.overwrite    = true;
pargs.autoSaveFlag = false;    
pargs.xyzp = fet_xy(Trial,sampleRate);
pargs.spk = copy(spk);

electrode = 0;



% Don't judge me
if pargs.SmoothingWeights(1)~=2,
    stag = [pargs.tag,'_SW',num2str(pargs.SmoothingWeights(1))];
else
    stag = pargs.tag;
end


for phase = 1:numel(binPhzc)
% CHECK existence of pfs object
    pargs.tag = ['allofield_theta_phase_',num2str(phase),stag];

    filepath = fullfile(Trial.spath, [Trial.filebase,'.Pfs.',pargs.tag,'.mat']);

    if exist(filepath,'file'),
        pfs{phase} = load(filepath).Pfs;
        if overwrite,
            pfs{phase}.purge_savefile();
        else,
            continue;
        end;% if
    end;% if exists file

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
        end;% if unit==1

        pargs.units  = units(unit);
        pfsArgs = struct2varargin(pargs);
        pfTemp = MTAApfs(pfTemp,pfsArgs{:});
        if unit==1,
            try,
                pfTemp.purge_savefile();
            end;
            pfTemp.save();        
        end;% if unit==1
        
    end;% for unit
    pfTemp.save();
    pfs{phase} = pfTemp;    
    pfTemp = Trial;
end;% for phase


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