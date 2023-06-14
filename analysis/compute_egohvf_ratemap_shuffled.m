function [pfs,bins] = compute_egohvf_ratemap_shuffled(Trial,units,varargin)
% function [pfs] = compute_egohb_ratemap(Trial,units,varargin)
% 
% Compute the mean firing rate of a neuron given spike position as the place field center 
% in the 2d coordinate system of the head, theta phase and head forward head velocity.
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

% SET theta phase bins
phzBinEdgs = linspace(0.5,2*pi-0.5,4);
phzBinCtrs = (phzBinEdgs(1:end-1)+phzBinEdgs(2:end))./2;

% SET head velocity bins
hvfBinEdgs = [-25,-5,5,25,80];
hvfBinCtrs = (hvfBinEdgs(1:end-1)+hvfBinEdgs(2:end))./2;

% INITIALLIZE output 
pfs  = cell([numel(phzBinCtrs),numel(hvfBinCtrs)]);
bins = {phzBinCtrs,hvfBinCtrs};

verbose = true;
if verbose,
    disp(['[status]        compute_egohvf_ratemap: processing trial: ',Trial.filebase]);
end

if isempty(units),  return;  end;

if overwrite
% COMPUTE head velocity in the frame of reference of the head
    hvfl = fet_href_HXY(Trial,sampleRate,false,'trb',4);
    hvflBinInds = discretize(hvfl(:,1),hvfBinEdgs);
% TRANSFORM Local Field Potential -> theta phase
    phz = load_theta_phase(Trial,xyz);
    phzBinInds = discretize(phz.data,phzBinEdgs);
% COMPUTE head basis
    hvec = xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    hvec = multiprod(hvec,...
                     [cos(headYawCorrection),-sin(headYawCorrection);...
                      sin(headYawCorrection),cos(headYawCorrection)],...
                     [2,3],...
                     [1,2]);
% GET theta state behaviors, minus rear, minus turn
    thetaState = resample(cast([Trial.stc{'walk+pause&theta'}],'TimeSeries'),xyz);
end%if overwrite
    
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
electrode = 0;

% Don't judge me
if pargs.SmoothingWeights(1)~=2,
    stag = ['_SW',num2str(pargs.SmoothingWeights(1))];
else
    stag = '';
end




for phase = 1:numel(phzBinCtrs)
    for hvf = 1:numel(hvfBinCtrs)    
        for iter = 1:nIter
        % CHECK existence of pfs object
        pargs.tag = ['egofield_theta_phase_hvf_',num2str(hvf),'_',num2str(phase),'_iter_',num2str(iter),stag];
        filepath = fullfile(Trial.spath, [Trial.filebase,'.Pfs.',pargs.tag,'.mat']);
        
        if exist(filepath,'file')
            pfs{phase} = load(filepath).Pfs;
            if overwrite,
                pfs{phase}.purge_savefile();
            else
                continue;
            end% if
        end% if exist
        

        for unit = 1:numel(units)
            if unit==1 
                pargs.spk = copy(spk);
                pargs.states = copy(thetaState);
                pargs.states.label = ['thetaPhz_',num2str(phase)];
                pargs.states.data( phzBinInds~=phase | hvfBinInds~=hvf ) = 0;
                cast(pargs.states,'TimePeriods');
                resInd = WithinRanges(pargs.spk.res,pargs.states.data);
                pargs.spk.res = pargs.spk.res(resInd);
                pargs.spk.clu = pargs.spk.clu(resInd);
            end% if unit==1
            [mxr,mxp] = pft.maxRate(units(unit));
            pfsCenterHR = MTADfet.encapsulate(Trial,                                           ...
                                              bsxfun(                                          ...
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
            end% if unit==1
        end% for unit
        pfTemp.save();
        pfs{phase,hvf} = pfTemp;    
        pfTemp = Trial;
        end%for iter
    end% for hvf
end% for phase
