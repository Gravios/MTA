function [pfs] = compute_egohbahvl_radMaze_ratemap(Trial,units,xyz,spk,pft,overwrite)

sampleRate = xyz.sampleRate;
binPhzs = linspace(0,2*pi,4);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
pfs = cell([1,numel(binPhzc)]);
verbose = true;

if verbose,
    disp(['[status]        compute_egohbahvl_ratemap: processing trial: ',Trial.filebase]);
end

if isempty(units),
    return;
end;% if

% LOAD VARIABLES -------------------------------------------------------------------------

% COMPUTE sign of the lateral velocity of the head
fhrvfl = fet_href_HXY(Trial,sampleRate,false,'trb',4);
hvlSign = MTADfet.encapsulate(Trial,sign(fhrvfl(:,2)),sampleRate,'headLatVel','hvl','l');
% COMPUTE head maze angle              -> partition feature
hma = fet_hma(Trial,sampleRate,false,[],xyz);

hdir = sq(xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2]));
hdir = atan2(hdir(:,2),hdir(:,1));

% COMPUTE head body angle               -> partiton feature
hba = fet_hba(Trial, sampleRate);
hbaSign = copy(hba);
hbaSign.data = sign(hbaSign.data);   
% TRANSFORM Local Field Potential       -> theta phase
phz = load_theta_phase(Trial, xyz);
% TRANFORM anteroposterior head vector  -> rotation matrix
rot = transform_vector_to_rotation_matrix( ...
        xyz,{'hcom','nose'}, Trial.meta.correction.headYaw);
% GET theta state behaviors, minus rear -> computational periods
thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);


% COMPUTE PLACEFIELDS --------------------------------------------------------------------
pfTemp = Trial;

pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.units        = units;
pargs.tag          = 'egofield';
pargs.binDims      = [20, 20, 0.6];
pargs.SmoothingWeights = [2.5, 2.5, 0.8];
pargs.halfsample   = false;
pargs.numIter      = 1;   
pargs.boundaryLimits = [-400,400;-400,400;-1.5,1.5];
pargs.states       = '';
pargs.overwrite    = true;
pargs.autoSaveFlag = false;    
electrode = 0;

for phase = 1:numel(binPhzc)

    % CHECK existence of pfs save file;
    tagbase = 'egofield_hbahvl_radMaze_theta_phase_';
    pargs.tag = [tagbase,num2str(phase)];
    filepath = fullfile(Trial.spath,...
                        [Trial.filebase,'.Pfs.',tagbase,num2str(phase),'.mat']);
    if exist(filepath,'file'),
        pfs{phase} = load(filepath).Pfs;
        if overwrite,
            pfs{phase}.purge_savefile();
        else,
            continue;
        end;% if
    end;% if exist
    
        
    for unit = 1:numel(units),

        % COMPUTE place fields
        [mxr,mxp] = pft.maxRate(units(unit));

        hpdiff = bsxfun(@circ_dist,hdir,atan2(mxp(2),mxp(1)));

        pargs.spk = copy(spk);
        pargs.spk.res = pargs.spk.res(pargs.spk.clu==units(unit));
        pargs.spk.clu = pargs.spk.clu(pargs.spk.clu==units(unit));
        pargs.states = copy(thetaState);
        pargs.states.label = ['thetaPhz_',num2str(phase)];        
        pargs.states.data(  (phz.data < binPhzs(phase) )                           ...
                            | (phz.data >= binPhzs(phase+1)) ... %& hvlSign.data ~= hbaSign.data                           ...
                            & ~(-pi/4 < hpdiff & hpdiff < pi/4 )) = 0;

        pargs.states.data = logical(pargs.states.data);
        %cast(pargs.states,'TimePeriods');
        %resInd = WithinRanges(pargs.spk.res, pargs.states.data);
        resInd = logical(pargs.states.data(pargs.spk.res));
        pargs.spk.res = pargs.spk.res(resInd);
        pargs.spk.clu = pargs.spk.clu(resInd);

        pargs.states.cast('TimePeriods');
        
        pfsCenterHR = MTADfet.encapsulate(Trial,                                       ...
                                          [multiprod(bsxfun(@minus,                    ...
                                                            mxp,                       ...
                                                            sq(xyz(:,'hcom',[1,2]))),  ...
                                                     rot,2,[2,3]),                     ...
                                           hba.data],                                  ...
                                          sampleRate,                                  ...
                                          'egocentric_placefield',                     ...
                                          'egopfs',                                    ...
                                          'p'                                          ...
                                          );
        pargs.xyzp = pfsCenterHR;
        pargs.units  = units(unit);
        pfsArgs = struct2varargin(pargs);
        pfTemp = MTAApfs(pfTemp,pfsArgs{:});

        % REMOVE old save file and save new file
        if unit==1,
            try
                pfTemp.purge_savefile();
            end
            pfTemp.save();        
        end% if 
    end;% for unit

    % SAVE results
    pfTemp.save();
    
    % ADD to output
    pfs{phase} = pfTemp;
    
    % Setup for next phase
    pfTemp = Trial;
    
end;% for phase
