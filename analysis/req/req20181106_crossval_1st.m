function [pfs,metaData] = req20181106_crossval_1st(varargin)
%
%function load_behavior_fields_theta_resolved(varargin)
% Input:
%    Trials
%    SessionListName
%    tag
%    units
%    sampleRate
%    stcMode
%    states
%    feature_fcn: function handel, provides 2 features of interest
%    overwriteFlag
%    
% Saved:         
%    tag
%    sessionListName
%    sampleRate
%    stcMode
%    states
%    feature_fcn
%    analysisHash
%
% Output:
%    pfs: either cellarray or MTAApfs object
%

% ANALYSIS DESCRIPTION -----------------------------------------------------------------------------
% PAG : Primary Analysis Goal
%    Deterimine the prefered theta phase for behavior state depedent rate changes associated with 
%    spatial tuning curves
%
% SAG : Secondary Analysis Goal 
% TAG : Tertiary Analysis Goal
%
% VAR - Analysis variables
%    spike time (res)
%    head position (xyz)
%    local field potential (lfp)
%
% TRN - Transformations
%    FFT of lfp -> gmHMM of mean power within 6-12Hz band -> thetaPeriods (tper).
%    HILBERT of 6-12 Hz bandpass filtered lfp -> thetaPhase
%    DISCRETIZE head position at spike times durng theta states -> theta spike occupancy map.
%    DISCRETIZE head position durng theta states -> theta position occupancy map.
%    ELEMENT wise division of theta spike occupancy map and theta position occpancy map -> spatial 
%      rate map.
%    SELECTION of spike times within the 90th percentile within fixed radius around center of 
%      spatial rate map -> place restricted selection of spike times during theta state.
%    CONVERT head markers positions from Cartesian to polar coordinate system.
%    
% MOD : Analysis Model 
%    meanFiringRate(thetaPhase, headPitch)
% CMD : Causal Model
%     {Sensory Information} EC3 -> CA1(theta[pi,2*pi])
%     {Auto Completion}     CA3 -> CA1(theta[0,pi])
%
% 
% PFSTAGS : hptp : head pitch, theta phase
%           hpef : head pitch, egocentric forward
%
% END ANALYSIS DESCRIPTION -------------------------------------------------------------------------



% GLOBALS ------------------------------------------------------------------------------------------
global MTA_PROJECT_PATH;
%---------------------------------------------------------------------------------------------------


%varargin = {};
% DEFARGS ------------------------------------------------------------------------------------------

pargs = struct('units',              [],                                   ...
               'states',             'theta',                              ...
               'overwrite',          false,                                ...
               'tag',                '',                                   ...
               'binDims',            [0.1,pi/8,0.1],                       ...
               'SmoothingWeights',   [1.8,0.6,1.8],                        ...
               'type',               'xyz',                                ...
               'spkShuffle',         false,                                ...
               'posShuffle',         false,                                ...
               'numIter',            1,                                    ...
               'xyzp',               [],                                   ...
               'boundaryLimits',     [-2,0.8;-pi,pi;-0.8,2],               ...
               'bootstrap',          false,                                ...
               'halfsample',         false,                                ...
               'compute_pfs',        @PlotPFCirc,                          ...
               'autoSaveFlag',       false,                                ...
               'spk',                []                                    ...
               );


defargs = struct('Trials',                        {{}},                                          ...
                 'sessionListName',               'MjgER2016',                                   ...
                 'tag',                           'hbpptbpFSCV1',                                ...
                 'units',                         [],                                            ...
                 'sampleRate',                    250,                                           ...
                 'stcMode',                       'msnn_ppsvd_raux',                             ...
                 'states',                        {{'theta','rear','hloc','hpause',              ...
                                                   'lloc','lpause','groom','sit'}},              ...
                 'feature_fcn',                   @fet_HB_pitchB,                                ...
                 'overwrite',                     false,                                         ...
                 'pargs',                         pargs                                          ...
);
[Trials,sessionListName,tag,units,sampleRate,stcMode,states,feature_fcn,overwrite,pargs] =       ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% SET Meta Data ------------------------------------------------------------------------------------
analysisHash = DataHash({sessionListName,tag,sampleRate,stcMode,states,func2str(feature_fcn)});
dataFilePath = fullfile(MTA_PROJECT_PATH,'analysis',[mfilename,'-meta-',analysisHash,'.mat']);
metaVars = {'sessionListName','tag','sampleRate','stcMode','states','feature_fcn','analysisHash'};
%---------------------------------------------------------------------------------------------------



% ANALYSIS MAIN ------------------------------------------------------------------------------------
if ~exist(dataFilePath,'file') || overwrite
% LOAD Trial data
% COLLATE unit info
    sessionList = get_session_list(sessionListName);
    Trials = af(@(s) MTATrial.validate(s), sessionList);
    units = cf(@(T)  select_placefields(T),  Trials); 
    units = req20180123_remove_bad_units(units);
    cluSessionMap = [];
    for u = 1:numel(units)
        cluSessionMap = cat(1,cluSessionMap,[u*ones([numel(units{u}),1]),units{u}(:)]);
    end


    % COMPUTE theta resovled behavior fields
    for t = 1:numel(sessionList),
        Trial = Trials{t}; 
        unitSubset = units{t};        
        subjectId = regexp(Trial.name,'^(\w*)-','tokens');
        subjectId = subjectId{1}{1};

% PFT : theta restricted place fields
        pft = pfs_2d_theta(Trial,unitSubset);        
% XYZ : marker positions        
        xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
% STCM : State Collection Matrix
        stcm = stc2mat(Trial.load('stc',stcMode),xyz,states);        
% LFP : Local Field Potential
        try,   lfp = load(Trial,'lfp',sessionList(t).thetaRefGeneral);
        catch, lfp = load(Trial,'lfp',sessionList(t).thetaRefGeneral);
        end
% $$$         try,   lfp = load(Trial,'lfp',sessionList(t).thetaRef);
% $$$         catch, lfp = load(Trial,'lfp',sessionList(t).thetaRef);
% $$$         end
% PHZ : LFP phase restricted to the theta band [6-12Hz]
        phz = lfp.phase([6,12]);
        phz.data = unwrap(phz.data);
        phz.resample(xyz);    
        phz.data = mod(phz.data+pi,2*pi)-pi;
        lfp.resample(xyz);    
% DRZ : directed rate zone
        drz = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);
% DDZ : directed distance zone        
        ddz = compute_ddz(Trial,unitSubset,pft,'sampleRate',sampleRate);        
% SPK : Unit spike times
        spk = Trial.load('spk',sampleRate,'',unitSubset,'deburst');    
% FET : Behavior space features
        fet = feval(feature_fcn,Trial,sampleRate,false,'trb');
% APER : theta periods without sit or groom
        for s = 2:6
            stcm(find(stcm(:,s),round(sum(nniz(stcm(:,s)))/2),'last'),s) = 0;
        end

        aper = stcm(:,1)==1 & any(stcm(:,2:6),2);

% PARGS : placefield computation arguments 
        
        pargs.units = unitSubset;
        pargs.spk = spk;
        
        pfs = MTAApfs(Trial,'tag',pargs.tag);
        pfs.purge_savefile();
        pfs = Trial;
        for unit = unitSubset,
            pargs.xyzp = copy(fet);
            pargs.xyzp.data = [fet(:,1),phz(:,1),fet(:,2)];
            %pargs.xyzp.data = [fet(:,1),phz(:,spk.map(spk.map(:,1)==unit,2)),fet(:,2)];            
            pargs.units = unit;
            pargs.tag   = tag;
            pargs.states = MTADepoch([],                                                   ...
                                     [],                                                   ...
                                     ThreshCross(aper & abs(drz(:,unit==unitSubset))<0.8   ...
                                                      & abs(ddz(:,unit==unitSubset))<250,  ...
                                                 0.5,0),                                   ...
                                     sampleRate,pargs.xyzp.sync.copy(),                    ...
                                     pargs.xyzp.origin,'TimePeriods','sts',[],'tdrz','d');
            pfsArgs = struct2varargin(pargs);
            pfs = MTAApfs(pfs,pfsArgs{:});    
            if unit==unitSubset(1),  pfs.save();  end
        end%for unit
        pfs.save();
        save(dataFilePath,metaVars{:});
    end%for t 
else
    if iscell(Trials),
        pfs = cf(@(t,u)  MTAApfs(t,u,'tag',tag), Trials,units);
    else
        pfs = MTAApfs(Trials,units,'tag',tag);
    end
end%if ~exist(dataFilePath,'file') || overwrite

metaData = load(dataFilePath);

