% ANALYSIS MAIN ------------------------------------------------------------------------------------

MjgER2016_load_data();
%  Variables:
%      Trials
%      units
%      cluSessionMap
%      pitchReferenceTrial
%      FigDir
%      sessionListName
%      sessionList
%      states
%      numStates
%      interpParPfsp
%      interpParDfs



% SET analysis parameters
sampleRate = 250;   % Hz
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};%,'ripple'};


% BATCH Process Trials
%for t = 1:23
for t = 17:23    
    %t = 20;    
    Trial = Trials{t}; 
    unitSubset = units{t};        
    subjectId = regexp(Trial.name,'^(\w*)-','tokens');
    subjectId = subjectId{1}{1};

    spk = Trial.load('spk',sampleRate,'',unitSubset,'deburst');

    pft = pfs_2d_theta(Trial,unitSubset);

    xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);

% STCM STATE Matrix
    stcm = stc2mat(Trial.load('stc','msnn_ppsvd_raux'),xyz,states);
    
% LFP - Local Field Potential
    try,   lfp = load(Trial,'lfp',sessionList(t).thetaRef);
    catch, lfp = load(Trial,'lfp',sessionList(t).thetaRef);
    end
 
% PHZ - LFP phase within theta band 
    phz = lfp.phase([6,12]);
    phz.data = unwrap(phz.data);
    phz.resample(xyz);    
    phz.data = mod(phz.data+pi,2*pi)-pi;
    lfp.resample(xyz);    

% DRZ - Directional Rate Zone
    [drz,~,drang] = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);
    [ddz] = compute_ddz(Trial,unitSubset,pft,'sampleRate',sampleRate);    
    [hrz] = compute_hrz(Trial,unitSubset,pft,'sampleRate',sampleRate);    

% APER - All good periods
    aper = stcm(:,1)==1 & ~any(stcm(:,[7,8]),2);
    
% PARGS - Default rate map arguments
    bstates = {1,2,[3,4],[5,6]};
    bstatesLabels = {'theta','rear','highBhv','lowBhv'};

for s = 1:numel(bstates)
    pargs = struct('units',              unitSubset,                           ...
                   'states',             bstates{s},                           ...
                   'overwrite',          false,                                ...
                   'tag',                ['HRZxTHP_',bstatesLabels{s}],        ...
                   'binDims',            [0.1,pi/8],                           ...
                   'SmoothingWeights',   [],                                   ...
                   'type',               'ht',                                 ...
                   'spkShuffle',         false,                                ...
                   'posShuffle',         false,                                ...
                   'numIter',            1,                                    ...
                   'xyzp',               [],                                   ...
                   'boundaryLimits',     [-1,1;-pi,pi],                        ...
                   'bootstrap',          false,                                ...
                   'halfsample',         false,                                ...
                   'compute_pfs',        @PlotPFCirc,                          ...
                   'autoSaveFlag',       false,                                ...
                   'spk',                spk                                   ...
                   );
    
    
    pfs = MTAApfs(Trial,'tag',pargs.tag);
    pfs.purge_savefile();
    pfs = Trial;
    for unit = unitSubset,
        [mxr,mxp] = pft.maxRate(unit);        
        pargs.xyzp = copy(xyz);
        pargs.xyzp.data = [hrz(:,unit==unitSubset),phz(:,spk.map(spk.map(:,1)==unit,2))];
        pargs.units = unit;
        pargs.states = MTADepoch([],                                                   ...
                                 [],                                                   ...
                                 ThreshCross(any(stcm(:,bstates{s}),2) & aper          ...
                                             & abs(ddz(:,unit==unitSubset))<250,       ...
                                             0.5,1),                                   ...
                                 sampleRate,pargs.xyzp.sync.copy(),                    ...
                                 pargs.xyzp.origin,'TimePeriods','sts',[],'tdrz','d');
        
        pfsArgs = struct2varargin(pargs);
        pfs = MTAApfs(pfs,pfsArgs{:});    
        if unit==unitSubset(1),
            pfs.save();
        end
    end
    pfs.save();
end

    
end

