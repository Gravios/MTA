%  req20190126
%      Tags: theta phase rate map conditioned on space
%      Status: active
%      Type: Analysis
%      Author: Justin Graboski
%      Final_Forms: NA
%      Project: MjgER2016: behavior-theta-modulation 
%      Description: Compilation of theta phase rate maps given positon and state
%      Protocol:
%      Figures: 

% ANALYSIS MAIN ------------------------------------------------------------------------------------

MjgER2016_load_data();
% help MjgER2016_load_data();

% SET analysis parameters
sampleRate = 250;   % Hz
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};%,'ripple'};

stl = {'theta','rear','high','low'};
sti = {[1],[2],[3,4],[5,6]};

% BATCH Process Trials

for t = 1:23
    %t = 18;    
    Trial = Trials{t}; 
    unitSubset = units{t};        
    subjectId = regexp(Trial.name,'^(\w*)-','tokens');
    subjectId = subjectId{1}{1};
    Trial.lfp.filename = [Trial.name,'.lfp'];

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
    [hrz,~,drang] = compute_hrz(Trial,unitSubset,pft,'sampleRate',sampleRate);
    [ddz] = compute_ddz(Trial,unitSubset,pft,'sampleRate',sampleRate);

    for s = 1:numel(stl)    
% APER - All good periods
        aper = stcm(:,1)==1 & any(stcm(:,sti{s}),2) & ~any(stcm(:,[7,8]),2);
% PARGS - Default rate map arguments
        pargs = struct('units',              unitSubset,                           ...
                       'states',             'theta',                              ...
                       'overwrite',          false,                                ...
                       'tag',                ['hztp-',stl{s}],                     ...
                       'binDims',            [0.1,pi/8],                           ...
                       'SmoothingWeights',   [1.2,0.8],                            ...
                       'type',               'ht',                                 ...
                       'spkShuffle',         false,                                ...
                       'posShuffle',         false,                                ...
                       'numIter',            1001,                                 ...
                       'xyzp',               [],                                   ...
                       'boundaryLimits',     [-1,1;-pi,pi],                        ...
                       'bootstrap',          false,                                ...
                       'halfsample',         true,                                 ...
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
                                     ThreshCross(aper                                      ...
                                                 & abs(ddz(:,unit==unitSubset))<250,       ...
                                                 0.5,1),                                   ...
                                     sampleRate,pargs.xyzp.sync.copy(),                    ...
                                     pargs.xyzp.origin,'TimePeriods','sts',[],'thrz','d');
            pfsArgs = struct2varargin(pargs);
            pfs = MTAApfs(pfs,pfsArgs{:});    
            if unit==unitSubset(1),
                pfs.save();
            end
        end
        pfs.save();
    end% for s
end% for t
