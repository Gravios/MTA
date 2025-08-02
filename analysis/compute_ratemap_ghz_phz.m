% req20190527
%    Tags: distance place field transform
%    Status: active
%    Type: analysis
%    Author: Justin Graboski
%    Final_Forms: NA
%    Project: MjgER2016: behavior-theta-modulation
%    Description: gaussian distance transform to compute linearized 2d trajectories rate maps
%    Protocol: 1. TRANSFORM distance to max normalized gaussian map
%              2. COMPUTE rate map 
%
%    Figures:

% DEFAULT VARS -------------------------------------------------------------------------------------
configure_default_args();
MjgER2016_load_data();

% SET analysis parameters
sampleRate = 250;   % Hz
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};%,'ripple'};

stl = {'theta','rear','high','low','hloc','hpause','lloc','lpause'};
sti = {[1],[2],[3,4],[5,6],3,4,5,6};

sigma = 150;
% ANALYSIS TESTING ---------------------------------------------------------------------------------


% ANALYSIS MAIN ------------------------------------------------------------------------------------


% BATCH Process Trials
for t = 1:numel(Trials)
    %t = 20;    
    Trial = Trials{t}; 
    unitSubset = units{t};        
    subjectId = regexp(Trial.name,'^(\w*)-','tokens');
    subjectId = subjectId{1}{1};
% LOAD spikes
    spk = Trial.load('spk',sampleRate,'',unitSubset,'deburst');
% LOAD theta placefield
    pft = pfs_2d_theta(Trial,unitSubset);
% LOAD marker positions
    xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
% LOAD State Matrix
    stcm = stc2mat(Trial.load('stc','msnn_ppsvd_raux'),xyz,states);
% LOAD theta phase
    phz = load_theta_phase(Trial,xyz,sessionList(t).thetaRefGeneral,phzCorrection(t));
% DRZ - Directional Rate Zone
    %[hrz,~,drang] = compute_hrz(Trial,unitSubset,pft,'sampleRate',sampleRate);
    [ddz] = compute_ddz(Trial,unitSubset,pft,'sampleRate',sampleRate);

    [ghz] = compute_ghz(Trial,unitSubset,pft,'sampleRate',sampleRate,'sigma',sigma);
% HEAD vector
    fxyz = filter(copy(xyz),'ButFilter',3,30,'low');
    hvec = fxyz(:,'nose',[1,2])-fxyz(:,'hcom',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    hvec =  multiprod([cos(headRotationCorrection{t}),-sin(headRotationCorrection{t});...
                       sin(headRotationCorrection{t}),cos(headRotationCorrection{t})],...
                      hvec,[1,2],[2,3]);
% REMOVE spikes which occur lateral to the head within 10cm
    pfhr = nan([size(xyz,1),numel(unitSubset),2]);
    for u = 1:numel(unitSubset)
        [mxr,mxp] = pft.maxRate(unitSubset(u));
        pfhr(:,u,:) = multiprod(bsxfun(@minus,mxp,sq(fxyz(:,'hcom',[1,2]))),hvec,2,[2,3]);
        ghz(abs(pfhr(:,u,2))>100,u) = nan;
    end
    
    for s = 1:numel(stl)    
% APER - All good periods
        aper = stcm(:,1)==1 & any(stcm(:,sti{s}),2) & ~any(stcm(:,[7,8]),2);
% PARGS - Default rate map arguments
        pargs = struct('units',              unitSubset,                           ...
                       'states',             'theta',                              ...
                       'overwrite',          true,                                 ...
                       'tag',                ['ddtp2-s',num2str(sigma),'-',stl{s}],...
                       'binDims',            [0.05,pi/8],                          ...
                       'SmoothingWeights',   [2.2,0.8],                            ...
                       'type',               'dt',                                 ...
                       'spkShuffle',         false,                                ...
                       'posShuffle',         false,                                ...
                       'numIter',            101,                                  ...
                       'xyzp',               [],                                   ...
                       'boundaryLimits',     [-1,1;0,2*pi],                        ...
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
            pargs.xyzp.data = [ghz(:,unit==unitSubset),phz.data];
            pargs.units = unit;
            pargs.states = MTADepoch([],                                                   ...
                                     [],                                                   ...
                                     ThreshCross(aper                                      ...
                                                 & abs(ddz(:,unit==unitSubset))<300,       ...
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

% $$$ figure,
% $$$ sp = tight_subplot(2,1,0,0.1);
% $$$ for u = units{t},
% $$$ axes(sp(1));plot(pfs,u,'mean','text',[],false);
% $$$ axes(sp(2));plot(pfs,u,'mean','text',[],false);
% $$$ title(num2str(u));
% $$$ waitforbuttonpress();
% $$$ end
% $$$ 
% $$$ 
% $$$ pftHZTPD    = cf(@(s) ...
% $$$                  cf(@(T,u) MTAApfs(T,u,'tag',['ddtp2-','s',num2str(sigma),'-',s]), Trials, units),...
% $$$                  stl);
% $$$ 
% $$$ t = 20;
% $$$ figure,
% $$$ sp = tight_subplot(2,8,0,0.1);
% $$$ sp = reshape(reshape(sp',8,2)',2,8);
% $$$ for u = units{t},
% $$$     for s = 1:8,
% $$$         axes(sp(s*2-1));plot(pftHZTPD{s}{t},u,'mean','text',[],false);
% $$$         title(stl{s});
% $$$         axes(sp(s*2));plot(pftHZTPD{s}{t},u,'mean','text',[],false);
% $$$         title(num2str(u));
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end

