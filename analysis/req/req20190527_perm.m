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

stl = {'rear','high','low',};
sti = {[2],[3,4],[5,6]};

sigma = 150;
numIter = 101;

overwrite = false;

% ANALYSIS TESTING ---------------------------------------------------------------------------------


% ANALYSIS MAIN ------------------------------------------------------------------------------------


% BATCH Process Trials

for t = 1:numel(Trials),%[24:28];
%for t = [29]
%for t =[24:28];
    %t = 20;    
    Trial = Trials{t}; 
    unitTrialSubset = units{t};        
% LOAD spikes
    spk = Trial.load('spk',sampleRate,'',unitTrialSubset,'deburst');
% LOAD theta placefield
    pft = pfs_2d_theta(Trial,unitTrialSubset);
% LOAD marker positions
    xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
% STCM STATE Matrix
    stcm = stc2mat(Trial.load('stc','msnn_ppsvd_raux'),xyz,states);
% LOAD theta phase
    phz = load_theta_phase(Trial,xyz,sessionList(t).thetaRefGeneral,phzCorrection(t));
% DRZ - Directional Rate Zone
    [hrz,~,drang] = compute_hrz(Trial,unitTrialSubset,pft,'sampleRate',sampleRate);
    [ddz] = compute_ddz(Trial,unitTrialSubset,pft,'sampleRate',sampleRate);
    [ghz] = compute_ghz(Trial,unitTrialSubset,pft,'sampleRate',sampleRate,'sigma',sigma);    
% HEAD vector
    fxyz = filter(copy(xyz),'ButFilter',3,30,'low');
    hvec = fxyz(:,'nose',[1,2])-fxyz(:,'hcom',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    hvec =  multiprod([cos(headRotationCorrection{t}),-sin(headRotationCorrection{t});...
                       sin(headRotationCorrection{t}),cos(headRotationCorrection{t})],...
                      hvec,[1,2],[2,3]);
% REMOVE spikes which occur lateral to the head within 10cm
    pfhr = nan([size(xyz,1),numel(unitTrialSubset),2]);
    for u = 1:numel(unitTrialSubset),%&pft.data.spar>0.15&pft.data.spar<0.3),
        [mxr,mxp] = pft.maxRate(unitTrialSubset(u));
        pfhr(:,u,:) = multiprod(bsxfun(@minus,mxp,sq(fxyz(:,'hcom',[1,2]))),hvec,2,[2,3]);
        ghz(abs(pfhr(:,u,2))>100,u) = nan;
    end
    
    for s = 1:numel(stl)    
        for k = 1:numel(stl)    
            if s==k, 
                continue, 
            end

% APER - All good periods
            aperS = stcm(:,1)==1 & any(stcm(:,sti{s}),2) & ~any(stcm(:,[7,8]),2);
            aperK = stcm(:,1)==1 & any(stcm(:,sti{k}),2) & ~any(stcm(:,[7,8]),2);
            aperBoth = (aperS | aperK);
            
            nS = sum(aperS);
            nK = sum(aperK);
            nB = sum(aperBoth);

            tBins = 1:round(0.1.*xyz.sampleRate):nB;
            tInd = discretize(1:nB,tBins);

            for n = 1:numIter
                subPeriod = false([nB,1]);
                subPeriod(ismember(tInd',randperm(numel(tBins)-1,round(nS/(0.1.*xyz.sampleRate)))')) = true;
                aper = false(size(aperBoth));
                aper(aperBoth) = subPeriod;
        
% PARGS - Default rate map arguments
                pargs = struct('units',              unitTrialSubset,                           ...
                               'states',             'theta',                              ...
                               'overwrite',          false,                                 ...
                               'tag',                ['ddtp2-s',num2str(sigma),'-',         ...
                                                      stl{s},'-',stl{k},'-',num2str(n)],   ...
                               'binDims',            [0.05,pi/8],                          ...
                               'SmoothingWeights',   [2.2,0.8],                            ...
                               'type',               'dt',                                 ...
                               'spkShuffle',         false,                                ...
                               'posShuffle',         false,                                ...
                               'numIter',            1,                                    ...
                               'xyzp',               [],                                   ...
                               'boundaryLimits',     [-1,1;0,2*pi],                        ...
                               'bootstrap',          false,                                ...
                               'halfsample',         false,                                ...
                               'compute_pfs',        @PlotPFCirc,                          ...
                               'autoSaveFlag',       false,                                ...
                               'spk',                spk                                   ...
                               );
                pfs = MTAApfs(Trial,'tag',pargs.tag);
                if all(ismember(pfs.data.clu,unitTrialSubset)) && ~overwrite
                    continue
                end
                
                pfs.purge_savefile();
                pfs = Trial;
                for unit = unitTrialSubset,
                    [mxr,mxp] = pft.maxRate(unit);        
                    pargs.xyzp = copy(xyz);
                    pargs.xyzp.data = [ghz(:,unit==unitTrialSubset),phz.data];
                    pargs.units = unit;
                    pargs.states = MTADepoch([],                                                   ...
                                             [],                                                   ...
                                             ThreshCross(aper                                      ...
                                                         & abs(ddz(:,unit==unitTrialSubset))<300,       ...
                                                         0.5,1),                                   ...
                                             sampleRate,pargs.xyzp.sync.copy(),                    ...
                                             pargs.xyzp.origin,'TimePeriods','sts',[],'thrz','d');
                    pfsArgs = struct2varargin(pargs);
                    pfs = MTAApfs(pfs,pfsArgs{:});    
                    if unit==unitTrialSubset(1),
                        pfs.save();
                    end
                end
                pfs.save();
            end%for n
        end%for k
    end%for s
end%for t


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
% $$$                  cf(@(T,u) MTAApfs(T,u,'tag',['ddtp-','s',num2str(sigma),'-',s]), Trials, units),...
% $$$                  stateLabels);
% $$$ 
% $$$ t = 5;
% $$$ figure,
% $$$ sp = tight_subplot(2,8,0,0.1);
% $$$ sp = reshape(reshape(sp',8,2)',2,8);
% $$$ for u = units{t},
% $$$     for s = 1:8,
% $$$         axes(sp(s*2-1));plot(pftHZTPD{s}{t},u,'mean','text',[],false);
% $$$         axes(sp(s*2));plot(pftHZTPD{s}{t},u,'mean','text',[],false);
% $$$         title(num2str(u));
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end

