



MjgER2016_load_data();
%  Variables:
%      sessionListName, sessionList
%      Trials, units, cluSessionMap, pitchReferenceTrial
%      states, numStates
%      FigDir, interpParPfsp, interpParDfs
%
%  Functions:
%      reshape_eigen_vector

unitCount = sum(cellfun(@numel,units));
phzBins = linpspace(-pi,pi,17);
spkRateFPhz = zeros(

for t = 1:23,
    Trial = Trials{t}; 
    unitSubset = units{t};

% LOAD theta place fields
    pft = pfs_2d_theta(Trial,unitSubset);

% LOAD xyz coordinates of subject markers;
    xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
    fxyz = filter(copy(xyz),'ButFilter',3,30,'low');
    
% COMPUTE polar coordinates of head relative to maze center;
    mazeCenterDist = sqrt(sum(xyz(:,'hcom',[1,2]).^2,3));
    mazeCenterAng = circ_dist(atan2(xyz(:,'hcom',2),xyz(:,'hcom',1)),...
                              atan2(diff(xyz(:,{'hcom','nose'},2),1,2),...
                                    diff(xyz(:,{'hcom','nose'},1),1,2)));
    hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    [drz,~,drang] = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);
    [ddz] = compute_ddz(Trial,unitSubset,pft,'sampleRate',sampleRate);
        
% LOAD State Collection;
    stc = Trial.load('stc','msnn_ppsvd_raux');
    stcm = stc2mat(stc,xyz,states);
    
% LOAD local field potential;
    % error sometimes due to readlink failure, cause: unknown
    try,  lfp = load(Trial,'lfp',sessionList(t).thetaRef);
    catch lfp = load(Trial,'lfp',sessionList(t).thetaRef);
    end    
    
% COMPUTE Theta phase;
    phz = lfp.phase([6,12]);
    phz.data = unwrap(phz.data);
    phz.resample(xyz);    
    phz.data = mod(phz.data+pi,2*pi)-pi;
    lfp.resample(xyz);    
    
% LOAD neuron cluster stats;
    nq = get(Trial.load('nq'),'nq');
    edist = nq.eDist(unitSubset)';

% LOAD spike clusters;
    spk = Trial.load('spk',sampleRate,'',unitSubset,'deburst');    
    
    [~,sind] = sort(pft.data.clu);
    si = subsref(pft.data.si(ismember(pft.data.clu,unitSubset)),substruct('()',{sind}));
    spar = subsref(pft.data.spar(ismember(pft.data.clu,unitSubset)),substruct('()',{sind}));    
    
    aper = stcm(:,1)==1 & ~any(stcm(:,[7,8]),2);
    for unit = unitSubset,
        [mxr,mxp] = pft.maxRate(unit);        
% SELECT only activity within placefield        
        state = MTADepoch([],                                                                   ...
                          [],                                                                   ...
                          ThreshCross(aper & abs(drz(:,unit==unitSubset))<0.8                   ...
                                           & abs(ddz(:,unit==unitSubset))<250,                  ...
                                      0.5,1),                                                   ...
                          sampleRate,pargs.xyzp.sync.copy(),                                    ...
                          pargs.xyzp.origin,'TimePeriods','sts',[],'tdrz','d');

        
% COMPUTE spike rate as a function of theta phase within placefield center
        res = spk(unit);
        res = res(WithinRanges(res,state.data));

        resPhzBinInds = discretize(phz(res,  spk.map(spk.map(:,1)==unit,2)),phzBins);
        thpPhzBinInds = discretize(phz(state,spk.map(spk.map(:,1)==unit,2)),phzBins);
        
        spkCount = accumarray(resPhzBinsInds,ones(size(res)),@sum);
        phzOcc   = accumarray(phzPhzBinsInds,ones(size(res)),@sum);
        
        
    end

end
