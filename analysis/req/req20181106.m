
% ANALYSIS DESCRIPTION -----------------------------------------------------------------------------

% SAG - Secondary Analysis Goal
% TAG - Tertiary Analysis Goal


% PAG : Primary Analysis Goal
%    Deterimine the prefered theta phase for behavior state depedent rate changes associated with spatial
%    tuning curves
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
%    ELEMENT wise division of theta spike occupancy map and theta position occpancy map -> spatial rate map.
%    SELECTION of spike times within the 90th percentile within fixed radius around center of spatial rate
%        map -> place restricted selection of spike times during theta state.
%    CONVERT head markers positions from Cartesian to polar coordinate system.
%    
% MOD : Analysis Model 
%    meanFiringRate(thetaPhase, headPitch)
% CMD : Causal Model
%     {Sensory Information} EC3 -> CA1(theta[pi,2*pi])
%     {Auto Completion}     CA3 -> CA1(theta[0,pi])
%

% END ANALYSIS DESCRIPTION -------------------------------------------------------------------------


% ANALYSIS MAIN ------------------------------------------------------------------------------------
global MTA_PROJECT_PATH;

MjgER2016_load_data();

% SET analysis parameters
sampleRate = 250;   % Hz
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};%,'ripple'};

thresholds.mazeCenterDist = 380;
thresholds.mazeCenterAng = pi/2;

dataFilePath = fullfile(MTA_PROJECT_PATH,'analysis',...
                        ['req20181106-data-',DataHash({[sessionList.sessionName],sampleRate,states}),'.mat']);



for t = 1:23;
    %t = 20;    
    Trial = Trials{t}; 
    unitSubset = units{t};        
    subjectId = regexp(Trial.name,'^(\w*)-','tokens');
    subjectId = subjectId{1}{1};

    switch subjectId,
      case 'jg05'        
        pft = MTAApfs(Trial,unitSubset,[],[],'CA1thetaCA3inputPhase');
      otherwise  
        pft = pfs_2d_theta(Trial,unitSubset);
    end

    stc = Trial.load('stc','msnn_ppsvd_raux');
    xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
    fxyz = filter(copy(xyz),'ButFilter',3,30,'low');
    mazeCenterDist = sqrt(sum(xyz(:,'hcom',[1,2]).^2,3));
    mazeCenterAng = circ_dist(atan2(xyz(:,'hcom',2),xyz(:,'hcom',1)),...
                              atan2(diff(xyz(:,{'hcom','nose'},2),1,2),...
                                    diff(xyz(:,{'hcom','nose'},1),1,2)));
    lfp = load(Trial,'lfp',sessionList(t).thetaRef);
    phz = lfp.phase([6,12]);
    phz.data = unwrap(phz.data);
    phz.resample(xyz);    
    phz.data = mod(phz.data+pi,2*pi)-pi;
    lfp.resample(xyz);    

    tvec = circshift(fxyz(:,'hcom',[1,2]),-50)-fxyz(:,'hcom',[1,2]);
    tvec = sq(bsxfun(@rdivide,tvec,sqrt(sum(tvec.^2,3))));
    tvec = cat(3,tvec,sq(tvec)*[0,-1;1,0]);

    tvecb = -circshift(fxyz(:,'hcom',[1,2]),50)-fxyz(:,'hcom',[1,2]);
    tvecb = sq(bsxfun(@rdivide,tvecb,sqrt(sum(tvecb.^2,3))));
    tvecb = cat(3,tvecb,sq(tvecb)*[0,-1;1,0]);
    
    hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    [drz,~,drang] = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);
    [ddz,~,drang] = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);
    
    vxy = vel(filter(copy(xyz),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
    vxy.data(vxy.data<1e-3) = 1e-3;
    vxy.data = log10(vxy.data);

    headAngle = atan2(diff(fxyz(:,{'hcom','nose'},2),1,2),diff(fxyz(:,{'hcom','nose'},1),1,2));
    vang = [circ_dist(circshift(headAngle,-25),headAngle),circ_dist(circshift(headAngle,-75),headAngle),...
            circ_dist(circshift(headAngle,-125),headAngle),circ_dist(circshift(headAngle,-200),headAngle)];
    vangb = [-circ_dist(circshift(headAngle,25),headAngle),-circ_dist(circshift(headAngle,75),headAngle),...
             -circ_dist(circshift(headAngle,125),headAngle),-circ_dist(circshift(headAngle,200),headAngle)];
    
    stcm = stc2mat(stc,xyz,states);

    nq = get(Trial.load('nq'),'nq');
    edist = nq.eDist(unitSubset)';
    
    spk = Trial.load('spk',sampleRate,'',unitSubset,'deburst');    
    
    [~,sind] = sort(pft.data.clu);
    si = subsref(pft.data.si(ismember(pft.data.clu,unitSubset)),substruct('()',{sind}));
    spar = subsref(pft.data.spar(ismember(pft.data.clu,unitSubset)),substruct('()',{sind}));    

    pch = fet_HB_pitch(Trial,sampleRate,false,'trb');
    pch.data  = pch.data(:,3);
    pch.name  = 'headPitch_thetaPhase';
    pch.label = 'fet_hptp';
    pch.key   = 'p';
    


    tper = stcm(:,1)==1 & ~(stcm(:,2)==2|stcm(:,7)==7|stcm(:,8)==8);

    pargs = struct('units',              unitSubset,                           ...
                   'states',             'theta',                              ...
                   'overwrite',          false,                                ...
                   'tag',                'hptp',                               ...
                   'binDims',            [ 0.1, 0.3],                          ...
                   'SmoothingWeights',   [1.5,1.5],                            ...
                   'type',               'xy',                                 ...
                   'spkShuffle',         false,                                ...
                   'posShuffle',         false,                                ...
                   'numIter',            1,                                    ...
                   'xyzp',               [],                                   ...
                   'boundaryLimits',     [-2,2;-pi,pi],                        ...
                   'bootstrap',          false,                                ...
                   'halfsample',         false,                                ...
                   'compute_pfs',        @PlotPFCirc,                          ...
                   'autoSaveFlag',       false                                 ...
                   );
    
    pfs = Trial;
    for unit = unitSubset,
        pargs.xyzp = copy(pch);
        pargs.xyzp.data = [pch.data,phz(:,spk.map(spk.map(:,1)==unit,2))];
        pargs.units = unit;
        pargs.states = MTADepoch([],                                                ...
                                 [],                                                ...
                                 ThreshCross(tper                                   ...
                                             & abs(drz(:,unit==unitSubset))<0.8   ...
                                             & abs(ddz(:,unit==unitSubset))<250,  ...
                                             0.5,1),                                ...
                                 sampleRate,xyzp.sync.copy(),                       ...
                                 xyzp.origin,'TimePeriods','sts',[],'tdrz','d');
        
        pfsArgs = struct2varargin(pargs);
        pfs = MTAApfs(pfs,pfsArgs{:});    
        if unit==unitSubset(1),
            pfs.save();
        end
    end
    pfs.save();
end

pfd = MTAApfs(Trial,'tag','hptp');

pft = pfs_2d_theta(Trial);

figure,
for u = unitSubset,
    clf();
    subplot(121);
    plot(pft,u,1,'colorbar',[],false);
    title(num2str(u));
    subplot(122);
    plot(pfd,u,1,'colorbar',[],false);
    ylabel('theta phase')
    xlabel('head pitch');
    waitforbuttonpress();    
end


rmap = pfd.data.rateMap;
rmap(~nniz(rmap(:)))=0;
[U,S,V] = svd(rmap',0);

figure,
for i = 1:10,
    subplot(2,5,i);
    imagesc(reshape(V(:,i),pfd.adata.binSizes')');
    axis('xy');
    caxis([-0.06,0.06]);
end


[~,V,FSr,VT] = erpPCA(rmap',10);

figure,
for i = 1:10,
    subplot(2,5,i);
    imagesc(reshape(V(:,i),pfd.adata.binSizes')');
    axis('xy');
    caxis([-0.03,0.03]);
end


% SET interp parameters                
interpParPfs = struct('bins',{{linspace(-500,500,50),...
                    linspace(-500,500,50),...
                    linspace(  -2,  2,50)}},...
                      'nanMaskThreshold', 0.1,...
                      'methodNanMap',     'cubic',...
                      'methodRateMap',    'cubic');
% --------------------------------------------------------------------------------------------------

% END ANALYSIS MAIN --------------------------------------------------------------------------------