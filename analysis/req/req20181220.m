global MTA_PROJECT_PATH;

MjgER2016_load_data();

% SET analysis parameters
sampleRate = 250;   % Hz
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};%,'ripple'};

thresholds.mazeCenterDist = 380;
thresholds.mazeCenterAng = pi/2;

dataFilePath = fullfile(MTA_PROJECT_PATH,'analysis',...
                        ['req20181220-data-',DataHash({[sessionList.sessionName],sampleRate,states}),'.mat']);


t = 20;    
Trial = Trials{t}; 
unitSubset = units{t};        
subjectId = regexp(Trial.name,'^(\w*)-','tokens');
subjectId = subjectId{1}{1};

pft = pfs_2d_theta(Trial,unitSubset);

stc = Trial.load('stc','msnn_ppsvd_raux');
xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
fxyz = filter(copy(xyz),'ButFilter',3,30,'low');
mazeCenterDist = sqrt(sum(xyz(:,'hcom',[1,2]).^2,3));
mazeCenterAng = circ_dist(atan2(xyz(:,'hcom',2),xyz(:,'hcom',1)),...
                          atan2(diff(xyz(:,{'hcom','nose'},2),1,2),...
                                diff(xyz(:,{'hcom','nose'},1),1,2)));
try
    lfp = load(Trial,'lfp',sessionList(t).thetaRef);
catch err
    lfp = load(Trial,'lfp',sessionList(t).thetaRef);
end

phz = lfp.phase([6,12]);
phz.data = unwrap(phz.data);
phz.resample(xyz);    
phz.data = mod(phz.data+pi,2*pi)-pi;
lfp.resample(xyz);    

hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
[drz,~,drang] = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);
[ddz] = compute_ddz(Trial,unitSubset,pft,'sampleRate',sampleRate);

stcm = stc2mat(stc,xyz,states);

nq = get(Trial.load('nq'),'nq');
edist = nq.eDist(unitSubset)';

spk = Trial.load('spk',sampleRate,'',unitSubset,'deburst');    

[~,sind] = sort(pft.data.clu);
si = subsref(pft.data.si(ismember(pft.data.clu,unitSubset)),substruct('()',{sind}));
spar = subsref(pft.data.spar(ismember(pft.data.clu,unitSubset)),substruct('()',{sind}));    


pchb = fet_HB_pitchB(Trial,sampleRate,false,'trb');        



aper = stcm(:,1)==1 & ~any(stcm(:,[7,8]),2);


pargs = struct('units',              unitSubset,                           ...
               'states',             'theta',                              ...
               'overwrite',          false,                                ...
               'tag',                'hbpptbpFS1v3',                       ...
               'binDims',            [0.1,pi/9,0.1],                       ...
               'SmoothingWeights',   [2.5,1.1,2.5],                        ...
               'type',               'xyz',                                ...
               'spkShuffle',         false,                                ...
               'posShuffle',         false,                                ...
               'numIter',            1,                                    ...
               'xyzp',               [],                                   ...
               'boundaryLimits',     [-1.8,1;-pi,pi;-0.8,2],               ...
               'bootstrap',          false,                                ...
               'halfsample',         false,                                ...
               'compute_pfs',        @PlotPFCirc,                          ...
               'autoSaveFlag',       false                                 ...
               );


pfs = MTAApfs(Trial,'tag',pargs.tag);
pfs.purge_savefile();
pfs = Trial;
for unit = unitSubset,
    [mxr,mxp] = pft.maxRate(unit);        
    pfsCenterHR = MTADfet.encapsulate(Trial,                                               ...
                                      multiprod(bsxfun(@minus,mxp,sq(fxyz(:,'hcom',[1,2]))),...
                                                hvec,2,[2,3]),                             ...
                                      sampleRate,                                          ...
                                      'placefield_center_referenced_to_head',              ...
                                      'pfsCenterHR',                                       ...
                                      'p'                                                  ...
                                      );
    pargs.xyzp = copy(pchb);
    pargs.xyzp.data = [pchb(:,1),phz(:,spk.map(spk.map(:,1)==unit,2)),pchb(:,2)];
    pargs.units = unit;
    pargs.states = MTADepoch([],                                                ...
                             [],                                                ...
                             ThreshCross(aper & abs(drz(:,unit==unitSubset))<0.8   ...
                                         & abs(ddz(:,unit==unitSubset))<250,  ...
                                         0.5,1),                           ...
                             sampleRate,pargs.xyzp.sync.copy(),                       ...
                             pargs.xyzp.origin,'TimePeriods','sts',[],'tdrz','d');
    
    pfsArgs = struct2varargin(pargs);
    pfs = MTAApfs(pfs,pfsArgs{:});    
    if unit==unitSubset(1),
        pfs.save();
    end
end
