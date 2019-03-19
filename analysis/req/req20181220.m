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
               'tag',                'hbpptbpTR1v5',                       ...
               'binDims',            [0.1,0.1],                            ...
               'SmoothingWeights',   [1.8,1.8],                            ...
               'type',               'xyz',                                ...
               'spkShuffle',         false,                                ...
               'posShuffle',         false,                                ...
               'numIter',            1,                                    ...
               'xyzp',               [],                                   ...
               'boundaryLimits',     [-1.8,1;-0.8,2],                      ...
               'bootstrap',          false,                                ...
               'halfsample',         false,                                ...
               'compute_pfs',        @PlotPFCirc,                          ...
               'autoSaveFlag',       false                                 ...
               );




ndim = numel(pargs.binDims);
binGrid = cell([1,ndim]);    
bins = cell([1,ndim]);
binc = cell([1,ndim]);
for f = 1:ndim,,
    bins{f} = pargs.boundaryLimits(f,1):pargs.binDims(f):pargs.boundaryLimits(f,2);
    binc{f} = (bins{f}(2:end)+bins{f}(1:end-1))./2;
end
[binGrid{:}] = ndgrid(binc{:});
binGrid = cat(ndim+1,binGrid{:});


ntt = zeros([cellfun(@numel,bins)-1,0]);
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
    pargs.xyzp.data = [pchb(:,1),pchb(:,2)];
    pargs.units = unit;
    pargs.states = MTADepoch([],                                                ...
                             [],                                                ...
                             ThreshCross(aper & abs(drz(:,unit==unitSubset))<0.8   ...
                                         & abs(ddz(:,unit==unitSubset))<250,  ...
                                         0.5,1),                           ...
                             sampleRate,pargs.xyzp.sync.copy(),                       ...
                             pargs.xyzp.origin,'TimePeriods','sts',[],'tdrz','d');


    tbinInds = MTADfet.encapsulate(Trial,...
                                   discretize(1:size(pargs.xyzp.data,1),1:sampleRate*3:size(pargs.xyzp.data,1)),...
                                   sampleRate,...
                                   'time bins',...
                                   'tbins',...
                                   't');
    
    pitches = pargs.xyzp(pargs.states,:);
    gridDist = sq(sqrt(sum((repmat(binGrid,[1,1,1,size(pitches,1)])...
                            -repmat(permute(pitches,[3,4,2,1]),[size(binGrid,1),size(binGrid,2),1,1])).^2,3)));

    [gds,gridDistTinds] = sort(gridDist,3);
    
    times = tbinInds(pargs.states);
    
    gds = mat2cell(gds,ones([1,28]),ones([1,28]),size(pitches,1));
    gridDistTinds = mat2cell(gridDistTinds,ones([1,28]),ones([1,28]),size(pitches,1));    
    
    tt = cf(@(dist,ind,time) time(nonzeros(double(sq(dist)<0.1).*sq(ind))),...
            gds,gridDistTinds,repmat({times},size(gds)));
    ntt = cat(3,ntt,cell2mat(cf(@(t) sum(diff(unique(t))>1), tt)));
% $$$     clf,imagesc(binc{:},tt'),axis('xy');
% $$$     caxis([3,10]);
% $$$     title(num2str(unit));
% $$$     waitforbuttonpress();
end


figure();
for u = 1:numel(unitSubset),
    clf,imagesc(binc{:},ntt(:,:,u)'),axis('xy');
    caxis([3,10]);
    title(num2str(unitSubset(u)));
    waitforbuttonpress();
end
