% req20180222 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Keywords: Placefields, drz, mds, t-sne
%  Description: Compute rate map of reduced dimesion subspace
%               see fet_rds.m
%  Bugs: NA


overwrite = false;
display = false;
sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';

states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',     ...
          'pause&theta','lpause&theta','hpause&theta',            ...
          'theta-groom-sit'};
statesCcg = {'loc','lloc','hloc','rear','pause','lpause','hpause',...
             'theta-groom-sit'};

FigDir = create_directory('/storage/gravio/figures/analysis/placefields_nonSpatialFeatures'); 

% LOAD Trials
% COMPUTE placefield unit stats
Trials  = af(@(S)  MTATrial.validate(S),   sessionList);
          cf(@(T)  T.load('nq'),           Trials);
numStates = numel(states);

units = cf(@(T)  select_placefields(T),  Trials);
        

tind = 1;
Trial = Trials{tind};

xyz = preproc_xyz(Trial,'trb');
pft = pfs_2d_theta(Trial,'overwrite',false);
pch = fet_HB_pitchB(Trial,[],[],[],pitchReferenceTrial);
drz = compute_drz(Trial,units{tind},pft);%,pfstats);
tper =[Trial.stc{'theta-groom-sit'}];
tper.resample(xyz);

% SET DBSTOP
dbstop in PlotPF.m at 78

pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units   = units{tind};
pargs.numIter = 1001;
pargs.halfsample = true;
pargs.tag            = 'DRZxHBPITCHxBPITCH_test';
pargs.boundaryLimits = [-pi/2,pi/2;-2,pi/2];
pargs.binDims        = [0.1,0.1];
pargs.SmoothingWeights = [2,2];
if overwrite,  
    pargs.overwrite = true;    
    pargs.xyzp = MTADxyz('data',pch.data,'sampleRate',xyz.sampleRate);
    drzState = {};
    for u = 1:numel(units{tind});
        dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                         xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
        drzState{u} = dper&tper;
        pargs.units  = units{tind}(u);
        pargs.states = drzState{u};
        pfsArgs = struct2varargin(pargs);
        MTAApfs(Trial,pfsArgs{:});
    end
end
pargs.states    = 'tdrz';
pargs.units     = units{tind};
pargs.overwrite = false;
pfsArgs = struct2varargin(pargs);
pfd{tind} = MTAApfs(Trial,pfsArgs{:});




% @ PlotPF
RateMap = SCount./SOcc;
RateMapOri = RateMap;

figure();
subplot(131);imagesc(RateMapOri')

rSmoother = Smoother;
rSmoother(rSmoother<1e-6) = 0;
rSmoother = rSmoother./sum(rSmoother(:));

rSOcc   = convn(Occupancy, rSmoother,'same');
rSCount = convn(SpikeCount,rSmoother,'same');


rRateMap = rSCount./rSOcc;
RateMapRed = rRateMap;
subplot(132);imagesc(RateMapRed')

RateMap(gtind) = SCount(gtind)./SOcc(gtind);
RateMap = RateMap(:);
RateMap(~gtind) = nan;
subplot(133);imagesc(reshape(RateMap,size(SCount))'),axis xy


pad_matrix = @(m,p) cat(2,...
                      zeros([size(m,1)+2*p(1),p(2)]),...
                      cat(1,zeros([p(1),size(m,2)]),m,zeros([p(1),size(m,2)])),...
                      zeros([size(m,1)+2*p(1),p(2)]));
rSmootherPadded  = pad_matrix(rSmoother, round(SmoothingWeights.*5));
OccupancyPadded  = pad_matrix(Occupancy, round(SmoothingWeights.*5));
SpikeCountPadded = pad_matrix(SpikeCount,round(SmoothingWeights.*5));
rSOccP   = convn(OccupancyPadded, rSmootherPadded,'same');
rSCountP = convn(SpikeCountPadded,rSmootherPadded,'same');

rRateMapP = rSCountP./rSOccP;
RateMapRedP = rRateMapP;
subplot(133);imagesc(RateMapRedP')


%subplot(133);imagesc([RateMapOri-RateMapRed]')


% t-SNE version 

Trial = MTATrial.validate('jg05-20120317.cof.all');
stc = Trial.load('stc','msnn_ppsvd_raux');

fet = fet_rds(Trial);

sfet = fet.copy;
sfet.resample(12);

states = {'loc','rear','pause'};
smat = sum(stc2mat(stc,sfet,states),2);
cmat = [0,0,1;...
        1,0,1;...
        0,1,1];


sind = logical(smat);
sfet = sfet(sind,:);
labels = cmat(smat(sind),:);
out = tsne(sfet,labels,2,3,80);


figure();
plot(out(:,1),out(:,2),'.');
xlim([-120,120]);
ylim([-120,120]);