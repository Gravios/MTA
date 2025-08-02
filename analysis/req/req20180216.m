% req20180216 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Keywords: Placefields, drz, behavioral transition
%  Description: Computing the mean spike theta phase as a function
%               of drz and behavior onset/offset
%  Bugs: NA


sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';
generateFigureParts = false;
marker = 'nose';


if generateFigureParts,
    FigDir = create_directory('/storage/gravio/figures/analysis/parts/placefields');
else,
    FigDir = create_directory('/storage/gravio/figures/analysis/placefields'); 
end        
        
Trials  = af(@(t)  MTATrial.validate(t),   sessionList);


statesTheta = {'loc&theta','lloc&theta','hloc&theta','rear&theta',     ...
               'pause&theta','lpause&theta','hpause&theta',            ...
               'theta-groom-sit'};
states      = {'loc','lloc','hloc','rear','pause','lpause','hpause',...
               'theta-groom-sit'};
numStates = numel(states);
% LOAD theta placefields
pft = cf(@(T,u)  pfs_2d_theta(T,u,'overwrite',false),  Trials,units);
% LOAD BHV fields and stats
[pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(Trials,[],[],false,true);


t = 20;

Trial = Trials{t};


xyz = preproc_xyz(Trial,'trb');                            % LOAD     Marker positions
xyz.filter('RectFilter');                                  % FILTER   marker positions
ang = create(MTADang,Trial,xyz);                           % COMPUTE  Intermarker spherical coordinates

units = select_placefields(Trial);                         % SELECT   subset of hippocampal neurons with high spatial information
pft = pfs_2d_theta(Trial,units);                           % COMPUTE  Expected Neuron firing rate given position in theta state
pfs = pfs_2d_states(Trial,units);
spk = create( Trial.spk.copy(), Trial, xyz.sampleRate,...  % LOAD     Spike identities and times
             'theta-groom-sit', units, 'deburst');         %                [ theta state; remove bursty spike ( isi < 10ms )]
drz = compute_drz(Trial,units,pft);                        % COMPUTE  Directional Rate Zone
ddz = compute_ddz(Trial,units,pft);                        % COMPUTE  Directional Distance Zone

lfp = Trial.load('lfp',sessionList(t).thetaRef);           % LOAD     Local Field Potential(LFP) of CA1pyr for each electrode shank
lfp.resample(xyz);                                         % RESAMPLE LFP to match xyz
phz = lfp.phase([6,12]);                                   % COMPUTE  LFP theta phase 

pch = fet_HB_pitchB(Trial);

vxy = xyz.vel({'spine_lower','hcom'},[1,2]);               % 
axy = vxy.copy();                                          % 
axy.data = circshift(axy.data,-1)-circshift(axy.data,1);   % 
vxy.data(vxy.data<1e-3) = 1e-4;                            % CLIP      lowerbound (<1e-3) of speed to 1e-4
vxy.data = log10(vxy.data);                                % TRANSFORM speed to log10 scale

spkll = create( Trial.spk.copy(), Trial, xyz.sampleRate,...  % LOAD     Spike identities and times
             'lloc&theta', [], 'deburst');         %                [ theta state; remove bursty spike ( isi < 10ms )]



binSize = 1;
halfBins = 60;
normalization = 'hz';
u = [79,129];
u = [79,119];
u = [151,119];
u = [79,85];
u = [79,151];
u = [119,85];
u = [129,111];

 
u = [151,85];
u = [151,79];
u = [151,119];

u = [119,79];
u = [102,85];
u = [102,73];
u = [73,85];

u = [103,85];
u = [31,139];
u = [103,119];

u = [80,35];
u = [116,35];
u = [116,80];
u = [80,134];
u = [116,134];
u = [35,134];

u = [63,104];

u = [35,68];  %HL,H
u = [68,80];  % H,L
u = [68,81];  % H,L
u = [80,81];  % L,L
u = [68,116]; %
u = [80,116];
u = [81,116];
u = [110,116]; %HL

% OTL
u = [83,104];

% CTL
u = [21,138]; % L,L
u = [21,31]; % L,HL

% OR
u = [72,79]; % LH,H
u = [79,111]; % H,HL
u = [79,119]; % H,HLR
u = [111,119]; % H,HLR
u = [25,79];  % HL,H
u = [25,111];  % HL,H
u = [25,72];  % HL,H

% OBR
u = [28,109];% LP, LL
u = [28,20];% LP, LLH
u = [20,109];% LLH, LL
u = [20,141];% LLH, LH
u = [109,141];% LL, LH
u = [20,140];% LLH, LLHP
u = [109,105]; %LLH, LL
u = [20,105]; %LLH, LL
u = [141,98]; %LLH, H
u = [105,98]; %LLH, H
u = [59,144]; %LL, HL
u = [59,141]; %LL, LH
u = [20,113]; %LL, LH
u = [59,105]; %LL, LL

uRes = spk(u(1));
uRes(ddz(uRes,u(1)==units)>300) = [];
iRes = spk(u(2));
iRes(ddz(iRes,u(1)==units)>300) = [];
% $$$ 
% $$$ uRes = spk(u(1));
% $$$ iRes = spk(u(2));
% $$$ 
[tccg,tbin] = CCG([uRes;iRes],...
                  [ones(size(uRes));2*ones(size(iRes))],...
                  binSize,halfBins,spk.sampleRate,[1,2],normalization);
figure,
drawnow();
set(gcf,'Position',[0,500,1000,400]);
subplot(231);plot(pft,u(1),'mean',true,[],false,0.99);title(num2str(u(1)))
subplot(232);plot(pft,u(2),'mean',true,[],false,0.99);title(num2str(u(2)))
subplot(234);bar(tbin,tccg(:,1,1));xlim(tbin([1,end]));ylim([0,18]);
subplot(235);bar(tbin,tccg(:,2,2));xlim(tbin([1,end]));ylim([0,18]);
subplot(236);bar(tbin,tccg(:,1,2));xlim(tbin([1,end]));ylim([0,18]);


figure,plot(drz(uRes,units==79),drz(uRes,units==119),'.');
figure,plot(drz(uRes,units==79),drz(uRes,units==85),'.');
figure,
subplot(121);plot(drz(uRes,units==79),drz(uRes,units==119),'.');
subplot(122);plot(drz(iRes,units==79),drz(iRes,units==119),'.');


% SET state
stsper = Trial.stc{'r'};
ststrans = xyz.copy('clear');
ststrans.data = nan([size(xyz,1),1]);
for period = stsper.data(2:end-1,2)',
    ststrans.data(period-60:period+60) = linspace(-1,1,121);
end

fet = ststrans;
%fet = pch;
%fet = vxy;
binx = linspace([-1,1,10]);                                % SET  bins for x axis drz 
biny = linspace([-1,1,10]);                                % SET  bins for y axis trans
%biny = linspace([-2,pi/2,20]);                             % SET  bins for y axis speed
%biny = linspace([-2,2,20]);                                % SET  bins for y axis speed
%biny = linspace([-0.5,0.5,10]);                            % SET  bins for y axis speed
mphzDVSize = [numel(binx),numel(biny)];                    % NOTE final map size
numIter = 1000;                                            % SET  count for bootstrap
figure,
for unit = units
    res = spk(unit);
    res(abs(ddz(res,unit==units))>250) = [];    
    subs = [discretize(drz(res,unit==units),binx),...
            discretize(fet(res,1),biny)];

    mphz = phz(res,spk.map(unit==spk.map(:,1),2));
    mphz(any(isnan(subs),2)) = [];
    subs(any(isnan(subs),2),:) = [];
    if length(subs)>3,
        
    mphzDV = accumarray(subs,mphz,mphzDVSize,@circ_mean,nan);    
    mphzDVOcc = accumarray(subs,ones([length(subs),1]),mphzDVSize,@sum);
    mphzDVBoot = nan([mphzDVSize,numIter]);
% $$$     for boot = 1:numIter,1,1,1
% $$$         ind = randi(round(0.5*length(subs)),[round(0.5*length(subs)),1]);
% $$$         mphzDVBoot(:,:,boot) = accumarray(subs(ind,:),mphz(ind),mphzDVSize,@circ_mean,nan);
% $$$     end
    clf();
`    subplot(1,4,[1,2]);plot(pft,unit,'mean',true,[],false,0.99);
    subplot(143);imagescnan({binx,biny,mphzDVOcc'});axis('xy');
    subplot(144);imagescnan({binx,biny,mphzDV'},[-pi,pi],'circular',true,[0,0,0],1,1,@hsv);axis('xy');
% $$$     subplot(144);imagesc(binx,biny,circ_mean(mphzDVBoot,[],3)');axis('xy');colormap('hsv');colorbar();
    waitforbuttonpress();
    end
end



% $$$ unit = 26;
% $$$ res = spk(unit);
% $$$ figure,
% $$$ %scatter(drz(res,unit==units),ang(res,'spine_middle','spine_upper',2),10,phz(res,spk.map(unit==spk.map(:,1),2)),'filled')
% $$$ scatter(drz(res,unit==units),vxy(res,1),20,phz(res,spk.map(unit==spk.map(:,1),2)),'filled')
% $$$ %scatter(drz(res,unit==units),axy(res,1),20,phz(res,spk.map(unit==spk.map(:,1),2)),'filled')
% $$$ colormap hsv



% pfs drz x sts transition

tags = {'DRZxSTSTRANS_ON','DRZxSTSTRANS_OFF'};
for s = 1:numel(tags)

    samples = round(xyz.sampleRate/2);
    stsper = Trial.stc{'r'};
    stsdiff = [[stsper.data(:,2);size(xyz,1)]-[0;stsper.data(:,1)]]';
    trans = stsper.data(:,s)';
    trans(stsdiff(s:end-2+s)<samples) = [];

    ststrans = xyz.copy('empty');
    ststrans.data = nan([size(xyz,1),1]);
    for transition = trans,
        ststrans.data(transition-samples:transition+samples) = linspace(-1,1,2*samples+1);
    end
    
    fet = ststrans.copy();

    pargs = get_default_args('MjgER2016','MTAApfs','struct');        
    pargs.tag              = tags{s};
    pargs.units            = units;
    pargs.numIter          = 1001;
    pargs.halfsample       = true;
    pargs.overwrite        = true;
    pargs.boundaryLimits   = [-1,1;-1,1];
    pargs.binDims          = [0.1,0.1];
    pargs.SmoothingWeights = [1.5,1.5];
    pargs.autoSaveFlag     = false;

    u = 1;        
    fet.data = [fet.data,drz(:,u)];
    pargs.xyzp   = fet;
    pargs.units  = units(u);
    pargs.states = MTADepoch([],[],ThreshCross(~isnan(ststrans.data),0.5,1),...% SELECT periods where drz
                             fet.sampleRate,fet.sync.copy(),fet.origin,...
                             'TimePeriods','sts',[],'tdrz','d');

    pfsArgs = struct2varargin(pargs);
    pfTemp = MTAApfs(Trial,pfsArgs{:});
    pfTemp.save();
    for u = 2:numel(units);
        pargs.units  = units(u);    
        fet           = ststrans.copy();
        fet.data     = [fet.data,drz(:,u)];
        pargs.xyzp   = fet;   
        pargs.states = MTADepoch([],[],ThreshCross(~isnan(ststrans.data),0.5,1),...
                                 fet.sampleRate,fet.sync.copy(),fet.origin,...
                                 'TimePeriods','sts',[],'tdrz','d');
        pfsArgs = struct2varargin(pargs);
        pfTemp  = MTAApfs(pfTemp,pfsArgs{:});
    end
    pfTemp.save();        
end



[accg,tbins] = autoccg(Trial);
pfd = {MTAApfs(Trial,'tag','DRZxSTSTRANS_ON'),MTAApfs(Trial,'tag','DRZxSTSTRANS_OFF'),MTAApfs(Trial,'tag','HBPITCHxBPITCH_v7')};


nx = numel(pfs)+1+1+1+1+1;
figure,
for u = 1:numel(units),
    maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                                 [pfs,pfd],repmat({units(u)},[1,numel(pfs)+1+1+1]))));
    subplot(1,nx,1);
    bar(tbins,accg(:,units(u)));axis tight;
    subplot(1,nx,2);        
    plot(pft,units(u),'mean',false,...
         [maxRate(pft,units(u),false,'prctile99',0.5)/2-0.2,...
          maxRate(pft,units(u),false,'prctile99',0.5)/2],true,0.5,false,[],@jet);
    title(num2str(units(u)));
    
    for s = 1:numel(pfs),
        subplot(1,nx,s+2);
        plot(pfs{s},units(u),1,false,[0,maxPfsRate],true,0.5,false,[],@jet);
        title(pfs{s}.parameters.states);
    end
    subplot(1,nx,nx-2);
    plot(pfd{1},units(u),1,false,[0,maxPfsRate],false,0.5,false,[],@jet);
    title(pfd{1}.tag);
    subplot(1,nx,nx-1);
    plot(pfd{2},units(u),1,false,[0,maxPfsRate],false,0.5,false,[],@jet);
    title(pfd{2}.tag);
    subplot(1,nx,nx);    
    plot(pfd{3},units(u),'mean',true,[0,maxPfsRate],false,0.5,false,[],@jet);
    title(pfd{3}.tag);
    cax = colorbar;
    set(cax,'Position',cax.Position+[0.05,0,0,0])
    
    waitforbuttonpress();
end
 