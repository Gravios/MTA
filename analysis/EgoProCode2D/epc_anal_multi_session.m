% Analysis of the lateral ego-centric bias of hippocampal predictive coding.
%  - selected trajectories which approach one side of the place
%    field with multiple head-body configurations in order to
%    compare the leftward vs rightward head body angles
%
%  - attempt to find 


configure_default_args();

sessionListName = 'BehaviorPlaceCode';
stateCollection = 'msnn_ppsvd_raux';
pitchReferenceTrial = 'Ed05-20140529.ont.all';
pfsState = 'theta-groom-sit-rear';
sampleRate = 250;

clear('states')
states.name = stateCollection;
states.labels =                                                               ...
    {{'theta-groom-sit', 'rear&theta',                                        ...
           'hloc&theta', 'hpause&theta',                                      ...
           'lloc&theta', 'lpause&theta'}};


global MTA_FIGURES_PATH
fig_dir = 'ego-multi-trl-drc';

% >>> Single Trial >>> --------------------------------------------------------

% $$$ trialIndex = 20;
% $$$ 
% $$$ Trial   = Trials     { trialIndex };
% $$$ units   = Units      { trialIndex };
% $$$ meta    = sessionList( trialIndex );
% $$$ Pft     = placeFieldsNoRear{ trialIndex };
% $$$ 
% $$$ bhvState = 'walk+turn+pause&theta';
% $$$ txyz = preproc_xyz(Trial,'trb',sampleRate);
% $$$ headYawCorrection = Trial.meta.correction.headYaw;
% $$$ headCenterCorrection = Trial.meta.correction.headCenter;
% $$$ % COMPUTE head basis
% $$$ hvec = txyz(:,'nose',[1,2])-txyz(:,'hcom',[1,2]);
% $$$ hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
% $$$ hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
% $$$ hvec = multiprod(hvec,...
% $$$                  [cos(headYawCorrection),-sin(headYawCorrection);...
% $$$                   sin(headYawCorrection), cos(headYawCorrection)],...
% $$$                  [2,3],...
% $$$                  [1,2]);
% $$$ hvfl = fet_href_HXY(Trial, sampleRate, [], 'trb', 2.4);
% $$$ hafl = circshift(hvfl.data,1)-hvfl.data;
% $$$ hba = fet_hba(Trial, sampleRate); % Head to Body Angle
% $$$ hav = fet_hbav(Trial, sampleRate);
% $$$ pch = fet_HB_pitch(Trial, sampleRate);
% $$$ pch.data = pch.data(:,3);
% $$$ phz = load_theta_phase(Trial, sampleRate);
% $$$ pyr = Trial.load('spk', sampleRate, bhvState, units{trialIndex}, 'deburst');
% $$$ fxyz = filter(txyz.copy(),'ButFilter',3,14,'low');
% $$$ vxy = vel(filter(txyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
% $$$ vxyz = vel(filter(txyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2,3]);
% $$$ lvxy = copy(vxy);
% $$$ lvxy.data(lvxy.data<=0.0001) = 0.0001;
% $$$ lvxy.data = log10(lvxy.data);


% >>> Bins --------------------------------------------------------------------

% $$$ bins.hfa.edges = linspace(-pi,pi,17);
% $$$ bins.hfa.centers = mean(cat(1,bins.hfa.edges(2:end),bins.hfa.edges(1:end-1)));
% $$$ bins.hfa.count = numel(bins.hfa.centers);
% $$$ bins.hma.edges = linspace(-pi,pi,17);
% $$$ bins.hma.centers = mean(cat(1,bins.hma.edges(2:end),bins.hma.edges(1:end-1)));
% $$$ bins.hma.count = numel(bins.hma.centers);
% $$$ bins.x.edges = linspace([-500,500,20]);
% $$$ bins.x.centers = mean(cat(1,bins.x.edges(2:end),bins.x.edges(1:end-1)));
% $$$ bins.x.count = numel(bins.x.centers);
% $$$ bins.y = bins.x;

% <<< Bins --------------------------------------------------------------------
%                                                       HB
% jg05-0310:  11    29    33    42    49    52    54    60    75    78    80
%bhvState = 'lbhv&theta';
%bhvState = 'walk+turn+pause&theta';
%bhvState = 'walk+turn&theta';
% $$$ bhvState = 'lbhv&theta';
% $$$ unit = 25;
% $$$ rmap = plot(Pft,unit);
% $$$ pfsRadius = sqrt(sum(rmap(:)>2)*20*20/pi);
% $$$ [mxr,mxp] = Pft.maxRate(unit);
% $$$ exy = [bsxfun(@plus,                                            ...
% $$$               multiprod(bsxfun(@minus,                          ...
% $$$                                mxp,                             ...
% $$$                                sq(txyz(:,'hcom',[1,2]))),  ...
% $$$                         hvec,2,[2,3]),                          ...
% $$$               0)];
% $$$ sts = [Trial.stc{bhvState}];
% $$$ sts.resample(txyz);
% $$$ xts = copy(sts);
% $$$ xts.cast('TimeSeries');
% $$$ xts.resample(txyz);
% $$$ xts.data = logical(xts.data);
% $$$ headAngle = sq(txyz(:,'nose',[1,2])-txyz(:,'hcom',[1,2]));
% $$$ headAngle = atan2(headAngle(:,2),headAngle(:,1));
% $$$ mazeAngle =  sq(txyz(:,'hcom',[1,2]));
% $$$ mazeAngle = atan2(mazeAngle(:,2),mazeAngle(:,1));     
% $$$ fieldAngle =  bsxfun(@minus,sq(txyz(:,'hcom',[1,2])),mxp);
% $$$ fieldAngle = atan2(fieldAngle(:,2),fieldAngle(:,1));
% $$$ hfa = circ_dist(headAngle,fieldAngle);
% $$$ pfa = bsxfun(@minus,sq(txyz(:,'hcom',[1,2])),mxp);
% $$$ pfa = atan2(pfa(:,2),pfa(:,1));
% $$$ usideleng = sqrt(sum(bsxfun(@minus, sq(txyz(:,'hcom',[1,2])),mxp).^2, 2));
% $$$ hma = circ_dist(headAngle, mazeAngle);
% $$$ figure()
% $$$ subplot(131)
% $$$     out = hist2([txyz(xts,'hcom',1),txyz(xts,'hcom',2)], ...
% $$$                   bins.x.edges, ...
% $$$                   bins.y.edges);
% $$$     imagesc(bins.x.centers, bins.x.centers, (out./sampleRate)');
% $$$     axis('xy');colormap(gca(),'jet');colorbar();hold('on');
% $$$     circle(mxp(1),mxp(2),pfsRadius,'r');
% $$$     set(gca(),'ColorScale','log');
% $$$ subplot(132)
% $$$     ind =  xts.data & within_ranges(usideleng,radius);
% $$$     hfaVpfaCount = hist2([hfa(ind), pfa(ind)], bins.hfa.edges,bins.hfa.edges);
% $$$     imagesc(bins.hfa.centers, bins.hfa.centers, (hfaVpfaCount./sampleRate)')
% $$$     axis('xy');colormap(gca(),'jet');colorbar()
% $$$     set(Lines(-pi/2,[],'k'),'LineWidth',2);
% $$$     set(Lines(pi/2,[],'k'),'LineWidth',2);
% $$$     set(gca(),'ColorScale','log');    
% $$$ subplot(133)
% $$$     pfa = fieldAngle;
% $$$     ind_L = within_ranges(hba,bins.hba.edges([1,2])) & ind;
% $$$     ind_R = within_ranges(hba,bins.hba.edges([3,4])) & ind;
% $$$     hfaInds = discretize(hfa(ind_L), bins.hfa.edges);
% $$$     pfaInds = discretize(pfa(ind_L), bins.hfa.edges);
% $$$     hbaLOcc = accumarray([hfaInds,pfaInds],hba(ind_L),[bins.hfa.count,bins.hfa.count],@numel);
% $$$     hfaInds = discretize(hfa(ind_R), bins.hfa.edges);
% $$$     pfaInds = discretize(pfa(ind_R), bins.hfa.edges);
% $$$     hbaROcc = accumarray([hfaInds,pfaInds],hba(ind_R),[bins.hfa.count,bins.hfa.count],@numel);
% $$$     imagescnan({bins.hfa.centers, bins.hfa.centers, (hbaROcc-hbaLOcc)'./sampleRate},[-1,1],'colorbarIsRequired',true,'colorMap',@jet)
% $$$     axis('xy');
% $$$     set(Lines(-pi/2,[],'k'),'LineWidth',2);
% $$$     set( Lines(pi/2,[],'k'),'LineWidth',2);
% $$$ hfaVpfaOcc = hfaVpfaCount./sampleRate;
% $$$ hbaRatio = (hbaROcc-hbaLOcc)'./sampleRate;
% $$$ set(gcf(),'Position',get(gcf(),'Position').*[1,1,3,1]);
% $$$ 
% $$$ %pfaRange =[-pi,-pi*3/4;3*pi/4,pi];
% $$$ %pfaRange =[-pi/4,pi/4];
% $$$ %pfaRange =[-pi/2,pi/2];
% $$$ pfaRange =[-pi,pi];
% $$$ %pfaRange =[-pi,0];
% $$$ %pfaRange =[0,pi/2];
% $$$ %pfaRange =[pi/2,2];
% $$$ 
% $$$ figure
% $$$ phzI = 1;
% $$$ phzOS = 1;
% $$$ pfsRadius = sqrt(sum(rmap(:)>2)*20*20/pi);
% $$$ for hbaI = 1:bins.hba.count;
% $$$     ures = pyr(unit,sts);
% $$$     ures = ures(within_ranges(phz(ures),bins.phz.edges([phzI,phzI+phzOS])));
% $$$     ind =   within_ranges(pfa,pfaRange) ...
% $$$           & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1]));
% $$$     ures = ures(ind(ures));
% $$$     exyInd = ind & xts.data;
% $$$     subplot2(3,1,hbaI,1); hold('on');
% $$$         %plot(exy(exyInd,2),exy(exyInd,1),'.','Color',[0.8,0.8,0.8]);
% $$$         scatter(exy(exyInd,2),exy(exyInd,1),2,hma(exyInd),'Filled');
% $$$         scatter( exy(ures,2), exy(ures,1), 30, hba(ures), 'Filled','MarkerEdgeColor','k');
% $$$         %plot( exy(ures,2), exy(ures,1),'.m','MarkerSize',10);
% $$$         circle(0,0,pfsRadius,'-k');
% $$$         Lines([],0,'k');Lines(0,[],'k');
% $$$         colorbar();colormap('hsv');caxis([-pi,pi]);xlim([-250,250]);ylim([-250,250]);daspect([1,1,1]);grid('on');
% $$$ end
% $$$ set(gcf(),'Position',get(gcf(),'Position').*[1,0,1,3]);
% $$$ 
% $$$  Double check the bhv state selection and compare the sector rates
% $$$ % to the all behavior-rear.
% $$$ 
% $$$ 
% $$$ radius = [25,150];
% $$$ %radius = [25,pfsRadius];
% $$$ hbaI = 1;
% $$$ phzI = 3;
% $$$ phzOS = 1;
% $$$ ures = pyr(unit,sts);
% $$$ ures = ures(within_ranges(phz(ures),bins.phz.edges([phzI,phzI+phzOS])));
% $$$ ind =   within_ranges(pfa,pfaRange) ...    
% $$$         & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
% $$$         & exy(:,2)<-25 ...
% $$$         & exy(:,1)>-100; 
% $$$ ures = ures(ind(ures));
% $$$ exyInd = ind & xts.data;
% $$$ rexy = [exy(exyInd,1),exy(exyInd,2)];
% $$$ sexy = exy(ures,:);
% $$$ scnt = sum(within_ranges(sqrt(sum(sexy.^2,2)),radius));
% $$$ pocc = sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))/sampleRate;
% $$$ scnt/pocc
% $$$ ures = pyr(unit,sts);
% $$$ ures = ures(within_ranges( phz(ures), bins.phz.edges([phzI,phzI+phzOS])));
% $$$ ind =   within_ranges(pfa,pfaRange) ...    
% $$$         & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
% $$$         & exy(:,2)>25 ...
% $$$         & exy(:,1)>-100; 
% $$$ ures = ures(ind(ures));
% $$$ exyInd = ind & xts.data;
% $$$ rexy = [exy(exyInd,1),exy(exyInd,2)];
% $$$ sexy = exy(ures,:);
% $$$ scnt = sum(within_ranges(sqrt(sum(sexy.^2,2)),radius));
% $$$ pocc = sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))/sampleRate;
% $$$ scnt/pocc
% $$$ 
% $$$ %radius = [0,100];
% $$$ radius = [25,pfsRadius];
% $$$ hbaI = 3;
% $$$ ures = pyr(unit,sts);
% $$$ ures = ures(within_ranges(phz(ures),bins.phz.edges([phzI,phzI+2])));
% $$$ ind = within_ranges(hma,[-pi,-pi*3/4; pi*3/4,pi]) ...
% $$$       & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
% $$$       & exy(:,2)>0; % E to C
% $$$ ures = ures(ind(ures)); 
% $$$ exyInd = ind & xts.data;
% $$$ rexy = [exy(exyInd,1),exy(exyInd,2)];
% $$$ sexy = exy(ures,:);
% $$$ sum(within_ranges(sqrt(sum(sexy.^2,2)),radius))/(sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))./sampleRate)

% <<< Single Trial <<< --------------------------------------------------------

% >>> Multi-Session Directional ego field decomposition >>> -------------------
% >>> Setup tracked unit sets >>> ---------------------------------------------

clear('pairs');
UnitSet = {};
UnitSet{1}.tid = 17
UnitSet{1}.tmap = { ...
    'jg05-20120309.cof.all'; ...
    'jg05-20120310.cof.all'; ...
    'jg05-20120311.cof.all'; ...
    'jg05-20120312.cof.all';  ...
};
UnitSet{1}.clm = [0,8];
UnitSet{1}.chn = 8;  UnitSet{1}.unt = [74, 57, 150, 124];
pairs = repmat(pairs,[9,1]);
UnitSet{2}.chn = 8;  UnitSet{2}.unt = [62, 50, 141, 105]; %*%
UnitSet{3}.chn = 8;  UnitSet{3}.unt = [89, 66, 144, 145]; % check 145
UnitSet{4}.chn = 8;  UnitSet{4}.unt = [80, 75, 142, 102]; %*%
UnitSet{5}.chn = 8;  UnitSet{5}.unt = [59, 49, 133, 128]; %*%
UnitSet{6}.chn = 8;  UnitSet{6}.unt = [61, 60, 134, 115];
UnitSet{7}.chn = 8;  UnitSet{7}.unt = [69, 56, 139, 114];
UnitSet{8}.chn = 8;  UnitSet{8}.unt = [70, 74, 149, 103];

pid = 2;
UnitSet{1}.tid = 10
UnitSet{1}.tmap = { ...
    'Ed10-20140814.cof.all'; ...
    'Ed10-20140815.cof.all'; ...
};
UnitSet{1}.clm = [0,8];
UnitSet{1}.chn = 4;  UnitSet{1}.unt = [124,94];
UnitSet{1} = repmat(pairs,[2,1]);

clear('pairs');                                              
UnitSet{1}.tid = 0;
UnitSet{1}.tmap = { ...
    'jg05-20120309.cof.all'; ...
    'jg05-20120310.cof.all'; ...
};
UnitSet{1}.clm = [0,8];
UnitSet{1}.chn = 8;  UnitSet{1}.unt = [20, 11];
pairs = repmat(pairs,[2,1]);

% <<< Setup tracked unit sets <<< ---------------------------------------------
% >>> Load Data >>> -----------------------------------------------------------
%for n = 1:numUnitSets
%Pft = cf(@(T)  pfs_2d_theta(T,[],'theta-groom-sit-rear','overwrite',true,'purge',true),  Trials);
%Xyz = cf(@(T)  preproc_xyz(T,'trb',sampleRate),            Trials);
%Bfs = cf(@(T,U)  compute_bhv_ratemaps(T,U), Trials,Units);
% <<< Load Data <<< -----------------------------------------------------------
% >>> Vars >>> ----------------------------------------------------------------
sampleRate = 250;
%bhvState = 'pause+loc&theta';
bhvState = 'theta-sit-groom-rear';
neuron_sets = get_neuron_sets(Trials);
mask_bhv =load('/storage/gravio/data/project/general/analysis/pfsHB_mask.mat');
% <<< Vars <<< ----------------------------------------------------------------
% >>> Clean Neuron Sets >>> ---------------------------------------------------

% Run multiple loops to remove:
%    - low firing rate cells
%    - silent cells
%    - Bad quality cells


bs = [];
gs = [];

hfig = figure();
for sid = bs%1:numel(neuron_sets)
    clf();
    neuron_sets{sid}
    trialIndex = cellfun(@(F) find_trial_index(Trials,F), ...
                         neuron_sets{sid}.Trial);
    for tid = 1:numel(trialIndex)
        subplot(1, numel(trialIndex), tid);
        plot(Pft{trialIndex(tid)}, neuron_sets{sid}.UnitID(tid), [],'text', ...
             [0,10]);
    end
    waitforbuttonpress();
    switch hfig.CurrentCharacter
      case '1'
        neuron_sets{sid}(1,:) = [];
      case '2'
        neuron_sets{sid}(2,:) = [];
      case '3'
        neuron_sets{sid}(3,:) = [];
      case '4'
        neuron_sets{sid}(4,:) = [];
      case '5'
        neuron_sets{sid}(5,:) = [];
      case '6'
        neuron_sets{sid}(6,:) = [];
      case '7'
        neuron_sets{sid}(7,:) = [];
      case 'b',
        bs(end+1) = sid;
      otherwise
        gs(end+1) = sid;
    end
    neuron_sets{sid}
end

bs = unique(bs);
gs = unique(gs);
% Remove the good sets from the bad sets once the unwanted cells are removed
for  g = gs
    bs(bs==g) = [];
end
% <<< Clean Neuron Sets <<< ---------------------------------------------------
for sid = gs
trialIndex = cellfun(@(F) find_trial_index(Trials,F), neuron_sets{sid}.Trial);
% >>> CHECK if all units of set exist

ux = find(~arrayfun(@(N,U) ...
                    ismember(N,U{1}), ...
                    neuron_sets{sid}.UnitID', ...
                    Units(trialIndex)));

if ~isempty(ux)
    for ui = ux(:)'
        tind = find_trial_index(Trials, neuron_sets{sid}.Trial{ui});
        Units{tind} = [ Units{tind}, neuron_sets{sid}.UnitID(ui)];
        compute_bhv_ratemaps(Trials{tind},neuron_sets{sid}.UnitID(ui));
        Bfs{tind} =  compute_bhv_ratemaps(Trials{tind}, Units{tind});
        plot(Bfs{t}, uid, [],'text',[0,6],true,'mazeMask',mask_bhv.mask);
    end
end

% <<< Check if all units of set exist
% >>> ACCUMULATE Trial Data >>> -----------------------------------------------

tTrial = {}; % Trials
tunits = {}; % Units
tpft = {};   % Placefields
txyz = {};   % Position
thba = {};   % Head-Body Angle
tphz = {};   % Theta Phase
tpyr = {};   % Spikes

for tid = 1:numel(trialIndex)
    tTrial{tid} = Trials{ trialIndex(tid) };
    tunits{tid} = Units { trialIndex(tid) };
    tpft{tid}   = Pft   { trialIndex(tid) };
    txyz{tid}   = Xyz   { trialIndex(tid) };

    sts{tid} = [tTrial{tid}.stc{bhvState}];
    sts{tid}.resample(txyz{tid});
    xts{tid} = copy(sts{tid});
    xts{tid}.cast('TimeSeries');
    xts{tid}.resample(txyz{tid});
    xts{tid}.data = logical(xts{tid}.data);

    thyc = tTrial{tid}.meta.correction.headYaw;
    headCenterCorrection = tTrial{tid}.meta.correction.headCenter;
    % COMPUTE head basis
    hvec = txyz{tid}(:,'nose',[1,2])-txyz{tid}(:,'hcom',[1,2]);
    hvec = sq(bsxfun(@rdivide, hvec, sqrt(sum(hvec.^2,3))));
    hvec = cat( 3, hvec, sq(hvec) * [ 0,-1;1,0]);
    thvec{tid} = multiprod(hvec,...
                           [cos(thyc),-sin(thyc);...
                            sin(thyc), cos(thyc)],...
                          [2,3],...
                          [1,2]);
    % LOAD feature - Head-Body Angle
    thba{tid} = fet_hba(tTrial{tid}, sampleRate);
    % LOAD feature - theta phase
    tphz{tid} = load_theta_phase(tTrial{tid}, txyz{tid});
    % LOAD data - spikes
    tpyr{tid} = tTrial{tid}.load('spk', sampleRate, sts{tid}, tunits{tid}, '');
end

% <<< Accumulate Trial Data <<< -----------------------------------------------
% >>> ACCUMULATE Ego-Fields >>> -----------------------------------------------

radius = [25,150];
eunit = neuronSets{sid}.UnitID;
ucnt = numel(eunit)

% >>> Vars >>> ----------------------------------------------------------------

exy={};
rmap={};       mazeOcc={};    thfaVpfaOcc={};
headAngle={};  mazeAngle={};  fieldAngle={};
thfa={};       tpfa={};       thma={};
usideleng={};  ind_L={};      ind_R={};
hfaInds_L=[];    pfaInds_L=[];
hfaInds_R=[];    pfaInds_R=[];
thba_L=[];     thba_R = [];

% <<< Vars < ----------------------------------------------------------------
% >>> GET place field positions >>> -----------------------------------------
mean_field_position = [];
for tid = 1:ucnt
    unit = eunit(tid);
    rmap{tid} = plot(tpft{tid}, unit);
    %pfsRadius = sqrt(sum(rmap(:)>2)*20*20/pi);
    [mxr,mxp] = tpft{tid}.maxRate(unit);
    mean_field_position = cat(1,mean_field_position,mxp);
end
% <<< GET place field positions <<< -----------------------------------------

for tid = 1:ucnt
    unit = eunit(tid);
    rmap{tid} = plot(tpft{tid}, unit);
    field_position =  mean_field_position(tid,:);
    head_position = sq(txyz{tid}(:,'hcom',[1,2]));
    % >>> Transform place field position to the head's frame of reference 
    exy{tid} = [                                                            ...
        bsxfun(                                                             ...
            @plus,                                                          ...
            multiprod(                                                      ...
                bsxfun( @minus, field_position, head_position),             ...
                thvec{tid},2,[2,3]),                                        ...
            0)];
    % <<< Transform place field position to the head's frame of reference <<< -
    % >>> Head angle w.r.t. arena 
    headAngle{tid} = sq(txyz{tid}(:,'nose',[1,2])-txyz{tid}(:,'hcom',[1,2]));
    headAngle{tid} = atan2(headAngle{tid}(:,2),headAngle{tid}(:,1));
    % <<< Head angle w.r.t. arena
    % >>> Field angle w.r.t. arena 
    fieldAngle{tid} =  bsxfun(@minus,sq(txyz{tid}(:,'hcom',[1,2])), mean_field_position(tid,:));
    fieldAngle{tid} = atan2(fieldAngle{tid}(:,2), fieldAngle{tid}(:,1));
    % <<< Field angle w.r.t. room
    % >>> Head to field angle 
    thfa{tid} = circ_dist(headAngle{tid}, fieldAngle{tid});
    % <<<
    %tpfa{tid} = bsxfun(@minus,sq(txyz{tid}(:,'hcom',[1,2])), mean_field_position(tid,:));
    %tpfa{tid} = atan2(tpfa{tid}(:,2),tpfa{tid}(:,1));

    usideleng{tid} = sqrt(sum(bsxfun(@minus, sq(txyz{tid}(:,'hcom',[1,2])),mean_field_position(tid,:)).^2, 2));

    mazeOcc{tid} = hist2([txyz{tid}(xts{tid},'hcom',1),...
                          txyz{tid}(xts{tid},'hcom',2)], ...
                          bins.x.edges, ...
                          bins.y.edges);

    ind =  xts{tid}.data & within_ranges(usideleng{tid},radius);
    thfaVpfaOcc{tid} = hist2([thfa{tid}(ind), ...
                              tpfa{tid}(ind)], ...
                              bins.hfa.edges, ...
                              bins.hfa.edges);

    ind_L{tid} = within_ranges(thba{tid}, bins.hba.edges([1,2])) & ind;
    ind_R{tid} = within_ranges(thba{tid}, bins.hba.edges([3,4])) & ind;

    mazeAngle{tid} =  sq(txyz{tid}(:,'hcom',[1,2]));
    mazeAngle{tid} = atan2(mazeAngle{tid}(:,2),mazeAngle{tid}(:,1));     
    thma{tid} = circ_dist(headAngle{tid}, mazeAngle{tid});

    hfaInds_L = cat(1, hfaInds_L, discretize(thfa{tid}(ind_L{tid}),       bins.hfa.edges));
    pfaInds_L = cat(1, pfaInds_L, discretize(fieldAngle{tid}(ind_L{tid}), bins.hfa.edges));
    thba_L  = cat(1, thba_L, thba{tid}(ind_L{tid}));

    hfaInds_R = cat(1, hfaInds_R, discretize(thfa{tid}(ind_R{tid}),       bins.hfa.edges));
    pfaInds_R = cat(1, pfaInds_R, discretize(fieldAngle{tid}(ind_R{tid}), bins.hfa.edges));
    thba_R  = cat(1, thba_R, thba{tid}(ind_R{tid}));
end
hbaLOcc = accumarray([hfaInds_L, pfaInds_L], ...
                     thba_L,...
                     [bins.hfa.count,bins.hfa.count],@numel);
hbaROcc = accumarray([hfaInds_R,pfaInds_R],...
                     thba_R,...
                     [bins.hfa.count,bins.hfa.count],@numel);

% <<< Accumulate Ego-Fields <<< -----------------------------------------------
% >>> PLOT Occupance by head-field angle >>> ----------------------------------
% $$$ hbaI = 1;
% $$$ phzI = 2;
% $$$ figure
% $$$ subplot(131);
% $$$     imagesc(bins.x.centers, bins.y.centers,log10(sum( cat(3, mazeOcc{:}), 3)'./sampleRate));
% $$$     colorbar(); colormap(gca(),'jet');
% $$$ subplot(132)
% $$$     imagesc(bins.hfa.centers, bins.hfa.centers, (sum(cat(3,thfaVpfaOcc{:}),3)./sampleRate)')
% $$$     axis('xy');colormap(gca(),'jet');colorbar()
% $$$     set(Lines(-pi/2,[],'k'),'LineWidth',2);
% $$$     set(Lines(pi/2,[],'k'),'LineWidth',2);
% $$$ subplot(133)
% $$$     imagescnan({bins.hfa.centers, bins.hfa.centers, (hbaROcc-hbaLOcc)'./sampleRate},[-1,1],'colorbarIsRequired',true,'colorMap',@jet)
% $$$     axis('xy');
% $$$     set(Lines(-pi/2,[],'k'),'LineWidth',2);
% $$$     set( Lines(pi/2,[],'k'),'LineWidth',2);
% $$$ % $$$ hfaVpfaOcc = hfaVpfaCount./sampleRate;
% $$$ % $$$ hbaRatio = (hbaROcc-hbaLOcc)'./sampleRate;
% $$$ set(gcf(),'Position',get(gcf(),'Position').*[1,1,3,1]);
% <<< Plot Occupance by head-field angle <<< ----------------------------------
% >>> SUMMARY figs
% $$$ 
% $$$ dermap = zeros([2,5,3]);
% $$$ figure();
% $$$ for hbaI = 1:bins.hba.count
% $$$     dermap(:,:,hbaI) = (ermap([4,5],:,hbaI)-ermap([2,1],:,hbaI));
% $$$     subplot2(3,1,hbaI,1);
% $$$     imagescnan({1:2, 1:5, dermap(:,:,hbaI)'}, ...
% $$$            [-4,4],                                   ...
% $$$            'colorbarIsRequired',true,                ...
% $$$            'colorMap',@jet)
% $$$ axis('xy');
% $$$ end
% $$$ set(gcf(),'Position',get(gcf(),'Position').*[1,0,1,3]);
% $$$ 
% $$$ 
% $$$ (ermap([4,5],2:5,1)-ermap([2,1],2:5,3))
% $$$ figure();
% $$$ dm1 = (ermap([4,5],2:5,1)-ermap([2,1],2:5,1));
% $$$ dm2 = (ermap([4,5],2:5,3)-ermap([2,1],2:5,3));
% $$$ ind = nniz(dm1(:)) & nniz(dm2(:));
% $$$     plot(dm1(ind),...
% $$$          dm2(ind) ,...
% $$$          '.');
% $$$     Lines([],0,'k');
% $$$     Lines(0,[],'k');
% $$$     mean(dm1(ind))
% $$$     mean(dm2(ind))
% $$$     
% $$$ 
% $$$ 
% $$$ 
% $$$ dermap = zeros([3,7,3]);
% $$$ figure();
% $$$ for hbaI = 1:bins.hba.count
% $$$     dermap(:,:,hbaI) = (ermap(5:7,:,hbaI)-ermap([3,2,1],:,hbaI));
% $$$     subplot2(3,1,hbaI,1);
% $$$ imagescnan({1:3, 1:7, dermap(:,:,hbaI)'}, ...
% $$$            [-2,2],                                   ...
% $$$            'colorbarIsRequired',true,                ...
% $$$            'colorMap',@jet)
% $$$ axis('xy');
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure();
% $$$ dm1 = (ermap(5:7,3:7,1)-ermap([3,2,1],3:7,1));
% $$$ dm2 = (ermap(5:7,3:7,3)-ermap([3,2,1],3:7,3));
% $$$ ind = nniz(dm1(:)) & nniz(dm2(:));
% $$$     plot(dm1(ind),...
% $$$          dm2(ind) ,...
% $$$          '.');
% $$$     Lines([],0,'k');
% $$$     Lines(0,[],'k');
% $$$     mean(dm1(ind))
% $$$     mean(dm2(ind))
% $$$     
% $$$ sum(nonzeros(dermap(1:2,4:6,3)))/6;
% $$$ sum(nonzeros(dermap(1:2,4:6,1)))/6;

% <<< summary figs
% >>> Mean Rate HBA  >>>
% $$$ 
% $$$ radius = [25,200];
% $$$ hbaI  = 1;
% $$$ phzI  = 3;
% $$$ phzOS = 1;
% $$$ sexy  = [];
% $$$ rexy  = [];
% $$$ phzBin = bins.phz.edges([phzI,phzI+phzOS]);
% $$$ hbaBin = bins.hba.edges([hbaI,hbaI+1]);
% $$$ for tid = 1:ucnt
% $$$ ures = tpyr{tid}(eunit(tid),sts{tid});
% $$$ ures = ures(within_ranges(tphz{tid}(ures), phzBin));
% $$$ ind =   within_ranges( tpfa{tid}, pfaRange) ...    
% $$$       & within_ranges( thba{tid}, hbaBin) ...
% $$$       & exy{tid}(:,2) < -25 ...
% $$$       & exy{tid}(:,1) > -100; 
% $$$ ures = ures( ind( ures));
% $$$ exyInd = ind & xts{tid}.data & within_ranges(tphz{tid}, phzBin);
% $$$ rexy = cat(1, rexy, [exy{tid}(exyInd,1), exy{tid}(exyInd,2)] );
% $$$ sexy = cat(1, sexy, exy{tid}(ures,:) );
% $$$ end
% $$$ sum(within_ranges(sqrt(sum(sexy.^2,2)),radius))/sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))*sampleRate
% $$$ sexy = [];
% $$$ rexy =[];
% $$$ for tid = 1:ucnt
% $$$ ures = tpyr{tid}(eunit(tid),sts{tid});
% $$$ ures = ures(within_ranges(tphz{tid}(ures), phzBin));
% $$$ ind =   within_ranges(tpfa{tid}, pfaRange) ...    
% $$$         & within_ranges(thba{tid}, hbaBin) ...
% $$$         & exy{tid}(:,2) > 25 ...
% $$$         & exy{tid}(:,1)> - 100; 
% $$$ ures = ures(ind(ures));
% $$$ exyInd = ind & xts{tid}.data & within_ranges(tphz{tid}, phzBin);
% $$$ rexy = cat(1,rexy,[exy{tid}(exyInd,1),exy{tid}(exyInd,2)]);
% $$$ sexy = cat(1,sexy,exy{tid}(ures,:));
% $$$ end
% $$$ sum(within_ranges(sqrt(sum(sexy.^2,2)),radius))/sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))*sampleRate

% $$$ ures = pyr(unit,sts);
% $$$ ures = ures(within_ranges(phz(ures),bins.phz.edges([phzI,phzI+phzOS])));
% $$$ ind =   within_ranges(pfa,pfaRange) ...    
% $$$         & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
% $$$         & exy(:,2)>25 ...
% $$$         & exy(:,1)>-100; 
% $$$ ures = ures(ind(ures));
% $$$ exyInd = ind & xts.data;
% $$$ rexy = [exy(exyInd,1),exy(exyInd,2)];
% $$$ sexy = exy(ures,:);
% $$$ sum(within_ranges(sqrt(sum(sexy.^2,2)),radius))/(sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))./sampleRate)

% <<< Mean Rate HBA  <<<
end%for sid = gs
% <<< Multi-Session Directional ego field decomposition <<< -------------------

% >>> Multi- session / field-entry-direction for ego fields >>> ---------------

% >>> VARS >>> ----------------------------------------------------------------
% >>>   Behavior States >>> ---------------------------------------------------
bhvState = 'pause&theta';
% <<<   Behavior States <<< ---------------------------------------------------
% >>>   Bins - Head to Field Angle >>> ----------------------------------------
hfa = {}
hfa{end+1} = [-pi,  pi];
hfa{end+1} = [-pi/4,pi/4];
hfa{end+1} = [0,   pi/2];
hfa{end+1} = [pi/4, 3*pi/4];
hfa{end+1} = [pi/2, pi];
hfa{end+1} = [-pi,-3*pi/4;3*pi/4,pi];
hfa{end+1} = [-pi,-pi/2];
hfa{end+1} = [-3*pi/4,-pi/4];
hfa{end+1} = [-pi/2,0];
binCnt = bins.ego.count;
% <<<   Bins - Head to Field Angle <<< ----------------------------------------
% >>>   ego_ratemaps >>> ------------------------------------------------------
ego_ratemaps = zeros(                                                       ...
    [binCnt,                                                                ...
     binCnt,                                                                ...
     numel(neuron_sets),                                                    ...
     bins.phz.count,                                                        ...
     bins.hba.count,                                                        ...
     numel(hfa)]                                                            ...
);
% <<<     ego_ratemaps     <<< ------------------------------------------------
% <<< VARS <<< ----------------------------------------------------------------

for sid = gs([11:60,62:end])
    
if neuron_sets{sid}.CellType == 0
    continue
end

trialIndex = cellfun(@(F) find_trial_index(Trials,F), neuron_sets{sid}.Trial);
% >>> CHECK if all units of set exist >>> -------------------------------------
ux = find(~arrayfun(@(N,U) ...
                    ismember(N,U{1}), ...
                    neuron_sets{sid}.UnitID', ...
                    Units(trialIndex)));

if ~isempty(ux)
    for ui = ux(:)'
        tind = find_trial_index(Trials, neuron_sets{sid}.Trial{ui});
        uind = neuron_sets{sid}.UnitID(ui);
        Units{tind} = [ Units{tind}, uind];
        pfs_2d_theta(Trials{tind},neuron_sets{sid}.UnitID(ui));
        Pft{tind} =  pfs_2d_theta(Trials{tind}, Units{tind});        
        compute_bhv_ratemaps(Trials{tind}, neuron_sets{sid}.UnitID(ui));
        Bfs{tind} =  compute_bhv_ratemaps(Trials{tind}, Units{tind});
        figure,plot(Bfs{tind}, uind, [],'text',[0,6],true,'mazeMask',mask_bhv.mask);
    end
end
% <<< Check if all units of set exist <<< -------------------------------------
% >>> GET place field positions >>> -------------------------------------------
mean_field_position = [];
for tid = trialIndex'
    u = find(tid==trialIndex);
    unit = neuron_sets{sid}.UnitID(u);
    [mxr,mxp] = Pft{tid}.maxRate(unit);
    mean_field_position = cat(1,mean_field_position,mxp);
end
% <<< GET place field positions <<< -------------------------------------------
% >>> COLLATE data over trials >>> --------------------------------------------

% >>>   Allocate Vars
tTrial = {}; % Trials
tunits = {}; % Units
tpft = {};   % Placefields
txyz = {};   % Position
thba = {};   % Head-Body Angle
tphz = {};   % Theta Phase
tpyr = {};   % Spikes
% <<<   Allocate Vars
for tid = 1:numel(trialIndex)
    % >>> Collate Trial Vars 
    tTrial{tid} = Trials{ trialIndex(tid) };
    tunits{tid} = Units { trialIndex(tid) };
    tpft{tid}   = Pft   { trialIndex(tid) };
    txyz{tid}   = Xyz   { trialIndex(tid) };
    % <<< Collate Trial Vars
    % >>> LOAD the periods of the active network and behavioral states 
    sts{tid} = [tTrial{tid}.stc{bhvState}];
    sts{tid}.resample(txyz{tid});
    xts{tid} = copy(sts{tid});
    xts{tid}.cast('TimeSeries');
    xts{tid}.resample(txyz{tid});
    xts{tid}.data = logical(xts{tid}.data);
    % <<< LOAD the periods of the active network and behavioral states
    % >>> COMPUTE head basis 
    thyc = tTrial{tid}.meta.correction.headYaw;
    %headCenterCorrection = tTrial{tid}.meta.correction.headCenter;
    hvec = txyz{tid}(:,'nose',[1,2])-txyz{tid}(:,'hcom',[1,2]);
    hvec = sq(bsxfun(@rdivide, hvec, sqrt(sum(hvec.^2,3))));
    hvec = cat( 3, hvec, sq(hvec) * [ 0,-1;1,0]);
    thvec{tid} = multiprod(hvec,...
                           [cos(thyc),-sin(thyc);...
                            sin(thyc), cos(thyc)],...
                          [2,3],...
                          [1,2]);
    % <<< COMPUTE head basis
    % >>> LOAD feature - Head-Body Angle 
    thba{tid} = fet_hba(tTrial{tid}, sampleRate);
    % <<< LOAD feature - Head-Body Angle
    % >>> LOAD theta phase 
    tphz{tid} = load_theta_phase(tTrial{tid}, txyz{tid});
    % <<< LOAD theta phase
    % >>> LOAD spikes by network and behavioral state (theta) 
    tpyr{tid} = tTrial{tid}.load('spk', sampleRate, sts{tid}, tunits{tid}, '');
    % <<< LOAD spikes by network and behavioral state (theta)
    % >>> Transform place field position to the head's frame of reference 
    field_position =  mean_field_position(tid,:);
    head_position = sq(txyz{tid}(:,'hcom',[1,2]));
    exy{tid} = [                                                            ...
        bsxfun(                                                             ...
            @plus,                                                          ...
            multiprod(                                                      ...
                bsxfun( @minus, field_position, head_position),             ...
                thvec{tid},2,[2,3]),                                        ...
            0)];
    % <<< Transform place field position to the head's frame of reference
    % >>> Head to field angle 

    thfa{tid} = bsxfun(@minus, head_position, field_position);
    thfa{tid} = atan2(thfa{tid}(:,2),thfa{tid}(:,1));

    % <<< Head to field angle 
end

% <<< COLLATE data over trials <<< --------------------------------------------
% >>> COMPUTE Ego Ratemaps >>> ------------------------------------------------
% >>> VAR ego_ratemap    >>> ----------------------------------------------

ego_ratemap = zeros( ...
    [binCnt,         ...
     binCnt,         ...
     1,              ...
     bins.phz.count, ...
     bins.hba.count, ...
     numel(hfa)] ...
);

% <<< VAR ego_ratemap    <<< ----------------------------------------------

filt_sigma = 1.5;

for hfaI = 1:numel(hfa)
    hfaBin = hfa{hfaI};
    for phzI = bins.phz.count
        ermap = zeros([binCnt, binCnt]);
        phzBin = bins.phz.edges([ phzI, phzI+1]);
        for hbaI = 1:bins.hba.count
            hbaBin = bins.hba.edges([ hbaI, hbaI+1]);
            % >>> Position occupancy given Head to Placefield angle >>> -------
            poccg = {};
            for tid = 1:numel(trialIndex)
                ind =   within_ranges(thfa{tid},        hfaBin)             ...
                        & within_ranges(thba{tid},      hbaBin)             ...
                        & within_ranges(tphz{tid}.data, phzBin);
                exyInd = ind & xts{tid}.data;
                pocci = discretize( [exy{tid}( exyInd, 2),                  ...
                                    exy{tid}( exyInd, 1)],                  ...
                                    bins.ego.edges);
                pocci = pocci(nniz(pocci),:);
                poccg{tid} = accumarray( pocci,                             ...
                                         ones([size(pocci,1),1]),           ...
                                         [binCnt, binCnt]);
            end
            % <<< Position occupancy given Head to Placefield angle <<< -------
            % >>> Spike occupancy given Head to Placefield angle >>> ----------

            soccg = {};
            for tid = 1:numel(trialIndex)
                unit = neuron_sets{sid}.UnitID(tid);                
                ind =   within_ranges( thfa{tid}, hfaBin)                 ...
                      & within_ranges( thba{tid}, hbaBin);
                % >>> SELECT spikes by behavioral state

                ures = tpyr{tid}( unit, sts{tid});

                % <<< SELECT spikes by behavioral state
                % >>> SELECT spikes by Theta Phase
                ures = ures( within_ranges( tphz{tid}(ures), phzBin));
                if numel(ures) > 1,
                    ures = ures(ind(ures));
                    % <<< SELECT spikes by Theta Phase
                    socci = discretize([exy{tid}( ures, 2),                     ...
                                        exy{tid}( ures, 1)],                    ...
                                       bins.ego.edges);
                    socci = socci(nniz(socci),:);
                    soccg{tid} = accumarray( socci,                             ...
                                             ones([size(socci,1),1]),           ...
                                             [binCnt, binCnt]);
                else
                    soccg{tid} = zeros([binCnt, binCnt]);
                end
                

            end

            % <<< Spike occupancy given Head to Placefield angle <<< ----------
            % >>> Sum Occupancies over Trials >>> -----------------------------
            soccg = sum(cat(3,soccg{:}),3);
            poccg = sum(cat(3,poccg{:}),3);
            % <<< Sum Occupancies over Trials <<< -----------------------------
            % >>> Smooth head and spike position occupancies >>> --------------

            fsoccg = imgaussfilt( soccg, filt_sigma);
            fpoccg = imgaussfilt( poccg, filt_sigma);
            ermap = fsoccg ./ (fpoccg / sampleRate);
            ermap(fpoccg/sampleRate<0.02) = nan;
            ego_ratemap(:,:,1,phzI,hbaI,hfaI) = ermap;

            % <<< Smooth head and spike position occupancies <<< --------------
        end
    end
end

% <<< Compute Ego Ratemaps <<< ------------------------------------------------
% >>> APPEND ratemaps to ego_ratemaps >>> -------------------------------------
ego_ratemaps(:,:,sid,:,:,:) = ego_ratemap;
% <<< Append ratemaps to ego_ratemaps <<< -------------------------------------
% >>> FIGURE - Ego Fields given HBA and HFA >>> -------------------------------

% >>> SETUP Figure >>> --------------------------------------------------------

global DELETE_CURRENT_AXES;
DELETE_CURRENT_AXES = false;

hfig = figure(666002);
figureFormat      = 'A4';
figureOrientation = 'portrait';
figureUnits       = 'centimeters';
subplotHeight     = 1.4;%cm
subplotWidth      = 1.4;%cm
subplotPadVert    = 0.1;%cm
subplotPadHorz    = 0.1;%cm
                        % SETUP 
setup_figure(hfig,                  ... 
             figureFormat,          ... 
             figureOrientation,     ... 
             figureUnits,           ...
             subplotWidth,          ...
             subplotHeight,         ...
             subplotPadHorz,        ...
             subplotPadVert         ...
             );

% <<< Setup Figure <<< --------------------------------------------------------
% >>> VARS >>> ----------------------------------------------------------------

phzI = 3;
phzBin = bins.phz.edges([ phzI, phzI+1]);
trlCnt = numel(trialIndex);
hbaCnt = bins.hba.count;
hfaCnt = numel(hfa);
ny = hbaCnt*2;
nx = hfaCnt;
clim = [0, max(nonzeros(ego_ratemap))];
xl = [-250,250];
yl = [-250,250];

% <<< VARS <<< ----------------------------------------------------------------

% >>> UPDATE global offsets >>> -----------------------------------------------
[ gyo, gxo ] = deal( 0,  0 );
% <<< GLOBAL OFFSETS <<< ------------------------------------------------------

% >>> Panel - Place and Behavior Field >>> ------------------------------------
% >>>   VARS >>> --------------------------------------------------------------

untCnt = size(neuron_sets{sid},1)
clim = cf(@(B)                                                              ...
          B.maxRate(neuron_sets{sid}.UnitID',                               ...
                    true,                                                   ...
                    1,                                                      ...
                    'mask', mask_bhv.mask),                                 ...
          Bfs(trialIndex));
clim = [0,max([clim{1}(1),clim{2}(2)])];

% <<<   VARS <<< --------------------------------------------------------------
for tid = trialIndex';    
    u = find(tid==trialIndex);
    uid = neuron_sets{sid}.UnitID(u);
    % >>> Place fields of each session >>> ------------------------------------
    % >>>    Axes Coordinates >>> ---------------------------------------------
    [yindex, yoffset] = deal( u+1,  (-u+1)/2 );
    [xindex, xoffset] = deal( 1,  0 );
    [   gyo, gxo    ] = deal( gyo,  gxo );
    [yscale, xscale ] = deal( 1,  1 );
    % <<<    Axes Coordinates    <<< ------------------------------------------
    % >>>    Axes Setup >>> ---------------------------------------------------
    sax = setup_axes(                                                       ...
        hfig,                                                               ...
        yindex, yoffset,                                                    ...
        xindex, xoffset,                                                    ...
        gyo,    gxo,                                                        ...
        yscale, xscale);
    % <<<    Axes Setup    <<< ------------------------------------------------
    % >>>    Plot Theta Place Field >>> ---------------------------------------

    plot(Pft{tid}, uid, [],'text',clim);

    % <<<    Plot Theta Place Field    <<< ------------------------------------
    % >>>    Annotations and Formating >>> ------------------------------------
    daspect(sax, [1,1,1]);
    sax.XTickLabel = {};
    sax.YTickLabel = {};
    if mod(u+1,2)==0,
        ylabel({neuron_sets{sid}.Trial{u},' '});
    else
        ylabel(neuron_sets{sid}.Trial(u));
    end
    % <<<    Annotations and Formating   <<< -----------------------------------
    % <<< Place fields of each session <<< ------------------------------------
    % >>> Behavior fields of each session >>> ---------------------------------
    % >>>    Axes Coordinates >>> ---------------------------------------------
    [yindex, yoffset] = deal( u+1, (-u+1)/2 );
    [xindex, xoffset] = deal( 2,  0 );
    [   gyo, gxo    ] = deal( gyo,  gxo );
    [yscale, xscale ] = deal( 1,  1 );
    % <<<    Axes Coordinates <<< ---------------------------------------------
    % >>>    Axes Setup >>> ---------------------------------------------------
    sax = setup_axes(                                                       ...
        hfig,                                                               ...
        yindex, yoffset,                                                    ...
        xindex, xoffset,                                                    ...
        gyo,    gxo,                                                        ...
        yscale, xscale);
    % <<<    Axes Setup    <<< ------------------------------------------------
    % >>>    Plot Behavior Field >>> ------------------------------------------

    plot(Bfs{tid}, uid, [],'text',clim,true,'mazeMask',mask_bhv.mask);

    % <<<    Plot Behavior Field    <<< ---------------------------------------    
    % >>>    Annotations and Formating >>> ------------------------------------
    daspect(sax, [1,1,1]);
    sax.XTickLabel = {};
    sax.YTickLabel = {};
    if tid == 1
        colorbar(sax)
        caxis(clim);
    end
    % <<<    Annotations and Formating   <<< -----------------------------------
    % <<< Behavior fields of each session <<< ---------------------------------
end
% <<< Panel - Place and Behavior Field <<< ------------------------------------

% >>> UPDATE global offsets >>> -----------------------------------------------
[ gyo, gxo ] = deal( 0,  5 );
% <<< GLOBAL OFFSETS <<< ------------------------------------------------------

% >>> PANEL - Subject above columns >>> ---------------------------------------
% >>>   LOAD the Rat Model >>> ------------------------------------------------
rat = load_patch_model('rat');
% <<<   LOAD the Rat Model <<< ------------------------------------------------
for hbaI = 1:bins.hba.count
    % >>>   Axes Coordinates >>> --------------------------------------------------
        [yindex, yoffset] = deal(    1,  0 );
        [xindex, xoffset] = deal( hbaI,0.5 );
        [   gyo, gxo    ] = deal(  gyo,  gxo );
        [yscale, xscale ] = deal(    1,  1 );
    % <<<   Axes Coordinates   <<< ----------------------------------------
    % >>>   Axes Setup >>> --------------------------------------------------------
        sax = setup_axes(                                                   ...
            hfig,                                                           ...
            yindex, yoffset,                                                ...
            xindex, xoffset,                                                ...
            gyo,    gxo,                                                    ...
            yscale, xscale);
        % <<<   Axes Setup   <<< ----------------------------------------------
    % >>>   PLOT Subject >>> ------------------------------------------------------
    subject = struct(rat);
    subject = update_subject_patch(subject, 'head',                         ...
                                   [], false,                               ...
                                   bins.hba.edges,                          ...
                                   bins.hba.centers);
    subject = update_subject_patch(subject, 'body',                         ...
                                   bins.hba.count+1-hbaI,  true,            ...
                                   bins.hba.edges,                          ...
                                   bins.hba.centers);
    patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
    patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
    patch(subject.body.overlay.vert{:}, [0.75,0.50,0.50],                   ...
          'FaceAlpha', 0.3,                                                 ...
          'EdgeColor',bins.hba.color(hbaI,:));
    line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
    line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);

% <<<   PLOT Subject <<< ------------------------------------------------------
    % >>>   Annotations and Formating >>> ---------------------------------
    xlim(sax, [-150, 150]);
    ylim(sax, [-150, 150]);
    Lines([],0,'k');
    Lines(0,[],'k');
    daspect(sax, [1,1,1]);
    sax.XTickLabel = {};
    sax.YTickLabel = {};
    % <<<   Annotations and Formating   <<< -------------------------------
end
% <<< PANEL - Subject above columns <<< ---------------------------------------
rmap = reshape(ego_ratemaps(:,:,sid,:,:,:),[],1);
rmap(~nniz(rmap)) = [];
elim = [0,prctile(nonzeros(rmap),99)];
if elim(2)==nan
    continue
end


for hfaI = 1:hfaCnt
    hfaBin = hfa{hfaI};
    % >>> PANEL - Direction of place cell entry >>> ---------------------------
    % >>>   Axes Coordinates >>> ----------------------------------------------
        [yindex, yoffset] = deal( hfaI+1,  0 );
        [xindex, xoffset] = deal( 1,  -1.4 );
        [   gyo, gxo    ] = deal( gyo,  gxo );
        [yscale, xscale ] = deal( 1,  1 );
    % <<<   Axes Coordinates <<< ----------------------------------------------
    % >>>   Axes Setup >>> ----------------------------------------------------
    sax = setup_axes(                                                       ...
        hfig,                                                               ...
        yindex, yoffset,                                                    ...
        xindex, xoffset,                                                    ...
        gyo,    gxo,                                                        ...
        yscale, xscale);
    % <<< Axes Setup <<< ------------------------------------------------------
    % >>>   PLOT the circle and arc >>> ---------------------------------------
    circle(0,0,1,'--k');
    theta = [];
    for hb = hfaBin'
        theta = [linspace([hb', 100]), theta];
    end
    patch([0 cos(theta) 0], ...
          [0 sin(theta) 0], ...
          'b',...
          'FaceAlpha', 0.5);
    % <<<   PLOT the circle and arc <<< ---------------------------------------
    % >>>   Annotations and Formating >>> -------------------------------------
    xlim(sax, [-1,1]);
    ylim(sax, [-1,1]);
    daspect(sax, [1,1,1]);
    sax.XTickLabel = {};
    sax.YTickLabel = {};
    % <<<   Annotations and Formating   <<< -------------------------------    
    % <<< PANEL - Direction of place cell entry <<< ---------------------------
    for hbaI = 1:hbaCnt
        hbaBin = bins.hba.edges([ hbaI, hbaI+1]);
        
        % >>> PANEL - EGO RATEMAP given head to maze angle >>> ----------------

        % >>>   Axes Coordinates >>> ------------------------------------------
        [yindex, yoffset] = deal( hfaI+1,  0 );
        [xindex, xoffset] = deal( hbaI,  0.5 );
        [   gyo, gxo    ] = deal( gyo,  gxo );
        [yscale, xscale ] = deal( 1,  1 );
        % <<<   Axes Coordinates   <<< ----------------------------------------
        % >>>   Axes Setup >>> ------------------------------------------------
        sax = setup_axes(                                                   ...
            hfig,                                                           ...
            yindex, yoffset,                                                ...
            xindex, xoffset,                                                ...
            gyo,    gxo,                                                    ...
            yscale, xscale);
        % <<<   Axes Setup   <<< ----------------------------------------------
        % >>>   Plot Ego Ratemap >>> ------------------------------------------
        imagescnan(                                                         ...
            {                                                               ...
                bins.ego.centers,                                           ...
                bins.ego.centers,                                           ...
                ego_ratemaps( :, :, sid, phzI, hbaI, hfaI)'                  ...
            },                                                              ...
            'colorLimits', elim,                                            ...
            'colorMap',    @jet                                             ...
            );        
        % <<<   Plot Ego Ratemap   <<< ----------------------------------------
        % >>>   Annotations and Formating >>> ---------------------------------
        Lines([],0,'w');
        Lines(0,[],'w');
        axis(sax, 'xy');
        xlim(sax, xl);
        ylim(sax, yl);
        daspect(sax, [1,1,1]);
        sax.XTickLabel = {};
        sax.YTickLabel = {};
        % <<<   Annotations and Formating   <<< -------------------------------

        % <<< PANEL - EGO RATEMAP given head to maze angle <<< ----------------
        
        % >>> PANEL - EGO SPIKES and OCCUPANCY scatter plot >>> ---------------

        % >>>   Axes Coordinates >>> ------------------------------------------
        [yindex, yoffset] = deal( hfaI+1, 0 );
        [xindex, xoffset] = deal( hbaI+3, 1 );
        [   gyo, gxo    ] = deal( gyo,  gxo );
        [yscale, xscale ] = deal( 1,  1 );
        % <<<   Axes Coordinates   <<< ----------------------------------------
        % >>>   Axes Setup >>> ------------------------------------------------
        sax = setup_axes(                                                   ...
            hfig,                                                           ...
            yindex, yoffset,                                                ...
            xindex, xoffset,                                                ...
            gyo,    gxo,                                                    ...
            yscale, xscale);
        % <<<   Axes Setup   <<< ----------------------------------------------
        % >>>   Plot Ego Occupancy >>> ----------------------------------------
        hold('on');
        for tid = 1:trlCnt
            ind =   within_ranges(thfa{tid},      hfaBin)                   ...
                  & within_ranges(thba{tid},      hbaBin)                   ...
                  & within_ranges(tphz{tid}.data, phzBin);
            exyInd = ind & xts{tid}.data;
            plot( exy{tid}(exyInd,2),                                       ...
                  exy{tid}(exyInd,1),                                       ...
                  '.',                                                      ...
                  'Color', [0.5,0.5,0.5]                                    ...
                  );
        end
        % <<<   Plot Ego Occupancy   <<< --------------------------------------
        % >>>   Plot Ego Spikes >>> -------------------------------------------

        for tid = 1:trlCnt
            unit = neuron_sets{sid}.UnitID(tid);
            ind =   within_ranges( thfa{tid}, hfaBin)                       ...
                  & within_ranges( thba{tid}, hbaBin);
            % SELECT spikes by state
            ures = tpyr{tid}( unit, sts{tid});
            % SELECT spikes by Theta Phase
            ures = ures( within_ranges( tphz{tid}(ures), phzBin));
            ures = ures(ind(ures));
            disp(numel(ures))
            plot( exy{tid}( ures, 2),                                       ...
                  exy{tid}( ures, 1),                                       ...
                  '.m',                                                     ...
                  'MarkerSize', 5                                           ...
            );
        end

        % <<<   Plot Ego Spikes   <<< -----------------------------------------
        % >>>   Annotations and Formating >>> ---------------------------------
        Lines([],0,'k');
        Lines(0,[],'k');
        xlim(sax, xl);
        ylim(sax, yl);
        daspect(sax, [1,1,1]);
        sax.XTickLabel = {};
        sax.YTickLabel = {};
        % <<<   Annotations and Formating   <<< -------------------------------

        % <<< PANEL - Ego spikes and occupancy scatter plot <<< ---------------
    end
end

% <<< FIGURE - Ego Fields given HBA and HFA <<< -------------------------------
% >>> SAVE Figure as pdf >>> --------------------------------------------------
fig_name = ['ego-Multi-Trl-Drc',                                            ...
            '_neuron_set-', num2str(sid,'%04.f'), '_',                      ...
            neuron_sets{sid}.Trial{1}];
create_directory(fullfile( MTA_FIGURES_PATH, fig_dir));
fig_path = fullfile( MTA_FIGURES_PATH, fig_dir, [fig_name, '.pdf']);
print(hfig, '-dpdf', fig_path);
% <<< SAVE Figure as pdf <<< --------------------------------------------------
end%for sid = gs

% <<< Multi- session / field-entry-direction for ego fields <<< ---------------
