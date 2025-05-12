% Analysis of the lateral ego-centric bias of hippocampal predictive coding.
%  - selected trajectories which approach one side of the place
%    field with multiple head-body configurations in order to
%    compare the leftward vs rightward head body angles
%
%  - attempt to find 

EgoProCode2D_load_data();

% >>> Single Trial >>> --------------------------------------------------------
trialIndex = 20;

Trial   = Trials     { trialIndex };
units   = Units      { trialIndex };
meta    = sessionList( trialIndex );
Pft     = placeFieldsNoRear{ trialIndex };

bhvState = 'walk+turn+pause&theta';
txyz = preproc_xyz(Trial,'trb',sampleRate);
headYawCorrection = Trial.meta.correction.headYaw;
headCenterCorrection = Trial.meta.correction.headCenter;
% COMPUTE head basis
hvec = txyz(:,'nose',[1,2])-txyz(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(headYawCorrection),-sin(headYawCorrection);...
                  sin(headYawCorrection), cos(headYawCorrection)],...
                 [2,3],...
                 [1,2]);
hvfl = fet_href_HXY(Trial, sampleRate, [], 'trb', 2.4);
hafl = circshift(hvfl.data,1)-hvfl.data;
hba = fet_hba(Trial, sampleRate); % Head to Body Angle
hav = fet_hbav(Trial, sampleRate);
pch = fet_HB_pitch(Trial, sampleRate);
pch.data = pch.data(:,3);
phz = load_theta_phase(Trial, sampleRate);
pyr = Trial.load('spk', sampleRate, bhvState, units{trialIndex}, 'deburst');
fxyz = filter(txyz.copy(),'ButFilter',3,14,'low');
vxy = vel(filter(txyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
vxyz = vel(filter(txyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2,3]);
lvxy = copy(vxy);
lvxy.data(lvxy.data<=0.0001) = 0.0001;
lvxy.data = log10(lvxy.data);


%%%<<< Bins --------------------------------------------------------------------

bins.hfa.edges = linspace(-pi,pi,17);
bins.hfa.centers = mean(cat(1,bins.hfa.edges(2:end),bins.hfa.edges(1:end-1)));
bins.hfa.count = numel(bins.hfa.centers);
bins.hma.edges = linspace(-pi,pi,17);
bins.hma.centers = mean(cat(1,bins.hma.edges(2:end),bins.hma.edges(1:end-1)));
bins.hma.count = numel(bins.hma.centers);
bins.x.edges = linspace([-500,500,20]);
bins.x.centers = mean(cat(1,bins.x.edges(2:end),bins.x.edges(1:end-1)));
bins.x.count = numel(bins.x.centers);
bins.y = bins.x;

%%%>>> Bins --------------------------------------------------------------------
%                                                       HB
% jg05-0310:  11    29    33    42    49    52    54    60    75    78    80
%bhvState = 'lbhv&theta';
%bhvState = 'walk+turn+pause&theta';
%bhvState = 'walk+turn&theta';
bhvState = 'lbhv&theta';
unit = 25;
rmap = plot(Pft,unit);
pfsRadius = sqrt(sum(rmap(:)>2)*20*20/pi);
[mxr,mxp] = Pft.maxRate(unit);
exy = [bsxfun(@plus,                                            ...
              multiprod(bsxfun(@minus,                          ...
                               mxp,                             ...
                               sq(txyz(:,'hcom',[1,2]))),  ...
                        hvec,2,[2,3]),                          ...
              0)];
sts = [Trial.stc{bhvState}];
sts.resample(txyz);
xts = copy(sts);
xts.cast('TimeSeries');
xts.resample(txyz);
xts.data = logical(xts.data);
headAngle = sq(txyz(:,'nose',[1,2])-txyz(:,'hcom',[1,2]));
headAngle = atan2(headAngle(:,2),headAngle(:,1));
mazeAngle =  sq(txyz(:,'hcom',[1,2]));
mazeAngle = atan2(mazeAngle(:,2),mazeAngle(:,1));     
fieldAngle =  bsxfun(@minus,sq(txyz(:,'hcom',[1,2])),mxp);
fieldAngle = atan2(fieldAngle(:,2),fieldAngle(:,1));
hfa = circ_dist(headAngle,fieldAngle);
pfa = bsxfun(@minus,sq(txyz(:,'hcom',[1,2])),mxp);
pfa = atan2(pfa(:,2),pfa(:,1));
usideleng = sqrt(sum(bsxfun(@minus, sq(txyz(:,'hcom',[1,2])),mxp).^2, 2));
hma = circ_dist(headAngle, mazeAngle);
figure()
subplot(131)
    out = hist2([txyz(xts,'hcom',1),txyz(xts,'hcom',2)], ...
                  bins.x.edges, ...
                  bins.y.edges);
    imagesc(bins.x.centers, bins.x.centers, (out./sampleRate)');
    axis('xy');colormap(gca(),'jet');colorbar();hold('on');
    circle(mxp(1),mxp(2),pfsRadius,'r');
    set(gca(),'ColorScale','log');
subplot(132)
    ind =  xts.data & within_ranges(usideleng,radius);
    hfaVpfaCount = hist2([hfa(ind), pfa(ind)], bins.hfa.edges,bins.hfa.edges);
    imagesc(bins.hfa.centers, bins.hfa.centers, (hfaVpfaCount./sampleRate)')
    axis('xy');colormap(gca(),'jet');colorbar()
    set(Lines(-pi/2,[],'k'),'LineWidth',2);
    set(Lines(pi/2,[],'k'),'LineWidth',2);
    set(gca(),'ColorScale','log');    
subplot(133)
    pfa = fieldAngle;
    ind_L = within_ranges(hba,bins.hba.edges([1,2])) & ind;
    ind_R = within_ranges(hba,bins.hba.edges([3,4])) & ind;
    hfaInds = discretize(hfa(ind_L), bins.hfa.edges);
    pfaInds = discretize(pfa(ind_L), bins.hfa.edges);
    hbaLOcc = accumarray([hfaInds,pfaInds],hba(ind_L),[bins.hfa.count,bins.hfa.count],@numel);
    hfaInds = discretize(hfa(ind_R), bins.hfa.edges);
    pfaInds = discretize(pfa(ind_R), bins.hfa.edges);
    hbaROcc = accumarray([hfaInds,pfaInds],hba(ind_R),[bins.hfa.count,bins.hfa.count],@numel);
    imagescnan({bins.hfa.centers, bins.hfa.centers, (hbaROcc-hbaLOcc)'./sampleRate},[-1,1],'colorbarIsRequired',true,'colorMap',@jet)
    axis('xy');
    set(Lines(-pi/2,[],'k'),'LineWidth',2);
    set( Lines(pi/2,[],'k'),'LineWidth',2);
hfaVpfaOcc = hfaVpfaCount./sampleRate;
hbaRatio = (hbaROcc-hbaLOcc)'./sampleRate;
set(gcf(),'Position',get(gcf(),'Position').*[1,1,3,1]);

%pfaRange =[-pi,-pi*3/4;3*pi/4,pi];
%pfaRange =[-pi/4,pi/4];
%pfaRange =[-pi/2,pi/2];
pfaRange =[-pi,pi];
%pfaRange =[-pi,0];
%pfaRange =[0,pi/2];
%pfaRange =[pi/2,2];

figure
phzI = 1;
phzOS = 1;
pfsRadius = sqrt(sum(rmap(:)>2)*20*20/pi);
for hbaI = 1:bins.hba.count;
    ures = pyr(unit,sts);
    ures = ures(within_ranges(phz(ures),bins.phz.edges([phzI,phzI+phzOS])));
    ind =   within_ranges(pfa,pfaRange) ...
          & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1]));
    ures = ures(ind(ures));
    exyInd = ind & xts.data;
    subplot2(3,1,hbaI,1); hold('on');
        %plot(exy(exyInd,2),exy(exyInd,1),'.','Color',[0.8,0.8,0.8]);
        scatter(exy(exyInd,2),exy(exyInd,1),2,hma(exyInd),'Filled');
        scatter( exy(ures,2), exy(ures,1), 30, hba(ures), 'Filled','MarkerEdgeColor','k');
        %plot( exy(ures,2), exy(ures,1),'.m','MarkerSize',10);
        circle(0,0,pfsRadius,'-k');
        Lines([],0,'k');Lines(0,[],'k');
        colorbar();colormap('hsv');caxis([-pi,pi]);xlim([-250,250]);ylim([-250,250]);daspect([1,1,1]);grid('on');
end
set(gcf(),'Position',get(gcf(),'Position').*[1,0,1,3]);

 Double check the bhv state selection and compare the sector rates
% to the all behavior-rear.

%%%<<< Mean Rate HBA L

radius = [25,150];
%radius = [25,pfsRadius];
hbaI = 1;
phzI = 3;
phzOS = 1;
ures = pyr(unit,sts);
ures = ures(within_ranges(phz(ures),bins.phz.edges([phzI,phzI+phzOS])));
ind =   within_ranges(pfa,pfaRange) ...    
        & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
        & exy(:,2)<-25 ...
        & exy(:,1)>-100; 
ures = ures(ind(ures));
exyInd = ind & xts.data;
rexy = [exy(exyInd,1),exy(exyInd,2)];
sexy = exy(ures,:);
scnt = sum(within_ranges(sqrt(sum(sexy.^2,2)),radius));
pocc = sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))/sampleRate;
scnt/pocc
ures = pyr(unit,sts);
ures = ures(within_ranges( phz(ures), bins.phz.edges([phzI,phzI+phzOS])));
ind =   within_ranges(pfa,pfaRange) ...    
        & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
        & exy(:,2)>25 ...
        & exy(:,1)>-100; 
ures = ures(ind(ures));
exyInd = ind & xts.data;
rexy = [exy(exyInd,1),exy(exyInd,2)];
sexy = exy(ures,:);
scnt = sum(within_ranges(sqrt(sum(sexy.^2,2)),radius));
pocc = sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))/sampleRate;
scnt/pocc

%radius = [0,100];
radius = [25,pfsRadius];
hbaI = 3;
ures = pyr(unit,sts);
ures = ures(within_ranges(phz(ures),bins.phz.edges([phzI,phzI+2])));
ind = within_ranges(hma,[-pi,-pi*3/4; pi*3/4,pi]) ...
      & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
      & exy(:,2)>0; % E to C
ures = ures(ind(ures)); 
exyInd = ind & xts.data;
rexy = [exy(exyInd,1),exy(exyInd,2)];
sexy = exy(ures,:);
sum(within_ranges(sqrt(sum(sexy.^2,2)),radius))/(sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))./sampleRate)
% <<< Single Trial <<< --------------------------------------------------------

% >>> Multi-Session Directional ego field decomposition >>> -------------------
% jg05-201203_10-11  49,152

% jg05-201203_(11)(12)
%              50, 21
%              53, 35
%              63, 31


% jg05-201203_(09)(10)(11)(12)
                                                         
clear('pairs');                                              
pairs(1).tid = 17
pairs(1).tmap = { ...
    'jg05-20120309.cof.all', ...
    'jg05-20120310.cof.all', ...
    'jg05-20120311.cof.all', ...
    'jg05-20120312.cof.all'  ...
};
pairs(1).clm = [0,8];
pairs(1).chn = 8;  pairs(1).unt = [74, 57, 150, 124];
pairs = repmat(pairs,[9,1]);
pairs(2).chn = 8;  pairs(2).unt = [62, 50, 141, 105]; %*%
pairs(3).chn = 8;  pairs(3).unt = [89, 66, 144, 145]; % check 145
pairs(4).chn = 8;  pairs(4).unt = [80, 75, 142, 102]; %*%
pairs(5).chn = 8;  pairs(5).unt = [59, 49, 133, 128]; %*%
pairs(6).chn = 8;  pairs(6).unt = [61, 60, 134, 115];
pairs(7).chn = 8;  pairs(7).unt = [69, 56, 139, 114];
pairs(8).chn = 8;  pairs(8).unt = [70, 74, 149, 103];


bhvState = 'pause+walk+turn&theta';
trialIndex = cellfun(@(F) find_trial_index(Trials,F), pairs(1).tmap);


% jg05-201203_(15)(16)(17)
%               6, 13, 29
%              27, 41, 50 ?
%              33, 30, 54
%              61, 48, 72 ?
%              63, 65, 69
%              77, 61, 63
% bhvState = 'pause+walk+turn&theta';
% trialIndex = [21,22,23];
Pft = cf(@(T)  pfs_2d_theta(T,[],'theta-groom-sit-rear'),  Trials);



tTrial = {};
tunits = {};
tpft = {};
txyz = {};
thyc = {};
tphz = {};
thba = {};
tpyr = {};
for tid = 1:numel(trialIndex)
    tTrial{tid} = Trials{ trialIndex(tid) };
    tunits{tid} = Units { trialIndex(tid) };
    tpft{tid}   = Pft   { trialIndex(tid) };
    txyz{tid}   = Xyz   { trialIndex(tid) };
    thyc{tid} = tTrial{tid}.meta.correction.headYaw;
    headCenterCorrection = tTrial{tid}.meta.correction.headCenter;
% COMPUTE head basis
    hvec = txyz{tid}(:,'nose',[1,2])-txyz{tid}(:,'hcom',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    thvec{tid} = multiprod(hvec,...
                 [cos(thyc{tid}),-sin(thyc{tid});...
                  sin(thyc{tid}), cos(thyc{tid})],...
                 [2,3],...
                 [1,2]);
    %hvfl = fet_href_HXY(Trial, sampleRate, [], 'trb', 2.4);
    %hafl = circshift(hvfl.data,1)-hvfl.data;
    thba{tid} = fet_hba(tTrial{tid}, sampleRate); % Head to Body Angle
    %hav = fet_hbav(tTrial{tid}, sampleRate);
    %pch = fet_HB_pitch(tTrial{tid},sampleRate);
    %pch.data = pch.data(:,3);
    tphz{tid} = load_theta_phase(tTrial{tid}, sampleRate);
    tpyr{tid} = tTrial{tid}.load('spk', sampleRate, bhvState, tunits{tid}, 'deburst');
end

radius = [25,150];
%radius = [25,pfsRadius];
% jg05-201203_(15)(16)(17)
%eunit = {6, 13, 29}; %NN pfaRange =[-pi/4,pi/4];
%eunit = {27, 41, 50}; %N (m, -pi/2:pi/2)
%eunit = {33, 30, 54}; %Y  [-pi,-3*pi/4;3*pi/4,pi];
%eunit = {61, 48, 72};
%eunit = {63, 65, 69}; N (m, 0:pi)
%eunit = {77, 61, 63}; N (m,-pi:0)



exy={};
rmap={};       mazeOcc={};    thfaVpfaOcc={};
sts={};        xts={};
headAngle={};  mazeAngle={};  fieldAngle={};
thfa={};       tpfa={};       thma={};
usideleng={};  ind_L={};      ind_R={};
hfaInds_L=[];    pfaInds_L=[];
hfaInds_R=[];    pfaInds_R=[];
thba_L=[];     thba_R = [];
meanFieldPosition = [];
for tid = 1:numel(trialIndex)
    unit = eunit{tid};
    rmap{tid} = plot(tpft{tid}, unit);
    %pfsRadius = sqrt(sum(rmap(:)>2)*20*20/pi);
    [mxr,mxp] = tpft{tid}.maxRate(unit);
    meanFieldPosition = cat(1,meanFieldPosition,mxp);
end
meanFieldPosition = mean(meanFieldPosition);

for tid = 1:numel(trialIndex)
    unit = eunit{tid};
    rmap{tid} = plot(tpft{tid}, unit);
    %pfsRadius = sqrt(sum(rmap(:)>2)*20*20/pi);
    exy{tid} = [bsxfun(@plus,                                            ...
              multiprod(bsxfun(@minus,                          ...
                               meanFieldPosition,                             ...
                               sq(txyz{tid}(:,'hcom',[1,2]))),  ...
                        thvec{tid},2,[2,3]),                          ...
              0)];
    sts{tid} = [tTrial{tid}.stc{bhvState}];
    sts{tid}.resample(txyz{tid});
    xts{tid} = copy(sts{tid});
    xts{tid}.cast('TimeSeries');
    xts{tid}.resample(txyz{tid});
    xts{tid}.data = logical(xts{tid}.data);
    headAngle{tid} = sq(txyz{tid}(:,'nose',[1,2])-txyz{tid}(:,'hcom',[1,2]));
    headAngle{tid} = atan2(headAngle{tid}(:,2),headAngle{tid}(:,1));
    fieldAngle{tid} =  bsxfun(@minus,sq(txyz{tid}(:,'hcom',[1,2])), meanFieldPosition);
    fieldAngle{tid} = atan2(fieldAngle{tid}(:,2), fieldAngle{tid}(:,1));
    thfa{tid} = circ_dist(headAngle{tid}, fieldAngle{tid});
    tpfa{tid} = bsxfun(@minus,sq(txyz{tid}(:,'hcom',[1,2])), meanFieldPosition);
    tpfa{tid} = atan2(tpfa{tid}(:,2),tpfa{tid}(:,1));
    usideleng{tid} = sqrt(sum(bsxfun(@minus, sq(txyz{tid}(:,'hcom',[1,2])),meanFieldPosition).^2, 2));
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


hbaI = 1;
phzI = 3;
figure
subplot(131);
    imagesc(bins.x.centers, bins.y.centers,log10(sum( cat(3, mazeOcc{:}), 3)'./sampleRate));
    colorbar(); colormap(gca(),'jet');
subplot(132)
    imagesc(bins.hfa.centers, bins.hfa.centers, (sum(cat(3,thfaVpfaOcc{:}),3)./sampleRate)')
    axis('xy');colormap(gca(),'jet');colorbar()
    set(Lines(-pi/2,[],'k'),'LineWidth',2);
    set(Lines(pi/2,[],'k'),'LineWidth',2);
subplot(133)
    imagescnan({bins.hfa.centers, bins.hfa.centers, (hbaROcc-hbaLOcc)'./sampleRate},[-1,1],'colorbarIsRequired',true,'colorMap',@jet)
    axis('xy');
    set(Lines(-pi/2,[],'k'),'LineWidth',2);
    set( Lines(pi/2,[],'k'),'LineWidth',2);
% $$$ hfaVpfaOcc = hfaVpfaCount./sampleRate;
% $$$ hbaRatio = (hbaROcc-hbaLOcc)'./sampleRate;
set(gcf(),'Position',get(gcf(),'Position').*[1,1,3,1]);


%pfaRange =[0,pi/2];
pfaRange =[-pi/3,pi/3];
%pfaRange =[-pi,-pi/2];
pfaRange =[-pi/4,pi/4];
%pfaRange =[pi/2,2];
%pfaRange =[pi/2,pi];
pfaRange = [-pi,-3*pi/4;3*pi/4,pi];
%pfaRange = [-pi,pi];
%pfaRange =[-pi/2,pi/2];
% pfaRange = [-pi,pi/2];
%pfaRange = [-pi/2,0];
%pfaRange = [-pi/2,pi/2];
%pfaRange = [-pi,pi/2];
%pfaRange = [0,pi/2];
%pfaRange = [pi/4,3*pi/2];

ptype = 1
figure
phzI = 3;
phzOS = 1;
%pfsRadius = sqrt(sum(rmap(:)>2)*20*20/pi);
sax = gobjects([1,3]);
sax(1) = subplot2(3,1,1,1); hold('on');sax(2) = subplot2(3,1,2,1); hold('on');sax(3) = subplot2(3,1,3,1); hold('on');
ermap = zeros([5,5,3]);
for hbaI = 1:bins.hba.count
    axes(sax(hbaI));
    poccg = {};
    for tid = 1:numel(trialIndex)
        ind =   within_ranges(tpfa{tid}, pfaRange) ...
                & within_ranges(thba{tid}, bins.hba.edges([hbaI,hbaI+1]));
        exyInd = ind & xts{tid}.data;
        pocci = discretize([exy{tid}(exyInd,2),exy{tid}(exyInd,1)],linspace(-250,250,6));
        pocci = pocci(nniz(pocci),:);
        poccg{tid} = accumarray(pocci,ones([size(pocci,1),1]),[5,5]); ...
        if ptype == 1
            plot(exy{tid}(exyInd,2),exy{tid}(exyInd,1),'.','Color', [0.5,0.5,0.5]);
        end
    end
    soccg = {};
    for tid = 1:numel(trialIndex)
        ind =   within_ranges(tpfa{tid}, pfaRange) ...
                & within_ranges(thba{tid}, bins.hba.edges([hbaI,hbaI+1]));
        ures = tpyr{tid}(eunit{tid}, sts{tid});
        ures = ures(within_ranges(tphz{tid}(ures), bins.phz.edges([phzI,phzI+phzOS])));
        ures = ures(ind(ures));
        socci = discretize([exy{tid}(ures,2),exy{tid}(ures,1)],linspace(-250,250,6));
        socci = socci(nniz(socci),:);
        soccg{tid} = accumarray(socci,ones([size(socci,1),1]),[5,5]);
        %scatter(exy{tid}(ures,2),   exy{tid}(ures,1), 30, thba{tid}(ures), ...
        %        'Filled','MarkerEdgeColor','k');
        if ptype == 1
            plot( exy{tid}(ures,2), exy{tid}(ures,1),'.m', 'MarkerSize',10);
            end
    end
    soccg = sum(cat(3,soccg{:}),3);
    poccg = sum(cat(3,poccg{:}),3);
    poccg(poccg/sampleRate<0.5) = nan;
    ermap(:,:,hbaI) = (soccg./(poccg/sampleRate));
    if ptype == 0
    imagescnan({1:5,1:5,ermap(:,:,hbaI)'},[0,6],'colorbarIsRequired',true);
    axis('xy');
    Lines([],0,'k');Lines(0,[],'k');
    axis('tight');
    else
        xlim([-250,250]); ylim([-250,250]);
    end
    daspect([1,1,1]); grid('on');    
end
set(gcf(),'Position',get(gcf(),'Position').*[1,0,1,3]);



dermap = zeros([2,5,3]);
figure();
for hbaI = 1:bins.hba.count
    dermap(:,:,hbaI) = (ermap([4,5],:,hbaI)-ermap([2,1],:,hbaI));
    subplot2(3,1,hbaI,1);
    imagescnan({1:2, 1:5, dermap(:,:,hbaI)'}, ...
           [-2,2],                                   ...
           'colorbarIsRequired',true,                ...
           'colorMap',@jet)
axis('xy');
end
set(gcf(),'Position',get(gcf(),'Position').*[1,0,1,3]);


(ermap([4,5],2:5,1)-ermap([2,1],2:5,3))
figure();
dm1 = (ermap([4,5],2:5,1)-ermap([2,1],2:5,1));
dm2 = (ermap([4,5],2:5,3)-ermap([2,1],2:5,3));
ind = nniz(dm1(:)) & nniz(dm2(:));
    plot(dm1(ind),...
         dm2(ind) ,...
         '.');
    Lines([],0,'k');
    Lines(0,[],'k');
    mean(dm1(ind))
    mean(dm2(ind))
    



dermap = zeros([3,7,3]);
figure();
for hbaI = 1:bins.hba.count
    dermap(:,:,hbaI) = (ermap(5:7,:,hbaI)-ermap([3,2,1],:,hbaI));
    subplot2(3,1,hbaI,1);
imagescnan({1:3, 1:7, dermap(:,:,hbaI)'}, ...
           [-2,2],                                   ...
           'colorbarIsRequired',true,                ...
           'colorMap',@jet)
axis('xy');
end




figure();
dm1 = (ermap(5:7,3:7,1)-ermap([3,2,1],3:7,1));
dm2 = (ermap(5:7,3:7,3)-ermap([3,2,1],3:7,3));
ind = nniz(dm1(:)) & nniz(dm2(:));
    plot(dm1(ind),...
         dm2(ind) ,...
         '.');
    Lines([],0,'k');
    Lines(0,[],'k');
    mean(dm1(ind))
    mean(dm2(ind))
    
sum(nonzeros(dermap(1:2,4:6,3)))/6;
sum(nonzeros(dermap(1:2,4:6,1)))/6;


%%%<<< Mean Rate HBA L
radius = [25,150];
hbaI  = 1;
phzI  = 3;
phzOS = 1;
sexy  = [];
rexy  = [];
for tid = 1:3
ures = tpyr{tid}(eunit{tid},sts{tid});
ures = ures(within_ranges(tphz{tid}(ures), bins.phz.edges([phzI,phzI+phzOS])));
ind =   within_ranges(tpfa{tid}, pfaRange) ...    
        & within_ranges(thba{tid}, bins.hba.edges([hbaI,hbaI+1])) ...
        & exy{tid}(:,2)<-25 ...
        & exy{tid}(:,1)>-30; 
ures = ures(ind(ures));
exyInd = ind & xts{tid}.data;
rexy = cat(1,rexy,[exy{tid}(exyInd,1),exy{tid}(exyInd,2)]);
sexy = cat(1,sexy,exy{tid}(ures,:));
end
sum(within_ranges(sqrt(sum(sexy.^2,2)),radius))/sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))*sampleRate
sexy = [];
rexy =[];
for tid = 1:3
ures = tpyr{tid}(eunit{tid},sts{tid});
ures = ures(within_ranges(tphz{tid}(ures), bins.phz.edges([phzI,phzI+phzOS])));
ind =   within_ranges(tpfa{tid}, pfaRange) ...    
        & within_ranges(thba{tid}, bins.hba.edges([hbaI,hbaI+1])) ...
        & exy{tid}(:,2)>25 ...
        & exy{tid}(:,1)>-30; 
ures = ures(ind(ures));
exyInd = ind & xts{tid}.data;
rexy = cat(1,rexy,[exy{tid}(exyInd,1),exy{tid}(exyInd,2)]);
sexy = cat(1,sexy,exy{tid}(ures,:));
end
sum(within_ranges(sqrt(sum(sexy.^2,2)),radius))/sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))*sampleRate



ures = pyr(unit,sts);
ures = ures(within_ranges(phz(ures),bins.phz.edges([phzI,phzI+phzOS])));
ind =   within_ranges(pfa,pfaRange) ...    
        & within_ranges(hba,bins.hba.edges([hbaI,hbaI+1])) ...
        & exy(:,2)>25 ...
        & exy(:,1)>-100; 
ures = ures(ind(ures));
exyInd = ind & xts.data;
rexy = [exy(exyInd,1),exy(exyInd,2)];
sexy = exy(ures,:);
sum(within_ranges(sqrt(sum(sexy.^2,2)),radius))/(sum(within_ranges(sqrt(sum(rexy.^2,2)),radius))./sampleRate)
% <<< Multi-Session Directional ego field decomposition <<< -------------------