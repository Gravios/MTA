
figure,
hold('on');
plot(pfsCenterHR(pargs.states,2),pfsCenterHR(pargs.states,1),'.');

xlim([-350,350]);
ylim([-350,350]);
daspect([1,1,1]);

tind = 17;
nbins = 8;



clear('pairs');                                              
pairs(1).tid = 17
pairs(1).clm = [0,8];
pairs(1).chn = 8;  pairs(1).unt = [74, 57, 150, 124];
pairs = repmat(pairs,[8,1]);
pairs(2).chn = 8;  pairs(2).unt = [62, 50, 141, 105]; %*%
pairs(3).chn = 8;  pairs(3).unt = [89, 66, 144, 145]; % check 145 MAYBE add
pairs(4).chn = 8;  pairs(4).unt = [80, 75, 142, 102]; %*%
pairs(5).chn = 8;  pairs(5).unt = [59, 49, 133, 128]; %*%
pairs(6).chn = 8;  pairs(6).unt = [61, 60, 134, 115]; % ADD group to ego units
pairs(7).chn = 8;  pairs(7).unt = [69, 56, 139, 114];
pairs(8).chn = 8;  pairs(8).unt = [70, 74, 149, 103];

unit = 89;
figure,plot(pft,unit,1,'text');


Trial = Trials{tind};
units = Units{tind};
pft = Pfnr{tind};
xyz = Xyz{tind};
spk = Spk{tind};
unit = units(1);

figure,plot(Pfnr{tind},unit,1,'text');
% COMPUTE anglular distance between the head and body
    headBodyAng = [xyz(:,'spine_upper',[1,2])-xyz(:,'bcom',[1,2]),...
                   xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2])];
    headBodyAng = sq(bsxfun(@rdivide,headBodyAng,sqrt(sum(headBodyAng.^2,3))));
    headBodyAng = cart2pol(headBodyAng(:,:,1),headBodyAng(:,:,2));
    headBodyAng = circ_dist(headBodyAng(:,2),headBodyAng(:,1));
    headBodyAng = MTADfet.encapsulate(Trial,...
                                      -(headBodyAng+Trial.meta.correction.headBody),...
                                      sampleRate,...
                                      'hba','hba','h');
% TRANSFORM Local Field Potential -> theta phase
    phz = load_theta_phase(Trial,xyz);
% COMPUTE head basis
    hvec = xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    hvec = multiprod(hvec,...
                     [cos(headYawCorrection),-sin(headYawCorrection);...
                      sin(headYawCorrection),cos(headYawCorrection)],...
                     [2,3],...
                     [1,2]);
    thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);
figure
clf();
for unit = units;
for hbaI = 1:3;
for phzI = 1:3;
[mxr,mxp] = pft.maxRate(unit);
pfsCenterHR = MTADfet.encapsulate(Trial,                                           ...
                                  bsxfun(                                         ...
                                      @plus,                                       ...
                                      multiprod(bsxfun(@minus,                     ...
                                                  mxp,                             ...
                                                  sq(xyz(:,'hcom',[1,2]))),        ...
                                                hvec,2,[2,3]),                        ...
                                      headCenterCorrection),                       ...
                                  sampleRate,                                      ...
                                  'egocentric_placefield',                         ...
                                  'egopfs',                                        ...
                                  'p'                                              ...
                                  );
pargs.spk = copy(spk);
pargs.states = copy(thetaState);
pargs.states.label = ['thetaPhz_',num2str(phzI)];
pargs.states.data(   phz(:,1) < binPhzs(phzI)             ...
                     | phz(:,1) >= binPhzs(phzI+1)          ... 
                     | headBodyAng.data < binHbas(hbaI)       ...
                     | headBodyAng.data >= binHbas(hbaI+1)) = 0;
cast(pargs.states,'TimePeriods');
resInd = WithinRanges( pargs.spk.res, pargs.states.data);
pargs.spk.res = pargs.spk.res(resInd);
pargs.spk.clu = pargs.spk.clu(resInd);
ind = pargs.spk(unit);
occ = hist2([pfsCenterHR(pargs.states,2),pfsCenterHR(pargs.states,1)],...
            linspace(-300,300,nbins),linspace(-300,300,nbins));
scc = hist2([pfsCenterHR(ind,2),pfsCenterHR(ind,1)], ...
            linspace(-300,300,nbins),linspace(-300,300,nbins));
subplot2(3,3,4-phzI,hbaI);
nind = occ<60;
scc(nind) = nan;
imagescnan({linspace(-300,300,nbins-1),linspace(-300,300,nbins-1),(scc./occ)'},'colorMap',@jet);
axis('xy');
colormap('jet');
caxis([0,0.08])
end
end
waitforbuttonpress();
end




% >>> Bins >>> ----------------------------------------------------------------
bins.x.name = 'x';
bins.x.description = 'x head position';
bins.x.edges = linspace(-500,500,16);
bins.x.centers = (bins.x.edges(1:end-1)+bins.x.edges(2:end))./2;
bins.x.count = numel(bins.x.centers);
bins.x.color = 'k';
bins.x.key = '';
bins.x.label = {};

bins.y.name = 'y';
bins.y.description = 'y head position';
bins.y.edges = linspace(-500,500,16);
bins.y.centers = (bins.y.edges(1:end-1)+bins.y.edges(2:end))./2;
bins.y.count = numel(bins.y.centers);
bins.y.color = 'k';
bins.y.key = '';
bins.y.label = {};

hdir = { ...    
    [-pi,  -pi/2], ...
    [-3*pi/4,  -pi/4], ...
    [-pi,  -pi/2], ...
    [-pi/2,     0], ...
    [-pi/4,  pi/4], ...        
    [0, pi/2], ...
    [pi/4,  3*pi/4], ...    
    [pi/2,   pi], ...
   [-3*pi/4,-pi, ; 3*pi/4, pi] , ...    
};

bins.hdr.edges = linspace(-pi, pi, 7);
bins.hdr.centers = (bins.hdr.edges(1:end-1)+bins.hdr.edges(2:end))./2;
bins.hdr.count = numel(bins.hdr.centers);
% <<< Bins <<< ----------------------------------------------------------------
sampleRate = 250;

pft = pfs_2d_theta(Trial);
figure,
plot(pft,41,1,'text');

xcc = zeros([bins.x.count,bins.y.count,bins.hdr.count]);
sxy = []; sph = []; shd = [];
%for unit = [68,41;19,20];
%for unit = [69,43;19,20];
for unit = [70,48;19,20];
Trial = Trials{unit(2)};
spk = Trial.load('spk',sampleRate,'theta',[41,43,48]);
stc = copy(Trial.stc);
xyz = preproc_xyz(Trial,'trb',sampleRate);
phz = load_theta_phase(Trial,xyz);
hir = atan2(xyz(:,'nose',2)-xyz(:,'hcom',2),xyz(:,'nose',1)-xyz(:,'hcom',1));
sper = [stc{'hloc+lloc+lpause&theta',sampleRate}];
uint = unit(1);
mres = spk(uint);
mres = mres(within_ranges(mres,sper.data));
mres2 = spk(uint);
mres2 = mres(within_ranges(mres,sper.data));

sxy = [sxy; sq(xyz(mres,'hcom',1:2))];
sph = [sph; phz(mres)];
shd = [shd; hir(mres)];

hir(hir==0) = nan;
hInd = discretize(SelectPeriods(hir,sper.data,'c',1,0), bins.hdr.edges);
xInd = discretize(xyz(sper,'hcom',1), bins.x.edges);
yInd = discretize(xyz(sper,'hcom',2), bins.y.edges);
xcc = xcc + accumarray([xInd,yInd,hInd],1, [bins.x.count,bins.y.count,bins.hdr.count],@sum)./sampleRate;
end





figure,
subplot(1,4,1);
imagescnan({bins.x.centers,...
            bins.x.centers,...
            circ_dist(ccm,2.5)'},...
           [-pi,pi],'circular',true,'colorMap',@hsv);
subplot(1,4,2);
imagescnan({bins.x.centers,...
            bins.x.centers,...
            ccv'},...
           [0,1],'linear',true,'colorMap',@jet);
subplot(1,4,3);
imagescnan({bins.x.centers,...
            bins.x.centers,...
            ccs'},...
           [0,pi/2],'linear',true,'colorMap',@jet);
subplot(1,4,4);
imagescnan({bins.x.centers,...
            bins.x.centers,...
            ccr'},...
           [0,0.5],'linear',true,'colorMap',@jet);

figure
for hI = 1:8
subplot(2,4,hI);
imagescnan({bins.x.centers,...
            bins.y.centers,...
            circ_dist(acm(:,:,hI),ccm)'},...
           [-pi,pi],'circular',true,'colorMap',@hsv);
axis('xy');
end


figure
for hI = 1:8
subplot(2,4,hI);
imagescnan({bins.x.centers,...
            bins.y.centers,...
            acs(:,:,hI)'},...
           [0,pi/2],'linear',true,'colorMap',@jet);
axis('xy');
end


figure
for hI = 1:8
subplot(2,4,hI);
imagescnan({bins.x.centers,...
            bins.y.centers,...
            acv(:,:,hI)'},...
           [0,1],'circular',true,'colorMap',@hsv);
axis('xy');
end

figure
for hI = 1:8
subplot(2,2,hI);
imagescnan({bins.x.centers,...
            bins.y.centers,...
            circ_dist(acm(:,:,hI),ccm)'},...
           [-pi/2,pi/2],'circular',true,'colorMap',@hsv);
axis('xy');
end


figure
for hI = 1:4
subplot(2,2,hI);
imagescnan({bins.x.centers,...
            bins.x.centers,...
            log10(acc(:,:,hI))'},...
           [],'linear',true,'colorMap',@jet);
axis('xy');
end




Pft = cf(@(T) pfs_2d_theta(T), Trials);
Pfs = cf(@(T) pfs_2d_states(T,'all',[],states.labels{1}), Trials);
figure,
u = 1; 
for unit = [10,11,12;...
            29,26,53;...
            11,25,41;];
subplot2(2,3,1,u);
    plot(Pft{unit(1)},unit(2),1,'text');
subplot2(2,3,2,u);
    plot(Pft{unit(1)},unit(3),1,'text');
u = u + 1;
end

figure,
u = 1; 
for unit = [10,11,12;...
            29,26,53;...
            11,25,41;];
for sts = 1:6
    subplot2(6,3,sts,u);
    plot(Pfs{unit(1)}{sts},unit(3),1,'text');
end

u = u + 1;
end





% >>> Accumulate vars >>> -----------------------------------------------------
sxy = []; sph = []; shd = [];
sft = 0;  ares = []; scl = [];
xcc = zeros([bins.x.count,bins.y.count,bins.hdr.count]);
%for unit = [68,41;19,20];
%for unit = [69,43;19,20];
% $$$ for unit = [19,20;...
% $$$             68,41;...
% $$$             70,48;];
for unit = [10,11,12;...
            29,26,53;...
            11,25,41;];
Trial = Trials{unit(1)};
spk = Trial.load('spk',sampleRate,'theta',unit(2:3));
stc = copy(Trial.stc);
xyz = preproc_xyz(Trial,'trb',sampleRate);
phz = load_theta_phase(Trial,xyz);
hir = atan2(xyz(:,'nose',2)-xyz(:,'hcom',2),xyz(:,'nose',1)-xyz(:,'hcom',1));
%sper = [stc{'hloc+lloc+lpause&theta',sampleRate}];
sper = [stc{'loc&theta',sampleRate}];
mres1 = spk(unit(2));
mres1 = mres1(within_ranges(mres1, sper.data));
mres2 = spk(unit(3));
mres2 = mres2(within_ranges(mres2, sper.data));
mres = [mres1; ...
        mres2];
sxy = [sxy; sq(xyz(mres,'hcom',1:2))];
sph = [sph; phz(mres)];
shd = [shd; hir(mres)];
scl = [scl; ones(size(mres1)); 2*ones(size(mres2))];
ares = [ ares ; mres + sft];
sft = sft+size(phz,1);
hir(hir==0) = nan;
hInd = discretize(SelectPeriods(hir,sper.data,'c',1,0), bins.hdr.edges);
xInd = discretize(xyz(sper,'hcom',1), bins.x.edges);
yInd = discretize(xyz(sper,'hcom',2), bins.y.edges);
xcc = accumarray([xInd,yInd,hInd],1, [bins.x.count,bins.y.count,bins.hdr.count],@sum)./sampleRate;

xInd = discretize(xyz(mres1,'hcom',1), bins.x.edges);
yInd = discretize(xyz(mres1,'hcom',2), bins.y.edges);
hInd = discretize(SelectPeriods(hir,sper.data,'c',1,0), bins.hdr.edges);
mres1 = SelectPeriods(mres1,sper.data,'d',1,1);
hInd = hInd(mres1);
scc1 = accumarray([xInd,yInd,hInd],1, [bins.x.count,bins.y.count,bins.hdr.count],@sum);

xInd = discretize(xyz(mres2,'hcom',1), bins.x.edges);
yInd = discretize(xyz(mres2,'hcom',2), bins.y.edges);
hInd = discretize(SelectPeriods(hir,sper.data,'c',1,0), bins.hdr.edges);
mres2 = SelectPeriods(mres2,sper.data,'d',1,1);
hInd = hInd(mres2);
scc2 = accumarray([xInd,yInd,hInd],1, [bins.x.count,bins.y.count,bins.hdr.count],@sum);



[W,H] = meshgrid(bins.x.centers, bins.x.centers);           
mazeMask = sqrt(W.^2 + H.^2) <480;

figure,
for hI = 1:bins.hdr.count
    subplot2(2,6,1,hI);
    imagesc(imgaussfilt(scc1(:,:,hI),0.75)./imgaussfilt(xcc(:,:,hI),0.75).*mazeMask)    
    colormap('jet');
    caxis([0,25]);
    axis('xy')
    subplot2(2,6,2,hI);
    imagesc(imgaussfilt(scc2(:,:,hI),0.75)./imgaussfilt(xcc(:,:,hI),0.75).*mazeMask)
    colormap('jet');
    caxis([0,15]);
    axis('xy')
end
end

% <<< Accumulate vars <<< -----------------------------------------------------

% >>> interneuron theta circ stats HDIR >>> -----------------------------------

acm1 =[];
acc1 = [];
acv1 = [];
acs1 = [];
acr1 = [];
acm2 =[];
acc2 = [];
acv2 = [];
acs2 = [];
acr2 = [];
for hI = 1:bins.hdr.count
    hdr = bins.hdr.edges([hI,hI+1]);
for xI = 1:bins.x.count
    xbnd = [bins.x.edges(xI),bins.x.edges(xI+1)];
    for yI =1:bins.y.count
        ybnd = [bins.y.edges(yI),bins.y.edges(yI+1)];
        ind = within_ranges(sxy(:,1), xbnd) ...
            & within_ranges(sxy(:,2), ybnd) ...
            & within_ranges(shd, hdr);
        cm1 = circ_mean(sph(ind&scl==1));
        cv1 = circ_var(sph(ind&scl==1));
        cs1 = circ_std(sph(ind&scl==1));
        cr1 = circ_r(sph(ind&scl==1));
        cm2 = circ_mean(sph(ind&scl==2));
        cv2 = circ_var(sph(ind&scl==2));
        cs2 = circ_std(sph(ind&scl==2));
        cr2 = circ_r(sph(ind&scl==2));
        if isempty(sph(ind))
            acm1(xI,yI,hI) = nan;
            acc1(xI,yI,hI) = nan;
            acv1(xI,yI,hI) = nan;
            acs1(xI,yI,hI) = nan;
            acm2(xI,yI,hI) = nan;
            acc2(xI,yI,hI) = nan;
            acv2(xI,yI,hI) = nan;
            acs2(xI,yI,hI) = nan;
        else
            acm1(xI,yI,hI) = cm1;
            acv1(xI,yI,hI) = cv1;
            acs1(xI,yI,hI) = cs1;
            acr1(xI,yI,hI) = cr1;
            acc1(xI,yI,hI) = sum(ind&scl==1);
            acm2(xI,yI,hI) = cm2;
            acv2(xI,yI,hI) = cv2;
            acs2(xI,yI,hI) = cs2;
            acr2(xI,yI,hI) = cr2;
            acc2(xI,yI,hI) = sum(ind&scl==2);
        end
        
    end
end
end

% <<< interneuron theta circ stats HDIR <<< -----------------------------------

% >>> interneuron theta circ stats ALL >>> ------------------------------------
ccm1 =[]; ccv1 =[]; ccs =[]; ccr1 = [];
ccm2 =[]; ccv2 =[]; ccs2 =[]; ccr2 = [];
for xI = 1:bins.x.count
    xbnd = [bins.x.edges(xI),bins.x.edges(xI+1)];
    for yI =1:bins.y.count
        ybnd = [bins.y.edges(yI),bins.y.edges(yI+1)];
        ind = within_ranges(sxy(:,1), xbnd) ...
            & within_ranges(sxy(:,2), ybnd);
        cm1 = circ_mean(sph(ind&scl==1));
        cv1 = circ_var(sph(ind&scl==1));
        cs1 = circ_std(sph(ind&scl==1));
        cr1 = circ_r(sph(ind&scl==1));
        cm2 = circ_mean(sph(ind&scl==2));
        cv2 = circ_var(sph(ind&scl==2));
        cs2 = circ_std(sph(ind&scl==2));
        cr2 = circ_r(sph(ind&scl==2));
        if isempty(sph(ind))
            ccm1(xI,yI) = nan;
            ccv1(xI,yI) = nan;
            ccs1(xI,yI) = nan;
            ccr1(xI,yI) = nan;            
            ccm2(xI,yI) = nan;
            ccv2(xI,yI) = nan;
            ccs2(xI,yI) = nan;
            ccr2(xI,yI) = nan;            
        else
            ccm1(xI,yI) = cm1;
            ccv1(xI,yI) = cv1;
            ccs1(xI,yI) = cs1;
            ccr1(xI,yI) = cr1;
            ccm2(xI,yI) = cm2;
            ccv2(xI,yI) = cv2;
            ccs2(xI,yI) = cs2;
            ccr2(xI,yI) = cr2;
        end
    end
end
% <<< interneuron theta circ stats ALL <<< ------------------------------------

% >>> interneuron ccg G x,y >>> -----------------------------------------------
acg =[]; acg1 =[]; acg2 =[];
for hI = 1:bins.hdr.count
    hdr = bins.hdr.edges([hI,hI+1]);
    for xI = 1:bins.x.count
        xbnd = [bins.x.edges(xI),bins.x.edges(xI+1)];
        for yI =1:bins.y.count
            ybnd = [bins.y.edges(yI),bins.y.edges(yI+1)];
            ind = within_ranges(sxy(:,1), xbnd) ...
                  & within_ranges(sxy(:,2), ybnd) ...
                  & within_ranges(shd, hdr);
            [mccg, tbins] = CCG(ares(ind),scl(ind),4,25,sampleRate,[1,2],'count');
            if ndims(mccg)<3
                acg(xI,yI,hI,:) = nan([51,1]);
                acg1(xI,yI,hI,:) = nan([51,1]);
                acg2(xI,yI,hI,:) = nan([51,1]);
            else            
                acg(xI,yI,hI,:) = mccg(:,1,2);
                acg1(xI,yI,hI,:) = mccg(:,1,1);
                acg2(xI,yI,hI,:) = mccg(:,2,2);
            end
        end
    end
end

for hI = 1:bins.hdr.count
    for xI = 1:bins.x.count
        for yI = 1:bins.y.count
            if xcc(xI,yI,hI)<0.5
                acg(xI,yI,hI,:) = nan([51,1]);
            end
        end
    end
end
% <<< interneuron ccg G x,y <<< -----------------------------------------------

% >>> diagnostic figs >>> -----------------------------------------------------

% $$$ figure
% $$$ for hI = 1:bins.hdr.count
% $$$ subplot(2,4,hI);
% $$$ imagescnan({bins.x.centers,...
% $$$             bins.y.centers,...
% $$$             acg(:,:,hI,51)'},...
% $$$            [],'linear',true,'colorMap',@jet);
% $$$ axis('xy');
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ figure
% $$$ for hI = 1:8
% $$$ subplot(2,4,hI);
% $$$ imagescnan({bins.x.centers,...
% $$$             bins.y.centers,...
% $$$             log10(sum(acg(:,:,hI,:),4))'},...
% $$$            [0,4.5],'linear',true,'colorMap',@jet);
% $$$ axis('xy');
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ currpoint = get(gca,'CurrentPoint');
% $$$ currpoint = currpoint(1,1:2,1);
% $$$ xI = find( bins.x.edges>currpoint(1), 1, 'first');
% $$$ yI = find( bins.y.edges>currpoint(2), 1, 'first');
% $$$ hI = 7;
% $$$ figure,bar(tbins,sq(acg(xI,yI,hI,:)));

% <<< diagnostic figs <<< -----------------------------------------------------

% >>> Main Figure >>> ---------------------------------------------------------
xcc(xcc<0.5) = nan;
scg = sum(acg,4);
ncg = bsxfun(@rdivide,acg,xcc);
ncg1 = bsxfun(@rdivide,acg1,xcc);
ncg2 = bsxfun(@rdivide,acg2,xcc);
xI = 7;
yI = 7;
ym = max(max(sq(ncg(xI,yI,:,:)),[],2));
ny = bins.hdr.count;
nx = 12;
ym = 35;
figure,
for hI = 1:ny
    subplot2(ny,nx,hI,1);
        circle(0,0,1,'--k');
        theta = [];
        hb = bins.hdr.edges([hI,hI+1]);
        theta = [linspace([hb([1,end]), 100]), theta];
        patch([0 cos(theta) 0], ...
              [0 sin(theta) 0], ...
              'b',...
              'FaceAlpha', 0.5);
        daspect([1,1,1]);
    subplot2(ny,nx,hI,2);
        imagescnan({bins.x.centers,...
                    bins.y.centers,...
                    scg(:,:,hI)'},...
                   [0,ym],'linear',true,'colorMap',@jet);
        axis('xy');
    subplot2(ny,nx,hI,3);
        bar(tbins,sq(ncg(xI,yI,hI,:))./xcc(xI,yI,hI));
         ylim([0,ym]);
    subplot2(ny,nx,hI,4);
        bar(tbins,sq(ncg1(xI,yI,hI,:))./xcc(xI,yI,hI));
        ylim([0,ym]);
    subplot2(ny,nx,hI,5);
        bar(tbins,sq(ncg2(xI,yI,hI,:))./xcc(xI,yI,hI));
        ylim([0,ym]);
    subplot2(ny,nx,hI,6);
        imagescnan({bins.x.centers,...
                    bins.y.centers,...
                    circ_dist(acm1(:,:,hI),2.5)'},...
                   [-pi,pi],'circular',true,'colorMap',@hsv);
        axis('xy');
    subplot2(ny,nx,hI,7);
        imagescnan({bins.x.centers,...
                    bins.y.centers,...
                    circ_dist(acm2(:,:,hI),0)'},...
                   [-pi,pi],'circular',true,'colorMap',@hsv);
        axis('xy');
    subplot2(ny,nx,hI,8);
        imagescnan({bins.x.centers,...
                    bins.y.centers,...
                    mean(ncg(:,:,hI,[20:30]),4)'},...
                   [],'linear',true,'colorMap',@jet);
        axis('xy');
        axis('xy');
    subplot2(ny,nx,hI,9);
        imagescnan({bins.x.centers,...
                    bins.y.centers,...
                    mean(ncg(:,:,hI,[15:25]),4,'omitnan')'-mean(ncg(:,:,hI,[27:36]),4,'omitnan')'},...
                   [-1,1],'linear',true,'colorMap',@jet);
        axis('xy');
    subplot2(ny,nx,hI,10);
        imagescnan({bins.x.centers,...
                    bins.y.centers,...
                    xcc(:,:,hI)'},...
                   [0,5],'linear',true,'colorMap',@jet);
        axis('xy');
    subplot2(ny,nx,hI,11);
        imagescnan({bins.x.centers(1:end-1),...
                    bins.y.centers(1:end-1),...
                    acr1(:,:,hI)'},...
                   [0,0.8],'linear',true,'colorMap',@jet);
        axis('xy');
    subplot2(ny,nx,hI,12);
        imagescnan({bins.x.centers(1:end-1),...
                    bins.y.centers(1:end-1),...
                    acr2(:,:,hI)'},...                    
                   [0,.8],'linear',true,'colorMap',@jet);
        axis('xy');
end
ForAllSubplots('grid on;');
% <<< Main Figure <<< ---------------------------------------------------------


currpoint = get(gca,'CurrentPoint');
currpoint = currpoint(1,1:2,1);
xI = find( bins.x.edges>currpoint(1), 1, 'first');
yI = find( bins.y.edges>currpoint(2), 1, 'first');

ym = max(max(sq(ncg(xI,yI,:,:)),[],2));
for hI = 1:ny
    subplot2(ny,nx,hI,3);
        bar(tbins,sq(ncg(xI,yI,hI,:)));
        ylim([0,ym]);
    subplot2(ny,nx,hI,4);
        bar(tbins,sq(ncg1(xI,yI,hI,:)));
        ylim([0,ym]);
    subplot2(ny,nx,hI,5);
        bar(tbins,sq(ncg2(xI,yI,hI,:)));
        ylim([0,ym]);
end

sacg = permute(RectFilter(permute(ncg,[4,1,2,3]),5,5),[2,3,4,1]);

tacg = nan([bins.x.count,bins.y.count,bins.hdr.count]);
racg = nan([bins.x.count,bins.y.count,bins.hdr.count]);
for hI = 1:bins.hdr.count
    for xI = 1:bins.x.count
        for yI = 1:bins.y.count
            [mr,mi] = max(sq(sacg(xI,yI,hI,:)));
            racg(xI,yI,hI) = mr;
            tacg(xI,yI,hI) = tbins(mi);
        end
    end
end

racg(tacg==-400) = nan;
tacg(tacg==-400) = nan;


figure,
for hI = 1:bins.hdr.count
subplot(6,1,hI);
bar(tbins,sq(sacg(8,8,hI,:)))
end

figure
for hI = 1:bins.hdr.count
subplot(6,1,hI);
    imagescnan({bins.x.centers,...
                    bins.y.centers,...
                    circ_dist(acm1(:,:,hI),acm2(:,:,hI))'},...
                   [-pi,pi],'circular',true,'colorMap',@hsv);
        axis('xy');
end


figure,
for xI = 1:bins.x.count
    for yI = 1:bins.y.count
subplot2(15,15,yI,xI);
hold('on');
plot(bins.hdr.centers,sq(tacg(xI,yI,:))./2,'+-')
plot(bins.hdr.centers,sq(racg(xI,yI,:)).*10,'+-r')
plot(bins.hdr.centers,mean(sq(ncg1(xI,yI,:,:)),2).*3,'+-g')
plot(bins.hdr.centers,mean(sq(ncg2(xI,yI,:,:)),2).*3,'+-c')
ylim([-50,50]);
xlim([-pi,pi]);
    end
end



figure,
for hI = 1:bins.hdr.count
subplot(6,1,hI);
plot(reshape(tacg(:,:,hI),[],1)+randn([size(reshape(tacg(:,:,hI),[],1)),1]),reshape(racg(:,:,hI),[],1),'.')
xlim([-200,200]);
end



figure,
for hI = 1:bins.hdr.count
    subplot2(ny,2,hI,1);
        imagescnan({bins.x.centers,...
                    bins.y.centers,...
                    racg(:,:,hI)'},...                    
                   [1,3],'linear',true,'colorMap',@jet);
        axis('xy');
    subplot2(ny,2,hI,2);
        imagescnan({bins.x.centers,...
                    bins.y.centers,...
                    tacg(:,:,hI)'},...                    
                   [-75,75],'linear',true,'colorMap',@jet);
        axis('xy');
end
    

figure,
subplot(121);
histcirc(sph(scl==1))
subplot(122);
histcirc(sph(scl==2))


% Uhg
Trial = Trials{20};
unit = [41,43,48];
spk = Trial.load('spk',sampleRate,'',unit);
stc = copy(Trial.stc);
xyz = preproc_xyz(Trial,'trb',sampleRate);
fxyz = filter(copy(xyz),'ButFilter',4,2.5,'low');
vxy = vel(fxyz,'hcom',[1,2]);
phz = load_theta_phase(Trial,xyz);
hir = copy(xyz);
hir.data = atan2(xyz(:,'nose',2)-xyz(:,'hcom',2),xyz(:,'nose',1)-xyz(:,'hcom',1));
%sper = [stc{'hloc+lloc+lpause&theta',sampleRate}];
sper = [stc{'loc&theta',sampleRate}];

ufr = Trial.load('ufr', phz, spk, unit, 0.005);



cfet = LoadFet('/storage/gravio/data/processed/nlx/jg05-20120312/jg05-20120312.fet.6',inf,[8,10,15]);

res = spk(unit(1));
cfet = LoadFet(fullfile(Trial.spath,[Trial.name,'.fet.6']),[],[Trial.spk.map(Trial.spk.map(:,1)==unit(1),3)]);
pres = res(within_ranges(phz(res),[0,2*pi]));
figure,
scatter3(xyz(pres,'hcom',1), xyz(pres,'hcom',2), cfet(pres,7),10,phz(pres) ,'Filled');
caxis([0,2*pi]);
colormap('hsv');



figure();
hold('on');
plot(double((ufr(sper,2)-eps)&(ufr(sper,3)-eps)).*(ufr(sper,2)-ufr(sper,3)))
plot(phz(sper,1))

yres1 = spk(43);
yres2 = spk(48);




ycg = [];
for i = 1:10000
[ycg(:,:,:,i),ybins] = CCG([yres1;yres2],[ones(size(yres1)); 2*ones(size(yres2))], 2, 25, sampleRate, [1,2],'count',[1,128]+64*i);
end

figure,bar(ybins,ycg(:,1,2));
figure,imagesc(sq(ycg(:,1,2,:)));

ylim([-0.5,1.5])

hxy = fet_H
axy = copy(vxy);
axy.data = [0;diff(vxy.data)];
hvfl = fet_hvfl(Trial,sampleRate);
figure,
plot(ufr(sper), hvfl(sper,1),'.');
pch = fet_HB_pitch(Trial,sampleRate);

figure,
hist2([ufr(sper), hvfl(sper,2)],50,50,'yprob');

figure,
hist2([ufr(sper), hvfl(sper,1)],50,50,'yprob');
caxis([0,0.15]);
colormap('jet');

figure,
hist2([ufr(sper,1), ufr(sper,2)],50,50,'yprob');
caxis([0,0.15]);
colormap('jet');

figure,
hist2([ufr(sper,1), ufr(sper,3)],50,50,'yprob');
caxis([0,0.15]);
colormap('jet');

figure,
hist2([ufr(sper,2), ufr(sper,3)],50,50,'yprob');
caxis([0,0.15]);
colormap('jet');



figure,
hist2([ufr(sper), hvfl(sper,2)],50,50,'yprob');
caxis([0,0.15]);
colormap('jet');

figure,
histcirc(phz(spk(unit(1))));

figure,
plot(phz(spk(unit)), hvfl(spk(unit),1), '.')

figure
plot(phz(spk(unit)), hir(spk(unit)), '.')

res = spk(unit(1);
res2 = spk(unit(2);
figure, scatter3(xyz(res,'hcom',1), xyz(res,'hcom',2), hir(res),10, ufr(res),'Filled'); colormap('jet');
zlim([0,pi/3]);


pres = res(within_ranges(phz(res),[pi,2*pi]));
figure,
scatter3(xyz(pres,'hcom',1), xyz(pres,'hcom',2), phz(pres),10,pch(pres,3) ,'Filled');
caxis([-1.8,.5]);
colormap('hsv');

figure,
hist2([pch(res,3), ufr(res)],50,50,'yprob');
caxis([0,0.15]);
colormap('jet');

figure,plot((ufr(sper,2)&ufr(sper,3)).*(ufr(sper,2)+ufr(sper,3)))
figure,plot(ufr(sper,2)&ufr(sper,3));



figure
scatter(ufr(sper,2), ufr(sper,3), 10, pch(sper,2), 'Filled');
colormap('jet');

