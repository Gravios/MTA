% MjgER2016_decode_xy
%
% Section 1. Loading data and decoding position
% Session 2. Plot 


MjgER2016_load_data();

% $$$ if ~exist('pfd','var'), [pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(Trials);  end
% $$$ numComp = size(eigVec{1},2);
% $$$ pfindex = 1;
% $$$ MjgER2016_load_bhv_erpPCA_scores();
% $$$ % output:
% $$$ %    fsrcz
% $$$ %    FSrC
% $$$ %    rmaps
% $$$ clear('pfd','tags','eigVec','eigVar','eigScore','validDims','zrmMean','zrmStd',...
% $$$       'clu','tlu','rind','D','LR','FSCFr','rsMean','rsStd','pfdShuffled','rmapsShuffledMean',...
% $$$       'rmapsShuffled','FSrM','FSrS','fsrsMean','fsrsStd','fsrsMean','fsrsStd');
 

% SET analysis parameters
sampleRate = 250;   % Hz
spikeWindow = 0.03; % ms [0.75,0.05,0.025,0.02,0.015]
mode = 'xy'; % alt vals: 'xy'
posteriorMaxThreshold = 0.001;

% $$$ 
% $$$ statesPfs = {'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta',...
% $$$           'lpause&theta'};
% $$$ 
% $$$ states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit','turn','walk'};%,'ripple'};
% $$$ 
% $$$ 
% $$$ thetaPerPeak = resample(cast([Trial.stc{'theta-groom-sit'}],'TimeSeries'),xyz);
% $$$ thetaPerTrgh = copy(thetaPerPeak);
% $$$ 
% $$$ thetaPerPeak.label = 'thetaPeak';
% $$$ thetaPerPeak.data(phz.data<0) = 0;
% $$$ thetaPerPeak = cast(thetaPerPeak,'TimePeriods');
% $$$ 
% $$$ thetaPerTrgh.label = 'thetaTrgh';
% $$$ thetaPerTrgh.data(phz.data>0) = 0;
% $$$ thetaPerTrgh = cast(thetaPerTrgh,'TimePeriods');



%-LOAD-DATA------------------------------------------------------------------------------%


% DECODE position from placefields and unit firing rate for all sessions
tind = [3:5,17:23];

Trials = Trials(tind);
units = units(tind);

%cf(@(t,u)  bhv_decode(t,sampleRate,[],u,mode,[],[],spikeWindow,true),  Trials(tind), units(tind));


% LOAD subject position object
xyz = cf(@(t)  resample(preproc_xyz(t,'trb'),sampleRate), Trials);

% TRIAL inds
dtind   = cf(@(x,t)  t.*ones([size(x,1),1]), xyz,num2cell(tind));    
dtind   = cat( 1, dtind{:}  );

% COMPUTE polar coordinates of horizontal position
dmcd = cf(@(x) sqrt(sum(x(:,'hcom',[1,2]).^2,3)),xyz);
dmcd = cat(1, dmcd{:});
dmca = cf(@(x) circ_dist(atan2(x(:,'hcom',2),x(:,'hcom',1)),...
                 atan2(diff(x(:,{'hcom','nose'},2),1,2),diff(x(:,{'hcom','nose'},1),1,2))),...
          xyz);
dmca = cat(1, dmca{:});

% LOAD lfp
cf(@(t) set(t.lfp,'filename',[t.name,'.lfp']), Trials);

lfp = cf(@(t,c) load(t,'lfp',c), Trials, num2cell([sessionList(tind).thetaRefGeneral]));
% COMPUTE theta LFP phase
dphz = cf(@(l) l.phase([6,12]), lfp);
cf(@(p) set(p,'data',unwrap(p.data)), dphz);
cf(@(p,x) p.resample(x), dphz,xyz);    
dphz = cf(@(p) mod(p.data+pi,2*pi)-pi, dphz);;
%clear('lfp');
dphz = cat(1,dphz{:});
dphz(ismember(dtind,[3,4,5])) = dphz(ismember(dtind,[3,4,5]))+2*0.78;
dphz(~ismember(dtind,[3,4,5])) = dphz(~ismember(dtind,[3,4,5]))+0.78;
dphz(dphz<0) = dphz(dphz<0)+2*pi;



% COMPUTE speed
% $$$ vxy = vel(filter(copy(xyz),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
% $$$ vxy.data(vxy.data<1e-3) = 1e-3;
% $$$ vxy.data = log10(vxy.data);
%fvxy = vel(filter(copy(xyz),'ButFilter',3,1.5,'low'),{'spine_lower','hcom'},[1,2]);




tRot = {0,0,0,0.17,0.17,0.17,0.17,0.17,0.17,0.17};
%pch  = cf(@(t)    fet_HB_pitchB(t,sampleRate),                                   Trials   );
hvec = cf(@(x)    x(:,'head_front',[1,2])-x(:,'head_back',[1,2]),                xyz      );
hvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                      hvec     );
hvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                     hvec     );
hvec = cf(@(h,r)  multiprod(h,[cos(r),-sin(r);sin(r),cos(r)],[2,3],[1,2]),       hvec, tRot);
tvec = cf(@(x)    circshift(x(:,'hcom',[1,2]),-1)-circshift(x(:,'hcom',[1,2]),1),xyz      );
tvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                      tvec     );
tvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                     tvec     );
stc  = cf(@(t)    t.load('stc','msnn_ppsvd_raux'),                               Trials       );

dstcm = cf(@(s,x)  stc2mat(s,x,states),                                           stc,xyz);
dstcm = cat(1, dstcm{:});

% HBP 
hbp = cf(@(p) p(:,1), pch);
hbp = cat(1,hbp{:});


% CONVERT MTAStateCollection into a state matrix 
stcm = stc2mat(stc,xyz,states);
% COMPUTE unit firing rate

% UNIT inclusion
duinc = cf(@(t,u,x) load(t,'ufr',x,[],u,spikeWindow,true,'gauss'), Trials,units,xyz);    
duinc = cf(@(u) sum(u.data>0.2,2), duinc);
duinc = cat(1,duinc{:});

% DECODE position from placefields and unit firing rate
overwrite = false;
[posEstCom,posEstMax,posEstSax,posteriorMax] = ...
    cf(@(t,u) bhv_decode(t,sampleRate,[],u,mode,[],[],spikeWindow,overwrite), Trials,units);

%-END-LOAD-DATA------------------------------------------------------------------------------%

% DIAGNOSTICS of req20191104


% HBA 
dhba = [];
for t = 1:numel(Trials),
    txyz = filter(copy(xyz{t}),'ButFilter',3,30,'low');    
    xycoor = cat(2,...
                 txyz(:,'hcom',[1,2])-txyz(:,'bcom',[1,2]),...
                 txyz(:,'nose',[1,2])-txyz(:,'hcom',[1,2]));
    tht = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    dhba = cat(1,dhba,circ_dist(tht(:,2),tht(:,1)));
end


% HRL 
dhrl = cf(@(t,x) transform_origin(t,filter(copy(x),'ButFilter',4,20),'hcom','nose',{'hcom','head_right'}), ...
          Trials, xyz);
dhrl = cf(@(h) real(h.roll), dhrl);
dhrl  = cat( 1, dhrl{:} );

% HVL 
% HVF 
hrv = cf(@(t) fet_href_HXY(t,sampleRate,false,'trb',[],4,0), Trials);
dhvf = cf(@(h) h(:,1), hrv);
dhvf = cat(1,dhvf{:});
dhvl = cf(@(h) h(:,2), hrv);
dhvl = cat(1,dhvl{:});

% DERROR 
dError = cf(@(e,x,h,t) ...
            [multiprod(e(:,[1,2])-sq(x(:,'hcom',[1,2])),h,2,[2,3]), ...
             multiprod(e(:,[1,2])-sq(x(:,'hcom',[1,2])),t,2,[2,3])], ...
            posEstSax, xyz, hvec, tvec);
dError = cat(1,dError{:});
dError(:,2) = dError(:,2)+16;

posteriorMax = cat(1,posteriorMax{:});

chrl = (dhrl-0.26*double(~ismember(dtind,[3,4,5])));
chba = -(dhba+0.2-0.4*double(~ismember(dtind,[3,4,5])));

%chvang = dhvang;
chvl = dhvl;
cError = dError;
a = chba;
chba  (a<0) = -chba    (a<0);
cError(a<0,2) = -cError(a<0,2);
cError(a<0,4) = -cError(a<0,4);
chvl  (a<0) = -chvl    (a<0);
chrl  (a<0) = -chrl    (a<0);
%chvang  ( a<0 ) = -chvang   ( a<0 );


hbaBinEdges = linspace(0,1.2,8);;
hbaBinCenters = (hbaBinEdges(2:end)+hbaBinEdges(1:end-1))./2;
hbaBinInd = discretize(chba,hbaBinEdges);

hvlBinEdges = [-60,-10,-2,2,10,60];
hvlBinCenters = (hvlBinEdges(2:end)+hvlBinEdges(1:end-1))./2;
hvlBinInd = discretize(chvl,hvlBinEdges);

hvfBinEdges = [-20,-10,-2,2,10,20,40,80];
hvfBinCenters = (hvfBinEdges(2:end)+hvfBinEdges(1:end-1))./2;
hvfBinInd = discretize(dhvf,hvfBinEdges);

hrlBinEdges = [-0.6,-0.3,-0.1,0.1,0.3,0.6];
hrlBinCenters = (hrlBinEdges(2:end)+hrlBinEdges(1:end-1))./2;
hrlBinInd = discretize(chrl,hrlBinEdges);


phzBinEdges = linspace(0,2*pi,26);
phzBinCenters = (phzBinEdges(2:end)+phzBinEdges(1:end-1))./2;
phzBinInd = discretize(dphz,phzBinEdges);

decBinEdges = linspace(-500,500,100);

ind = duinc>=1&duinc<=6 ...
      & posteriorMax>0.002 ...
      & dstcm(:,1)         ...   
      & any(dstcm(:,[9]),2) ...
      & dmcd < thresholds.mazeCenterDist;
%      & (dmcd < thresholds.mazeCenterDist  |  abs(dmca) < thresholds.mazeCenterAng);

txhv = hbaBinInd(ind);
tyhv = hvlBinInd(ind);
%tyhv = hrlBinInd(ind);
tphz = dphz(ind);
tpbi = phzBinInd(ind);
tdec = cError(ind,:);

figure();
for x = 1:numel(hbaBinCenters),
    for y = 1:numel(hvlBinCenters),
        dind = txhv==x & tyhv==y;    
        subplot2(numel(hvlBinCenters),numel(hbaBinCenters),y,x);
        hist2([tdec(dind,2),tphz(dind)],decBinEdges,phzBinEdges);
        for p = 1:numel(phzBinCenters),
            mj(p,x,y) = median(tdec(dind&tpbi==p,2),'omitnan');
        end
        
    end
end
ForAllSubplots('xlim([-300,300]);');
ForAllSubplots('Lines(0,[],''k'');');

figure,
for p = 1:25,
    subplot(1,25,p);
    imagesc(sq(mj(p,:,:))')
    axis('xy');
    caxis([-100,100]);
end


txhv = hbaBinInd(ind);
tyhv = hvfBinInd(ind);
tphz = dphz(ind);
tdec = cError(ind,:);
figure();
for x = 1:numel(hbaBinCenters),
    for y = 1:numel(hvfBinCenters),
        dind = txhv==x&tyhv==y;    
        subplot2(numel(hvfBinCenters),numel(hbaBinCenters),y,x);
        hist2([tdec(dind,1),tphz(dind)], decBinEdges,phzBinEdges);
    end
end
ForAllSubplots('xlim([-300,300]);');


% START SECTION ------------------------------------------------------------------------------------
% Section 2: phase precession (drz X theta phase) color coded as circular distance between head 
%            direction and vector to placefield center.
% CURRENT : there exists high rate place fields which demonstrate both distance and angle coding in the
%           the phase of theta.
% TO DO   : create a list of high rate units which show this phenomena


% OMIT points where rat faceing away from the maze center near the boundary. 
thresholds.mazeCenterDist = 380;
thresholds.mazeCenterAng = pi/2;
cind =    mazeCenterDist < thresholds.mazeCenterDist  |  abs(mazeCenterAng) < thresholds.mazeCenterAng ;

% PLOT Place field along side the theta phase as a function of drz. Points colored as heading
%      to place field
figure();
for unit = unitSubset,
    subplot(331),
    plot(pft,unit,'mean',true,[],true);
    title(['unit: ',num2str(unit)]);
    res = spk(unit);
    sres = res(WithinRanges(res,[stc{'x&t',sampleRate}.data]));
    sres = sres(abs(ddz(sres,unit==unitSubset))<250 & cind(sres));
    
    subplot(332),
    scatter([drz(sres,unit==unitSubset);drz(sres,unit==unitSubset)],...
         [phz(sres);phz(sres)+2*pi],10,abs([drang(sres,unit==unitSubset);drang(sres,unit==unitSubset)]), ...
            'filled');
    colorbar();    
    sres = res(WithinRanges(res,[stc{'x&t',sampleRate}.data]));
    sres = sres(abs(ddz(sres,unit==unitSubset))<250);    
    subplot(333),
    scatter([hrz(sres,unit==unitSubset);hrz(sres,unit==unitSubset)],...
         [phz(sres);phz(sres)+2*pi],10,abs([hrang(sres,unit==unitSubset);hrang(sres,unit==unitSubset)]),'filled');
    colorbar();    
    sres = res(WithinRanges(res,[stc{'p&t',sampleRate}.data]));
    sres = sres(abs(ddz(sres,unit==unitSubset))<250);    
    subplot(335),
    scatter([drz(sres,unit==unitSubset);drz(sres,unit==unitSubset)],...
         [phz(sres);phz(sres)+2*pi],10,abs([drang(sres,unit==unitSubset);drang(sres,unit==unitSubset)]),'filled');
    colorbar();
    sres = res(WithinRanges(res,[stc{'p&t',sampleRate}.data]));
    sres = sres(abs(ddz(sres,unit==unitSubset))<250);    
    subplot(336),
    scatter([hrz(sres,unit==unitSubset);hrz(sres,unit==unitSubset)],...
         [phz(sres);phz(sres)+2*pi],10,abs([hrang(sres,unit==unitSubset);hrang(sres,unit==unitSubset)]), ...
            'filled');
    colorbar();    
    sres = res(WithinRanges(res,[stc{'r&t',sampleRate}.data]));
    sres = sres(abs(ddz(sres,unit==unitSubset))<250);    
    subplot(338),
    scatter([drz(sres,unit==unitSubset);drz(sres,unit==unitSubset)],...
         [phz(sres);phz(sres)+2*pi],10,[fet(sres,2);fet(sres,2)],'filled');
    colorbar();
    sres = res(WithinRanges(res,[stc{'r&t',sampleRate}.data]));
    sres = sres(abs(ddz(sres,unit==unitSubset))<250);    
    subplot(339),
    scatter([hrz(sres,unit==unitSubset);hrz(sres,unit==unitSubset)],...
         [phz(sres);phz(sres)+2*pi],10,[fet(sres,2);fet(sres,2)], ...
            'filled');
    colorbar();
    colormap('jet');
    waitforbuttonpress();
end

% section 2
% END SECTION --------------------------------------------------------------------------------------



% DIAGNOSTICS 
% $$$ overwrite = false;
% $$$ [decPosC_long,decPosM_long,desPosS_long,posteriorMax_long] = ...
% $$$     bhv_decode(Trial,sampleRate,unitSubset,mode,[],[],0.75,overwrite);
% $$$ [decPosC_short,decPosM_short,desPosS_short,posteriorMax_short] = ...
% $$$     bhv_decode(Trial,sampleRate,unitSubset,mode,[],[],0.02,overwrite);
% $$$ % REVIEW decoded output
% $$$ ff = sq(xyz(:,'head_back',[1,2]));
% $$$ nf = size(ff,2);
% $$$ figure,
% $$$ for s = 1:nf
% $$$     subplot(nf+1,1,s);
% $$$     hold('on');
% $$$     plot(ff(:,s))
% $$$     plot(posEstCom(:,s));
% $$$     plot(posEstMax(:,s));
% $$$     plot(posEstSax(:,s));    
% $$$ end
% $$$ subplot(nf+1,1,nf+1);
% $$$ plot(unitInclusion)
% $$$ linkaxes(findobj(gcf,'Type','Axes'),'x');
% $$$ % EXPLORE the decoded positions' angle to head position 
% $$$ ang = create(MTADang,Trial,filter(copy(xyz),'RectFilter',3,3));
% $$$ ind = ':';
% $$$ decError = [multiprod(posEstSax(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3])];
% $$$ decError(~(unitInclusion>=2  & posteriorMax>0.002)) = nan;
% $$$ hxy = fet_href_H(Trial,sampleRate);
% $$$ figure,
% $$$ subplot(3,1,[1,2]);
% $$$ hold('on');
% $$$ %plot(circ_dist(circshift(ang(:,'head_back','head_front',1),-1),ang(:,'head_back','head_front',1)).*ang.sampleRate);
% $$$ plot(hxy(:,2)./10);
% $$$ %plot(diff(RectFilter(circ_dist(circshift(ang(:,'head_back','head_front',1),-1),ang(:,'head_back','head_front',1)),3,3).*ang.sampleRate));
% $$$ plot(atan2(decError(:,2),decError(:,1)),'.r');
% $$$ plot(phz(:,1),'c')
% $$$ Lines([],0,'k');
% $$$ Lines([],-pi,'k');
% $$$ Lines([],pi,'k');
% $$$ Lines([],pi/2,'b');
% $$$ Lines([],-pi/2,'b');
% $$$ ylim([-10,10])
% $$$ subplot(3,1,[3]);
% $$$ plotSTC(stc,ang.sampleRate,[],fliplr({'rear','loc','pause'}),fliplr('rbk'));
% $$$ ylim([0,4]);
% $$$ linkaxes(findobj(gcf,'Type','Axes'),'x')



% START section ------------------------------------------------------------------------------------
% Section 3: Theta phase as a function of decoded spatial representation projected onto 2d 
%            orthonormal basis aligned to the head, partitioned by behavioral state.

thresholds.unitInclusion = 2;
thresholds.posteriorMax = 0.002;
thresholds.mazeCenterDist = 380;
thresholds.mazeCenterAng = pi/2;

cind =   unitInclusion  >=  thresholds.unitInclusion         ...
       & posteriorMax   >   thresholds.posteriorMax          ...
       & (  mazeCenterDist      <  thresholds.mazeCenterDist ...
         |  abs(mazeCenterAng)  <  thresholds.mazeCenterAng  );


normalize = @(out) bsxfun(@rdivide,out,max(out,[],2)); tagNormalization = 'ColumnMax';
normalize = @(out) out; tagNormalization = 'None';
hfig = figure();
pause(0.2);
hfig.Position =[0, 0, 775, 1400];
pause(0.2);
ns = 6;
for s = 1:ns; %states
sind =      stcm(:,1)   ...
         &  any(stcm(:,s),2);
ind = cind&sind;

decError = multiprod(posEstSax(ind,[1,2])-sq(mean(xyz(ind,'hcom',[1,2]),2)),hvec(ind,:,:),2,[2,3]);
%decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3]);
% Forward decoding projection onto head
subplot2(ns,3,s,1);
out = hist2([[decError(:,1);decError(:,1)],...
             [phz(ind,1);phz(ind,1)+2*pi]],...
            linspace(-300,300,50), ...
            linspace(-pi,3*pi,30));
imagesc(linspace(-300,300,50),linspace(-pi,pi*3,30),normalize(out)');axis('xy');
title('forward projection')
xlabel('mm');
ylabel({states{s},'theta phase'});
% Lateral decoding projection onto head
subplot2(ns,3,s,2);
out = hist2([[decError(:,2);decError(:,2)],...
             [phz(ind,1);phz(ind,1)+2*pi]],...
            linspace(-300,300,50), ...
            linspace(-pi,3*pi,30));
imagesc(linspace(-300,300,50),linspace(-pi,pi*3,30),normalize(out)');axis('xy');
title('lateral projection')
xlabel('mm');
ylabel('theta phase');
% Lateral decoding projection onto head
decAng = atan2(decError(:,2),decError(:,1));
subplot2(ns,3,s,3);
out = hist2([[decAng(:,1);decAng(:,1)],...
             [phz(ind,1);phz(ind,1)+2*pi]],...
            linspace(-pi,pi,50), ...
            linspace(-pi,3*pi,30));
imagesc(linspace(-pi,pi,50),linspace(-pi,3*pi,30),normalize(out)');axis('xy');
title('angle of projection')
xlabel('rad');
ylabel('theta phase');
end

hax = findobj(hfig,'Type','Axes');
af(@(ax) set(ax,'Units','centimeters'), hax);
af(@(ax) set(ax,'Position',[ax.Position(1:2),2.5,2.5]), hax);

fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
xlim([0,hfig.Position(3)]);
ylim([0,hfig.Position(4)]);
text( 100,hfig.Position(4)-60,{['Theta phase as a function of projected decoded position on to head frame of reference.'],...
          ['Normalization: ',tagNormalization],...
          ['Session: ',Trial.filebase],...
          ['Constrants:'],...
          ['unitInclusion >= ',num2str(thresholds.unitInclusion) ],...
          ['posteriorMax  >  ',num2str(thresholds.posteriorMax)  ],...
          ['(mazeCenterDist < ',num2str(thresholds.mazeCenterDist),...
           ' | abs(mazeCenterAng) < ',num2str(thresholds.mazeCenterAng)]});

print(hfig,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_decodedXthetaPhaseXstates_',Trial.filebase,'_Norm-',tagNormalization,'.eps']]);
print(hfig,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
      ['dc_decodedXthetaPhaseXstates_',Trial.filebase,'_Norm-',tagNormalization,'.png']]);


% END S3 -------------------------------------------------------------------------------------------



% S4 -----------------------------------------------------------------------------------------------
% Section 4: Theta phase as a function of decoded spatial representation projected onto 2d 
%            orthonormal basis aligned to the head, partitioned by speed, and omitting rearing 
%            periods.

% ANALYSIS -----------------------------------------------------------------------------------------
% CIRCULAR LINEAR correlation between theta phase (LFP{pyr,5-12Hz}) and projection
%   of the decoded position on longitudinal axis of the head

% DEFINITIONS ----
% DPTR : decoded head position in head trajectory reference frame (DPTR) 
% DPHR : decoded head position in head reference frame (DPHR) 
% TP   : Theta phase


cind =   unitInclusion  >=  thresholds.unitInclusion         ...
       & posteriorMax   >   thresholds.posteriorMax          ...
       & (  mazeCenterDist      <  thresholds.mazeCenterDist ...
         |  abs(mazeCenterAng)  <  thresholds.mazeCenterAng  );

sind =      stcm(:,1)   ...   
         & ~stcm(:,2)   ...
         & ~stcm(:,8)   ...          
         & ~stcm(:,7);         

sum(cind&sind)

decPos = posEstCom(:,[1,2]);
%decPos = posEstSax(:,[1,2]);

tags = {'body','head'};
speedInd = 2,
for speedInd = 1:2,

    
    tag = tags{speedInd};

% COMPUTE the circular-linear correlation of DPTR vs TP
    % SETUP data struct to collect model parameters across speed partions
    clear('modelTraj'); 
    modelTraj.partitions = prctile(vxy(sind,speedInd),linspace(0,100,7));
    % VISUALDOC partitioning of head speed space
% $$$     figure,
% $$$     hist(vxy(sind,speedInd),1000)
% $$$     Lines(modelTraj.partitions,[],'r');
    % SET fiting resolution and range
    fitRes = 10000;
    fitRange = linspace(-2*pi,2*pi,fitRes);
    % FIT circular linear model for each speed partition 
    for v = 1:numel(modelTraj.partitions)-1,
        ind = cind  &  sind      ...
              &  vxy(:,speedInd)>modelTraj.partitions(v)       ...
              &  vxy(:,speedInd)<modelTraj.partitions(v+1);
        % COMPUTE projection of the decoded position on longitudinal axis of the head
        decError = multiprod(decPos(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
        decPhase = phz(ind,1);
        % ASSIGN inputs
        lin = decError(:,1);
        circ = decPhase;
        % TRANSFORM linear and circular components into cos and sin parts
        cosPart = sum(cos(bsxfun(@minus,circ,2*pi*lin*fitRange)),1);
        sinPart = sum(sin(bsxfun(@minus,circ,2*pi*lin*fitRange)),1);
        R = sqrt((cosPart./length(circ)).^2 + ...
                 (sinPart./length(circ)).^2 );
        [lmi,lmv] = LocalMinima(-R',round(fitRes/10),0,1);
        % SAVE model parameters
        modelTraj.R(v) = R(lmi); % R: fit quality 
        modelTraj.parameters(v,:) = [2*pi*fitRange(lmi),... %Regression Parm : Slope
                            atan2(sinPart(lmi),cosPart(lmi))];%: Offset
                                                              % COMPUTE circular-linear correlation coefficient
        linC = mod(abs(modelTraj.parameters(v,1))*lin,2*pi);
        circMean = atan2(sum(sin(circ)),sum(cos(circ)));
        linCMean = atan2(sum(sin(linC)),sum(cos(linC)));
        modelTraj.rho(v) = sum(sin(circ-circMean).*sin(linC-linCMean))...
            ./sqrt(sum(sin(circ-circMean).^2).*sum(sin(linC-linCMean).^2));
    end

% COMPUTE the circular-linear correlation of DPHR vs TP
    clear('modelHead'); 
    modelHead.partitions = modelTraj.partitions;
    % SET fiting resolution and range
    fitRes = 10000;
    fitRange = linspace(-pi,pi,fitRes);
    % FIT circular linear model for each speed partition 
    for v = 1:numel(modelHead.partitions)-1,
        ind = cind  &  sind      ...
              &  vxy(:,speedInd)>modelHead.partitions(v)       ...
              &  vxy(:,speedInd)<modelHead.partitions(v+1);
        % COMPUTE projection of the decoded position on longitudinal axis of the head
        decError = multiprod(decPos(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3]);
        decPhase = phz(ind,1);
        % ASSIGN inputs
        lin = decError(:,1);
        circ = decPhase;
        % TRANSFORM linear and circular components into cos and sin parts
        cosPart = sum(cos(bsxfun(@minus,circ,2*pi*bsxfun(@times,fitRange,lin))),1);
        sinPart = sum(sin(bsxfun(@minus,circ,2*pi*bsxfun(@times,fitRange,lin))),1);
        R = sqrt((cosPart./length(circ)).^2+...
                 (sinPart./length(circ)).^2 );
        [lmi,lmv] = LocalMinima(-R',round(fitRes/10),0,1);
        % SAVE model parameters
        modelHead.R(v) = R(lmi); % R: fit quality 
        modelHead.parameters(v,:) = [2*pi*fitRange(lmi)',... %Regression Parm : Slope
                            atan2(sinPart(lmi),cosPart(lmi))'];%: Offset

        linC = mod(abs(modelHead.parameters(v,1))*lin,2*pi);
        circMean = atan2(sum(sin(circ)),sum(cos(circ)));
        linCMean = atan2(sum(sin(linC)),sum(cos(linC)));
        modelHead.rho(v) = sum(sin(circ-circMean).*sin(linC-linCMean))...
            ./sqrt(sum(sin(circ-circMean).^2).*sum(sin(linC-linCMean).^2));
    end



% VISUALIZE circular-linear relation ship between over range of speeds
    hfig = figure();
    % STATE - locomotion
    nPart = numel(modelTraj.partitions)-1;

    clf();
    hfig.Units = 'centimeters';
    hfig.Position = [0,0,30,5*nPart];
    for v = 1:nPart
        y = nPart-v+1;
        % select data with in speed partition v
        ind = cind  &  sind  ...
              &  vxy(:,speedInd)>modelTraj.partitions(v)   ...
              &  vxy(:,speedInd)<modelTraj.partitions(v+1);
        % REFERENCE trajectory coordinate system
        %decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
        decError = multiprod(decPos(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
        
        % Anterior-Posterior axis
        subplot2(nPart,6,y,1);
        % JPDF phase X error
        hist2([[decError(:,1);decError(:,1)],...
               [phz(ind,1);phz(ind,1)+pi*2]],...
              linspace(-300,300,40),...
              linspace(-pi,pi*3,30)); 
        line([-300,300],polyval(modelTraj.parameters(v,:),[-300,300])-2*pi,'Color','m');    
        line([-300,300],polyval(modelTraj.parameters(v,:),[-300,300]),     'Color','m');
        line([-300,300],polyval(modelTraj.parameters(v,:),[-300,300])+2*pi,'Color','m');
        title(['rho: ',num2str(modelTraj.rho(v))]);
        xlabel('mm');
        ylabel('theta phase');
        Lines(0,[],'k');
        % medial-lateral axis
        subplot2(nPart,6,y,2);
        hist2([[decError(:,2);decError(:,2)],...
               [phz(ind,1);phz(ind,1)+pi*2]],...
              linspace(-300,300,40),...
              linspace(-pi,pi*3,30)); 
        xlabel('mm');
        ylabel('theta phase');
        Lines(0,[],'k');
        % direction
        subplot2(nPart,6,y,3);
        hist2([[atan2([decError(:,2);decError(:,2)],...
                      [decError(:,1);decError(:,1)])],...
               [phz(ind,1);phz(ind,1)+pi*2]],...
              linspace(-pi,pi,40),...
              linspace(-pi,pi*3,30)); 
        xlabel('yaw');
        ylabel('theta phase');
        Lines(0,[],'k');
        % REFERENCE head coordinate system
        %decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3]);    
        decError = multiprod(decPos(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3]);    
        subplot2(nPart,6,y,4);
        hist2([[decError(:,1);decError(:,1)],...
               [phz(ind,1);phz(ind,1)+pi*2]],...
              linspace(-300,300,40),...
              linspace(-pi,pi*3,30)); 
        line([-300,300],polyval(modelHead.parameters(v,:),[-300,300])-2*pi,'Color','m');
        line([-300,300],polyval(modelHead.parameters(v,:),[-300,300]),'Color','m');
        line([-300,300],polyval(modelHead.parameters(v,:),[-300,300])+2*pi,'Color','m');
        
        title(['rho: ',num2str(modelHead.rho(v))]);
        xlabel('mm');
        ylabel('theta phase');
        Lines(0,[],'k');
        % lateral
        subplot2(nPart,6,y,5);
        hist2([[decError(:,2);decError(:,2)],...
               [phz(ind,1);phz(ind,1)+pi*2]],...
              linspace(-300,300,40),...
              linspace(-pi,pi*3,30)); 
        xlabel('mm');
        ylabel('theta phase');
        Lines(0,[],'k');
        % yaw
        subplot2(nPart,6,y,6);
        hist2([atan2([decError(:,2);decError(:,2)],...
                     [decError(:,1);decError(:,1)]),...
               [phz(ind,1);phz(ind,1)+pi*2]],...
              linspace(-pi,pi,40),...
              linspace(-pi,pi*3,30));
        Lines(0,[],'k');
        xlabel('yaw');
        ylabel('theta phase');
    end


    ForAllSubplots('hax=gca;hax.Units=''centimeters'';hax.Position(3:4)=[3,2];')

    af(@(hax) set(hax,'Position',hax.Position+[0,-2,0,0]), findobj(gcf,'Type','Axes'));

    axPos = cell2mat([get(findobj(gcf,'Type','Axes'),'Position')]);
    axYPos = unique(axPos(:,2));

    fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
    xlim([0,hfig.Position(3)]);
    ylim([0,hfig.Position(4)]);
    Lines(hfig.Position(3)/2+0.2,[],'k');

    Lines([],axYPos+3,'k');

    for v = 1:nPart, 
        text(1,axYPos(nPart-v+1)+1,...
             {[tag,' speed:'],...
              num2str(modelTraj.partitions(nPart-v+2)),...
              '',...
              num2str(modelTraj.partitions(nPart-v+1))},...
             'Rotation',0);
    end

    ht = text(hfig.Position(3).*0.3,hfig.Position(4)*.9,...
              'Head Movement Basis', 'HorizontalAlignment','center');
    ht = text(hfig.Position(3).*0.7,hfig.Position(4)*.9,...
              'Head Direction Basis','HorizontalAlignment','center');

    
    text(1,hfig.Position(4)*.95,...
         {[Trial.filebase,': locomation and pause'],...
          ['unitInclusion >= ',num2str(thresholds.unitInclusion) ],...
          ['posteriorMax  >  ',num2str(thresholds.posteriorMax)  ],...
          ['JPDF of projected decoded positions and theta phase,'],...
          ['partition over equal partitions of ' tag ' speed'    ],...
          ['(mazeCenterDist < ',num2str(thresholds.mazeCenterDist),...
           ' | abs(mazeCenterAng) < ',num2str(thresholds.mazeCenterAng)]});

    print(hfig,'-depsc2',...
          ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
           ['dc_decodedXthetaPhaseX',tag,'Speed_',mode,'_',Trial.filebase,'.eps']]);
    print(hfig,'-dpng',...
          ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
           ['dc_decodedXthetaPhaseX',tag,'Speed_',mode,'_',Trial.filebase,'.png']]);

% SAVE vars
    save(fullfile(Trial.spath,...
                  [Trial.filebase,'_dc_decodedXthetaPhaseX',tag,'Speed_',mode,'.mat']),...
         'modelTraj',...
         'modelHead');
end


% Section 4
% END SECTION ------------------------------------------------------------------------------------





% ANALYSIS ----------------------------------------------------------------------------------------%

% DECODED values

uphz = phz.copy;
uphz.data = unwrap(phz.data);
uphz.data = discretize(uphz.data,-pi:pi:uphz(end));

cind = unitInclusion>=2 & posteriorMax>0.002;
sind =     stcm(:,1)   ...   
         & ~stcm(:,2)              ...
         & ~stcm(:,8)            ...          
         & ~stcm(:,7);            
sum(cind&sind)
ind = cind&sind&phz.data>-pi&phz.data<0;


%decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3]);    

cm = accumarray(uphz.data(ind),,[max(uphz.data),1],@circ_mean);
cs = accumarray(uphz.data(ind),exp(i*atan2(decError(:,2),decError(:,1))),[max(uphz.data),1],@mean);
cp = accumarray(uphz.data(ind),atan2(decError(:,2),decError(:,1)),[max(uphz.data),1],@PPC);
cn = accumarray(uphz.data(ind),1,[max(uphz.data),1],@sum);

figure,
%scatter(cm,abs(cs),5,cn,'filled');
scatter(cm(cn>3),cp(cn>3),5,cn(cn>3),'filled');
colormap('jet');
colorbar;
caxis([3,10]);

decError = multiprod(posEstCom(:,[1,2])-sq(xyz(:,'hcom',[1,2])),hvec(:,:,:),2,[2,3]);    
figure();
hold('on');
plot(atan2(decError(:,2),decError(:,1)))
plot(phz.data);
plot(unitInclusion/10-10)
plot(decError(:,2),decError(:,1))
%%

cind = unitInclusion>=4 & posteriorMax>0.001;% & stcm(:,7)~=7;%
 
sum(cind)


%decError = sqrt(sum((sq(xyz(ind,'hcom',[1,2]))-posEstCom(ind,[1,2])).^2,2)).*sign(decError);


hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [2,7,50,30];
% STATE - locomotion
ind = cind & stcm(:,1)==1&(stcm(:,5)==5|stcm(:,3)==3)&stcm(:,2)~=2;
% frontal
subplot(261),
decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
%decError = multiprod(posEstMax(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
hist2([[phz(ind,1);phz(ind,1)+pi*2],[decError(:,1);decError(:,1)]],...
      linspace(-pi,pi*3,30),linspace(-400,400,75)); 
title('position error')
ylabel('mm');
xlabel('theta phase');

% lateral
subplot(262),
hist2([[phz(ind,1);phz(ind,1)+pi*2],[decError(:,2);decError(:,2)]],...
      linspace(-pi,pi*3,30),linspace(-400,400,75)); 
title('position error')
ylabel('mm');
xlabel('theta phase');

% direction
subplot(263),
hist2([[phz(ind,1);phz(ind,1)+pi*2],atan2([decError(:,2);decError(:,2)],...
                                          [decError(:,1);decError(:,1)])],...
      linspace(-pi,pi*3,30),linspace(-pi,pi,75)); 
title('position error')
ylabel('yaw');
xlabel('theta phase');



% STATE - locomotion 
% CONDITION - head direction
% frontal
subplot(264),
decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3]);
%decError = multiprod(posEstMax(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
hist2([[phz(ind,1);phz(ind,1)+pi*2],[decError(:,1);decError(:,1)]],...
      linspace(-pi,pi*3,30),linspace(-400,400,75)); 
title('position error')
ylabel('mm');
xlabel('theta phase');

% lateral
subplot(265),
hist2([[phz(ind,1);phz(ind,1)+pi*2],[decError(:,2);decError(:,2)]],...
      linspace(-pi,pi*3,30),linspace(-400,400,75)); 
title('position error')
ylabel('mm');
xlabel('theta phase');

% direction
subplot(266),
hist2([[phz(ind,1);phz(ind,1)+pi*2],atan2([decError(:,2);decError(:,2)],...
                                          [decError(:,1);decError(:,1)])],...
      linspace(-pi,pi*3,30),linspace(-pi,pi,75)); 
title('position error')
ylabel('yaw');
xlabel('theta phase');


% STATE -Pause
ind = cind & stcm(:,1)==1&(stcm(:,4)==4|stcm(:,6)==6)&stcm(:,2)~=2;
% frontal 
decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
subplot(267),
hist2([[phz(ind,1);phz(ind,1)+pi*2],[decError(:,1);decError(:,1)]],...
      linspace(-pi,pi*3,30),linspace(-400,400,75)); 
title('position error')
ylabel('mm');
xlabel('theta phase');

% lateral
subplot(268),
hist2([[phz(ind,1);phz(ind,1)+pi*2],[decError(:,2);decError(:,2)]],...
      linspace(-pi,pi*3,30),linspace(-400,400,75)); 
title('position error')
ylabel('mm');
xlabel('theta phase');

% direction
subplot(269),
hist2([[phz(ind,1);phz(ind,1)+pi*2],atan2([decError(:,2);decError(:,2)],...
                                          [decError(:,1);decError(:,1)])],...
      linspace(-pi,pi*3,30),linspace(-pi,pi,75)); 
title('position error')
ylabel('yaw');
xlabel('theta phase');


% STATE - locomotion 
% CONDITION - head direction
% frontal
subplot(2,6,10),
decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3]);
hist2([[phz(ind,1);phz(ind,1)+pi*2],[decError(:,1);decError(:,1)]],...
      linspace(-pi,pi*3,30),linspace(-400,400,75)); 
title('position error')
ylabel('mm');
xlabel('theta phase');

% lateral
subplot(2,6,11),
hist2([[phz(ind,1);phz(ind,1)+pi*2],[decError(:,2);decError(:,2)]],...
      linspace(-pi,pi*3,30),linspace(-400,400,75)); 
title('lateral error')
ylabel('mm');
xlabel('theta phase');

% direction
subplot(2,6,12),
hist2([[phz(ind,1);phz(ind,1)+pi*2],atan2([decError(:,2);decError(:,2)],...
                                          [decError(:,1);decError(:,1)])],...
      linspace(-pi,pi*3,30),linspace(-pi,pi,75)); 
title('position error')
ylabel('yaw');
xlabel('theta phase');

ForAllSubplots('caxis([0,255])');
%ForAllSubplots('caxis([0,55])');
ForAllSubplots('hax=gca;hax.Units=''centimeters'';hax.Position(3:4)=[4,8];')

fax = axes('Position',[0,0,1,1],'Visible','off');

Lines(0.5,[],'k');
Lines([],0.5,'k');
xlim([0,1]);
ylim([0,1]);

ht = text(0.25,0.9, 'Head Movement Basis', 'HorizontalAlignment','center');
ht = text(0.75,0.9, 'Head Direction Basis','HorizontalAlignment','center');
ht = text(0.05,0.75,'Locomotion',          'HorizontalAlignment','center','Rotation',90);
ht = text(0.05,0.25,'Pause',               'HorizontalAlignment','center','Rotation',90);

print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_headingXthetaPhase_',mode,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_headingXthetaPhase_',mode,'_',Trial.filebase,'.png']]);



figure,plot(mind)

figure,psot(ep
decError = sqrt(sum((sq(xyz(ind,'hcom',[1,2]))-posEstCom(ind,[1,2])).^2,2));


figure();
imagesc(phasevals,shiftvals/40,eposTPRError)
axis('xy');




shiftvals = -round(sampleRate*2):round(sampleRate/50):round(sampleRate*2);
% SET time bins
slidingTimeWindowSize = round(sampleRate*1);
timevals = 1:size(xyz,1);
timevals(~nniz(xyz)) = [];
timevals = timevals(1:slidingTimeWindowSize:numel(timevals));
timeBinInds  = discretize([1:size(xyz,1)]',timevals);
badTimeBinInds = find(diff(timevals)>slidingTimeWindowSize)+1;
timeBinInds(ismember(timeBinInds,badTimeBinInds)) = 0;
% SET phase bins
phasevals = linspace(-pi,pi,13);
phaseBinInds = discretize(phz.data,phasevals);
% SET speed bins
speedvals = linspace(0,1.8,7);
speedBinInds = discretize(vxy(:,2),speedvals);

eposTPRErrorXY = nan([numel(timevals)-1,numel(shiftvals),numel(phasevals)-1,numel(speedvals)-1]);

eposTPRErrorXYTrj = nan([numel(timevals)-1,numel(shiftvals),numel(phasevals)-1,numel(speedvals)-1,2]);

meanTPRSpeed = nan([numel(timevals)-1,numel(phasevals)-1,numel(speedvals)-1]);
meanTPRPostMax = nan([numel(timevals)-1,numel(phasevals)-1,numel(speedvals)-1]);
meanTPRPhase = nan([numel(timevals)-1,numel(phasevals)-1,numel(speedvals)-1]);

for time = 1:numel(timevals)-2,
    tic
    disp([num2str(time),' of ',num2str(numel(timevals)-3)])

    timeWindow = time-2:time+2;
    
    timeWindowIndex = ismember(timeBinInds,timeWindow);

    txyz = sq(xyz(timeWindowIndex,'nose',[1,2]));
    tvxy = vxy(timeWindowIndex,2);
    ttvec = tvec(timeWindowIndex,:,:);    
    tposteriorMax = posteriorMax(timeWindowIndex);
    tphz = phz(timeWindowIndex);    
    
    tind =   unitInclusion(timeWindowIndex)>=3  ...
             & stcm(timeWindowIndex,1)==1       ...
             & stcm(timeWindowIndex,2)~=2       ...
             & timeBinInds((timeWindowIndex))==time...
             & tposteriorMax>0.002;

    tempPhaseBinInds = phaseBinInds(timeWindowIndex);
    tempSpeedBinInds = speedBinInds(timeWindowIndex);

    eposThetaPhaseRes = nan([size(tind,1),2]);
    eposThetaPhaseRes(timeBinInds(timeWindowIndex)==time,:) = posEstCom(timeBinInds==time,1:2);

    for speed = 1:numel(speedvals)-1,   
        for phase = 1:numel(phasevals)-1,   
            ind = tind  &  tempPhaseBinInds==phase  &  tempSpeedBinInds==speed;
            
            if sum(ind)>3
                meanTPRSpeed(time,phase,speed) = mean(tvxy(ind));
                meanTPRPostMax(time,phase,speed) = mean(tposteriorMax(ind));
                meanTPRPhase(time,phase,speed) = circ_mean(tphz(ind));
                tepos = eposThetaPhaseRes;
                tepos(~ind,:) = nan;


                                    
                for shift = 1:numel(shiftvals)
                    eposTPRErrorXYTrj(time,shift,phase,speed,:) = ...
                       mean(multiprod(circshift(tepos(:,1:2),shiftvals(shift))-txyz,ttvec,2,[2,3]),'omitnan');
                    
                    eposTPRErrorXY(time,shift,phase,speed) = ...
                        mean(sqrt(sum((txyz-circshift(tepos,shiftvals(shift))).^2,2)),'omitnan');
                end
            end
        end
    end
    toc
end



% $$$ figure();
% $$$ for speed = 1:numel(speedvals)-1,
% $$$ subplot2(2,numel(speedvals)-1,1,speed);    
% $$$     imagesc(phasevals,shiftvals/sampleRate,eposTPRErrorXY(:,:,speed));axis('xy');
% $$$ caxis([100,300])
% $$$ subplot2(2,numel(speedvals)-1,2,speed);
% $$$     imagesc(phasevals,shiftvals/sampleRate0,1./(bsxfun(@rdivide,eposTPRErrorXY(:,:,speed),min( ...
% $$$         eposTPRErrorXY(:,:,speed)))+eps));axis('xy');
% $$$     title(num2str(mean(10.^speedvals(speed:speed+1))))
% $$$ end



derror = eposTPRErrorXY;
%derror = abs(eposTPRErrorXYTrj(:,1));
% $$$ derror = eposTPRErrorBP;
% $$$ derror = eposTPRErrorHP;

mind = [];
minv = [];
%for time = 1:numel(timevals)-1,
for t = 1:size(derror,1)
    for speed = 1:numel(speedvals)-1,
        %[~,mind(t,speed,:)] = max(RectFilter(RectFilter(1./(bsxfun(@rdivide,sq(derror(t,:,:,speed)),min(sq(derror(t,:,:,speed))))+eps)',3,1)',3,1));
        [minv(t,speed,:),mind(t,speed,:)] = min(RectFilter(RectFilter(sq(derror(t,:,:,speed))',3,1)',3,1));
    end
end
mind(mind==0) = nan;
minv(minv==0) = nan;


shiftTimeBins = shiftvals(2:end)./sampleRate;

figure;
for s = 1:numel(speedvals)-1,   
    for p = 1:numel(phasevals)-1,   
        subplot2(numel(phasevals)-1,numel(speedvals)-1,p,s);

        ind = mind(:,s,p)~=1&mind(:,s,p)~=201;

        %hist2([shiftTimeBins(mind(ind,s,p))',minv(ind,s,p)],linspace(-1,2,30),linspace(40,400,6));
% $$$         hist2([shiftTimeBins(mind(ind,s,p))',eposTPRErrorXY(ind,100,p,s,1)],...
% $$$               linspace(-2,2,40),linspace(40,300,8));
% $$$         scatter([-2;-2;shiftTimeBins(mind(ind,s,p))'],[0;0;eposTPRErrorXY(ind,100,p,s,1)],...
% $$$               5,[0;400;minv(ind,s,p)],'filled');
% $$$         Lines(mean(shiftTimeBins(mind(ind,s,p)),'omitnan'),[],'m');
% $$$         xlim([-2,2]);
% $$$         ylim([0,250]);
% $$$         Lines(0,[],'k');

% $$$         bar(shiftTimeBins,histc(shiftTimeBins(mind(ind,s,p)),shiftTimeBins),'histc');
% $$$         Lines(mean(shiftTimeBins(mind(ind,s,p)),'omitnan'),[],'g');
% $$$         Lines(0,[],'m');
% $$$         xlim(shiftvals([1,end])./sampleRate);
% $$$         ylim([0,15]);
        

        scatter(log10(eposTPRErrorXY([ind;false;false],100,p,s)),log10(minv(ind,s,p)),10, ...
                shiftTimeBins(mind(ind,s,p)),'filled');
        xlim([1,2.7]);
        ylim([1,2.7]);
        
% $$$         scatter(log10([1;1;eposTPRErrorXY([ind;false;false],100,p,s,1)]),...
% $$$                 [1;1;log10(minv(ind,s,p))],5, ...
% $$$                 shiftTimeBins([1;150;mind(ind,s,p)]),'filled');
% $$$         line([0,3],[0,3],'Color','k');
% $$$         xlim([0,3]);
% $$$         ylim([0,3]);
        
% $$$         scatter(([0;0;eposTPRErrorXYTrj([ind;false;false],50,p,s,1)]),...
% $$$                 [1;1;log10(minv(ind,s,p))],5, ...
% $$$                 shiftTimeBins([1;150;mind(ind,s,p)]),'filled');
% $$$         xlim([-400,400]);
% $$$         ylim([1,2.7]);
        
        if s==1,ylabel(num2str(mean(phasevals([p:p+1]))));end
        if p==(numel(phasevals)-1),xlabel(num2str(speedvals([s:s+1])));end
            
    end
end
colormap('jet');
ForAllSubplots('caxis([-1,1])');
linkaxes(findobj(gcf,'Type','Axes'),'xy')




% $$$ 
% $$$ cm = reshape(permute(repmat(jet(numel(shiftvals)-1),[1,1,20]),[1,3,2]),[],3);
% $$$ figure();
% $$$ scatter(reshape(repmat(phasevals(1:end-1),[numel(shiftvals)-1,1]),[],1),shiftvals(reshape(mind,[],1))/40,20,cm,'filled');

figure();
imagesc(phasevals,speedvals,(mind-40)/40)
axis('xy');
caxis([-0.4,0.4]);

imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorXY(:,:,18),min(eposTPRErrorXY(:,:,18)))+eps));axis('xy');

figure();
subplot(131);
imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorXY,min(eposTPRErrorXY))+eps));axis('xy');
subplot(132);
imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorHP,min(eposTPRErrorHP))+eps));axis('xy');
subplot(133);
imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorBP,min(eposTPRErrorBP))+eps));axis('xy');







% OLD 
% $$$ 
% $$$ figure();
% $$$ for speed = 1:numel(speedvals)-1,
% $$$ subplot2(2,numel(speedvals)-1,1,speed);    
% $$$     imagesc(phasevals,shiftvals/40,eposTPRErrorXY(:,:,speed));axis('xy');
% $$$ caxis([100,300])
% $$$ subplot2(2,numel(speedvals)-1,2,speed);
% $$$     imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorXY(:,:,speed),min( ...
% $$$         eposTPRErrorXY(:,:,speed)))+eps));axis('xy');
% $$$     title(num2str(mean(10.^speedvals(speed:speed+1))))
% $$$ end
% $$$ 
% $$$ 
% $$$ mind = [];
% $$$ for speed = 1:numel(speedvals)-1,
% $$$ [~,mind(speed,:)] = max(RectFilter(RectFilter(1./(bsxfun(@rdivide,eposTPRErrorXY(:,:,speed),min(eposTPRErrorXY(:,:,speed)))+eps)',3,1)',3,1));
% $$$ end
% $$$ % $$$ 
% $$$ % $$$ cm = reshape(permute(repmat(jet(numel(shiftvals)-1),[1,1,20]),[1,3,2]),[],3);
% $$$ % $$$ figure();
% $$$ % $$$ scatter(reshape(repmat(phasevals(1:end-1),[numel(shiftvals)-1,1]),[],1),shiftvals(reshape(mind,[],1))/40,20,cm,'filled');
% $$$ 
% $$$ figure();
% $$$ imagesc(phasevals,speedvals,(mind-40)/40)
% $$$ axis('xy');
% $$$ caxis([-0.4,0.4]);
% $$$ 
% $$$ imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorXY(:,:,18),min(eposTPRErrorXY(:,:,18)))+eps));axis('xy');
% $$$ 
% $$$ figure();
% $$$ subplot(131);
% $$$ imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorXY,min(eposTPRErrorXY))+eps));axis('xy');
% $$$ subplot(132);
% $$$ imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorHP,min(eposTPRErrorHP))+eps));axis('xy');
% $$$ subplot(133);
% $$$ imagesc(phasevals,shiftvals/40,1./(bsxfun(@rdivide,eposTPRErrorBP,min(eposTPRErrorBP))+eps));axis('xy');






    


% CREATE video of decoding
% $$$ xbinsReal = discretize(xyz(:,'nose',1),binsE{1});
% $$$ ybinsReal = discretize(xyz(:,'nose',2),binsE{2});
% $$$ 
% $$$ xbinsReal = discretize(tpos(:,1),binsE{1});
% $$$ ybinsReal = discretize(tpos(:,2),binsE{2});
% $$$ hbinsReal = discretize(tpos(:,3),binsE{3});
% $$$ bbinsReal = discretize(tpos(:,4),binsE{4});
% $$$ 
% $$$ 
% $$$ [C,H] = bhv_contours();
% $$$ 
% $$$ H{1}.ZData(sqrt(sum(cat(3,H{1}.XData,H{1}.YData).^2,3))<0.2)=0;
% $$$ 
% $$$ 
% $$$ vidObj = VideoWriter('/storage/share/Projects/BehaviorPlaceCode/decode/posterior_example_immobile.avi','Uncompressed AVI');
% $$$ open(vidObj);
% $$$ 
% $$$ hfig = figure();
% $$$ hfig.Units = 'centimeters';
% $$$ pause(0.1);
% $$$ hfig.Position(3:4) = [30,20];
% $$$ pause(0.1);
% $$$ axPosition = axes('Units','centimeters','Position',[2,4,6,6]);
% $$$ hold(axPosition,'on');
% $$$ daspect(axPosition,[1,1,1]);
% $$$ axPosture = axes('Units','centimeters','Position',[2,12,6,6]);
% $$$ hold(axPosture,'on');
% $$$ daspect(axPosition,[1,1,1]);
% $$$ axBg = axes('Position',[0 0 1 1],'Visible','off');    
% $$$ 
% $$$ axPosX = axes('Units',         'centimeters',...
% $$$               'Position',      [13,16,14,3],...
% $$$               'NextPlot',      'add',...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on',...
% $$$               'XTickLabel',    {});
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,1)/10,'m','LineWidth',2);
% $$$ plot((1:tind)/decodingSampleRate,xyz(1:tind,'nose',1)/10,'g','LineWidth',2);
% $$$ xlim([310,375]);
% $$$ ylim([-50,50]);
% $$$ ylabel('X (cm)');
% $$$ haxLinesTime(1) = animatedline([310,310],[-50,50],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ axPosY = axes('Units',         'centimeters',...
% $$$               'Position',      [13,13,14,3],...
% $$$               'NextPlot',      'add',...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on',...
% $$$               'XTickLabel',    {});
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,2)/10,'m','LineWidth',2);
% $$$ plot((1:tind)/decodingSampleRate,xyz(1:tind,'nose',2)/10,'g','LineWidth',2);
% $$$ xlim([310,375]);
% $$$ ylim([-50,50]);
% $$$ ylabel('Y (cm)');
% $$$ haxLinesTime(2) = animatedline([310,310],[-50,50],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ 
% $$$ axPosY = axes('Units',         'centimeters',...
% $$$               'Position',      [13,10,14,3],...
% $$$               'NextPlot',      'add',...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on',...
% $$$               'XTickLabel',    {});
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,3),'m','LineWidth',2);
% $$$ plot((1:tind)/decodingSampleRate,fet(1:tind,1),'g','LineWidth',2);
% $$$ xlim([310,375]);
% $$$ ylim([-1.5,1.5]);
% $$$ ylabel({'Head','Pitch (rad)'});
% $$$ haxLinesTime(3) = animatedline([310,310],[-1.5,1.5],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ 
% $$$ axPosY = axes('Units',         'centimeters',...
% $$$               'Position',      [13,7,14,3],...
% $$$               'NextPlot',      'add',...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on',...
% $$$               'XTickLabel',    {});
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,4),'m','LineWidth',2);
% $$$ plot((1:tind)/decodingSampleRate,fet(1:tind,2),'g','LineWidth',2);
% $$$ xlim([310,375]);
% $$$ ylim([-1.5,1.5]);
% $$$ ylabel({'Body','Pitch (rad)'});
% $$$ haxLinesTime(4) = animatedline([310,310],[-1.5,1.5],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ axStc  = axes('Units',         'centimeters',...
% $$$               'Position',      [13, 4,14,3],...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on');
% $$$ plotSTC(stc,1,[],fliplr({'rear','loc','pause'}),fliplr('rbk'));
% $$$ axStc.YTickLabel = {'pause','loc','rear'};
% $$$ xlim([310,375]);
% $$$ xlabel('Time (s)');
% $$$ haxLinesTime(5) = animatedline([310,310],[1,4],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ limitsBound = [-50,50;-50,50;-1.5,1.5;-1.5,1.5,;1,4];
% $$$ 
% $$$ for t = 3100:1:3750
% $$$ 
% $$$     axes(axPosition);    
% $$$     if t == 3100,    
% $$$         imagesc(linspace(-50,50,100),linspace(-50,50,100),zeros([100,100]));
% $$$         axPositionCirc = circle(0,0,42,'--r');        
% $$$     end
% $$$ 
% $$$     if ~isnan(hbinsReal(t))&~isnan(bbinsReal(t)),
% $$$         axPositionIm = imagesc(binsE{1}/decodingSampleRate,binsE{2}/decodingSampleRate,sq(tE(:,:,hbinsReal(t),bbinsReal(t),t))');
% $$$         axPositionScXyz = scatter(xyz(t,'nose',1)/decodingSampleRate,xyz(t,'nose',2)/decodingSampleRate,40,'g','filled');
% $$$         axPositionScTpos = scatter(tpos(t,1)/decodingSampleRate,tpos(t,2)/decodingSampleRate,40,'m','filled');
% $$$     end    
% $$$     uistack(axPositionCirc,'top')
% $$$ 
% $$$     if t == 3100,        
% $$$ 
% $$$     end
% $$$     
% $$$     xlabel('X Position (cm)');
% $$$     ylabel('Y Position (cm)');
% $$$     title('Decoded Position');
% $$$     xlim([-50,50]);    
% $$$     ylim([-50,50]);    
% $$$     daspect([1,1,1]);
% $$$ 
% $$$     axes(axPosture);
% $$$     if t == 3100,    
% $$$         imagesc(linspace(-1.65,1.4,100),linspace(-1,1.75,100),zeros([100,100]));
% $$$     end
% $$$     
% $$$     if ~isnan(xbinsReal(t))&~isnan(ybinsReal(t)),
% $$$         axPostureIm = imagesc(binsE{3},binsE{4},sq(tE(xbinsReal(t),ybinsReal(t),:,:,t))');
% $$$         axPostureScFet = scatter(fet(t,1),fet(t,2),40,'g','filled');
% $$$         axPostureScTpos = scatter(tpos(t,3),tpos(t,4),40,'m','filled');
% $$$     end
% $$$     xlabel('Head-Body Pitch (rad)');
% $$$     ylabel('Body Pitch (rad)')
% $$$     title('Decoded Posture');
% $$$     xlim([-1.65,1.4]);
% $$$     ylim([-1,1.75]);
% $$$     
% $$$     if t == 3100,
% $$$         HC = cf(@(h) copyobj(h,gca), H);
% $$$         h = plot(0,0,'w');
% $$$         h.Visible = 'off';
% $$$         h = plot(0,0,'c');
% $$$         h.Visible = 'off';
% $$$         h = plot(0,0,'r');
% $$$         h.Visible = 'off';        
% $$$         lax = legend({'real','decode','low','high','rear'},'Location','eastoutside');
% $$$         lax.Position = lax.Position+[0.1,0,0,0];
% $$$         lax.Color = [0.6,0.6,0.6];
% $$$     end
% $$$     cf(@(h) uistack(h,'top'), HC);
% $$$     
% $$$ 
% $$$     
% $$$     axes(axBg);
% $$$     cla(axBg);
% $$$     text(0.14,0.1,['time: ',num2str(t/decodingSampleRate),' s']);
% $$$     
% $$$     af(@(h) clearpoints(h), haxLinesTime);
% $$$     for h = 1:numel(haxLinesTime),
% $$$          addpoints(haxLinesTime(h),[t/decodingSampleRate,t/decodingSampleRate],limitsBound(h,:));
% $$$     end
% $$$ 
% $$$     drawnow('update');
% $$$ 
% $$$     writeVideo(vidObj,getframe());
% $$$     
% $$$     % Clear axes
% $$$     delete([axPositionIm,axPositionScXyz,axPositionScTpos]);    
% $$$     delete([axPostureIm,axPostureScFet,axPostureScTpos]);
% $$$     
% $$$ end
% $$$ close(vidObj);

% $$$ 
% $$$ % DIAGNOSTIC plot of decoded and real positions in physical and behavioral spaces
% $$$ figure();
% $$$ subplot(511);
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,1)/10);
% $$$ hold('on');
% $$$ plot((1:tind)/decodingSampleRate,xyz(1:tind,'nose',1)/10);
% $$$ title('x pos');
% $$$ subplot(512);
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,2)/10);
% $$$ hold('on');
% $$$ plot((1:tind)/decodingSampleRate,xyz(1:tind,'nose',2)/10);
% $$$ title('y pos');
% $$$ subplot(513);
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,3));
% $$$ hold('on');
% $$$ plot((1:tind)/decodingSampleRate,fet(1:tind,1));
% $$$ title('head pitch');
% $$$ subplot(514);
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,4));
% $$$ hold('on');
% $$$ plot((1:tind)/decodingSampleRate,fet(1:tind,2));
% $$$ title('body pitch');
% $$$ subplot(515);
% $$$ plotSTC(stc,1,[],fliplr({'rear','hloc','hpause','lloc','lpause','theta'}),fliplr('rbcgkm'));
% $$$ xlabel('Time (s)');
% $$$ linkaxes(get(gcf,'Children'),'x');
% $$$ linkaxes(get(gcf,'Children'),'x');


% ENCAPSULATE decoded position in MTAData object

rpos = MTADfet.encapsulate(Trial,                                                             ... % MTATrial object
                           tpos,                                                              ... % Data
                           decodingSampleRate,                                                ... % Sample rate
                           ['bayesiandecoded_',pfs.parameters.type,'_',pfs.parameters.states],... % Name
                           ['bd',pfs.parameters.type],                                        ... % Label
                           'b'                                                                ... % key
);

ind = [stc{'a&t-m-s',xyz.sampleRate}];


% FIGURE - jpdf of physical and behavioral spaces
figure,
subplot(121);
out = hist2(rpos(ind,[1,2])/10,linspace([-50,50,50]),linspace([-50,50,50]));
imagesc(linspace([-50,50,50]),linspace([-50,50,50]),out'./sum(out(:)));
axis('xy');
caxis([0,max(caxis)/2]);
daspect([1,1,1]);
title({'Decoded Position',Trial.filebase});
xlabel('cm');
ylabel('cm');
if strcmp(pfstype,'xyhb');
    subplot(122);
    out = hist2(rpos(ind,[3,4]),linspace(-1.5,1,50),linspace(-1,1.8,50));
    imagesc(linspace(-1.5,1,50),linspace(-1,1.8,50),out'./sum(out(:)));
    axis('xy');
    caxis([0,max(caxis)/2]);
    daspect([1,1,1]);
end
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_jpdfs_pos_a_bhv_',pfstype,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_jpdfs_pos_a_bhv_',pfstype,'_',Trial.filebase,'.png']]);



meanPopRate = sum(ufr(ind,:),2);
meanPopIncl = sum(ufr(ind,:)>1,2);

rposError = [sqrt((sq(xyz(ind,'nose',[1,2]))-rpos(ind,[1,2])).^2)];
spatialError = sqrt(sum((sq(xyz(ind,'nose',[1,2]))-rpos(ind,[1,2])).^2,2));
if strcmp(pfstype,'xyhb'),
    rposError = cat(2,rposError,sqrt((sq(fet(ind,[1,2]))-rpos(ind,[3,4])).^2));  
    bhvError = sqrt(sum((sq(fet(ind,[1,2]))-rpos(ind,[3,4])).^2,2));
else
    bhvError = [];
    bhvErrorCondPos = [];
end






% FIGURE - jpdf of phisical and behavioral errors
if strcmp(pfstype,'xyhb'),
figure();
hist2([spatialError/10,bhvError],linspace(0,20,25),linspace(0,1.4,25));
title({'Spatial Vs Behavioral Error',[Trial.filebase,' ',num2str(numel(unitSubset)),' units']});
xlabel('cm');
ylabel('rad');
hax = gca();
hax.Units = 'centimeters';
hax.Position = [hax.Position([1,2]),4,4];
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_jpdf_pos_x_bhv_',pfstype,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_jpdf_pos_x_bhv_',pfstype,'_',Trial.filebase,'.png']]);
end


% COMPUTE mean error for decoded components individually
shuffle = @(x) circshift(x,randi(size(x,1)));
ind = [stc{'a&t-m-s',xyz.sampleRate}];
meanPopRate = sum(ufr(ind,:),2);
txss = sq(xyz(ind,'nose',[1,2]));
tfss = fet(ind,[1,2]);
trps = rpos(ind,:);
meanError = mean(rposError);
% COMPUTE shuffled error for decoded components
meanErrorShuffled = nan([1000,numel(pfsBins)]);
bhvErrorShuffled = [];
mbes = [];
for i = 1:1000,
    sptErrorShuffled = sqrt((shuffle(txss)-trps(:,[1,2])).^2);
    if strcmp(pfstype,'xyhb'),
        bhvErrorShuffled = sqrt((shuffle(tfss)-trps(:,[3,4])).^2);
        mbes = mean(bhvErrorShuffled(meanPopRate>0.6,:));
    end
    meanErrorShuffled(i,:) = [mean(sptErrorShuffled(meanPopRate>0.6,:)), mbes];
end

% COMPUTE zscore of
for e =1:numel(meanError),
zsMSError(e) = (meanError(e)-mean(meanErrorShuffled(:,e)))/std(meanErrorShuffled(:,e));
end



% FIGURE - mean error conditiond on position
posBins = linspace(-50,50,25);
xBins = discretize(xyz(ind,'nose',1)./10,posBins);
yBins = discretize(xyz(ind,'nose',2)./10,posBins);
aind = xBins>0&yBins>0;
figure();
subplot(121);
spatialErrorCondPos = accumarray([xBins(aind),yBins(aind)],spatialError(aind)/10,[numel(posBins),numel(posBins)],@median);
spatialErrorCondPos(spatialErrorCondPos==0) = nan;
imagescnan({posBins,posBins,spatialErrorCondPos'},[0,20],'linear',true,'colorMap',@jet);axis('xy');
title({'Mean Decoded Spatial Error','Conditioned on Position'});
if strcmp(pfstype,'xyhb'),
    subplot(122);
    bhvErrorCondPos = accumarray([xBins(aind),yBins(aind)],bhvError(aind),[numel(posBins),numel(posBins)],@median);
    bhvErrorCondPos(bhvErrorCondPos==0) = nan;
    imagescnan({posBins,posBins,bhvErrorCondPos'},[0,1],'linear',true,'colorMap',@jet);axis('xy');
    title({'Mean Decoded Pitch Error','Conditioned on Position'});
end
hax = findobj(gcf,'Type','Axes');
af(@(h) set(h,'Units','centimeters'),  hax);
hax(1).Position = [hax(1).Position([1,2]),4,4];
hax(2).Position = [hax(2).Position([1,2]),0.5,4];
if strcmp(pfstype,'xyhb'),
    hax(3).Position = [hax(3).Position([1,2]),4,4];
    hax(4).Position = [hax(4).Position([1,2]),0.5,4];
end
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_pos_a_bhv_Cond_pos',pfstype,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_pos_a_bhv_Cond_pos',pfstype,'_',Trial.filebase,'.png']]);



% FIGURE - Independence of errors
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1,1,14,16];

subplot(423);hist2(mud([rposError(:,[1,2])]),50,50);title('x vs y');
if strcmp(pfstype,'xyhb'),
    subplot(421);hist2(mud([spatialError,bhvError]),50,50);title('xy vs hb');
    subplot(424);hist2(mud([rposError(:,[3,4])]),50,50);title('h vs b');
    subplot(425);hist2(mud([rposError(:,[1,3])]),50,50);title('x vs h');
    subplot(426);hist2(mud([rposError(:,[2,4])]),50,50);title('y vs b');
    subplot(427);hist2(mud([rposError(:,[1,4])]),50,50);title('x vs b');
    subplot(428);hist2(mud([rposError(:,[2,3])]),50,50);title('y vs h');
end
hax = findobj(gcf,'Type','Axes');
af(@(h) set(h,'Units','centimeters'),  hax);
af(@(h) set(h,'Position',[h.Position(1:2),2,2]),  hax);
axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
text(0.5,0.9,{Trial.filebase,['units: ',num2str(numel(unitSubset))]});
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_independence_',pfstype,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_independence_',pfstype,'_',Trial.filebase,'.png']]);





save(fullfile(MTA_PROJECT_PATH,'analysis',['bhv_decode_',pfstype,'_',Trial.filebase,'.mat']),...
     'ind',                         ...
     'rposError',                   ...
     'spatialError',                ...
     'bhvError',                    ...
     'spatialErrorCondPos',         ...     
     'bhvErrorCondPos',             ...
     'posBins',                     ...
     'meanError',                   ...
     'meanErrorShuffled',           ...
     'zsMSError',                   ...
     'meanPopRate',                 ...
     'meanPopIncl',                 ...
     'unitSubset'                   ...
);

end


clear('ds');
trialSubset = [1:4,6,7,17,18,20:23];
pfstypes = {'xy','xyhb'};
for trialIndex =  1:numel(trialSubset);
    Trial = Trials{trialSubset(trialIndex)}; 
    for typeIndex = 1:numel(pfstypes)
        ds(trialIndex,typeIndex) = load(fullfile(MTA_PROJECT_PATH,'analysis',...
                              ['bhv_decode_',pfstypes{typeIndex},'_',Trial.filebase,'.mat']));
    end
end

me2d = reshape([ds(:,1).meanError],[2,12])';
me4d = reshape([ds(:,2).meanError],[4,12])';

me2d = [ds(:,1).spatialError']';
me4d = reshape([ds(:,2).meanError],[4,12])';



me2dx = cell2mat(af(@(s) median(s.rposError(:,1)), ds(:,1)));
me4dx = cell2mat(af(@(s) median(s.rposError(:,1)), ds(:,2)));

me2dy = cell2mat(af(@(s) median(s.rposError(:,2)), ds(:,1)));
me4dy = cell2mat(af(@(s) median(s.rposError(:,2)), ds(:,2)));

me2s = cell2mat(af(@(s) median(s.spatialError(:,1)), ds(:,1)));
me4s = cell2mat(af(@(s) median(s.spatialError(:,1)), ds(:,2)));

me4p = cell2mat(af(@(s) median(s.bhvError(:,1)), ds(:,2)));

nunits = cell2mat(af(@(s) numel(s.unitSubset(:)), ds(:,1)));

figure();
subplot
hold('on');
scatter(nunits,me2dx,20,'b','filled')
scatter(nunits,me4dx,20,'r','filled')

figure();
hold('on');
scatter(nunits,me2dy,20,'b','filled')
scatter(nunits,me4dy,20,'r','filled')



figure();
subplot(121);
hold('on');
scatter(nunits,me2s/10,20,'b','filled');
scatter(nunits,me4s/10,20,'r','filled');
xlabel('Number of Units');
ylabel('Distance (cm)');
title('Median Spatial Error');
legend({'2d place fields','4d place fields'},'location','eastoutside');
grid('on')
set(gca,'XTick',[0,25,50,75,100]);
set(gca,'YTick',[4,5,6,7,8,9,10]);

subplot(122);
scatter(nunits,me4p,20,'r','filled')
xlabel('Number of Units');
ylabel('Distance (rad)');
title('Median Behavioral Error');
legend({'4d place fields'},'location','eastoutside');
grid('on');
set(gca,'XTick',[0,25,50,75,100]);

hax = findobj(gcf,'Type','Axes');
af(@(h) set(h,'Units','centimeters'),  hax);
af(@(h) set(h,'Position',[h.Position(1:2),3,3]),  hax);

print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_regress_nunits_x_error.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_regress_nunits_x_error.png']]);




% 20180725 ----------------------------------------------
trialIndex = 20;
Trial = Trials{trialIndex};
stc = Trial.load('stc','msnn_ppsvd_raux');
unitSubset = units{trialIndex};

xyz = preproc_xyz(Trial,'trb');
xyz.resample(250);

vxy = vel(filter(copy(xyz),'ButFilter',3,2.5,'low'),{'spine_lower','head_back'});

states = {'rear','hloc','hpause','lloc','lpause','groom','sit','theta','spw'};

drz = xyz.copy();
drz.data = compute_drz(Trial,unitSubset,pfs,'feature',xyz.copy());

lfp = Trial.load('lfp',[68,72,76,82]);

phz = phase(resample(copy(lfp),xyz),[5,12]);

spk = Trial.load('spk',xyz.sampleRate,[],unitSubset);

ufrAll = Trial.load('ufr',xyz,[],[],0.04);

unitsInt = select_units(Trial,'int');
ufrInt = Trial.load('ufr',xyz,[],unitsInt,0.04);

int = Trial.load('spk',xyz.sampleRate,'',unitsInt);

figure
for unit = 1:numel(unitsInt),
    subplotfit(unit,20);
    rose(phz(int(unitsInt(unit)),1));
end



figure,imagesc(log10(sq(ys(:,21,:)))');

ufrIntObj = ufrInt.copy();
ufrIntObj.data = mean(ufrInt(:,[5,6,10]),2);
uiphz = ufrIntObj.phase([5,12]);

ufrAllObj = ufrAll.copy();
ufrAllObj.data = mean(ufrAll.data,2);


figure,
plot(mean(ufrInt(:,[5,6,10]),2));
hold('on');
plot(phz(:,2)*10+50);
plot(uiphz(:)*10+50);


ufr = Trial.ufr.copy;
ufr = ufr.create(Trial,xyz,'',unitSubset,0.04,true,'gauss');

xyn = xyz.copy();
xyn.data = sq(xyz(:,'nose',[1,2]));
xyl = xyn.copy();;
xyl.resample(lfp);

specArgs = struct('nFFT',2^11,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,40]);
   
[ys,fs,ts] = fet_spec(Trial,lfp,[],[],[],specArgs);

lys = ys.copy;
lys.data = mean(log10(ys(:,fs<15,2)),2);
lys.resample(xyz);

lrs = ys.copy;
lrs.data = mean(log10(ys(:,5<fs&fs<12,3)),2)./mean(log10(ys(:,5<fs&fs<12,4)),2);
lrs.resample(xyz);

uind = sum(ufr.data-eps,2)~=0&nniz(xyz);

tE = decode_bayesian_poisson(smap,(ufr(:,:)'));%,'bins',bins,'baseline',0.001);

%tE(:,:,uind) = tE;
%tE(:,:,~uind) = 0;



% ESTIMATE positions from posterior
epos = nan([size(tE,ndims(tE)),ndims(tE)-1]);
binsE = {};
for i = 1:numel(pfs.adata.bins),  binsE{i} = pfsBins{i}(grows{i});  end
gbins = cell([1,numel(pfsBins)]); 
[gbins{:}] = ndgrid(binsE{:});
gbins = cat(numel(gbins)+1,gbins{:});
ss = substruct('()',[repmat({':'},[1,ndims(gbins)-1]),{1}]);

for tind = 1:size(tE,ndims(tE)),
    ss.subs{end} = tind;
    epos(tind,:) = sum(reshape(gbins.*repmat(subsref(tE,ss),[ones([1,ndims(gbins)-1]),...
                        size(gbins,ndims(gbins))]),...
                               [],size(gbins,ndims(gbins))));
end

% $$$ tE = RectFilter(tE,3,1);
% $$$ tE = permute(RectFilter(permute(tE,[2,1,3]),3,1),[2,1,3]);

tpos = nan([size(tE,ndims(tE)),ndims(tE)-1]);
apos = nan([size(tE,ndims(tE)),1]);
for tind = 1:size(tE,ndims(tE)),
    tbin = LocalMinima2(-tE(:,:,tind),0,30);
    if ~isempty(tbin),
        tpos(tind,:) = [gpfsBins{1}(tbin(1)),gpfsBins{2}(tbin(2))];
        apos(tind)   = tE(tbin(1),tbin(2),tind);
    end
end



stcm = stc2mat(stc,xyz,states);

% $$$ timeSpec = [1:size(ys,1)]./ys.sampleRate;
% $$$ cmapLims = [1,2.5;1,2.5;1,3;1,3.5];

figure();

sp = [];
for s = 1:4,
    sp(s+1) = subplot(6,1,s+1);
    imagesc(ts,fs,log10(ys(:,:,s))');
    axis('xy');
    caxis(cmapLims(s,:));
    colormap(gca,'jet');
end
sp(6) = subplot(6,1,6);
haxSTS = plotSTC(Trial.stc,1,[],fliplr(states),fliplr('rbcgkmy'));
haxSTS.YTickLabel = fliplr(states);
xlabel('Time (s)');

linkaxes(get(gcf,'Children'),'x');

sp(1) = subplot(611);
hold('on');
position = animatedline();
position.Marker = '*';
position.MarkerSize = 20;
position.MarkerFaceColor  = 'm';
position.MarkerEdgeColor  = 'm';

predicted = animatedline();
predicted.Marker = '*';
predicted.MarkerSize = 20;
predicted.MarkerFaceColor  = 'm';
predicted.MarkerEdgeColor  = 'm';

for i = 5266163:2:5286163,
    if apos(i)<0.001, continue;end
    
    posterior = imagesc(pfsBins{1}(3:end-2),...
                pfsBins{2}(3:end-2),...
                tE(:,:,i)');
    
    position.clearpoints();
    position.addpoints(xyn(i,1),xyn(i,2),10);
    predicted.clearpoints();
    predicted.addpoints(tpos(i,1),tpos(i,2),10);
    daspect([1,1,1]);
    
    t = NearestNeighbour(timeSpec,stper(i));

    xlim(sp(2),[t-5,t+5]);
    drawnow();
    
    delete(posterior);
end


unitInclusion = sum(logical(ufr(:,:)-eps),2);
statePeriods = any(stcm==7|stcm==7,2);
statePeriods = any(stcm==6|stcm==6,2);

%positionError = sqrt(sum((tpos-xyls).^2,2));
positionError = sqrt(sum((epos-sq(xyz(:,'nose',[1,2]))).^2,2));
lowFreqMeanPower =lys(:);
thetaRatioRadLM =lrs(:);




ind = unitInclusion>=1&round(tpos(:,1)/5)~=0&round(tpos(:,2)/5)~=0&statePeriods&nniz(tpos)&apos>0.001;
figure,
subplot(121);
hist2([lowFreqMeanPower(ind),positionError(ind)],linspace(1,2.5,50),linspace(0,800,50));
colormap('jet');
caxis([0,500])
subplot(122);
hist2([thetaRatioRadLM(ind),positionError(ind)],linspace(0.6,1.1,50),linspace(0,800,50));
colormap('jet');
caxis([0,500])



%% segmentation of immobility sequences ------------------------------------------------------------

mufr = unity(mean(ufr.data,2));



hfig = figure();
hfig.CurrentCharacter = 's';
sp = gobjects([1,2]);;
sp(1) = subplot(221);
plot([1:numel(mufr)]./xyz.sampleRate,mufr);
hold('on');
plot([1:numel(mufr)]./xyz.sampleRate,apos);
plot([1:numel(mufr)]./xyz.sampleRate,unity(ufrAllObj.data));
sp(2) = subplot(223);
haxSTS = plotSTC(stc,1,[],fliplr(states),fliplr('rbcgkmy'));
haxSTS.YTickLabel = fliplr(states);
xlabel('Time (s)');
linkaxes(get(gcf,'Children'),'x');


sp(3) = subplot(222);
while ~strcmp(hfig.CurrentCharacter,'q')
    if diff(sp(1).XLim).*xyz.sampleRate<10000,
        cla(sp(3));
        segInd = round(sp(1).XLim.*xyz.sampleRate);
        segInd(segInd<1) = 1;
        segInd(segInd>numel(apos)) = numel(apos);    
        segInd = segInd(1):segInd(2);
        segInd = segInd(apos(segInd)>0.1);
        if ~isempty(segInd),
            axes(sp(3));
            imagesc(gpfsBins{:},sum(tE(:,:,segInd),3)');
            hold('on');
            plot(tpos(segInd,1),tpos(segInd,2),'-g');
            plot(epos(segInd,1),epos(segInd,2),'-c');            
            scatter(tpos(segInd(1),1),tpos(segInd(1),2),10,'m','filled');
            plot(xyz(segInd(1):segInd(end),'nose',1),xyz(segInd(1):segInd(end),'nose',2),'.','MarkerSize',10);
            axis('xy');
        end
        
        drawnow();
        waitforbuttonpress();
    else
        waitforbuttonpress();
        continue    
    end
end


% Pfeiffer, Foster 2013
% Criteria for spw decoding events
% normalized mean Population must > 1 Hz
% no sequence may jump more than 10cm
%
% output vars
% mean observed speed of each segment
% mean speed of segment


tper = [stc{'theta',xyz.sampleRate}];
sequencePeriods = ThreshCross(mufr,1,5);
% REMOVE theta periods
sequencePeriods(WithinRanges(mean(sequencePeriods,2),tper.data),:) = [];


seq.sampleRate = xyz.sampleRate;
seq.periods = sequencePeriods;

seq.duration = nan([size(seq.periods,1),1]);
seq.obsDistance = nan([size(seq.periods,1),1]);
seq.decDistance = nan([size(seq.periods,1),1]);
seq.obsDiffMean = nan([size(seq.periods,1),1]);
seq.decDiffMean = nan([size(seq.periods,1),1]);
seq.obsDiffStd = nan([size(seq.periods,1),1]);
seq.decDiffStd = nan([size(seq.periods,1),1]);
seq.obsDiffMax = nan([size(seq.periods,1),1]);
seq.decDiffMax = nan([size(seq.periods,1),1]);
seq.meanError = nan([size(seq.periods,1),1]);

seq.obsPathDistMean = nan([size(seq.periods,1),1]);
seq.decPathDistMean = nan([size(seq.periods,1),1]);
seq.obsPathDistStd = nan([size(seq.periods,1),1]);
seq.decPathDistStd = nan([size(seq.periods,1),1]);
seq.obsPathAngMean = nan([size(seq.periods,1),1]);
seq.decPathAngMean = nan([size(seq.periods,1),1]);
seq.obsPathAngStd = nan([size(seq.periods,1),1]);
seq.decPathAngStd = nan([size(seq.periods,1),1]);

seq.obsPathAngPPC = nan([size(seq.periods,1),1]);
seq.decPathAngPPC = nan([size(seq.periods,1),1]);


 
dvxy = vxy.data;
dxyz = xyz(:,'nose',[1,2]);
depos = permute(epos,[1,3,2]);
for p = 1:size(seq.periods,1),
    segInd = seq.periods(p,1):seq.periods(p,2);
    segInd = segInd(apos(segInd)>0.1);

    if numel(segInd)>5,
        seq.duration(p) = diff(segInd([1,end]));
        seq.obsDistance(p) = sqrt(sum(diff(dxyz(segInd([1,end]),1,:)).^2,3));
        seq.decDistance(p) = sqrt(sum(diff(depos(segInd([1,end]),1,:)).^2,3));

        mdiffXYZ = repmat(dxyz(segInd,1,:),[1,numel(segInd),1])-repmat(permute(dxyz(segInd,1,:),[2,1,3]),[numel(segInd),1,1]);    
        mdiffXYZ = repmat(dxyz(segInd,1,:),[1,numel(segInd),1])-repmat(permute(dxyz(segInd,1,:),[2,1,3]),[numel(segInd),1,1]);    
        mdistXYZ = sqrt(sum(mdiffXYZ.^2,3));

        mdiffEpos = repmat(depos(segInd,1,:),[1,numel(segInd),1])-repmat(permute(depos(segInd,1,:),[2,1,3]),[numel(segInd),1,1]);    
        mdiffEpos = repmat(depos(segInd,1,:),[1,numel(segInd),1])-repmat(permute(depos(segInd,1,:),[2,1,3]),[numel(segInd),1,1]);    
        mdistEpos = sqrt(sum(mdiffEpos.^2,3));    

        seq.obsPathDistMean(p) = mean(diag(mdistXYZ,1));
        seq.decPathDistMean(p) = mean(diag(mdistEpos,1));

        seq.obsPathDistStd(p) = std(diag(mdistXYZ,1));
        seq.decPathDistStd(p) = std(diag(mdistEpos,1));
        

        obsAng = atan2(diag(mdiffXYZ(:,:,1),1),diag(mdiffXYZ(:,:,2),1));
        decAng = atan2(diag(mdiffEpos(:,:,1),1),diag(mdiffEpos(:,:,2),1));        
        
        seq.obsPathAngMean(p) = circ_mean(circ_dist(obsAng(1:end-1),obsAng(2:end)));
        seq.decPathAngMean(p) = circ_mean(circ_dist(decAng(1:end-1),decAng(2:end)));
        
        seq.obsPathAngStd(p) = circ_std(circ_dist(obsAng(1:end-1),obsAng(2:end)));
        seq.decPathAngStd(p) = circ_std(circ_dist(decAng(1:end-1),decAng(2:end)));

        seq.obsPathAngPPC(p) = PPC(circ_dist(obsAng(1:end-1),obsAng(2:end)));
        seq.decPathAngPPC(p) = PPC(circ_dist(decAng(1:end-1),decAng(2:end)));
        
        
        seq.obsDiffMean(p) = mean(mdistXYZ(mdistXYZ(:)~=0  ));
        seq.decDiffMean(p) = mean(mdistEpos(mdistEpos(:)~=0));
        
        seq.obsDiffStd(p) = std(mdistXYZ(mdistXYZ(:)~=0  ));
        seq.decDiffStd(p) = std(mdistEpos(mdistEpos(:)~=0));

        if sum(mdistXYZ(:)~=0)~=0,
            seq.obsDiffMax(p) = max(mdistXYZ(mdistXYZ(:)~=0  ));
        else
            seq.obsDiffMax(p) = 0;
        end
        seq.decDiffMax(p) = max(mdistEpos(mdistEpos(:)~=0));
        
        seq.meanError(p)   = mean(sqrt(sum([depos(segInd,1,:)-dxyz(segInd,1,:)].^2,3)));
    end
    
end


figure,
figure();plot(seq.obsPathAngPPC,log10(seq.obsDiffMean),'.')

figure();
hold('on');
plot(seq.obsPathAngPPC,log10(seq.obsPathDistMean),'.')
plot(seq.decPathAngPPC,log10(seq.decPathDistMean),'.')

figure();
subplot(121);
sp = plot(seq.obsPathAngPPC,seq.obsPathAngMean,'.');



subplot(122);
plot(seq.decPathAngPPC,seq.decPathAngMean,'.')



hfig = figure();
sp = gobjects([0,1]);
sp(end+1) = subplot(221);
hold('on');
scatter(seq.decPathAngPPC,log10(seq.decPathDistMean),8,seq.duration,'filled')
sp(end+1) = subplot(222);
hold('on');
scatter(seq.obsPathAngPPC,log10(seq.obsPathDistMean),8,seq.duration,'filled')
linkaxes(get(gcf,'Children'),'xy');
drawnow();
c = gobjects([1,numel(sp)]);;
xy = [0,0];
while isempty(hfig.CurrentCharacter)||hfig.CurrentCharacter~='q',
    waitforbuttonpress();
% REMOVE old unit marker
    delete(c);     
% GET current subplot index    
    axind = find(arrayfun(@isequal,sp,repmat([gca],size(sp))));
% GET xy position of currsor on left mouse button down within current subplot    
    xy = sp(axind).CurrentPoint(1,1:2);
% GET axes Data
    axData = [sp(axind).Children.XData',sp(axind).Children.YData'];
% FIND closest point to currsor on left mouse button down
    [~,mind]=min(sqrt(sum(bsxfun(@minus,xy,axData).^2,2)));
% PLOT Posterior
    segInd = seq.periods(mind,1):seq.periods(mind,2);
    segInd = segInd(apos(segInd)>0.1);
    subplot(223);
    imagesc(gpfsBins{:},sum(tE(:,:,segInd),3)');
    caxis([0,0.75]);
    colormap('hot');
    axis('xy');
    hold('on');    
    plot(epos(segInd,1),epos(segInd,2),'m','LineWidth',2);
    plot(epos(segInd(1),1),epos(segInd(1),2),'c*','MarkerSize',10);    
    plot(xyn(segInd,1), xyn(segInd,2),'g*','MarkerSize',3);
    quiver(xyn(segInd,1), xyn(segInd,2),...
             cos(ang(segInd,'head_back','head_front',1))*100,...
             sin(ang(segInd,'head_back','head_front',1))*100,0,'Color',[1,1,1]);
    axes(sp(axind));
% HIGHLIGHT selected unit on current subplot
    c(axind) = circle(axData(mind,1),axData(mind,2),0.25,'g');
    for a = find(~ismember(1:numel(sp),axind))
        axes(sp(a));        
        c(a) = circle(sp(a).Children.XData(mind),sp(a).Children.YData(mind),0.25,'k');
    end
end




figure();
plot(seq.obsPathAngPPC,log10(seq.decDiff),'.')
plot(seq.decPathAngPPC,log10(seq.decDiffMax),'.')

figure,
hist2(log10(abs([seq.decMeanDiff,seq.obsDistance])+1e-4),100,100)
Lines([],log10(200),'m');

% Sequence analysis 

% Types of sequences
% THETA 
% PROSPECTIVE 
% SPW 
% 
% Theta Sequence: 
% The hippocampus represents the environment as a squence of action potentials of place cells arranged
% temoporally by the spatial order of the placefields with respect to the rats trajectory. Each sequence
% can last between 
% 
% Filter the sequences based on continuity 





% effect of acceleration sign on decoding head project
eds = linspace(-300,300,500);
phzBins = linspace(-pi,pi,8);
figure();
for m = 1:2,
    for p = 1:numel(phzBins)-1,
        subplot2(numel(phzBins)-1,2,p,m);
        hold('on');
        ind = nniz(decError) & (sign([0;diff(fvxy(:,m))])==-1) & fvxy(:,m)>5 ...
              & stcm(:,1) & ~stcm(:,2) ...
              & phz(:,1)>phzBins(p) & phz(:,1)<phzBins(p+1);
        hax = bar(eds,histc(decError(ind,1),eds),'histc');
        Lines(median(decError(ind,1),'omitnan'),[],'c');
        hax.EdgeColor = 'none';
        hax.FaceColor = 'c';
        hax.EdgeAlpha = 0.3;
        hax.FaceAlpha = 0.3;
        ind = nniz(decError) & (sign([0;diff(fvxy(:,m))])==1) &  fvxy(:,m)>5 ...
              & stcm(:,1) & ~stcm(:,2) ...
              & phz(:,1)>phzBins(p) & phz(:,1)<phzBins(p+1);
        hax = bar(eds,histc(decError(ind,1),eds),'histc');
        Lines(median(decError(ind,1),'omitnan'),[],'r');
        hax.EdgeColor = 'none';
        hax.FaceColor = 'r';
        hax.EdgeAlpha = 0.3;
        hax.FaceAlpha = 0.3;
    end
end



% DIAGNOSTIC figures ------------------------------------------------------------------------------

% $$$ figure,
% $$$ for u = 1:numel(pfs.data.clu),
% $$$ clf();
% $$$ srmap = pfs.plot(unitSubset(u),'mean',false,[],false,0.25,false,interpParPfsNdim,@jet,mazeMask);
% $$$ subplot2(8,6,1,[1,2]);
% $$$ plot(pft,unitSubset(u),'mean',true,[],false);
% $$$ title(num2str(u));
% $$$ rmax = pft.maxRate(unitSubset(u));
% $$$ for i = 0:5,
% $$$ for j = 0:6,    
% $$$     subplot2(8,6,j+2,i+1);
% $$$     imagescnan({pfsBins{1},pfsBins{2},srmap(:,:,i*3+20,j*3+15)'},[0,rmax],[],true,'colorMap',@jet);axis('xy');
% $$$ end
% $$$ end
% $$$ waitforbuttonpress();
% $$$ end
% $$$ mrt = pft.maxRate(unitSubset);



% TEST N-dimensional interpolation and masking of rate map
% $$$ 
% $$$ %srmap = srmap.*mask;
% $$$ figure,
% $$$ rmax = pft.maxRate(unitSubset(u));
% $$$ subplot(121);
% $$$ imagescnan({pfsBins{1},pfsBins{2},sq(srmap(:,:,25,25))'},[0,rmax],[],true,'colorMap',@jet);axis('xy');
% $$$ subplot(122);
% $$$ imagescnan({pfsBins{3},pfsBins{4},sq(srmap(30,20,:,:))'},[0,rmax],[],true,'colorMap',@jet);axis('xy');



