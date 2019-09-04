% MjgER2016_figure6
%
% DECODING position and behavior based on spatio-behavioral ratemaps
% 
% 
% 

%global AP
global MTA_PROJECT_PATH

MjgER2016_load_data();
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};
spikeWindow = 0.3;
sampleRate = 30;
sampleRateTPD = 250;
phzBins = -pi:pi/4:pi;
phzBinCenters = sum([phzBins(1:end-1);phzBins(2:end)])'./2
smoothingWeights = [250.^2,250.^2, 0.4.^2, 0.4.^2];
interpParPfs = struct('bins',{{linspace(-500,500,50),...
                               linspace(-500,500,50),...
                               linspace(  -2,  0.8,50),...
                               linspace(  -0.8,2,  50)}},...
                      'nanMaskThreshold', 0.1,...
                      'methodNanMap',     'linear',...
                      'methodRateMap',    'linear');



% SECTION 1 ----------------------------------------------------------------------------------------
%
%  Behavior rate maps
%
%      Generate the behavior rate maps which are computed from the data within the placefield during
%  theta periods
%
%  Joint physical-behavioral rate maps
%
%      Generate the spatial-behavioral rate maps which are computed from the data within the 
%  placefield during theta periods
%


% SET analysis args
MjgER2016_figure6_args('section 1');

% GET behavior field rate map
bfrm        = cf(@(t,u)   compute_bhv_ratemaps(t,u),                 Trials, units);
bfrmShuff   = cf(@(t,u)   compute_bhv_ratemaps_shuffled(t,u),        Trials, units);

% GET bhv ratemap erpPCA
[eigVecs, eigScrs, eigVars, unitSubsets, validDims, zrmMean, zrmStd] = ...
    compute_bhv_ratemaps_erpPCA(bfrm, units);

% GET bhv ratemap erpPCA scores
[fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS,fsrsMean,fsrsStd,rmapsShuffledMean,rmapsShuffled] = ...
    compute_bhv_ratemaps_erpPCA_scores(Trials,units,bfrm,bfrmShuff,eigVecs,validDims,unitSubsets,false);

% MAKE bhv ratemap mask
maskbhv = false(bfrm{1}.adata.binSizes');
maskbhv(validDims) = true;


% SELECT trial subset
tind = [3:5,17:23];
Trials = Trials(tind);
units  = units(tind);

stc  = cf(@(t)                                                                    ...
          t.load('stc','msnn_ppsvd_raux'),                                        ...
          Trials);

% GET xyhb placefields
dcxyhb.pfs  = cf(@(t,u)  compute_xyhb_ratemaps(t,u),  Trials, units);
% MAKE xy ratemap mask
maskcirc = create_tensor_mask(dcxyhb.pfs{1}.adata.bins(1:2));



% MAKE xyhb ratemap mask
% $$$ mask = repmat(maskcirc,[1,1,size(maskbhv)]).*repmat(permute(maskbhv,[3,4,1,2]),[size(maskcirc),1,1]);
%save(fullfile(MTA_PROJECT_PATH,'analysis','pfsXYHB_mask.mat'),'mask','-v7.3');
% LOAD xyhb ratemap mask
load(fullfile(MTA_PROJECT_PATH,'analysis','pfsXYHB_mask.mat'));


%%%<<< DIAGNOSTIC FIGURE
% $$$ rmp = plot(pfsXYHB{1},61,1,[],[],false);
% $$$ % CHECK mask is correct
% $$$ nx = pfsXYHB{1}.adata.binSizes(1);
% $$$ ny = pfsXYHB{1}.adata.binSizes(2);
% $$$ figure,
% $$$ nx = size(mask,1);
% $$$ ny = size(mask,2);
% $$$ sp = reshape(tight_subplot(ny,nx,0,0)',[ny,nx])';
% $$$ for x = 1:nx,
% $$$     for y = 1:ny,
% $$$         axes(sp(x,y));
% $$$         imagescnan(sq(mask(x,y,:,:))'.*sq(rmp(x,y,:,:))'); 
% $$$         axis('xy')
% $$$     end
% $$$ end
%%%>>>



%%%<<< DECODE xy 300ms window
dcxy.tag = 'xy_sr30_sw300_HighRes';
dcxy.sampleRate = 30;
dcxy.mask = maskcirc;
dcxy.smoothingWeights = [250.^2,250.^2];

% GET subject position objects
dcxy.xyz = cf(@(t)                                                              ...
              preproc_xyz( t, 'trb', sampleRate),                               ...
              Trials);

% GET spike time objects
dcxy.spk = cf(@(t,x,u)                                                          ...
              create( copy(t.spk), t, x.sampleRate, '', u, 'deburst'),          ...
              Trials, dcxy.xyz, units);

% GET unit firing rate objects
dcxy.ufr = cf(@(t,x,u,s)                                                        ...
              t.load('ufr',x,s,u,spikeWindow,true,'gauss'),                     ...
              Trials, dcxy.xyz, units, dcxy.spk);

% GET placefield ratemap objects
dcxy.pfs = cf(@(t,u)                                                            ...
              compute_ratemaps(t,u),                                            ...
              Trials, units);

% GET unit inclusion 
dcxy.uinc =  cf(@(u)                                                            ...
                sum(u.data>0.2,2),                                              ...
                dcxy.ufr);

% DECODE position from unit firing rate
[dcxy.com, dcxy.max, dcxy.sax, dcxy.post] =                                     ...
        cf(@(t,u,r,p)                                                           ...
           decode_ufr(t, u,                                                     ...
                      dcxy.sampleRate,                                          ...
                      r, p, [],                                                 ...
                      dcxy.mask,                                                ...
                      dcxy.smoothingWeights,                                    ...
                      'tag',dcxy.tag),                                          ...
           Trials, units, dcxy.ufr, dcxy.pfs);

% GET state matrix synced with dcxy.xyz objects
dcxy.stcm = cf(@(s,x)                                                           ...
               stc2mat(s,x,states),                                             ...
               stc, dcxy.xyz);
%%%>>>


%%%<<< DECODE xyhb 300ms window 
% LOAD decoded position based on low sampleRate and full theta phase
dcxyhb.tag = 'xy_sr30_HighRes';
dcxyhb.sampleRate = dcxy.sampleRate;
dcxyhb.mask = mask;
dcxyhb.smoothingWeights = [250.^2,250.^2, 0.4.^2, 0.4.^2];
dcxyhb.stcm = dcxy.stcm;
[dcxyhb.com, dcxyhb.max, dcxyhb.sax, dcxyhb.post] =                             ...
        cf(@(t,u,r,p)                                                           ...        
           decode_ufr(t, u,                                                     ...
                      sampleRate,                                               ...
                      r, p, [],                                                 ...
                      dcxyhb.mask, dcxyhb.smoothingWeights,                     ...
                      'tag',dcxyhb.tag),                                        ...
           Trials,units,dcxy.ufr,dcxyhb.pfs);
%%%>>>


dcxy.error.xy   = sqrt(sum((  dcxy.com{t}(:,[1,2])-sq(dcxy.xyz{t}(:,'nose',[1,2]))).^2,2));
dcxyhb.error.xy = sqrt(sum((dcxyhb.com{t}(:,[1,2])-sq(dcxy.xyz{t}(:,'nose',[1,2]))).^2,2));

t = 7;
ind =   dcxy.stcm{t}(:,1)==1 ...
      & any(dcxy.stcm{t}(:,3:6),2) ...
      & dcxy.uinc{t}>1 ...
      & dcxy.post{t}>0.0005 ...
      & dcxyhb.post{t}>0.0005;

figure,
    subplot(211);
        histogram(dcxy.error.xy(ind),...
                  0:10:500);
    subplot(212);
        histogram(dcxyhb.error.xy(ind),...
                  0:10:500);

        
[median(dcxy.error.xy(ind)),median(dcxyhb.error.xy(ind))]
[std(nonzeros(  dcxy.error.xy(ind).*double(0.5 < rand([sum(ind),1])))),...
 std(nonzeros(dcxyhb.error.xy(ind).*double(0.5 < rand([sum(ind),1]))))]

[p,h] = ranksum(  dcxy.error.xy(ind),...
                dcxyhb.error.xy(ind))

uincBin = 1:26;
xyeBin = 0:10:600;
uincInd  = discretize(dcxy.uinc{t}(ind),uincBin);
xyeInd   = discretize(  dcxy.error.xy(ind),xyeBin);
xyhbeInd = discretize(dcxyhb.error.xy(ind),xyeBin);

nind = nniz(xyeInd) & nniz(xyhbeInd) & nniz(uincInd);

figure();
    subplot(211);
        imagesc(xyeBin,                                                         ...
                uincBin,                                                        ...
                accumarray([xyeInd(nind),uincInd(nind)],                        ...
                           1,                                                   ...
                           [numel(xyeBin)-1,numel(uincBin)-1],                  ...
                           @sum)'                                               ...
                );
        axis('xy');
    subplot(212);
        imagesc(xyeBin,                                                         ...
                uincBin,                                                        ...
                accumarray([xyhbeInd(nind),uincInd(nind)],                        ...
                           1,                                                   ...
                           [numel(xyeBin)-1,numel(uincBin)-1],                  ...
                           @sum)'                                               ...
                );
        axis('xy');


dcxyhb.fet = cf(@(t)  fet_HB_pitchB(t,dcxyhb.sampleRate),  Trials);




% LOAD decoded position based on high sampleRate and binned by theta phase
decEstComTPD     = cell([1,numel(tind)]);
decEstMaxTPD     = cell([1,numel(tind)]);
decEstSaxTPD     = cell([1,numel(tind)]);
posteriorMaxTPD  = cell([1,numel(tind)]);
unitInclusionTPD = cell([1,numel(tind)]);
xyz  = cf(@(t)     resample(preproc_xyz(t,'trb'),sampleRateTPD),      Trials(tind));
for t = tind;
    Trials{t}.lfp.filename = [Trials{t}.name,'.lfp'];    
    tn = find(t==tind);
    try,   lfp = load(Trials{t},'lfp',sessionList(t).thetaRefGeneral);
    catch, lfp = load(Trials{t},'lfp',sessionList(t).thetaRefGeneral);
    end
    % PHZ : LFP phase restricted to the theta band [6-12Hz]
    phz = lfp.phase([6,12]);
    phzBinInds = discretize(phz.data,phzBins);
    spk  = create(copy(Trials{t}.spk),Trials{t},lfp.sampleRate,'',units{t},'deburst');
    for p = 1:numel(phzBins)-1
        spkt = copy(spk);
        spkt.clu = spkt.clu(phzBinInds(spkt.res)==p);
        spkt.res = round(spkt.res(phzBinInds(spkt.res)==p)./lfp.sampleRate.*sampleRateTPD);    
        spkt.sampleRate = sampleRateTPD;
        ufr = Trials{t}.load('ufr',xyz{tn},spkt,units{t},spikeWindow,true,'gauss');
        ufr.data =  ufr.data.*6;
        unitInclusionTPD{tn}(:,p) = sum(ufr.data>0.2,2);
        tag = ['xyhb_thpD',num2str(sampleRateTPD),'_',num2str(p),'o',num2str(numel(phzBins)-1)];
        [decEstComTPD{tn}(:,:,p),decEstMaxTPD{tn}(:,:,p),...
         decEstSaxTPD{tn}(:,:,p),posteriorMaxTPD{tn}(:,p)] = ...
            decode_ufr(Trials{t},units{t},sampleRateTPD,ufr,pfsXYHB{tn},...
                       [],mask,smoothingWeights,'tag',tag,'overwrite',false);
    end
end




% COMPILE egocentric variables 
fet  = cf(@(t)    fet_HB_pitchB(t,sampleRateTPD),                                    Trials(tind));
hvec = cf(@(x)    x(:,'head_front',[1,2])-x(:,'head_back',[1,2]),                    xyz          );
hvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                          hvec         );
hvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                         hvec         );
tvec = cf(@(x)    circshift(x(:,'hcom',[1,2]),-1)-circshift(x(:,'hcom',[1,2]),1),    xyz          );
tvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                          tvec         );
tvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                         tvec         );
stc  = cf(@(t)    t.load('stc','msnn_ppsvd_raux'),                                   Trials(tind));
stcm = cf(@(s,x)  stc2mat(s,x,states),                                               stc,   xyz   );

% COMPUTE error
dErr = decEstComTPD;
%dErr = decEstSaxTPD;
%dErr = decEstMaxTPD;
dErr = cf(@(e,x,h,f)  cat(2,...
                          sq(multiprod(permute(bsxfun(@minus,...
                                                  e(:,[1,2],:),...
                                                  sq(x(:,'hcom',[1,2]))),...
                                           [1,2,4,3]),...
                                   h(:,:,:),2,[2,3])),...
                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
          dErr,xyz,hvec,fet);



%%%<<< LOAD example data
t = 20;
Trial = Trials{t};
unitSubset = units{t};
exyz = resample(Trial.load('xyz','trb'),30);
%efet = fet_HB_pitchB(Trial,sampleRate);
efet = resample(copy(fet{tind==t}),sampleRate);
ets = [1:size(exyz,1)]./sampleRate;
estcm = stc2mat(Trial.load('stc'),exyz,states);
eind =   posteriorMax{t==tind} > 0.001 ...
       & unitInclusion{t==tind} > 3     ...
       & estcm(:,1)==1;
enanmask = ones([size(exyz,1),1]);
enanmask(~eind)=nan;

%%%>>>

% STARTFIG -----------------------------------------------------------------------------------------


eTS = [290,345];
%eTS = [530,600];
%eTS = [1000,1045];
eInds = round(eTS.*sampleRate)+1;
eInds = [eInds(1):eInds(2)]';

[hfig,fig,fax,sax] = set_figure_layout(figure(666006),'A4','portrait',[],1.5,1.5,0,0.2);




%%%<<< PLOT trajectory example in xy space
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(1,-fig.subplot.height, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*2,                        ...
                              fig.subplot.height*2],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(exyz(eInds,5,1),...
     exyz(eInds,5,2),...
     '-k','LineWidth',1);
plot(sq(decEstCom{t==tind}(eInds,1)).*enanmask(eInds),...
     sq(decEstCom{t==tind}(eInds,2)).*enanmask(eInds),...
     '-r','LineWidth',1);
xlim([-500,500]);
ylim([-500,500]);
sax(end).XTick = [];
sax(end).YTick = [];
%%%>>>


%%%<<< PLOT trajectory example in bhv space
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(3, -fig.subplot.height, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*2,                        ...
                              fig.subplot.height*2],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(efet(eInds,1),...
     efet(eInds,2),...
     '-k','LineWidth',1);
plot(sq(decEstCom{t==tind}(eInds,3)).*enanmask(eInds),...
     sq(decEstCom{t==tind}(eInds,4)).*enanmask(eInds),...
     '-r','LineWidth',1);
xlim(sax(end),[-1.8,0.6]);
ylim(sax(end),[-0.6,1.8]);
sax(end).XTick = [];
sax(end).YTick = [];
%%%>>>


%%%<<< PLOT x vs EST x
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(1, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*5,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(eInds),...
     exyz(eInds,5,1),...
     '-k','LineWidth',1);
plot(ets(eInds),...
     sq(decEstCom{t==tind}(eInds,1)).*enanmask(eInds),...
     '-r','LineWidth',1);
sax(end).XTick = [];
sax(end).YTick = [];
%%%>>>


%%%<<< PLOT y vs EST y
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*5,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(eInds),...
     exyz(eInds,5,2),...
     '-k','LineWidth',1);
plot(ets(eInds),...
     sq(decEstCom{t==tind}(eInds,2)).*enanmask(eInds),...
     '-r','LineWidth',1);
sax(end).XTick = [];
sax(end).YTick = [];
%%%>>>


%%%<<< PLOT hp vs EST hp
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(3, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*5,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(eInds),...
     efet(eInds,1),...
     '-k','LineWidth',1);
plot(ets(eInds),...
     sq(decEstCom{t==tind}(eInds,3)).*enanmask(eInds),...
     '-r','LineWidth',1);
sax(end).XTick = [];
sax(end).YTick = [];
ylim([-1.8,0.6]);
%%%>>>


%%%<<< PLOT bp vs EST bp
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(4, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*5,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(eInds),...
     efet(eInds,2),...
     '-k','LineWidth',1);
plot(ets(eInds),...
     sq(decEstCom{t==tind}(eInds,4)).*enanmask(eInds),...
     '-r','LineWidth',1);
sax(end).YTick = [];
ylim([-0.6,1.8]);
%%%>>>


%%%<<< COMPUTE error for TDP decoding
stid = [3,4,5,6];
indTDP  = cf(@(p,u,s)                                           ...
             all(p>0.0005,2)                                    ...
               & sum(double(u>=1),2)>6                          ...
               & s(:,1)==1                                      ...
               & any(logical(s(:,stid)),2),                     ...
             posteriorMaxTPD,unitInclusionTPD,stcm);
dcTDPfrontal = cell([1,10]);
dcTDPlateral = cell([1,10]);
for tn = 1:10,
    for j = 1:8;
        dcTDPfrontal{tn} = cat(2,                                                                ...
                               dcTDPfrontal{tn},                                                 ...
                               histc(sq(dErr{tn}(indTDP{tn},1,j)),linspace(-500,500,100)));
        dcTDPlateral{tn} = cat(2,                                                                ...
                               dcTDPlateral{tn},                                                 ...
                               histc(sq(dErr{tn}(indTDP{tn},2,j)),linspace(-500,500,100)));
    end
end
%%%>>>


%%%<<< PLOT frontal error for TDP decoding
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(5, -0.8, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
imagesc(linspace(-500,500,250),                                 ...
        circ_rad2ang([phzBinCenters;phzBinCenters+2*pi]),       ...
        repmat(sum(cat(3,dcTDPfrontal{:}),3),1,2)')
axis('xy');
xlim([-200,200]);
ylim([-60,420]);            
% MASK areas outside single cycle
patch([-300,-300,300,300],...
      [ 360, 420,420,360],...
      [0.2,0.2,0.2],...
      'EdgeAlpha',0,...
      'FaceAlpha',0.4);
patch([-300,-300,300,300],...
      [ -60,   0,  0,-60],...
      [0.2,0.2,0.2],...
      'EdgeAlpha',0,...
      'FaceAlpha',0.4);
sax(end).XTick = [-150,0,150];
sax(end).XTickLabels = [-15,0,15];
%%%>>>

%%%<<< PLOT lateral error as function of theta phase
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(5, -0.8, 2, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
imagesc(linspace(-500,500,250),                                 ...
        circ_rad2ang([phzBinCenters;phzBinCenters+2*pi]),       ...
        repmat(sum(cat(3,dcTDPlateral{:}),3),1,2)')
axis('xy');
xlim([-200,200]);
ylim([-60,420]);            
% MASK areas outside single cycle
patch([-300,-300,300,300],...
      [ 360, 420,420,360],...
      [0.2,0.2,0.2],...
      'EdgeAlpha',0,...
      'FaceAlpha',0.4);
patch([-300,-300,300,300],...
      [ -60,   0,  0,-60],...
      [0.2,0.2,0.2],...
      'EdgeAlpha',0,...
      'FaceAlpha',0.4);
sax(end).XTick = [-150,0,150];
sax(end).XTickLabels = [-15,0,15];
%%%>>>



%%%<<< supfig
figure();
for tn = 1:10
    subplot2(2,10,1,tn);
    imagesc(linspace(-500,500,250),                             ...
            [phzBinCenters;phzBinCenters+2*pi],                 ...
            repmat(dcTDPfrontal{tn},1,2)');
    axis('xy');
    xlim([-300,300]);    
    
    subplot2(2,10,2,tn);
    imagesc(linspace(-500,500,250),1:8,repmat(dcTDPlateral{tn},1,2)')
    axis('xy');
    xlim([-300,300]);
end
%%%>>>



% STOPFIG ------------------------------------------------------------------------------------------
    

ind  = cf(@(p,u,s) all(p>0.00005,2) & all(u>=2,2) & s(:,1)==1 & any(logical(s(:,stid)),2),...
          posteriorMaxTPD,unitInclusionTPD,stcm);


xD = linspace(-300,300,30);
yD = linspace(-300,300,30);
hD = linspace(-1.5,1.5,30);
bD = linspace(-1.5,1.5,30);

hCntXY = cf(@(e,i)  cat(3,...
                        hist2(e(i,[1,2],1),xD,yD),...
                        hist2(e(i,[1,2],2),xD,yD),...
                        hist2(e(i,[1,2],3),xD,yD),...
                        hist2(e(i,[1,2],4),xD,yD),...
                        hist2(e(i,[1,2],5),xD,yD),...
                        hist2(e(i,[1,2],6),xD,yD)),...
                        dErr,ind);
hCntXY = sum(cat(4,hCntXY{:}),4);

hCntHB = cf(@(e,i)  cat(3,...
                        hist2(e(i,[3,4],1),hD,bD),...
                        hist2(e(i,[3,4],2),hD,bD),...
                        hist2(e(i,[3,4],3),hD,bD),...
                        hist2(e(i,[3,4],4),hD,bD),...
                        hist2(e(i,[3,4],5),hD,bD),...
                        hist2(e(i,[3,4],6),hD,bD)),...
                        dErr,ind);
hCntHB = sum(cat(4,hCntHB{:}),4);

figure,
spp = gobjects([0,1]);
spb = gobjects([0,1]);
for p = 1:numel(phzBins)-1,
    spp(end+1) = subplot2(numel(phzBins)-1,2,p,1);    imagesc(xD,yD,hCntXY(:,:,p)');   axis('xy');
    spb(end+1) = subplot2(numel(phzBins)-1,2,p,2);    imagesc(hD,bD,hCntHB(:,:,p)');   axis('xy');
end
af(@(a) caxis(a,[0,120]),spp);
af(@(a) caxis(a,[0,400]),spb);



t = 3
ind{t} = ':';
figure,
subplot(4,1,1);plot([xyz{t}(ind{t},5,1),sq(dErr{t}(ind{t},1,:))]);    
subplot(4,1,2);plot([xyz{t}(ind{t},5,2),sq(dErr{t}(ind{t},2,:))]);        
subplot(4,1,3);plot([fet{t}(ind{t},1),  sq(dErr{t}(ind{t},3,:))]);
subplot(4,1,4);plot([fet{t}(ind{t},2),  sq(dErr{t}(ind{t},4,:))]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');




tind = 20;
Trial = Trials{tind};
unitSubset = units{tind};
xyz = preproc_xyz(Trial,'trb');
xyz.resample(sampleRate);
fet = fet_HB_pitchB(Trial,sampleRate);
hvec = xyz(:,'head_front',[1,2])-xyz(:,'head_back',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
% CONVERT MTAStateCollection into a state matrix 
stc = Trial.load('stc','msnn_ppsvd_raux');
stcm = stc2mat(stc,xyz,states);



dfet = decEstComTPD{3};
%dfet = decEstSaxTPD{1};
decError = zeros([size(xyz{t},1),size(decEstSax{j},2),numel(phzBins)-1]);
for p = 1:numel(phzBins)-1,
    decError(:,[1,2],p) = [multiprod(dfet(:,[1,2],p)-sq(xyz{t}(:,'hcom',[1,2])),hvec{t}(:,:,:),2,[2,3])];
end            
decError(:,[3,4],:) = [bsxfun(@minus,fet{t}(:,1),dfet(:,3,:)),bsxfun(@minus,fet{t}(:,2),dfet(:,4,:))];


j = 3
ind =   all(posteriorMaxTPD{j}>0.0001,2)...
      & all(unitInclusionTPD{j}>=4,2) ...
      & stcm(:,1)==1 & any(logical(stcm(:,[3,4,5,6])),2);
figure,
spp = gobjects([0,1]);
spb = gobjects([0,1]);
for p = 1:numel(phzBins)-1,
    spp(end+1) = subplot2(numel(phzBins)-1,2,p,1);
    hist2(decError(ind,[1,2],p),linspace(-300,300,50),linspace(-400,400,50));
    spb(end+1) = subplot2(numel(phzBins)-1,2,p,2);
    hist2(decError(ind,[3,4],p),linspace(-1.5,1.5,50),linspace(-1.5,1.5,50));
end
af(@(a) caxis(a,[0,50]),spp);
af(@(a) caxis(a,[0,150]),spb);



figure,
ind = all(posteriorMaxTPD{j}>0.000001,2)      & all(unitInclusionTPD{t}>=4,2);
subplot(4,1,1);plot([xyz{t}(ind,5,1),sq(dfet(ind,1,:))]);    
subplot(4,1,2);plot([xyz{t}(ind,5,2),sq(dfet(ind,2,:))]);        
subplot(4,1,3);plot([fet{t}(ind,1),  sq(dfet(ind,3,:))]);
subplot(4,1,4);plot([fet{t}(ind,2),  sq(dfet(ind,4,:))]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');



figure,
ind = posteriorMax{j}>0.0001;
subplot(5,1,1);plot([xyz(ind,5,1),sq(decEstCom{j}(ind,1))]);    
subplot(5,1,2);plot([xyz(ind,5,2),sq(decEstCom{j}(ind,2))]);        
subplot(5,1,3);plot([fet(ind,1),  sq(decEstCom{j}(ind,3))]);
subplot(5,1,4);plot([fet(ind,2),  sq(decEstCom{j}(ind,4))]);
subplot(5,1,5);plot([unitInclusion{j}(ind)]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');
