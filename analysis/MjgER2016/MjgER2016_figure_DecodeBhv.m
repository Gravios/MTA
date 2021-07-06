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
%MjgER2016_figure6_args('section 1');
configure_default_args();

% GET behavior field rate map
bfrm        = cf(@(t,u)   compute_bhv_ratemaps(t,u),                 Trials, units);
bfrmShuff   = cf(@(t,u)   compute_bhv_ratemaps_shuffled(t,u),        Trials, units);

% GET bhv ratemap erpPCA
[eigVecs, eigScrs, eigVars, unitSubsets, validDims, zrmMean, zrmStd] = ...
    compute_bhv_ratemaps_erpPCA(bfrm, units);

cluSessionMapSubset = cluSessionMap(unitSubsets,:);

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
numTrials = numel(Trials);

stc  = cf(@(t)                                                                    ...
          t.load('stc','msnn_ppsvd_raux'),                                        ...
          Trials);

% GET xyhb placefields
dc4.pfs  = cf(@(t,u)  compute_xyhb_ratemaps(t,u),  Trials, units);
% MAKE xy ratemap mask
maskcirc = create_tensor_mask(dc4.pfs{1}.adata.bins(1:2));



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

dc2.tag = 'xy_sr30_sw300_HighRes';
dc2.sampleRate = 30;
dc2.mask = maskcirc;
dc2.smoothingWeights = [250.^2,250.^2];

% GET subject position objects
dc2.xyz = cf(@(t)                                                               ...
              preproc_xyz( t, 'trb', sampleRate),                               ...
              Trials);

% GET spike time objects
dc2.spk = cf(@(t,x,u)                                                           ...
              create( copy(t.spk), t, x.sampleRate, '', u, 'deburst'),          ...
              Trials, dc2.xyz, units);

% GET unit firing rate objects
dc2.ufr = cf(@(t,x,u,s)                                                         ...
              t.load('ufr',x,s,u,spikeWindow,'gauss',true),                     ...
              Trials, dc2.xyz, units, dc2.spk);

% GET placefield ratemap objects
dc2.pfs = cf(@(t,u)                                                             ...
              compute_ratemaps(t,u),                                            ...
              Trials, units);

% GET unit inclusion 
dc2.uinc =  cf(@(u)                                                             ...
                sum(u.data>0.2,2),                                              ...
                dc2.ufr);

% $$$ % GENERATE behavioral class vector
% $$$ dc2.ubhv = {};
% $$$ for t = 1:numel(tind)
% $$$     for u = units{t}
% $$$         dc2.ubhv{t==tind} = zeros([size(dc2.ufr{t==tind}),3]);
% $$$         usind = find(ismember(cluSessionMapSubset,[t,u],'rows'));
% $$$         if ~isempty(usind)
% $$$             dc2.ubhv{t==tind}(:,u,:) = repmat(fsrcz(usind,1:3),size(dc2.ubhv{t==tind},1),1);
% $$$         end
% $$$     end
% $$$ end
% $$$ 


% DECODE position from unit firing rate
[dc2.com, dc2.max, dc2.sax, dc2.post] =                                         ...
        cf(@(t,u,r,p)                                                           ...
           decode_ufr(t, u,                                                     ...
                      dc2.sampleRate,                                           ...
                      r, p, [],                                                 ...
                      dc2.mask,                                                 ...
                      dc2.smoothingWeights,                                     ...
                      'tag',dc2.tag),                                           ...
           Trials, units, dc2.ufr, dc2.pfs);

% GET state matrix synced with dc2.xyz objects
dc2.stcm = cf(@(s,x)                                                            ...
               stc2mat(s,x,states),                                             ...
               stc, dc2.xyz);

%%%>>>


%%%<<< DECODE xyhb 300ms window 
% LOAD decoded position based on low sampleRate and full theta phase
dc4.tag = 'xy_sr30_HighRes';
dc4.sampleRate = dc2.sampleRate;
dc4.mask = mask;
dc4.smoothingWeights = [250.^2,250.^2, 0.4.^2, 0.4.^2];
dc4.stcm = dc2.stcm;

[dc4.com, dc4.max, dc4.sax, dc4.post] =                                         ...
        cf(@(t,u,r,p)                                                           ...        
           decode_ufr(t, u,                                                     ...
                      sampleRate,                                               ...
                      r, p, [],                                                 ...
                      dc4.mask, dc4.smoothingWeights,                           ...
                      'tag',dc4.tag),                                           ...
           Trials,units,dc2.ufr,dc4.pfs);
%%%>>>


%%%<<< COMPUTE errors
dc4.xyz = dc2.xyz;
dc4.fet = cf(@(t)  fet_HB_pitchB(t,dc4.sampleRate),  Trials);

dc2.error.xy   = cf(@(c,x)                                                      ...
                     sqrt( sum(( c(:,[1,2])-sq( x(:,'nose',[1,2]))).^2,2)),     ...
                     dc2.com, dc2.xyz);
dc4.error.xy = cf(@(c,x)                                                        ...
                     sqrt( sum(( c(:,[1,2])-sq( x(:,'nose',[1,2]))).^2,2)),     ...
                     dc4.com, dc2.xyz);
dc4.error.hp = cf(@(c,f)                                                        ...
                     c(:,3)-sq(f(:,1)),                                         ...
                     dc4.com, dc4.fet);
dc4.error.bp = cf(@(c,f)                                                        ...
                     c(:,4)-sq(f(:,2)),                                         ...
                     dc4.com,dc4.fet);
dc4.error.hb = cf(@(c,f)                                                        ...
                  sqrt( sum((c(:,[3,4])-f(:,[1,2])).^2,2)),                     ...
                  dc4.com,dc4.fet);
%%%>>>



uincBin = 0:26;
xyeBin = 0:2.5:300;
hpeBin = -1:0.05:1;
bpeBin = -1:0.05:1;
hbeBin = 0:0.01:1;
xyBin = -500:50:500;
hpBin = -2:0.1:0.5;
bpBin = -0.5:0.1:2;

% $$$ ddz        = cf(@(t,u,p)                                        ...
% $$$                 compute_ddz(t,u,p,'sampleRate',dc2.sampleRate), ...
% $$$                 Trials, units, dc2.pfs);
% $$$ 
% $$$ % for each trial
% $$$     %for each time point:
% $$$     %    identify the units within 20cm radius
% $$$     %    compute unit composition
% $$$     
% $$$ 
% $$$ cFSrC = bsxfun(@plus,FSrC(:,[1:3]),[0.5,1.5,1.5]);
% $$$ figure,
% $$$ plot3(cFSrC(:,1),cFSrC(:,2),cFSrC(:,3),'.');
% $$$ hold('on');
% $$$ line([0,3],[0,3],[0,3],'Color','r')
% $$$ 
% $$$ rvec = [1/sqrt(numel(cvec)).*ones(size(cvec))]';
% $$$ 
% $$$ zang = repmat({[]},size(dc2.xyz));
% $$$ for t = 1:numel(tind),
% $$$     ti = tind(t);
% $$$     for s = find(ind{t})',
% $$$         usub = units{t}(abs(ddz{t}(s,:))<200);
% $$$         zvec = cFSrC( ismember( cluSessionMapSubset,               ...
% $$$                                 [repmat(ti,[numel(usub),1]),usub'],...
% $$$                                 'rows'),                           ...
% $$$                       1:3);
% $$$         if size(zvec,1) < 3,
% $$$             zang{t}(end+1) = 0;
% $$$             continue;
% $$$         else,
% $$$             zvec = sum(bsxfun(@rdivide,zvec,sqrt(sum(zvec.^2,2))));
% $$$             zang{t}(end+1) = acos((zvec * rvec)./sqrt( sum( zvec.^2, 2)));            
% $$$         end
% $$$     end
% $$$ end



%%%<<< BIN errors
uincInd    = cf(@(u) discretize(u,uincBin),                    dc2.uinc);
xyInd      = cf(@(x) discretize(sq(x(:,'hcom',[1,2])), xyBin), dc2.xyz );
xyInd      = cf(@(x) discretize(sq(x(:,'hcom',[1,2])), xyBin), dc2.xyz );
hpInd      = cf(@(h) discretize(h(:,1), hpBin), dc4.fet);
bpInd      = cf(@(b) discretize(b(:,2), bpBin), dc4.fet);

% XY - xy error
dc2.xyEInd = cf(@(e)  discretize(e, xyeBin),  dc2.error.xy);
dc4.xyEInd = cf(@(e)  discretize(e, xyeBin),  dc4.error.xy);
% XYHB - hp and bp error
dc4.hpEInd = cf(@(e)  discretize(e, hpeBin),  dc4.error.hp);
dc4.bpEInd = cf(@(e)  discretize(e, bpeBin),  dc4.error.bp);
dc4.hbEInd = cf(@(e)  discretize(e, hbeBin),  dc4.error.hb);

%%%>>>

%%%<<< JPDF errors xy hb

ind = cf(@(ei1,ei2,hpi,bpi,hbi,u,s)                              ...
         s(:,1)==1 & any(s(:,3:6),2)                           ...
           & nniz(ei1) & nniz(ei2)                             ...
           & nniz(hpi) & nniz(bpi)                             ...
           & nniz(hbi) & nniz(u),                              ...
         dc2.xyEInd, dc4.xyEInd,                               ...
           dc4.hpEInd, dc4.bpEInd,                             ...
           dc4.hbEInd, uincInd, dc2.stcm);

JPDF_d2xyErr_Uinc =                                             ...
    cf(@(x,u,i)                                                 ...
       accumarray([x(i),u(i)],                                  ...
                  1,                                            ...
                  [numel(xyeBin)-1,numel(uincBin)-1],           ...
                  @sum,nan),                                        ...
       dc2.xyEInd, uincInd, ind);

JPDF_d4xyErr_Uinc =                                             ...
    cf(@(x,u,i)                                                 ...
       accumarray([x(i),u(i)],                                  ...
                  1,                                            ...
                  [numel(xyeBin)-1,numel(uincBin)-1],           ...
                  @sum,nan),                                        ...
       dc4.xyEInd, uincInd, ind);


JPDF_d4hpErr_Uinc =                                             ...
    cf(@(hp,u,i)                                                ...
       accumarray([hp(i),u(i)],                                 ...
               1,                                               ...
               [numel(hpeBin)-1,numel(uincBin)-1],              ...
               @sum,nan),                                           ...
       dc4.bpEInd, uincInd, ind);


JPDF_d4bpErr_Uinc =                                             ...
    cf(@(bp,u,i)                                                ...
       accumarray([bp(i),u(i)],                                 ...
               1,                                               ...
               [numel(bpeBin)-1,numel(uincBin)-1],              ...
               @sum,nan),...
       dc4.hpEInd, uincInd, ind);

JPDF_d4hbErr_Uinc =                                             ...
    cf(@(hb,u,i)                                                ...
       accumarray([hb(i),u(i)],                                 ...
               1,                                               ...
               [numel(hbeBin)-1,numel(uincBin)-1],              ...
               @sum,nan),...
       dc4.hbEInd, uincInd, ind);

% $$$ CONEXP_d4hbErr_Uinc_post =                                      ...
% $$$     cf(@(hb,u,i)                                                ...
% $$$        accumarray([hb(i),u(i)],                                 ...
% $$$                dc4.post,                                        ...
% $$$                [numel(bpeBin)-1,numel(uincBin)-1],              ...
% $$$                @sum,nan),...
% $$$        dc4.hbEInd, uincInd, nind,ind);



xyeBinCenters  = (xyeBin(1:end-1)+xyeBin(2:end))./2;
uincBinCenters = (uincBin(1:end-1)+uincBin(2:end))./2;
hpeBinCenters  = (hpeBin(1:end-1)+hpeBin(2:end))./2;
bpeBinCenters  = (bpeBin(1:end-1)+bpeBin(2:end))./2;
hbeBinCenters  = (hbeBin(1:end-1)+hbeBin(2:end))./2;

y = 1;
figure();
colormap('jet');
    subplot2(5,2,y,1);
        imagesc(xyeBinCenters./10,                              ...
                uincBinCenters,                                 ...
                sum(cat(3, JPDF_d2xyErr_Uinc{:}),3,'omitnan')');
        axis('xy');
        y = y+1;
    subplot2(5,2,y,1);
        imagesc(xyeBinCenters./10,                              ...
                uincBinCenters,                                 ...
                sum(cat(3, JPDF_d4xyErr_Uinc{:}),3,'omitnan')');
        axis('xy');
        y = y+1;                
    subplot2(5,2,y,1);
        imagesc(hbeBinCenters,                                  ...
                uincBinCenters,                                 ...
                sum(cat(3, JPDF_d4hbErr_Uinc{:}),3,'omitnan')');
        axis('xy');
        y = y+1;        
    subplot2(5,2,y,1);
        imagesc(hpeBinCenters,                                  ...
                uincBinCenters,                                 ...
                sum(cat(3, JPDF_d4hpErr_Uinc{:}),3,'omitnan')');
        axis('xy');
        y = y+1;
    subplot2(5,2,y,1);
        imagesc(bpeBinCenters,                                  ...
                uincBinCenters,                                 ...
                sum(cat(3, JPDF_d4bpErr_Uinc{:}),3,'omitnan')');
        axis('xy');

% $$$ t = 9;
% $$$ figure,subplot(211);plot(dc4.fet{t}(:,1)),hold on,plot(dc4.com{t}(:,3)); 
% $$$       subplot(212);plot(dc4.fet{t}(:,2)),hold on,plot(dc4.com{t}(:,4));
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');

%%%>>>
        
        
        
CONEXP_xy_hbErr =                                               ...
    cf(@(xy,bp,u,i)                                             ...
       accumarray([xy(i,1),xy(i,2)],                            ...
                  hbeBin(bp(i)),                                ...
                  [numel(xyBin)-1,numel(xyBin)-1],              ...
                  @mean,                                        ...
                  NaN),                                         ...
       xyInd, dc4.hpEInd, uincInd, ind);

CONSTD_xy_hbErr =                                               ...
    cf(@(xy,bp,u,i)                                             ...
       accumarray([xy(i,1),xy(i,2)],                            ...
                  hbeBin(bp(i)),                                ...
                  [numel(xyBin)-1,numel(xyBin)-1],              ...
                   @std,                                        ...
                  NaN),                                         ...
       xyInd, dc4.hpEInd, uincInd, ind);
             
figure();
for s = 1:numTrials,
    subplot2(2,numTrials,1,s);
    imagescnan({linspace(-500,500,numel(xyBin)-1),              ...
                linspace(-500,500,numel(xyBin)-1),              ...
                CONEXP_xy_hbErr{s}'},                           ...
               [-0.75,0.75],                                    ...
               'linear',                                        ...
               true,                                            ...
               [0.2,0.2,0.2],                                   ...
               [],[],@jet);
    
    subplot2(2,numTrials,2,s);
    imagescnan({linspace(-500,500,numel(xyBin)-1),              ...
                linspace(-500,500,numel(xyBin)-1),              ...
                CONSTD_xy_hbErr{s}'},                           ...
               [0,0.5],                                         ...
               'linear',                                        ...
               true,                                            ...
               [0.2,0.2,0.2],                                   ...
               [],[],@jet);
end




nind = cf(@(ei1,ei2,hpi,bpi,hbi,hp,bp,u,s)                              ...
         s(:,1)==1 & any(s(:,2:6),2)                           ...
           & nniz(ei1) & nniz(ei2)                             ...
           & nniz(hpi) & nniz(bpi)                             ...
           & nniz(hbi) & nniz(u)...
           & nniz(hp)  & nniz(bp),                              ...
         dc2.xyEInd, dc4.xyEInd,                               ...
           dc4.hpEInd, dc4.bpEInd,                             ...
           dc4.hbEInd, hpInd, bpInd, ...
           uincInd, dc2.stcm);

CONEXP_hb_xyErr =                                               ...
    cf(@(hp,bp,xye,u,i)                                         ...
       accumarray([hp(i), bp(i)],                               ...
                  xyeBin(xye(i)),                               ...
                  [numel(hpBin)-1,numel(bpBin)-1],              ...
                  @mean,                                        ...
                  NaN),                                         ...
       hpInd, bpInd, dc4.xyEInd, uincInd, nind);
CONSTD_hb_xyErr =                                               ...
    cf(@(hp,bp,xye,u,i)                                         ...
       accumarray([hp(i), bp(i)],                               ...
                  xyeBin(xye(i)),                               ...
                  [numel(hpBin)-1,numel(bpBin)-1],              ...
                  @std,                                         ...
                  NaN),                                         ...
       hpInd, bpInd, dc4.xyEInd, uincInd, nind);

figure();
for s = 1:numTrials,
    subplot2(2,numTrials,1,s);
    imagescnan({linspace([hpBin([1,end]),numel(hpBin)-1]),      ...
                linspace([bpBin([1,end]),numel(bpBin)-1]),      ...
                CONEXP_hb_xyErr{s}'},                           ...
               [0,200],                                         ...
               'linear',                                        ...
               true,                                            ...
               [0.2,0.2,0.2],                                   ...
               [],[],@jet);
    axis('xy');
    
    subplot2(2,numTrials,2,s);
    imagescnan({linspace([hpBin([1,end]),numel(hpBin)-1]),      ...
                linspace([bpBin([1,end]),numel(bpBin)-1]),      ...
                CONSTD_hb_xyErr{s}'},                           ...
               [0,100],                                         ...
               'linear',                                        ...
               true,                                            ...
               [0.2,0.2,0.2],                                   ...
               [],[],@jet);
    axis('xy');
end









% DCTPD 
% LOAD decoded position based on high sampleRate and binned by theta phase
load(fullfile(MTA_PROJECT_PATH,'analysis','pfsXYHB_mask.mat'));
%dct.phzBins = -pi:pi/4:pi;
dct.phzBins = 0:pi/4:2*pi;
dct.states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};
dct.sampleRate = 250;
dct.spikeWindow = 0.3;
dct.com  = cell([1,numTrials]);
dct.max  = cell([1,numTrials]);
dct.sax  = cell([1,numTrials]);
dct.post = cell([1,numTrials]);
dct.uinc = cell([1,numTrials]);
dct.ubhv = cell([1,numTrials]);
dct.xyz  = cf(@(t)    resample(preproc_xyz(t,'trb'),dct.sampleRate),      Trials);
dct.stc  = cf(@(t)    t.load('stc','msnn_ppsvd_raux'),                    Trials);
dct.stcm = cf(@(s,x)  stc2mat(s,x,dct.states),                   dct.stc,dct.xyz);
dct.pfs  = cf(@(t,u)  compute_xyhb_ratemaps(t,u),                  Trials, units);

for t = 1:numTrials
    
    phz = load_theta_phase(Trials{t},                           ...
                           Trials{t}.lfp.sampleRate,            ...
                           sessionList(tind(t)).thetaRefGeneral,...
                           0);
    %phzCorrection(tind(t)));
    phzBinInds = discretize(phz.data,dct.phzBins);
    
    spk  = create(copy(Trials{t}.spk),Trials{t},lfp.sampleRate,'',units{t},'deburst');
    
    for p = 1:numel(dct.phzBins)-1
        % LOAD phase specific spikes
        spkt = copy(spk);
        spkt.clu = spkt.clu(phzBinInds(spkt.res)==p);
        spkt.res = round(spkt.res(phzBinInds(spkt.res)==p)./Trials{t}.lfp.sampleRate.*dct.sampleRate)+1;
        spkt.sampleRate = dct.sampleRate;
        
        % GENERATE unit firing rate object
        dct.ufr{t}{p} = Trials{t}.load('ufr',             ... load target
                                       dct.xyz{t},        ... reference object
                                       spkt,              ... spike container
                                       units{t},          ... units
                                       dct.spikeWindow,   ... spike windoww
                                       'gauss',           ... mode
                                       true);               % overwrite
        dct.ufr{t}{p}.data =  dct.ufr{t}{p}.data.*6;
        dct.uinc{t}(:,p) = sum(dct.ufr{t}{p}.data>0.2,2);
        
        tag = ['xyhb_thpD',num2str(dct.sampleRate),'_',num2str(p),'o',num2str(numel(phzBins)-1)];
        %tag = ['xyhb_phzcrt_thpD',num2str(dct.sampleRate),'_',num2str(p),'o',num2str(numel(phzBins)-1)];        
        % DECODE xyhb
        [dct.com{t}(:,:,p),dct.max{t}(:,:,p),dct.sax{t}(:,:,p),dct.post{t}(:,p)] = ...
            decode_ufr(Trials{t},units{t},dct.sampleRate,dct.ufr{t}{p},dct.pfs{t},  ...
                       [],mask,smoothingWeights,'tag',tag,'overwrite',false);
    end
end


for t = 1:numTrials
    dct.ubhv{t} = zeros([size(dct.ufr{t}{1}),3]);
    for u = 1:numel(units{t})
        usind = find(ismember(cluSessionMapSubset,[tind(t),units{t}(u)],'rows'));
        if ~isempty(usind) 
           dct.ubhv{t}(:,u,:) = permute(repmat(fsrcz(usind,1:3),size(dct.ubhv{t},1),1),[1,3,2]);
        end
    end
end



% only because recomputing with corrected phases would take forever
dct.phzCorrection = [2,2,2,1,1,1,1,1,1,1];
for t = 1:numTrials,
    dct.ufr{t} = circshift(dct.ufr{t},dct.phzCorrection(t));
    dct.com{t} = circshift(dct.com{t},dct.phzCorrection(t),3);
    dct.max{t} = circshift(dct.max{t},dct.phzCorrection(t),3);
    dct.sax{t} = circshift(dct.sax{t},dct.phzCorrection(t),3);
    dct.post{t} = circshift(dct.post{t},dct.phzCorrection(t),2);
    dct.uinc{t} = circshift(dct.uinc{t},dct.phzCorrection(t),2);
end


dct.ufrp = cell([1,numTrials]);
for t = 1:numTrials,
    for p = 1:numel(dct.phzBins)-1
        if p == 1,
            dct.ufrp{t} = copy(dct.ufr{t}{p});
        else
            dct.ufrp{t}.data = cat(3,dct.ufrp{t}.data, dct.ufr{t}{p}.data);
        end
    end
end



dct.hRot = {0,0,0,0.17,0.17,0.17,0.17,0.17,0.17,0.17};
dct.hvec = cf(@(x)    x(:,'nose',[1,2])-x(:,'hcom',[1,2]),                           dct.xyz      );
dct.hvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                      dct.hvec     );
dct.hvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                     dct.hvec     );
dct.hvec = cf(@(h,r)  multiprod(h,[cos(r),-sin(r);sin(r),cos(r)],[2,3],[1,2]),       dct.hvec, dct.hRot);
dct.tvec = cf(@(x)    circshift(x(:,'hcom',[1,2]),-round(250*0.1))-circshift(x(:,'hcom',[1,2]),round(250*0.1)), dct.xyz   );
dct.tvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                      dct.tvec     );
dct.tvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                     dct.tvec     );
dct.ts   = cf(@(x)    [1:size(x,1)]./x.sampleRate,                                   dct.xyz      );
dct.fet  = cf(@(t)    fet_HB_pitchB(t,dct.sampleRate),                               Trials       );

% COMPILE egocentric variables 

% $$$ dct.hvec = cf(@(x)    x(:,'head_front',[1,2])-x(:,'head_back',[1,2]),                    dct.xyz      );
% $$$ dct.hvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                          dct.hvec     );
% $$$ dct.hvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                         dct.hvec     );
% $$$ dct.tvec = cf(@(x)    circshift(x(:,'hcom',[1,2]),-1)-circshift(x(:,'hcom',[1,2]),1),    dct.xyz      );
% $$$ dct.tvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                          dct.tvec     );
% $$$ dct.tvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                         dct.tvec     );

% COMPUTE error

% project 
dct.ego.com = cf(@(e,x,h,f)  cat(2,...
                          sq(multiprod(permute(bsxfun(@minus,...
                                                  e(:,[1,2],:),...
                                                  sq(x(:,'hcom',[1,2]))),...
                                           [1,2,4,3]),...
                                   h(:,:,:),2,[2,3])),...
                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
          dct.com,dct.xyz,dct.hvec,dct.fet);
dct.ego.sax = cf(@(e,x,h,f)  cat(2,...
                          sq(multiprod(permute(bsxfun(@minus,...
                                                  e(:,[1,2],:),...
                                                  sq(x(:,'hcom',[1,2]))),...
                                           [1,2,4,3]),...
                                   h(:,:,:),2,[2,3])),...
                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
          dct.sax,dct.xyz,dct.hvec,dct.fet);
dct.ego.max = cf(@(e,x,h,f)  cat(2,...
                          sq(multiprod(permute(bsxfun(@minus,...
                                                  e(:,[1,2],:),...
                                                  sq(x(:,'hcom',[1,2]))),...
                                           [1,2,4,3]),...
                                   h(:,:,:),2,[2,3])),...
                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
          dct.max,dct.xyz,dct.hvec,dct.fet);




%%%<<< COMPUTE error for TDP decoding per state

stid = [1,2,3,4,5,6];
dcTDPfrontal = repmat({cell([1,numTrials])},[1,numel(stid)]);
dcTDPlateral = repmat({cell([1,numTrials])},[1,numel(stid)]);
indTDP = cell([1,numel(stid)]);
errorBinEdges = linspace(-500,500,100);      
errorBinCenters = mean(GetSegs(errorBinEdges,1:numel(errorBinEdges)-1,2));
for sts = 1:numel(stid),
    indTDP{sts}  = cf(@(p,u,s,b,r)                                       ...
                      all(p>0.0005,2)                                    ...
                      & sum(double(u>=1),2)>6                            ...
                      & sum(double(u>=1),2)>6                            ...                      
                      & s(:,1)==1                                        ...
                      & any(logical(s(:,stid(sts))),2)                   ...
                      & ~any(logical(s(:,[7,8])),2),                     ...
                      dct.post,dct.uinc,dct.stcm,dct.ubhv,dct.ufr);
    for t = 1:numel(Trials),
        for j = 1:numel(phzBins)-1;
            dcTDPfrontal{sts}{t} =                                             ...
                cat(2,                                                          ...
                    dcTDPfrontal{sts}{t},                                      ...
                    histcounts(sq(dErr{t}(indTDP{sts}{t},1,j)),errorBinEdges)');
            dcTDPlateral{sts}{t} =                                             ...
                cat(2,                                                          ...
                    dcTDPlateral{sts}{t},                                      ...
                    histcounts(sq(dErr{t}(indTDP{sts}{t},2,j)),errorBinEdges)');
        end
    end
end

% $$$ testInd = all(dct.post{7}>0.0005,2)                                    ...
% $$$           & sum(double(dct.uinc{7}>=1),2)>6                            ...
% $$$     & dct.stcm{7}(:,1)==1                                        ...
% $$$     & any(logical(dct.stcm{7}(:,stid(sts))),2)                   ...
% $$$     & ~any(logical(dct.stcm{7}(:,[7,8])),2);



%%%>>>


%%%<<< COMPUTE error for TDP decoding per state
dErr = dct.ego.sax;

pbound =[-1.35,1.35];
stid = [1,2,3,4,5,6];
dcTDPfrontal = repmat({cell([1,numTrials])},[1,numel(stid)]);
dcTDPlateral = repmat({cell([1,numTrials])},[1,numel(stid)]);
dcTDPhp = repmat({cell([1,numTrials])},[1,numel(stid)]);
dcTDPbp = repmat({cell([1,numTrials])},[1,numel(stid)]);
indTDP = cell([1,numel(stid)]);
errorBinEdges = linspace(-500,500,100);      
errorBinCenters = mean(GetSegs(errorBinEdges,1:numel(errorBinEdges)-1,2));
perrorBinEdges = linspace([pbound,50]);      
perrorBinCenters = mean(GetSegs(perrorBinEdges,1:numel(perrorBinEdges)-1,2));

for sts = 1:numel(stid),
    indTDP{sts}  = cf(@(p,u,s)                                           ...
                      all(p>0.0005,2)                                    ...
                      & sum(double(u>=3),2)>6                            ...
                      & s(:,1)==1                                        ...
                      & any(logical(s(:,stid(sts))),2)                   ...
                      & ~any(logical(s(:,[7,8])),2),                     ...
                      dct.post,dct.uinc,dct.stcm);
    for t = 1:numel(Trials),
        for j = 1:numel(phzBins)-1;
            dcTDPfrontal{sts}{t} =                                               ...
                cat(2,                                                           ...
                    dcTDPfrontal{sts}{t},                                        ...
                    histcounts(sq(dErr{t}(indTDP{sts}{t},1,j)),errorBinEdges)');
            dcTDPlateral{sts}{t} =                                               ...
                cat(2,                                                           ...
                    dcTDPlateral{sts}{t},                                        ...
                    histcounts(sq(dErr{t}(indTDP{sts}{t},2,j)),errorBinEdges)');
            dcTDPhp{sts}{t} =                                               ...
                cat(2,                                                           ...
                    dcTDPhp{sts}{t},                                        ...
                    histcounts(sq(dErr{t}(indTDP{sts}{t},3,j)),perrorBinEdges)');
            dcTDPbp{sts}{t} =                                               ...
                cat(2,                                                           ...
                    dcTDPbp{sts}{t},                                        ...
                    histcounts(sq(dErr{t}(indTDP{sts}{t},4,j)),perrorBinEdges)');
        end
    end
end

%%%>>>



%%%<<< COMPUTATION : find ridge values for each state's dcTDPfrontal error as function of theta phz

dcTpdLongSum = {};
dcTpdLatSum = {};
dcTpdhpSum = {};
dcTpdbpSum = {};
dcTpdInd = [];
dcTpdWmn = [];
dcTpdMax = [];
for sts = 1:numel(stid),
    dcTpdLongSum{sts} = RectFilter(sum(cat(3,dcTDPfrontal{sts}{:}),3),5,3);
    dcTpdLongSum{sts} = bsxfun(@rdivide,dcTpdLongSum{sts},sum(dcTpdLongSum{sts}));
    [~, dcTpdInd(end+1,:)] = max(dcTpdLongSum{sts});

    dcTpdLatSum{sts} = RectFilter(sum(cat(3,dcTDPlateral{sts}{:}),3),5,3);
    dcTpdLatSum{sts} = bsxfun(@rdivide,dcTpdLatSum{sts},sum(dcTpdLatSum{sts}));
    
    dcTpdMax(end+1,:) = errorBinCenters(dcTpdInd(end,:)');    
    dcTpdWmn(end+1,:) = sum(bsxfun(@times,dcTpdLongSum{sts},errorBinCenters'))';
    
    dcTpdhpSum{sts} = RectFilter(sum(cat(3,dcTDPhp{sts}{:}),3),5,3);
    dcTpdhpSum{sts} = bsxfun(@rdivide,dcTpdhpSum{sts},sum(dcTpdhpSum{sts}));
    
    dcTpdbpSum{sts} = RectFilter(sum(cat(3,dcTDPbp{sts}{:}),3),5,3);
    dcTpdbpSum{sts} = bsxfun(@rdivide,dcTpdbpSum{sts},sum(dcTpdhpSum{sts}));
end




%%%>>>



%%%<<< LOAD example data
t = 7;
Trial = Trials{t};
unitSubset = units{t};
exyz = resample(Trial.load('xyz','trb'),30);
efet = fet_HB_pitchB(Trial,sampleRate);
resample(efet,sampleRate);
ets = [1:size(exyz,1)]./sampleRate;
estcm = stc2mat(Trial.load('stc'),exyz,states);
eind =   dc4.post{t} > 0.002 ...
       & dc2.uinc{t} > 3     ...
       & estcm(:,1)==1;
enanmask = ones([size(exyz,1),1]);
enanmask(~eind)=nan;


%%%>>>

%%%<<< MAIN FIGURE portrait ------------------------------------------------------------------------


eTS = [290,345];
%eTS = [530,600];
%eTS = [1000,1045];
eInds = round(eTS.*sampleRate)+1;
eInds = [eInds(1):eInds(2)]';
eIndsSub = eInds(1300:1600);

%startfig

t = 7;
[hfig,fig,fax,sax] = set_figure_layout(figure(666006),'A4','portrait',[],2,2,0.2,0.2);

errType = 'com';

%%%<<< PLOT example trajectory within physical space
% ADJUST subplot coordinates
t = 7; 
[yind, yOffSet, xind, xOffSet] = deal(1,-fig.subplot.height, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*2,                        ...
                              fig.subplot.height*2],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold(sax(end),'on');
imagesc(dc4.pfs{1}.adata.bins{1:2},~maskcirc');
colormap(sax(end),'gray')
caxis(sax(end),[-2,1])

circle(0,0,480,'-k');
plot(exyz(eIndsSub,5,1),...
     exyz(eIndsSub,5,2),...
     '-k','LineWidth',1);
plot(sq(dc4.(errType){t}(eIndsSub,1)).*enanmask(eIndsSub),...
     sq(dc4.(errType){t}(eIndsSub,2)).*enanmask(eIndsSub),...
     '-c','LineWidth',1);
xlim([-500,500]);
ylim([-500,500]);
sax(end).XTick = [];
sax(end).YTick = [];
%%%>>>




%%%<<< PLOT example trajectory within behavior space
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
hold(sax(end),'on');
imagesc(dc4.pfs{1}.adata.bins{3:4},~maskbhv');
colormap(sax(end),'gray')
caxis(sax(end),[-2,1])

axis('xy');
plot(efet(eIndsSub,1),...
     efet(eIndsSub,2),...
     '-k','LineWidth',1);
plot(sq(dc4.(errType){t}(eIndsSub,3)).*enanmask(eIndsSub),...
     sq(dc4.(errType){t}(eIndsSub,4)).*enanmask(eIndsSub),...
     '-c','LineWidth',1);
xlim(sax(end),[-1.8,0.6]);
ylim(sax(end),[-0.6,1.8]);
sax(end).XTick = [];
sax(end).YTick = [];
%%%>>>




%%%<<< PLOT timeseries of xcoord vs decoded xcoord
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
     sq(dc4.(errType){t}(eInds,1)).*enanmask(eInds),...
     '-r','LineWidth',1);
plot(ets(eIndsSub),...
     sq(dc4.(errType){t}(eIndsSub,1)).*enanmask(eIndsSub),...
     '-c','LineWidth',1);
sax(end).XTick = [];
sax(end).YTick = [];
%%%>>>




%%%<<< PLOT timeseries of ycoord vs decoded ycoord
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
     sq(dc4.(errType){t}(eInds,2)).*enanmask(eInds),...
     '-r','LineWidth',1);
plot(ets(eIndsSub),...
     sq(dc4.(errType){t}(eIndsSub,2)).*enanmask(eIndsSub),...
     '-c','LineWidth',1);
sax(end).XTick = [];
sax(end).YTick = [];
%%%>>>




%%%<<< PLOT timeseries of head pitch vs decoded head pitch
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
     sq(dc4.(errType){t}(eInds,3)).*enanmask(eInds),...
     '-r','LineWidth',1);
plot(ets(eIndsSub),...
     sq(dc4.(errType){t}(eIndsSub,3)).*enanmask(eIndsSub),...
     '-c','LineWidth',1);
sax(end).XTick = [];
sax(end).YTick = [];
ylim([-1.6,0.5]);
%%%>>>




%%%<<< PLOT timeseries of body pitch vs decoded body pitch
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
     sq(dc4.(errType){t}(eInds,4)).*enanmask(eInds),...
     '-r','LineWidth',1);
plot(ets(eIndsSub),...
     sq(dc4.(errType){t}(eIndsSub,4)).*enanmask(eIndsSub),...
     '-c','LineWidth',1);

sax(end).YTick = [];
ylim([-0.5,1.6]);
%%%>>>




%%%<<< PLOT longitudinal error for TDP decoding
% SET constants
xcoords =linspace(-500,500,250);
ycoords = circ_rad2ang([phzBinCenters;phzBinCenters+2*pi]);
% ADJUST subplot coordinates
xOffSet = -0.2;
for sts = 1:numel(stid),
% ADJUST subplot coordinates    
    [yind, yOffSet, xind, xOffSet] = deal(5, -0.8, sts, xOffSet+0.2);
% CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
% PLOT JPDF of longitudinal error vs theta phase
    imagesc(xcoords,                                 ...
            ycoords,       ...
            repmat(RectFilter(dcTpdLongSum{sts},5,3),1,2)')
    axis('xy');
    xl = [-200,200];
    xlim(xl);
    ylim([-60,420]);            
% MASK areas outside single cycle
    patch([xl(1),xl(1),xl(2),xl(2)],...
          [ 360, 420,420,360],...
          [0.2,0.2,0.2],...
          'EdgeAlpha',0,...
          'FaceAlpha',0.4);
    patch([xl(1),xl(1),xl(2),xl(2)],...    
          [ -60,   0,  0,-60],...
          [0.2,0.2,0.2],...
          'EdgeAlpha',0,...
          'FaceAlpha',0.4);
    sax(end).XTick = [-100,0,100];
    sax(end).XTickLabels = [-15,0,15];
    sax(end).YTick = [];    
    sax(end).XTickLabels = [];    
    Lines(0,[],'w',[],'LineWidth',1);
    plot(repmat(dcTpdMax(sts,:),[1,2])', ycoords,'k','LineWidth',2);
    %plot(repmat(dcTpdMax(sts,:),[1,2])', ycoords,'k','LineWidth',2);
    %plot(repmat(dcTpdWmn(1,:),[1,2])', ycoords,'k','LineWidth',2);    
    %plot(repmat(dcTpdWmn(sts,:),[1,2])', ycoords,'m','LineWidth',2);
    title(states{sts});
    colormap(sax(end),'jet');    
    
    if sts == 1,
        ylabel(sax(end),{'RCP Error'});
    end
    
end
% ADJUST subplot coordinates        
[yind, yOffSet, xind, xOffSet] = deal(5, -0.8, sts+1, xOffSet);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,      ...
                              fig.page.ypos(yind)+yOffSet,      ...
                              fig.subplot.width/1.5,              ...
                              fig.subplot.height],              ...
                  'FontSize', 8,                                ...
                  'LineWidth',1);
hold(sax(end),'on');
% ADD phase "bar"
plot(-cos(circ_ang2rad(-60:420)),-60:420,'-k','LineWidth',2);

sax(end).YAxisLocation = 'right';
sax(end).Visible = 'off';                                    
sax(end).YTick = [];
sax(end).YTickLabels = [];
sax(end).XTick = [];
ylim([-60,420]);
xlim([-1.5,1.5]);


text( 0.25, 360, '360', 'HorizontalAlignment','center',     ...
         'FontSize', 8, 'VerticalAlignment',  'middle');
text(-0.25, 180, '180', 'HorizontalAlignment','center',     ...
         'FontSize', 8, 'VerticalAlignment',  'middle');                    
text( 0.25,   0,   '0', 'HorizontalAlignment','center',     ...
         'FontSize', 8, 'VerticalAlignment',  'middle');                    

%%%>>>




%%%<<< PLOT lateral error for TDP decoding
% SET constants
xcoords =linspace(-500,500,250);
ycoords = circ_rad2ang([phzBinCenters;phzBinCenters+2*pi]);
% ADJUST subplot coordinates
xOffSet = -0.2;
for sts = 1:numel(stid),
% ADJUST subplot coordinates    
    [yind, yOffSet, xind, xOffSet] = deal(6, -0.8, sts, xOffSet+0.2);
% CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
% PLOT JPDF of lateral error vs theta phase
    imagesc(xcoords,                                 ...
            ycoords,       ...
            repmat(RectFilter(dcTpdLatSum{sts},5,3),1,2)')
    axis('xy');
    xl = [-200,200];
    xlim(xl);
    ylim([-60,420]);            
% MASK areas outside single cycle
    patch([xl(1),xl(1),xl(2),xl(2)],...
          [ 360, 420,420,360],...
          [0.2,0.2,0.2],...
          'EdgeAlpha',0,...
          'FaceAlpha',0.4);
    patch([xl(1),xl(1),xl(2),xl(2)],...    
          [ -60,   0,  0,-60],...
          [0.2,0.2,0.2],...
          'EdgeAlpha',0,...
          'FaceAlpha',0.4);
    sax(end).XTick = [-100,0,100];
    sax(end).XTickLabels = [-15,0,15];
    sax(end).YTick = [];    
    sax(end).XTickLabels = [];    
    Lines(0,[],'w',[],'LineWidth',1);
    
    %title(states{sts});
    colormap(sax(end),'jet');    
    
    if sts == 1,
        ylabel(sax(end),{'Lat Error'});
    end
    
end
% ADJUST subplot coordinates        
[yind, yOffSet, xind, xOffSet] = deal(6, -0.8, sts+1, xOffSet);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,      ...
                              fig.page.ypos(yind)+yOffSet,      ...
                              fig.subplot.width/3,              ...
                              fig.subplot.height],              ...
                  'FontSize', 8,                                ...
                  'LineWidth',1);
hold(sax(end),'on');
% ADD phase "bar"
plot(-cos(circ_ang2rad(-60:420)),-60:420,'-k','LineWidth',2);


sax(end).YAxisLocation = 'right';
sax(end).Visible = 'off';                                    
sax(end).YTick = [];
sax(end).YTickLabels = [];
sax(end).XTick = [];
ylim([-60,420]);
xlim([-1.5,1.5]);

%%%>>>




%%%<<< PLOT head pitch error
% SET constants
xcoords =linspace([-1.25,1.25,50]);
ycoords = circ_rad2ang([phzBinCenters;phzBinCenters+2*pi]);
% ADJUST subplot coordinates
xOffSet = -0.2;
for sts = 1:numel(stid),
% ADJUST subplot coordinates    
    [yind, yOffSet, xind, xOffSet] = deal(7, -0.8, sts, xOffSet+0.2);
% CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
% PLOT JPDF of longitudinal error vs theta phase
    imagesc(xcoords,                                 ...
            ycoords,       ...
            repmat(RectFilter(dcTpdhpSum{sts},5,3),1,2)')
    axis('xy');
    xl = [pbound];
    xlim(xl);
    ylim([-60,420]);            
% MASK areas outside single cycle
    patch([xl(1),xl(1),xl(2),xl(2)],...
          [ 360, 420,420,360],...
          [0.2,0.2,0.2],...
          'EdgeAlpha',0,...
          'FaceAlpha',0.4);
    patch([xl(1),xl(1),xl(2),xl(2)],...    
          [ -60,   0,  0,-60],...
          [0.2,0.2,0.2],...
          'EdgeAlpha',0,...
          'FaceAlpha',0.4);
    sax(end).XTick = [];
    sax(end).XTickLabels = [];
    sax(end).YTick = [];    
    sax(end).XTickLabels = [];    
    Lines(0,[],'w',[],'LineWidth',1);
    %plot(repmat(dcTpdMax(sts,:),[1,2])', ycoords,'k','LineWidth',2);
    %plot(repmat(dcTpdMax(sts,:),[1,2])', ycoords,'k','LineWidth',2);
    %plot(repmat(dcTpdWmn(1,:),[1,2])', ycoords,'k','LineWidth',2);    
    %plot(repmat(dcTpdWmn(sts,:),[1,2])', ycoords,'m','LineWidth',2);
    %title(states{sts});
    colormap(sax(end),'jet');    
    if sts == 1,
        %ylabel(sax(end),{'Head','Pitch'});
        ylabel(sax(end),{'HP Error'});        
    end
    
end
% ADJUST subplot coordinates        
[yind, yOffSet, xind, xOffSet] = deal(7, -0.8, sts+1, xOffSet);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,      ...
                              fig.page.ypos(yind)+yOffSet,      ...
                              fig.subplot.width/3,              ...
                              fig.subplot.height],              ...
                  'FontSize', 8,                                ...
                  'LineWidth',1);
hold(sax(end),'on');
% ADD phase "bar"
plot(-cos(circ_ang2rad(-60:420)),-60:420,'-k','LineWidth',2);

sax(end).YAxisLocation = 'right';
sax(end).Visible = 'off';                                    
sax(end).YTick = [];
sax(end).YTickLabels = [];
sax(end).XTick = [];
ylim([-60,420]);
xlim([-1.25,1.25]);



%%%>>>



%%%<<< PLOT body pitch error
% SET constants

xcoords =linspace([pbound,50]);
ycoords = circ_rad2ang([phzBinCenters;phzBinCenters+2*pi]);
% ADJUST subplot coordinates
xOffSet = -0.2;
for sts = 1:numel(stid),
% ADJUST subplot coordinates    
    [yind, yOffSet, xind, xOffSet] = deal(8, -0.8, sts, xOffSet+0.2);
% CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
% PLOT JPDF of longitudinal error vs theta phase
    imagesc(xcoords,                                 ...
            ycoords,       ...
            repmat(RectFilter(dcTpdbpSum{sts},5,3),1,2)')
    axis('xy');
    xl = [pbound];
    xlim(xl);
    ylim([-60,420]);            
% MASK areas outside single cycle
    patch([xl(1),xl(1),xl(2),xl(2)],...
          [ 360, 420,420,360],...
          [0.2,0.2,0.2],...
          'EdgeAlpha',0,...
          'FaceAlpha',0.4);
    patch([xl(1),xl(1),xl(2),xl(2)],...    
          [ -60,   0,  0,-60],...
          [0.2,0.2,0.2],...
          'EdgeAlpha',0,...
          'FaceAlpha',0.4);
    sax(end).XTick = [];
    sax(end).XTickLabels = [];
    sax(end).YTick = [];    
    sax(end).XTickLabels = [];    
    Lines(0,[],'w',[],'LineWidth',1);
    %plot(repmat(dcTpdMax(sts,:),[1,2])', ycoords,'k','LineWidth',2);
    %plot(repmat(dcTpdMax(sts,:),[1,2])', ycoords,'k','LineWidth',2);
    %plot(repmat(dcTpdWmn(1,:),[1,2])', ycoords,'k','LineWidth',2);    
    %plot(repmat(dcTpdWmn(sts,:),[1,2])', ycoords,'m','LineWidth',2);
    %title(states{sts});
    colormap(sax(end),'jet');    
    if sts == 1,
        %ylabel(sax(end),{'Body','Pitch'});
        ylabel(sax(end),{'BP Error'});        
    end
end

% ADJUST subplot coordinates        
[yind, yOffSet, xind, xOffSet] = deal(8, -0.8, sts+1, xOffSet);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,      ...
                              fig.page.ypos(yind)+yOffSet,      ...
                              fig.subplot.width/3,              ...
                              fig.subplot.height],              ...
                  'FontSize', 8,                                ...
                  'LineWidth',1);
hold(sax(end),'on');
% ADD phase "bar"
plot(-cos(circ_ang2rad(-60:420)),-60:420,'-k','LineWidth',2);

sax(end).YAxisLocation = 'right';
sax(end).Visible = 'off';                                    
sax(end).YTick = [];
sax(end).YTickLabels = [];
sax(end).XTick = [];
ylim([-60,420]);
xlim(pbound);

%%%>>>

% stopfig


%>>>
% END MAIN FIGURE ----------------------------------------------------------------------------------



%%%<<< supfig
figure();
for t = 1:numTrials
    subplot2(numTrials,2,t,1);
    imagesc(linspace(-500,500,250),                             ...
            [phzBinCenters;phzBinCenters+2*pi],                 ...
            repmat(dcTDPfrontal{1}{t},1,2)');
    axis('xy');
    xlim([-300,300]);    
    
    subplot2(numTrials,2,t,2);    
    imagesc(linspace(-500,500,250),1:8,repmat(dcTDPlateral{1}{t},1,2)')
    axis('xy');
    xlim([-300,300]);
end
%%%>>>



ind  = cf(@(p,u,s) all(p>0.00005,2) & all(u>=2,2) & s(:,1)==1 & any(logical(s(:,stid)),2),...
          dct.post,dct.uinc,stcm);


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



dfet = dct.com{3};
%dfet = dct.sax{1};
decError = zeros([size(xyz{t},1),size(decEstSax{j},2),numel(phzBins)-1]);
for p = 1:numel(phzBins)-1,
    decError(:,[1,2],p) = [multiprod(dfet(:,[1,2],p)-sq(xyz{t}(:,'hcom',[1,2])),hvec{t}(:,:,:),2,[2,3])];
end            
decError(:,[3,4],:) = [bsxfun(@minus,fet{t}(:,1),dfet(:,3,:)),bsxfun(@minus,fet{t}(:,2),dfet(:,4,:))];


j = 3
ind =   all(dct.post{j}>0.0001,2)...
      & all(dct.uinc{j}>=4,2) ...
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
ind = all(dct.post{j}>0.000001,2)      & all(dct.uinc{t}>=4,2);
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
