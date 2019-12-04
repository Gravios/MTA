% MjgER2016_figure6
%
% DECODING position and behavior based on spatio-behavioral ratemaps
% 

% Describe the problem:
%     Hippocampal place cells' firing rate variance is greater than expected. (Fenton)
%     Previous studies lack the spatio temporal resolution required to separate
%        and examine the effects of sensori-motor state on firing rate variance. ( other studies)
%     
%     Extra firing rate variance is explained largely by behaivoral state (Figure 3)
%     
%     Does the 
%     
% 
% Describe
% Research 
% Function
% A
% Rip apart
% Model
%
%
% Decode{ window:300ms, space:xy }
% Decode{ window:300ms, space:xyhb }
% Decode{ window:300ms, space:xyhb, condition:thetaPhz[8]}
%
% A.1) example traj xy
%  .2) example traj hb
%
% B.1) example timeseries x
%  .2) example timeseries y
%  .3) example timeseries h
%  .4) example timeseries b
%
% C.0) JPDF longitudual head projection versus theta phaze
%  .1) theta
%  .2) rear
%  .3) hloc
%  .4) hpause
%  .5) lloc
%  .6) lpause
%
% D.0) 

%global AP
global MTA_PROJECT_PATH

%%%<<< LOAD data and SET analysis parameters

MjgER2016_load_data();
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};
spikeWindow = 0.3;
sampleRate = 30;
sampleRateTPD = 250;
phzBins = -pi:pi/4:pi;
phzBinCenters = sum([phzBins(1:end-1);phzBins(2:end)])'./2;
smoothingWeights = [250.^2,250.^2, 0.4.^2, 0.4.^2];
interpParPfs = struct('bins',{{linspace(-500,500,50),...
                               linspace(-500,500,50),...
                               linspace(  -2,  0.8,50),...
                               linspace(  -0.8,2,  50)}},...
                      'nanMaskThreshold', 0.1,...
                      'methodNanMap',     'linear',...
                      'methodRateMap',    'linear');
pbound =[-2.25,2.25];
perrorBinEdges = linspace([pbound,100]);      
perrorBinCenters = mean(GetSegs(perrorBinEdges,1:numel(perrorBinEdges)-1,2));

errorBinEdges = linspace(-500,500,100);      
errorBinCenters = mean(GetSegs(errorBinEdges,1:numel(errorBinEdges)-1,2));

uincBin = 0:2:40;
xyeBin = 0:10:400;
hpeBin = linspace([pbound,100]);
bpeBin = linspace([pbound,100]);
hbeBin = linspace(0,2.25,100);
xyBin = -500:50:500;
hpBin = -2:0.1:0.5;
bpBin = linspace(-0.5,2,50);



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

cluSessionMapSubset = cluSessionMap(unitSubsets,:);

% GET bhv ratemap erpPCA scores
% $$$ [fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS,fsrsMean,fsrsStd,rmapsShuffledMean,rmapsShuffled] = ...
% $$$     compute_bhv_ratemaps_erpPCA_scores(Trials,units,bfrm,bfrmShuff,eigVecs,validDims,unitSubsets,false);

% MAKE xyhb ratemap mask
% $$$ mask = repmat(maskcirc,[1,1,size(maskbhv)]).*repmat(permute(maskbhv,[3,4,1,2]),[size(maskcirc),1,1]);
%save(fullfile(MTA_PROJECT_PATH,'analysis','pfsXYHB_mask.mat'),'mask','-v7.3');
% LOAD xyhb ratemap mask
load(fullfile(MTA_PROJECT_PATH,'analysis','pfsXYHB_mask.mat'));


% SELECT trial subset
tind = [3:5,17:23];
Trials = Trials(tind);
units  = units(tind);
numTrials = numel(Trials);

stc  = cf(@(t)  t.load('stc','msnn_ppsvd_raux'),  Trials);

%%%>>>

%%%<<< DC2 DECODE xy 300ms window

dc2.tag = 'xy_sr30_sw300_HighRes';
dc2.sampleRate = 30;
dc2.smoothingWeights = [250.^2,250.^2];

%Trials{t}.lfp.filename = [Trials{t}.name,'.lfp'];    

dc2.phz = cf(@(t,c) ...
             phase( load( copy( t.lfp ), t, c), [6, 12]), ...
             Trials, num2cell([sessionList(tind).thetaRefGeneral]));
cf(@(p)   set(p,'data',unwrap(p.data)),         dc2.phz        );
cf(@(p,x) resample(p,x),                        dc2.phz,dc2.xyz);
cf(@(p)   set(p,'data',mod(p.data+pi,2*pi)-pi), dc2.phz        );
      


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
              t.load('ufr',x,s,u,spikeWindow,true,'gauss'),                     ...
              Trials, dc2.xyz, units, dc2.spk);



% GET placefield ratemap objects
dc2.pfs = cf(@(t,u)                                                             ...
              compute_ratemaps(t,u),                                            ...
              Trials, units);

% MAKE xy ratemap mask
maskcirc = create_tensor_mask(dc2.pfs{1}.adata.bins(1:2));
dc2.mask = maskcirc;

% GET unit inclusion 
dc2.uinc =  cf(@(u)                                                             ...
                sum(u.data>0.2,2),                                              ...
                dc2.ufr);

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
% SET time stamps
dc2.ts   = cf(@(x)                                                              ...
              [1:size(x,1)]./x.sampleRate,                                      ...
              dc2.xyz);

dc2.tRot = {0,0,0,0.17,0.17,0.17,0.17,0.17,0.17,0.17};
dc2.hvec = cf(@(x)    x(:,'head_front',[1,2])-x(:,'head_back',[1,2]),                dc2.xyz      );
dc2.hvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                      dc2.hvec     );
dc2.hvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                     dc2.hvec     );
dc2.hvec = cf(@(h,r)  multiprod(h,[cos(r),-sin(r);sin(r),cos(r)],[2,3],[1,2]),       dc2.hvec, dc2.tRot);
dc2.tvec = cf(@(x)    circshift(x(:,'hcom',[1,2]),-1)-circshift(x(:,'hcom',[1,2]),1),dc2.xyz      );
dc2.tvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                      dc2.tvec     );
dc2.tvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                     dc2.tvec     );

%%%>>>

%%%<<< DC4 DECODE xyhb 300ms window 

% MAKE bhv ratemap mask
maskbhv = false(bfrm{1}.adata.binSizes');
maskbhv(validDims) = true;

% GET xyhb placefields
dc4.pfs  = cf(@(t,u)  compute_xyhb_ratemaps(t,u),  Trials, units);

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

dc4.fet = cf(@(t)  fet_HB_pitchB(t,dc4.sampleRate),  Trials);

%%%>>>

%%%<<< COMPUTE decoding stats
dc2.error.xy   = cf(@(c,x)                                                      ...
                     sqrt( sum(( c(:,[1,2])-sq( x(:,'nose',[1,2]))).^2,2)),     ...
                     dc2.com, dc2.xyz);

dc2.error.ego = cf(@(e,x,h)                                                              ...
                   sq(multiprod(permute(bsxfun(@minus,                                   ...
                                               e(:,[1,2],:),                             ...
                                               sq(x(:,'hcom',[1,2]))),                   ...
                                        [1,2,4,3]),                                      ...
                                h(:,:,:),2,[2,3])),                                      ...
               dc2.com, dc2.xyz, dc2.hvec);



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


%%%<<< vectorize dc2 data

d2phz = cf(@(e) e(:,1), dc2.phz);
d2phz = cat(1,d2phz{:});
d2lon = cf(@(e) e(:,1), dc2.error.ego);
d2lon = cat(1,d2lon{:});
d2lat = cf(@(e) e(:,2), dc2.error.ego);
d2lat = cat(1,d2lat{:});

d2ufz = cf(@(e) e(:,1), dc2.ufz);
d2ufz = cat(1, d2ufz{:});

d2uinc = cat(1,dc2.uinc{:});
d2post = cat(1,dc2.post{:});
d2stcm = cat(1,dc2.stcm{:});



dc2.phzBinEdges = -pi:pi/8:pi;
dc2.phzBinCenters = mean(cat(1,dc2.phzBinEdges(2:end),dc2.phzBinEdges(1:end-1)));
d2phzInd = discretize(d2phz,dc2.phzBinCenters);
d2lonInd = discretize(d2lon,ferrorBinEdges{1});
d2latInd = discretize(d2lat,ferrorBinEdges{2});




d2ufz = cf(@(s,x) ...
         accumarray(s.res(s.res>0&s.res<size(x,1)),1,[size(x,1),1],@sum),...
         dc2.spk, dc2.xyz);
d2ufz = cat(1,d2ufz{:});

d2ufz = circ_mean(circshift(GetSegs(d2phz,1:size(d2phz,1),11,nan),-6,2),...
                circshift(GetSegs(d2ufz,1:size(d2ufz,1),11,nan),-6,2),...
                1)';


d2ufzInd = discretize(d2ufz,dc2.phzBinCenters);
%%%>>>


%%%<<< SUPFIG , egocentric projections of the decoded position as a function of 

%               windowed unit coactivation


uceglon = cf(@(e,u,s) ...
          hist2([e(any(logical(s(:,3:6)),2),1), u(any(logical(s(:,3:6)),2))],...
                errorBinEdges,1:40), ...
          dc2.error.ego, dc2.uinc, dc2.stcm);

uceglat = cf(@(e,u,s) ...
          hist2([e(any(logical(s(:,3:6)),2),2), u(any(logical(s(:,3:6)),2))],...
                linspace(-500,500,100),1:40), ...
             dc2.error.ego, dc2.uinc, dc2.stcm);

nuceglon = bsxfun(@rdivide,sum(cat(3,uceglon{:}),3),sum(sum(cat(3,uceglon{:}),3)));
nuceglat = bsxfun(@rdivide,sum(cat(3,uceglat{:}),3),sum(sum(cat(3,uceglat{:}),3)));

nuceglon = sum(cat(3,uceglon{:}),3);
nuceglat = sum(cat(3,uceglat{:}),3);

qtls = [0.25,0.50,0.75];
for u = 1:30
    [x, index] = unique(cumsum(nuceglon(:,u)));
    qntlNuLon(:,u) = interp1(x,ferrorBinCenters{1}(index),qtls);
    [x, index] = unique(cumsum(nuceglat(:,u)));
    qntlNuLat(:,u) = interp1(x,ferrorBinCenters{1}(index),qtls);
end

% STARTFIG
mins = [];
ubins = 1:3:25;
figure();
for u = 1:25,
    d2ind =   logical(d2stcm(:,1))                             ...
        & any(logical(d2stcm(:,[3,4,5,6])),2)                  ...
        & ~any(logical(d2stcm(:,[7,8])),2)               ...
        & d2post>=0.001                                     ...
        & ismember(d2uinc,[ubins(u):ubins(u)+2])  ...
        & nniz(d2lonInd) ...
        & nniz(d2latInd) ...
        & nniz(d2phzInd);


    subplot2(4,9,2,u);    
        out = accumarray([d2lonInd(d2ind),d2latInd(d2ind)],1,[numel(ferrorBinCenters{1}),numel(ferrorBinCenters{2})],@sum);
        [~,mins(u,1)] = max(reshape(imgaussfilt(out,[3,3]),[],1));
        [mins(u,1),mins(u,2)] = ind2sub([numel(ferrorBinCenters{1}),numel(ferrorBinCenters{2})],mins(u,1));
        imagesc(ferrorBinCenters{1},ferrorBinCenters{2},imgaussfilt(out,[3,3])');        
        axis('xy');
        Lines(0,[],'k');
        Lines([],0,'k');
        ylim([-250,250]);
        xlim([-250,250]);                

    subplot2(4,9,3,u);
        out = accumarray([d2lonInd(d2ind),d2latInd(d2ind)],...
                         d2ufz(d2ind),...
                         [numel(ferrorBinCenters{1}),numel(ferrorBinCenters{2})],...
                         @circ_mean);    
        out(out==0) = nan;
        imagescnan({ferrorBinCenters{1},ferrorBinCenters{2},out'},[-pi,pi],'circular',true,[0.8,0.8,0.8],'colorMap',@hsv);    
        Lines(0,[],'k');
        Lines([],0,'k');
        ylim([-250,250]);
        xlim([-250,250]);        
        
    subplot2(4,9,4,u);
        out = accumarray([d2lonInd(d2ind),d2latInd(d2ind)],...
                         d2ufz(d2ind),...
                         [numel(ferrorBinCenters{1}),numel(ferrorBinCenters{2})],...
                         @circ_var);
        imagescnan({ferrorBinCenters{1},ferrorBinCenters{2},out'},[0,1],'circular',true,[0.8,0.8,0.8],'colorMap',@jet);
        Lines(0,[],'k');
        Lines([],0,'k');
        ylim([-250,250]);
        xlim([-250,250]);                
    
    axis('xy');
end

    subplot2(4,9,1,1)
        d2ind =   logical(d2stcm(:,1))                             ...
        & any(logical(d2stcm(:,[3,4,5,6])),2)                  ...
        & ~any(logical(d2stcm(:,[7,8])),2)               ...
        & d2post>=0.001                                     ...
        & nniz(d2lonInd) ...
        & nniz(d2latInd) ...
        & nniz(d2phzInd);
        out = accumarray([d2lonInd(d2ind),d2phzInd(d2ind)],1,[numel(ferrorBinCenters{1}),numel(dc2.phzBinCenters)],@sum);
        imagesc(ferrorBinCenters{1},dc2.phzBinCenters,out');        
        axis('xy');
        Lines(0,[],'k');
        Lines([],0,'k');
        xlim([-250,250]);
        
    subplot2(4,9,1,2)
        out = accumarray([d2latInd(d2ind),d2phzInd(d2ind)],1,[numel(ferrorBinCenters{2}),numel(dc2.phzBinCenters)],@sum);
        imagesc(ferrorBinCenters{2},dc2.phzBinCenters,out');        
        axis('xy');
        Lines(0,[],'k');
        Lines([],0,'k');
        xlim([-250,250]);

    subplot2(4,9,1,5)
        plot(ferrorBinCenters{1}(mins(:,1)),1.5:3:26.5);
    subplot2(4,9,1,6)
        plot(ferrorBinCenters{2}(mins(:,2)),1.5:3:26.5);    

% ENDFIG
        
        
sax = gobjects([1,0]);
figure();
    sax(end+1) = subplot(211); hold('on');
        imagesc(linspace(-500,500,99),1:39,(nuceglon)');
        [mp,mpind] = max(nuceglon);
        %plot(errorBinCenters(mpind),1:39,'-k');
        plot(qntlNuLon(1,:),1:30,'-m');
        plot(qntlNuLon(2,:),1:30,'-g');        
        plot(qntlNuLon(3,:),1:30,'-m');
        caxis([0,0.055]);

    sax(end+1) = subplot(212); hold('on');
        imagesc(linspace(-500,500,99),1:39,(nuceglat)');
        [mp,mpind] = max(nuceglat);
        %plot(errorBinCenters(mpind),1:39,'-k');
        plot(qntlNuLat(1,:),1:30,'-m');
        plot(qntlNuLat(2,:),1:30,'-g');        
        plot(qntlNuLat(3,:),1:30,'-m');
        caxis([0,0.065]);
        
af(@(h) colormap(h,'jet'), sax);
af(@(h) axis(h,'xy'),      sax);
af(@(h) axis(h,'tight'),   sax);
af(@(h) xlim(h,[-250,250]),sax);

%%%>>>

%%%<<< SUPFIG , compare dc2 xy error to dc4 xy error

ind = cf(@(s,u,p,phb)                                            ...
         s(:,1)==1                                               ...
           & any(s(:,3:6),2)                                     ...
           & u>1                                                 ...
           & p>0.005                                            ...
           & phb>0.0005,                                         ...
         dc2.stcm,dc2.uinc,dc2.post,dc4.post);

out = cf(@(e4,e2,i) ...
         hist2(mud([e4(i),e2(i)]),50,50),...
         dc4.error.xy,dc2.error.xy,ind);

figure();
    for s = 1:numTrials+1,
        subplot(1,11,trl),
        if s < numTrials,
            imagesc(out{trl}');
            axis('xy');
            line([0,800],[0,800],'Color','r');    
            title(Trials{trl}.filebase);
        else
            imagesc(sum(cat(3,out{:}),3)');
            axis('xy');
            caxis([0,1000]);
            title('All Sessions');
        end
    end

%%%>>>
    
%%%<<< COMPUTE behavior space distances

ind = cf(@(s,u,p,phb)                                            ...
         s(:,1)==1                                               ...
           & any(s(:,3:6),2)                                     ...
           & u>1                                                 ...
           & p>0.0005                                            ...
           & phb>0.0005,                                         ...
         dc2.stcm,dc2.uinc,dc2.post,dc4.post);



% COMPUTE behavior space distances
ddz        = cf(@(t,u,p)                                        ...
                compute_ddz(t,u,p,'sampleRate',dc2.sampleRate), ...
                Trials, units, dc2.pfs);

% for each trial
    %for each time point:
    %    identify the units within 20cm radius
    %    compute unit composition
    

cFSrC = bsxfun(@plus,FSrC(:,[1:3]),[0.5,1.5,1.5]);
figure,
plot3(cFSrC(:,1),cFSrC(:,2),cFSrC(:,3),'.');
hold('on');
line([0,3],[0,3],[0,3],'Color','r')

rvec = [1/sqrt(numel(cvec)).*ones(size(cvec))]';


% COMPUTE distance between placefield center and maze center
% OR 
% COMPUTE behavior space distances
zang = repmat({[]},size(dc2.xyz));
for t = 1:numel(tind),
    ti = tind(t);
    for s = find(ind{t})',
        usub = units{t}(abs(ddz{t}(s,:))<200);
        zvec = cFSrC( ismember( cluSessionMapSubset,               ...
                                [repmat(ti,[numel(usub),1]),usub'],...
                                'rows'),                           ...
                      1:3);
        if size(zvec,1) < 3,
            zang{t}(end+1) = 0;
            continue;
        else,
            zvec = sum(bsxfun(@rdivide,zvec,sqrt(sum(zvec.^2,2))));
            zang{t}(end+1) = acos((zvec * rvec)./sqrt( sum( zvec.^2, 2)));            
        end
    end
end

%%%>>>

%%%<<< COMPUTE ERROR: DC4

uincInd    = cf(@(u) discretize(u,uincBin),                    dc2.uinc);
xyInd      = cf(@(x) discretize(sq(x(:,'hcom',[1,2])), xyBin), dc2.xyz );
xyInd      = cf(@(x) discretize(sq(x(:,'hcom',[1,2])), xyBin), dc2.xyz );
dc4.hpInd      = cf(@(h) discretize(h(:,1), hpBin), dc4.fet);
dc4.bpInd      = cf(@(b) discretize(b(:,2), bpBin), dc4.fet);

% XY - xy error
dc2.xyEInd = cf(@(e)  discretize(e, xyeBin),  dc2.error.xy);
dc4.xyEInd = cf(@(e)  discretize(e, xyeBin),  dc4.error.xy);
% XYHB - hp and bp error
dc4.hpEInd = cf(@(e)  discretize(e, hpeBin),  dc4.error.hp);
dc4.bpEInd = cf(@(e)  discretize(e, bpeBin),  dc4.error.bp);
dc4.hbEInd = cf(@(e)  discretize(e, hbeBin),  dc4.error.hb);
dc2.uincInd = cf(@(u) discretize(u,uincBin),                    dc2.uinc);

%%%>>>

%%%<<< COMPUTE JPDF : DC2, DC4 decoding error vs unitInc 

out = {};
oe = {};
for sts = 1:6,
    jind = cf(@(ei1,ei2,hpi,bpi,hbi,u,s,bi,hi,po)                     ...
             s(:,1)==1  & any(s(:,sts),2)                          ...
               & nniz(ei1) & nniz(ei2)                             ...
               & nniz(hpi) & nniz(bpi)                             ...
               & nniz(hbi) & nniz(u)                               ...
               & nniz(bi)  & nniz(hi) & po>0.01,                      ...
             dc2.xyEInd, dc4.xyEInd,                               ...
               dc4.hpEInd, dc4.bpEInd,                             ...
               dc4.hbEInd, dc2.uincInd, dc2.stcm,...
               dc4.bpInd,  dc4.hpInd,dc4.post);
% $$$     for trl =  1:10,
% $$$         out{sts}{trl} = hist2([dc4.fet{trl}(jind{trl},2),...
% $$$                             dc2.uinc{trl}(jind{trl})],...
% $$$                               -0.5:0.02:1.8,...
% $$$                               1:40);
% $$$     end
    % mean xy error as function of 
    bpa = cf(@(b,i) b(i), dc4.bpInd  ,jind);
    uia = cf(@(b,i) b(i), dc2.uincInd,jind);
    exa = cf(@(b,i) b(i), dc4.error.xy,jind);    
    oe{sts} = accumarray([cat(1,bpa{:}), cat(1,uia{:})],...
                              cat(1,exa{:}),...
                              [numel(bpBin)-1,numel(uincBin)-1],...
                              @mean,nan);
    
end

out = cf(@(o) sum(cat(3,o{:}),3), out);


figure,
for sts = 1:6,
    subplot(6,1,sts);
    imagesc(bsxfun(@rdivide,out{sts},sum(out{sts}))');
    axis('xy');
end

figure,
for sts = 1:6,
    subplot(6,1,sts);
    imagescnan({bpBin(1:end-1),uincBin(1:end-1),oe{sts}'},[0,150],[],true,'colorMap',@jet);
    axis('xy');
    colormap('jet');
end




    

JPDF_d2xyErr_Uinc = {};
JPDF_d4xyErr_Uinc = {};
JPDF_d4hpErr_Uinc = {};
JPDF_d4bpErr_Uinc = {};
JPDF_d4hbErr_Uinc = {};
JPDF_d4xyErr_d4hbErr = {};
for sts = 1:6
    jind = cf(@(ei1,ei2,hpi,bpi,hbi,u,s,p)                           ...
             s(:,1)==1 & any(s(:,sts),2)                           ...
               & nniz(ei1) & nniz(ei2)                             ...
               & nniz(hpi) & nniz(bpi)                             ...
               & nniz(hbi) & nniz(u) & p>0.01,                              ...
             dc2.xyEInd, dc4.xyEInd,                               ...
               dc4.hpEInd, dc4.bpEInd,                             ...
               dc4.hbEInd, uincInd, dc2.stcm,dc4.post);

    JPDF_d2xyErr_Uinc{sts} =                                       ...
        cf(@(x,u,i)                                                ...
           accumarray([x(i),u(i)],                                 ...
                      1,                                           ...
                      [numel(xyeBin)-1,numel(uincBin)-1],          ...
                      @sum,nan),                                   ...
           dc2.xyEInd, uincInd, jind);

    JPDF_d4xyErr_Uinc{sts} =                                       ...
        cf(@(x,u,i)                                                ...
           accumarray([x(i),u(i)],                                 ...
                      1,                                           ...
                      [numel(xyeBin)-1,numel(uincBin)-1],          ...
                      @sum,nan),                                   ...
           dc4.xyEInd, uincInd, jind);
    
    JPDF_d4hpErr_Uinc{sts} =                                       ...
        cf(@(hp,u,i)                                               ...
           accumarray([hp(i),u(i)],                                ...
                      1,                                           ...
                      [numel(hpeBin)-1,numel(uincBin)-1],          ...
                      @sum,nan),                                   ...
           dc4.bpEInd, uincInd, jind);

    JPDF_d4bpErr_Uinc{sts} =                                       ...
        cf(@(bp,u,i)                                               ...
           accumarray([bp(i),u(i)],                                ...
                      1,                                           ...
                      [numel(bpeBin)-1,numel(uincBin)-1],          ...
                      @sum,nan),                                   ...
           dc4.hpEInd, uincInd, jind);
    
    JPDF_d4hbErr_Uinc{sts} =                                       ...
        cf(@(hb,u,i)                                               ...
           accumarray([hb(i),u(i)],                                ...
                      1,                                           ...
                   [numel(hbeBin)-1,numel(uincBin)-1],             ...
               @sum,nan),                                          ...
       dc4.hbEInd, uincInd, jind);

    JPDF_d4xyErr_d4hbErr{sts} =                                    ...
        cf(@(x,h,i)                                                ...
           accumarray([x(i),h(i)],                                 ...
                      1,                                           ...
                      [numel(xyeBin)-1,numel(hbeBin)-1],           ...
                      @sum,nan),                                   ...
           dc4.xyEInd, dc4.hbEInd, jind);
end

%%%>>>

%%%<<< PLOTING JPDF : DC2, DC4 decoding error vs unitInc

y = 1;

JPDF_d4xyErr_d4hbErr_gfltr = imgaussfilt(sum(cat(3, JPDF_d4xyErr_d4hbErr{:}),3,'omitnan'),2);
JPDF_d4xyErr_d4hbErr_gfltr = JPDF_d4xyErr_d4hbErr_gfltr./sum(nonzeros())
figure();
    colormap('jet');
    imagesc(xyeBinCenters./10,                              ...
            hbeBinCenters,                                 ...
            );
        axis('xy');




xyeBinCenters  = (xyeBin(1:end-1)+xyeBin(2:end))./2;
uincBinCenters = (uincBin(1:end-1)+uincBin(2:end))./2;
hpeBinCenters  = (hpeBin(1:end-1)+hpeBin(2:end))./2;
bpeBinCenters  = (bpeBin(1:end-1)+bpeBin(2:end))./2;
hbeBinCenters  = (hbeBin(1:end-1)+hbeBin(2:end))./2;



figure();
colormap('jet');
for sts = 1:6,
    y = 1;
    subplot2(5,6,y,sts);
        imagesc(xyeBinCenters./10,                              ...
                uincBinCenters,                                 ...
                log10(sum(cat(3, JPDF_d2xyErr_Uinc{sts}{:}),3,'omitnan')+1)')

% $$$         imagesc(xyeBinCenters./10,                              ...
% $$$                 uincBinCenters,                                 ...
% $$$                 bsxfun(@rdivide,...
% $$$                        sum(cat(3, JPDF_d2xyErr_Uinc{sts}{:}),3,'omitnan'),...
% $$$                        sum(sum(cat(3, JPDF_d2xyErr_Uinc{sts}{:}),3,'omitnan')))');
        axis('xy');
% $$$         caxis([0,0.1]);        
        y = y+1;
    subplot2(5,6,y,sts);
        imagesc(xyeBinCenters./10,                              ...
                uincBinCenters,                                 ...
                       log10(sum(cat(3, JPDF_d4xyErr_Uinc{sts}{:}),3,'omitnan')+1)');
    
% $$$         imagesc(xyeBinCenters./10,                              ...
% $$$                 uincBinCenters,                                 ...
% $$$                 bsxfun(@rdivide,...
% $$$                        sum(cat(3, JPDF_d4xyErr_Uinc{sts}{:}),3,'omitnan'),...
% $$$                        sum(sum(cat(3, JPDF_d4xyErr_Uinc{sts}{:}),3,'omitnan')))');
        axis('xy');
% $$$         caxis([0,0.1]);
        y = y+1;                
    subplot2(5,6,y,sts);
        imagesc(hbeBinCenters,                                  ...
                uincBinCenters,                                 ...
                log10(sum(cat(3, JPDF_d4hbErr_Uinc{sts}{:}),3,'omitnan')+1)');
    
% $$$         imagesc(hbeBinCenters,                                  ...
% $$$                 uincBinCenters,                                 ...
% $$$                 bsxfun(@rdivide,...
% $$$                        sum(cat(3, JPDF_d4hbErr_Uinc{sts}{:}),3,'omitnan'),...
% $$$                        sum(sum(cat(3, JPDF_d4hbErr_Uinc{sts}{:}),3,'omitnan')))');
        axis('xy');
% $$$         caxis([0,0.01])
        y = y+1;        
    subplot2(5,6,y,sts);
        imagesc(hpeBinCenters,                                  ...
                uincBinCenters,                                 ...
                log10(sum(cat(3, JPDF_d4hpErr_Uinc{sts}{:}),3,'omitnan')+1)');

% $$$         imagesc(hpeBinCenters,                                  ...
% $$$                 uincBinCenters,                                 ...
% $$$                 bsxfun(@rdivide,...
% $$$                        sum(cat(3, JPDF_d4hpErr_Uinc{sts}{:}),3,'omitnan'),...
% $$$                        sum(sum(cat(3, JPDF_d4hpErr_Uinc{sts}{:}),3,'omitnan')))');
        axis('xy');
% $$$         caxis([0,0.1])
        y = y+1;
    subplot2(5,6,y,sts);
        imagesc(bpeBinCenters,                                  ...
                uincBinCenters,                                 ...
                log10(sum(cat(3, JPDF_d4bpErr_Uinc{sts}{:}),3,'omitnan')+1)');
    
% $$$         imagesc(bpeBinCenters,                                  ...
% $$$                 uincBinCenters,                                 ...
% $$$                 bsxfun(@rdivide,...
% $$$                        sum(cat(3, JPDF_d4bpErr_Uinc{sts}{:}),3,'omitnan'),...
% $$$                        sum(sum(cat(3, JPDF_d4bpErr_Uinc{sts}{:}),3,'omitnan')))');
        axis('xy');
% $$$         caxis([0,0.1])
end

% $$$ t = 9;
% $$$ figure,subplot(211);plot(dc4.fet{t}(:,1)),hold on,plot(dc4.com{t}(:,3)); 
% $$$       subplot(212);plot(dc4.fet{t}(:,2)),hold on,plot(dc4.com{t}(:,4));
% $$$ linkaxes(findobj(gcf(),'Type','Axes'),'x');
       
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

%%%>>>

%%%<<< DCT DECODE xyhb|[phz] 300ms window : 

% decoded position based on high sampleRate and binned by theta phase
dct.sampleRate = sampleRateTPD;
dct.com  = cell([1,numTrials]);
dct.max  = cell([1,numTrials]);
dct.som  = cell([1,numTrials]);
dct.post = cell([1,numTrials]);
dct.uinc = cell([1,numTrials]);
dct.pufr = cell([1,numTrials]);
dct.xyz  = cf(@(t)     resample(preproc_xyz(t,'trb'),dct.sampleRate),      Trials);
for t = 1:numTrials
    Trials{t}.lfp.filename = [Trials{t}.name,'.lfp'];    

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
        spkt.res = round(spkt.res(phzBinInds(spkt.res)==p)./lfp.sampleRate.*dct.sampleRate);
        spkt.sampleRate = dct.sampleRate;
        ufr = Trials{t}.load('ufr',dct.xyz{t},spkt,units{t},spikeWindow,true,'gauss');
        dct.pufr{t}(:,p) = sum(ufr.data,2,'omitnan');
        ufr.data =  ufr.data.*6;
        dct.uinc{t}(:,p) = sum(ufr.data>0.2,2);
        tag = ['xyhb_thpD',num2str(dct.sampleRate),'_',num2str(p),'o',num2str(numel(phzBins)-1)];
        [dct.com{t}(:,:,p),dct.max{t}(:,:,p),...
         dct.sax{t}(:,:,p),dct.post{t}(:,p)] = ...
            decode_ufr(Trials{t},units{t},dct.sampleRate,ufr,dc4.pfs{t},...
                       [],mask,smoothingWeights,'tag',tag,'overwrite',false);
    end
end

%%%>>>

%%%<<< COMPUTE DCT PROJECTION : projection of decoded position onto egocentric coordinates
dct.tRot = {0,0,0,0.17,0.17,0.17,0.17,0.17,0.17,0.17};
dct.fet  = cf(@(t)    fet_HB_pitchB(t,dct.sampleRate),                               Trials       );
dct.hvec = cf(@(x)    x(:,'head_front',[1,2])-x(:,'head_back',[1,2]),                dct.xyz      );
dct.hvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                      dct.hvec     );
dct.hvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                     dct.hvec     );
dct.hvec = cf(@(h,r)  multiprod(h,[cos(r),-sin(r);sin(r),cos(r)],[2,3],[1,2]),       dct.hvec, dct.tRot);
dct.tvec = cf(@(x)    circshift(x(:,'hcom',[1,2]),-1)-circshift(x(:,'hcom',[1,2]),1),dct.xyz      );
dct.tvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                      dct.tvec     );
dct.tvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                     dct.tvec     );
dct.stc  = cf(@(t)    t.load('stc','msnn_ppsvd_raux'),                               Trials       );
dct.stcm = cf(@(s,x)  stc2mat(s,x,states),                                           dct.stc,dct.xyz);
dct.ts   = cf(@(x)    [1:size(x,1)]./x.sampleRate,                                   dct.xyz      );

for t = 1:numel(dct.xyz),
    dct.hvang{t} = filter(copy(dct.xyz{t}),'ButFilter',3,1.5,'low');
% $$$     xycoor = cat(2,...
% $$$                  dct.hvang{t}(:,'hcom',[1,2])-dct.hvang{t}(:,'bcom',[1,2]),...
% $$$                  dct.hvang{t}(:,'nose',[1,2])-dct.hvang{t}(:,'hcom',[1,2]));
% $$$     xycoor = cat(2,...
% $$$                  dct.hbang{t}(:,'spine_upper',[1,2])-dct.hbang{t}(:,'pelvis_root',[1,2]),...
% $$$                  dct.hbang{t}(:,'head_front',[1,2])-dct.hbang{t}(:,'head_back',[1,2]));
    xycoor = cat(2,...
                 dct.hvang{t}(:,'spine_upper',[1,2])-dct.hvang{t}(:,'spine_lower',[1,2]),...
                 dct.hvang{t}(:,'head_front',[1,2])-dct.hvang{t}(:,'head_back',[1,2]));
    dct.hvang{t}.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
% $$$     dct.hvang{t}.data = circ_dist(circshift(dct.hvang{t}.data(:,2),-5),...
% $$$                                   circshift(dct.hvang{t}.data(:,2),+5));
    dct.hvang{t}.data = circ_dist(circshift(circ_dist(dct.hvang{t}.data(:,2),dct.hvang{t}.data(:,1)),-5),...
                                  circshift(circ_dist(dct.hvang{t}.data(:,2),dct.hvang{t}.data(:,1)),+5));

    
    dct.hbang{t} = filter(copy(dct.xyz{t}),'ButFilter',3,15,'low');    
% $$$     xycoor = cat(2,...
% $$$                  dct.hbang{t}(:,'hcom',[1,2])-dct.hbang{t}(:,'bcom',[1,2]),...
% $$$                  dct.hbang{t}(:,'nose',[1,2])-dct.hbang{t}(:,'hcom',[1,2]));
% $$$     xycoor = cat(2,...
% $$$                  dct.hbang{t}(:,'spine_upper',[1,2])-dct.hbang{t}(:,'pelvis_root',[1,2]),...
% $$$                  dct.hbang{t}(:,'head_front',[1,2])-dct.hbang{t}(:,'head_back',[1,2]));
    xycoor = cat(2,...
                 dct.hbang{t}(:,'spine_upper',[1,2])-dct.hbang{t}(:,'spine_lower',[1,2]),...
                 dct.hbang{t}(:,'head_front',[1,2])-dct.hbang{t}(:,'head_back',[1,2]));
    dct.hbang{t}.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    dct.hbang{t}.data = circ_dist(dct.hbang{t}.data(:,2),dct.hbang{t}.data(:,1));
end

for t = 1:numel(dct.xyz),
    dct.hdist{t} = filter(copy(dct.xyz{t}),'ButFilter',3,5,'low');
    xycoor =    dct.xyz{t}(:,'hcom',[1,2]);
    [~,dct.hdist{t}.data] = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
end


for t = 1:numTrials
    hang = transform_origin(Trials{t},dct.xyz{t});
    dct.roll{t} = real(hang.roll);
end

dct.bref = cf(@(t) fet_bref(t,dct.sampleRate), Trials);
for t = 1:numel(dct.bref),
    dct.brefd{t} = filter(copy(dct.bref{t}),'ButFilter',4,[1.5,6],'bandpass');
end


for t = 1:numel(dct.xyz),
    dct.xyvel{t} = vel(filter(copy(dct.xyz{t}),'ButFilter',3,3),'hcom',[1,2]);
end

dct.hrvfl = cf(@(t) filter(fet_href_H(t,dct.sampleRate),'ButFilter',3,25,'low'), Trials);


% COMPUTE error
dct.error = dct.com;
%dErr = dct.sax;
%dErr = dct.max;
dct.error = cf(@(e,x,h,f)  cat(2,                                                                ...
                           sq(multiprod(permute(bsxfun(@minus,                                   ...
                                                       e(:,[1,2],:),                             ...
                                                       sq(x(:,'hcom',[1,2]))),                   ...
                                                [1,2,4,3]),...
                                        h(:,:,:),2,[2,3])),...
                           [bsxfun(@minus,f(:,1),e(:,3,:)),...
                            bsxfun(@minus,f(:,2),e(:,4,:))]), ...
               dct.error,dct.xyz,dct.hvec,dct.fet);
%%%>>>

%%%<<< COMPUTE error for TDP decoding per state

stid = [1,2,3,4,5,6];
dct.lon = repmat({cell([1,numTrials])},[1,numel(stid)]);
dct.lat = repmat({cell([1,numTrials])},[1,numel(stid)]);
dct.hp = repmat({cell([1,numTrials])},[1,numel(stid)]);
dct.bp = repmat({cell([1,numTrials])},[1,numel(stid)]);
dct.hpRnd = repmat({cell([1,numTrials])},[1,numel(stid)]);
dct.bpRnd = repmat({cell([1,numTrials])},[1,numel(stid)]);
dct.ind = cell([1,numel(stid)]);
%norm = 'probability';
norm = 'count';
for sts = 1:numel(stid),
    dct.ind{sts}  = cf(@(p,u,s)                                 ...
                      all(p>0.001,2)                           ...
                      & sum(double(u>=3),2)>6                   ...
                      & s(:,1)==1                               ...
                      & any(logical(s(:,stid(sts))),2)          ...
                      & ~any(logical(s(:,[7,8])),2),            ...
                        dct.post,dct.uinc,dct.stcm);
    for t = 1:numel(Trials),
        for j = 1:numel(phzBins)-1;
            dct.lon{sts}{t} =                                   ...
                cat(2,                                          ...
                    dct.lon{sts}{t},                            ...
                    histcounts(sq(dct.error{t}(dct.ind{sts}{t},1,j)), ...
                               errorBinEdges,...
                               'Normalization',norm)');
            dct.lat{sts}{t} =                                   ...
                cat(2,                                          ...
                    dct.lat{sts}{t},                            ...
                    histcounts(sq(dct.error{t}(dct.ind{sts}{t},2,j)),      ...
                               errorBinEdges,...
                               'Normalization',norm)');                    
            dct.hp{sts}{t} =                                               ...
                cat(2,                                                     ...
                    dct.hp{sts}{t},                                        ...
                    histcounts(sq(dct.error{t}(dct.ind{sts}{t},3,j)),      ...
                               perrorBinEdges,...
                               'Normalization',norm)');                               
            dct.bp{sts}{t} =                                               ...
                cat(2,                                                     ...
                    dct.bp{sts}{t},                                        ...
                    histcounts(sq(dct.error{t}(dct.ind{sts}{t},4,j)),      ...
                               perrorBinEdges,...
                               'Normalization',norm)');

            dct.hpRnd{sts}{t} =                                            ...
                cat(2,                                                     ...
                    dct.hpRnd{sts}{t},                                     ...
                    histcounts(sq(dct.fet{t}(dct.ind{sts}{t},1)) - ...
                               circshift(sq(dct.fet{t}(dct.ind{sts}{t},1)),5000*j),      ...
                               perrorBinEdges,...
                               'Normalization',norm)');
            dct.bpRnd{sts}{t} =                                            ...
                cat(2,                                                     ...
                    dct.bpRnd{sts}{t},                                     ...
                    histcounts(sq(dct.fet{t}(dct.ind{sts}{t},2)) - ...
                               circshift(sq(dct.fet{t}(dct.ind{sts}{t},2)),5000*j),      ...
                               perrorBinEdges,...
                               'Normalization',norm)');
            
        end
            
    end
end

%%%>>>

%%%<<< SUPFIG States x Trials

figure();
if exist('sp','var'),delete(sp);end;
sp = gobjects([numel(stid),numTrials]);
for sts = 1:numel(stid),
    for t = 1:numTrials,
        sp(sts,t) = subplot2(numel(stid),numTrials+1,sts,t);
        imagesc(errorBinCenters,...
                [phzBinCenters;phzBinCenters+2*pi],...
                repmat(dct.lon{sts}{t},[1,2])');
        axis    (sp(sts,t),'xy');
        colormap(sp(sts,t),'jet');
        xlim    (sp(sts,t),[-200,200]);        
    end
end
for sts = 1:numel(stid),
    sp(sts,t+1) = subplot2(numel(stid),numTrials+1,sts,numTrials+1);
    imagesc(errorBinCenters,...
            [phzBinCenters;phzBinCenters+2*pi],...
            repmat(sum(cat(3,dct.lon{sts}{:}),3),[1,2])');
    axis    (sp(sts,t+1),'xy');
    colormap(sp(sts,t+1),'jet');
    xlim    (sp(sts,t+1),[-200,200]);        
end

figure();
if exist('spt','var'),delete(spt);end;
spt = gobjects([numel(stid),numTrials]);
for sts = 1:numel(stid),
    for t = 1:numTrials,
        spt(sts,t) = subplot2(numel(stid),numTrials,sts,t);
        imagesc(errorBinCenters,...
                [phzBinCenters;phzBinCenters+2*pi],...
                repmat(dct.lat{sts}{t},[1,2])');
        axis    (spt(sts,t),'xy');
        colormap(spt(sts,t),'jet');
        xlim    (spt(sts,t),[-200,200]);
    end
end
for sts = 1:numel(stid),
    spt(sts,t+1) = subplot2(numel(stid),numTrials+1,sts,numTrials+1);
    imagesc(errorBinCenters,...
            [phzBinCenters;phzBinCenters+2*pi],...
            repmat(sum(cat(3,dct.lat{sts}{:}),3),[1,2])');
    axis    (spt(sts,t+1),'xy');
    colormap(spt(sts,t+1),'jet');
    xlim    (spt(sts,t+1),[-200,200]);        
end



figure();
if exist('sph','var'),delete(sph);end;
sph = gobjects([numel(stid),numTrials]);
for sts = 1:numel(stid),
    for t = 1:numTrials,
        sph(sts,t) = subplot2(numel(stid),numTrials+1,sts,t);
        imagesc(perrorBinCenters,...
                [phzBinCenters;phzBinCenters+2.*pi],...
                repmat(dct.hp{sts}{t},[1,2])');
        axis    (sph(sts,t),'xy');
        colormap(sph(sts,t),'jet');
    end
end
for sts = 1:numel(stid),
    sph(sts,t+1) = subplot2(numel(stid),numTrials+1,sts,numTrials+1);
    imagesc(perrorBinCenters,...
            [phzBinCenters;phzBinCenters+2*pi],...
            repmat(sum(cat(3,dct.hp{sts}{:}),3),[1,2])');
    axis    (sph(sts,t+1),'xy');
    colormap(sph(sts,t+1),'jet');
end

figure();
if exist('spb','var'),delete(spb);end;
spb = gobjects([numel(stid),numTrials]);
for sts = 1:numel(stid),
    for t = 1:numTrials,
        spb(sts,t) = subplot2(numel(stid),numTrials+1,sts,t);
        imagesc(perrorBinCenters,...
                [phzBinCenters;phzBinCenters+2.*pi],...
                repmat(dct.bp{sts}{t},[1,2])');
        axis    (spb(sts,t),'xy');
        colormap(spb(sts,t),'jet');
    end
end
for sts = 1:numel(stid),
    spb(sts,t+1) = subplot2(numel(stid),numTrials+1,sts,numTrials+1);
    imagesc(perrorBinCenters,...
            [phzBinCenters;phzBinCenters+2*pi],...
            repmat(sum(cat(3,dct.bp{sts}{:}),3),[1,2])');
    axis    (spb(sts,t+1),'xy');
    colormap(spb(sts,t+1),'jet');
end



figure();
if exist('spbr','var'),delete(spbr);end;
spbr = gobjects([numel(stid),numTrials]);
for sts = 1:numel(stid),
    for t = 1:numTrials,
        spbr(sts,t) = subplot2(numel(stid),numTrials+1,sts,t);
        imagesc(perrorBinCenters,...
                [phzBinCenters;phzBinCenters+2.*pi],...
                repmat(dct.bpRnd{sts}{t},[1,2])');
        axis    (spbr(sts,t),'xy');
        colormap(spbr(sts,t),'jet');
    end
end
for sts = 1:numel(stid),
    spbr(sts,t+1) = subplot2(numel(stid),numTrials+1,sts,numTrials+1);
    imagesc(perrorBinCenters,...
            [phzBinCenters;phzBinCenters+2*pi],...
            repmat(sum(cat(3,dct.bpRnd{sts}{:}),3),[1,2])');
    axis    (spbr(sts,t+1),'xy');
    colormap(spbr(sts,t+1),'jet');
end


figure();
if exist('sphr','var'),delete(sphr);end;
sphr = gobjects([numel(stid),numTrials]);
for sts = 1:numel(stid),
    for t = 1:numTrials,
        sphr(sts,t) = subplot2(numel(stid),numTrials+1,sts,t);
        imagesc(perrorBinCenters,...
                [phzBinCenters;phzBinCenters+2.*pi],...
                repmat(dct.hpRnd{sts}{t},[1,2])');
        axis    (sphr(sts,t),'xy');
        colormap(sphr(sts,t),'jet');
    end
end
for sts = 1:numel(stid),
    sphr(sts,t+1) = subplot2(numel(stid),numTrials+1,sts,numTrials+1);
    imagesc(perrorBinCenters,...
            [phzBinCenters;phzBinCenters+2*pi],...
            repmat(sum(cat(3,dct.hpRnd{sts}{:}),3),[1,2])');
    axis    (sphr(sts,t+1),'xy');
    colormap(sphr(sts,t+1),'jet');
end



%%%>>>

%%%<<< VECTORIZE dct data

dErrlon = cf(@(e) reshape( sq(e(:,1,:))', [], 1), dct.error);
dErrlat = cf(@(e) reshape( sq(e(:,2,:))', [], 1), dct.error);
dErrhp  = cf(@(e) reshape( sq(e(:,3,:))', [], 1), dct.error);
dErrbp  = cf(@(e) reshape( sq(e(:,4,:))', [], 1), dct.error);

dbref   = cf(@(e) reshape( repmat( e.data(:,17),[1,8])', [], 1), dct.bref);
dbreff   = cf(@(e) reshape( repmat( e.data(:,17),[1,8])', [], 1), dct.brefd);
dbrefd  = cf(@(e) reshape( repmat( [0;diff(e.data(:,17))],[1,8])', [], 1), dct.brefd);

dhroll  = cf(@(e) reshape( repmat( e     ,[1,8])', [], 1), dct.roll);
dhbang  = cf(@(e) reshape( repmat( e.data,[1,8])', [], 1), dct.hbang);
dhvang  = cf(@(e) reshape( repmat( e.data,[1,8])', [], 1), dct.hvang);
dxyvel  = cf(@(e) reshape( repmat( e.data,[1,8])', [], 1), dct.xyvel);
dhrvf  = cf(@(e) reshape( repmat( e(:,1),[1,8])', [], 1), dct.hrvfl);
dhrvl  = cf(@(e) reshape( repmat( e(:,2),[1,8])', [], 1), dct.hrvfl);
dhdist  = cf(@(e) reshape( repmat( e(:,1),[1,8])', [], 1), dct.hdist);

dstcm   = cf(@(s) reshape( permute( repmat( s, [1,1,8]),[2,3,1]),  8,[])', dct.stcm);
dphz    = cf(@(e) reshape( repmat( phzBinCenters',[size(e,1),1])', [], 1), dct.error);
duinc   = cf(@(u) reshape( u', [], 1), dct.uinc);
dpufr   = cf(@(u) reshape( u', [], 1), dct.pufr);
dpost   = cf(@(p) reshape( p', [], 1), dct.post);
duincI  = cf(@(u) reshape( repmat( sum(double(u>=3),2)>6,[1,8])',[],1), dct.uinc);
dpostI  = cf(@(p) reshape( repmat( all(p>0.0001,2),[1,8])',[],1), dct.post);

dtind  = cf(@(p,t)  repmat(t.*ones([size(p,1),1]),[8,1]), dct.post,num2cell(tind));

dErrlon = cat( 1, dErrlon{:} );
dErrlat = cat( 1, dErrlat{:} );
dErrhp  = cat( 1, dErrhp{:} );
dErrbp  = cat( 1, dErrbp{:} );

dbref   = cat( 1, dbref{:}  );
dbreff  = cat( 1, dbreff{:} );
dbrefd  = cat( 1, dbrefd{:} );

dhroll  = cat( 1, dhroll{:} );
dhbang  = cat( 1, dhbang{:} );
dhvang  = cat( 1, dhvang{:} );
dxyvel  = cat( 1, dxyvel{:} );
dhrvf   = cat( 1, dhrvf{:} );
dhrvl   = cat( 1, dhrvl{:} );
dhdist   = cat( 1, dhdist{:} );

dstcm   = cat( 1, dstcm{:}  );
dphz    = cat( 1, dphz{:}   );
duinc   = cat( 1, duinc{:}  );
dpufr   = cat( 1, dpufr{:}  );
dpost   = cat( 1, dpost{:}  );
duincI  = cat( 1, duincI{:} );
dpostI  = cat( 1, dpostI{:} );
dtind   = cat( 1, dtind{:}  );





tdphz = dphz;


for p = 1:8,
    tdphz(ismember(dtind,[3,4,5])&phzBinCenters(p)==dphz) = ...
        phzBinCenters((p+2).*double(p<=6)+mod(p+2,8).*double(p>6));
end

for p = 1:8,
    tdphz(~ismember(dtind,[3,4,5])&phzBinCenters(p)==dphz) = ...
        phzBinCenters((p+1).*double(p<=7)+mod(p+1,8).*double(p>7));
end
dphz = tdphz;


dphz(~ismember(dtind,[3,4,5])) = dphz(~ismember(dtind,[3,4,5]))+diff(phzBinCenters(1:2));
dphz(dphz>pi) = dphz(dphz>pi)-2.*pi;


save('/storage/gravio/data/project/general/analysis/MjgER2016_decodeBhv.mat',...
     'dErrlon', 'dErrlat', 'dErrhp', 'dErrbp', 'dstcm', 'dphz', 'duinc', 'dpost',...
     'duincI', 'dpostI', 'dtind','-v7.3');
     
load('/storage/gravio/data/project/general/analysis/MjgER2016_decodeBhv.mat');

%%%>>>

%%%<<< COLLECT : ripple associated stats

ts = cf(@(t) ([1:7:round(diff(t.lfp.sync.sync([1,end]),1,2).*t.lfp.sampleRate)]+26)'./t.lfp.sampleRate, Trials);

% set rippleDetectionChannels 

[rpts,rppw] = cf(@(trl,rchans)                                                            ...
    LocalMinima(                                                                          ...
     -mean(                                                                               ...
      nunity(sq(                                                                          ... 
       sum(                                                                               ...
        segs(                                                                             ...
         filter(trl.load('lfp',rchans),                                                   ...
                'ButFilter',                                                              ...
                3,                                                                        ...
                [180,280],                                                                ...
                'bandpass'),                                                              ...
         1:7:round(diff(trl.lfp.sync.sync([1,end]),1,2).*trl.lfp.sampleRate),             ...
         51,                                                                              ...
         nan).^2)),                                                                       ...
       [],[],[],[5,95]),                                                                  ...
      2),                                                                                 ...
     5,                                                                                   ...
    -3),                                                                                  ...
    Trials,                                                                               ...
    mat2cell(cat(1,sessionList(tind).rippleDetectionChannels)',8,ones(size(tind)))        ...
);

%
[xx,xi,yi] = cf(@(dts,rts) NearestNeighbour(dts,rts), dct.ts, rpts);

[xx,xi,yi] = cf(@(dts,rts) NearestNeighbour(dts,rts), dc2.ts, rpts);

 
figure();
hold('on');
cf(@(u,s,r,x) plot(u(x)+randn([numel(x),1]),r,'.'), dc2.uinc, dc2.stcm, rppw, xi);
 
[usi,usv] = cf(@(u) LocalMinima(-u,20,-16), dc2.uinc);
 


rpow = cf(                                                                                ...
    @(trl,rchans)                                                                         ...
    mean(                                                                                 ...
      nunity(sq(                                                                          ...
       sqrt(                                                                              ... 
        sum(                                                                              ...
         segs(                                                                            ...
          filter(trl.load('lfp',rchans),                                                  ...
                'ButFilter',                                                              ...
                 3,                                                                       ...
                 [180,280],                                                               ...
                'bandpass'),                                                              ...
          1:7:round(diff(trl.lfp.sync.sync([1,end]),1,2).*trl.lfp.sampleRate),            ...
          51,                                                                             ...
          nan).^2))),                                                                     ...
       [],[],[],[5,95]),                                                                  ...
     2),                                                                                  ...
    Trials,                                                                               ...
    mat2cell(cat(1,sessionList(tind).rippleDetectionChannels)',8,ones(size(tind)))        ...
);
 

[xx,xi,yi] = cf(@(dts,rts) NearestNeighbour(dts,rts), ts, dc2.ts);

figure();
hold('on');
cf(@(u,s,r,x,o) plot((u+randn([numel(u),1])).*double(any(logical(s(:,1:6)),2)),(r(x)+o).*double(any(logical(s(:,1:6)),2)),'.'), dc2.uinc, dc2.stcm, rpow, xi,num2cell(100:100:1000));

 
sclr = jet(6);
figure();
hold('on');
cf(@(u,s,r,x,o) ...
   scatter((u+randn([numel(u),1])).*double(any(logical(s(:,1:6)),2)),...
           (r(x)+o).*double(any(logical(s(:,1:6)),2)),...
           5,...
           sclr(sum(s(any(logical(s(:,1:6)),2),3:6),2),:)), ...
   dc2.uinc, dc2.stcm, rpow, xi,num2cell(100:100:1000));

 
sclr = jet(6);
o = 100:100:1000;
outT = {}
outN = {} 
 
figure();
for t = 1:numel(dc2.uinc),
    s = dc2.stcm{t};
    u = dc2.uinc{t};
    for sts = 2:8,
        subplot2(8,2,sts,1);
        hold('on');
        ind = logical(s(:,sts)) & logical(s(:,1));
% $$$         scatter(u(ind)+randn([sum(ind),1])/5,...
% $$$                 log10(rpow{t}(xi{t}(ind))+3),...
% $$$                 5,...
% $$$                 sclr(1,:),...
% $$$                 'filled');
        outT{sts}{t} =hist2([u(ind),...
               log10(rpow{t}(xi{t}(ind))+3)],...
              1:30,...
              linspace(-0.5,1.5,40));
        
        subplot2(8,2,sts,2);
        hold('on');
        ind = logical(s(:,sts)) & ~logical(s(:,1));
% $$$         scatter(u(ind)+randn([sum(ind),1])/5,...
% $$$                log10(rpow{t}(xi{t}(ind))+3),...
% $$$                 5,...
% $$$                 sclr(1,:),...
% $$$                 'filled');
        outN{sts}{t} = hist2([u(ind),...
               log10(rpow{t}(xi{t}(ind))+3)],...
              1:30,...
              linspace(-0.5,1.5,40));
        
    end
end 
linkaxes(findobj(gcf(),'Type','Axes'),'xy');


figure()
    for sts = 2:8,
        subplot2(8,2,sts,1);    
        imagesc(1:30,linspace(-0.5,1.5,40),sum(cat(3,outT{sts}{:}),3)');
        axis('xy');
        subplot2(8,2,sts,2);    
        imagesc(1:30,linspace(-0.5,1.5,40),sum(cat(3,outN{sts}{:}),3)');
        axis('xy');        
    end
linkaxes(findobj(gcf(),'Type','Axes'),'xy');

 
figure();
hold('on');
cf(@(u,s,r,x,o) plot((u+randn([numel(u),1])),(r(x)+o),'.'), dc2.uinc, dc2.stcm, rpow, xi,num2cell(100:100:1000));

%%%>>>

% LNO vs Phz
%%%<<< MUTUAL INFORMATION between longitudinal decoding offset and theta phase

Ixy  = zeros([100,6,6]);
Ixyr = zeros([100,6,6]);
for i = 1:100,
    disp(['i: ',num2str(i)]);
    for sts = 1:6,
        ind =   logical(dstcm(:,1))                             ...
                & logical(dstcm(:,sts))                         ...
                & ~any(logical(dstcm(:,[7,8])),2)               ...
                & dpostI                                        ...
                & duincI;        
        ind(ind) = randn([sum(ind),1]) > 2;        
        si = sum(ind);
        for sto = 1:6,    
            tic
            indo =   logical(dstcm(:,1))                        ...
                     & logical(dstcm(:,sto))                    ...
                     & ~any(logical(dstcm(:,[7,8])),2)          ...
                     & dpostI                                   ...
                     & duincI;        
            indo(indo) = randn([sum(indo),1]) > 2;
            so = sum(indo);
            
            indb = find(ind | indo);
            indb = randsample(indb,sum([si,so].*double([sts>sto,sts<=sto])));
            
            px  = histcounts(dErrlon(indb),                     ...
                             errorBinEdges,                     ...
                             'Normalization','probability');
            py  = histcounts(dphz(indb),                        ...
                             phzBins,                           ...
                             'Normalization','probability');
            pxyt = histcounts2(dErrlon(indb),                    ...
                              dphz(indb),                       ...
                              errorBinEdges,                    ...
                              phzBins,                          ...
                              'Normalization','probability');
            Ixy(i,sts,sto) = sum(nonzeros(pxy.*log(pxy./bsxfun(@times,px',py))),'omitnan');
            
            rndPhz = dphz(indb)+randn([numel(indb),1])*2;
            rndPhz( rndPhz < -pi ) = rndPhz( rndPhz < -pi ) + 2.*pi;
            rndPhz( rndPhz >  pi ) = rndPhz( rndPhz >  pi ) - 2.*pi;
            py  = histcounts(rndPhz,                            ...
                             phzBins,                           ...
                             'Normalization','probability');
            pxy = histcounts2(dErrlon(indb),                    ...
                              rndPhz,                           ...
                              errorBinEdges,                    ...
                              phzBins,                          ...
                              'Normalization','probability');
            Ixyr(i,sts,sto) = sum(nonzeros(pxy.*log(pxy./bsxfun(@times,px',py))),'omitnan');
            toc
        end
    end
end
% $$$ save('/storage/gravio/data/project/general/analysis/MjgER2016_decodeBhv_Ixy.mat',...
% $$$      'Ixy','Ixyr');
load('/storage/gravio/data/project/general/analysis/MjgER2016_decodeBhv_Ixy.mat');
     
figure,hold('on');
histogram(log10(Ixyr(:,1)),linspace(-4,-1,100))
hist(log10(Ixy(:,1)),linspace(-4,-1,100))

%%%>>>

% HPO vs PHZ
%%%<<< MUTUAL INFORMATION between head pitch decoding offset and theta phase

Ihp  = zeros([100,6,6]);
Ihpr = zeros([100,6,6]);
for i = 1:100,
    disp(['i: ',num2str(i)]);
    for sts = 1:6,
        ind =   logical(dstcm(:,1))                             ...
                & logical(dstcm(:,sts))                         ...
                & ~any(logical(dstcm(:,[7,8])),2)               ...
                & dpostI                                        ...
                & duincI;        
        ind(ind) = randn([sum(ind),1]) > 1;        
        si = sum(ind);
        for sto = 1:6,    
            tic
            indo =   logical(dstcm(:,1))                        ...
                     & logical(dstcm(:,sto))                    ...
                     & ~any(logical(dstcm(:,[7,8])),2)          ...
                     & dpostI                                   ...
                     & duincI;        
            indo(indo) = randn([sum(indo),1]) > 1;
            so = sum(indo);
            
            indb = find(ind | indo);
            indb = randsample(indb,sum([si,so].*double([sts>sto,sts<=sto])));
            
            px  = histcounts(dErrhp(indb),                     ...
                             perrorBinEdges,                    ...
                             'Normalization','probability');
            py  = histcounts(dphz(indb),                        ...
                             phzBins,                           ...
                             'Normalization','probability');
            pxy = histcounts2(dErrhp(indb),                    ...
                              dphz(indb),                       ...
                              perrorBinEdges,                   ...
                              phzBins,                          ...
                              'Normalization','probability');
            Ihp(i,sts,sto) = sum(nonzeros(pxy.*log(pxy./bsxfun(@times,px',py))),'omitnan');
            
            rndPhz = dphz(indb)+randn([numel(indb),1])*2;
            rndPhz( rndPhz < -pi ) = rndPhz( rndPhz < -pi ) + 2.*pi;
            rndPhz( rndPhz >  pi ) = rndPhz( rndPhz >  pi ) - 2.*pi;
            py  = histcounts(rndPhz,                            ...
                             phzBins,                           ...
                             'Normalization','probability');
            pxy = histcounts2(dErrhp(indb),                    ...
                              rndPhz,                           ...
                              perrorBinEdges,                    ...
                              phzBins,                          ...
                              'Normalization','probability');
            Ihpr(i,sts,sto) = sum(nonzeros(pxy.*log(pxy./bsxfun(@times,px',py))),'omitnan');
            toc
        end
    end
end

% $$$ save('/storage/gravio/data/project/general/analysis/MjgER2016_decodeBhv_Ihp.mat',...
% $$$      'Ihp','Ihpr');
load('/storage/gravio/data/project/general/analysis/MjgER2016_decodeBhv_Ihp.mat');

%%%>>>

% BPO vs PHZ
%%%<<< MUTUAL INFORMATION between head pitch decoding offset and theta phase

Ibp  = zeros([100,6,6]);
Ibpr = zeros([100,6,6]);
for i = 1:100,
    disp(['i: ',num2str(i)]);
    for sts = 1:6,
        ind =   logical(dstcm(:,1))                             ...
                & logical(dstcm(:,sts))                         ...
                & ~any(logical(dstcm(:,[7,8])),2)               ...
                & dpostI                                        ...
                & duincI;        
        ind(ind) = randn([sum(ind),1]) > 1;        
        si = sum(ind);
        for sto = 1:6,    
            tic
            indo =   logical(dstcm(:,1))                        ...
                     & logical(dstcm(:,sto))                    ...
                     & ~any(logical(dstcm(:,[7,8])),2)          ...
                     & dpostI                                   ...
                     & duincI;        
            indo(indo) = randn([sum(indo),1]) > 1;
            so = sum(indo);
            
            indb = find(ind | indo);
            indb = randsample(indb,sum([si,so].*double([sts>sto,sts<=sto])));
            
            px  = histcounts(dErrbp(indb),                     ...
                             perrorBinEdges,                    ...
                             'Normalization','probability');
            py  = histcounts(dphz(indb),                        ...
                             phzBins,                           ...
                             'Normalization','probability');
            pxy = histcounts2(dErrbp(indb),                    ...
                              dphz(indb),                       ...
                              perrorBinEdges,                   ...
                              phzBins,                          ...
                              'Normalization','probability');
            Ibp(i,sts,sto) = sum(nonzeros(pxy.*log(pxy./bsxfun(@times,px',py))),'omitnan');
            
            rndPhz = dphz(indb)+randn([numel(indb),1])*2;
            rndPhz( rndPhz < -pi ) = rndPhz( rndPhz < -pi ) + 2.*pi;
            rndPhz( rndPhz >  pi ) = rndPhz( rndPhz >  pi ) - 2.*pi;
            py  = histcounts(rndPhz,                            ...
                             phzBins,                           ...
                             'Normalization','probability');
            pxy = histcounts2(dErrbp(indb),                    ...
                              rndPhz,                           ...
                              perrorBinEdges,                    ...
                              phzBins,                          ...
                              'Normalization','probability');
            Ibpr(i,sts,sto) = sum(nonzeros(pxy.*log(pxy./bsxfun(@times,px',py))),'omitnan');
            toc
        end
    end
end
%save('/storage/gravio/data/project/general/analysis/MjgER2016_decodeBhv_Ibp.mat',...
%     'Ibp','Ibpr');
load('/storage/gravio/data/project/general/analysis/MjgER2016_decodeBhv_Ibp.mat');

% matlab -nojvm -nodisplay -nosplash 

%%%>>>

%%%<<< RESAMPLE : thetaXerror distributions 

% $$$ indTDP{sts}  = cf(@(p,u,s)                                           ...
% $$$                   all(p>0.0005,2)                                    ...
% $$$                   & sum(double(u>=3),2)>6                            ...
% $$$                   & s(:,1)==1                                        ...
% $$$                   & any(logical(s(:,stid(sts))),2)                   ...
% $$$                   & ~any(logical(s(:,[7,8])),2),                     ...
% $$$                   posteriorMaxTPD,unitInclusionTPD,stcm);


numIter = 100;
ferrorBinEdges = {errorBinEdges,errorBinEdges,perrorBinEdges,perrorBinEdges};
ferrorBinCenters = {errorBinCenters,errorBinCenters,perrorBinCenters,perrorBinCenters};
ferr = {dErrlon,dErrlat,dErrhp,dErrbp};
ttid = 1:7;

pxp = zeros( [ numel( errorBinCenters ), ...
               numel( phzBinCenters   ), ...
               numIter,                  ...
               numel(stid),              ...
               numel(ttid)+1,            ...
               numel(ferrorBinEdges)     ...
             ]                           ...
           );

smMask = false([size(dphz,1)./numel(phzBinCenters),1]);
smMask(1:32:end) = true;
smMask = reshape(repmat(smMask,[1,numel(phzBinCenters)])',[],1);



shifts = 0:8:2^8;
for sts = 1:6,
    ind =   logical(dstcm(:,1))                             ...
            & logical(dstcm(:,sts))                         ...
            & ~any(logical(dstcm(:,[7,8])),2)               ...
            & dpostI                                        ...
            & duincI;        
    for t = [1:numel(ttid),numel(ttid)+1],
        if t <= numel(ttid)
            tid = tind(ttid(t));
        end
        
        for i = 1:100,
            tic        
            disp(['s: ',num2str(sts),'  t: ',num2str(t),'  i: ',num2str(i)]);
            if t == numel(ttid)+1,
                indb = ind & circshift(smMask,randsample(shifts,1));            
            else,
                indb = ind & circshift(smMask,randsample(shifts,1)) & ( dtind == tid );
            end
            indb(indb) = randn([sum(indb),1]) > 0;
            for f = 1:4
                pxp(:,:,i,sts,t,f) = ...
                    histcounts2(ferr{f}(indb),                  ...
                                dphz(indb),                     ...
                                ferrorBinEdges{f},              ...
                                phzBins,                        ...
                                'Normalization','probability');
            end
            toc
        end
    end
end

npxp = bsxfun(@rdivide,pxp,sum(pxp));

qtls = [0.25,0.50,0.75];
% QNTL( qtl, phz, itr, sts, trl, fet)
qntl = zeros([3,size(npxp,2),size(npxp,3),size(npxp,4),size(npxp,5),size(npxp,6)]);
for p = 1:size(npxp,2),
    for i = 1:size(npxp,3),
        for s =1:size(npxp,4),
            for t = 1:size(npxp,5),
                for f = 1:size(npxp,6),
                    [x, index] = unique(cumsum(npxp(:,p,i,s,t,f)));
                    qntl(:,p,i,s,t,f) = interp1(x,ferrorBinCenters{f}(index),qtls);
                end
            end
        end
    end
end

% MPX( qtl, phz, itr, sts, trl, fet)

mpxp = mean(npxp,3,'omitnan');

%%%>>>

%%%<<< SUPFIG : DHPR head referenced projection of decoded position 

lnstyle = {'--m','-k','--m'};

clims = {[0,0.05],[0,0.05],[0,0.075],[0,0.1]};

%pboundaries = {[-200,300],[-200,200],[-0.75,0.75],[-0.75,0.75]};
pboundaries = {[-200,300],[-200,200],[-1,1],[-1,1]};


for sts = 1:6;
figure,
hax = gobjects([8,4]);
for f = 1:4
    for trl = 1:8,
        hax(trl,1) = subplot2(4,8,f,trl);
        hold(hax(trl,1),'on');
        imagesc(ferrorBinCenters{f},...
                [phzBinCenters;phzBinCenters+2.*pi],...
                repmat(mpxp(:,:,1,sts,trl,f),[1,2,1])');
        hax(trl,1).YTick = [-pi,0,pi,2*pi,3.*pi];
        hax(trl,1).YTickLabels = {'-\pi','0','\pi','2\pi','3\pi'};
        axis    (hax(trl,1),'xy');
        colormap(hax(trl,1),'jet');
        caxis(clims{f});

        ylim([-pi,3.*pi]);
        for q = 1:3,
            plot(repmat(sq(mean(qntl(q,:,:,sts,trl,f),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],lnstyle{q});
        end

        if trl == 8
            for t = 1:8;
                for q = 1:3,
                    plot(repmat(sq(mean(qntl(q,:,:,sts,t,f),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],lnstyle{q});
                end
            end
        end
        xlim(pboundaries{f}([1,end]));        
    end
end
suptitle(states{sts});
end




figure,
hax = gobjects([6,4]);
trl = 8;
for sts = 1:6;
    for f = 1:4
        hax(sts,f) = subplot2(4,6,f,sts);
        hold(hax(sts,f),'on');
        imagesc(ferrorBinCenters{f},...
                [phzBinCenters;phzBinCenters+2.*pi],...
                repmat(mpxp(:,:,1,sts,trl,f),[1,2,1])');
        hax(sts,f).YTick       = [-pi,0,pi,2*pi,3.*pi];
        hax(sts,f).YTickLabels = {'-\pi','0','\pi','2\pi','3\pi'};
        axis    (hax(sts,f),'xy');
        colormap(hax(sts,f),'jet');
        caxis(clims{f});
        ylim([-pi,3.*pi]);
        for t = 1:8;
            for q = 1:3,
                plot(repmat(sq(mean(qntl(q,:,:,sts,t,f),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],lnstyle{q});
            end
        end
        if f ==1, title(states{sts});end
        xlim(pboundaries{f}([1,end]));
    end
end





for sts = 1:8,
hax(sts,2) = subplot2(4,8,2,sts);
hold(hax(sts,2),'on');
imagesc(errorBinCenters,...
        [phzBinCenters;phzBinCenters+2.*pi],...
        repmat(mplt(:,:,5,sts),[1,2,1])');
hax(sts,2).YTick = [-pi,0,pi,2*pi,3.*pi];
hax(sts,2).YTickLabels = {'-\pi','0','\pi','2\pi','3\pi'};
axis    (hax(sts,2),'xy');
colormap(hax(sts,2),'jet')
caxis([0,0.06])
xlim([-200,300]);
ylim([-pi,3.*pi]);
plot(repmat(sq(mean(q25lt(:,:,:,5,8),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
plot(repmat(sq(mean(q50lt(:,:,:,5,8),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'-k');
plot(repmat(sq(mean(q75lt(:,:,:,5,8),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
if sts == 1
    for t = 1:8;
    plot(repmat(sq(mean(q25lt(:,:,:,5,t),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
    plot(repmat(sq(mean(q50lt(:,:,:,5,t),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'-k');
    plot(repmat(sq(mean(q75lt(:,:,:,5,t),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
end
end
end

for sts = 1:8,
hax(sts,3) = subplot2(4,8,3,sts);
hold(hax(sts,3),'on');
imagesc(perrorBinCenters,...
        [phzBinCenters;phzBinCenters+2.*pi],...
        repmat(mphp(:,:,5,sts),[1,2,1])');
hax(sts,3).YTick = [-pi,0,pi,2*pi,3.*pi];
hax(sts,3).YTickLabels = {'-\pi','0','\pi','2\pi','3\pi'};
axis    (hax(sts,3),'xy');
colormap(hax(sts,3),'jet')
caxis([0,0.06])
xlim([pbound]);
ylim([-pi,3.*pi]);
plot(repmat(sq(mean(q25hp(:,:,:,5,8),3,'omitnan')),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
plot(repmat(sq(mean(q50hp(:,:,:,5,8),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'-k');
plot(repmat(sq(mean(q75hp(:,:,:,5,8),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
if sts == 1
    for t = 1:8;
    plot(repmat(sq(mean(q25hp(:,:,:,5,t),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
    plot(repmat(sq(mean(q50hp(:,:,:,5,t),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'-k');
    plot(repmat(sq(mean(q75hp(:,:,:,5,t),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
end
end
end


for sts = 1:8,
hax(sts,4) = subplot2(4,8,4,sts);
hold(hax(sts,4),'on');
imagesc(perrorBinCenters,...
        [phzBinCenters;phzBinCenters+2.*pi],...
        repmat(mpbp(:,:,5,sts),[1,2,1])');
hax(sts,4).YTick = [-pi,0,pi,2*pi,3.*pi];
hax(sts,4).YTickLabels = {'-\pi','0','\pi','2\pi','3\pi'};
axis    (hax(sts,4),'xy');
colormap(hax(sts,4),'jet')
caxis([0,0.06])
xlim([pbound]);
ylim([-pi,3.*pi]);
plot(repmat(sq(mean(q25bp(:,:,:,5,8),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
plot(repmat(sq(mean(q50bp(:,:,:,5,8),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'-k');
plot(repmat(sq(mean(q75bp(:,:,:,5,8),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
if sts == 1
    for t = 1:8;
    plot(repmat(sq(mean(q25bp(:,:,:,5,t),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
    plot(repmat(sq(mean(q50bp(:,:,:,5,t),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'-k');
    plot(repmat(sq(mean(q75bp(:,:,:,5,t),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
end
end
end



figure();
hax = gobjects([6,4]);
for sts = 1:6,
hax(sts,1) = subplot2(6,4,sts,1);
hold(hax(sts,1),'on');
imagesc(errorBinCenters,...
        [phzBinCenters;phzBinCenters+2.*pi],...
        repmat(mpxy(:,:,sts,8),[1,2,1])');
hax(sts,1).YTick = [-pi,0,pi,2*pi,3.*pi];
hax(sts,1).YTickLabels = {'-\pi','0','\pi','2\pi','3\pi'};
axis    (hax(sts,1),'xy');
colormap(hax(sts,1),'jet')
caxis([0,0.045])
xlim([-200,300]);
ylim([-pi,3.*pi]);
plot(repmat(sq(mean(q25xy(:,:,:,sts,1),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
plot(repmat(sq(mean(q25xy(:,:,:,sts,1),3)),[1,2])- ...
     3.*repmat(sq(std(q25xy(:,:,:,sts,1),[],3)),[1,2]),...
     [phzBinCenters;phzBinCenters+2.*pi],'--m');
plot(repmat(sq(mean(q25xy(:,:,:,sts,1),3)),[1,2])+ ...
     3.*repmat(sq(std(q25xy(:,:,:,sts),[],3)),[1,2]),...
     [phzBinCenters;phzBinCenters+2.*pi],'--m');
plot(repmat(sq(mean(q50xy(:,:,:,sts,1),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'-k');
plot(repmat(sq(mean(q50xy(:,:,:,sts,1),3)),[1,2])- ...
     3.*repmat(sq(std(q50xy(:,:,:,sts,1),[],3)),[1,2]),...
     [phzBinCenters;phzBinCenters+2.*pi],'--g');
plot(repmat(sq(mean(q50xy(:,:,:,sts,1),3)),[1,2])+ ...
     3.*repmat(sq(std(q50xy(:,:,:,sts,1),[],3)),[1,2]),...
     [phzBinCenters;phzBinCenters+2.*pi],'--g');
plot(repmat(sq(mean(q75xy(:,:,:,sts,1),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
plot(repmat(sq(mean(q75xy(:,:,:,sts,1),3)),[1,2])- ...
     3.*repmat(sq(std(q75xy(:,:,:,sts,1),[],3)),[1,2]),...
     [phzBinCenters;phzBinCenters+2.*pi],'--m');
plot(repmat(sq(mean(q75xy(:,:,:,sts,1),3)),[1,2])+ ...
     3.*repmat(sq(std(q75xy(:,:,:,sts,1),[],3)),[1,2]),...
     [phzBinCenters;phzBinCenters+2.*pi],'--m');
end

for sts = 1:6,
hax(sts,1) = subplot2(6,4,sts,2);
hold(hax(sts,1),'on');
imagesc(errorBinCenters,...
        [phzBinCenters;phzBinCenters+2.*pi],...
        repmat(mplt(:,:,sts),[1,2,1])');
hax(sts,1).YTick = [-pi,0,pi,2*pi,3.*pi];
hax(sts,1).YTickLabels = {'-\pi','0','\pi','2\pi','3\pi'};
axis    (hax(sts,1),'xy');
colormap(hax(sts,1),'jet')
caxis([0,0.045])
xlim([-200,200]);
ylim([-pi,3.*pi]);
plot(repmat(sq(mean(q25lt(:,:,:,sts),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
plot(repmat(sq(mean(q50lt(:,:,:,sts),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'-k');
plot(repmat(sq(mean(q75lt(:,:,:,sts),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
end

for sts = 1:6,
hax(sts,2) = subplot2(6,4,sts,3);
hold(hax(sts,2),'on');
imagesc(perrorBinCenters,...
        [phzBinCenters;phzBinCenters+2.*pi],...
        repmat(mphp(:,:,sts),[1,2,1])');
hax(sts,2).YTick = [-pi,0,pi,2*pi,3.*pi];
hax(sts,2).YTickLabels = {'-\pi','0','\pi','2\pi','3\pi'};
axis(hax(sts,2),'xy');
colormap(hax(sts,2),'jet')
caxis([0,0.075])
xlim([-1.25,1.25]);
ylim([-pi,3.*pi]);
plot(repmat(sq(mean(q25hp(:,:,:,sts),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
plot(repmat(sq(mean(q50hp(:,:,:,sts),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'-k');
plot(repmat(sq(mean(q75hp(:,:,:,sts),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
end


for sts = 1:6,
hax(sts,3) = subplot2(6,4,sts,4);
hold(hax(sts,3),'on');
imagesc(perrorBinCenters,...
        [phzBinCenters;phzBinCenters+2.*pi],...
        repmat(mpbp(:,:,sts),[1,2,1])');
hax(sts,3).YTick = [-pi,0,pi,2*pi,3.*pi];
hax(sts,3).YTickLabels = {'-\pi','0','\pi','2\pi','3\pi'};
axis(hax(sts,3),'xy');
colormap(hax(sts,3),'jet');
caxis([0,0.1]);
xlim([-1.25,1.25]);
ylim([-pi,3.*pi]);
plot(repmat(sq(mean(q25bp(:,:,:,sts),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
plot(repmat(sq(mean(q50bp(:,:,:,sts),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'-k');
plot(repmat(sq(mean(q75bp(:,:,:,sts),3)),[1,2]),[phzBinCenters;phzBinCenters+2.*pi],'--k');
end

%%%>>>

%%%<<< COMPUTE : unit coactivation as function of theta phase 

stsGrpId = {2,[3,5],[4,6],1};
pxz = nan([numel(ferrorBinCenters{1}),...
           numel(phzBinCenters),...
           100,...
           numel(stsGrpId),...
           numel(ttid),...
           numel(ferrorBinCenters),...
           2]);
for sts = 4:numel(stsGrpId),
    ind =   logical(dstcm(:,1))                             ...
            & any(logical(dstcm(:,stsGrpId{sts})),2)        ...
            & ~any(logical(dstcm(:,[7,8])),2)               ...
            & dpostI                                        ...
            & duincI;        
    for t = [1:numel(ttid),numel(ttid)+1],
        if t <= numel(ttid)
            tid = tind(ttid(t));
        end
        for i = 1:100,
            if t == numel(ttid)+1,                
                indb = find(ind & circshift(smMask,randsample(1:32,1)) & duinc < 7);
            else
                indb = find(ind & circshift(smMask,randsample(1:32,1)) & duinc < 7 &  dtind == tid);
            end
            indb = randsample(indb,round(numel(indb)/2));
            for f = 1:4,
                pxz(:,:,i,sts,t,f,1) = ...
                    histcounts2(ferr{f}(indb),                  ...
                                dphz(indb),                     ...
                                ferrorBinEdges{f},              ...
                                phzBins,                        ...
                                'Normalization','probability');
            end
            
            if t == numel(ttid)+1,                                
                indb = find(ind & circshift(smMask,randsample(1:32,1)) & duinc >= 7);
            else
                indb = find(ind & circshift(smMask,randsample(1:32,1)) & duinc >= 7 & dtind == tid);
            end
            indb = randsample(indb,round(numel(indb)/2));            
            for f = 1:4,
                pxz(:,:,i,sts,t,f,2) = ...
                    histcounts2(ferr{f}(indb),                  ...
                                dphz(indb),                     ...
                                ferrorBinEdges{f},              ...
                                phzBins,                        ...
                                'Normalization','probability');
            end

        end
    end
end

%%%>>>

%%%<<< SUPFIG : unit coativation as function of theta phase

% MPX( qtl, phz, uincPart, sts, trl, fet)
figure();
for sts = 1:3,
    subplot2(6,2,sts,1);
        imagesc(errorBinCenters,...
                [phzBinCenters;phzBinCenters+2.*pi],...
                repmat(mean(pxz(:,:,:,sts,8,1,1),3),[1,2])');
        axis('xy');
    subplot2(6,2,sts,2);
        imagesc(errorBinCenters,...
                [phzBinCenters;phzBinCenters+2.*pi],...
                repmat(mean(pxz(:,:,:,sts,8,1,2),3),[1,2])');
        axis('xy');
end
colormap('jet');
ForAllSubplots('Lines(0,[],''k'');');
ForAllSubplots('Lines([],2*pi,''k'');');
ForAllSubplots('Lines([],pi,''k'');');
ForAllSubplots('Lines([],0,''k'');');
ForAllSubplots('caxis([0,0.01]);');
ForAllSubplots('caxis([0,0.01]);');
ForAllSubplots('xlim([-300,300]);');

% SMOOTH jpdf
spxz = imgaussfilt(repmat(pxz,[1,3,1,1,1,1,1]),[1,0.5]);
spxz = spxz(:,[1:numel(phzBinCenters)]+numel(phzBinCenters),:,:,:,:,:);
% NORMALIZE jpdf
spxz = reshape(bsxfun(@rdivide,...
                      reshape(spxz,[],1,size(spxz,3),...
                                   size(spxz,4),size(spxz,5),size(spxz,6),size(spxz,7)),...
                      sum(sum(spxz))),...
               size(pxz));

% GET Contour 
coords = cell([1,2]);
[coords{:}] = ndgrid(ferrorBinCenters{1},[phzBinCenters;phzBinCenters+2.*pi]);


plevels = [0.50];
cclr= 'rgk';
figure();
    for sts = 1:3
    subplot2(1,2,1,1); 
        mspxz = repmat(mean(spxz(:,:,:,sts,8,1,1),3),[1,2]);
        if sts == 1,
            imagesc(errorBinCenters,...
                    [phzBinCenters;phzBinCenters+2.*pi],...
                    repmat(mean(spxz(:,:,:,size(spxz,4),8,1,1),3),[1,2])');
        end
        hold('on');

        %plevels = [0.25,0.50,0.75];

        clevels = zeros(size(plevels));;
        maxIter = 10;
        for p = 1:numel(plevels),
            step = 0.0001;
            thresh = 0.0001;        
            iter = 1;
            while round(sum(mspxz(mspxz(:)>(thresh))).*100)~=plevels(p).*100 && iter<maxIter,
                while sum(mspxz(mspxz(:)>(thresh)))>plevels(p),
                    thresh = thresh+step;
                end
                step = step./2;
                while sum(mspxz(mspxz(:)>(thresh)))<plevels(p),
                    thresh = thresh-step;            
                end
                step = step./2;
                iter = iter+1;
            end
            clevels(p) = thresh;
        end
        
        if numel(clevels) == 1,
            clevels = repmat(clevels,[1,2]);
        end
        
        [~,ctax] = contour(coords{:},mspxz,clevels);
        ctax.LineWidth = 1;
        ctax.LineColor = cclr(sts);

        ylim(circ_ang2rad([-60,420]));
        axis('xy');

        xl = [-200,200];
        xlim(xl);
        % MASK areas outside single cycle
        patch([xl(1),xl(1),xl(2),xl(2)],...
              circ_ang2rad([ 360, 420,420,360]),...
              [0.2,0.2,0.2],...
              'EdgeAlpha',0,...
              'FaceAlpha',0.4);
        patch([xl(1),xl(1),xl(2),xl(2)],...    
              circ_ang2rad([ -60,   0,  0,-60]),...
              [0.2,0.2,0.2],...
              'EdgeAlpha',0,...
              'FaceAlpha',0.4);
        sax(end).XTick = [-100,0,100];
        sax(end).XTickLabels = [-15,0,15];
        sax(end).YTick = [];    
        sax(end).XTickLabels = [];    
        Lines(0,[],'w',[],'LineWidth',1);
        

    subplot2(1,2,1,2); 
        mspxz = repmat(mean(spxz(:,:,:,sts,8,1,2),3),[1,2]);
        if sts == 1,
            imagesc(errorBinCenters,...
                    [phzBinCenters;phzBinCenters+2.*pi],...
                    repmat(mean(spxz(:,:,:,size(spxz,4),8,1,2),3),[1,2])');
        end
        hold('on');
        
        clevels = zeros(size(plevels));
        maxIter = 10;        
        for p = 1:numel(plevels),
            step = 0.0001;
            thresh = 0.0001;        
            iter = 1;
            while round(sum(mspxz(mspxz(:)>(thresh))).*100)~=plevels(p).*100 && iter<maxIter,
                while sum(mspxz(mspxz(:)>(thresh)))>plevels(p),
                    thresh = thresh+step;
                end
                step = step./2;
                while sum(mspxz(mspxz(:)>(thresh)))<plevels(p),
                    thresh = thresh-step;            
                end
                step = step./2;
                iter = iter+1;
            end
            clevels(p) = thresh;
        end
        
        if numel(clevels) == 1,
            clevels = repmat(clevels,[1,2]);
        end
        
        [~,ctax] = contour(coords{:},mspxz,clevels);
        ctax.LineWidth = 1;
        ctax.LineColor = cclr(sts);
        
        ylim(circ_ang2rad([-60,420]));
        axis('xy');
        xl = [-200,200];
        xlim(xl);
        % MASK areas outside single cycle
        patch([xl(1),xl(1),xl(2),xl(2)],...
              circ_ang2rad([ 360, 420,420,360]),...
              [0.2,0.2,0.2],...
              'EdgeAlpha',0,...
              'FaceAlpha',0.4);
        patch([xl(1),xl(1),xl(2),xl(2)],...    
              circ_ang2rad([ -60,   0,  0,-60]),...
              [0.2,0.2,0.2],...
              'EdgeAlpha',0,...
              'FaceAlpha',0.4);
        sax(end).XTick = [-100,0,100];
        sax(end).XTickLabels = [-15,0,15];
        sax(end).YTick = [];    
        sax(end).XTickLabels = [];    
        Lines(0,[],'w',[],'LineWidth',1);
    end
% figure
        
% KL divergence


plot(repmat(sq(mean(q25xy(:,:,:,sts),3)),[1,2]),ycoords,'--k');
plot(repmat(sq(mean(q50xy(:,:,:,sts),3)),[1,2]),ycoords,'-k');
plot(repmat(sq(mean(q75xy(:,:,:,sts),3)),[1,2]),ycoords,'--k');

%plot(repmat(dcTpdMax(sts,:),[1,2])', ycoords,'k','LineWidth',2);
%plot(repmat(dcTpdMax(sts,:),[1,2])', ycoords,'k','LineWidth',2);
%plot(repmat(dcTpdWmn(1,:),[1,2])', ycoords,'k','LineWidth',2);    
%plot(repmat(dcTpdWmn(sts,:),[1,2])', ycoords,'m','LineWidth',2);





figure,
for sts = 1:6,
subplot2(6,2,sts,1);
imagesc(errorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,1,sts,8,2),[1,2])');
axis('xy')
Lines(0,[],'k');
subplot2(6,2,sts,2);
imagesc(errorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,2,sts,8,2),[1,2])');
axis('xy')
Lines(0,[],'k');
end
colormap('jet');


figure,
for sts = 1:6,
subplot2(6,2,sts,1);
imagesc(perrorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,1,sts,8,4),[1,2])');
axis('xy')
subplot2(6,2,sts,2);
imagesc(perrorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,2,sts,8,4),[1,2])');
axis('xy')
end

figure,
for sts = 1:6,
subplot2(6,2,sts,1);
imagesc(perrorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,1,sts,8,3),[1,2])');
axis('xy')
subplot2(6,2,sts,2);
imagesc(perrorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,2,sts,8,3),[1,2])');
axis('xy')
end


figure,
for sts = 1:6,
    subplot(6,1,sts);
    ind =   logical(dstcm(:,1))                             ...
            & logical(dstcm(:,sts))                         ...
            & ~any(logical(dstcm(:,[7,8])),2)               ...
            & dpostI                                        ...
            & duincI;        
    out = hist2([duinc(ind),dphz(ind)],1:16,phzBins);
    imagesc(1:16,[phzBinCenters;phzBinCenters+2*pi],repmat(bsxfun(@rdivide,out,sum(out)),[1,2])');
    axis('xy');
end

%%%>>>

%%%<<< COMPUTATION : find ridge values for each state's dcTDPfrontal error as function of theta phz

dcTpdLongSum = {};
dcTpdhpSum = {};
dcTpdbpSum = {};
dcTpdInd = [];
dcTpdWmn = [];
dcTpdMax = [];
for sts = 1:numel(stid),
    dcTpdLongSum{sts} = RectFilter(sum(cat(3,dcTDPfrontal{sts}{:}),3),5,3);
    dcTpdLongSum{sts} = bsxfun(@rdivide,dcTpdLongSum{sts},sum(dcTpdLongSum{sts}));
    [~, dcTpdInd(end+1,:)] = max(dcTpdLongSum{sts});
    dcTpdMax(end+1,:) = errorBinCenters(dcTpdInd(end,:)');    
    dcTpdWmn(end+1,:) = sum(bsxfun(@times,dcTpdLongSum{sts},errorBinCenters'))';
    
    dcTpdhpSum{sts} = RectFilter(sum(cat(3,dcTDPhp{sts}{:}),3),5,3);
    dcTpdhpSum{sts} = bsxfun(@rdivide,dcTpdhpSum{sts},sum(dcTpdhpSum{sts}));
    
    dcTpdbpSum{sts} = RectFilter(sum(cat(3,dcTDPbp{sts}{:}),3),5,3);
    dcTpdbpSum{sts} = bsxfun(@rdivide,dcTpdbpSum{sts},sum(dcTpdhpSum{sts}));
end

%%%>>>

%%%<<< LOAD example data

t = 20;
Trial = Trials{t==tind};
unitSubset = units{t==tind};
exyz = resample(Trial.load('xyz','trb'),30);
%efet = fet_HB_pitchB(Trial,sampleRate);
efet = resample(copy(dc4.fet{tind==t}),sampleRate);
ets = [1:size(exyz,1)]./sampleRate;
estcm = stc2mat(Trial.load('stc'),exyz,states);
eind =   dc4.post{t==tind} > 0.001 ...
       & dc2.uinc{t==tind} > 3     ...
       & estcm(:,1)==1;
enanmask = ones([size(exyz,1),1]);
enanmask(~eind)=nan;

%%%>>>

%%%<<< MAIN FIGURE

% STARTFIG -----------------------------------------------------------------------------------------


eTS = [290,345];
%eTS = [530,600];
%eTS = [1000,1045];
eInds = round(eTS.*sampleRate)+1;
eInds = [eInds(1):eInds(2)]';
eIndsSub = eInds(1300:1600);


[hfig,fig,fax,sax] = set_figure_layout(figure(666007),'A4','portrait',[],1.5,1.5,0,0.2);


%%%<<< PLOT example trajectory within physical space

% ADJUST subplot coordinates
t = 20; 
[yind, yOffSet, xind, xOffSet] = deal(1,-fig.subplot.height-fig.subplot.verticalPadding, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2+ fig.subplot.verticalPadding, ...
                              fig.subplot.height*2+fig.subplot.verticalPadding],...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
imagesc(dc4.pfs{1}.adata.bins{1:2},~maskcirc');
colormap(sax(end),'gray')
caxis(sax(end),[-2,1])

circle(0,0,480,'-k');
plot(exyz(eIndsSub,5,1),...
     exyz(eIndsSub,5,2),...
     '-k','LineWidth',1);
plot(sq(dc4.com{t==tind}(eIndsSub,1)).*enanmask(eIndsSub),...
     sq(dc4.com{t==tind}(eIndsSub,2)).*enanmask(eIndsSub),...
     '-c','LineWidth',1);
xlim([-500,500]);
ylim([-500,500]);
sax(end).XTick = [];
sax(end).YTick = [];
axes(fax);
rectangle('Position',sax(end).Position,'LineWidth',1);

%%%>>>


%%%<<< PLOT example trajectory within behavior space
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(3, -fig.subplot.height-fig.subplot.verticalPadding, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2+fig.subplot.verticalPadding,  ...
                              fig.subplot.height*2+fig.subplot.verticalPadding],...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
imagesc(dc4.pfs{1}.adata.bins{3:4},~maskbhv');
colormap(sax(end),'gray')
caxis(sax(end),[-2,1])

axis('xy');
plot(efet(eIndsSub,1),...
     efet(eIndsSub,2),...
     '-k','LineWidth',1);
plot(sq(dc4.com{t==tind}(eIndsSub,3)).*enanmask(eIndsSub),...
     sq(dc4.com{t==tind}(eIndsSub,4)).*enanmask(eIndsSub),...
     '-c','LineWidth',1);
xlim(sax(end),[-1.8,0.6]);
ylim(sax(end),[-0.6,1.8]);
sax(end).XTick = [];
sax(end).YTick = [];
axes(fax);
rectangle('Position',sax(end).Position,'LineWidth',1);
%%%>>>


%%%<<< PLOT timeseries of xcoord vs decoded xcoord

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(1, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*7,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(eInds),...
     exyz(eInds,5,1),...
     '-k','LineWidth',1);
plot(ets(eInds),...
     sq(dc4.com{t==tind}(eInds,1)).*enanmask(eInds),...
     '-r','LineWidth',1);
plot(ets(eIndsSub),...
     sq(dc4.com{t==tind}(eIndsSub,1)).*enanmask(eIndsSub),...
     '-c','LineWidth',1);
sax(end).XTick = [];
sax(end).YTick = [];
% ADD bounding box        
axes(fax);
rectangle('Position',sax(end).Position,'LineWidth',1);

%%%>>>


%%%<<< PLOT timeseries of ycoord vs decoded ycoord
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*7,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(eInds),...
     exyz(eInds,5,2),...
     '-k','LineWidth',1);
plot(ets(eInds),...
     sq(dc4.com{t==tind}(eInds,2)).*enanmask(eInds),...
     '-r','LineWidth',1);
plot(ets(eIndsSub),...
     sq(dc4.com{t==tind}(eIndsSub,2)).*enanmask(eIndsSub),...
     '-c','LineWidth',1);
sax(end).XTick = [];
sax(end).YTick = [];
% ADD bounding box        
axes(fax);
rectangle('Position',sax(end).Position,'LineWidth',1);
%%%>>>


%%%<<< PLOT timeseries of head pitch vs decoded head pitch
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(3, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*7,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(eInds),...
     efet(eInds,1),...
     '-k','LineWidth',1);
plot(ets(eInds),...
     sq(dc4.com{t==tind}(eInds,3)).*enanmask(eInds),...
     '-r','LineWidth',1);
plot(ets(eIndsSub),...
     sq(dc4.com{t==tind}(eIndsSub,3)).*enanmask(eIndsSub),...
     '-c','LineWidth',1);
sax(end).XTick = [];
sax(end).YTick = [];
ylim([-1.6,0.5]);
% ADD bounding box        
axes(fax);
rectangle('Position',sax(end).Position,'LineWidth',1);
%%%>>>


%%%<<< PLOT timeseries of body pitch vs decoded body pitch

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(4, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*7,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(eInds),...
     efet(eInds,2),...
     '-k','LineWidth',1);
plot(ets(eInds),...
     sq(dc4.com{t==tind}(eInds,4)).*enanmask(eInds),...
     '-r','LineWidth',1);
plot(ets(eIndsSub),...
     sq(dc4.com{t==tind}(eIndsSub,4)).*enanmask(eIndsSub),...
     '-c','LineWidth',1);

sax(end).YTick = [];
ylim([-0.5,1.6]);
% ADD bounding box        
axes(fax);
rectangle('Position',sax(end).Position,'LineWidth',1);

%%%>>>


%%%<<< PLOT spatial vs behaivoral error
% sss

%%%>>>


%%%<<< PLOT BLOCK : theta phase * decoding error
xs = 1;
%%%>>>


%%%<<< PLOT longitudinal error for TDP decoding

% SET constants
xcoords =linspace(-500,500,250);
ycoords = circ_rad2ang([phzBinCenters;phzBinCenters+2*pi]);
for sts = 1:numel(stid),
% ADJUST subplot coordinates    
    [yind, yOffSet, xind, xOffSet] = deal(4+sts, -0.8, xs, 0);

% CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width.*2+0.2,                        ...
                                  fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
% PLOT JPDF of longitudinal error vs theta phase
    imagesc(xcoords,                                 ...
            ycoords,                                 ...
            repmat(mpxy(:,:,sts),[1,2,1])');
          % repmat(RectFilter(dcTpdLongSum{sts},5,3),1,2)')    
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
    plot(repmat(sq(mean(q25xy(:,:,:,sts),3)),[1,2]),ycoords,'--k');
    plot(repmat(sq(mean(q50xy(:,:,:,sts),3)),[1,2]),ycoords,'-k');
    plot(repmat(sq(mean(q75xy(:,:,:,sts),3)),[1,2]),ycoords,'--k');
    
    %plot(repmat(dcTpdMax(sts,:),[1,2])', ycoords,'k','LineWidth',2);
    %plot(repmat(dcTpdMax(sts,:),[1,2])', ycoords,'k','LineWidth',2);
    %plot(repmat(dcTpdWmn(1,:),[1,2])', ycoords,'k','LineWidth',2);    
    %plot(repmat(dcTpdWmn(sts,:),[1,2])', ycoords,'m','LineWidth',2);

    colormap(sax(end),'jet');    
    caxis([0,0.045]);
    ylabel(sax(end),states{sts});
    if sts == 1,
        title(sax(end),{'RCP Error'});
    end
end

% PLOT sptial scale bar 
axes(fax);
line([sax(end).Position(1)+sax(end).Position(3)./4,sax(end).Position(1)+sax(end).Position(3).*3./4],...
     [1,1].*(sax(end).Position([2])-0.2), 'Color', 'k', 'LineWidth', 1);
text(sax(end).Position(1)+sax(end).Position(3)/2,               ...
     sax(end).Position([2])-0.4,                                ...
     '20 cm',                                                   ...
     'HorizontalAlignment','center');

%%%>>>


%%%<<< PLOT lateral error for TDP decoding

% SET constants
xcoords =linspace(-500,500,250);
ycoords = circ_rad2ang([phzBinCenters;phzBinCenters+2*pi]);
for sts = 1:numel(stid),
% ADJUST subplot coordinates    
    [yind, yOffSet, xind, xOffSet] = deal(4+sts, -0.8, xs+2, 0.4);

% CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width.*2+0.2,                        ...
                                  fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
% PLOT JPDF of longitudinal error vs theta phase
    imagesc(xcoords,                                 ...
            ycoords,                                 ...
            repmat(mplt(:,:,sts),[1,2,1])');
          % repmat(RectFilter(dcTpdLongSum{sts},5,3),1,2)')    
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
    plot(repmat(sq(mean(q25lt(:,:,:,sts),3)),[1,2]),ycoords,'--k');
    plot(repmat(sq(mean(q50lt(:,:,:,sts),3)),[1,2]),ycoords,'-k');
    plot(repmat(sq(mean(q75lt(:,:,:,sts),3)),[1,2]),ycoords,'--k');
    
    %plot(repmat(dcTpdMax(sts,:),[1,2])', ycoords,'k','LineWidth',2);
    %plot(repmat(dcTpdMax(sts,:),[1,2])', ycoords,'k','LineWidth',2);
    %plot(repmat(dcTpdWmn(1,:),[1,2])', ycoords,'k','LineWidth',2);    
    %plot(repmat(dcTpdWmn(sts,:),[1,2])', ycoords,'m','LineWidth',2);

    colormap(sax(end),'jet');    
    caxis([0,0.045]);
end

% PLOT sptial scale bar 
axes(fax);
line([sax(end).Position(1)+sax(end).Position(3)./4,sax(end).Position(1)+sax(end).Position(3).*3./4],...
     [1,1].*(sax(end).Position([2])-0.2), 'Color', 'k', 'LineWidth', 1);
text(sax(end).Position(1)+sax(end).Position(3)/2,               ...
     sax(end).Position([2])-0.4,                                ...
     '20 cm',                                                   ...
     'HorizontalAlignment','center');

%%%>>>


%%%<<< PLOT head pitch error

% SET constants
xcoords =linspace([-1.25,1.25,50]);
ycoords = circ_rad2ang([phzBinCenters;phzBinCenters+2*pi]);
% ADJUST subplot coordinates
for sts = 1:numel(stid),
% ADJUST subplot coordinates    
    [yind, yOffSet, xind, xOffSet] = deal(4+sts, -0.8, xs+4, 0.8);
% CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width.*2+0.2,                 ...
                                  fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
% PLOT JPDF of longitudinal error vs theta phase
    imagesc(xcoords,                                 ...
            ycoords,       ...
            repmat(mphp(:,:,sts),[1,2,1])')
          % repmat(RectFilter(dcTpdhpSum{sts},5,3),1,2)')    
    axis('xy');
    xl = [-1.25,1.25];
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
    plot(repmat(sq(mean(q25hp(:,:,:,sts),3)),[1,2]),ycoords,'--k');
    plot(repmat(sq(mean(q50hp(:,:,:,sts),3)),[1,2]),ycoords,'-k');
    plot(repmat(sq(mean(q75hp(:,:,:,sts),3)),[1,2]),ycoords,'--k');
    
    colormap(sax(end),'jet');    
    caxis(sax(end),[0,0.075]);
    if sts == 1,
        %ylabel(sax(end),{'Head','Pitch'});
        title(sax(end),{'HP Error'});        
    end
    
end

% PLOT scale bar 
axes(fax);
line([sax(end).Position(1)+sax(end).Position(3).*1.5./5,        ...
      sax(end).Position(1)+sax(end).Position(3).*3.5./5],       ...
      [1,1].*(sax(end).Position([2])-0.2),                      ...
      'Color', 'k',                                             ...
      'LineWidth', 1);
text(sax(end).Position(1)+sax(end).Position(3)/2,               ...
     sax(end).Position([2])-0.5,                                ...
     '1 rad',                                                   ...
     'HorizontalAlignment','center');

%%%>>>


%%%<<< PLOT body pitch error

% SET constants

xcoords =linspace([pbound,50]);
ycoords = circ_rad2ang([phzBinCenters;phzBinCenters+2*pi]);
for sts = 1:numel(stid),
% ADJUST subplot coordinates    
    [yind, yOffSet, xind, xOffSet] = deal(4+sts, -0.8, xs+6, 1.2);
% CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width.*2+0.2,                        ...
                                  fig.subplot.height],                      ...
                      'FontSize', 8,                                        ...
                      'LineWidth',1);
    hold(sax(end),'on');
% PLOT JPDF of longitudinal error vs theta phase
    imagesc(xcoords,                                 ...
            ycoords,       ...
            repmat(mpbp(:,:,sts),[1,2,1])');
           %repmat(RectFilter(dcTpdbpSum{sts},5,3),1,2)')
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
          'Edgealpha',0,...
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
    plot(repmat(sq(mean(q25bp(:,:,:,sts),3)),[1,2]),ycoords,'--k');
    plot(repmat(sq(mean(q50bp(:,:,:,sts),3)),[1,2]),ycoords,'-k');
    plot(repmat(sq(mean(q75bp(:,:,:,sts),3)),[1,2]),ycoords,'--k');
    
    colormap(sax(end),'jet');    
    caxis([0,0.1]);
    if sts == 1,
        %ylabel(sax(end),{'Body','Pitch'});
        title(sax(end),{'BP Error'});        
% ADJUST subplot coordinates        
        [yind, yOffSet, xind, xOffSet] = deal(4+sts, -0.8, xs+8, 1.4);
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
    else
% ADJUST subplot coordinates        
        [yind, yOffSet, xind, xOffSet] = deal(4+sts, -0.8, xs+8, 1.4);
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
    end
    
end
% PLOT scale bar 
axes(fax);
line([sax(end-1).Position(1)+sax(end-1).Position(3).*1.5./5,        ...
      sax(end-1).Position(1)+sax(end-1).Position(3).*3.5./5],       ...
      [1,1].*(sax(end).Position(2)-0.2),                      ...
      'Color', 'k',                                             ...
      'LineWidth', 1);
text(sax(end-1).Position(1)+sax(end-1).Position(3)/2,               ...
     sax(end-1).Position(2)-0.5,                                ...
     '1 rad',                                                   ...
     'HorizontalAlignment','center');

%%%>>>

% STOPFIG ------------------------------------------------------------------------------------------

%%%>>>



