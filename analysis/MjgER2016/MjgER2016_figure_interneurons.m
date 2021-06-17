% MjgER2016_figure_interneurons




%%%<<< LOAD DATA -----------------------------------------------------------------------------------


configure_default_args();

global AP
% compute_bhv_ratemaps_shuffled --------------------------------------------------------------------
AP.compute_bhv_ratemaps_shuffled =                                                               ...
    struct('get_featureSet',            @fet_HB_pitchB,                                          ...
           'sampleRate',                16,                                                      ...
           'pfsArgs',                   struct('states',           activeState,                  ...
                                               'binDims',          [0.1,0.1],                    ...
                                               'SmoothingWeights', [1.8,1.8],                    ...
                                               'numIter',          1001,                         ...
                                               'posShuffle',       true,                         ...
                                               'boundaryLimits',   [-2,0.8;-0.8,2],              ...
                                               'halfsample',       true),                        ...
           'threshRate',                1,                                                       ...
           'threshDist',                1000                                                     ...
           );
%---------------------------------------------------------------------------------------------------
n
MjgER2016_load_data();

pargs = struct('numIter'   , 1000,...
               'halfsample', true,...
               'posShuffle', true);

bfs   = cf(@(t,u)  compute_bhv_ratemaps(t,u,'overwrite',false),          Trials, units);
pftin = cf(@(t,u)  pfs_2d_theta(t,u),                                    Trials, unitsInts);
pftis = cf(@(t,u)  pfs_2d_theta(t,u,'pfsArgsOverride',pargs),            Trials, unitsInts);
pfsin = cf(@(t,u)  pfs_2d_states(t,u),                                   Trials, unitsInts);
pfsis = cf(@(t,u)  pfs_2d_states(t,u, 'pfsArgsOverride',pargs),          Trials, unitsInts);
bfsin = cf(@(t,u)  compute_bhv_ratemaps(t,u,'overwrite',false),          Trials, unitsInts);
bfsis = cf(@(t,u)  compute_bhv_ratemaps_shuffled(t,u,'overwrite',false), Trials, unitsInts);


cf(@(t) t.load('nq'), Trials);
ampSymInt = cf(@(t,u) t.nq.AmpSym(u), Trials,unitsInts);
ampSymInt = cat(1,ampSymInt{:});
ampSymPyr = cf(@(t,u) t.nq.AmpSym(u), Trials,units);
ampSymPyr = cat(1,ampSymPyr{:});

spkWidthRInt = cf(@(t,u) t.nq.SpkWidthR(u), Trials,unitsInts);
spkWidthRInt = cat(1,spkWidthRInt{:});
spkWidthRPyr = cf(@(t,u) t.nq.SpkWidthR(u), Trials,units);
spkWidthRPyr = cat(1,spkWidthRPyr{:});

centerMaxInt = cf(@(t,u) t.nq.CenterMax(u), Trials,unitsInts);
centerMaxInt = cat(1,centerMaxInt{:});
centerMaxPyr = cf(@(t,u) t.nq.CenterMax(u), Trials,units);
centerMaxPyr = cat(1,centerMaxPyr{:});

figure,
hold('on');
scatter(spkWidthRPyr,ampSymPyr,10,log10(abs(centerMaxPyr)),'filled');
scatter(spkWidthRInt,ampSymInt,10,log10(abs(centerMaxInt)),'filled');
colormap('jet');
colorbar();


figure,
hold('on');
plot(spkWidthRPyr,ampSymPyr,'.');
plot(spkWidthRInt,ampSymInt,'.r');


% reorder state fields and concat with theta state
pfs = cf(@(t,s) cat(2,{t},{s{[4,3,7,2,6]}}), pftin, pfsin);


[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfs, units, [], [], false);

overwrite = false;
% SPACE X BHV
tag = 'interneurons_xyhb_2020';
tag = 'interneurons_xyhb_2020_final';    
tag = 'interneurons_xyhb_2020_loc';      

tag = 'interneurons_xyhb_2020_theta';    
pfi = cf(@(t,u) req20201117(t,u,tag,false,false,overwrite), Trials,unitsInts);
tag = 'interneurons_xyhb_2020_loc_rear'; 
pfiLR = cf(@(t,u) req20201117(t,u,tag,false,false,overwrite), Trials,unitsInts);
tag = 'interneurons_xyhb_2020_pause_rear'; 
pfiPR = cf(@(t,u) req20201117(t,u,tag,false,false,overwrite), Trials,unitsInts);

%%%>>>



%%%<<< COMPUTE SI from behavior ratemaps -----------------------------------------------------------
[rmapNB,cmap] = decapsulate_and_concatenate_mtaapfs(bfsin,unitsInts);
[rmapSB] = decapsulate_and_concatenate_mtaapfs(bfsis,unitsInts);
% MASK 
maskBhv = validDims;

rmapNB(~maskBhv,:) = nan;
rmapSB(~maskBhv,:,:) = nan;

% sinfNB -> spatial information of normal behavior field
sinfNB = sum(bsxfun(@times,                                                                 ...
                    1./sum(double(bsxfun(@times,validDims,~isnan(rmapNB))),'omitnan')       ...
                 ,     bsxfun(@rdivide,rmapNB,mean(rmapNB,'omitnan'))                       ... 
                 .*log2(bsxfun(@rdivide,rmapNB,mean(rmapNB,'omitnan')))),                   ...
             'omitnan')';

% sinfSB -> spatial information of shuffled behavior field
sinfSB = sq(sum(bsxfun(@times,                                                              ...
                       1./sum(double(bsxfun(@times,validDims,~isnan(rmapSB))),'omitnan')    ...
                       ,        bsxfun(@rdivide,rmapSB,mean(rmapSB,'omitnan'))              ... 
                         .*log2(bsxfun(@rdivide,rmapSB,mean(rmapSB,'omitnan')))),'omitnan'));
%%%>>>



%%%<<< COMPUTE SI from place ratemaps --------------------------------------------------------------
% rmapNP -> rate map normal place field
rmapNP = decapsulate_and_concatenate_mtaapfs(pftin,unitsInts);
rmapSP = decapsulate_and_concatenate_mtaapfs(pftis,unitsInts);
% MASK 
width = pftin{1}.adata.binSizes(1);
height = pftin{1}.adata.binSizes(2);
radius = round(width/2)-find(pftin{1}.adata.bins{1}<=-450,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
maskPlace = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);

rmapNP(~maskPlace(:),:,:) = nan;
rmapSP(~maskPlace(:),:,:) = nan;

% sinfNP -> spatial information of normal place field
sinfNP = sum(bsxfun(@times,                                                                 ...
                    1./sum(double(bsxfun(@times,maskPlace(:),~isnan(rmapNP))),'omitnan')    ...
                    ,bsxfun(@rdivide,rmapNP,mean(rmapNP,'omitnan'))                         ... 
                     .*log2(bsxfun(@rdivide,rmapNP,mean(rmapNP,'omitnan')))),               ...
             'omitnan')';

% sinfSP -> spatial information of shuffled place field
sinfSP = sq(sum(bsxfun(@times,...
                       1./sum(double(bsxfun(@times,maskPlace(:),~isnan(rmapSP))),'omitnan') ...
                       ,bsxfun(@rdivide,rmapSP,mean(rmapSP,'omitnan'))                      ... 
                .*log2(bsxfun(@rdivide,rmapSP,mean(rmapSP,'omitnan')))),'omitnan'));
%%%>>>



%%%<<< COMPUTE mean rate from ratemaps -------------------------------------------------------------

% TODO Interneuron bhv score ratio (rear/walk) vs theta phase preference
rmapNT = decapsulate_and_concatenate_mtaapfs(cf(@(p) p{1}, pfs),unitsInts);
rmapNR = decapsulate_and_concatenate_mtaapfs(cf(@(p) p{2}, pfs),unitsInts);
rmapNH = decapsulate_and_concatenate_mtaapfs(cf(@(p) p{3}, pfs),unitsInts);
rmapNL = decapsulate_and_concatenate_mtaapfs(cf(@(p) p{5}, pfs),unitsInts);


mrateNR = mean(rmapNR,'omitnan')';
mrateNH = mean(rmapNH,'omitnan')';
mrateNL = mean(rmapNL,'omitnan')';
mrateNB = mean(rmapNB,'omitnan')';
mrateNT = mean(rmapNT,'omitnan')';

%%%>>>



%%%<<< COMPUTE zscore of SI ------------------------------------------------------------------------

% ZSCORE 
zsinfNP = (sinfNP-mean(sinfSP,2))./std(sinfSP,[],2);
zsinfNB = (sinfNB-mean(sinfSB,2))./std(sinfSB,[],2);

% $$$ 
% $$$ figure();
% $$$ subplot(121);
% $$$ hold('on');
% $$$ plot(zsinfNP,zsinfNB,'.');
% $$$ line([-20,100],[-20,100]);
% $$$ circle(0,0,3.1)
% $$$ xlim([-10,100]);
% $$$ ylim([-10,100]);
% $$$ subplot(122);
% $$$ plot(log10(sinfNP),log10(sinfNB),'.');
% $$$ xlim([-4,0]);
% $$$ ylim([-4,0]);
% $$$ line([-4,0],[-4,0]);


%%%>>>



%%%<<< DECOMPOSE factors of BhvPos ratemaps --------------------------------------------------------

%%%<<< DECAPSULATE rate mapsSPACE X BHV

[rmaps,cluSessionMap] = decapsulate_and_concatenate_mtaapfs(pfi,unitsInts);
[rmapsLR,~] = decapsulate_and_concatenate_mtaapfs(pfiLR,unitsInts);
[rmapsPR,~] = decapsulate_and_concatenate_mtaapfs(pfiPR,unitsInts);

bhvMask = false(size(validDims));
bhvMask(validDims) = true;
bhvMask = reshape_eigen_vector(bhvMask,bfs)';
bhvLims = [-1.6, 0.6; ...
           -0.5, 1.7];

bins = pfi{1}.adata.bins;
binDims = pfi{1}.adata.binSizes';

% RESHAPE rmaps into [ xPosition, yPosition, headPitch, bodyPitch, unit ]
rmapa = nan([binDims,size(rmaps,2)]);
rmapaLR = nan([binDims,size(rmaps,2)]);
rmapaPR = nan([binDims,size(rmaps,2)]);
for u = 1:size(rmaps,2),rmapa(:,:,:,:,u) = reshape(rmaps(:,u),binDims); end
for u = 1:size(rmapsLR,2),rmapaLR(:,:,:,:,u) = reshape(rmapsLR(:,u),binDims); end
for u = 1:size(rmapsPR,2),rmapaPR(:,:,:,:,u) = reshape(rmapsPR(:,u),binDims); end

%%%>>>
[bmaps,cluSessionMap] = decapsulate_and_concatenate_mtaapfs(bfsin,unitsInts);

%%%<<< DIAGNOSTIC PLOT - bhv rate map
% PLOT example
% $$$ mind = find(ismember(cluSessionMap,[20,106],'rows'))
cluSessionMap(mind,:)
rmapLR = rmapaLR(:,:,:,:,mind);
rmapPR = rmapaPR(:,:,:,:,mind);
hfig = figure();
for x = 2:5,
    for y = 2:5,
        subplot2(6,6,7-y,x);
        pax = pcolor(bins{3:4},sq(rmapLR(x,y,:,:)-rmapPR(x,y,:,:))');
        pax.EdgeColor = 'none';
        axis('xy');
        caxis([-10,10]);
        colormap('jet');
    end
end
cax = colorbar();
cax.Position(1) = 0.01;
drawnow();
pause(0.1);
cax.Position(1) = sum(hfig.CurrentAxes.Position([1,3]))+cax.Position(1);
%%%>>>

% bmaps = reshape(rmapa(3,4,:,:,:),[],size(rmapa,length(binDims)+1));



%%%<<< COLLECT inner spatial bins -> bmaps
% $$$ bmaps = [];
% $$$ for k = 2:5,
% $$$     for j = 2:5,
% $$$         bmaps(:,:,(k-2)*4+(j-1)) = reshape(rmapa(k,j,:,:,:),[],size(rmapa,length(binDims)+1));
% $$$     end
% $$$ end
% $$$ bmaps(:,:,[1,4,9,16]) =[];
% $$$ bmaps = reshape(bmaps,[size(bmaps,1),prod(size(bmaps))/size(bmaps,1)]);
%%%>>>

% COMPRESS bmaps to only valid behavior space elements
vmaps = bmaps(validDims,:);
%vmaps(:,sum(isnan(vmaps))>0)= 0;
vmaps(isnan(vmaps(:)))=0;

% COMPUTE factors
[LU,LR,FSr,VT] = erpPCA(vmaps',5);

% PLOT first 5 eigen vectors
figure()
for v = 1:5,
    eigVec = nan([binDims(end-1:end)]);
    eigVec(validDims) = LU(:,v);
    subplot(1,6,v);
    pax = pcolor(bins{end-1:end},eigVec');
    pax.EdgeColor = 'none';
    axis('xy');
    colorbar();
end
subplot(1,6,6);
plot(VT(1:5,4),'-+')

% $$$ figure();
% $$$ scatter(FSr(:,1),FSr(:,2),10,FSr(:,3),'Filled');
% $$$ 
% $$$ figure();
% $$$ scatter(FSr(:,1),FSr(:,2),10,mrateNT','Filled');

%%%>>>



%%%<<< COMPUTE intra units bhv ratemap rank correlation --------------------------------------------

% COLLECT inner spatial bins
imaps = [];
for k = 2:5,
    for j = 2:5,
        imaps(:,:,k-1,j-1) = bsxfun(@times,reshape(rmapa(k,j,:,:,:),[], size(rmapa,length(binDims)+1)),maskBhv(:));
    end
end

nimap = [];
dcoor = [];
d = 1;
for k = 1:4,
    for j = 1:4,
        if (k==1&&j==1) || (k==1&&j==4)|| (k==4&&j==4)|| (k==4&&j==1),
            continue
        end
        nimap(:,:,d) = imaps(:,:,k,j);
        dcoor(d,:) = [j,k];
        d=d+1;
    end
end

rho = [];
for u = 1:size(imaps,2),
    d = 1;    
    for j = 1:11
        for k= j+1:12
            x1 = nimap(:,u,k);
            x2 = nimap(:,u,j);
            ind = nniz(x1) & nniz(x2) & validDims;
            rho(u,d) = corr(nimap(ind,u,j),nimap(ind,u,k),'type','Spearman');
            d = d+1;
        end
    end
end

% $$$ rhoRank = [];
% $$$ for u = 1:size(imaps,2),
% $$$     d = 1;    
% $$$     for j = 1:11
% $$$         for k= j+1:12
% $$$             x1 = nimap(:,u,k);
% $$$             x2 = nimap(:,u,j);
% $$$             ind = nniz(x1) & nniz(x2) & validDims;
% $$$             [~,rank1] = sort(nimap(ind,u,j));
% $$$             [~,rank2] = sort(nimap(ind,u,k));
% $$$             rhoRank(u,d) = corr(rank1,rank2,'type','Pearson');
% $$$             d = d+1;
% $$$         end
% $$$     end
% $$$ end
% $$$ [~,rank1] = sort(nimap(ind,u,j));
% $$$ [~,rank2] = sort(nimap(ind,u,k));
% $$$ figure,plot(rank1,rank2,'.');
% $$$ figure,plot(nimap(ind,u,j),nimap(ind,u,k),'.');


dst = [];
d = 1;    
for j = 1:11
    for k= j+1:12
        dst(d) = sqrt(sum((dcoor(j,:)-dcoor(k,:)).^2));
        d = d+1;
    end
end

%rhoDR =[];
%rhoDRMean = [];
rhoDRPoly = [];
for u = 1:size(imaps,2),
    %rhoDR(u) = corr(dst',rho(u,:)');
    %rhoDRMean(u) = mean(rho(u,:));
    rhoDRPoly(u,:) =polyfit(dst',rho(u,:)',1);
end

dsti = discretize(dst,[0,1.5,2.5,3.5]);
thr = [0:0.1:0.9];
rhoDRratio = [];
for u = 1:225,
    for b = 1:3,
        for t =1:numel(thr)
        rhoDRratio(u,b,t) = sum(rho(u,dsti==b)>thr(t))./sum(rho(u,dsti==b)>-1);
    end
end

    
        

figure
u = 1
plot(dst(:),rho(165,:),'.')

% $$$ 
% $$$ rho = [];
% $$$ dcor = [];
% $$$ for u = 1:size(rmapa,5),
% $$$     d = 1;
% $$$     for k = 1:4,
% $$$         for j = 1:4,
% $$$             if k==j || (k==1&&j==1) || (k==1&&j==4)|| (k==4&&j==4)|| (k==4&&j==1),
% $$$                 continue
% $$$             end
% $$$             x1 = imaps(:,u,k,j);
% $$$             x2 = imaps(:,u,j,k);
% $$$             ind = nniz(x1) & nniz(x2) & validDims;
% $$$             rho(u,d) = corr(imaps(ind,u,k,j), imaps(ind,u,j,k), 'type', 'Spearman');
% $$$             %nimap(:,:,d) = 
% $$$             d = d+1;
% $$$         end
% $$$     end
% $$$ end

%%%>>>



%%%<<< PLOT examples of bhv rate correlations ------------------------------------------------------

% INTERACTIVE figure of inter map rank correlations ------------------------------------------------
sax = gobjects([1,5]);
hfig = figure()
sax(1) = subplot2(2,3,1,1);
hold('on');
%    pax = scatter(rhoDR, zsinfNP, 10, mrateNT, 'Filled');
%    pax = scatter(rhoDR, zsinfNP, 10, tppRUN', 'Filled');
    pax = scatter( zsinfNP',zsinfNB', 20, rhoDR', 'Filled');    
%pax = scatter(rhoDR, zsinfNP, 10, tppRUN', 'Filled');    
%    plot(rho(mind,:),dst(mind,:),'.');
%$plot(rho(mind,:),dst(mind,:),'.');
%    pax = scatter(mean(rho,2), std(rho,[],2), 10, mrateNT, 'Filled');    
    circAx = plot(0,0,'or','MarkerSize',20);
    %colormap(gca(),'hsv');
    colormap(gca(),'jet');    
    colorbar();
    %caxis(gca(),[0,10]);    
    caxis(gca(),[-1,0]);        
    %caxis(gca(),[-pi,pi]);    
    %Lines([],3.1,'k');
    %xlabel('Mean Spearman Corr');
    ylabel('BHV Information z-score');    
    xlabel('POS Information z-score');
sax(5) = subplot2(2,3,2,1);
hold('on');
    pax2 = scatter(rhoDR,std(rho,[],2),10,tppRUN','Filled');
    circAx2 = plot(0,0,'or','MarkerSize',20);
    colormap(gca(),'hsv');
    caxis(gca(),[-pi,pi]);    
    
while true
    [~,mind] = min(sqrt(sum(bsxfun(@minus,[pax.XData',pax.YData'],sax(1).CurrentPoint(1,[1,2])).^2,2)));    
    rmap = rmapa(:,:,:,:,mind);
    rmapEX = [];
    for x = 2:5,
        for y = 2:5,
            if (x==2&&y==2) || (x==2&&y==5)|| (x==5&&y==5)|| (x==5&&y==2)
                continue
            end
            rmapEX((28*(x-2)+1):(28*(x-2)+28),(28*(y-2)+1):(28*(y-2)+28)) = sq(rmap(x,y,:,:))'.*double(maskBhv);
        end
    end
    clim = [min(nonzeros(rmapEX(:))),max(nonzeros(rmapEX(:)))];
    
delete(sax(2));
sax(2) = subplot2(2,3,1,2);
    trlunt = cluSessionMap(mind,:);
    plot(pftin{trlunt(1)},trlunt(2),1,'colorbar',clim,'colorMap',@jet);
    grid('on');    
    sax(2).XTick = [-400,-200,0,200,400];
    sax(2).YTick = [-400,-200,0,200,400];    
    sax(2).GridColor = 'm';
    circAx.XData = pax.XData(mind);
    circAx.YData = pax.YData(mind);
    
delete(sax(3));
sax(3) = subplot2(2,3,1,3);
    imagesc(rmapEX);
    colormap(gca(),'jet');
    colorbar();
    caxis(clim);
    axis('xy');
    title(num2str(trlunt,'Trial:%i Unit:%i'));
sax(3) = subplot2(2,3,2,3);    
    plot(dst(:),rho(mind,:),'.')
    ylim([-1,1]);

    circAx2.XData = pax2.XData(mind);
    circAx2.YData = pax2.YData(mind);
    
drawnow()
waitforbuttonpress();
end    

%%%>>> 



%%%<<< COMPUTE interneuron theta phase preference --------------------------------------------------

% COMPUTE theta phase preference for run and rem states
nUnits = sum(cellfun(@numel,unitsInts));
nBins = 36;
% $$$ tphRUN = nan([nUnits,nBins]);
% $$$ tppRUN = nan([nUnits,1]);
% $$$ tprRUN = nan([nUnits,1]);
tptRUN = nan([nUnits,1]);

for t = 1:numel(Trials)
    Trial = Trials{t};    
    disp(Trial.name);
% LOAD theta phase
    phz = load_theta_phase(Trial,[],sessionList(t).thetaRefGeneral,phzCorrection(t));
% LOAD spike groups
    spk = Trial.load('spk', phz.sampleRate, 'theta-groom-sit', unitsInts{t});
% COLLECT unit theta phase stats
    for u = find(cluSessionMap(:,1) == t)',
        res = spk(cluSessionMap(u,2));               
        if numel(res)>3
            tphz = phz(res);
% $$$             tphRUN(u,:) = histcounts(tphz,linspace(0,2*pi,nBins+1));
% $$$             tppRUN(u,1) = circ_mean (tphz);
% $$$             tprRUN(u,1) = circ_r    (tphz);
            tptRUN(u,1) = circ_rtest(tphz);
        end%if
    end%for u
end%for t
%%%>>>



%%%<<< COMPUTE speed vs rate corr ------------------------------------------------------------------

%vrCorr = nan([nUnits,1]);
vrCorrLog = nan([nUnits,1]);
for t = 1:numel(Trials)
    Trial = Trials{t};    
    disp(Trial.name);
    xyz = preproc_xyz(Trial,'trb');
    vxy = vel(filter(copy(xyz),'ButFilter',4,1.5),'hcom',[1,2]);
    vxy.data(vxy.data<0.01) = 0.01;
% LOAD spike groups
    spk = Trial.load('spk', xyz.sampleRate, 'theta-groom-sit', unitsInts{t});
    ufr = Trial.load('ufr', xyz,spk,unitsInts{t},0.3);
% COLLECT unit theta phase stats
    uset = find(cluSessionMap(:,1) == t);
    for u = uset',
        res = spk(cluSessionMap(u,2));
        res(res>size(ufr,1)) = [];
        if numel(res)>3
            vrCorrLog(u,1) = corr(log10(vxy(res,1)),ufr(res,uset==u));            
            %vrCorr(u,1) = corr(vxy(res,1),ufr(res,uset==u));
        end%if
    end%for u
end%for t

%%%>>>



%%%<<< DECOMPOSE (erpPCA) interneuron autocorrelogram ----------------------------------------------

halfBins = 100;
accgScale = nan([nUnits,2*halfBins+1]);
accgHz = nan([nUnits,2*halfBins+1]);
for t = 1:numel(Trials)
    Trial = Trials{t};    
    disp(Trial.name);
% LOAD spike groups
    spk = Trial.load('spk', Trial.lfp.sampleRate, 'theta-groom-sit', unitsInts{t});
    for u = find(cluSessionMap(:,1) == t)',    
        res = spk(cluSessionMap(u,2));
        [accgHz(u,:),tbin] = CCG(res,u,1,halfBins,Trial.lfp.sampleRate,u,'hz');
        %[accgScale(u,:),tbin] = CCG(res,u,1,halfBins,Trial.lfp.sampleRate,u,'scale');        
    end
end


%[LUa,LRa,FSra,VTa] = erpPCA(accgScale(:,:));
[LUa,LRa,FSra,VTa] = erpPCA(RectFilter(accgHz',3,3)');

figure,
imagesc(RectFilter(accgHz(:,:)',3,3)')

figure,
for v = 1:10,
    subplot(11,1,v);
    plot(LRa(:,v));
end
subplot(11,1,v+1);
plot(VTa(1:10,4),'+-');


figure
for v = 1:5,
    subplot(2,6,v);
    plot(LRa(:,v));
    title(['Factor: ',num2str(v)]);
    subplot(2,6,v+6);
    plot(FSra(:,v),tprRUN,'.');
    xlabel('score');
    ylabel('tpr');
end
subplot(2,6,v+1);
plot(VTa(1:5,4),'+-');
xlabel('Factor');
ylabel('exp var');

figure();
ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]);
scatter(FSr(ind,1),FSr(ind,2),20,vrCorr(ind),'filled');
colormap('jet');
colorbar();

f = 2;
figure
subplot(221);
    ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN>0.15;
    scatter(tppRUN(ind)+2*pi*double(tppRUN(ind)<0),FSra(ind,f),20,tprRUN(ind),'Filled');
    xlim([0,2*pi]);
    ylim([-3,3]);
    grid('on');
subplot(223);
    ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN<=0.15;
    scatter(tppRUN(ind)+2*pi*double(tppRUN(ind)<0),FSra(ind,f),20,tprRUN(ind),'Filled');
    xlim([0,2*pi]);
    ylim([-3,3]);    
    grid('on');
subplot(222);
    ind = ~ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN>0.15;
    scatter(tppRUN(ind)+2*pi*double(tppRUN(ind)<0),FSra(ind,f),20,tprRUN(ind),'Filled');
    xlim([0,2*pi]);
    ylim([-3,3]);    
    grid('on');
subplot(224);
    ind = ~ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN<=0.15;
    scatter(tppRUN(ind)+2*pi*double(tppRUN(ind)<0),FSra(ind,f),20,tprRUN(ind),'Filled');
    xlim([0,2*pi]);
    ylim([-3,3]);        
    grid('on');

%%%>>>



%%%<<< MAIN FIGURE portrait ------------------------------------------------------------------------

% bhvState spatial ratemap examples (placefields)
% place-bhehavior ratemaps          (space binned behaviorRatemaps)
% spatial and behavioral variation 
% Factors and explained variance

%%%<<< setup

[hfig,fig,fax,sax] = set_figure_layout(figure(666001),'A4','portrait',[],1.5,1.5,0.1,0.2);
expUnitsPP = {[ 3],[124];...
              [18],[ 59];...
              [4],[ 60];...
              [22],[ 49];...
              [25],[ 7]};

eupInd = [2,1;...
          3,1;...
          5,1];
eclr = 'rmg';
mrk  = '^so';

numUnits = sum(cellfun(@numel,expUnitsPP(:,2)));
nanColor = [0.25,0.25,0.25];
labels = {{'Theta','Behavior'},{'Theta','Arena'},'Rear','H Loc','H Pause','L Loc','L Pause'};
maskBhv = reshape(validDims,bfsin{1}.adata.binSizes');
axOpts = {'Units',                 'centimeters',...
          'FontSize',              8,            ...
          'LineWidth',             1,            ...
          'PlotBoxAspectRatioMode','manual'};

%%%>>>


%%%<<< Panel - Example Interneuorns by Behavior State

[yind, yOffSet, xind, xOffSet] = deal(1, 0, 1, 0);
for tind = 1:size(expUnitsPP,1),
    t = expUnitsPP{tind,1}(1);
    for uind = 1:numel(expUnitsPP{tind,2}),
        unit = expUnitsPP{tind,2}(uind);
        
        % GET maximum firing rate of all states
        maxPfsRate = max(cell2mat(cf(@(p,u,m)                                   ...
                                  p.maxRate(u,true&numel(m)==1,1,[],[],[],m),   ...
                                  [pfs{t},{bfsin{t}}],                          ...
                                  [repmat({unit},[1,numel(pfs{t})+1])],         ...
                                  [repmat({true},[1,numel(pfs{t})]),{maskBhv(:)}])));
        minPfsRate = min(cell2mat(cf(@(p,u,m)                                   ...
                                  minRate(p,u,true&numel(m)==1,1,[],[],[],m),   ...
                                  [pfs{t}],                          ...
                                  [repmat({unit},[1,numel(pfs{t})])],         ...
                                  [repmat({true},[1,numel(pfs{t})])])));

        
% BHV fields MTAApfs
        [yind, yOffSet, xind, xOffSet] = deal(yind, yOffSet, 1, 0);
        sax(end+1) = axes(axOpts{:},                                            ...
                          'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                      fig.page.ypos(yind)+yOffSet,              ...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height]);
        hold(sax(end),'on');
        plot(bfsin{t},unit,1,false,[minPfsRate,maxPfsRate],false,0.5,false,[],@jet,maskBhv,nanColor);            
        %plot(pfs{t}{s},unit,1,false,[0,maxPfsRate],true,0.5,false,interpParPfs,@jet,[],nanColor);
        sax(end).YTickLabel = {};
        sax(end).XTickLabel = {};
        box(sax(end),'on');
        
        xlim([-1.7,0.5]);
        ylim([-0.5,1.7]);
        
        exampleUnitInd = find(yind==eupInd(:,1));
        if ~isempty(exampleUnitInd),
            scatter(0.25,1.5,30,eclr(exampleUnitInd),mrk(exampleUnitInd),'filled'),
        end

% ADD Behavior space scale bars                
        if yind == 1
            title(labels{xind},'FontSize',8,'FontWeight','normal');
            line([-1.65,-1.2],[1.65,1.65],'Color','w','LineWidth',1);
            line([-1.65,-1.65],[1.2,1.65],'Color','w','LineWidth',1);        
        end
        
            
        end 
        if yind == 3
            ylabel(sax(end),'Body Pitch');
        end
        if yind == numUnits
            xlabel(sax(end),'Head Pitch');
        end
        
        
        for s = 1:numel(pfs{t}),
% PLACEFIELDS MTAApfs
            [yind, yOffSet, xind, xOffSet] = deal(yind, yOffSet, s+1, 0.5+0.25*double(s~=1));
            sax(end+1) = axes(axOpts{:},                                            ...
                              'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                          fig.page.ypos(yind)+yOffSet,              ...
                                          fig.subplot.width,                        ...
                                          fig.subplot.height]);

            plot(pfs{t}{s},unit,1,false,[minPfsRate,maxPfsRate],true,0.5,false,[],@jet,[],nanColor);            
            %plot(pfs{t}{s},unit,1,false,[0,maxPfsRate],true,0.5,false,interpParPfs,@jet,[],nanColor);
            sax(end).YTickLabel = {};
            sax(end).XTickLabel = {};
            
            if yind == 1, 
                title(labels{xind});
            end 
% ADD Arena scale bars        
            if yind == 1 && s ==1,
                line([-490,-250],[490,490],'Color','w','LineWidth',3);
                line([-490,-490],[250,490],'Color','w','LineWidth',3);        
            end
            
            

            if s==6,
                text(sax(end),...
                     492,                                  ... xpos
                     -392,                                  ... ypos
                     num2str(round(maxPfsRate)),           ... rate
                                'FontSize', 8,             ...
                                   'Color', [1,1,1],       ...
                     'HorizontalAlignment', 'right');
                text(sax(end),...
                     -472,                                 ... xpos
                     -392,                                  ... ypos
                     num2str(round(minPfsRate)),           ... rate
                                'FontSize', 8,             ...
                                   'Color', [1,1,1],       ...
                     'HorizontalAlignment', 'left');
            end

            % ADD arena scale bar
% $$$             if  yind == numUnits && s == 1,            
% $$$                 line(fax,...
% $$$                      [sax(end).Position(1),sax(end).Position(1)+sax(end).Position(3)./5],...
% $$$                      [sax(end).Position(2),sax(end).Position(2)]-0.1,...
% $$$                      'LineWidth',2,...
% $$$                      'Color','k');
% $$$                 %line(fax,...
% $$$                 %     [sax(end).Position(1),sax(end).Position(1)]-0.1,...
% $$$                 %     [sax(end).Position(2),sax(end).Position(2)+sax(end).Position(4)./5],...
% $$$                 %     'LineWidth',1,...
% $$$                 %     'Color','k');
% $$$                 text(fax,...
% $$$                      sax(end).Position(1),...
% $$$                      sax(end).Position(2)-0.3,...
% $$$                      '20 cm',...
% $$$                      'HorizontalAlignment','left');
% $$$             end
            
            % ADD colorbar
            if  yind == numUnits && s == numStates,
                cax = colorbar(sax(end),'southoutside');
                colormap(cax,'jet');
                cax.Units = 'centimeters';
                %cax.Position = [sum(sax(end).Position([1,3]))+0.15,sax(end).Position(2),0.15,sax(end).Position(4)];
 
                cax.Position = [sum(sax(end).Position([1])),...
                                sax(end).Position(2)-0.2,...
                                sax(end).Position(3),...
                                0.15];
                cax.XTick = [0,1];
                cax.XTickLabel = {'min','Max'};
                cax.Label.String = 'Hz';
                cax.Label.Units = 'centimeters';
                cax.Label.HorizontalAlignment = 'center';
                cax.Label.FontSize = 8;
                cax.Label.Position(1) = sax(end).Position(3)/2;
                cax.Label.Position(2) = cax.Label.Position(2)+0.1;
            end
        end% for s
        
% THETA phase histograms
        [yind, yOffSet, xind, xOffSet] = deal(yind, yOffSet, s+2, 1.3);
        sax(end+1) = axes(axOpts{:},                                            ...
                          'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                      fig.page.ypos(yind)+yOffSet,              ...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height]);
        bar(linspace(0,360,36),tphRUN(ismember(cluSessionMap,[t,unit],'rows'),:),'EdgeColor','none');
        sax(end).YTick = [];
        sax(end).XTick = [];
        if yind==1,
            title({'Theta','Phase'});
        end
        box(sax(end),'off');

% ACCG 
        [yind, yOffSet, xind, xOffSet] = deal(yind, yOffSet, s+3, 1.5);
        sax(end+1) = axes(axOpts{:},                                            ...
                          'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                      fig.page.ypos(yind)+yOffSet,              ...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height]);
        bar(tbin,accgHz(ismember(cluSessionMap,[t,unit],'rows'),:),'EdgeColor','none');
        sax(end).YTick = [];
        sax(end).XTick = [];
        xlim([-20,20]);
        box(sax(end),'off');
        if yind==1,
            title('ACCG');
        end
        
        if yind == 5
            xlabel('time');
            sax(end).XTick = [-20,0,20];            
            % Theta label
        end
        
        if yind == 3        
            sax(end).YAxisLocation = 'right';
            ylabel('Spike Count');
        end

        
        yind = yind+1;                    
    end
end

% DRAW line and arrows for theta to behavioral state decomposition
p = 9*(yind-1)-2;
% HORIZONTAL line
line(fax,                                                                                                ...
     [sax(end-p).Position(1)+fig.subplot.width-0.2,sax(end-(p-5)).Position(1)+fig.subplot.width/2+0.07], ...
      repmat(sax(end-p).Position(2)+sax(end-p).Position(4)+0.575,[1,2]),                                 ...
     'Color','k',                                                                                        ...
     'LineWidth',1);
% ARROWHEADS 
for a = 1:5,
    patch(fax,                                                                        ...
          repmat(sax(end-(p-a)).Position(1)+fig.subplot.width/2,[1,3])+[-0.07,0.07,0],...
          repmat(sum(sax(end-p).Position([2,4]))+0.475,[1,3])+[0.1,0.1,-0.09],        ...
          'k');
end

% THETA wave for theta phase histograms
[yind, yOffSet, xind, xOffSet] = deal(yind, 0, s+2, 1.3);
sax(end+1) = axes('Units','centimeters',                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,      ...
                    fig.page.ypos(yind)+yOffSet+fig.subplot.height/1.5,      ...
                    fig.subplot.width,            ...
                    fig.subplot.height/3],            ...
                  'FontSize', 8,                                ...
                  'LineWidth',1);
hold(sax(end),'on');
% ADD phase "bar"
plot(0:360,cos(circ_ang2rad(0:360)),'-k','LineWidth',2);
sax(end).Visible = 'off';                                    
sax(end).YTick = [];
sax(end).YTickLabels = [];
sax(end).XTick = [];
xlim([0,360]);
ylim([-1.5,1.5]);

text(360,-0.25, '2\pi', 'HorizontalAlignment','center',     ...
     'FontSize', 8, 'VerticalAlignment',  'middle');
text(180, 0.25, '\pi', 'HorizontalAlignment','center',      ...
     'FontSize', 8, 'VerticalAlignment',  'middle');                    
text(0,  -0.25, '0', 'HorizontalAlignment','center',        ...
     'FontSize', 8, 'VerticalAlignment',  'middle');                    


%%%>>>


%%%<<< GRAPHIC - erpPCA 

% $$$ sectionXind = 1;
% $$$ sectionYind = 5;
% $$$ 
% $$$ [yind, yOffSet, xind, xOffSet] = deal(sectionYind, -3.25, sectionXind, -0.25);
% $$$ sax(end+1) = axes(axOpts{:},                                            ...
% $$$                   'Position',[fig.page.xpos(xind)+xOffSet,              ...
% $$$                               fig.page.ypos(yind)+yOffSet,              ...
% $$$                               fig.subplot.width*1.5,                        ...
% $$$                               fig.subplot.height*1.5]);
% $$$ k = 1;
% $$$ hold('on');
% $$$ for tind = 1:size(expUnitsPP,1),
% $$$     t = expUnitsPP{tind,1}(1);
% $$$     for uind = 1:numel(expUnitsPP{tind,2}),
% $$$         unit = expUnitsPP{tind,2}(uind);
% $$$         brm = plot(bfsin{t},unit,1,false,[0,maxPfsRate],false,0.5,false);
% $$$         brm(~validDims) = nan;
% $$$         pax = pcolor(bins{end-1:end},brm');
% $$$         pax.ZData = -k.*ones(size(pax.ZData));
% $$$         set(pax,'EdgeColor','none');
% $$$         line(bins{end-1}([1,1]),bins{end}([1,end]),[-k;-k],'Color','k');
% $$$         line(bins{end-1}([1,end]),bins{end}([end,end]),[-k;-k],'Color','k');
% $$$         line(bins{end-1}([1,end]),bins{end}([1,1]),[-k;-k],'Color','k');    
% $$$         line(bins{end-1}([end,end]),bins{end}([1,end]),[-k;-k],'Color','k');
% $$$         k = k + 1;
% $$$     end
% $$$ end
% $$$ set(sax(end),'Visible',true);
% $$$ colormap(sax(end),'jet');
% $$$ view([20,20]);
% $$$ zlabel('Varimax PCA');
% $$$ sax(end).XTick = [];
% $$$ sax(end).YTick = [];
% $$$ sax(end).ZTick = [];
% $$$ 
% $$$ patch(fax,                                                                                    ...
% $$$       repmat(sax(end).Position(1)+fig.subplot.width/2,[1,3])+[-0.1,0.1,0]+0.2,                        ...
% $$$       repmat(sax(end).Position(2)-0.475,[1,3])+[0.1,0.1,-0.09],        ...
% $$$       'k');
% $$$ patch(fax,                                                                                    ...
% $$$       repmat(sax(end).Position(1)+fig.subplot.width/2,[1,3])+[-0.1,0.1,0]+0.2,                        ...
% $$$       repmat(sum(sax(end).Position([2,4]))+0.2,[1,3])+[0.1,0.1,-0.09],        ...
% $$$       'k');
% $$$           

%%%>>>



%%%<<< Panel - bhvRmap factor decomposition 
% PLOT first 3 eigen vectors
sectionXind = 2;
sectionYind = 7;
for v = 1:3,
% BHV fields MTAApfs
    [yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0.5, sectionXind+(v-1), 0.2.*v);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height]);
    
    eigVec = nan(binDims(end-1:end));
    eigVec(validDims) = LU(:,v);
    set(pcolor(bins{end-1:end},eigVec'),'EdgeColor','none');
    
    axis('xy');
    title(['Factor ',num2str(v)]);
    
    clear_axes_labels(sax(end));
    
    apply_colorbar(sax(end),'southoutside','jet');
    
    xlim([-1.7,0.5]);
    ylim([-0.5,1.7]);
    sax(end).Color = nanColor;
    
end
[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0.5, 1, 0);
sax(end+1) = axes(axOpts{:},                                            ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height]);
plot(VT(1:5,4),'-+')
xlim([0,6]);
xlabel('Factor');
ylabel('Exp Var');

%%%>>>




%%%<<< Panel - Asc vs Dsc factor score

sectionXind = 2;
sectionYind = 8;
clear('statsFSrTPP');
for f = 1:3,
    [yind, yOffSet, xind, xOffSet] = deal(sectionYind, -0.5, sectionXind+(f-1), 0.2.*f);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height]);
    %factor 2 boxplot phz ascending vs descending

    FSrTPP = {};
    mthr = 0.2;
    FSrTPP{1} = FSr(ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN>mthr & tppRUN > 0.4, f);    
    FSrTPP{2} = FSr(ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN>mthr & tppRUN <-0.4, f);
    FSrTPP{3} = FSr(ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN<=mthr & tppRUN > 0.4, f);    
    FSrTPP{4} = FSr(ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN<=mthr & tppRUN <-0.4, f);
    for c = 1:2,
        [statsFSrTPP(f,c).h,statsFSrTPP(f,c).p,statsFSrTPP(f,c).ci,statsFSrTPP(f,c).stats] = ...
            ttest2(FSrTPP{[1:2]+2.*(c==2)});
    end    
% $$$     FSrGRP = af(@(f,g) g.*ones(size(f{1})), FSrTPP,1:2);
% $$$     FSrGRP = cat(1,FSrGRP{:});
% $$$     FSrTPP = cat(1,FSrTPP{:});

    violin(FSrTPP,...
           'xlabels',{'mDsc','mAsc','Dsc','Asc'},...
           'facecolor',[0,1,1;0,1,1;1,0.5,0;1,0.5,0],...
           'plotmean', false,...
           'plotlegend',false);
        
    ylim(ylim+[0,1.5]);
    grid('on');
    hold('on');
    for c = 1:2    
        line(sax(end).XTick([1:2]+2*(c==2)),sax(end).YTick(end).*[1,1]-0.5,'Color','k')    
        if statsFSrTPP(f,c).h
            text(mean(sax(end).XTick([1:2]+2*(c==2))),...
                 sax(end).YTick(end),...
                 repmat('*',[1,floor(abs(log10(statsFSrTPP(f,c).p)))]),...
                 'HorizontalAlignment','center');
        else
            text(mean(sax(end).XTick([1:2]+2*(c==2))),...
                 sax(end).YTick(end)+0.1,...
                 'N.S.',...
                 'HorizontalAlignment','center');
        end
    end
    
    %title(['Factor ',num2str(f)]);
    if f == 2,
        title('Score vs TPP')
    end    
    
    if f ~= 3
        sax(end).XTickLabel = {};
    end
    %plot( sax(end).XTick(1).*ones([sum(FSrGRP==1),1]), FSrTPP(FSrGRP==1),'.');
    
end

%%%>>>



%%%<<< Panel - Theta Phase Preference (TPP)
sectionXind = 6;
sectionYind = 8;

[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0, sectionXind, -0.25);
sax(end+1) = axes(axOpts{:},                                            ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width.*2,                        ...
                              fig.subplot.height.*2]);
ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]);
pax = polarplot(tppRUN(ind),tprRUN(ind),'.');
hold('on');
polarplot([0,circ_mean(tppRUN(ind))],[0,circ_r(tppRUN(ind))],'-r');
rlim([0,0.7]);
pax.Parent.ThetaTick = [0,90,180,270];
pax.Parent.ThetaTickLabel = {'0','90','180','270'};
pax.Parent.GridColor = [0.1,0.1,0.1];
pax.Parent.GridAlphaMode = 'manual';
pax.Parent.GridAlpha = 0.3;
title({'TPP CA1'})


for e = 1:3,    
    trlUntPair = [expUnitsPP{eupInd(e,1),1},expUnitsPP{eupInd(e,1),2}];
    mind = find(ismember(cluSessionMap,trlUntPair,'rows'));
    polarscatter(tppRUN(mind),tprRUN(mind),30,eclr(e),mrk(e),'Filled');
end

% $$$ ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]);
% $$$ sectionXind = 7;
% $$$ sectionYind = 5;
% $$$ [yind, yOffSet, xind, xOffSet] = deal(sectionYind, -2.75, sectionXind, 0.5);
% $$$ sax(end+1) = axes(axOpts{:},                                            ...
% $$$                   'Position',[fig.page.xpos(xind)+xOffSet,              ...
% $$$                               fig.page.ypos(yind)+yOffSet,              ...
% $$$                               fig.subplot.width,                        ...
% $$$                               fig.subplot.height]);
% $$$ ind = ~ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]);
% $$$ polarplot(tppRUN(ind),tprRUN(ind),'.')
% $$$ hold('on');
% $$$ polarplot([0,circ_mean(tppRUN(ind))],[0,circ_r(tppRUN(ind))],'-r')
% $$$ title({'TPP CA3'})
%%%>>>



%%%<<< Panel - zscrBhvInfo vs zscrPosInfo
sectionXind = 9;
sectionYind = 8;

[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0, sectionXind, 0);
sax(end+1) = axes(axOpts{:},                                            ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width.*2,                        ...
                              fig.subplot.height.*2]);
hold(sax(end),'on');
plot(log10(zsinfNP),log10(zsinfNB),'.')
title({'PosInfo vs','BhvInfo'})
xlabel('log10(z-score)');
ylabel('log10(z-score)');
xlim(xlim+[-0.1,0.1]);
ylim(ylim+[-0.1,0.1]+0.2);
grid('on');
for e = 1:3,
    trlUntPair = [expUnitsPP{eupInd(e,1),1},expUnitsPP{eupInd(e,1),2}];
    mind = find(ismember(cluSessionMap,trlUntPair,'rows'));
    scatter(log10(zsinfNP(mind)),log10(zsinfNB(mind)),30,eclr(e),mrk(e),'Filled');
end
box(sax(end),'on');

%%%>>>



%%%<<< Panel - Example Interneurons by Behavior States and Maze Postition

sectionXind = 1;
sectionYind = 12;
% eupInd -> first examples ( expUnitsPP )
xOffSet = -1;
for  e = 1:size(eupInd,1),
    [yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0.5, sectionXind+e*3-3, 0);
    sax(end+1) = axes(axOpts{:},                                                         ...
                      'Position', [fig.page.xpos(xind)+xOffSet,                          ...
                                   fig.page.ypos(yind)+yOffSet,                          ...
                                   fig.subplot.width*3+fig.subplot.horizontalPadding,    ...
                                   fig.subplot.height*3+fig.subplot.verticalPadding]);
    hold('on');

    bhvSpcWidth = 22;
    bhvSpcHeight = 22;
    reducedMaskBhv = maskBhv(4:25,4:25);

    trlUntPair = [expUnitsPP{eupInd(e,1),1},expUnitsPP{eupInd(e,1),2}];

    mind = find(ismember(cluSessionMap,trlUntPair,'rows'));
    rmap = rmapa(:,:,:,:,mind);
    rmapEX = nan([bhvSpcWidth*4,bhvSpcHeight*4]);
    for x = 2:5,
        for y = 2:5,
            if (x==2&&y==2) || (x==2&&y==5)|| (x==5&&y==5)|| (x==5&&y==2)
                continue
            end
            trmap = sq(rmap(x,y,4:25,4:25));
            trmap(~reducedMaskBhv) = nan;
            rmapEX(( bhvSpcWidth*(x-2)+1):(bhvSpcWidth*(x-2) +bhvSpcWidth), ...
                   (bhvSpcHeight*(y-2)+1):(bhvSpcHeight*(y-2)+bhvSpcHeight)) = trmap;
        end
    end
    clim = [min(nonzeros(rmapEX(:))), max(nonzeros(rmapEX(:)))];

    pax = pcolor(rmapEX');
    pax.EdgeColor='none';

    axis('xy');
    Lines([22:22:22*4]+0.5,[],'k');
    Lines([],[22:22:22*4]+0.5,'k');
    xlim([0,22*4+0.5]);
    ylim([0,22*4+0.5]);
    
    scatter(11.5,78,40,eclr(e),mrk(e),'Filled');
    
    clear_axes_labels(sax(end));
    
    apply_colorbar(sax(end),'southoutside','jet');

    sax(end).Color = nanColor;    
    box(sax(end),'on');
% $$$     axes(fax);
% $$$     rectangle('Position',sax(end).Position,'LineWidth',1);
end



%%%>>>

%%%<<< Panel - example map corr vs distance
sectionXind = 1;
sectionYind = 15;
[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0.75, sectionXind, 0);
sax(end+1) = axes(axOpts{:},                                                         ...
                  'Position', [fig.page.xpos(xind)+xOffSet,                          ...
                               fig.page.ypos(yind)+yOffSet,                          ...
                               fig.subplot.width.*2,                                 ...
                               fig.subplot.height.*2]);
hold('on');
for e = 1:3,
    trlUntPair = [expUnitsPP{eupInd(e,1),1},expUnitsPP{eupInd(e,1),2}];
    mind = find(ismember(cluSessionMap,trlUntPair,'rows'));
    scatter(dst+randn([1,66])./30,                                     ... % x: distance
            rho(mind,:),                                               ... % y: map correlation
            10,                                                        ... % marker size
            eclr(e),                                                   ... % marker color
            mrk(e),                                                    ... % marker shape
            'Filled');
    xlim([0.95,3.1]);
    line(xlim,polyval(rhoDRPoly(mind,:),xlim),'LineWidth',1,'color',eclr(e));
    drawnow();
end

sax(end).XTickLabel = {20,40,60};
xlim([0.8,3.4]);

ylabel('map correlation');
xlabel('map distance (cm)');
grid('on');        

%%%>>>


%%%<<< Panel - cdf dist
sectionXind = 4;
sectionYind = 15;
[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0.75, sectionXind, 0);
sax(end+1) = axes(axOpts{:},                                                         ...
                  'Position', [fig.page.xpos(xind)+xOffSet,                          ...
                               fig.page.ypos(yind)+yOffSet,                          ...
                               fig.subplot.width.*2,                                 ...
                               fig.subplot.height.*2]);
hold('on');
ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]);
scatter(diff(rhoDRratio(ind,[2,1],3),1,2)+randn([sum(ind),1])./scl,    ... % x: near corr ratio diff
        diff(rhoDRratio(ind,[3,1],3),1,2)+randn([sum(ind),1])./scl,    ... % y: far  corr ratio diff
        5,                                                     ... % marker size
        rhoDR(ind),                                                    ... % marker color  rhoDRratio(:,1,3)
        'filled');
line([-0.2,1],[-0.2,1],'color','k');

for e = 1:3,
    trlUntPair = [expUnitsPP{eupInd(e,1),1},expUnitsPP{eupInd(e,1),2}];
    mind = find(ismember(cluSessionMap,trlUntPair,'rows'));
    scatter(diff(rhoDRratio(mind,[2,1],3),1,2)+randn([1,1])./scl,      ... %x
            diff(rhoDRratio(mind,[3,1],3),1,2)+randn([1,1])./scl,      ... %y
            30,                                                        ... marker size
            eclr(e),                                                   ... marker color
            mrk(e),                                                    ... marker shape
            'Filled');
end
xlim([-0.2,1]);
ylim([-0.2,1]);
apply_colorbar(sax(end),'eastoutside','jet');
box(sax(end),'on');
% $$$ axes(fax);
% $$$ rectangle('Position',sax(end).Position,'LineWidth',1);

%%%>>>








% PRESERVE figure layout
preserve_figure_layout(hfig);

set_figure_layout_mode(hfig,'absolute');
set_figure_layout_mode(hfig,'normalized');
set_figure_layout_mode(hfig,'fixedaspect');

reset_figure_layout(hfig);

%%%>>>



%%%<<< MAIN FIGURE landscape -----------------------------------------------------------------------

% bhvState spatial ratemap examples (placefields)
% place-bhehavior ratemaps          (space binned behaviorRatemaps)
% spatial and behavioral variation 
% Factors and explained variance

%%%<<< setup

[hfig,fig,fax,sax] = set_figure_layout(figure(666002),'A4','landscape',[],1.5,1.5,0.1,0.2);
expUnitsPP = {[ 3],[124];...
              [18],[ 59];...
              [4],[ 60];...
              [22],[ 49];...
              [25],[ 7]};

eupInd = [5,1;...
          3,1;...
          2,1];
eclr = 'rmg';
mrk  = '^so';

numUnits = sum(cellfun(@numel,expUnitsPP(:,2)));
nanColor = [0.25,0.25,0.25];
labels = {{'Theta','Behavior'},{'Theta','Arena'},'Rear','H Loc','H Pause','L Loc','L Pause'};
maskBhv = reshape(validDims,bfsin{1}.adata.binSizes');
axOpts = {'Units',                 'centimeters',...
          'FontSize',              8,            ...
          'LineWidth',             1,            ...
          'PlotBoxAspectRatioMode','manual'};

%%%>>>


%%%<<< Panel - Example Interneuorns by Behavior State

[yind, yOffSet, xind, xOffSet] = deal(1, 0, 1, 0);
for tind = 1:size(expUnitsPP,1),
    t = expUnitsPP{tind,1}(1);
    for uind = 1:numel(expUnitsPP{tind,2}),
        unit = expUnitsPP{tind,2}(uind);
        
        % GET maximum firing rate of all states
        maxPfsRate = max(cell2mat(cf(@(p,u,m)                                   ...
                                  p.maxRate(u,true&numel(m)==1,1,[],[],[],m),   ...
                                  [pfs{t},{bfsin{t}}],                          ...
                                  [repmat({unit},[1,numel(pfs{t})+1])],         ...
                                  [repmat({true},[1,numel(pfs{t})]),{maskBhv(:)}])));
        minPfsRate = min(cell2mat(cf(@(p,u,m)                                   ...
                                  minRate(p,u,true&numel(m)==1,1,[],[],[],m),   ...
                                  [pfs{t}],                          ...
                                  [repmat({unit},[1,numel(pfs{t})])],         ...
                                  [repmat({true},[1,numel(pfs{t})])])));

        
% BHV fields MTAApfs
        [yind, yOffSet, xind, xOffSet] = deal(yind, yOffSet, 1, 0);
        sax(end+1) = axes(axOpts{:},                                            ...
                          'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                      fig.page.ypos(yind)+yOffSet,              ...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height]);
        hold(sax(end),'on');
        plot(bfsin{t},unit,1,false,[minPfsRate,maxPfsRate],false,0.5,false,[],@jet,maskBhv,nanColor);            
        %plot(pfs{t}{s},unit,1,false,[0,maxPfsRate],true,0.5,false,interpParPfs,@jet,[],nanColor);
        sax(end).YTickLabel = {};
        sax(end).XTickLabel = {};
        box(sax(end),'on');
        
        xlim([-1.7,0.5]);
        ylim([-0.5,1.7]);
        
        exampleUnitInd = find(yind==eupInd(:,1));
        if ~isempty(exampleUnitInd),
            scatter(0.25,1.5,30,eclr(exampleUnitInd),mrk(exampleUnitInd),'filled'),
        end

% ADD Behavior space scale bars                
        if yind == 1
            title(labels{xind},'FontSize',8,'FontWeight','normal');
            line([-1.65,-1.2],[1.65,1.65],'Color','w','LineWidth',1);
            line([-1.65,-1.65],[1.2,1.65],'Color','w','LineWidth',1);        
        end
        
            
        end 
        if yind == 3
            ylabel(sax(end),'Body Pitch');
        end
        if yind == numUnits
            xlabel(sax(end),'Head Pitch');
        end
        
        
        for s = 1:numel(pfs{t}),
% PLACEFIELDS MTAApfs
            [yind, yOffSet, xind, xOffSet] = deal(yind, yOffSet, s+1, 0.5+0.25*double(s~=1));
            sax(end+1) = axes(axOpts{:},                                            ...
                              'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                          fig.page.ypos(yind)+yOffSet,              ...
                                          fig.subplot.width,                        ...
                                          fig.subplot.height]);

            plot(pfs{t}{s},unit,1,false,[minPfsRate,maxPfsRate],true,0.5,false,[],@jet,[],nanColor);            
            %plot(pfs{t}{s},unit,1,false,[0,maxPfsRate],true,0.5,false,interpParPfs,@jet,[],nanColor);
            sax(end).YTickLabel = {};
            sax(end).XTickLabel = {};
            
            if yind == 1, 
                title(labels{xind});
            end 
% ADD Arena scale bars        
            if yind == 1 && s ==1,
                line([-490,-250],[490,490],'Color','w','LineWidth',3);
                line([-490,-490],[250,490],'Color','w','LineWidth',3);        
            end
            
% ADD spatial grid to explain later 4D tuning maps
            if yind == numUnits && s ==1,            
                gridXY = bins{1}(2:6)-diff(bins{1}(1:2))/2;
                for x = gridXY'
                    for y = gridXY'
                        line([x,x],[x,y],'LineWidth',1,'Color','w');
                        line([x,y],[x,x],'LineWidth',1,'Color','w');
                    end
                end
                gridSax = sax(end);
            end
            
            

            if s==6,
                text(sax(end),...
                     492,                                  ... xpos
                     -392,                                  ... ypos
                     num2str(round(maxPfsRate)),           ... rate
                                'FontSize', 8,             ...
                                   'Color', [1,1,1],       ...
                     'HorizontalAlignment', 'right');
                text(sax(end),...
                     -472,                                 ... xpos
                     -392,                                  ... ypos
                     num2str(round(minPfsRate)),           ... rate
                                'FontSize', 8,             ...
                                   'Color', [1,1,1],       ...
                     'HorizontalAlignment', 'left');
            end

            % ADD arena scale bar
% $$$             if  yind == numUnits && s == 1,            
% $$$                 line(fax,...
% $$$                      [sax(end).Position(1),sax(end).Position(1)+sax(end).Position(3)./5],...
% $$$                      [sax(end).Position(2),sax(end).Position(2)]-0.1,...
% $$$                      'LineWidth',2,...
% $$$                      'Color','k');
% $$$                 %line(fax,...
% $$$                 %     [sax(end).Position(1),sax(end).Position(1)]-0.1,...
% $$$                 %     [sax(end).Position(2),sax(end).Position(2)+sax(end).Position(4)./5],...
% $$$                 %     'LineWidth',1,...
% $$$                 %     'Color','k');
% $$$                 text(fax,...
% $$$                      sax(end).Position(1),...
% $$$                      sax(end).Position(2)-0.3,...
% $$$                      '20 cm',...
% $$$                      'HorizontalAlignment','left');
% $$$             end
            
            % ADD colorbar
            if  yind == numUnits && s == numStates,
                cax = colorbar(sax(end),'southoutside');
                colormap(cax,'jet');
                cax.Units = 'centimeters';
                %cax.Position = [sum(sax(end).Position([1,3]))+0.15,sax(end).Position(2),0.15,sax(end).Position(4)];
 
                cax.Position = [sum(sax(end).Position([1])),...
                                sax(end).Position(2)-0.2,...
                                sax(end).Position(3),...
                                0.15];
                cax.XTick = [0,1];
                cax.XTickLabel = {'min','Max'};
                cax.Label.String = 'Hz';
                cax.Label.Units = 'centimeters';
                cax.Label.HorizontalAlignment = 'center';
                cax.Label.FontSize = 8;
                cax.Label.Position(1) = sax(end).Position(3)/2;
                cax.Label.Position(2) = cax.Label.Position(2)+0.1;
            end
        end% for s
        
% THETA phase histograms
        [yind, yOffSet, xind, xOffSet] = deal(yind, yOffSet, s+2, 1.3);
        sax(end+1) = axes(axOpts{:},                                            ...
                          'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                      fig.page.ypos(yind)+yOffSet,              ...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height]);
        bar(linspace(0,360,36),tphRUN(ismember(cluSessionMap,[t,unit],'rows'),:),'EdgeColor','none');
        sax(end).YTick = [];
        sax(end).XTick = [];
        if yind==1,
            title({'Theta','Phase'});
        end
        box(sax(end),'off');

% ACCG 
        [yind, yOffSet, xind, xOffSet] = deal(yind, yOffSet, s+3, 1.5);
        sax(end+1) = axes(axOpts{:},                                            ...
                          'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                      fig.page.ypos(yind)+yOffSet,              ...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height]);
        bar(tbin,accgHz(ismember(cluSessionMap,[t,unit],'rows'),:),'EdgeColor','none');
        sax(end).YTick = [];
        sax(end).XTick = [];
        xlim([-20,20]);
        box(sax(end),'off');
        if yind==1,
            title('ACCG');
        end
        
        if yind == 5
            xlabel('time');
            sax(end).XTick = [-20,0,20];            
            % Theta label
        end
        
        if yind == 3        
            sax(end).YAxisLocation = 'right';
            ylabel('Spike Count');
        end

        
        yind = yind+1;                    
    end
end

% DRAW line and arrows for theta to behavioral state decomposition
p = 9*(yind-1)-2;
% HORIZONTAL line
line(fax,                                                                                                ...
     [sax(end-p).Position(1)+fig.subplot.width-0.2,sax(end-(p-5)).Position(1)+fig.subplot.width/2+0.07], ...
      repmat(sax(end-p).Position(2)+sax(end-p).Position(4)+0.575,[1,2]),                                 ...
     'Color','k',                                                                                        ...
     'LineWidth',1);
% ARROWHEADS 
for a = 1:5,
    patch(fax,                                                                        ...
          repmat(sax(end-(p-a)).Position(1)+fig.subplot.width/2,[1,3])+[-0.07,0.07,0],...
          repmat(sum(sax(end-p).Position([2,4]))+0.475,[1,3])+[0.1,0.1,-0.09],        ...
          'k');
end

% THETA wave for theta phase histograms
[yind, yOffSet, xind, xOffSet] = deal(yind, 0, s+2, 1.3);
sax(end+1) = axes('Units','centimeters',                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,      ...
                    fig.page.ypos(yind)+yOffSet+fig.subplot.height/1.5,      ...
                    fig.subplot.width,            ...
                    fig.subplot.height/3],            ...
                  'FontSize', 8,                                ...
                  'LineWidth',1);
hold(sax(end),'on');
% ADD phase "bar"
plot(0:360,cos(circ_ang2rad(0:360)),'-k','LineWidth',2);
sax(end).Visible = 'off';                                    
sax(end).YTick = [];
sax(end).YTickLabels = [];
sax(end).XTick = [];
xlim([0,360]);
ylim([-1.5,1.5]);

text(360,-0.25, '2\pi', 'HorizontalAlignment','center',     ...
     'FontSize', 8, 'VerticalAlignment',  'middle');
text(180, 0.25, '\pi', 'HorizontalAlignment','center',      ...
     'FontSize', 8, 'VerticalAlignment',  'middle');                    
text(0,  -0.25, '0', 'HorizontalAlignment','center',        ...
     'FontSize', 8, 'VerticalAlignment',  'middle');                    


%%%>>>


%%%<<< GRAPHIC - erpPCA 

% $$$ sectionXind = 1;
% $$$ sectionYind = 5;
% $$$ 
% $$$ [yind, yOffSet, xind, xOffSet] = deal(sectionYind, -3.25, sectionXind, -0.25);
% $$$ sax(end+1) = axes(axOpts{:},                                            ...
% $$$                   'Position',[fig.page.xpos(xind)+xOffSet,              ...
% $$$                               fig.page.ypos(yind)+yOffSet,              ...
% $$$                               fig.subplot.width*1.5,                        ...
% $$$                               fig.subplot.height*1.5]);
% $$$ k = 1;
% $$$ hold('on');
% $$$ for tind = 1:size(expUnitsPP,1),
% $$$     t = expUnitsPP{tind,1}(1);
% $$$     for uind = 1:numel(expUnitsPP{tind,2}),
% $$$         unit = expUnitsPP{tind,2}(uind);
% $$$         brm = plot(bfsin{t},unit,1,false,[0,maxPfsRate],false,0.5,false);
% $$$         brm(~validDims) = nan;
% $$$         pax = pcolor(bins{end-1:end},brm');
% $$$         pax.ZData = -k.*ones(size(pax.ZData));
% $$$         set(pax,'EdgeColor','none');
% $$$         line(bins{end-1}([1,1]),bins{end}([1,end]),[-k;-k],'Color','k');
% $$$         line(bins{end-1}([1,end]),bins{end}([end,end]),[-k;-k],'Color','k');
% $$$         line(bins{end-1}([1,end]),bins{end}([1,1]),[-k;-k],'Color','k');    
% $$$         line(bins{end-1}([end,end]),bins{end}([1,end]),[-k;-k],'Color','k');
% $$$         k = k + 1;
% $$$     end
% $$$ end
% $$$ set(sax(end),'Visible',true);
% $$$ colormap(sax(end),'jet');
% $$$ view([20,20]);
% $$$ zlabel('Varimax PCA');
% $$$ sax(end).XTick = [];
% $$$ sax(end).YTick = [];
% $$$ sax(end).ZTick = [];
% $$$ 
% $$$ patch(fax,                                                                                    ...
% $$$       repmat(sax(end).Position(1)+fig.subplot.width/2,[1,3])+[-0.1,0.1,0]+0.2,                        ...
% $$$       repmat(sax(end).Position(2)-0.475,[1,3])+[0.1,0.1,-0.09],        ...
% $$$       'k');
% $$$ patch(fax,                                                                                    ...
% $$$       repmat(sax(end).Position(1)+fig.subplot.width/2,[1,3])+[-0.1,0.1,0]+0.2,                        ...
% $$$       repmat(sum(sax(end).Position([2,4]))+0.2,[1,3])+[0.1,0.1,-0.09],        ...
% $$$       'k');
% $$$           

%%%>>>



%%%<<< Panel - bhvRmap factor decomposition 
% PLOT first 3 eigen vectors
sectionXind = 13;
sectionYind = 1;
for v = 1:3,
% BHV fields MTAApfs
    [yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0.5, sectionXind+(v-1), 0.2.*v+0.25);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height]);
    
    eigVec = nan(binDims(end-1:end));
    eigVec(validDims) = LU(:,v);
    set(pcolor(bins{end-1:end},eigVec'),'EdgeColor','none');
    
    axis('xy');
    title(['Factor ',num2str(v)]);
    
    clear_axes_labels(sax(end));
    
    apply_colorbar(sax(end),'southoutside','jet');
    
    xlim([-1.7,0.5]);
    ylim([-0.5,1.7]);
    sax(end).Color = nanColor;
    
end
[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0.5, 12, 0);
sax(end+1) = axes(axOpts{:},                                            ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height]);
plot(VT(1:5,4),'-+')
xlim([0,6]);
xlabel('Factor');
ylabel('Exp Var');

%%%>>>



%%%<<< Panel - Asc vs Dsc factor score

sectionXind = 13;
sectionYind = 2;
clear('statsFSrTPP');
for f = 1:3,
    [yind, yOffSet, xind, xOffSet] = deal(sectionYind, -0.5, sectionXind+(f-1), 0.2.*f+0.25);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height]);
    %factor 2 boxplot phz ascending vs descending

    FSrTPP = {};
    mthr = 0.2;
    FSrTPP{1} = FSr(ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN>mthr & tppRUN > 0.4, f);    
    FSrTPP{2} = FSr(ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN>mthr & tppRUN <-0.4, f);
    FSrTPP{3} = FSr(ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN<=mthr & tppRUN > 0.4, f);    
    FSrTPP{4} = FSr(ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN<=mthr & tppRUN <-0.4, f);
    for c = 1:2,
        [statsFSrTPP(f,c).h,statsFSrTPP(f,c).p,statsFSrTPP(f,c).ci,statsFSrTPP(f,c).stats] = ...
            ttest2(FSrTPP{[1:2]+2.*(c==2)});
    end    
% $$$     FSrGRP = af(@(f,g) g.*ones(size(f{1})), FSrTPP,1:2);
% $$$     FSrGRP = cat(1,FSrGRP{:});
% $$$     FSrTPP = cat(1,FSrTPP{:});

    violin(FSrTPP,...
           'xlabels',{'mDsc','mAsc','Dsc','Asc'},...
           'facecolor',[0,1,1;0,1,1;1,0.5,0;1,0.5,0],...
           'plotmean', false,...
           'plotlegend',false);
        
    ylim(ylim+[0,1.5]);
    grid('on');
    hold('on');
    for c = 1:2    
        line(sax(end).XTick([1:2]+2*(c==2)),sax(end).YTick(end).*[1,1]-0.5,'Color','k')    
        if statsFSrTPP(f,c).h
            text(mean(sax(end).XTick([1:2]+2*(c==2))),...
                 sax(end).YTick(end),...
                 repmat('*',[1,floor(abs(log10(statsFSrTPP(f,c).p)))]),...
                 'HorizontalAlignment','center');
        else
            text(mean(sax(end).XTick([1:2]+2*(c==2))),...
                 sax(end).YTick(end)+0.1,...
                 'N.S.',...
                 'HorizontalAlignment','center');
        end
    end
    
    %title(['Factor ',num2str(f)]);
    if f == 2,
        title('Score vs TPP')
    end    
    
    if f ~= 3
        sax(end).XTickLabel = {};
    end
    %plot( sax(end).XTick(1).*ones([sum(FSrGRP==1),1]), FSrTPP(FSrGRP==1),'.');
    
end

%%%>>>



%%%<<< Panel - Theta Phase Preference (TPP)
sectionXind = 12;
sectionYind = 5;

[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0, sectionXind, -0.25);
sax(end+1) = axes(axOpts{:},                                            ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width.*2,                        ...
                              fig.subplot.height.*2]);
ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]);
pax = polarplot(tppRUN(ind),tprRUN(ind),'.');
hold('on');
polarplot([0,circ_mean(tppRUN(ind))],[0,circ_r(tppRUN(ind))],'-r');
rlim([0,0.7]);
pax.Parent.ThetaTick = [0,90,180,270];
pax.Parent.ThetaTickLabel = {'0','90','180','270'};
pax.Parent.GridColor = [0.1,0.1,0.1];
pax.Parent.GridAlphaMode = 'manual';
pax.Parent.GridAlpha = 0.3;
title({'TPP CA1'})


for e = 1:3,    
    trlUntPair = [expUnitsPP{eupInd(e,1),1},expUnitsPP{eupInd(e,1),2}];
    mind = find(ismember(cluSessionMap,trlUntPair,'rows'));
    polarscatter(tppRUN(mind),tprRUN(mind),30,eclr(e),mrk(e),'Filled');
end

% $$$ ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]);
% $$$ sectionXind = 7;
% $$$ sectionYind = 5;
% $$$ [yind, yOffSet, xind, xOffSet] = deal(sectionYind, -2.75, sectionXind, 0.5);
% $$$ sax(end+1) = axes(axOpts{:},                                            ...
% $$$                   'Position',[fig.page.xpos(xind)+xOffSet,              ...
% $$$                               fig.page.ypos(yind)+yOffSet,              ...
% $$$                               fig.subplot.width,                        ...
% $$$                               fig.subplot.height]);
% $$$ ind = ~ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]);
% $$$ polarplot(tppRUN(ind),tprRUN(ind),'.')
% $$$ hold('on');
% $$$ polarplot([0,circ_mean(tppRUN(ind))],[0,circ_r(tppRUN(ind))],'-r')
% $$$ title({'TPP CA3'})
%%%>>>



%%%<<< Panel - zscrBhvInfo vs zscrPosInfo
sectionXind = 14;
sectionYind = 5;

[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0, sectionXind, 0.75);
sax(end+1) = axes(axOpts{:},                                            ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width.*2,                        ...
                              fig.subplot.height.*2]);
hold(sax(end),'on');
plot(log10(zsinfNP),log10(zsinfNB),'.')
title({'PosInfo vs','BhvInfo'})
xlabel('log10(z-score)');
ylabel('log10(z-score)');
xlim(xlim+[-0.1,0.1]);
ylim(ylim+[-0.1,0.1]+0.2);
grid('on');
for e = 1:3,
    trlUntPair = [expUnitsPP{eupInd(e,1),1},expUnitsPP{eupInd(e,1),2}];
    mind = find(ismember(cluSessionMap,trlUntPair,'rows'));
    scatter(log10(zsinfNP(mind)),log10(zsinfNB(mind)),30,eclr(e),mrk(e),'Filled');
end
box(sax(end),'on');

%%%>>>



%%%<<< Panel - Example Interneurons by Behavior States and Maze Postition

sectionXind = 1;
sectionYind = 9;
% eupInd -> first examples ( expUnitsPP )
xOffSet = -1;
for  e = 1:size(eupInd,1),
    [yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0, sectionXind+e*3-3, 0);
    sax(end+1) = axes(axOpts{:},                                                         ...
                      'Position', [fig.page.xpos(xind)+xOffSet,                          ...
                                   fig.page.ypos(yind)+yOffSet,                          ...
                                   fig.subplot.width*3+fig.subplot.horizontalPadding,    ...
                                   fig.subplot.height*3+fig.subplot.verticalPadding]);
    hold('on');

    bhvSpcWidth = 22;
    bhvSpcHeight = 22;
    reducedMaskBhv = maskBhv(4:25,4:25);

    trlUntPair = [expUnitsPP{eupInd(e,1),1},expUnitsPP{eupInd(e,1),2}];

    mind = find(ismember(cluSessionMap,trlUntPair,'rows'));
    rmap = rmapa(:,:,:,:,mind);
    rmapEX = nan([bhvSpcWidth*4,bhvSpcHeight*4]);
    for x = 2:5,
        for y = 2:5,
            if (x==2&&y==2) || (x==2&&y==5)|| (x==5&&y==5)|| (x==5&&y==2)
                continue
            end
            trmap = sq(rmap(x,y,4:25,4:25));
            trmap(~reducedMaskBhv) = nan;
            rmapEX(( bhvSpcWidth*(x-2)+1):(bhvSpcWidth*(x-2) +bhvSpcWidth), ...
                   (bhvSpcHeight*(y-2)+1):(bhvSpcHeight*(y-2)+bhvSpcHeight)) = trmap;
        end
    end
    clim = [min(nonzeros(rmapEX(:))), max(nonzeros(rmapEX(:)))];

    pax = pcolor(rmapEX');
    pax.EdgeColor='none';

    axis('xy');
    Lines([22:22:22*4]+0.5,[],'k');
    Lines([],[22:22:22*4]+0.5,'k');
    xlim([0,22*4+0.5]);
    ylim([0,22*4+0.5]);
    
    scatter(11.5,78,40,eclr(e),mrk(e),'Filled');
    
    clear_axes_labels(sax(end));
    
    apply_colorbar(sax(end),'southoutside','jet');

    sax(end).Color = nanColor;    
    box(sax(end),'on');

% Connect theta pfs with grid to 4D example
    if e == 1
        axes(fax)
        line([sax(end).Position(1),gridSax(end).Position(1)+0.2],...
             [sum(sax(end).Position([2,4])),gridSax.Position(2)+0.2],...
             'LineWidth',2,...
             'Color','k');
        line([sum(sax(end).Position([1,3])),sum(gridSax.Position([1,3]))-0.2],...
             [sum(sax(end).Position([2,4])),gridSax.Position(2)+0.2],...
             'LineWidth',2,...
             'Color','k');
    end
end



%%%>>>



%%%<<< Panel - example map corr vs distance
sectionXind = 11;
sectionYind = 9;
[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0, sectionXind, 0);
sax(end+1) = axes(axOpts{:},                                                         ...
                  'Position', [fig.page.xpos(xind)+xOffSet,                          ...
                               fig.page.ypos(yind)+yOffSet,                          ...
                               fig.subplot.width.*2,                                 ...
                               fig.subplot.height.*2]);
hold('on');
for e = 1:3,
    trlUntPair = [expUnitsPP{eupInd(e,1),1},expUnitsPP{eupInd(e,1),2}];
    mind = find(ismember(cluSessionMap,trlUntPair,'rows'));
    scatter(dst+randn([1,66])./30,                                     ... % x: distance
            rho(mind,:),                                               ... % y: map correlation
            10,                                                        ... % marker size
            eclr(e),                                                   ... % marker color
            mrk(e),                                                    ... % marker shape
            'Filled');
    xlim([0.95,3.1]);
    line(xlim,polyval(rhoDRPoly(mind,:),xlim),'LineWidth',1,'color',eclr(e));
    drawnow();
end

sax(end).XTickLabel = {20,40,60};
xlim([0.8,3.4]);

ylabel('map correlation');
xlabel('map distance (cm)');
grid('on');        

%%%>>>



%%%<<< Panel - Unkown title
sectionXind = 14;
sectionYind = 9;
[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0, sectionXind, 0.5);
sax(end+1) = axes(axOpts{:},                                           ...
                  'Position', [fig.page.xpos(xind)+xOffSet,            ...
                               fig.page.ypos(yind)+yOffSet,            ...
                               fig.subplot.width.*2,                   ...
                               fig.subplot.height.*2]);
hold('on');
ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]);
scatter(diff(rhoDRratio(ind,[2,1],3),1,2)+randn([sum(ind),1])./scl,    ... % x: near corr ratio diff
        diff(rhoDRratio(ind,[3,1],3),1,2)+randn([sum(ind),1])./scl,    ... % y: far  corr ratio diff
        5,                                                             ... % marker size
        rhoDR(ind),                                                    ... % marker color  rhoDRratio(:,1,3)
        'filled');
line([-0.2,1],[-0.2,1],'color','k');

for e = 1:3,
    trlUntPair = [expUnitsPP{eupInd(e,1),1},expUnitsPP{eupInd(e,1),2}];
    mind = find(ismember(cluSessionMap,trlUntPair,'rows'));
    scatter(diff(rhoDRratio(mind,[2,1],3),1,2)+randn([1,1])./scl,      ... %x
            diff(rhoDRratio(mind,[3,1],3),1,2)+randn([1,1])./scl,      ... %y
            30,                                                        ... marker size
            eclr(e),                                                   ... marker color
            mrk(e),                                                    ... marker shape
            'Filled');
end
xlim([-0.2,1]);
ylim([-0.2,1]);
apply_colorbar(sax(end),'eastoutside','jet');
box(sax(end),'on');


pdSix = fitdist(diff(rhoDRratio(ind,[3,1],3),1,2)+randn([sum(ind),1])./scl,'Kernel','BandWidth',0.05);
x = linspace([ylim,40]);
ySix = pdf(pdSix,x);
ySix = ySix/sum(ySix);

hold('on');
plot(1-ySix.*5,x,'k-','LineWidth',1)
fill([1,1-ySix.*5,1],x([1,1:end,end]),'c');


% $$$ axes(fax);
% $$$ rectangle('Position',sax(end).Position,'LineWidth',1);

%%%>>>




% PRESERVE figure layout
% $$$ preserve_figure_layout(hfig);
% $$$ set_figure_layout_mode(hfig,'absolute');
% $$$ set_figure_layout_mode(hfig,'normalized');
% $$$ set_figure_layout_mode(hfig,'fixedaspect');
% $$$ reset_figure_layout(hfig);

% END MAIN FIGURE landscape
%%%>>>



%%%<<< SUPFIG interneuron intermap correlation change with distance ------------------------------- 


[hfig,fig,fax,sax] = set_figure_layout(figure(666004),'A4','portrait',[],1.5,1.5,1,1);
axOpts = {'Units',                 'centimeters',...
          'FontSize',              8,            ...
          'LineWidth',             1,            ...
          'PlotBoxAspectRatioMode','manual'};


%%%<<< Panel - place and behavior rate maps
supExampleUnits = [115,134,194];
gridSax = gobjects([1,3]);
for e = 1:3;
    tid  = cluSessionMap(supExampleUnits(e),1);
    uid  = cluSessionMap(supExampleUnits(e),2);
% PLOT theta behaivor field
    [yind, yOffSet, xind, xOffSet] = deal(1, 1, 1+(e-1)*2, 1);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height]);
    maxPfsRate = pfs{tid}{1}.maxRate(uid,true,1,[],[],[],true);
    minPfsRate = minRate(pfs{tid}{1},uid,true,1,[],[],[],true);

    hold(sax(end),'on');
    plot(bfsin{tid},uid,1,false,[minPfsRate,maxPfsRate],false,0.5,false,[],@jet,maskBhv,nanColor);
    sax(end).YTickLabel = {};
    sax(end).XTickLabel = {};
    sax(end).Color = nanColor;
    box(sax(end),'on');
    xlim([-1.7,0.5]);
    ylim([-0.5,1.7]);

% PLOT theta place field
    [yind, yOffSet, xind, xOffSet] = deal(1, 1, 2+(e-1)*2, 1);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                        fig.page.ypos(yind)+yOffSet,              ...
                        fig.subplot.width,                        ...
                        fig.subplot.height]);

    plot(pfs{tid}{1},uid,1,'',[minPfsRate,maxPfsRate],'nanColor',nanColor);
    sax(end).YTickLabel = {};
    sax(end).XTickLabel = {};
    sax(end).Color = nanColor;    
    xlim([-550,550]);
    ylim([-550,550]);
    box(sax(end),'on');
    %gridSax(e) = sax(end);
% ADD min and max rates to pfs plot
    text(sax(end),                             ... ax
         492,                                  ... xpos
         460,                                  ... ypos
         num2str(round(maxPfsRate)),           ... rate
        'FontSize', 8,                         ...
        'Color', [1,1,1],                      ...
        'HorizontalAlignment', 'right');
    text(sax(end),                             ... ax
         -472,                                 ... xpos
         460,                                  ... ypos
         num2str(round(minPfsRate)),           ... rate
        'FontSize', 8,                         ...
        'Color', [1,1,1],                      ...
        'HorizontalAlignment', 'left');
        
% ADD spatial grid to explain later 4D tuning maps
        if e == 1
            gridXY = bins{1}(2:6)-diff(bins{1}(1:2))/2;
            for x = gridXY'
                for y = gridXY'
                    line([x,x],[x,y],'LineWidth',1,'Color','w');
                    line([x,y],[x,x],'LineWidth',1,'Color','w');
                end
            end
            gridSax = sax(end);
        end
% ADD colorbar
    if e == 3
        % ADD colorbar
        cax = colorbar(sax(end),'northoutside');
        colormap(cax,'jet');
        cax.Units = 'centimeters';
        %cax.Position = [sum(sax(end).Position([1,3]))+0.15,sax(end).Position(2),0.15,sax(end).Position(4)];
        
        cax.Position = [sax(end).Position([1]),...
                        sum(sax(end).Position([2,4]))+0.2,...
                        sax(end).Position(3),...
                        0.15];
        cax.XTick = [0,1];
        cax.XTickLabel = {'min','Max'};
        cax.Label.String = 'Hz';
        cax.Label.Units = 'centimeters';
        cax.Label.HorizontalAlignment = 'center';
        cax.Label.FontSize = 8;
        cax.Label.Position(1) = sax(end).Position(3)/2;
        cax.Label.Position(2) = cax.Label.Position(2)+0.1;
    end
            
end
%%%>>>


%%%<<< Panel - Example Interneurons by Behavior States and Maze Postition

sectionXind = 1;
sectionYind = 3;
for  e = 1:numel(supExampleUnits),
    [yind, yOffSet, xind, xOffSet] = deal(sectionYind, 1, sectionXind+e*2-2, 1);
    sax(end+1) = axes(axOpts{:},                                                         ...
                      'Position', [fig.page.xpos(xind)+xOffSet,                          ...
                                   fig.page.ypos(yind)+yOffSet,                          ...
                                   fig.subplot.width*2+fig.subplot.horizontalPadding,    ...
                                   fig.subplot.height*2+fig.subplot.verticalPadding]);
    hold('on');

    bhvSpcWidth = 22;
    bhvSpcHeight = 22;
    reducedMaskBhv = maskBhv(4:25,4:25);
    trlUntPair  = cluSessionMap(supExampleUnits(e),:);

    mind = find(ismember(cluSessionMap,trlUntPair,'rows'));
    rmap = rmapa(:,:,:,:,mind);
    rmapEX = nan([bhvSpcWidth*4,bhvSpcHeight*4]);
    for x = 2:5,
        for y = 2:5,
            if (x==2&&y==2) || (x==2&&y==5)|| (x==5&&y==5)|| (x==5&&y==2)
                continue
            end
            trmap = sq(rmap(x,y,4:25,4:25));
            trmap(~reducedMaskBhv) = nan;
            rmapEX(( bhvSpcWidth*(x-2)+1):(bhvSpcWidth*(x-2) +bhvSpcWidth), ...
                   (bhvSpcHeight*(y-2)+1):(bhvSpcHeight*(y-2)+bhvSpcHeight)) = trmap;
        end
    end
    clim = [min(nonzeros(rmapEX(:))), max(nonzeros(rmapEX(:)))];

    pax = pcolor(rmapEX');
    pax.EdgeColor='none';

    axis('xy');
    Lines([22:22:22*4]+0.5,[],'k');
    Lines([],[22:22:22*4]+0.5,'k');
    xlim([0,22*4+0.5]);
    ylim([0,22*4+0.5]);

    if e == 1
        hold('on');
        plot([33,55],[77,77],'-',...
             'LineWidth',2,...
             'MarkerEdgeColor',[0.929 0.694 0.125],...
             'MarkerFaceColor',[0.929 0.694 0.125]);
        plot([33,77],[77,55],'-',...
             'LineWidth',2,...             
             'MarkerEdgeColor',[0.85 0.325 0.098],...
             'MarkerFaceColor',[0.85 0.325 0.098]);
        plot([33,55],[77,11],'-',...
             'LineWidth',2,...             
             'MarkerEdgeColor',[0.929 0.694 0.125],...
             'MarkerFaceColor',[0.929 0.694 0.125]);
    end
    
% ADD marker
    scatter(11.5,78,40,eclr(e),mrk(e),'Filled');
    
    clear_axes_labels(sax(end));
    
    apply_colorbar(sax(end),'southoutside','jet');

    sax(end).Color = nanColor;    
    box(sax(end),'on');

    % Connect theta pfs with grid to 4D example
    if e == 1
    axes(fax)
    line([sax(end).Position(1),gridSax(e).Position(1)+0.2],...
         [sum(sax(end).Position([2,4])),gridSax(e).Position(2)+0.2],...
         'LineWidth',2,...
         'Color','k');
    line([sum(sax(end).Position([1,3])),sum(gridSax.Position([1,3]))-0.2],...
         [sum(sax(end).Position([2,4])),gridSax(e).Position(2)+0.2],...
         'LineWidth',2,...
         'Color','k');
    end

end



%%%>>>



%%%<<< Panel - Rate correlation

sectionXind = 1;
sectionYind = 5;
for  e = 1:numel(supExampleUnits),
    [yind, yOffSet, xind, xOffSet] = deal(sectionYind, 0, sectionXind+e*2-2, 1);
    sax(end+1) = axes(axOpts{:},                                                         ...
                      'Position', [fig.page.xpos(xind)+xOffSet,                          ...
                                   fig.page.ypos(yind)+yOffSet,                          ...
                                   fig.subplot.width*2+fig.subplot.horizontalPadding,    ...
                                   fig.subplot.height*2+fig.subplot.verticalPadding]);
    hold('on');
    x1 = nimap(:,u,1);
    hold('on');
    for m = [2,6,12];
        x2 = nimap(:,supExampleUnits(e),m);        
        ind = nniz(x1) & nniz(x2) & validDims;        
        plot(nimap(ind,supExampleUnits(e),1),nimap(ind,supExampleUnits(e),m),'.');
    end
    
    if e ==2 
        title('Correlation between behavior rate maps');
    end

    if e == 1
        legend({'20cm','40cm','60cm'},'Location','southeast');
        xlabel('rate Hz (Rate Map 1)');
        ylabel('rate Hz (Rate Map 2)');        
    end
end    
    

%%%>>>


%%%<<< Panel - Rate correlation

sectionXind = 1;
sectionYind = 7;
dstEdges = [0.5,1.5,2.5,3.5];
dstCenters = [1,2,3];
dsti = discretize(dst,dstEdges);
dstColors = [0.7,0.7,0.7;...
             0.8,0.8,0.8;...
             0.9,0.9,0.9];

for  e = 1:numel(supExampleUnits),
    [yind, yOffSet, xind, xOffSet] = deal(sectionYind, -1, sectionXind+e*2-2, 1);
    sax(end+1) = axes(axOpts{:},                                                         ...
                      'Position', [fig.page.xpos(xind)+xOffSet,                          ...
                                   fig.page.ypos(yind)+yOffSet,                          ...
                                   fig.subplot.width*2+fig.subplot.horizontalPadding,    ...
                                   fig.subplot.height*2+fig.subplot.verticalPadding]);
    hold('on');
    ylim([-0.5,1]);
    for m = 1:numel(dstCenters)
        patch([dstEdges(m),dstEdges(m),dstEdges(m+1),dstEdges(m+1)],...
              [min(ylim),max(ylim),max(ylim),min(ylim)],...
              dstColors(m,:),...
              'EdgeColor','none');
    end
    
    plot(dst,rho(supExampleUnits(e),:),'.');
    plot(dstCenters,...
         [mean(rho(supExampleUnits(e),dsti==1)),...
          mean(rho(supExampleUnits(e),dsti==2)),...
          mean(rho(supExampleUnits(e),dsti==3))],...
         '-+');
    if e ==2 
        title('Behavior Rate Map Correlation vs Distance ');
    end
    if e == 1
        xlabel('map distance (cm)');
        ylabel('map Correlation');
    end
    
    line(xlim(),[0.3,0.3],'Color','m');
    plot(dstCenters,...
         rhoDRratio(supExampleUnits(e),:,4),...
         '-+m');
    
    xlim([0.5,3.5]);    
    sax(end).XTick = dstCenters;
    sax(end).XTickLabel = {'20','40','60'};
    scatter(3,0.9,40,eclr(e),mrk(e),'Filled');
end    
    

%%%>>>

%%%<<< Panel - correlation proportion
sectionXind = 1;
sectionYind = 10;

    [yind, yOffSet, xind, xOffSet] = deal(sectionYind, 1, sectionXind, 1);
    sax(end+1) = axes(axOpts{:},                                                         ...
                      'Position', [fig.page.xpos(xind)+xOffSet,                          ...
                                   fig.page.ypos(yind)+yOffSet,                          ...
                                   fig.subplot.width*2+fig.subplot.horizontalPadding,    ...
                                   fig.subplot.height*2+fig.subplot.verticalPadding]);
    hold('on');

    plot(rhoDRPoly(:,1),rhoDR','.');
    for e = 1:3,
        scatter(rhoDRPoly(supExampleUnits(e),1),...
                rhoDR(supExampleUnits(e)),...
                20,...
                eclr(e),...
                mrk(e),...
                'Filled');    
    end
    xlabel('MapCorr-Dist Slope');
    ylabel('MapCorr-Dist Correlation');
    xlim([-0.45,0]);
%%%>>>

%%%<<< Panel - correlation proportion


sectionXind = 3;
sectionYind = 10;
[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 1, sectionXind, 1+0.1);
sax(end+1) = axes(axOpts{:},                                                            ...
                  'Position', [fig.page.xpos(xind)+xOffSet,                             ...
                               fig.page.ypos(yind)+yOffSet,                             ...
                               fig.subplot.width,      ...
                               fig.subplot.height.*2+fig.subplot.horizontalPadding]);
hold('on');
sax(end).PlotBoxAspectRatioMode = 'auto';
out = hist(diff(rhoDRratio(:,[2,1],:),1,2),binsDRe);
imagesc(thr(1:7),binsDRc,out(:,1:7));
line([0.3,0.3],ylim,'Color','m');
ylim([-0.2,1]);
%ylabel('Change in proportion samples above correlation threshold')    
% $$$ xlabel('Thresh');

title('\Delta[40,20]');

[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 1, sectionXind+1, 1-0.1);
sax(end+1) = axes(axOpts{:},                                                            ...
                  'Position', [fig.page.xpos(xind)+xOffSet,                             ...
                               fig.page.ypos(yind)+yOffSet,                             ...
                               fig.subplot.width,      ...
                               fig.subplot.height.*2+fig.subplot.horizontalPadding]);
hold('on');
sax(end).PlotBoxAspectRatioMode = 'auto';
sax(end).YTickLabel = {};
out = hist(diff(rhoDRratio(:,[3,1],:),1,2),binsDRe);
imagesc(thr(1:7),binsDRc,out(:,1:7));
line([0.3,0.3],ylim,'Color','m');
ylim([-0.2,1]);
title('\Delta[60,20]');
xlabel('Correlation Threshold');

% binned bhv ratemap spearman correlation 
% Correlation of spatially independent behavior-ratemaps

%%%>>>

%%%<<< Panel - Unknow Title
sectionXind = 5;
sectionYind = 10;
[yind, yOffSet, xind, xOffSet] = deal(sectionYind, 1, sectionXind, 1);
sax(end+1) = axes(axOpts{:},                                                            ...
                  'Position', [fig.page.xpos(xind)+xOffSet,                             ...
                               fig.page.ypos(yind)+yOffSet,                             ...
                               fig.subplot.width.*2+fig.subplot.horizontalPadding,      ...
                               fig.subplot.height.*2+fig.subplot.horizontalPadding]);
hold('on');
%ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]);
ind = true([225,1]);;
scatter(diff(rhoDRratio(ind,[2,1],3),1,2)+randn([sum(ind),1])./scl,    ... % x: near corr ratio diff
        diff(rhoDRratio(ind,[3,1],3),1,2)+randn([sum(ind),1])./scl,    ... % y: far  corr ratio diff
        5,                                                             ... % marker size
        rhoDR(ind),                                                    ... % marker color  rhoDRratio(:,1,3)
        'filled');
line([-0.2,1],[-0.2,1],'color','k');

for e = 1:3,
    mind = supExampleUnits(e);
    scatter(diff(rhoDRratio(mind,[2,1],3),1,2)+randn([1,1])./scl,      ... %x
            diff(rhoDRratio(mind,[3,1],3),1,2)+randn([1,1])./scl,      ... %y
            30,                                                        ... marker size
            eclr(e),                                                   ... marker color
            mrk(e),                                                    ... marker shape
            'Filled');
end
xlim([-0.2,1]);
ylim([-0.2,1]);
apply_colorbar(sax(end),'eastoutside','jet');
box(sax(end),'on');


pdSix = fitdist(diff(rhoDRratio(ind,[3,1],3),1,2)+randn([sum(ind),1])./scl,'Kernel','BandWidth',0.05);
x = linspace([ylim,40]);
ySix = pdf(pdSix,x);
ySix = ySix/sum(ySix);

hold('on');
plot(1-ySix.*5,x,'k-','LineWidth',1)
fill([1,1-ySix.*5,1],x([1,1:end,end]),'c');


% $$$ axes(fax);
% $$$ rectangle('Position',sax(end).Position,'LineWidth',1);

%%%>>>



%%%>>>


%%%<<< SUPFIG interneuron details ------------------------------------------------------------------

ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]);
figure,
histogram(rhoDR(ind));
%plot(tppRUN(ind),rhoDR(ind),'.');
ylabel('count');
xlabel('corr');
hold('on');
eupInd = [4,1;...
          2,1;...
          5,1];
for  e = 1:size(eupInd,1),
    trlUntPair = [expUnitsPP{eupInd(e,1),1},expUnitsPP{eupInd(e,1),2}];
end

u = 131;
Lines(rhoDR(u),[],'r');
text(rhoDR(u),45,num2str(u));

[~,sind] = sort(rhoDR,'descend');

hfigg = figure();
sp = gobjects([0,1]);
sp(1) = subplot(261);
histogram(rhoDR(ind));
sp(2) = subplot(262);
sp(3) = subplot(263);
hold('on');
plot(rhoDRMean,rhoDR,'.');
sp(4) = subplot(264);
sp(5) = subplot(265);

lax = [];
tax = []
u = sind(1);
while u ~= -1
    delete(lax);
    delete(tax);
    delete(pax);    
    cla(sp(2));
axes(sp(2));

% BHVPOS fields
mind = u;
rmap = rmapa(:,:,:,:,mind);
rmapEX = nan([bhvSpcWidth*4,bhvSpcHeight*4]);
for x = 2:5,
    for y = 2:5,
        if (x==2&&y==2) || (x==2&&y==5)|| (x==5&&y==5)|| (x==5&&y==2)
            continue
        end
        trmap = sq(rmap(x,y,4:25,4:25));
        trmap(~reducedMaskBhv) = nan;
        rmapEX(( bhvSpcWidth*(x-2)+1):(bhvSpcWidth*(x-2) +bhvSpcWidth), ...
               (bhvSpcHeight*(y-2)+1):(bhvSpcHeight*(y-2)+bhvSpcHeight)) = trmap;
    end
end
clim = [min(nonzeros(rmapEX(:))), max(nonzeros(rmapEX(:)))];
pax = pcolor(rmapEX');
pax.EdgeColor='none';
colormap(gca(),'jet');
cax = colorbar();
caxis(clim);
axis('xy');
Lines([22:22:22*4]+0.5,[],'k');
Lines([],[22:22:22*4]+0.5,'k');
xlim([0,22*4+0.5]);
ylim([0,22*4+0.5]);
clear_axes_labels(sax(end));
colorbar(sp(2));
colormap('jet');

% Place fields
for s = 1:6
    subplot(2,6,s+6);
    plot(pfs{cluSessionMap(u,1)}{s},cluSessionMap(u,2),1,'text');
end     

% Change line positions in histogram
axes(sp(1));
    lax = Lines(rhoDR(u),[],'r');
    tax = text(rhoDR(u),45,num2str(u));

axes(sp(3));
    pax = plot(rhoDRMean(u),rhoDR(u),'.r','MarkerSize',20);
    
axes(sp(4));
    bar(tbin,accgHz(u,:));
    axis('tight');
    xlim([-30,30]);
    
axes(sp(5));
        bar(linspace(0,360,36),tphRUN(u,:));
    
u = figure_controls(hfigg,u,sind);
end

134
137
162
125
129
126
32
194
25

189
117
120
115

%%%>>>



%%%<<< SUPFIG THETA phase preference vs Factor Scores ----------------------------------------------

[hfig,fig,fax,sax] = set_figure_layout(figure(666002),'A4','portrait',[],3,3,2,2);
axOpts = {'Units',                 'centimeters',...
          'FontSize',              8,            ...
          'LineWidth',             1,            ...
          'PlotBoxAspectRatioMode','manual'};

for f = 1:3;
% PLOT factors
    [yind, yOffSet, xind, xOffSet] = deal(1, 2, f, 2-f);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height]);
        eigVec = nan(binDims(end-1:end));
        eigVec(validDims) = LU(:,f);
        pax = pcolor(bins{end-1:end},eigVec');
        pax.EdgeColor = 'none';
        sax(end).XTick = [];
        sax(end).YTick = [];        
        axis('xy');
        colormap(hfig.CurrentAxes,'parula');                
        ylabel(colorbar(),'A.U.');
        title(num2str(f,'Factor %d'));
        xlabel('Head Pitch');
        ylabel('Body Pitch');
        daspect([1,1,1]);
% PLOT thetaPhasePreference vs Factor Score (STRONGLY theta modulated)
    [yind, yOffSet, xind, xOffSet] = deal(2, 2, f, 1);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height]);
        ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN>0.2;%& vrCorr>-0.08;
        %scatter(circ_rad2ang(tppRUN(ind,1)+2*pi*double(tppRUN(ind,1)<0)),FSr(ind,2),20,vrCorr(ind),'filled');
        scatter(circ_rad2ang(tppRUN(ind,1)+2*pi*double(tppRUN(ind,1)<0)),FSr(ind,f),15,vrCorrLog(ind),'filled');
        ylim([-3,3]);
        xlim([0,360]);
        hfig.CurrentAxes.XTick = [0,90,180,270,360];        
        xlabel({'mean \theta phase'});
        ylabel(['Factor ',num2str(f),' score']);
        colormap(hfig.CurrentAxes,'jet');                
        if f == 3
            cax = colorbar();
            ylabel(cax,'speed-rate corr');
            cax.Position(1) = cax.Position(1) + .001;
            drawnow();            
            cax.Position(1) = sum(hfig.CurrentAxes.Position([1,3]))+0.01;
        end
        caxis([-0.4,0.6]);
        [rhoCL,pval] = circ_corrcl(tppRUN(ind,1),FSr(ind,f));
        title({num2str(rhoCL,'Rho: %0.4g'),...
               num2str(pval,'P-Val: %10.4e')});
        grid('on');
% PLOT thetaPhasePreference vs Factor Score (WEAKLY theta modulated)
    [yind, yOffSet, xind, xOffSet] = deal(4, 1, f, 1);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height]);
        ind = ismember(cluSessionMap(:,1),[3,4,5,8:12,17:25]) & tprRUN<=0.2;%& vrCorr>-0.08;
        %scatter(circ_rad2ang(tppRUN(ind,1)+2*pi*double(tppRUN(ind,1)<0)),FSr(ind,2),20,vrCorr(ind),'filled');
        scatter(circ_rad2ang(tppRUN(ind,1)+2*pi*double(tppRUN(ind,1)<0)),FSr(ind,f),15,vrCorrLog(ind),'filled');
        ylim([-3,3]);
        xlim([0,360]);
        hfig.CurrentAxes.XTick = [0,90,180,270,360];
        xlabel({'mean \theta phase'});
        ylabel(['Factor ',num2str(f),' score']);
        colormap(hfig.CurrentAxes,'jet');        
        if f == 3
            cax = colorbar();
            ylabel(cax,'speed-rate corr');
            cax.Position(1) = cax.Position(1) + .001;
            drawnow();
            cax.Position(1) = sum(hfig.CurrentAxes.Position([1,3]))+0.01;
        end
        caxis([-0.4,0.6]);
        [rhoCL,pval] = circ_corrcl(tppRUN(ind,1),FSr(ind,f));
        title({num2str(rhoCL,'Rho: %0.4g'),...
               num2str(pval,'P-Val: %10.4e')});
        grid('on');
% PLOT thetaPhaseHistogram foreach Interneuron sorted by factor score (STRONGLY theta modulated)
    [yind, yOffSet, xind, xOffSet] = deal(3, 2, f, 1);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height]);
        [~,sind] = sort(FSr(:,f));
        stphRUN = tphRUN(sind,:);
        stprRUN = tprRUN(sind);
        scluSessionMap = cluSessionMap(sind,:);
        ind = ismember(scluSessionMap(:,1),[3,4,5,8:12,17:25]) & stprRUN>0.2;%& vrCorr>-0.08;
        imagesc(linspace(0,360,36),1:sum(ind),bsxfun(@rdivide,stphRUN(ind,:),max(stphRUN(ind,:),[],2)));
        axis('xy');
        hfig.CurrentAxes.XTick = [0,90,180,270,360];        
        ylabel('Unit');
        xlabel('Normalize Count');
        title({'\theta Phase Histogram',['sorted by factor ',num2str(f),' score']});
        colormap(hfig.CurrentAxes,'parula')        
% PLOT thetaPhaseHistogram foreach Interneuron sorted by factor score (WEAKLY theta modulated)
    [yind, yOffSet, xind, xOffSet] = deal(5, 1, f, 1);
    sax(end+1) = axes(axOpts{:},                                            ...
                      'Position',[fig.page.xpos(xind)+xOffSet,              ...
                                  fig.page.ypos(yind)+yOffSet,              ...
                                  fig.subplot.width,                        ...
                                  fig.subplot.height]);
        [~,sind] = sort(FSr(:,f));
        stphRUN = tphRUN(sind,:);
        stprRUN = tprRUN(sind);
        scluSessionMap = cluSessionMap(sind,:);
        ind = ismember(scluSessionMap(:,1),[3,4,5,8:12,17:25]) & stprRUN<0.2;%& vrCorr>-0.08;
        imagesc(linspace(0,360,36),1:sum(ind),bsxfun(@rdivide,stphRUN(ind,:),max(stphRUN(ind,:),[],2)));        
        axis('xy');
        hfig.CurrentAxes.XTick = [0,90,180,270,360];
        ylabel('Unit');
        xlabel('Normalize Count');
        title({'\theta Phase Histogram',['Sorted by factor ',num2str(f),' score']})
        colormap(hfig.CurrentAxes,'parula')
end
[yind, yOffSet, xind, xOffSet] = deal(1, 2, 4, -1.5);
sax(end+1) = axes(axOpts{:},                                            ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width,                        ...
                              fig.subplot.height]);
plot(VT(1:5,4),'+-');
xlim([0,6]);
xlabel('Factors');
ylabel('Exp Var');
sax(end).Position([2,3]) = sax(1).Position([2,3])+[0.5,0];
sax(end).Position([4]) = 2;

%hfig.CurrentAxes = fax;

text(fax,2,26.5,'A','FontSize',18)
text(fax,2,22.5,'B','FontSize',18)
text(fax,2,11.5,'C','FontSize',18)



preserve_figure_layout(hfig);
set_figure_layout_mode(hfig,'absolute');
set_figure_layout_mode(hfig,'normalized');
set_figure_layout_mode(hfig,'fixedaspect');
reset_figure_layout(hfig);

%%%>>>



%%%<<< SUPFIG Spatial Vs Behavioral Information ----------------------------------------------------

[hfig,fig,fax,sax] = set_figure_layout(figure(666003),'A4','portrait',[],2,2,2,2);
axOpts = {'Units',                 'centimeters',...
          'FontSize',              8,            ...
          'LineWidth',             1,            ...
          'PlotBoxAspectRatioMode','manual'};

% SPATIAL and BEHAVIORAL information
[yind, yOffSet, xind, xOffSet] = deal(1, 0, 1, 0);
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                          ...
                              fig.page.ypos(yind)+yOffSet,                          ...
                              fig.subplot.width*2,    ...
                              fig.subplot.height*2],    ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
scatter(log10(sinfNP),log10(sinfNB),10,mrateNT,'Filled');
colormap('jet');
colorbar();
title({'Spatial Vs Behavioral','Information'});
xlabel({'Spatial', 'log10(bits/spk)'});
ylabel({'Behavioral', 'log10(bits/spk)'});
daspect([1,1,1]);
xlim([-4,0]);
ylim([-4,0]);
line([-4,0],[-4,0]);


[yind, yOffSet, xind, xOffSet] = deal(1, 0, 3,0);
sax(end+1) = axes('Units','centimeters',                                            ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                          ...
                              fig.page.ypos(yind)+yOffSet,                          ...
                              fig.subplot.width*2,    ...
                              fig.subplot.height*2],    ...
                  'FontSize', 8,                                                    ...
                  'LineWidth',1);
hold('on');
scatter(zsinfNP,zsinfNB,10,mrateNT,'Filled');
colormap('jet');
colorbar();
title({'Spatial Vs Behavioral','Information'});
xlabel({'Spatial','z-score'});
ylabel({'Behaviorial', 'z-score'});
line([-20,100],[-20,100]);
daspect([1,1,1])
circle(0,0,3.1)
xlim([-10,100]);
ylim([-10,100]);

%%%>>>



%%%<<< supplemental panel - inter bhv-map(pos-mapCorre) by distance

xOffSet = -1;
for  e = 1:size(eupInd,1),
    [yind, yOffSet, xind, xOffSet] = deal(sectionYind+3, -1, e*2-1, xOffSet+1);
    sax(end+1) = axes(   'Units', 'centimeters',                                         ...
                         'Position', [fig.page.xpos(xind)+xOffSet,                          ...
                        fig.page.ypos(yind)+yOffSet,                          ...
                        fig.subplot.width*2+fig.subplot.horizontalPadding,    ...
                        fig.subplot.height],    ...
                         'FontSize',  8,                                                    ...
                         'LineWidth',  1);
    hold('on');
    trlUntPair = [expUnitsPP{eupInd(e,1),1},expUnitsPP{eupInd(e,1),2}];    
    mind = find(ismember(cluSessionMap,trlUntPair,'rows'));    
    % Change dist from bins to mm
    plot(dst(:),rho(mind,:),'.')
    ylim([-1,1]);
end

%%%>>>



%%%<<< TODO compute the theta phase preference for each interneuron



rFSr = reshape(FSr,16,255,5);
mrateNRdT = mean(rmapNR-rmapNT,'omitnan')';

figure();
hold('on');
plot([tppRUN;tppRUN+2*pi],[mrateNB;mrateNB],'.');
plot([tppRUN;tppRUN+2*pi],[mrateNR;mrateNR],'.r');

figure();
hold('on');
plot([tppRUN;tppRUN+2*pi],[mrateNB;mrateNB],'.');
plot([tppRUN;tppRUN+2*pi],[mrateNR;mrateNR],'.r');


ind = tprRUN>0.20;
thetaVec = linspace(-2*pi,4*pi,200);
figure();
hold('on');
subplot(321);
    hold('on');    
    scatter([tppRUN(~ind);tppRUN(~ind)+2*pi],...
            [mrateNRdT(~ind);mrateNRdT(~ind)],...%[mrateNH(~ind)-mrateNL(~ind);mrateNH(~ind)-mrateNL(~ind)],...
            10,...
            [mrateNT(~ind);mrateNT(~ind)],...
            'Filled');
% $$$     scatter([tppRUN(~ind);tppRUN(~ind)+2*pi],...
% $$$             [mrateNRdT(~ind);mrateNRdT(~ind)],...%[mrateNH(~ind)-mrateNL(~ind);mrateNH(~ind)-mrateNL(~ind)],...
% $$$             10,...
% $$$             [mrateNT(~ind);mrateNT(~ind)],...
% $$$             'Filled');
    plot(thetaVec,3*cos(thetaVec),'c','LineWidth',2);
    grid('on');
    xlim([-pi,3*pi]);
    ylim([-7,7]);
    ylabel ({'Rate difference between','theta and rear state'});
    xlabel ('Theta Phase')
subplot(322);
    hold('on');
    scatter([tppRUN(ind);tppRUN(ind)+2*pi],...
            [mrateNRdT(ind);mrateNRdT(ind)],...%[mrateNH(ind)-mrateNL(ind);mrateNH(ind)-mrateNL(ind)],...
            10,...
            [mrateNT(ind);mrateNT(ind)],...
            'Filled');
    plot(thetaVec,3*cos(thetaVec),'c','LineWidth',2);
    grid('on');    
    xlim([-pi,3*pi]);    
    ylim([-7,7]);
    xlabel ('Theta Phase')
subplot(323);
    hold('on');    
    scatter([tppRUN(~ind);tppRUN(~ind)+2*pi],...
            [zsinfNB(~ind);zsinfNB(~ind)],...%[mrateNH(~ind)-mrateNL(~ind);mrateNH(~ind)-mrateNL(~ind)],...
            10,...
            [mrateNT(~ind);mrateNT(~ind)],...
            'Filled');
% $$$     scatter([tppRUN(~ind);tppRUN(~ind)+2*pi],...
% $$$             [mrateNRdT(~ind);mrateNRdT(~ind)],...%[mrateNH(~ind)-mrateNL(~ind);mrateNH(~ind)-mrateNL(~ind)],...
% $$$             10,...
% $$$             [mrateNT(~ind);mrateNT(~ind)],...
% $$$             'Filled');
    plot(thetaVec,3*cos(thetaVec),'c','LineWidth',2);
    grid('on');
    xlim([-pi,3*pi]);
    %    ylim([-7,7]);
    ylabel ({'bhv zscore'});
    xlabel ('Theta Phase')
subplot(324);
    hold('on');
    scatter([tppRUN(ind);tppRUN(ind)+2*pi],...
            [zsinfNB(ind);zsinfNB(ind)],...%[mrateNH(ind)-mrateNL(ind);mrateNH(ind)-mrateNL(ind)],...
            10,...
            [mrateNT(ind);mrateNT(ind)],...
            'Filled');
    plot(thetaVec,3*cos(thetaVec),'c','LineWidth',2);
    grid('on');    
    xlim([-pi,3*pi]);    
    %    ylim([-7,7]);
    xlabel ('Theta Phase')
subplot(325);
    hold('on');    
    scatter([tppRUN(~ind);tppRUN(~ind)+2*pi],...
            [zsinfNP(~ind);zsinfNP(~ind)],...%[mrateNH(~ind)-mrateNL(~ind);mrateNH(~ind)-mrateNL(~ind)],...
            10,...
            [mrateNT(~ind);mrateNT(~ind)],...
            'Filled');
% $$$     scatter([tppRUN(~ind);tppRUN(~ind)+2*pi],...
% $$$             [mrateNRdT(~ind);mrateNRdT(~ind)],...%[mrateNH(~ind)-mrateNL(~ind);mrateNH(~ind)-mrateNL(~ind)],...
% $$$             10,...
% $$$             [mrateNT(~ind);mrateNT(~ind)],...
% $$$             'Filled');
    plot(thetaVec,3*cos(thetaVec),'c','LineWidth',2);
    grid('on');
    xlim([-pi,3*pi]);
    %    ylim([-7,7]);
    ylabel ({'pos zscore'});
    xlabel ('Theta Phase')
subplot(326);
    hold('on');
    scatter([tppRUN(ind);tppRUN(ind)+2*pi],...
            [zsinfNP(ind);zsinfNP(ind)],...%[mrateNH(ind)-mrateNL(ind);mrateNH(ind)-mrateNL(ind)],...
            10,...
            [mrateNT(ind);mrateNT(ind)],...
            'Filled');
    plot(thetaVec,3*cos(thetaVec),'c','LineWidth',2);
    grid('on');    
    xlim([-pi,3*pi]);    
    %    ylim([-7,7]);
    xlabel ('Theta Phase')
colormap('jet');
cax = colorbar();


figure();
plot([tppRUN;tppRUN+2*pi],[tppREM;tppREM+2*pi],'.')

% TODO scatter bhv x pos information colored by mean rate
figure,
scatter(log10(sinfNP),log10(sinfNB),15,mrateNT,'Filled');
colormap('jet');
colorbar();

% TODO bhv and pos information CDF
figure,
subplot(221);    plot(mrateNT, zsinfNP,'.');
subplot(222);    plot(mrateNB, zsinfNB,'.');
subplot(223);    plot(mrateNT, log10(sinfNP),'.');
subplot(224);    plot(mrateNB, log10(sinfNB),'.');

% TODO bin bhv x pos information z-score along the slope or by the angle


%spk gamma phz lock pyr and int

% body position given head kinematics
% Given the vector head_com to head_front, what is the head_com to spine_upper vector
% convert frame of reference from room to head for all markers 

% x_body_now = [s,t]*x_head_embedded +a*INTEGRAL(HAV)+b*INTEGRAL(HVL)+c*INTEGRAL(HVF)


xyz = preproc_xyz(Trials{20},'trb');;
ang = create(MTADang(

href = fet_href(Trials{20},[],[],'trb');


ind = Trials{20}.stc{'n'};
mind = Trials{20}.stc{'m'};
wind = Trials{20}.stc{'w'};
hfig = figure();
for m = 1:5,
    sax = subplot2(5,1,m,1);
    hold(sax,'on');
    plot(href(:,m),href(:,m+5),'.');
    plot(href(mind,m),href(mind,m+5),'.m');
    plot(href(ind,m),href(ind,m+5),'.g');    
    plot(href(wind,m),href(wind,m+5),'.k');    
end


% A. Examples
% B.

% TODO COMPUTE the phase resolved 

%%%>>>



%%%<<< TESTING  ------------------------------------------------------------------------------------

ind = tprRUN>0.25;
thetaVec = linspace(-2*pi,4*pi,200);
[yind, yOffSet, xind, xOffSet] = deal(sectionYind+2, -2, 5, 1);
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                    fig.page.ypos(yind)+yOffSet,              ...
                    fig.subplot.width.*3+fig.subplot.horizontalPadding*2,...
                    fig.subplot.height.*1.5],  ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
scatter([tppRUN(~ind)';tppRUN(~ind)'+2*pi],...
        [mrateNRdT(~ind);mrateNRdT(~ind)],...%[mrateNH(~ind)-mrateNL(~ind);mrateNH(~ind)-mrateNL(~ind)],...
        10,...
        [mrateNT(~ind);mrateNT(~ind)],...
        'Filled');
plot(thetaVec,3*cos(thetaVec),'Color',[0.75,0.75,0.75],'LineWidth',2);
grid('on');
xlim([-pi,3*pi]);
ylim([-7,7]);
%ylabel ({'Rate Difference','theta-rear'});
xlabel ('Theta Phase')

[yind, yOffSet, xind, xOffSet] = deal(sectionYind+4, -2, 5, 1);
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width.*3+fig.subplot.horizontalPadding*2,...
                              fig.subplot.height.*1.5],  ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(thetaVec,3*cos(thetaVec),'Color',[0.75,0.75,0.75],'LineWidth',2);
scatter([tppRUN(ind)';tppRUN(ind)'+2*pi],...
        [mrateNRdT(ind);mrateNRdT(ind)],...%[mrateNH(ind)-mrateNL(ind);mrateNH(ind)-mrateNL(ind)],...
        10,...
        [mrateNT(ind);mrateNT(ind)],...
        'Filled');
grid('on');    
xlim([-pi,3*pi]);    
ylim([-7,7]);
xlabel ('Theta Phase')

figure()
hold('on');
grid('on');
ind = ~~mrateNT;
plot(thetaVec,3*cos(thetaVec),'Color',[0.75,0.75,0.75],'LineWidth',2);
scatter([tppRUN(ind)';tppRUN(ind)'+2*pi],...
        [mrateNRdT(ind);mrateNRdT(ind)],...%[mrateNH(ind)-mrateNL(ind);mrateNH(ind)-mrateNL(ind)],...
        10,...
        [tprRUN(ind)';tprRUN(ind)'],...
        'Filled');
colormap('jet');
colorbar();
xlim([-pi,3*pi]);    
ylim([-7,7]);

%%%>>>



%%%<<< OLD parts -----------------------------------------------------------------------------------

%%%<<< Panel - ThetaPhasePreference vs factor score (or) rear/non-rear state rate diff

ind = tprRUN>0.25;
thetaVec = linspace(-2*pi,4*pi,200);
[yind, yOffSet, xind, xOffSet] = deal(sectionYind+2, -2, 5, 1);
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                    fig.page.ypos(yind)+yOffSet,              ...
                    fig.subplot.width.*3+fig.subplot.horizontalPadding*2,...
                    fig.subplot.height.*1.5],  ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
scatter([tppRUN(~ind);tppRUN(~ind)+2*pi],...
        [mrateNRdT(~ind);mrateNRdT(~ind)],...%[mrateNH(~ind)-mrateNL(~ind);mrateNH(~ind)-mrateNL(~ind)],...
        10,...
        [mrateNT(~ind);mrateNT(~ind)],...
        'Filled');
plot(thetaVec,3*cos(thetaVec),'Color',[0.75,0.75,0.75],'LineWidth',2);
grid('on');
xlim([-pi,3*pi]);
ylim([-7,7]);
%ylabel ({'Rate Difference','theta-rear'});
xlabel ('Theta Phase')

[yind, yOffSet, xind, xOffSet] = deal(sectionYind+4, -2, 5, 1);
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width.*3+fig.subplot.horizontalPadding*2,...
                              fig.subplot.height.*1.5],  ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(thetaVec,3*cos(thetaVec),'Color',[0.75,0.75,0.75],'LineWidth',2);
scatter([tppRUN(ind);tppRUN(ind)+2*pi],...
        [mrateNRdT(ind);mrateNRdT(ind)],...%[mrateNH(ind)-mrateNL(ind);mrateNH(ind)-mrateNL(ind)],...
        10,...
        [mrateNT(ind);mrateNT(ind)],...
        'Filled');
grid('on');    
xlim([-pi,3*pi]);    
ylim([-7,7]);
xlabel ('Theta Phase')

%%%>>>

%%%>>>



