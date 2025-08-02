% req20210705
% 
% XYHB 30sps 300ms

global MTA_PROJECT_PATH

configure_default_args();

MjgER2016_load_data();

sampleRate = 30;

%%%<<< XYHB Decoding -------------------------------------------------------------------------------

states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};

tind = [3:5,17:25];
cind = 1;
dc = {};
for t = tind;
    Trial = Trials{t};
    tunits = units{t};
    thetaRefChan = sessionList(t).thetaRefGeneral;

    halfSpkWindow = 0.15;
    ufrWindow = 0.15;

    smoothingWeights = [800.^2,800.^2, 1.2.^2, 1.2.^2];
    %smoothingWeights = [250.^2,250.^2, 0.4.^2, 0.4.^2];

    xyz = preproc_xyz(Trial,'trb');
    xyz.resample(sampleRate);

    pfs = compute_xyhb_ratemaps( Trial, tunits);

    ds = load(fullfile(MTA_PROJECT_PATH,'analysis','pfsXYHB_mask.mat'));
    mask = ds.mask;

    spk = Trial.load('spk', sampleRate, '', tunits, 'deburst');
    ufr = Trial.load('ufr', xyz,spk,tunits,ufrWindow,'boxcar',true);

    dc{cind} = decode_ufr_boxcar(Trial, ...
                                 tunits, ...
                                 sampleRate, ...
                                 ufr, ...
                                 pfs, ...
                                 mask,...
                                 halfSpkWindow,...
                                 smoothingWeights,...
                                 'overwrite',false);

    dc{cind}.pfs = pfs;

    dc{cind}.xyz  = copy(xyz);
    dc{cind}.xyz.data  = xyz.data(dc{cind}.ind,:,:);

    dc{cind}.hRot = headRotation{t};
    dc{cind}.vxy  = filter(copy(xyz),'ButFilter',4,1.5,'low');;
    dc{cind}.vxy  = dc{cind}.vxy.vel('hcom',[1,2]);
    dc{cind}.vxy.data  = dc{cind}.vxy(dc{cind}.ind,1);
    dc{cind}.lvxy = copy(dc{cind}.vxy);
    dc{cind}.lvxy.data(dc{cind}.lvxy.data<=0) = 0.001;
    dc{cind}.lvxy.data = log10(dc{cind}.lvxy.data);

    dc{cind}.hvec = dc{cind}.xyz(:,'nose',[1,2])-dc{cind}.xyz(:,'hcom',[1,2]);
    dc{cind}.hvec = sq(bsxfun(@rdivide,dc{cind}.hvec,sqrt(sum(dc{cind}.hvec.^2,3))));
    dc{cind}.hvec = cat(3,dc{cind}.hvec,sq(dc{cind}.hvec)*[0,-1;1,0]);
    dc{cind}.hvec = multiprod(dc{cind}.hvec,                           ...
                              [cos(dc{cind}.hRot),-sin(dc{cind}.hRot); ...
                        sin(dc{cind}.hRot), cos(dc{cind}.hRot)],...
                              [2,3],                                   ...
                              [1,2]);

    dc{cind}.tvec = circshift(xyz(:,'hcom',[1,2]),-round(250*0.1)) - ...
        circshift(xyz(:,'hcom',[1,2]),round(250*0.1));
    dc{cind}.tvec = sq(bsxfun(@rdivide,dc{cind}.tvec,sqrt(sum(dc{cind}.tvec.^2,3))));
    dc{cind}.tvec = cat(3,dc{cind}.tvec,sq(dc{cind}.tvec)*[0,-1;1,0]);
    dc{cind}.tvec = dc{cind}.tvec(dc{cind}.ind,:,:);

    dc{cind}.fet  = fet_HB_pitchB(Trial,dc{cind}.sampleRate);
    dc{cind}.fet.data  = dc{cind}.fet.data(dc{cind}.ind,:);

    dc{cind}.hdist = filter(copy(xyz),'ButFilter',3,5,'low');
    xycoor =    xyz(:,'hcom',[1,2]);
    [~,dc{cind}.hdist.data] = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    dc{cind}.hdist.data = dc{cind}.hdist.data(dc{cind}.ind);

    % COM 
    dc{cind}.ecom = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.com},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.ecom = dc{cind}.ecom{1};
    % SAX 
    dc{cind}.esax = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.sax},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.esax = dc{cind}.esax{1};
    % MAX 
    dc{cind}.emax = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.max},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.emax = dc{cind}.emax{1};
    % LOM 
    dc{cind}.elom = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.lom},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.elom = dc{cind}.elom{1};
    % LAX 
    dc{cind}.elax = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.lax},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.elax = dc{cind}.elax{1};

    % HVANG
    dc{cind}.hvang = filter(copy(xyz),'ButFilter',4,2,'low');
    xycoor = cat(2,...
                 dc{cind}.hvang(:,'spine_upper',[1,2])-dc{cind}.hvang(:,'bcom',[1,2]),...
                 dc{cind}.hvang(:,'nose',[1,2])-dc{cind}.hvang(:,'hcom',[1,2]));
    dc{cind}.hvang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    % Positive: CCW (Left)     Negative: CW (Right)
    dc{cind}.hvang.data = circ_dist(circshift(dc{cind}.hvang.data(:,2),-10),...
                                    circshift(dc{cind}.hvang.data(:,2),+10));
    dc{cind}.hvang = dc{cind}.hvang.data(dc{cind}.ind);

    
    
% COMPUTE corrected hbang
    dc{cind}.hbang = filter(copy(xyz),'ButFilter',3,9,'low');    
    xycoor = cat(2,...
                 dc{cind}.hbang(:,'spine_upper',[1,2])-dc{cind}.hbang(:,'bcom',[1,2]),...
                 dc{cind}.hbang(:,'nose',[1,2])-dc{cind}.hbang(:,'hcom',[1,2]));
    dc{cind}.hbang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    dc{cind}.hbang.data = circ_dist(dc{cind}.hbang.data(:,2),dc{cind}.hbang.data(:,1));
    dc{cind}.hbang = dc{cind}.hbang.data(dc{cind}.ind)+0.2+-0.4*double(cind>4);

% GENERATE state collection matrix    
    dc{cind}.stcm = stc2mat(Trial.stc,xyz,states);
    dc{cind}.stcm = dc{cind}.stcm(dc{cind}.ind,:);

% COMPUTE head and body relative velocity vectors
    hvfl = fet_href_HXY(Trial,sampleRate,[],'trb');
    bvfl = fet_bref_BXY(Trial,sampleRate,[],'trb');
    dc{cind}.hvfl = hvfl(dc{cind}.ind,:);
    dc{cind}.bvfl = bvfl(dc{cind}.ind,:);    
    
    cind  = cind +1;
end                     

ind = dc{cind}.stcm(:,1) == 1 & (dc{cind}.stcm(:,2) == 2 |dc{cind}.stcm(:,3) == 3 | ...
                                 dc{cind}.stcm(:,4) == 4 |dc{cind}.stcm(:,5) == 5 | ...
                                 dc{cind}.stcm(:,6) == 6);
figure();
hist2([dc{cind}.elax(ind,3),dc{cind}.esax(ind,3)],linspace(-1.25,1.25,30),linspace(-1.25,1.25,30),[],'mud')

figure();
hist2([dc{cind}.elax(ind,4),dc{cind}.esax(ind,4)],linspace(-1.25,1.25,30),linspace(-1.25,1.25,30),[],'mud')

ind = dc{cind}.stcm(:,1) == 1 & (dc{cind}.stcm(:,2) == 2);
ind = dc{cind}.stcm(:,1) == 1 & (dc{cind}.stcm(:,3) == 3);
ind = dc{cind}.stcm(:,1) == 1 & (dc{cind}.stcm(:,5) == 5);

% $$$  |dc{cind}.stcm(:,3) == 3 | ...
% $$$                                  dc{cind}.stcm(:,4) == 4 |dc{cind}.stcm(:,5) == 5 | ...
% $$$                                  dc{cind}.stcm(:,6) == 6
figure();
subplot(221)
histogram([dc{cind}.esax(ind,3)],linspace(-1.25,1.25,30));
ylim([0,2000]);
subplot(222)
histogram([dc{cind}.elax(ind,3)],linspace(-1.25,1.25,30));
ylim([0,2000]);
subplot(223)
histogram([dc{cind}.esax(ind,4)],linspace(-1.25,1.25,30));
ylim([0,5000]);
subplot(224)
histogram([dc{cind}.elax(ind,4)],linspace(-1.25,1.25,30));
ylim([0,5000]);

std(dc{cind}.esax(ind,4))
std(dc{cind}.elax(ind,4))
mean(dc{cind}.esax(ind,4))
mean(dc{cind}.elax(ind,4))

hberrSax = sqrt(sum(dc{cind}.esax(:,3:4).^2,2));
hberrLax = sqrt(sum(dc{cind}.elax(:,3:4).^2,2));

[mean(hberrSax(ind)),std(hberrSax(ind))]
[mean(hberrLax(ind)),std(hberrLax(ind))]


%%%>>>



%%%<<< CONVERT cell to array -----------------------------------------------------------------------
cind = 1;
dca = dc{cind};
dca.sessionId = tind(cind).*ones(size(dc{cind}.ind));

%dca.xyz = dca.xyz.data;
dca.pos = sq(dca.xyz(:,'hcom',:));
dca.fet = dca.fet.data;
dca.vxy = dca.vxy.data;
dca.lvxy = dca.lvxy.data;
dca.hdist = dca.hdist.data;

dca.uinc = dc{1}.uinc;
maxUinc = max(cell2array(cf(@(d) size(d.uinc,2),dc)));
dca.uinc = cat(2,dc{1}.uinc,nan(abs([[0,maxUinc]-size(dc{1}.uinc)])));

for cind = 2:numel(dc),
    dca.sessionId = cat(1,dca.sessionId,tind(cind).*ones(size(dc{cind}.ind)));
    dca.ind = cat(1,dca.ind,dc{cind}.ind);
    dca.max = cat(1,dca.max,dc{cind}.max);
    dca.com = cat(1,dca.com,dc{cind}.com);    
    dca.sax = cat(1,dca.sax,dc{cind}.sax);        
    dca.lom = cat(1,dca.com,dc{cind}.lom);    
    dca.lax = cat(1,dca.sax,dc{cind}.lax);        
    dca.post = cat(1,dca.post,dc{cind}.post);
    dca.ucnt = cat(1,dca.ucnt,dc{cind}.ucnt);    
    dca.uinc = cat(1,dca.uinc,cat(2,dc{cind}.uinc,nan(abs([[0,maxUinc]-size(dc{cind}.uinc)]))));    
    %dca.xyz = cat(1,dca.xyz,dc{cind}.xyz.data);    
    dca.pos = cat(1,dca.pos,sq(dc{cind}.xyz(:,'hcom',:)));
    dca.vxy = cat(1,dca.vxy,dc{cind}.vxy.data);        
    dca.lvxy = cat(1,dca.lvxy,dc{cind}.lvxy.data);            

    dca.hvec = cat(1,dca.hvec,dc{cind}.hvec);    
    dca.tvec = cat(1,dca.tvec,dc{cind}.tvec);        
    dca.fet = cat(1,dca.fet,dc{cind}.fet.data);

    dca.hvang = cat(1,dca.hvang,dc{cind}.hvang);    
    dca.hbang = cat(1,dca.hbang,dc{cind}.hbang);        
    dca.hdist = cat(1,dca.hdist,dc{cind}.hdist.data);        
    
    dca.ecom  = cat(1,dca.ecom,dc{cind}.ecom);        
    dca.esax  = cat(1,dca.esax,dc{cind}.esax);        
    dca.emax  = cat(1,dca.emax,dc{cind}.emax);

    dca.elom  = cat(1,dca.elom,dc{cind}.elom);        
    dca.elax  = cat(1,dca.elax,dc{cind}.elax);        

    dca.stcm = cat(1,dca.stcm,dc{cind}.stcm);    
    dca.hvfl = cat(1,dca.hvfl,dc{cind}.hvfl);        
    dca.bvfl = cat(1,dca.bvfl,dc{cind}.bvfl);            
    
    dca.hRot  = cat(1,dca.hRot,dc{cind}.hRot);
    dca.pfs   = cat(1,dca.pfs,dc{cind}.pfs);
end


dca.posi = discretize(dca.pos,linspace(-500,500,21));
dca.feti(:,1) = discretize(dca.fet(:,1),linspace(-2,0.8,29));
dca.feti(:,2) = discretize(dca.fet(:,2),linspace(-0.8,2,29));


%%%>>>


%%%<<< Compare COM vs LOM  
ind =   nniz(dca.feti(:,1:2)) ...
      & dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2|dca.stcm(:,3)==3|dca.stcm(:,4)==4|dca.stcm(:,5)==5|dca.stcm(:,6)==6) ...
      & dca.ucnt >= 2;


figure,
subplot(221);
hist2([dca.ecom(ind,2),dca.ecom(ind,1)], linspace(-300,300,100), linspace(-300,300,100));
subplot(222);
hist2([dca.ecom(ind,3),dca.ecom(ind,4)], linspace(-1.5,1.5,100), linspace(-1.5,1.5,100));
subplot(223);
hist2([dca.elax(ind,2),dca.elax(ind,1)], linspace(-300,300,100), linspace(-300,300,100));
subplot(224);
hist2([dca.elax(ind,3),dca.elax(ind,4)], linspace(-1.5,1.5,100), linspace(-1.5,1.5,100));

figure();
subplot(221);
histogram(dca.elax(ind,3), linspace(-1.5,1.5,100));
Lines(0,[],'r');
subplot(223);
histogram(dca.esax(ind,3), linspace(-1.5,1.5,100));
Lines(0,[],'r');
subplot(222);
histogram(dca.elax(ind,4), linspace(-1.5,1.5,100));
Lines(0,[],'r');
subplot(224);
histogram(dca.esax(ind,4), linspace(-1.5,1.5,100));
Lines(0,[],'r');


[std(dca.ecom(ind,1)),std(dca.esax(ind,1)),std(dca.elom(ind,1)),std(dca.elax(ind,1))]
[mean(dca.ecom(ind,1)),mean(dca.esax(ind,1)),mpean(dca.elom(ind,1)),mean(dca.elax(ind,1))]

[std(dca.ecom(ind,2)),std(dca.esax(ind,2)),std(dca.elom(ind,2)),std(dca.elax(ind,2))]
[mean(dca.ecom(ind,2)),mean(dca.esax(ind,2)),mean(dca.elom(ind,2)),mean(dca.elax(ind,2))]

[std(dca.ecom(ind,3)),std(dca.esax(ind,3)),std(dca.elom(ind,3)),std(dca.elax(ind,3))]
[mean(dca.ecom(ind,3)),mean(dca.esax(ind,3)),mean(dca.elom(ind,3)),mean(dca.elax(ind,3))]

[std(dca.ecom(ind,4)),std(dca.esax(ind,4)),std(dca.elom(ind,4)),std(dca.elax(ind,4))]
[mean(dca.ecom(ind,4)),mean(dca.esax(ind,4)),mean(dca.elom(ind,4)),mean(dca.elax(ind,4))]

%%%>>>


%%%<<< PLOT MUD XY vs HB error

bfsn = cf(@(t,u)  compute_bhv_ratemaps(t,u),          Trials, units);
bfss = cf(@(t,u)  compute_bhv_ratemaps_shuffled(t,u), Trials, units);

[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfsn, units, [], [], false);
numComp = size(eigVecs,2);
bhvMask = false(size(validDims));
bhvMask(validDims) = true;
bhvMask = reshape_eigen_vector(bhvMask,bfsn)';

pftn = cf(@(t,u) pfs_2d_theta(t,u), Trials, units);

maskPos = create_tensor_mask(pftn{1}.adata.bins);
posInfoN = {};
for t = 1:numel(bfsn)
for u = 1:numel(units{t})
rmap = pftn{t}.data.rateMap(logical(maskPos(:)),pftn{t}.data.clu==units{t}(u));
MRate = mean(rmap,'omitnan');
posInfoN{t}(u) = sum(1/ngBins.*(rmap./MRate).*log2(rmap./MRate),'omitnan');
end
end


ngBins = sum(validDims);

bhvInfoN = {};
for t = 1:numel(bfsn)
for u = 1:numel(units{t})
rmap = bfsn{t}.data.rateMap(validDims,bfsn{t}.data.clu==units{t}(u));
MRate = mean(rmap,'omitnan');
bhvInfoN{t}(u) = sum(1/ngBins.*(rmap./MRate).*log2(rmap./MRate),'omitnan');
end
end

bhvInfoS = {};
for t = 1:numel(bfsn)
for u = 1:numel(units{t})
rmap = sq(bfss{t}.data.rateMap(validDims,bfsn{t}.data.clu==units{t}(u),:));
MRate = mean(rmap,'omitnan');
bhvInfoS{t}(u,:) = sum(1/ngBins.*(rmap./MRate).*log2(rmap./MRate),'omitnan');
end
end


% $$$ u = 60;
% $$$ (bhvInfoN{20}(66)-mean(bhvInfoS{20}(66,:),2))./std(bhvInfoS{20}(66,:),[],2)

dca.BhvInfoM = nan(size(dca.ind));
dca.BhvInfoS = nan(size(dca.ind));
dca.PosInfoM = nan(size(dca.ind));
dca.PosInfoS = nan(size(dca.ind));
for ind = 1:size(dca.uinc,1)
% $$$     nzd = nonzeros(posInfoN{dca.sessionId(ind)}.*dca.uinc(ind,~isnan(dca.uinc(ind,:))));
% $$$     dca.PosInfoM(ind,1) = mean(nzd);    
% $$$     dca.PosInfoS(ind,1) = std(nzd);        
    nzd = nonzeros(bhvInfoN{dca.sessionId(ind)}.*dca.uinc(ind,~isnan(dca.uinc(ind,:))));    
    dca.BhvInfoM(ind,1) = mean(nzd);
    dca.BhvInfoS(ind,1) = std(nzd);    
end

ind = dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2 | dca.stcm(:,3)==3 | dca.stcm(:,4)==4 | dca.stcm(:,5)==5 | dca.stcm(:,6)==6);

nEdges = 11;
nBins = nEdges-1;
nXTicks = 2;
nYTicks = 2;

[uPosErr,xPosErr] = MakeUniformDistr(sqrt(sum(bsxfun(@plus,dca.elax(ind,1:2), [-35,0]).^2,2))/10);
[uBhvErr,xBhvErr] = MakeUniformDistr(sqrt(sum(dca.elax(ind,3:4).^2,2)));
uPosErrInd = discretize(uPosErr,linspace([xPosErr([1,end])',nEdges]));
uBhvErrInd = discretize(uBhvErr,linspace([xBhvErr([1,end])',nEdges]));

[uPosErr,xPosErr] = MakeUniformDistr(sqrt(sum(bsxfun(@plus,dca.elax(ind,1:2), [-35,0]).^2,2))/10);
[uBhvErr,xBhvErr] = MakeUniformDistr(sqrt(sum(dca.elax(ind,3:4).^2,2)));

uPosErrInd = discretize(uPosErr,linspace([xPosErr([1,end])',31]));
uBhvErrInd = discretize(uBhvErr,linspace([xBhvErr([1,end])',31]));



xticks = round(xPosErr(round([1:size(xPosErr,1)/(nBins/nXTicks):size(xPosErr,1),size(xPosErr,1)])),0);
yticks = round(xBhvErr(round([1:size(xBhvErr,1)/(nBins/nYTicks):size(xBhvErr,1),size(xBhvErr,1)])),2);

% ESAX and ELAX
figure,
[uPosErr,xPosErr] = MakeUniformDistr(sqrt(sum(bsxfun(@plus,dca.esax(ind,1:2), [-35,0]).^2,2))/10);
[uBhvErr,xBhvErr] = MakeUniformDistr(sqrt(sum(dca.esax(ind,3:4).^2,2)));
uPosErrInd = discretize(uPosErr,linspace([xPosErr([1,end])',nEdges]));
uBhvErrInd = discretize(uBhvErr,linspace([xBhvErr([1,end])',nEdges]));
subplot(221);
    uposInfoM = dca.PosInfoM(ind);
    nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(uposInfoM);
    out = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],uposInfoM(nind),[nBins,nBins],@mean);
    imagesc(out')
    set(gca(),'XTick',0:nXTicks:nBins);
    set(gca(),'YTick',0:nYTicks:nBins);
    set(gca(),'XTickLabel',xticks);
    set(gca(),'YTickLabel',yticks);
    axis('xy');
    colormap('jet');
    caxis([0.54,0.81]);    
subplot(222);
    ubhvInfoM = dca.BhvInfoM(ind);
    nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(ubhvInfoM);
    out = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],ubhvInfoM(nind),[nBins,nBins],@mean);
    imagesc(out')
    set(gca(),'XTick',0:nXTicks:nBins);
    set(gca(),'YTick',0:nYTicks:nBins);
    set(gca(),'XTickLabel',xticks);
    set(gca(),'YTickLabel',yticks);
    axis('xy');
    colormap('jet');
    caxis([0.17,0.3])

[uPosErr,xPosErr] = MakeUniformDistr(sqrt(sum(bsxfun(@plus,dca.elax(ind,1:2), [-35,0]).^2,2))/10);
[uBhvErr,xBhvErr] = MakeUniformDistr(sqrt(sum(dca.elax(ind,3:4).^2,2)));
uPosErrInd = discretize(uPosErr,linspace([xPosErr([1,end])',nEdges]));
uBhvErrInd = discretize(uBhvErr,linspace([xBhvErr([1,end])',nEdges]));
subplot(223);
    uposInfoM = dca.PosInfoM(ind);
    nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(uposInfoM);
    outP = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],uposInfoM(nind),[nBins,nBins],@mean);
    imagesc(outP')
    set(gca(),'XTick',0:nXTicks:nBins);
    set(gca(),'YTick',0:nYTicks:nBins);
    set(gca(),'XTickLabel',xticks);
    set(gca(),'YTickLabel',yticks);
    axis('xy');
    colormap('jet');
    caxis([0.54,0.81]);    
subplot(224);
    ubhvInfoM = dca.BhvInfoM(ind);
    nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(ubhvInfoM);
    outB = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],ubhvInfoM(nind),[nBins,nBins],@mean);
    imagesc(outB')
    set(gca(),'XTick',0:nXTicks:nBins);
    set(gca(),'YTick',0:nYTicks:nBins);
    set(gca(),'XTickLabel',xticks);
    set(gca(),'YTickLabel',yticks);
    axis('xy');
    colormap('jet');
    caxis([0.17,0.3])


% ECOM and ELOM
figure,
[uPosErr,xPosErr] = MakeUniformDistr(sqrt(sum(bsxfun(@plus,dca.ecom(ind,1:2), [-35,0]).^2,2))/10);
[uBhvErr,xBhvErr] = MakeUniformDistr(sqrt(sum(dca.ecom(ind,3:4).^2,2)));
uPosErrInd = discretize(uPosErr,linspace([xPosErr([1,end])',nEdges]));
uBhvErrInd = discretize(uBhvErr,linspace([xBhvErr([1,end])',nEdges]));
subplot(221);
    uposInfoM = dca.PosInfoM(ind);
    nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(uposInfoM);
    out = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],uposInfoM(nind),[nBins,nBins],@mean);
    imagesc(out')
    set(gca(),'XTick',0:nXTicks:nBins);
    set(gca(),'YTick',0:nYTicks:nBins);
    set(gca(),'XTickLabel',xticks);
    set(gca(),'YTickLabel',yticks);
    axis('xy');
    colormap('jet');
    caxis([2.75,4.75]);
subplot(222);
    ubhvInfoM = dca.BhvInfoM(ind);
    nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(ubhvInfoM);
    out = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],ubhvInfoM(nind),[nBins,nBins],@mean);
    imagesc(out')
    set(gca(),'XTick',0:nXTicks:nBins);
    set(gca(),'YTick',0:nYTicks:nBins);
    set(gca(),'XTickLabel',xticks);
    set(gca(),'YTickLabel',yticks);
    axis('xy');
    colormap('jet');
    caxis([0.17,0.3])

[uPosErr,xPosErr] = MakeUniformDistr(sqrt(sum(bsxfun(@plus,dca.elom(ind,1:2), [-35,0]).^2,2))/10);
[uBhvErr,xBhvErr] = MakeUniformDistr(sqrt(sum(dca.elom(ind,3:4).^2,2)));
uPosErrInd = discretize(uPosErr,linspace([xPosErr([1,end])',nEdges]));
uBhvErrInd = discretize(uBhvErr,linspace([xBhvErr([1,end])',nEdges]));
subplot(223);
    uposInfoM = dca.PosInfoM(ind);
    nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(uposInfoM);
    outP = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],uposInfoM(nind),[nBins,nBins],@mean);
    imagesc(outP')
    set(gca(),'XTick',0:nXTicks:nBins);
    set(gca(),'YTick',0:nYTicks:nBins);
    set(gca(),'XTickLabel',xticks);
    set(gca(),'YTickLabel',yticks);
    axis('xy');
    colormap('jet');
    caxis([3,4.75]);    
subplot(224);
    ubhvInfoM = dca.BhvInfoM(ind);
    nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(ubhvInfoM);
    outB = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],ubhvInfoM(nind),[nBins,nBins],@mean);
    imagesc(outB')
    set(gca(),'XTick',0:nXTicks:nBins);
    set(gca(),'YTick',0:nYTicks:nBins);
    set(gca(),'XTickLabel',xticks);
    set(gca(),'YTickLabel',yticks);
    axis('xy');
    colormap('jet');
    caxis([0.17,0.3])
    

    
ucntm = dca.ucnt(ind);
nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(ucntm);
out = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],...
                 ucntm(nind),...
                 [nBins,nBins],...
                 @mean);
figure, imagesc(out');axis('xy')    
    set(gca(),'XTick',0:nXTicks:nBins);
    set(gca(),'YTick',0:nYTicks:nBins);
    set(gca(),'XTickLabel',xticks);
    set(gca(),'YTickLabel',yticks);
    axis('xy');
    colormap('jet');
  
    
figure();
    imagesc(sqrt(outP*outB)');
    axis('xy');
    set(gca(),'XTick',0:nXTicks:nBins);
    set(gca(),'YTick',0:nYTicks:nBins);
    set(gca(),'XTickLabel',xticks);
    set(gca(),'YTickLabel',yticks);
    colormap('jet');
    colorbar();



%%%>>>
