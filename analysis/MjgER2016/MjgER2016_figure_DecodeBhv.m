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

% DEF selection of trials by index
tInds = [3:5,17:23];

% DEF behavior field rate map
bfrm        = cf(@(t,u)   compute_bhv_ratemaps(t,u),                    Trials(tInds), units(tInds));
%bfrmShuff   = cf(@(t,u)   compute_bhv_ratemaps_shuffled(t,u),           Trials(tInds), units(tInds));

% COMPUTE bhv ratemap erpPCA
[eigVecs, eigScrs, eigVars, unitSubsets, validDims, zrmMean, zrmStd] = ...
    compute_bhv_ratemaps_erpPCA(bfrm, units(tInds));

% COMPUTE bhv ratemap erpPCA scores
% $$$ [fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS,fsrsMean,fsrsStd,rmapsShuffledMean,rmapsShuffled] = ...
% $$$     compute_bhv_ratemaps_erpPCA_scores(Trials(tInds),units(tInds),bfrm,bfrmShuff,eigVecs,validDims,unitSubsets);

% GENERATE bhv ratemap mask
maskbhv = false(bfrm{1}.adata.binSizes');
maskbhv(validDims) = true;

% GENERATE xyhb placefields
pfsXYHB  = cf(@(t,u)  compute_xyhb_ratemaps(t,u),  Trials(tInds), units(tInds));
% GENERATE xy ratemap mask
maskcirc = create_tensor_mask(pfsXYHB{1}.adata.bins(1:2));

% COMPLETE mask
% $$$ mask = repmat(maskcirc,[1,1,size(maskbhv)]).*repmat(permute(maskbhv,[3,4,1,2]),[size(maskcirc),1,1]);
%save(fullfile(MTA_PROJECT_PATH,'analysis','pfsXYHB_mask.mat'),'mask','-v7.3');
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



% LOAD unit firing rate
% $$$ xyz  =  cf(@(t) resample(t.load('xyz'),sampleRate), Trials(tInds));
% $$$ ufr  =  cf(@(t,x,u) t.load('ufr',x,[],u,spikeWindow,true,'gauss'), ...
% $$$            Trials(tInds), xyz, units(tInds));



decEstCom     = cell([1,numel(tInds)]);
decEstMax     = cell([1,numel(tInds)]);
decEstSax     = cell([1,numel(tInds)]);
posteriorMax  = cell([1,numel(tInds)]);
unitInclusion = cell([1,numel(tInds)]);

decEstComTPD     = cell([1,numel(tInds)]);
decEstMaxTPD     = cell([1,numel(tInds)]);
decEstSaxTPD     = cell([1,numel(tInds)]);
posteriorMaxTPD  = cell([1,numel(tInds)]);
unitInclusionTPD = cell([1,numel(tInds)]);


for t = tInds;
    Trial = Trials{t};
    unitSubset = units{t};
    
    xyz = preproc_xyz(Trial,'trb');
    xyz.resample(sampleRate);
    spk  = create(copy(Trial.spk),Trial,xyz.sampleRate,'',unitSubset,'deburst');
    ufr = Trial.load('ufr',xyz,spk,unitSubset,spikeWindow,true,'gauss');
    unitInclusion{t==tInds} = sum(ufr.data>0.2,2);    
    tag = 'xy_sr30_HighRes';
    %tag = '';
    [decEstCom{t==tInds},decEstMax{t==tInds},decEstSax{t==tInds},posteriorMax{t==tInds}] = ...
        decode_ufr(Trial,                                                                  ...
                   unitSubset,                                                             ...
                   sampleRate,                                                             ...
                   ufr,pfsXYHB{t==tInds},[],mask,smoothingWeights,'tag',tag);
end


t = 20;
Trial = Trials{t};
unitSubset = units{t};
xyz = resample(Trial.load('xyz','trb'),30);
fet = fet_HB_pitchB(Trial,sampleRate);
ts = [1:size(xyz,1)]./sampleRate;
stcm = stc2mat(Trial.load('stc'),xyz,states);


figure,
ind =    posteriorMax{t==tInds} > 0.001 ...
      & unitInclusion{t==tInds} > 3     ...
      & stcm(:,1)==1;
nanmask = ones([size(xyz,1),1]);
nanmask(~ind)=nan;
% X vs X_est
subplot(4,1,1);
    hold('on');
    plot(ts,xyz(:,5,1),'-k','LineWidth',1);
    plot(ts,sq(decEstCom{t==tInds}(:,1)).*nanmask,'-r','LineWidth',1);
% X vs X_est
subplot(4,1,2);
    hold('on');
    plot(ts,xyz(:,5,2),'-k','LineWidth',1);
    plot(ts,decEstCom{t==tInds}(:,2).*nanmask,'-r','LineWidth',1);
subplot(4,1,3);
    hold('on');
    plot(ts,fet(:,1),'-k','LineWidth',1);
    plot(ts,decEstCom{t==tInds}(:,3).*nanmask,'-r','LineWidth',1);
    ylim([-1.5,0.5]);    
subplot(4,1,4);
    hold('on');
    plot(ts, fet(:,2),'-k','LineWidth',1);
    plot(ts, decEstCom{t==tInds}(:,4).*nanmask,'-r','LineWidth',1);
    ylim([-0.5,1.5]);
linkaxes(findobj(gcf(),'Type','Axes'),'x');

xlim([290,345]);
xlim([530,600]);
xlim([1000,1045]);


% $$$ decError = zeros([size(xyz,1),size(decEstSax,2)]);
% $$$ decError(:,[1,2]) = [multiprod(decEstSax(:,[1,2])-sq(xyz(:,'hcom',[1,2])),hvec(:,:,:),2,[2,3])];
% $$$ decError(:,[3,4],:) = [fet(:,1)-decEstSaxTPD(:,3),fet(:,2)-decEstSax(:,4)];
% $$$ 
% $$$ 
% $$$ ind = all(posteriorMax>0.001,2)&stcm(:,1)==1&~any(logical(stcm(:,[7,8,9])),2)&any(logical(stcm(:,[2,3,5])),2);
% $$$ figure,
% $$$ subplot2(1,2,1,1);
% $$$ hist2(decError(ind,[1,2]),linspace(-300,300,50),linspace(-300,300,50));
% $$$ subplot2(1,2,1,2);
% $$$ hist2(decError(ind,[3,4]),linspace(-1.5,1.5,50),linspace(-1.5,1.5,50));
% $$$ 

sampleRateTPD = 250;
decEstComTPD     = cell([1,numel(tInds)]);
decEstMaxTPD     = cell([1,numel(tInds)]);
decEstSaxTPD     = cell([1,numel(tInds)]);
posteriorMaxTPD  = cell([1,numel(tInds)]);
unitInclusionTPD = cell([1,numel(tInds)]);
xyz       = cf(@(t)     resample(preproc_xyz(t,'trb'),sampleRateTPD),      Trials(tInds));



phzBins = -pi:pi/4:pi;

for t = tInds;
    Trials{t}.lfp.filename = [Trials{t}.name,'.lfp'];    
    tn = find(t==tInds);
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


%[yind, yOffSet, xind, xOffSet] = deal(0, 0, 0, 0);




% COMPILE egocentric variables 
fet  = cf(@(t)    fet_HB_pitchB(t,sampleRateTPD),                                    Trials(tInds));
hvec = cf(@(x)    x(:,'head_front',[1,2])-x(:,'head_back',[1,2]),                    xyz          );
hvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                          hvec         );
hvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                         hvec         );
tvec = cf(@(x)    circshift(x(:,'hcom',[1,2]),-1)-circshift(x(:,'hcom',[1,2]),1),    xyz          );
tvec = cf(@(h)    sq(bsxfun(@rdivide,h,sqrt(sum(h.^2,3)))),                          tvec         );
tvec = cf(@(h)    cat(3,h,sq(h)*[0,-1;1,0]),                                         tvec         );
stc  = cf(@(t)    t.load('stc','msnn_ppsvd_raux'),                                   Trials(tInds));
stcm = cf(@(s,x)  stc2mat(s,x,states),                                               stc,   xyz   );


[hfig,figOpts,xpos,ypos,fax] = set_figure_layout(figure(666006),'A4','portrait');



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

ind  = cf(@(p,u,s) all(p>0.0005,2) & sum(double(u>=1),2)>6 & s(:,1)==1 & any(logical(s(:,stid)),2),...
          posteriorMaxTPD,unitInclusionTPD,stcm);

dcTDPfrontal = cell([1,10]);
dcTDPlateral = cell([1,10]);
for tn = 1:10,
    for j = 1:8;
        dcTDPfrontal{tn} = cat(2,dcTDPfrontal{tn},histc(sq(dErr{tn}(ind{tn},1,j)),linspace(-500,500,100)));
        dcTDPlateral{tn} = cat(2,dcTDPlateral{tn},histc(sq(dErr{tn}(ind{tn},2,j)),linspace(-500,500,100)));        
    end
end

figure,
for tn = 1:10,
    subplot2(2,10,1,tn);
    imagesc(linspace(-500,500,250),1:8,repmat(dcTDPfrontal{tn},1,2)')
    axis('xy');
    xlim([-300,300]);    
    subplot2(2,10,2,tn);
    imagesc(linspace(-500,500,250),1:8,repmat(dcTDPlateral{tn},1,2)')
    axis('xy');
    xlim([-300,300]);
end

figure,
subplot(121);
imagesc(linspace(-500,500,250),1:8,repmat(sum(cat(3,dcTDPfrontal{:}),3),1,2)')
axis('xy');
xlim([-300,300]);
subplot(122);
imagesc(linspace(-500,500,250),1:8,repmat(sum(cat(3,dcTDPlateral{:}),3),1,2)')
axis('xy');
xlim([-300,300]);


out = [];
for j = 1:8;
    out = cat(2,out,histc(sq(dErr{1}(ind{1},3,j)),linspace(-2,2,100)));
end
figure,
imagesc(linspace(-2,2,100),1:8,out')
axis('xy');
    
stid = [3,5];
stid = [4,6];
stid = [3,4,5,6];
stid = [2];
stid = [3];
stid = [5];

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
