sstates% MjgER2016 Figure4
%
% HP:= Head Pitch
% BP:= Body Pitch
% PFD:= behavior field restricted to DRZ[-0.5,0.5]
%
% Subplots:
%    A. phase precession examples across states
%        1. auto correlogram of unit
%        2. PFD behavioral space HPxBP
%        3. PFD behavioral space BSxHS
%        4. PFS theta
%    B. Place field examples
%        1. auto correlogram of unit
%        2. PFD behavioral space HPxBP
%        3. PFD behavioral space BSxHS
%        4. PFS theta
%        5. PFS rear & theta
%        6. PFS hloc & theta
%        7. PFS lloc & theta
%        8. PFS hpause & theta
%        9. PFS lpause & theta
%    C. Eigenvectors of PFD decomposition of HPxBP space
%    D. t-SNE mapping of fscores within HPxBP of first 3 eigenvectors
%    E. ECDF of zscores


MjgER2016_load_data();
%  MjgER2016_load_data:
%
%  Variables:
%      Trials
%      units
%      cluSessionMap
%      pitchReferenceTrial
%      FigDir
%      sessionListName
%      sessionList
%      states
%      numStates
%      interpParPfsp
%      interpParDfs
%
%  Functions:
%      reshape_eigen_vector


if ~exist('pfd','var'), ...
        [pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = ...
        req20180123_ver5(Trials,[],'13','loadPFDFlag',false);  
end



% $$$ 
% $$$ numComp = 3;
% $$$ fpc  = cell([1,numComp]);
% $$$ for i = 1:numComp,
% $$$     fpc{i} = nan(size(validDims{1}));
% $$$     fpc{i}(validDims{1}) = eigVec{1}(:,i);
% $$$ end
% $$$ fpcMinMax = [min(cellfun(@min,fpc)),max(cellfun(@max,fpc))];
% $$$ figure();
% $$$ yind = yind + 1;
% $$$ for i = 1:3,
% $$$     subplot(1,3,i)    
% $$$     imagescnan({pfd{1}.adata.bins{:},abs(reshape_eigen_vector(fpc{i},pfd(1,1)))},...
% $$$                fpcMinMax,'linear',false,nanColor,1,1);                % PRINT eigenvectors
% $$$     axis('xy');
% $$$     axis('tight');
% $$$     ylabel(['F',num2str(i)])
% $$$     xlim(pfd{1}.adata.bins{1}([1,end]));
% $$$     xlim([-2,nonzeros(xlim.*[0,1])-0.2])            
% $$$     ylim(pfd{1}.adata.bins{2}([1,end]));
% $$$ end


MjgER2016_load_bhv_erpPCA_scores();
%     rmaps    - matrix[D x S](numeric); rate maps corresponding to the valid eigenvector dims
%     FSrC     - matrix[U x V](Numeric); fscores
%     fsrcz    - matrix[U x V](Numeric); normalized fscores
clear(varNonEss{:},varNonAux{:});
sessionUnitCnt = cellfun(@numel,units);
unitCnt = sum(sessionUnitCnt);
statesLabels = {'theta','rear','H Loc','H Pause','L Loc','L Pause'};
uind = ismember(cluSessionMap,cluSessionMap(unitSubsets{1},:),'rows');
sampleRate = 250;



% ACCUMULATE phase precession stats ----------------------------------------------------------------
overwrite = false;
numIter = 100;
parmDRZall     = zeros([unitCnt,2,numIter+1,numStates]);
phzStatsDRZall = zeros([unitCnt,2,numIter+1,numStates]);
rDRZall        = zeros([unitCnt,1,numIter+1,numStates]);
rhoDRZall      = zeros([unitCnt,1,numIter+1,numStates]);
parmHRZall     = zeros([unitCnt,2,numIter+1,numStates]);
phzStatsHRZall = zeros([unitCnt,2,numIter+1,numStates]);
rHRZall        = zeros([unitCnt,1,numIter+1,numStates]);
rhoHRZall      = zeros([unitCnt,1,numIter+1,numStates]);

for t = 1:numel(Trials),

    Trial = Trials{t};    
    unitSubset = units{t};
    fprintf('Processing Trial: %s\nComputational time: ',Trial.filebase)        
    
% LOAD marker position from motion capture data
    xyz = preproc_xyz(Trial,'trb');
    xyz.resample(sampleRate);
    
% LOAD local field potential (lfp)
% RESAMPLE lfp to xyz sample rate
% COMPUTE lfp phase in theta band (6-12 Hz)
    try, lfp = Trial.load('lfp',sessionList(t).thetaRef);
    catch, lfp = Trial.load('lfp',sessionList(t).thetaRef);
    end
    phz = lfp.phase([5,13]);    
    phz.data = unwrap(phz.data);
    phz.resample(xyz);    
    phz.data = mod(phz.data+pi,2*pi)-pi;
    lfp.resample(xyz);

% LOAD theta state placefields
    pft = pfs_2d_theta(Trial);
% COMPUTE placefield centered rate scaled distance metric with approach/depart signature
    hrz = compute_hrz(Trial,unitSubset,pft,'sampleRate',sampleRate);
    ddz = compute_ddz(Trial,unitSubset,pft,'sampleRate',sampleRate);
    drz = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);

% COMPUTE head frame of reference vectors for all time points
    fxyz = filter(copy(xyz),'ButFilter',3,30,'low');
    hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    

% NAN sample points where spikes occured lateral to the head
    pfhr = nan([size(xyz,1),numel(unitSubset),2]);
    for u = 1:numel(unitSubset),%&pft.data.spar>0.15&pft.data.spar<0.3),
        [mxr,mxp] = pft.maxRate(unitSubset(u));
        pfhr(:,u,:) = multiprod(bsxfun(@minus,mxp,sq(fxyz(:,'hcom',[1,2]))),hvec,2,[2,3]);
        hrz(abs(pfhr(:,u,2))>100,u) = nan;
    end
    
% ISOLATE approaching and retreating trajectories
    drzp = drz;    drzp(drzp<0)=nan;
    ddzp = ddz;    ddzp(ddzp<0)=nan;    
    drzn = drz;    drzn(drzn>0)=nan;
    ddzn = ddz;    ddzn(ddzn>0)=nan;    
    hrzp = hrz;    hrzp(hrzp<0)=nan;
    hdzp = ddz;    hdzp(hdzp<0)=nan;    
    hrzn = hrz;    hrzn(hrzn>0)=nan;
    hdzn = ddz;    hdzn(hdzn>0)=nan;        

% LOCATE unit row in group stats matrix
    if t == 1,  ind = 1:sessionUnitCnt(1);
    else,       ind = (sum(sessionUnitCnt(1:(t-1)))+1):sum(sessionUnitCnt(1:t));
    end

% COMPUTE phase precession stats for all states    
    for s = 1:numStates,
        spkpp = Trial.spk.copy();
        spkpp.create(Trial,xyz.sampleRate,states{s},unitSubset,'deburst');

        [parmDRZall(ind,:,:,s), phzStatsDRZall(ind,:,:,s), rDRZall(ind,:,:,s), rhoDRZall(ind,:,:,s)] = ...
            MjgER2016_phasePrecession(Trial,drz,ddz,phz,spkpp,unitSubset,[],[],numIter,...
                                      ['DRZ-',states{s},'-',num2str(numIter)],overwrite);

% $$$         [parmDRZpos(ind,:,s), phzStatsDRZpos(ind,:,s), rDRZpos(ind,:,s), rhoDRZpos(ind,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,drzp,ddzp,phz,spkpp,unitSubset);
% $$$         [parmDRZneg(ind,:,s), phzStatsDRZneg(ind,:,s), rDRZneg(ind,:,s), rhoDRZneg(ind,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,drzn,ddzn,phz,spkpp,unitSubset);
 
        [parmHRZall(ind,:,:,s), phzStatsHRZall(ind,:,:,s), rHRZall(ind,:,:,s), rhoHRZall(ind,:,:,s)] = ...
            MjgER2016_phasePrecession(Trial,hrz,ddz,phz,spkpp,unitSubset,[],[],numIter,...
                                      ['HRZ-',states{s},'-',num2str(numIter)],overwrite);        
% $$$         [parmHRZpos(ind,:,s), phzStatsHRZpos(ind,:,s), rHRZpos(ind,:,s), rhoHRZpos(ind,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,hrzp,hdzp,phz,spkpp,unitSubset);
% $$$         [parmHRZneg(ind,:,s), phzStatsHRZneg(ind,:,s), rHRZneg(ind,:,s), rhoHRZneg(ind,:,s)] = ...
% $$$             MjgER2016_phasePrecession(Trial,hrzn,hdzn,phz,spkpp,unitSubset);
    end
        
    
end



overwrite = false;
numIter = 100;

parmDRZall     = zeros([unitCnt,2,numIter+1,numStates]);
phzStatsDRZall = zeros([unitCnt,2,numIter+1,numStates]);
rDRZall        = zeros([unitCnt,1,numIter+1,numStates]);
rhoDRZall      = zeros([unitCnt,1,numIter+1,numStates]);
parmHRZall     = zeros([unitCnt,2,numIter+1,numStates]);
phzStatsHRZall = zeros([unitCnt,2,numIter+1,numStates]);
rHRZall        = zeros([unitCnt,1,numIter+1,numStates]);
rhoHRZall      = zeros([unitCnt,1,numIter+1,numStates]);

for t = 1:numel(Trials),
    Trial = Trials{t};    
    unitSubset = units{t};
    fprintf('Processing Trial: %s\nComputational time: ',Trial.filebase)        
    if t == 1,
        ind = 1:sessionUnitCnt(1);
    else,
        ind = (sum(sessionUnitCnt(1:(t-1)))+1):sum(sessionUnitCnt(1:t));
    end

    for s = 1:numStates,
        [parmDRZall(ind,:,:,s), phzStatsDRZall(ind,:,:,s), rDRZall(ind,:,:,s), rhoDRZall(ind,:,:,s)] = ...
            MjgER2016_phasePrecession(Trial,[],[],[],[],unitSubset,[],[],numIter,...
                                      ['DRZ-',states{s},'-',num2str(numIter)],overwrite);
        [parmHRZall(ind,:,:,s), phzStatsHRZall(ind,:,:,s), rHRZall(ind,:,:,s), rhoHRZall(ind,:,:,s)] = ...
            MjgER2016_phasePrecession(Trial,[],[],[],[],unitSubset,[],[],numIter,...
                                      ['HRZ-',states{s},'-',num2str(numIter)],overwrite);        
    end
end


% END ACCUMULATE phase precession stats ------------------------------------------------------------


figure,plot(parmHRZall(:,1,1,3),parmHRZall(:,2,1,3),'.');
figure,plot(parmHRZall(:,1,1,5),parmHRZall(:,2,1,5),'.');


figure,hold('on');
plot(parmHRZall(:,1,1,5),parmHRZall(:,2,1,5),'.');
plot(parmHRZall(:,1,1,3),parmHRZall(:,2,1,3),'.r');



figure,
subplot(221),
hist((rhoHRZall(:,1,1,3)-mean(rhoHRZall(:,1,2:end,3),3))./std(rhoHRZall(:,1,2:end,3),[],3),30);
%sum((rhoHRZall(:,1,1,3)-mean(rhoHRZall(:,1,2:end,3),3))./std(rhoHRZall(:,1,2:end,3),[],3)<-2)
subplot(223),
plot((rhoHRZall(:,1,1,3)-mean(rhoHRZall(:,1,2:end,3),3))./std(rhoHRZall(:,1,2:end,3),[],3), ...
rhoHRZall(:,1,1,3),'.')
title('Zscore vs Rho')
subplot(224),plot(mean(rhoHRZall(:,1,2:end,3),3),rhoHRZall(:,1,1,3),'.')
title('E[Rho_{shf}] vs Rho')



% LOAD EXAMPLE DATA --------------------------------------------------------------------------------

pfds = req20180123_ver5(Trials,[],'13');

exampleUnits = [74,83,79,59,103];

t = 20;
Trial = Trials{t};    
unitSubset = units{t};
pft = pfs_2d_theta(Trial);
pfd = req20180123_ver5(Trial,[],'13');



if ~exist('spkhrz','var')|isempty(spkhrz), req20180621(); end

csthresh = 2;
csoffset = 1;
cluSessionMapSubset = cluSessionMap(unitSubsets{1},:);
cluSessionMapSubset_L = cluSessionMapSubset(fsrcz(:,1)>2.1&fsrcz(:,2)<(csthresh-csoffset)&fsrcz(:,3)<(csthresh-csoffset),:);
cluSessionMapSubset_H = cluSessionMapSubset(fsrcz(:,3)>2.1&fsrcz(:,1)<(csthresh-csoffset)&fsrcz(:,2)<(csthresh-csoffset),:);
cluSessionMapSubset_R = cluSessionMapSubset(fsrcz(:,2)>2.1&fsrcz(:,1)<(csthresh-csoffset)&fsrcz(:,3)<(csthresh-csoffset),:);
cluSessionMapSubset_N = cluSessionMapSubset(fsrcz(:,2)<csthresh&fsrcz(:,1)<csthresh&fsrcz(:,3)<csthresh,:);
cluSessionMapSubset_C = cluSessionMapSubset(~ismember(cluSessionMapSubset,...
                                                  [cluSessionMapSubset_L;...
                                                   cluSessionMapSubset_H;...
                                                   cluSessionMapSubset_R;...
                                                   cluSessionMapSubset_N],'rows'),:);
%cluSessionMapSubset = cluSessionMap(~ismember(unitSubsets{1},:));


pfdMaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfds(:,1),units');
pfdMaps = cat(2,pfdMaps{:});

indmatHighLow = zeros(pfds{1}.adata.binSizes');
indmatHighLow(pfds{1}.adata.bins{1}<-0.25,pfds{1}.adata.bins{1}<0.5) = 1;
indmatHighLow(pfds{1}.adata.bins{1}>-0.25,pfds{1}.adata.bins{1}<0.5) = 2;
indmatHighLow(:,pfds{1}.adata.bins{1}>0.5) = 3;

hlScore = (mean(pfdMaps(indmatHighLow(:)==1,:),'omitnan')-mean(pfdMaps(indmatHighLow(:)==2,:),'omitnan'))'./ ...
          (mean(pfdMaps(indmatHighLow(:)==2,:),'omitnan')+mean(pfdMaps(indmatHighLow(:)==1,:),'omitnan'))';
rScore =  mean(pfdMaps(indmatHighLow(:)==3,:),'omitnan');
rScore = (mean(pfdMaps(indmatHighLow(:)==3,:),'omitnan')-mean(mean(pfdMaps(ismember(indmatHighLow(:),[1,2]),:),'omitnan')))'./ ...
         (mean(pfdMaps(indmatHighLow(:)==3,:),'omitnan')+mean(mean(pfdMaps(ismember(indmatHighLow(:),[1,2]),:),'omitnan')))';

% FIGURE START -------------------------------------------------------------------------------------
cond_round = @(rate) max([round(rate,0),round(rate,1)].*[rate>=10,rate<10]);
nanColor = [0.15,0.15,0.15];

pageWidth  = 21.0;
pageHeight = 29.7;

pwidth = 2;
pheight = 2;

xpad = 0.0;
ypad = 0.1;

xpos = 3.5:(pwidth+xpad):pageWidth;
ypos = fliplr(0.5:(pheight+ypad):pageHeight-3.7);

% SET figure opts
hfig = figure(666001);
hfig.Units = 'centimeters';
hfig.Position = [1, 1, pageWidth,pageHeight];
hfig.PaperPositionMode = 'auto';

clf();

sp = gobjects([1,0]);
fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
xlim([0,hfig.Position(3)]);
ylim([0,hfig.Position(4)]);



% PLOT rate maps and phase precession examples
for u = 1:numel(exampleUnits),
    maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                             [{pft},pfd],repmat({exampleUnits(u)},[1,1+numel(pfd)]))));
% PLOT theta example    
    yind = u;
    xind = 1;    
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    plot(pft,exampleUnits(u),'mean',false,[],true,0.5,false,interpParPfs,@jet,[],nanColor);
    text(-490,-380,num2str(cond_round(pft.maxRate(exampleUnits(u)))),'FontSize',10,'Color',[1,1,1]);
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    axes(fax);
    rectangle('Position',sp(end).Position,'LineWidth',1);

% PLOT Behavior field restricted to theta placefield center
    xind = 2;
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    plot(pfd{1},exampleUnits(u),'mean',false,[],false,0.5,false,interpParDfs,@jet,[],nanColor);
    text(-1.7,-0.45,num2str(cond_round(maxPfsRate)),'FontSize',10,'Color',[1,1,1]);
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    xlim([-1.8,0.5]);
    ylim([-0.75,1.8]);
    axes(fax);
    rectangle('Position',sp(end).Position,'LineWidth',1);
    
    
% PLOT phase precession for each state
    for s = 1:numel(states)-1
        xind = 2+s;
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(xind),ypos(yind),pwidth,pheight],...
                         'FontSize', 8,...
                         'LineWidth',1);
        
        uind =  ismember(spkmap,[t,exampleUnits(u)],'rows') ...
                & logical(spkstc(:,1)) ...
                & logical(spkstc(:,statesInds(s))) ...
                & abs(spkego(:,2))<100;
        
        hold('on');
        plot([spkhrz(uind);spkhrz(uind)],...
             [circ_rad2ang(spkphz(uind));circ_rad2ang(spkphz(uind))+360],...
             '.','MarkerSize',5);
        xlim([-1,1]);
        ylim([-180,540]);
        
        if u == 1
            title(statesLabels{s+1});
        end

        if s ~= numel(states)-1,
            sp(end).XTickLabels = {};
            sp(end).YTickLabels = {};
        else,
            sp(end).YAxisLocation = 'right';
            sp(end).YTick = [0,180,360,540];
        end
        
        uExInd = ismember(cluSessionMap,[t,exampleUnits(u)],'rows');
        
        plot([-1,1],circ_rad2ang(P_H(uExInd,s+1,1)*[-1,1]+P_H(uExInd,s+1,2)),'-r','LineWidth',1)
        plot([-1,1],circ_rad2ang(P_H(uExInd,s+1,1)*[-1,1]+P_H(uExInd,s+1,2))+360,'-r','LineWidth',1)

        axes(fax);
        rectangle('Position',sp(end).Position,'LineWidth',1);
                              
        
% $$$         plot([0,1],circ_rad2ang(PP_H(uExInd,s+1,1)*[0,1]+PP_H(uExInd,s+1,2)),'-r','LineWidth',1)
% $$$         plot([0,1],circ_rad2ang(PP_H(uExInd,s+1,1)*[0,1]+PP_H(uExInd,s+1,2))+360,'-r','LineWidth',1)
% $$$ 
% $$$         plot([-1,0],circ_rad2ang(PN_H(uExInd,s+1,1)*[-1,0]+PN_H(uExInd,s+1,2)),'-g','LineWidth',1)
% $$$         plot([-1,0],circ_rad2ang(PN_H(uExInd,s+1,1)*[-1,0]+PN_H(uExInd,s+1,2))+360,'-g','LineWidth',1)

        %text(-(['HRZ rho: ',num2str(round(rho_H(uExInd,s+1),2))]);
        
    end        
end%for u
% END plot phase precession examples


% $$$ ind = all(spkstc(:,[1,7]),2)&ismember(spkmap(:,1),17:23)&ismember(spkmap,cluSessionMapSubset_L,'rows');
% $$$ ind = all(spkstc(:,[1,5]),2)&ismember(spkmap(:,1),17:23)&ismember(spkmap,cluSessionMapSubset_HLoc,'rows');
% $$$ ind = all(spkstc(:,[1,4]),2)&ismember(spkmap(:,1),17:23)&ismember(spkmap,cluSessionMapSubset_Rear,'rows');
% $$$ %ind = all(spkstc(:,[1,2]),2)&ismember(spkmap(:,1),17:23)&ismember(spkmap,[18,11],'rows');
% $$$ 
% $$$ figure,plot(spkdrz(ind),spkphz(ind),'.')
% $$$ 
% $$$ figure,plot(spkhrz(ind),spkphz(ind),'.')
% $$$ figure,
% $$$ statesInds = [4,5,6,7,8];
% $$$ statesIndComp = {[5,6,7,8],7,8,5,6};
% $$$ statesUSubs  = {cluSessionMapSubset_R, ...
% $$$                 cluSessionMapSubset_H, ...
% $$$                 cluSessionMapSubset_H, ...
% $$$                 cluSessionMapSubset_L, ...
% $$$                 cluSessionMapSubset_L};
% $$$ for s = 1:numel(statesInds);
% $$$     subplot2(2,5,1,s);
% $$$     ind = all(spkstc(:,[1,statesInds(s)]),2)&ismember(spkmap(:,1),17:23)& ...
% $$$           ismember(spkmap,cluSessionMapSubset_LLoc,'rows');
% $$$     plot(spkhrz(ind),spkphz(ind),'.');
% $$$     subplot2(2,5,2,s);
% $$$     ind = spkstc(:,1)==1&any(spkstc(:,statesIndComp{s}),2)...
% $$$           &ismember(spkmap(:,1),17:23)& ...
% $$$           ismember(spkmap,statesUSubs{s},'rows');
% $$$     plot(spkhrz(ind),spkphz(ind),'.');
% $$$ end

statesInds = [4,5,6,7,8];
statesIndComp = {[5,6,7,8],7,8,5,6};

statesUSubs  = {cluSessionMapSubset_R, ...
                cluSessionMapSubset_H, ...
                cluSessionMapSubset_H, ...
                cluSessionMapSubset_L, ...
                cluSessionMapSubset_L};

% $$$ statesUSubs  = {cluSessionMap(rScore>-0.5  & RmaxN>0.25,:), ...
% $$$                 cluSessionMap(hlScore<-0.1 & rScore<-0.2 & RmaxN>0.25,:), ...
% $$$                 cluSessionMap(hlScore<-0.1 & rScore<-0.2 & RmaxN>0.25,:), ...
% $$$                 cluSessionMap(hlScore>0.1  & rScore<-0.2 & RmaxN>0.25,:), ...
% $$$                 cluSessionMap(hlScore>0.1  & rScore<-0.2 & RmaxN>0.25,:)};


% CA1 
% JPDF HRZ vs PHZ : all units 
sesIds = [3:5,8:12,17:23];
tind = spkstc(:,1)==1 & ismember(spkmap(:,1),sesIds);
% JPDF of plot VS spatial position for all units
for s = 1:numel(statesInds);
    xind = 2+s;
% CREATE subplot axes
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+1)-pheight/2,pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
% INDEX data
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows') & abs(spkego(:,2)<100);
%    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMap(genericScores{s},:),'rows');
%    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMap(RmaxN(:,s+1)>0.25,:),'rows');
% PLOT jpdf of data
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    %title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
    title(statesLabels{s+1});    
    %caxis([0,150]);
    if s == 1, ylabel({'CA1'}); end    
    %if s == 1, ylabel({'Behavior','Specific Units'}); end
    if s == numel(statesInds),
        sp(end).YAxisLocation = 'right';
        sp(end).YTick = [0,pi,2*pi,3*pi];
        sp(end).YTickLabels = [0,180,360,540];
        sp(end).XTickLabels = {};        
    else        
        sp(end).XTickLabels = {};
        sp(end).YTickLabels = {};
    end    
    Lines([],pi,'k');    
end


% CA3 
% JPDF HRZ vs PHZ : all units 
sesIds = [1,2,6,7,13:16];
tind = spkstc(:,1)==1 & ismember(spkmap(:,1),sesIds);
% JPDF of plot VS spatial position for all units
for s = 1:numel(statesInds);
    xind = 2+s;
% CREATE subplot axes
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+2)-pheight/2,pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
% INDEX data
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows') & abs(spkego(:,2)<100);
%    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMap(genericScores{s},:),'rows');
%    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMap(RmaxN(:,s+1)>0.25,:),'rows');
% PLOT jpdf of data
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    %title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
    %title(statesLabels{s});
    %caxis([0,150]);
    if s == 1, ylabel({'CA2&3'}); end
    if s == numel(statesInds),
        sp(end).YAxisLocation = 'right';
        sp(end).YTick = [0,pi,2*pi,3*pi];
        sp(end).YTickLabels = [0,180,360,540];
        sp(end).XTickLabels = {};        
    else        
        sp(end).XTickLabels = {};
        sp(end).YTickLabels = {};
    end    
    Lines([],pi,'k');    
end


% PLOT phase precession group stats
% 1. slopes of phase precession in units preferred behavior
% 2. slopes of phase precession in units non-preferred behavior
% COMPUTE vars for indexing

% CA1
ppZscores = sq(bsxfun(@minus,rhoHRZall(:,1,1,:),mean(rhoHRZall(:,1,2:end,:),3))./std(rhoHRZall(:,1,2:end,:),[],3));
sesIds = [3:5,8:12,17:23];
zthresh = -2;
rthresh = 0.2;
statesIndCompPhz = [2,5,6,3,4];
% ACCUMULATE slope within units' preferred behavior
prefBhvId = [];
prefBhvSlopes = [];
prefBhvPhases = [];
nprefBhvId = [];
nprefBhvSlopes = [];
nprefBhvPhases = [];
for s = 1:5;    
    ind =   ismember(cluSessionMap,statesUSubs{s},'rows') ...
          & ismember(cluSessionMap(:,1),sesIds)           ...
          & rHRZall(:,1,1,[s+1])   > rthresh;
    %          & ppZscores(:,[s+1]) < zthresh                  ...    
    prefBhvId = cat(1,prefBhvId,s.*ones([sum(ind),1]));    
    prefBhvSlopes = cat(1,prefBhvSlopes,parmHRZall(ind,1,1,s+1));        
    prefBhvPhases = cat(1,prefBhvSlopes,parmHRZall(ind,2,1,s+1));
    nprefBhvId = cat(1,nprefBhvId,(statesIndCompPhz(s)-1).*ones([sum(ind),1]));
    nprefBhvSlopes = cat(1,nprefBhvSlopes,parmHRZall(ind,1,1,statesIndCompPhz(s)));
    nprefBhvPhases = cat(1,nprefBhvSlopes,parmHRZall(ind,2,1,statesIndCompPhz(s)));    
end            
% CREATE subplot axes
xind = 1;    
yind = 6;
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind)-pheight/2,pwidth*1.9,pheight],...
                 'FontSize', 8,...
                 'LineWidth',1);
boxplot([prefBhvSlopes;nprefBhvSlopes],      ...
        [prefBhvId*2-1;nprefBhvId*2],              ...
        'plotstyle',    'traditional',       ...
        'boxstyle',     'filled',            ...
        'colors',       'rrbbccggmm',             ...
        'symbol',       '.',                 ...
        'datalim',      [-10,7.5],           ...
        'labels',       {});
grid('on');


% CA3
xind = 1;    
yind = 7;
sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind)-pheight/2,pwidth*1.9,pheight],...
                 'FontSize', 8,...
                 'LineWidth',1);
sesIds = [1,2,6,7,13:16];
% ACCUMULATE slope within units' preferred behavior
prefBhvId = [];
prefBhvSlopes = [];
prefBhvPhases = [];
nprefBhvId = [];
nprefBhvSlopes = [];
nprefBhvPhases = [];
for s = 1:5;    
    ind =   ismember(cluSessionMap,statesUSubs{s},'rows') ...
          & ismember(cluSessionMap(:,1),sesIds)           ...
          & rHRZall(:,1,1,[s+1])   > rthresh;
    %          & ppZscores(:,[s+1]) < zthresh                  ...    
    prefBhvId = cat(1,prefBhvId,s.*ones([sum(ind),1]));    
    prefBhvSlopes = cat(1,prefBhvSlopes,parmHRZall(ind,1,1,s+1));        
    prefBhvPhases = cat(1,prefBhvSlopes,parmHRZall(ind,2,1,s+1));
    nprefBhvId = cat(1,nprefBhvId,(statesIndCompPhz(s)-1).*ones([sum(ind),1]));
    nprefBhvSlopes = cat(1,nprefBhvSlopes,parmHRZall(ind,1,1,statesIndCompPhz(s)));
    nprefBhvPhases = cat(1,nprefBhvSlopes,parmHRZall(ind,2,1,statesIndCompPhz(s)));    
end            
boxplot([prefBhvSlopes;nprefBhvSlopes],      ...
        [prefBhvId*2-1;nprefBhvId*2],        ...
        'plotstyle',    'traditional',       ...
        'boxstyle',     'filled',            ...
        'colors',       'rrbbccggmm',        ...
        'symbol',       '.',                 ...
        'datalim',      [-10,7.5],           ...
        'labelorientation','inline',       ...
        'labels',       reshape([statesLabels(2:end);repmat({''},[1,numel(statesLabels)-1])],[],1)');        
%        'labels',       reshape(repmat(statesLabels(2:end),[2,1]),[],1)');
grid('on');
sp(end).Position =[xpos(xind),ypos(yind)-pheight/2,pwidth*1.9,pheight];


% END FIGURE 4 -------------------------------------------------------------------------------------






% SUPPLEMENTARY figures of figure 4
% 4.1 examples of phase precession a comparison of DRZ vs HRZ
%     4.1.1 [pfs,bfs,[drz for multiple behaviors],[hrz for multiple behaviors]]
%     4.1.2 R value of HRZ circular-linear correlation vs R value of DRZ circular-linear correlation 
%     4.1.2 Population of Units' Rho of HRZxTHP vs Rho of DRZxTHP
%
% DRZ: Trajectory directed rate zone (Huxter 2008)
% HRZ: Head directed rate zone
% THP: Theta phase
%




% SUP FIG 4.1 --------------------------------------------------------------------------------------




s = 2;
    ind =   ismember(cluSessionMap,statesUSubs{s},'rows') ...
          & ismember(cluSessionMap(:,1),sesIds)           ...
          & rHRZall(:,1,1,[s+1])   > rthresh;

figure();
plot(parmHRZall(ind,1,1,3),parmHRZall(ind,1,1,5),'.');
daspect([1,1,1]);
grid('on');
hold('on');
line([-20,20],[-20,20]);



ppZscoresHRZ = sq(bsxfun(@minus,rhoHRZall(:,1,1,:),mean(rhoHRZall(:,1,2:end,:),3))./std(rhoHRZall(:,1,2:end,:),[],3));
ppZscoresDRZ = sq(bsxfun(@minus,rhoDRZall(:,1,1,:),mean(rhoDRZall(:,1,2:end,:),3))./std(rhoDRZall(:,1,2:end,:),[],3));

figure,
zthresh = -2;
sesIds = [3:5,8:12,17:23];
%sesIds = [1,2,6,7,13:16];
for s = 1:6
    ind = ppZscoresDRZ(:,s) < zthresh & ppZscoresHRZ(:,s) < zthresh & ismember(cluSessionMap(:,1),sesIds);    
    subplot(2,6,s);
    hold('on');
    scatter(rHRZall(ind,1,1,s),rDRZall(ind,1,1,s),15,parmHRZall(ind,1,1,s),'filled');
    caxis([-5,5]);
    line([0,1],[0,1]);
    res = multiprod([rHRZall(ind,1,1,s),rDRZall(ind,1,1,s)],[cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4)],[2],[1,2]);
    resbins = 0.0:0.2:1;
    resbinsInds = discretize(res(:,1),rbins);
    resSTest = accumarray(resbinsInds(nniz(resbinsInds)),res(nniz(resbinsInds),2),[numel(resbins-1),1],@signtest,nan);
    title({statesLabels{s},['Sign Test p=',num2str(signtest(res(:,2)))]});
    ylabel('R_{DRZ}');
    xlabel('R_{HRZ}');
    xlim([0,0.5]);
    ylim([0,0.5]);
    colormap('jet');

    subplot(2,6,s+6);
    hold('on');
    plot(rhoHRZall(ind,1,1,s),rhoDRZall(ind,1,1,s),'.');
    line([-1,1],[-1,1]);    
    res = multiprod([rhoHRZall(ind,1,1,s),rhoDRZall(ind,1,1,s)],[cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4)],[2],[1,2]);
    resBins = -1:0.2:0;
    resBinsInds = discretize(res(:,1),resBins);
    resSTest = accumarray(resBinsInds(nniz(resBinsInds)),res(nniz(resBinsInds),2),[numel(resBins-1),1],@signtest,nan);
    resMeans = multiprod([accumarray(resBinsInds(nniz(resBinsInds)),res(nniz(resBinsInds),1),[numel(resBins-1),1],@mean,nan),...
                          accumarray(resBinsInds(nniz(resBinsInds)),res(nniz(resBinsInds),2),[numel(resBins-1),1],@mean,nan)],...
                         [cos(-pi/4),-sin(-pi/4);sin(-pi/4),cos(-pi/4)],[2],[1,2]);
    plot(resMeans(:,1),resMeans(:,2),'-+');
    title({statesLabels{s},['Sign Test p=',num2str(signtest(res(:,2)))]});
    ylabel('rho_{DRZ}');
    xlabel('rho_{HRZ}');
    xlim([-0.7,0]);
    ylim([-0.7,0]);
end






% TEST the significance between the r values  THPxDRZ vs THPxHRZ
% xy THP x DRZ
% xz THP x HRZ
% yz HRZ x DRZ
% $$$ diff <- xy-xz
% $$$ determin = 1-xy*xy - xz*xz - yz*yz + 2*xy*xz*yz;
% $$$ av=(xy+xz)/2
% $$$ cube= (1-yz)*(1-yz)*(1-yz)
% $$$ t2 = diff * sqrt((n-1)*(1+yz)/(((2*(n-1)/(n-3))*determin+av*av*cube)))
% $$$ p <- pt(abs(t2),n-3,lower.tail=FALSE)    #changed to n-3 12/15/18
% $$$ if(twotailed) p <- 2*p
% $$$     value <- list(test="test of difference between two correlated  correlations",t=t2,p=p,Call=cl)



% NON PREFERRED bhv
% $$$ 
% $$$ % PP slopes vs R
% $$$ uind = ismember(cluSessionMap,cluSessionMap(unitSubsets{1},:),'rows');
% $$$ for s = 1:numel(statesInds);
% $$$     xind = 2+s;
% $$$ % CREATE subplot axes
% $$$     sp(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[xpos(xind),ypos(yind+3)-pheight/2,pwidth,pheight],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$     plot(RmaxN_H(uind,s+1),PN_H(uind,s+1,1),'.')
% $$$     %scatter(RmaxN_H(uind,s),PN_H(uind,s,1),4,bsxfun(@minus,FSrC(uind(unitSubsets{1}),[2,3,1]),[0,0.2,0.2]),'Filled');
% $$$     xlim([0,1]);
% $$$     ylim([-15,15]);
% $$$     Lines([],0,'k');
% $$$     if s==1, ylabel({'entering','field'}); end    
% $$$ end
% $$$ 
% $$$ for s = 1:numel(statesInds);
% $$$     xind = 2+s;
% $$$ % CREATE subplot axes
% $$$     sp(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[xpos(xind),ypos(yind+4)-pheight/2,pwidth,pheight],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$     %scatter(RmaxP_H(uind,s),PP_H(uind,s,1),4,bsxfun(@minus,FSrC(uind(unitSubsets{1}),[2,3,1]),[0,0.2,0.2]),'Filled');    
% $$$     plot(RmaxP_H(uind,s+1),PP_H(uind,s+1,1),'.')
% $$$     xlim([0,1]);
% $$$     ylim([-15,15]);
% $$$     Lines([],0,'k');
% $$$     if s==1, ylabel({'exiting','field'}); end
% $$$ end



% $$$ sesIds = [8,9,10,11,12,17:23];
% $$$ tind = spkstc(:,1)==1 & ismember(spkmap(:,1),sesIds);
% $$$ hrzBins = linspace(-1,1,21);
% $$$ % SUBDIVIDED by behavioral selectivity
% $$$ for s = 1:numel(statesInds);
% $$$     xind = 2+s;
% $$$     sp(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[xpos(xind),ypos(yind+5)-pheight/2,pwidth,pheight],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$     ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows');          
% $$$     hold('on');
% $$$ % $$$     plot(linspace(hrzBins(1),hrzBins(end),numel(hrzBins)-1),...
% $$$ % $$$          accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_mean)','+');
% $$$     errorbar(linspace(hrzBins(1),hrzBins(end),numel(hrzBins)-1)',...
% $$$              accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_mean)',...
% $$$              accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_var)'./2);    
% $$$     sp(end).XTickLabels = {};
% $$$     sp(end).YTickLabels = {};
% $$$     title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
% $$$     Lines([],pi,'k');    
% $$$ end
% $$$ 
% $$$ sesIds = [8,9,10,11,12,17:23];
% $$$ tind = spkstc(:,1)==1 & ismember(spkmap(:,1),sesIds);
% $$$ hrzBins = linspace(-1,1,21);
% $$$ % SUBDIVIDED by behavioral selectivity
% $$$ for s = 1:numel(statesInds);
% $$$     xind = 2+s;
% $$$     sp(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[xpos(xind),ypos(yind+5)-pheight/2,pwidth,pheight],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$     ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows');          
% $$$     hold('on');
% $$$ % $$$     plot(linspace(hrzBins(1),hrzBins(end),numel(hrzBins)-1),...
% $$$ % $$$          accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_mean)','+');
% $$$     errorbar(linspace(hrzBins(1),hrzBins(end),numel(hrzBins)-1)',...
% $$$              accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_mean)',...
% $$$              accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_var)'./2);    
% $$$     sp(end).XTickLabels = {};
% $$$     sp(end).YTickLabels = {};
% $$$     title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
% $$$     Lines([],pi,'k');    
% $$$ end





    
% $$$     sp(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[xpos(xind),ypos(yind+2)-pheight/2,pwidth,pheight],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$     ind = tind & any(spkstc(:,statesIndComp{s}),2) & ismember(spkmap,statesUSubs{s},'rows');
% $$$     hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
% $$$           linspace(-1,1,21),...
% $$$           linspace(-pi,3*pi,21));    
% $$$     sp(end).XTickLabels = {};
% $$$     sp(end).YTickLabels = {};
% $$$     %caxis([0,150]); 
% $$$     Lines([],pi,'k');    
% $$$     if s == 1, ylabel({'Complementary','behavior'}); end    
% $$$     if s == numel(statesInds),
% $$$         sp(end).YAxisLocation = 'right';
% $$$         sp(end).YTick = [0,pi,2*pi,3*pi];
% $$$         sp(end).YTickLabels = [0,180,360,540];
% $$$     else        
% $$$         sp(end).XTickLabels = {};
% $$$         sp(end).YTickLabels = {};
% $$$     end
% $$$     
% $$$     
% $$$     sp(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[xpos(xind),ypos(yind+4),pwidth,pheight],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$     ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMapSubset_N,'rows');
% $$$     hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
% $$$           linspace(-1,1,21),...
% $$$           linspace(-pi,3*pi,21));    
% $$$     if s == numel(statesInds),
% $$$         sp(end).YAxisLocation = 'right';
% $$$         sp(end).YTick = [0,pi,2*pi,3*pi];
% $$$         sp(end).YTickLabels = [0,180,360,540];
% $$$     else        
% $$$         sp(end).XTickLabels = {};
% $$$         sp(end).YTickLabels = {};
% $$$     end
% $$$     Lines([],pi,'k');
% $$$     title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);    
% $$$     if s == 1, ylabel({'Behavior','non-selective'}); end        
% $$$     %caxis([0,150]);
% $$$ end




% CA1 
% JPDF HRZ vs PHZ : selective behaviors
sesIds = [8,9,10,11,12,17:23];
tind = spkstc(:,1)==1 & ismember(spkmap(:,1),sesIds);
figure
for s = 1:numel(statesInds);
% PREFFERED state
    xind = 2+s;
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+1)-pheight/2,pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows') & abs(spkego(:,2)<100);
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
    if s == 1, ylabel({'Behavior','Specific Units'}); end
    if s == numel(statesInds),
        sp(end).YAxisLocation = 'right';
        sp(end).YTick = [0,pi,2*pi,3*pi];
        sp(end).YTickLabels = [0,180,360,540];
        sp(end).XTickLabels = {};        
    else        
        sp(end).XTickLabels = {};
        sp(end).YTickLabels = {};
    end    
    Lines([],pi,'k');    

% NONPREFFERED state    
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+2)-pheight/2,pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesIndComp{s}),2) & ismember(spkmap,statesUSubs{s},'rows') & abs(spkego(:,2)<100);
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));    
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    Lines([],pi,'k');    
    if s == 1, ylabel({'Complementary','behavior'}); end    
    if s == numel(statesInds),
        sp(end).YAxisLocation = 'right';
        sp(end).YTick = [0,pi,2*pi,3*pi];
        sp(end).YTickLabels = [0,180,360,540];
    else        
        sp(end).XTickLabels = {};
        sp(end).YTickLabels = {};
    end
    
% THEREST state        
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+4),pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMapSubset_N,'rows') & abs(spkego(:,2)<100);
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));    
    if s == numel(statesInds),
        sp(end).YAxisLocation = 'right';
        sp(end).YTick = [0,pi,2*pi,3*pi];
        sp(end).YTickLabels = [0,180,360,540];
    else        
        sp(end).XTickLabels = {};
        sp(end).YTickLabels = {};
    end
    Lines([],pi,'k');
    title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);    
    if s == 1, ylabel({'Behavior','non-selective'}); end        
end








statesUSubs  = {cluSessionMapSubset_R, ...
                cluSessionMapSubset_H, ...
                cluSessionMapSubset_H, ...
                cluSessionMapSubset_L, ...
                cluSessionMapSubset_L};

statesUSubs  = {cluSessionMap(rScore>-0.5 ,:), ...
                cluSessionMap(hlScore<-0.1 & rScore<-0.2,:), ...
                cluSessionMap(hlScore<-0.1 & rScore<-0.2,:), ...
                cluSessionMap(hlScore>0.1  & rScore<-0.2,:), ...
                cluSessionMap(hlScore>0.1  & rScore<-0.2,:)};

ppZscores = sq(bsxfun(@minus,rhoHRZall(:,1,1,:),mean(rhoHRZall(:,1,2:end,:),3))./std(rhoHRZall(:,1,2:end,:),[],3));
sesIds = [8,9,10,11,12,17:23];
zthresh = -2;
rthresh = 0.2;
statesIndCompPhz = [4,5,6,3,4];

ppSlopeBins = linspace(-10,10,100);



s = 1;
figure,
subplot(131);plot(rHRZall(:,1,1,s),rhoHRZall(:,1,1,s),'.');
subplot(132);plot(rHRZall(:,1,1,s),ppZscores(:,s),'.');
subplot(133);scatter(rHRZall(:,1,1,s),parmHRZall(:,1,1,s),10,ppZscores(:,s),'Filled');
caxis([-5,5])
colormap(gca,'jet');

figure();plot3(rHRZall(:,1,1,1),ppZscores(:,1),rhoHRZall(:,1,1,1),'.');

% $$$ uind = ismember(cluSessionMap,cluSessionMap(unitSubsets{1},:),'rows') & ismember(cluSessionMap(:,1),[17:23]);
figure,
sp = gobjects([0,1]);
sb = gobjects([0,1]);
for s = 1:5;
    sp(end+1) = subplot(5,5,s);
    ind = ismember(cluSessionMap,statesUSubs{s},'rows') ...
          & ismember(cluSessionMap(:,1),sesIds) ...
          & (any(ppZscores(:,[s+1])<zthresh,2))...          
          & any(rHRZall(:,1,1,[s+1])>rthresh,4);% ...          
    %& (any(ppZscores(:,[s+1])<zthresh,2))...

%          & (any(ppZscores(:,[s+1,statesIndCompPhz(s)])<zthresh,2));% ...    
    %| any(rHRZall(:,1,1,[s+1,statesIndCompPhz(s)])<rthresh,4));
    indR = ind & all(reshape(WithinRanges(parmHRZall(:,1,1,[s+1,statesIndCompPhz(s)]),[-10,10]),[],2),2);
                                                                                                 %indR = ismember(cluSessionMap,statesUSubs{s},'rows') & any(abs(ppZscores(:,[s+1]))>zthresh,2) & any(sq(rHRZall(:,1,1,[s+1]))>0.2,2);
% $$$     ind = ismember(cluSessionMap,statesUSubs{s},'rows') & any(abs(ppZscores(:,[3,5]))>2,2);%|sq(rHRZall(:,1,1,[3,5]))>0.2,2);
% $$$     indR = ismember(cluSessionMap,statesUSubs{s},'rows') & any(abs(ppZscores(:,[3,5]))>2,2) & any(sq(rHRZall(:,1,1,[3,5]))>0.2,2);
% PREFERRED bhv
    plot(parmHRZall(ind,1,1,s+1),...
         parmHRZall(ind,2,1,s+1)+double(parmHRZall(ind,2,1,s+1)<0).*2.*pi,'.');
    title(statesLabels{s+1});
% NON PREFERRED bhv
    sp(end+1) = subplot(5,5,s+5);
    plot(parmHRZall(ind,1,1,statesIndCompPhz(s)),...
         parmHRZall(ind,2,1,statesIndCompPhz(s))+double(parmHRZall(ind,2,1,statesIndCompPhz(s))<0).*2.*pi,'.');    

% MARGINAL pp slope
    sb(end+1) = subplot(5,5,s+10); 
    hold('on');
    hax = bar(ppSlopeBins,histc(parmHRZall(ind,1,1,s+1),ppSlopeBins),'histc');
    hax.FaceColor = 'b';
    hax.EdgeColor = 'b';
    hax.FaceAlpha = 0.4;
    hax.EdgeAlpha = 0.4;
    ylim([0,20]);
    %Lines(mean(nonzeros(parmHRZall(ind,1,1,s+1).*double(WithinRanges(parmHRZall(ind,1,1,s+1),[-7.5,5]))),'omitnan'),[],'r');

    hax = bar(ppSlopeBins,histc(parmHRZall(ind,1,1,statesIndCompPhz(s)),ppSlopeBins),'histc');
    hax.FaceColor = 'r';
    hax.EdgeColor = 'r';
    hax.FaceAlpha = 0.4;
    hax.EdgeAlpha = 0.4;
    ylim([0,20]);
    %Lines(mean(nonzeros(parmHRZall(ind,1,1,statesIndCompPhz(s)).*double(WithinRanges(parmHRZall(ind,1,1,statesIndCompPhz(s)),[-7.5,5]))),'omitnan'),[],'b');    


% TRANSITION 
    sp(end+1) = subplot(5,5,s+15);
    hold('on');
    plot([parmHRZall(indR,1,1,statesIndCompPhz(s)),parmHRZall(indR,1,1,s+1)]',...
         [parmHRZall(indR,2,1,statesIndCompPhz(s))+double(parmHRZall(indR,2,1,statesIndCompPhz(s))<0).*2.*pi,...
          parmHRZall(indR,2,1,s+1)+double(parmHRZall(indR,2,1,s+1)<0).*2.*pi]','b');    
    plot(parmHRZall(ind,1,1,s+1),...
         parmHRZall(ind,2,1,s+1)+double(parmHRZall(ind,2,1,s+1)<0).*2.*pi,'.r')
% $$$     quiver(parmHRZall(ind,1,1,s+1),...
% $$$            parmHRZall(ind,2,1,s+1),...
% $$$            parmHRZall(ind,1,1,statesIndCompPhz(s))-parmHRZall(ind,1,1,s+1),...
% $$$            parmHRZall(ind,2,1,statesIndCompPhz(s))-parmHRZall(ind,2,1,s+1),...
% $$$            0);
    subplot(5,5,s+20);    
    plot([parmHRZall(indR,1,1,statesIndCompPhz(s))-parmHRZall(indR,1,1,s+1)],...
         [parmHRZall(indR,2,1,statesIndCompPhz(s))+double(parmHRZall(indR,2,1,statesIndCompPhz(s))<0).*2.*pi]-...
          [parmHRZall(indR,2,1,s+1)+double(parmHRZall(indR,2,1,s+1)<0).*2.*pi],'.b');
    xlim([-10,10]);
    ylim([-5,5]);

end
ForAllSubplots('grid(''on'');');
af(@(s) ylim(s,[0,2*pi]), sp);
af(@(s) set(s,'YTick',[0:pi/4:2*pi]),sp);
af(@(s) set(s,'XTick',[-10:2.5:10]),[sp,sb]);
af(@(s) xlim(s,[-7.5,5]), [sp,sb]);
linkaxes(sp,'xy');


figure,
plot(parmHRZall(:,1,1,6),parmHRZall(:,2,1,3),'.');







% $$$ ind = spkstc(:,1)==1&spkstc(:,1)~=4&any(spkstc(:,[5,6,7,8]),2)...
% $$$       &ismember(spkmap(:,1),17:23)& ...
% $$$       ismember(spkmap,statesUSubs{s},'rows');
% $$$ %vbins = discretize(spkvxy(ind,3))
% $$$ phzBins = linspace(-pi,3*pi,41);
% $$$ phzBinInds = discretize(reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1),phzBins);
% $$$ hrzBins = linspace(-1,1,21);
% $$$ hrzBinInds = discretize(repmat(spkhrz(ind),[2,1]),hrzBins);
% $$$ figure,imagesc(accumarray([hrzBinInds,phzBinInds],repmat(spkvxy(ind,3),[2,1]),[numel(hrzBins),numel(phzBins)],@mean)');







% OTHER FIG ---------------------------------


figure();
binDims = [30,20];

binsTrans = linspace(-pi/2,pi/2,binDims(2));
binsTransInd = discretize(spkpch(:,2),binsTrans);
labelTrans = 'Body Pitch (rad)';
binsDrz = linspace(-1,1,binDims(1));
binsDrzInd  =  discretize(spkdrz,binsDrz);
labelDRZ = 'DRZ';
saveLabel = 'drz_x_bpitch';
 

binsTrans = linspace(-10,10,binDims(2));
binsTransInd = discretize(spksvd(:,2),binsTrans);
labelTrans = 'svdpc2';
binsDrz = linspace(-1,1,binDims(1));
binsDrzInd  =  discretize(spkdrz,binsDrz);
labelDRZ = 'DRZ';
saveLabel = 'drz_x_fsvd';


binsTrans = linspace(-1,0.6,binDims(2));
binsTransInd = discretize(spkpch(:,1),binsTrans);
labelTrans = 'Head Pitch (rad)';
binsDrz = linspace(-1,1,binDims(1));
binsDrzInd  =  discretize(spkdrz,binsDrz);
labelDRZ = 'DRZ';
saveLabel = 'drz_x_hpitch';

binsTrans = linspace(20,300,binDims(2));
binsTransInd = discretize(spkhgt,binsTrans);
labelTrans = 'Height (mm)';
binsDrz = linspace(-1,1,binDims(1));
binsDrzInd  =  discretize(spkdrz,binsDrz);
labelDRZ = 'DRZ';
saveLabel = 'drz_x_height';

clf();
ny = 4;
for y = 1:ny,
    switch y
      case 1
      % DRZ X Head Speed
        binsTrans = linspace(-0.5,1.9,binDims(2));
        binsTransInd = discretize(log10(spkvxy(:,2)),binsTrans);
        labelTrans = 'Head Speed log10(cm/s)';
        binsDrz = linspace(-1,1,binDims(1));
        binsDrzInd  =  discretize(spkdrz,binsDrz);
        %binsDrzInd  =  discretize(spkdrz,binsDrz);
        labelDRZ = 'DRZ';
        aveLabel = 'drz_x_hspeed';
      case 2
        % HRZ X Head Speed
        binsTrans = linspace(-0.5,1.9,binDims(2));
        binsTransInd = discretize(log10(spkvxy(:,2)),binsTrans);
        labelTrans = 'Head Speed log10(cm/s)';
        binsDrz = linspace(-1,1,binDims(1));
        binsDrzInd  =  discretize(spkhrz,binsDrz);
        labelDRZ = 'HRZ';
        saveLabel = 'hrz_x_hspeed';
      case 3
        % DRZ X Body Speed
        binsTrans = linspace(-0.5,1.9,binDims(2));
        binsTransInd = discretize(log10(spkvxy(:,1)),binsTrans);
        labelTrans = 'Body Speed log10(cm/s)';
        binsDrz = linspace(-1,1,binDims(1));
        binsDrzInd  =  discretize(spkdrz,binsDrz);
        labelDRZ = 'DRZ';
        saveLabel = 'drz_x_bspeed';
      case 4
        % HRZ X Body Speed
        binsTrans = linspace(-0.5,1.9,binDims(2));
        binsTransInd = discretize(log10(spkvxy(:,1)),binsTrans);
        labelTrans = 'Body Speed log10(cm/s)';
        binsDrz = linspace(-1,1,binDims(1));
        binsDrzInd  =  discretize(spkhrz,binsDrz);
        labelDRZ = 'HRZ';
        saveLabel = 'hrz_x_bspeed';
    end


% $$$ sind = 9;
% $$$ binsTrans = linspace(-1,1,binDims(2));
% $$$ binsTransInd =  discretize(spktrans{1,sind},binsTrans);
% $$$ labelTrans = [states{sind},' state onset'];
% $$$ binsDrz = linspace(-1,1,binDims(1));
% $$$ binsDrzInd  =  discretize(spkdrz,binsDrz);
% $$$ labelDRZ = 'DRZ';
% $$$ saveLabel = ['drz_x_',states{sind},'ON'];



%sigUnitsBhv = true([size(FSrC,1),1]);
%sigUnitsBhv = any(FSrC(:,[1,3])>=-0.5,2);
%sigUnitsBhv = any(FSrC(:,[1,3])>=-0,2);
%sigUnitsBhv = any(fsrcz(:,[1,3])>=2,2);
%sigUnitsBhv = any(FSrC(:,[2])<0,2);
%sigUnitsBhv = any(FSrC(:,[2])<=-0,2);
%sigUnitsBhv = true([size(fsrcz,1),1]);
%sigUnitsBhv = all(fsrcz(:,[1,3])>-2,2)&all(fsrcz(:,[1,3])<2,2)&any(fsrcz(:,[2])<2,2);
%sigUnitsBhv = any(fsrcz(:,[1])>1,2);
%sigUnitsBhv = any(fsrcz(:,[2])<0,2);

ind = nniz(binsTransInd)&nniz(binsDrzInd)           ...
& spkstc(:,1)                                     ... theta
&~spkstc(:,9)                                     ... not groom
&~spkstc(:,10)                                    ... not sit
&~spkstc(:,4);                                     ... not rear
%& ismember(spkmap,cluSessionSubset(sigUnitsBhv,:),'rows');
A = accumarray([binsDrzInd(ind),binsTransInd(ind)],spkphz(ind),[numel(binsDrz)-1,numel(binsTrans)-1],@circ_mean);
S = accumarray([binsDrzInd(ind),binsTransInd(ind)],spkphz(ind),[numel(binsDrz)-1,numel(binsTrans)-1],@circ_std);
C = accumarray([binsDrzInd(ind),binsTransInd(ind)],ones([sum(ind),1]),[numel(binsDrz)-1,numel(binsTrans)-1],@sum);


%A(A<0) = A(A<0)+2*pi;

subplot2(ny,3,y,1);imagesc(binsDrz,binsTrans,A');axis('xy');colormap(gca,'hsv');
xlabel(labelDRZ);ylabel(labelTrans);
cax = colorbar();ylabel(cax,'Mean Theta Phase');
subplot2(ny,3,y,2);imagesc(binsDrz,binsTrans,S');colorbar();axis('xy');colormap(gca,'default');
xlabel(labelDRZ);ylabel(labelTrans);
cax = colorbar();ylabel(cax,'STD Theta Phase');
subplot2(ny,3,y,3);imagesc(binsDrz,binsTrans,C');colorbar();axis('xy');colormap(gca,'default');
xlabel(labelDRZ);ylabel(labelTrans);
cax = colorbar();ylabel(cax,'Count');

end


hax = findobj(gcf,'Type','Axes');
af(@(h) set(h,'Units','centimeters'),  hax);
af(@(h) set(h,'Position',[h.Position(1:2),2.5,2.5]),  hax);

print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
       ['pp_pop_drz_hrz.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
       ['pp_pop_drz_hrz.png']]);




















% $$$ 
% $$$ % PLOT Slope histogram of units
% $$$ figure,hist(PN(uind,1,1),100)
% $$$ 
% $$$ s = 3;
% $$$ figure,plot(PN(uind,s,1),exp(-phzStatsN(uind,s,2).^2./2),'.')
% $$$ 
% $$$ 
% $$$ s = 3
% $$$ figure,
% $$$ subplot(121);
% $$$ plot(P(uind,s,1),Rmax(uind,s),'.')
% $$$ ylim([0,1]);
% $$$ xlim([-10,10]);
% $$$ subplot(122);
% $$$ plot(P_H(uind,s,1),Rmax_H(uind,s),'.')
% $$$ ylim([0,1]);
% $$$ xlim([-10,10]);
% $$$ 
% $$$ 
% $$$ 
% $$$ s = 4
% $$$ figure,
% $$$ subplot(121);
% $$$ plot(PN(uind,s,1),RmaxN(uind,s),'.')
% $$$ ylim([0,1]);
% $$$ xlim([-10,10]);
% $$$ subplot(122);
% $$$ plot(PN_H(uind,s,1),RmaxN_H(uind,s),'.')
% $$$ ylim([0,1]);
% $$$ xlim([-10,10]);
% $$$ 
% $$$ 
% $$$ s = 3;
% $$$ figure()
% $$$ plot(Rmax_H(uind,s),Rmax_H(uind,s+2),'.');
% $$$ xlim([0,1]);
% $$$ ylim([0,1]);
% $$$ line([0,1],[0,1]);
% $$$ 
% $$$ 
% $$$ s = 3;
% $$$ figure,
% $$$ subplot(121);
% $$$ plot(PN(uind,s,1),rhoN(uind,s),'.'),ylim([0,1]);
% $$$ ylim([-1,1]);
% $$$ xlim([-10,10]);
% $$$ subplot(122);
% $$$ plot(PN_H(uind,s,1),rhoN_H(uind,s),'.')
% $$$ ylim([-1,1]);
% $$$ xlim([-10,10]);
% $$$ 
% $$$ 
% $$$ 
% $$$ uind = ismember(cluSessionMap,cluSessionMap(unitSubsets{1},:),'rows');
% $$$ uind = ismember(cluSessionMap,cluSessionMap(unitSubsets{1},:),'rows') & ismember(cluSessionMap(:,1),[17:23]);
% $$$ 
% $$$ Rmax_H
% $$$ 
% $$$ s = 3;
% $$$ figure,
% $$$ subplot(121);
% $$$ plot(P_H(uind,s,1),P_H(uind,s+2,1),'.');
% $$$ xlabel(['slope, ' statesLabels{s}]);
% $$$ ylabel(['slope, ' statesLabels{s+2}]);
% $$$ ylim([-7.5,5]);
% $$$ xlim([-7.5,5]);
% $$$ line([-7.5,5],[-7.5,5]);
% $$$ subplot(122);
% $$$ plot(P_H(uind,s,2),P_H(uind,s+2,2),'.');
% $$$ ylim([-pi,pi]);
% $$$ xlim([-pi,pi]);
% $$$ line([-pi,pi],[-pi,pi]);
% $$$ 
% $$$ s = 3;
% $$$ figure,
% $$$ subplot(121);
% $$$ plot(PN_H(uind,s,1),PN_H(uind,s+2,1),'.');
% $$$ ylim([-7.5,5]);
% $$$ xlim([-7.5,5]);
% $$$ line([-7.5,5],[-7.5,5]);
% $$$ subplot(122);
% $$$ plot(PN_H(uind,s,2),PN_H(uind,s+2,2),'.');
% $$$ ylim([-pi,pi]);
% $$$ xlim([-pi,pi]);
% $$$ line([-pi,pi],[-pi,pi]);
