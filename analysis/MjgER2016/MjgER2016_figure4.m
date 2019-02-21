% MjgER2016 Figure4
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

statesLabels = {'theta','rear',{'H Loc'},{'H Pause'},{'L Loc'},{'L Pause'}};

uind = ismember(cluSessionMap,cluSessionMap(unitSubsets{1},:),'rows');


sampleRate = 250;

P  = {};   phzStats  = {};   Rmax  = {}; rho  = {};
PP = {};   phzStatsP = {};   RmaxP = {}; rhoP = {};
PN = {};   phzStatsN = {};   RmaxN = {}; rhoN = {};

P_H  = {};   phzStats_H  = {};   Rmax_H  = {}; rho_H= {};
PP_H = {};   phzStatsP_H = {};   RmaxP_H = {}; rhoP_H= {};
PN_H = {};   phzStatsN_H = {};   RmaxN_H = {}; rhoN_H= {};

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

    pft = pfs_2d_theta(Trial);
    hrz = compute_hrz(Trial,unitSubset,pft,'sampleRate',sampleRate);
    ddz = compute_ddz(Trial,unitSubset,pft,'sampleRate',sampleRate);
    drz = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);

    fxyz = filter(copy(xyz),'ButFilter',3,30,'low');
    hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    

    pfhr = nan([size(xyz,1),numel(unitSubset),2]);
    for u = 1:numel(unitSubset),%&pft.data.spar>0.15&pft.data.spar<0.3),
        [mxr,mxp] = pft.maxRate(unitSubset(u));
        pfhr(:,u,:) = multiprod(bsxfun(@minus,mxp,sq(fxyz(:,'hcom',[1,2]))),hvec,2,[2,3]);
        hrz(abs(pfhr(:,u,2))>100,u) = nan;
    end
    
    

    drzp = drz;    drzp(drzp<0)=nan;
    ddzp = ddz;    ddzp(ddzp<0)=nan;    
    drzn = drz;    drzn(drzn>0)=nan;
    ddzn = ddz;    ddzn(ddzn>0)=nan;    
    hrzp = hrz;    hrzp(hrzp<0)=nan;
    hdzp = ddz;    hdzp(hdzp<0)=nan;    
    hrzn = hrz;    hrzn(hrzn>0)=nan;
    hdzn = ddz;    hdzn(hdzn>0)=nan;        

    if t == 1,
        ind = 1:sessionUnitCnt(1);
    else,
        ind = (sum(sessionUnitCnt(1:(t-1)))+1):sum(sessionUnitCnt(1:t));
    end

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
        
% $$$         % DRZ phase precession
% $$$         [P{t,s},  phzStats{t,s},  Rmax{t,s},  rho{t,s}] = ...
% $$$             MjgER2016_phasePrecession(Trial,drz,ddz,phz,spkpp,unitSubset);
% $$$         [PP{t,s}, phzStatsP{t,s},RmaxP{t,s}, rhoP{t,s}] = ...
% $$$             MjgER2016_phasePrecession(Trial,drzp,ddzp,phz,spkpp,unitSubset);
% $$$         [PN{t,s}, phzStatsN{t,s},RmaxN{t,s}, rhoN{t,s}] = ...
% $$$             MjgER2016_phasePrecession(Trial,drzn,ddzn,phz,spkpp,unitSubset);
% $$$ 
% $$$         % HRZ phase precession
% $$$         [P_H{t,s},   phzStats_H{t,s}, Rmax_H{t,s}, rho_H{t,s}] = ...
% $$$             MjgER2016_phasePrecession(Trial,hrz,ddz,phz,spkpp,unitSubset);
% $$$         [PP_H{t,s}, phzStatsP_H{t,s},RmaxP_H{t,s}, rhoP_H{t,s}] = ...
% $$$             MjgER2016_phasePrecession(Trial,hrzp,hdzp,phz,spkpp,unitSubset);
% $$$         [PN_H{t,s}, phzStatsN_H{t,s},RmaxN_H{t,s}, rhoN_H{t,s}] = ...
% $$$             MjgER2016_phasePrecession(Trial,hrzn,hdzn,phz,spkpp,unitSubset);
% $$$     end
    
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

% DRZ 
% $$$ P = cf(@(p) permute(p,[1,3,2]), P);   P = cell2mat(P);
% $$$ PP = cf(@(p) permute(p,[1,3,2]), PP); PP = cell2mat(PP);
% $$$ PN = cf(@(p) permute(p,[1,3,2]), PN); PN = cell2mat(PN);
% $$$ phzStats = cf(@(p) permute(p,[1,3,2]), phzStats);    phzStats = cell2mat(phzStats);
% $$$ phzStatsP = cf(@(p) permute(p,[1,3,2]), phzStatsP);  phzStatsP = cell2mat(phzStatsP);
% $$$ phzStatsN = cf(@(p) permute(p,[1,3,2]), phzStatsN);  phzStatsN = cell2mat(phzStatsN);
% $$$ Rmax = cf(@(p) permute(p,[1,3,2]), Rmax);   Rmax = cell2mat(Rmax);
% $$$ RmaxP = cf(@(p) permute(p,[1,3,2]), RmaxP); RmaxP = cell2mat(RmaxP);
% $$$ RmaxN = cf(@(p) permute(p,[1,3,2]), RmaxN); RmaxN = cell2mat(RmaxN);
% $$$ rho = cf(@(p) permute(p,[1,3,2]), rho);     rho = cell2mat(rho);
% $$$ rhoP = cf(@(p) permute(p,[1,3,2]), rhoP);   rhoP = cell2mat(rhoP);
% $$$ rhoN = cf(@(p) permute(p,[1,3,2]), rhoN);   rhoN = cell2mat(rhoN);
% $$$ 
% $$$ % HRZ 
% $$$ P_H = cf(@(p) permute(p,[1,3,2]), P_H);   P_H = cell2mat(P_H);
% $$$ PP_H = cf(@(p) permute(p,[1,3,2]), PP_H); PP_H = cell2mat(PP_H);
% $$$ PN_H = cf(@(p) permute(p,[1,3,2]), PN_H); PN_H = cell2mat(PN_H);
% $$$ phzStats_H = cf(@(p) permute(p,[1,3,2]), phzStats_H);   phzStats = cell2mat(phzStats_H);
% $$$ phzStatsP_H = cf(@(p) permute(p,[1,3,2]), phzStatsP_H); phzStatsP = cell2mat(phzStatsP_H);
% $$$ phzStatsN_H = cf(@(p) permute(p,[1,3,2]), phzStatsN_H); phzStatsN = cell2mat(phzStatsN_H);
% $$$ Rmax_H = cf(@(p) permute(p,[1,3,2]), Rmax_H);   Rmax_H = cell2mat(Rmax_H);
% $$$ RmaxP_H = cf(@(p) permute(p,[1,3,2]), RmaxP_H); RmaxP_H = cell2mat(RmaxP_H);
% $$$ RmaxN_H = cf(@(p) permute(p,[1,3,2]), RmaxN_H); RmaxN_H = cell2mat(RmaxN_H);
% $$$ rho_H = cf(@(p) permute(p,[1,3,2]), rho_H);     rho_H = cell2mat(rho_H);
% $$$ rhoP_H = cf(@(p) permute(p,[1,3,2]), rhoP_H);   rhoP_H = cell2mat(rhoP_H);
% $$$ rhoN_H = cf(@(p) permute(p,[1,3,2]), rhoN_H);   rhoN_H = cell2mat(rhoN_H);
% $$$ 
% exp(-s^2./2)


%figure,plot(parmHRZall(:,1,1,6),parmHRZall(:,2,1,6),'.');


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

exampleUnits = [74,83,79,59,103];

t = 20;
Trial = Trials{t};    
unitSubset = units{t};

pft = pfs_2d_theta(Trial);
pfd = req20180123_ver5(Trial,[],'13');
pfds = req20180123_ver5(Trials,[],'13');


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

genericScores = {rScore>-0.5  & RmaxN>0.25, ...
                 hlScore<-0.1 & rScore<-0.2 & RmaxN>0.25, ...
                 hlScore<-0.1 & rScore<-0.2 & RmaxN>0.25, ...
                 hlScore>0.1  & rScore<-0.2 & RmaxN>0.25, ...
                 hlScore>0.1  & rScore<-0.2 & RmaxN>0.25};

sesIds = [17:23];

%sesIds = [8,9,10,11,12];

% JPDF HRZ vs PHZ : all units
sesIds = [8,9,10,11,12,17:23];
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
%    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows');          
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMap(genericScores{s},:),'rows');
%    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMap(RmaxN(:,s+1)>0.25,:),'rows');
% PLOT jpdf of data
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
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
%    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows');          
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMap(genericScores{s},:),'rows');
%    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMap(RmaxN(:,s+1)>0.25,:),'rows');
% PLOT jpdf of data
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
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


% PP slopes vs R
uind = ismember(cluSessionMap,cluSessionMap(unitSubsets{1},:),'rows');
for s = 1:numel(statesInds);
    xind = 2+s;
% CREATE subplot axes
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+3)-pheight/2,pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    plot(RmaxN_H(uind,s+1),PN_H(uind,s+1,1),'.')
    %scatter(RmaxN_H(uind,s),PN_H(uind,s,1),4,bsxfun(@minus,FSrC(uind(unitSubsets{1}),[2,3,1]),[0,0.2,0.2]),'Filled');
    xlim([0,1]);
    ylim([-15,15]);
    Lines([],0,'k');
    if s==1, ylabel({'entering','field'}); end    
end

for s = 1:numel(statesInds);
    xind = 2+s;
% CREATE subplot axes
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+4)-pheight/2,pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    %scatter(RmaxP_H(uind,s),PP_H(uind,s,1),4,bsxfun(@minus,FSrC(uind(unitSubsets{1}),[2,3,1]),[0,0.2,0.2]),'Filled');    
    plot(RmaxP_H(uind,s+1),PP_H(uind,s+1,1),'.')
    xlim([0,1]);
    ylim([-15,15]);
    Lines([],0,'k');
    if s==1, ylabel({'exiting','field'}); end
end







sesIds = [8,9,10,11,12,17:23];
tind = spkstc(:,1)==1 & ismember(spkmap(:,1),sesIds);
hrzBins = linspace(-1,1,21);
% SUBDIVIDED by behavioral selectivity
for s = 1:numel(statesInds);
    xind = 2+s;
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+5)-pheight/2,pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows');          
    hold('on');
% $$$     plot(linspace(hrzBins(1),hrzBins(end),numel(hrzBins)-1),...
% $$$          accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_mean)','+');
    errorbar(linspace(hrzBins(1),hrzBins(end),numel(hrzBins)-1)',...
             accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_mean)',...
             accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_var)'./2);    
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
    Lines([],pi,'k');    
end

sesIds = [8,9,10,11,12,17:23];
tind = spkstc(:,1)==1 & ismember(spkmap(:,1),sesIds);
hrzBins = linspace(-1,1,21);
% SUBDIVIDED by behavioral selectivity
for s = 1:numel(statesInds);
    xind = 2+s;
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+5)-pheight/2,pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows');          
    hold('on');
% $$$     plot(linspace(hrzBins(1),hrzBins(end),numel(hrzBins)-1),...
% $$$          accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_mean)','+');
    errorbar(linspace(hrzBins(1),hrzBins(end),numel(hrzBins)-1)',...
             accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_mean)',...
             accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_var)'./2);    
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
    Lines([],pi,'k');    
end


    
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+2)-pheight/2,pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesIndComp{s}),2) & ismember(spkmap,statesUSubs{s},'rows');
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));    
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    %caxis([0,150]); 
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
    
    
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+4),pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMapSubset_N,'rows');
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
    %caxis([0,150]);
end




% selective behaviors
sesIds = [8,9,10,11,12,17:23];
tind = spkstc(:,1)==1 & ismember(spkmap(:,1),sesIds);
figure
for s = 1:numel(statesInds);
    xind = 2+s;
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+1)-pheight/2,pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows');          
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
    %caxis([0,150]);
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
    
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+2)-pheight/2,pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesIndComp{s}),2) & ismember(spkmap,statesUSubs{s},'rows');
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));    
    sp(end).XTickLabels = {};
    sp(end).YTickLabels = {};
    %caxis([0,150]); 
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
    
    
    sp(end+1) = axes('Units','centimeters',...
                     'Position',[xpos(xind),ypos(yind+4),pwidth,pheight],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMapSubset_N,'rows');
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
    %caxis([0,150]);
end










s = 2;
P(ismember(spkmap,statesUSubs{s},'rows');



ind = spkstc(:,1)==1&spkstc(:,1)~=4&any(spkstc(:,[5,6,7,8]),2)...
      &ismember(spkmap(:,1),17:23)& ...
      ismember(spkmap,statesUSubs{s},'rows');

%vbins = discretize(spkvxy(ind,3))
phzBins = linspace(-pi,3*pi,41);
phzBinInds = discretize(reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1),phzBins);

hrzBins = linspace(-1,1,21);
hrzBinInds = discretize(repmat(spkhrz(ind),[2,1]),hrzBins);


figure,imagesc(accumarray([hrzBinInds,phzBinInds],repmat(spkvxy(ind,3),[2,1]),[numel(hrzBins),numel(phzBins)],@mean)');







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
