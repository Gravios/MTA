% MjgER2016 Figure4
%
% HP:= Head Pitch
% BP:= Body Pitch
% PFD:= behavior field restricted to DRZ[-0.5,0.5]
%
%% ACCUMULATE phase precession stats ----------------------------------------------------------------
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

global MTA_PROJECT_PATH

MjgER2016_load_data();
MjgER2016_general_args('section 1');

%%%<<< SET args
overwrite = false;
stateLabels = {'theta','rear','high','low','hloc','hpause','lloc','lpause'};
statesLabels = {'theta','rear','H Loc','H Pause','L Loc','L Pause'};
sampleRate = 250;
sigma = 150;
sessionUnitCnt = cellfun(@numel,units);
unitCnt = sum(sessionUnitCnt);
phzOrder = [9:16,1:8];
%CA1 units
sesIds = [3:5,8:12,17:23];
%CA3 units
%sesIds = [1:2,6:7,13:16];
sesInds = ismember(cluSessionMap(:,1),sesIds);
statesPfss = {'loc&theta','lloc&theta','hloc&theta','rear&theta',         ...
              'pause&theta','lpause&theta','hpause&theta'};
statesPfssInd = [4,3,2];
stsColor = 'rgb';
%%%>>>

%%%<<< LOAD data

% DEF place field rate maps
pfts        = cf(@(t,u)  pfs_2d_theta(t,u),                                       Trials, units);
pfss        = cf(@(t,u)  pfs_2d_states(t,u),                                      Trials, units);
% DEF behavior field rate maps
pfbs        = cf(@(t,u)  compute_bhv_ratemaps(t,u),                               Trials, units);
pfbsShuff   = cf(@(t,u)  compute_bhv_ratemaps_shuffled(t,u),                      Trials, units);
% DEF Conditional Expectation: req20190527(): phase VS gdz (gaussian distance zone)
pftHZTPD    = cf(@(s) ...
                 cf(@(T,u) MTAApfs(T,u,'tag',['ddtp-','s',num2str(sigma),'-',s]), Trials, units),...
                 stateLabels);
% LOAD bhv ratemap erpPCA
[eigVecs,eigScrs,eigVars,unitSubset,validDims] = compute_bhv_ratemaps_erpPCA(pfbs,units,'overwrite',overwrite);
cluSessionMapSubset = cluSessionMap(unitSubset,:);
% LOAD bhv ratemap erpPCA Scores
[fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS,fsrsMean,fsrsStd] = ...
    compute_bhv_ratemaps_erpPCA_scores(Trials,units,pfbs,pfbsShuff,eigVecs,validDims,unitSubset,overwrite);


% GENERAL variables
phaseBinCenters = pftHZTPD{1}{1}.adata.bins{1}(phzOrder)+2*pi ...
                  .*double(pftHZTPD{1}{1}.adata.bins{1}(phzOrder)<0);

uind = ismember(cluSessionMap,cluSessionMap(unitSubset,:),'rows');

trjCntFilePath = fullfile(MTA_PROJECT_PATH,'analysis',...
                          ['unit_uniqueTrajectoryCount-',DataHash({cf(@(t)t.filebase,Trials),units}),'.mat'])
if ~exist(trjCntFilePath,'file'),
    [tcount,scount,tper,states] = cf(@(t,u) compute_unit_uniqueTrajectoryCount(t,u), Trials,units);
    [tcount,scount,tper,states] = deal(cat(1,tcount{:}),cat(1,scount{:}),cat(1,tper{:}),cat(1,states{:}));
    save(trjCntFilePath,'tcount','scount','tper','states');
else
    load(trjCntFilePath);
end

%[cluSessionMap(cluSessionMap(:,1)==19,2),diff(sq(tper(cluSessionMap(:,1)==19,1,:)),1,2)./60]
dbins = pftHZTPD{1}{1}.adata.bins{1};
pbins = circ_rad2ang([pftHZTPD{1}{1}.adata.bins{2};pftHZTPD{1}{1}.adata.bins{2}+2*pi]);
sclr = 'rgb';
sclrm = eye(3);
colorMap = @cool;
%%%>>>

%req20181220

% LOAD phase precession statistics
MjgER2016_figure_thetaPhasePrecession_stats();
MjgER2016_figure_thetaPhasePrecession_drzCorrectedPhasePreference();
MjgER2016_figure_thetaPhasePrecession_correctedPhasePreferenceDecomposition_statePermutations();

%%%<<< DIAGNOSTIC figures
figure,plot(parmHRZall(:,1,1,3),parmHRZall(:,2,1,3),'.');
figure,plot(parmHRZall(:,1,1,5),parmHRZall(:,2,1,5),'.');
figure,
    hold('on');
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
%%%>>>


% LOAD EXAMPLE DATA --------------------------------------------------------------------------------

exampleUnits = [20,74;...
                20,83;...
                20,79;...
                20,59;...
                20,103];

expUnitsPP = {[20],[79,25]};

if ~exist('spkhrz','var')|isempty(spkhrz), req20180621(); end

% INDEX cluSessionMap for each spike
[~,spkCluSessionMapInd] = max(all(bsxfun(@eq,spkmap,permute(cluSessionMap,[3,2,1])),2)==true,[],3);
spkppr = sq(rHRZall(spkCluSessionMapInd,1,1,:));
spkppr(isnan(spkppr)) = 0;

csthresh = 2.1;
csoffset = 1;
cluSessionMapSubset_L = cluSessionMapSubset(  fsrcz(:,1)>csthresh                                ...
                                            & fsrcz(:,2)<(csthresh-csoffset)                     ...
                                            & fsrcz(:,3)<(csthresh-csoffset),:);
cluSessionMapSubset_H = cluSessionMapSubset(  fsrcz(:,3)>csthresh                                ...
                                            & fsrcz(:,1)<(csthresh-csoffset)                     ...
                                            & fsrcz(:,2)<(csthresh-csoffset),:);
cluSessionMapSubset_R = cluSessionMapSubset(  fsrcz(:,2)>csthresh                                ...
                                            & fsrcz(:,1)<(csthresh-csoffset)                     ...
                                            & fsrcz(:,3)<(csthresh-csoffset),:);
cluSessionMapSubset_N = cluSessionMapSubset(  fsrcz(:,2)<csthresh                                ...
                                            & fsrcz(:,1)<csthresh                                ...
                                            & fsrcz(:,3)<csthresh,:);
cluSessionMapSubset_C = cluSessionMapSubset(~ismember(cluSessionMapSubset,...
                                                  [cluSessionMapSubset_L;...
                                                   cluSessionMapSubset_H;...
                                                   cluSessionMapSubset_R;...
                                                   cluSessionMapSubset_N],'rows'),:);
%cluSessionMapSubset = cluSessionMap(~ismember(unitSubset,:));


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


pfdMaps = cf(@(p,u) mean(p{1}.data.rateMap(:,ismember(p{1}.data.clu,u),:),3,'omitnan'), cf(@(p) pfbs(p), expUnitsPP(:,1)),expUnitsPP(:,2));
pfdMaps = cat(2,pfdMaps{:});

% $$$ indmatHighLow = zeros(pfbs{1}.adata.binSizes');
% $$$ indmatHighLow(pfbs{1}.adata.bins{1}<-0.25,pfbs{1}.adata.bins{1}<0.5) = 1;
% $$$ indmatHighLow(pfbs{1}.adata.bins{1}>-0.25,pfbs{1}.adata.bins{1}<0.5) = 2;
% $$$ indmatHighLow(:,pfbs{1}.adata.bins{1}>0.5) = 3;
% $$$ 
% $$$ hlScore = (mean(pfdMaps(indmatHighLow(:)==1,:),'omitnan')-mean(pfdMaps(indmatHighLow(:)==2,:),'omitnan'))'./ ...
% $$$           (mean(pfdMaps(indmatHighLow(:)==2,:),'omitnan')+mean(pfdMaps(indmatHighLow(:)==1,:),'omitnan'))';
% $$$ rScore =  mean(pfdMaps(indmatHighLow(:)==3,:),'omitnan');
% $$$ rScore = (mean(pfdMaps(indmatHighLow(:)==3,:),'omitnan')-mean(mean(pfdMaps(ismember(indmatHighLow(:),[1,2]),:),'omitnan')))'./ ...
% $$$          (mean(pfdMaps(indmatHighLow(:)==3,:),'omitnan')+mean(mean(pfdMaps(ismember(indmatHighLow(:),[1,2]),:),'omitnan')))';

%% FIGURE START ------------------------------------------------------------------------------------
%% STARTFIG ----------------------------------------------------------------------------------------
cond_round = @(rate) max([round(rate,0),round(rate,1)].*[rate>=10,rate<10]);
nanColor = [0.15,0.15,0.15];



% SET figure opts
[hfig,fig,fax,sax] = set_figure_layout(figure(666001),'A4','portrait',[],1.5,1.5);




%%%<<< PLOT rate maps and phase precession examples
for tind = 1:size(expUnitsPP,1),
    t = expUnitsPP{tind,1}(1);
    for uind = 1:numel(expUnitsPP{tind,2}),
        unit = expUnitsPP{tind,2}(uind);
        maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                                     [pfts(t),pfbs(t)],repmat({unit},[1,1+numel(pfbs{t})]))));
        %%%<<< PLOT theta example    
        yind = uind;
        xind = 1;    
        sax(end+1) = axes('Units','centimeters',...
                          'Position',[fig.page.xpos(xind),...
                            fig.page.ypos(yind),fig.subplot.width,fig.subplot.height],...
                          'FontSize', 8,...
                          'LineWidth',1);
        plot(pfts{t},unit,'mean',false,[],true,0.5,false,interpParPfs,colorMap,[],nanColor);
        text(-480,-380,num2str(cond_round(pfts{t}.maxRate(unit,true,'interpPar',interpParPfs))),'FontSize',10,'Color',[1,1,1]);
        sax(end).XTickLabels = {};
        sax(end).YTickLabels = {};
        axes(fax);
        rectangle('Position',sax(end).Position,'LineWidth',1);
        %%%>>>

        %%%<<< PLOT Behavior field restricted to theta placefield center
        xind = 2;
        sax(end+1) = axes('Units','centimeters',                               ...
                          'Position',[fig.page.xpos(xind),                      ...
                            fig.page.ypos(yind),                      ...
                            fig.subplot.width,                        ...
                            fig.subplot.height],                      ...
                          'FontSize', 8,                                        ...
                          'LineWidth',1);
        %plot(pfbs{t},unit,'mean',false,[],false,0.5,false,interpParDfs,@jet,reshape(validDims,pfbs{t}.adata.binSizes'),nanColor);
        plot(pfbs{t},unit,'mean',[],[],false,0.5,false,...
             interpParDfs,colorMap,reshape(validDims,pfbs{t}.adata.binSizes'),nanColor);    
        text(-1.7,-0.45,num2str(cond_round(pfbs{t}.maxRate(unit,false,'mask',validDims))),'FontSize',10,'Color',[1,1,1]);
        %text(-1.7,-0.45,num2str(cond_round(maxPfsRate)),'FontSize',10,'Color',[1,1,1]);
        sax(end).XTickLabels = {};
        sax(end).YTickLabels = {};
        xlim([-1.8,0.5]);
        ylim([-0.75,1.8]);
        axes(fax);
        rectangle('Position',sax(end).Position,'LineWidth',1);
        %%%>>>
        
        %%%<<< PLOT phase precession for each state
        for s = 1:numel(states)-1
            xind = 2+s;
            sax(end+1) = axes('Units','centimeters',...
                              'Position',[fig.page.xpos(xind),...
                                fig.page.ypos(yind),...
                                fig.subplot.width,...
                                fig.subplot.height],...
                              'FontSize', 8,...
                              'LineWidth',1);
            
            uind =  ismember(spkmap,[t,unit],'rows') ...
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
                sax(end).XTickLabels = {};
                sax(end).YTickLabels = {};
            else,
                sax(end).YAxisLocation = 'right';
                sax(end).YTick = [0,180,360,540];
            end
            
            uExInd = ismember(cluSessionMap,[t,unit],'rows');
            
            plot([-1,1],circ_rad2ang(parmHRZall(uExInd,1,1,s+1)*[-1,1]+parmHRZall(uExInd,2,1,s+1)),'-r','LineWidth',1)
            plot([-1,1],circ_rad2ang(parmHRZall(uExInd,1,1,s+1)*[-1,1]+parmHRZall(uExInd,2,1,s+1))+360,'-r','LineWidth',1)

            axes(fax);
            rectangle('Position',sax(end).Position,'LineWidth',1);
        end
        %%%>>>
    end%for u
end%for t
% END plot phase precession examples ---------------------------------------------------------------
%%%>>>

%%%<<< JPDF HRZ vs PHZ : all spikes
% CA1 ----------------------------------------------------------------------------------------------
% JPDF HRZ vs PHZ : all spikes
tind = spkstc(:,1)==1 & ismember(spkmap(:,1),sesIds);
% JPDF of plot VS spatial position for all units
for s = 1:numel(statesInds);
    xind = 2+s;
% CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',...
                     'Position',[fig.page.xpos(xind),...
                                 fig.page.ypos(yind+1)-fig.subplot.height/2,...
                                 fig.subplot.width,...
                                 fig.subplot.height],...
                     'FontSize', 8,...
                     'LineWidth',1);
% INDEX data
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows') ...
          & abs(spkego(:,2)<100) & spkppr(:,s+1) > 0.2;
% PLOT jpdf of data
    hist2([repmat(spkghz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));
    sax(end).XTickLabels = {};
    sax(end).YTickLabels = {};
    title(statesLabels{s+1});    
    if s == 1, ylabel({'Selective'}); end    
    %if s == 1, ylabel({'Behavior','Specific Units'}); end
    if s == numel(statesInds),
        sax(end).YAxisLocation = 'right';
        sax(end).YTick = [0,pi,2*pi,3*pi];
        sax(end).YTickLabels = [0,180,360,540];
        sax(end).XTickLabels = {};        
    else        
        sax(end).XTickLabels = {};
        sax(end).YTickLabels = {};
    end    
    Lines([],pi,'k');    
end
%%%>>>

%%%<<< JPDF HRZ vs PHZ : Behaviorally non selective neurons
% CA1 ----------------------------------------------------------------------------------------------
% RESTRICT anaysis to Sessions with CA1 units and theta activity
tind = spkstc(:,1)==1 & ismember(spkmap(:,1),sesIds);
% JPDF of plot VS spatial position for all units
for s = 1:numel(statesInds);
    xind = 2+s;
% CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',...
                     'Position',[fig.page.xpos(xind),fig.page.ypos(yind+2)-fig.subplot.height/2,fig.subplot.width,fig.subplot.height],...
                     'FontSize', 8,...
                     'LineWidth',1);
% INDEX data
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMapSubset_N,'rows') & ...
          abs(spkego(:,2)<100) & spkppr(:,s+1) > 0.2;
% PLOT jpdf of data
    hist2([repmat(spkghz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));
    sax(end).XTickLabels = {};
    sax(end).YTickLabels = {};
    if s == 1, ylabel({'Conical'}); end    
    %if s == 1, ylabel({'Behavior','Specific Units'}); end
    if s == numel(statesInds),
        sax(end).YAxisLocation = 'right';
        sax(end).YTick = [0,pi,2*pi,3*pi];
        sax(end).YTickLabels = [0,180,360,540];
        sax(end).XTickLabels = {};        
    else        
        sax(end).XTickLabels = {};
        sax(end).YTickLabels = {};
    end    
    Lines([],pi,'k');    
end
%%%>>>

%%%<<< REMOVED : too few samples for phase precession analysis in CA2/3
% $$$ % CA3 
% $$$ % JPDF HRZ vs PHZ : all units 
% $$$ sesIds = [1,2,6,7,13:16];
% $$$ tind = spkstc(:,1)==1 & ismember(spkmap(:,1),sesIds);
% $$$ % JPDF of plot VS spatial position for all units
% $$$ for s = 1:numel(statesInds);
% $$$     xind = 2+s;
% $$$ % CREATE subplot axes
% $$$     sp(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[fig.page.xpos(xind),fig.page.ypos(yind+2)-fig.subplot.height/2,fig.subplot.width,fig.subplot.height],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$ % INDEX data
% $$$     ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows') & abs(spkego(:,2)<100);
% $$$ %    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMap(genericScores{s},:),'rows');
% $$$ %    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMap(RmaxN(:,s+1)>0.25,:),'rows');
% $$$ % PLOT jpdf of data
% $$$     hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
% $$$           linspace(-1,1,21),...
% $$$           linspace(-pi,3*pi,21));
% $$$     sax(end).XTickLabels = {};
% $$$     sax(end).YTickLabels = {};
% $$$     %title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
% $$$     %title(statesLabels{s});
% $$$     %caxis([0,150]);
% $$$     if s == 1, ylabel({'CA2&3'}); end
% $$$     if s == numel(statesInds),
% $$$         sax(end).YAxisLocation = 'right';
% $$$         sax(end).YTick = [0,pi,2*pi,3*pi];
% $$$         sax(end).YTickLabels = [0,180,360,540];
% $$$         sax(end).XTickLabels = {};        
% $$$     else        
% $$$         sax(end).XTickLabels = {};
% $$$         sax(end).YTickLabels = {};
% $$$     end    
% $$$     Lines([],pi,'k');    
% $$$ end
%%%>>>

%%%<<< PLOT phase precession group stats
% PLOT phase precession group stats ----------------------------------------------------------------
% 1. slopes of phase precession in units preferred behavior
% 2. slopes of phase precession in units non-preferred behavior
% COMPUTE vars for indexing

%%%<<< CA1 Behaviorally selective population phase precession statistics
% CA1
%ppZscores = sq(bsxfun(@minus,...
%                      rhoHRZall(:,1,1,:),...
%                      mean(rhoHRZall(:,1,2:end,:),3))./std(rhoHRZall(:,1,2:end,:),[],3));
%sesIds = [3:5,8:12,17:23];
sesIds = [8:12,17:23];
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
    prefBhvId     = cat(1,prefBhvId,s.*ones([sum(ind),1]));    
    prefBhvSlopes = cat(1,prefBhvSlopes,parmHRZall(ind,1,1,s+1));        
    prefBhvPhases = cat(1,prefBhvSlopes,parmHRZall(ind,2,1,s+1));
    nprefBhvId     = cat(1,nprefBhvId,(statesIndCompPhz(s)-1).*ones([sum(ind),1]));
    nprefBhvSlopes = cat(1,nprefBhvSlopes,parmHRZall(ind,1,1,statesIndCompPhz(s)));
    nprefBhvPhases = cat(1,nprefBhvSlopes,parmHRZall(ind,2,1,statesIndCompPhz(s)));    
end
% CREATE subplot axes
xind = 1;    
yind = 3;
sax(end+1) = axes('Units','centimeters',...
                 'Position',[fig.page.xpos(xind),...
                             fig.page.ypos(yind)-fig.subplot.height/2,...
                             fig.subplot.width*1.9,...
                             fig.subplot.height],...
                 'FontSize', 8,...
                 'LineWidth',1);
boxplot([prefBhvSlopes;nprefBhvSlopes],      ...
        [prefBhvId*2-1;nprefBhvId*2],              ...
        'plotstyle',    'traditional',       ...
        'boxstyle',     'filled',            ...
        'colors',       'rrggmmbbcc',             ...
        'symbol',       '.',                 ...
        'datalim',      [-10,7.5],           ...
        'labels',       {});
grid('on');
%%%>>>

%%%<<< CA1 Non selective population phase precession statistics
%ppZscores = sq(bsxfun(@minus,rhoHRZall(:,1,1,:),mean(rhoHRZall(:,1,2:end,:),3))./std(rhoHRZall(:,1,2:end,:),[],3));
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
    ind =   ismember(cluSessionMap,cluSessionMapSubset_N,'rows') ...
          & ismember(cluSessionMap(:,1),sesIds)           ...
          & rHRZall(:,1,1,[s+1]) > rthresh;
    prefBhvId     = cat(1,prefBhvId,s.*ones([sum(ind),1]));    
    prefBhvSlopes = cat(1,prefBhvSlopes,parmHRZall(ind,1,1,s+1));        
    prefBhvPhases = cat(1,prefBhvSlopes,parmHRZall(ind,2,1,s+1));
    nprefBhvId     = cat(1,nprefBhvId,(statesIndCompPhz(s)-1).*ones([sum(ind),1]));
    nprefBhvSlopes = cat(1,nprefBhvSlopes,parmHRZall(ind,1,1,statesIndCompPhz(s)));
    nprefBhvPhases = cat(1,nprefBhvSlopes,parmHRZall(ind,2,1,statesIndCompPhz(s)));    
end            
% CREATE subplot axes
xind = 1;    
yind = 4;
sax(end+1) = axes('Units','centimeters',                                                         ...
                 'Position',[fig.page.xpos(xind),                                                ...
                             fig.page.ypos(yind)-fig.subplot.height/2,                           ...
                             fig.subplot.width*1.9,                                              ...
                             fig.subplot.height],                                                ...
                 'FontSize', 8,                                                                  ...
                 'LineWidth',1);
boxplot([prefBhvSlopes;nprefBhvSlopes],      ...
        [prefBhvId*2-1;nprefBhvId*2],              ...
        'plotstyle',    'traditional',       ...
        'boxstyle',     'filled',            ...
        'colors',       'rrggmmbbcc',             ...
        'symbol',       '.',                 ...
        'datalim',      [-10,7.5],           ...
        'labelorientation','inline',       ...
        'labels',       reshape([statesLabels(2:end);repmat({''},[1,numel(statesLabels)-1])],[],1)');        
grid('on');
sax(end).Position =[fig.page.xpos(xind),fig.page.ypos(yind)-fig.subplot.height/2,fig.subplot.width*1.9,fig.subplot.height];
%%%>>>
%%%>>>

%%%<<< GHZ VS PHZ ratemap examples 

statesIndsGPE = [2,3,4];
expUnitsGPE = {[20],[21];...
               [21],[14];...               
               [ 3],[158];
               [21],[22];...
               [22],[58];...
               [20],[103]};

markers = '^vposhp';
tppUnits = {[5], [4];  ...  REAR , high @ peak    
            [21],[22]; ...  HIGH , low  @ peak 
            [22],[61]; ...  HIGH , low  @ peak 
            [18],[18]; ...  LOW  , high @ peak
            [22],[58]; ...  LOW  , high @ trough
            [20],[103];...  LOW
            [19],[15]; ...  HIGH
           };

yind = 5;                
ucnt = 1 
for tind = 1:size(tppUnits,1)
    t = tppUnits{tind,1}(1);
    for uid = 1:numel(tppUnits{tind,2});
        yind = yind+1;                
        u = tppUnits{tind,2}(uid);
        uind =  ismember(cluSessionMapSubset,[t,u],'rows');

        %%%<<< PLOT place field
        xind = 1;
        sax(end+1) = axes('Units','centimeters',                                             ...
                          'Position',[fig.page.xpos(xind),                                   ...
                                      fig.page.ypos(yind)-fig.subplot.height/2,              ...
                                      fig.subplot.width,                                     ...
                                      fig.subplot.height],                                   ...
                          'FontSize', 8,                                                     ...
                          'LineWidth',1);
        hold(sax(end),'on');
        plot(pfts{t},u,'mean',false,[],true,0.5,false,interpParPfs,colorMap,[],nanColor);
        text(-480,-380,num2str(cond_round(pfts{t}.maxRate(u,true,'interpPar',interpParPfs))),...
             'FontSize',10,'Color',[1,1,1]);
        sax(end).XTickLabels = {};
        sax(end).YTickLabels = {};
        % plot state placefield centers 
        pfssCntrMrate = zeros([numel(statesPfssInd),1]);
        pfssCntrPos = zeros([numel(statesPfssInd),2]);
        for s = 1:numel(statesPfssInd),
            [pfssCntrMrate(s),pfssCntrPos(s,:)] = maxRate(pfss{t}{statesPfssInd(s)},...
                                                          u,true,'interpPar',interpParPfs);
            if s~=1 & any(ismember(round(pfssCntrPos(1:s-1,:)),round(pfssCntrPos(s,:)),'rows')),
                shift = 10;
            else
                shift = 0;
            end
            
            if pfssCntrMrate(s)>2,
                hax = scatter(pfssCntrPos(s,1)+shift,pfssCntrPos(s,2),15,sclr(s));
            end
            
        end

        axis(sax(end),'tight');
        axes(fax);
        rectangle('Position',sax(end).Position,'LineWidth',1);
        %%%>>>

        %%%<<< PLOT behavior field
        xind = 2;
        sax(end+1) = axes('Units','centimeters',                                ...
                          'Position',[fig.page.xpos(xind),                      ...
                                      fig.page.ypos(yind)-fig.subplot.height/2, ...
                                      fig.subplot.width,                        ...
                                      fig.subplot.height],                      ...
                          'FontSize', 8,                                        ...
                          'LineWidth',1);
        plot(pfbs{t},u,'mean',[],[],false,0.5,false,...
             interpParDfs,@cool,reshape(validDims,pfbs{t}.adata.binSizes'),nanColor);    
        text(-1.7,-0.45,num2str(cond_round(pfbs{t}.maxRate(u,false,'mask',validDims))),...
             'FontSize',10,'Color',[1,1,1]);
        %text(-1.7,-0.45,num2str(cond_round(maxPfsRate)),'FontSize',10,'Color',[1,1,1]);
        sax(end).XTickLabels = {};
        sax(end).YTickLabels = {};
        xlim([-1.8,0.5]);
        ylim([-0.75,1.8]);
        axes(fax);
        rectangle('Position',sax(end).Position,'LineWidth',1);
        %%%>>>
        
        %%%<<< PLOT drzphz field
        mrateHZTPD = max(cell2mat(cf(@(p) p{t}.maxRate(u,false), pftHZTPD(statesIndsGPE(:)))));
        for s = 1:numel(statesIndsGPE),
% SET horizontal offset    
            xind = 2+s;    
% CREATE subplot axes
            sax(end+1) = axes('Units','centimeters',                                             ...
                              'Position',[fig.page.xpos(xind),                                   ...
                                          fig.page.ypos(yind)-fig.subplot.height/2,              ...
                                          fig.subplot.width,                                     ...
                                          fig.subplot.height/2],                                 ...
                              'FontSize', 8,                                                     ...
                              'LineWidth',1);
            
            hold('on');                            
            plot(pftHZTPD{statesIndsGPE(s)}{t},u,'mean','',[0,mrateHZTPD],false,[],[],[],colorMap,[],nanColor);
            text(-.9,-1.8,num2str(cond_round(mrateHZTPD)),...
                 'FontSize',10,'Color',[1,1,1]);
            if mean(maxPhzRate(:,uind,s)) > 4,
                plot(mprPos(:,uind,s)',            ... x = drz 
                     pftHZTPD{1}{1}.adata.bins{2}',... y = phase
                     ['-k'],                       ... LineColor
                     'LineWidth',1);                 % LineWidth
                plot(meanDTPos(:,uind,1,s)',         ... x = drz
                     pftHZTPD{1}{1}.adata.bins{2}',... y = phase
                     ['-w'],                       ... LineColor
                     'LineWidth',1);                 % LineWidth
                               
                ylim([-pi,pi]);
            end

% $$$             [hl,hp] = boundedline(mprPos(:,uind,s)',...
% $$$                         pftHZTPD{1}{1}.adata.bins{2}',...
% $$$                         0,...
% $$$                         ...% mprPosBSStd(:,ismember(cluSessionMapSubset,[t,u],'rows'),s)'.*1.96./10,...
% $$$                         ['-',sclr(s)],...
% $$$                         'alpha',...                        
% $$$                         'transparency',0.3,...
% $$$                         'orientation','horiz');
% $$$             hl.LineWidth = 1;
% $$$             hp.EdgeColor = 'none';
            Lines(0,[],'w');
            
            sax(end+1) = axes('Units','centimeters',                                             ...
                              'Position',[fig.page.xpos(xind),                                   ...
                                          fig.page.ypos(yind),                                   ...
                                          fig.subplot.width,                                     ...
                                          fig.subplot.height/2],                                 ...
                              'FontSize', 8,                                                     ...
                              'LineWidth',1);
            hold('on');
            plot(pftHZTPD{statesIndsGPE(s)}{t},u,'mean','',[0,mrateHZTPD],false,[],[],[],colorMap,[],nanColor);
            if mean(maxPhzRate(:,uind,s)) > 4,
                plot(mprPos(:,uind,s)',            ... x = drz
                     pftHZTPD{1}{1}.adata.bins{2}',... y = phase
                     ['-k'],                       ... LineColor
                     'LineWidth',1);                 % LineWidth
                plot(meanDTPos(:,uind,1,s)',       ... x = drz
                     pftHZTPD{1}{1}.adata.bins{2}',... y = phase
                     ['-w'],                       ... LineColor
                     'LineWidth',1);                 % LineWidth
                ylim([-pi,pi]);                               
            end
            
% $$$             [hl,hp] = boundedline(mprPos(:,uind,s)',...
% $$$                         pftHZTPD{1}{1}.adata.bins{2}',...
% $$$                         0,...
% $$$                         ... %mprPosBSStd(:,uind,s)'.*1.96./10,...
% $$$                         ['-',sclr(s)],...
% $$$                         'alpha',...
% $$$                         'transparency',0.3,...
% $$$                         'orientation','horiz');
% $$$             hl.LineWidth = 1;            
% $$$             hp.EdgeColor = 'none';
            Lines(0,[],'w',1);            
            
            if s == numel(statesIndsGPE),
                sax(end).YAxisLocation = 'right';
                sax(end).YTick = [0,pi,2*pi,3*pi];
                sax(end).YTickLabels = [0,180,360,540];
                sax(end).XTickLabels = {};   
            else        
                sax(end).XTickLabels = {};
                sax(end).YTickLabels = {};
            end    
            sax(end).Color = 'none';
            sax(end-1).Color = 'none';            

            

            axes(fax);
            rectangle('Position',[sax(end-1).Position.*[1,1,1,2]],'LineWidth',1);
            
        end%for s
        %%%>>>
        
        %%%<<< PLOT drz corrected theta phase preference
        sax(end+1) = axes('Units','centimeters',                                             ...
                          'Position',[fig.page.xpos(xind+1),                                 ...
                                      fig.page.ypos(yind)-fig.subplot.height/2,              ...
                                      fig.subplot.width,                                     ...
                                      fig.subplot.height],                                   ...
                          'FontSize', 8,                                                     ...
                          'LineWidth',1);
        hold('on');
        for s = 1:numel(statesIndsGPE),
            [hl,hp] = boundedline(repmat(maxPhzRate(:,uind,s),[2,1])',...
                        pbins',...
                        repmat(phzRateStd(:,uind,s),[2,1])',...
                        ['-',sclr(s)],...
                        'alpha',...
                        'transparency',0.3,...
                        'orientation','horiz');
            hl.LineWidth = 1;
            hp.EdgeColor = 'none';
            [hl,hp] = boundedline(repmat(phzRateMean(:,uind,s),[2,1])',...
                        pbins',...
                        repmat(phzRateStd(:,uind,s),[2,1])'.*1.96./10,...
                        ['-',sclr(s)],...
                        'transparency',1,...
                        'orientation','horiz');
            hl.LineWidth = 1;
            hp.EdgeColor = 'none';
        end
        plot(30,360,markers(ucnt),'MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10);
        sax(end).YAxisLocation = 'right';
        sax(end).XTickLabels = {};
        sax(end).YTick = [0,180,360,540];
        sax(end).YTickLabels = [0,180,360,540];
        
        ylim(pbins([1,end]));
        %xlim([0,max(xlim)]);
        %xlim([0,mrateHZTPD+8]);
        xlim([0,25]);
        %%%>>>
        
        ucnt = 1 + ucnt;
    end%for uid
end%for tind
%%%>>>

%%%<<< PLOT phase eigenvectors 
xind = 7;
yind = 6;
sax(end+1) = axes('Units','centimeters',                                             ...
                  'Position',[fig.page.xpos(xind)+fig.subplot.width/2,               ...
                              fig.page.ypos(yind)-fig.subplot.height/2,              ...
                              fig.subplot.width,                                     ...
                              fig.subplot.height],                                   ...
                  'FontSize', 8,                                                     ...
                  'LineWidth',1);
hold('on');
for s = 1:3,
    plot(repmat(W(:,s),2,1),pbins,'LineWidth',2);
end
ylim([-180,540]);
sax(end).YAxisLocation = 'right';
sax(end).XTickLabels = {};
sax(end).YTick = [0,180,360,540];
sax(end).YTickLabels = [0,180,360,540];


yind = 8;
sax(end+1) = axes('Units','centimeters',                                             ...
                  'Position',[fig.page.xpos(xind)+fig.subplot.width/2,               ...
                              fig.page.ypos(yind)-fig.subplot.height/2,              ...
                              fig.subplot.width*2,                                   ...
                              fig.subplot.height*2],                                 ...
                  'FontSize', 8,                                                     ...
                  'LineWidth',1);
hold('on');
%%%>>>

%%%<<< PLOT rate diff vs phase diff
figure,hold('on');
ucnt = 1;
% EXAMPLES 
for tind = 1:size(tppUnits,1)
    t = tppUnits{tind,1}(1);
    dispOpts = {['-',markers(ucnt)],                              ...
                'MarkerFaceColor','m',                            ...
                'MarkerEdgeColor','m',                            ...
                'MarkerSize',10};
    for uid = 1:numel(tppUnits{tind,2});
        u = tppUnits{tind,2}(uid);
        ind = ismember(cluSessionMapSubset,[t,u],'rows');% example unit index        
        plot(prmRateSrt(ind,1)'-prmRateSrt(ind,2)',       ...
             -circ_dist(prmPhzAngSrt(ind,1)',                     ...
                        prmPhzAngSrt(ind,2)'),                    ...
         dispOpts{:});
    end
    ucnt = ucnt+1;
end
% POPULATION 
smarkers = '^sv';
for s = 1:3
    ind = dsi(:,1)==s & validPosOccRlx & sesInds(unitSubset) & prmMaxRateSrt(:,2) > 4;
    scatter(prmRateSrt(ind,1)'-prmRateSrt(ind,2)',            ...
            -circ_dist(prmPhzAngSrt(ind,1)',                          ...
                       prmPhzAngSrt(ind,2)'),                         ...
            15,prmPhzAngSrtShift(ind),smarkers(s),'filled');
end    
Lines([],0,'k');
Lines(0,[],'k');    
colormap(hsv(3));
caxis([0,2*pi]);
xlim([0,15]);
Lines(4.5,[],'k');
% MEAN phz diff
rbins = linspace(0,15,8);
rdfBinInd = discretize(prmRateSrt(:,1)-prmRateSrt(:,2),rbins);
nind = nniz(rdfBinInd) & phzBinInd==2 & validPosOccRlx & sesInds(unitSubset) & prmMaxRateSrt(:,2) > 4;
out = accumarray(rdfBinInd(nind),                     ... Subs
                 circ_dist(prmPhzAngSrt(nind,1),  ... Vals
                           prmPhzAngSrt(nind,2)), ...
                 [numel(rbins)-1,1],                  ... Size
                 @median);%                               Fun
plot(diff(rbins)./2+rbins(1:end-1),out,'-+');
%%%>>>


figure();
for s = 1:3;
    subplot(1,3,s);
    hold('on');
    for j = nonzeros(double(s~=[1:3]).*[1:3])',
        ind =  dsi(:,1)==s & dsi(:,2)==j & validPosOccRlx & sesInds(unitSubset) & prmMaxRateSrt(:,2) > 4 & validUnits;
        plot(prmPhzAngSrt(ind,1)+double(prmPhzAngSrt(ind,1)<0).*2*pi,...
                prmPhzAngSrt(ind,2)+double(prmPhzAngSrt(ind,2)<0).*2.*pi,smarkers(j));
    end
    line([0,2*pi],[0,2*pi])
end




phzBinInd = discretize(prmPhzAngSrtShift,linspace(0,2*pi,4));
rdfBinInd = discretize(prmRateSrt(:,1)-prmRateSrt(:,2),[0,4.5,30]);
ssc = 'mc';
phzBinsR = linspace(-pi,pi,13);
figure();hold('on');
for r = 1:2,
    ind = phzBinInd==2 & rdfBinInd==r & validPosOccRlx & sesInds(unitSubset) & prmMaxRateSrt(:,2) > 4;
    hax = bar(phzBinsR,                                                       ...
              histc(circ_dist(prmPhzAngSrt(ind,1)',                       ...
                              prmPhzAngSrt(ind,2)'),phzBinsR),            ...
              'histc');
    hax.FaceColor = ssc(r);
    hax.EdgeColor = ssc(r);
    hax.FaceAlpha = 0.4;
    hax.EdgeAlpha = 0;
    Lines(mean(circ_dist(prmPhzAngSrt(ind,1)',prmPhzAngSrt(ind,2)')),[],ssc(r));
end    

figure();
s = 3;
phzBinsR = linspace(-pi,pi,17);
ind = dsi(:,1)==s &phzBinInd==2 & validPosOccRlx & sesInds(unitSubset) & prmMaxRateSrt(:,2) > 4;
hist(circ_dist(prmPhzAngSrt(ind,1)',prmPhzAngSrt(ind,2)'),phzBinsR);

figure();
phzBinsR = linspace(-pi,pi,17);
ind = validPosOccRlx & sesInds(unitSubset) & prmMaxRateSrt(:,2) > 4;
hist(circ_dist(prmPhzAngSrt(ind,1)',prmPhzAngSrt(ind,2)'),phzBinsR);


sesIds = [3:5,8:12,17:23];
sesInds = ismember(cluSessionMap(:,1),sesIds);

figure,hist(prmPhzRlnSrt(ind,1),0:0.01:0.5)



% ENDFIG




%% NNMF analysis figures ---------------------------------------------------------------------------

% PLOT fscr diff of Dominate vs subdominate state
yind = 8;
sax(end+1) = axes('Units','centimeters',                                             ...
                  'Position',[fig.page.xpos(xind)+fig.subplot.width/2,               ...
                              fig.page.ypos(yind)-fig.subplot.height/2,              ...
                              fig.subplot.width*2,                                   ...
                              fig.subplot.height*2],                                 ...
                  'FontSize', 8,                                                     ...
                  'LineWidth',1);
hold('on');

yind = 10;
sax(end+1) = axes('Units','centimeters',                                             ...
                  'Position',[fig.page.xpos(xind)+fig.subplot.width/2,               ...
                              fig.page.ypos(yind)-fig.subplot.height/2,              ...
                              fig.subplot.width*2,                                   ...
                              fig.subplot.height*2],                                 ...
                  'FontSize', 8,                                                     ...
                  'LineWidth',1);
hold('on');


% PLOT the differences between phase features of dominate ves subdominate 
validPosOccRlx =  all(sq(sum(any(isnan(pftTPZDa(10:30,:,:,1,3:4)),2)))<5,2);
validPosOccRr =  all(sq(sum(any(isnan(pftTPZDa(10:30,:,:,1,2)),2)))<5,2);

for s = 1:3
    ind = dsi(:,1)==s & validPosOccRlx & sesInds(unitSubset);
    dfscr = bsxfun(@minus,fscrCPPDAll(ind,~ismember(1:size(fscrCPPDAll,2),s),:),fscrCPPDAll(ind,s,:));
    dispOpts = {'.',...
                'MarkerFaceColor',stsColor(s),...
                'MarkerEdgeColor',stsColor(s),...
                'MarkerSize',3};
    for j = 1:size(dfscr,1),
        if dsi(j,2)>dsi(j,3)
            axes(sax(end-1));
            plot(dfscr(j,2,1),dfscr(j,2,2),dispOpts{:});
            plot(dfscr(j,1,1),dfscr(j,1,2),dispOpts{:});
            axes(sax(end));                
            plot(dfscr(j,2,1),dfscr(j,2,3),dispOpts{:});
            plot(dfscr(j,1,1),dfscr(j,1,3),dispOpts{:});
        else
            axes(sax(end-1));            
            plot(dfscr(j,1,1),dfscr(j,1,2),dispOpts{:});
            plot(dfscr(j,2,1),dfscr(j,2,2),dispOpts{:});
            axes(sax(end));            
            plot(dfscr(j,1,1),dfscr(j,1,3),dispOpts{:});
            plot(dfscr(j,2,1),dfscr(j,2,3),dispOpts{:});
        end
    end
end
axes(sax(end));            
Lines([],0,'k','-',1);
Lines(0,[],'k','-',1);
axes(sax(end-1));            
Lines([],0,'k','-',1);
Lines(0,[],'k','-',1);
            
grid(sax(end),'on');
grid(sax(end-1),'on');
ylim(sax(end),[-0.09,0.045]);
xlim(sax(end),[-0.12,0.045]);
ylim(sax(end-1),[-0.09,0.045]);
xlim(sax(end-1),[-0.12,0.045]);
sax(end-1).XTickLabels = {};
sax(end-1).YTickLabels = {};
sax(end).XTickLabels = {};
sax(end).YTickLabels = {};

tppUnits = {[20],[21];...
               [21],[14];...               
               [21],[22];...
               [22],[58];...
               [20],[103]};

%cla(sax(end))
%cla(sax(end-1))
ucnt = 1 
for tind = 1:size(tppUnits,1)
    t = tppUnits{tind,1}(1);
    dispOpts = {['-',markers(ucnt)],...
                'MarkerFaceColor','m',...
                'MarkerEdgeColor','m',...
                'MarkerSize',5};
    
    for uid = 1:numel(tppUnits{tind,2});
        u = tppUnits{tind,2}(uid);
        ind = ismember(cluSessionMapSubset,[t,u],'rows');% example unit index
        dfscr = bsxfun(@minus,...
                       fscrCPPDAll(ind,~ismember(1:size(fscrCPPDAll,2),dsi(ind,1)),:),...
                       fscrCPPDAll(ind,dsi(ind,1),:));
        if dsi(ind,2)>dsi(ind,3)
            axes(sax(end-1));
            plot(dfscr(1,:,1),dfscr(1,:,2),dispOpts{:});
            %plot(dfscr(1,2,1),dfscr(1,2,2),dispOpts{:});
            %plot(dfscr(1,1,1),dfscr(1,1,2),dispOpts{:});
            axes(sax(end));                
            plot(dfscr(1,:,1),dfscr(1,:,3),dispOpts{:});
            %plot(dfscr(1,2,1),dfscr(1,2,3),dispOpts{:});
            %plot(dfscr(1,1,1),dfscr(1,1,3),dispOpts{:});
        else
            axes(sax(end-1));            
            plot(dfscr(1,:,1),dfscr(1,:,2),dispOpts{:});
            %plot(dfscr(1,1,1),dfscr(1,1,2),dispOpts{:});
            %plot(dfscr(1,2,1),dfscr(1,2,2),dispOpts{:});
            axes(sax(end));            
            plot(dfscr(1,:,1),dfscr(1,:,3),dispOpts{:});
            %plot(dfscr(1,1,1),dfscr(1,1,3),dispOpts{:});
            %plot(dfscr(1,2,1),dfscr(1,2,3),dispOpts{:});
        end
        ucnt = 1 + ucnt;
    end        
end

        





%cla(sax(end))
%cla(sax(end-1))





% $$$ sp = tight_subplot(2,8,0,0.1);
% $$$ sp = reshape(reshape(sp',8,2)',2,8);
% $$$ for u = units{t},
% $$$     for s = 1:8,
% $$$         axes(sp(s*2-1));plot(pftHZTPD{s}{t},u,'mean','text',[],false);
% $$$         axes(sp(s*2));plot(pftHZTPD{s}{t},u,'mean','text',[],false);
% $$$         title(num2str(u));
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end




% $$$ % CA3
% $$$ xind = 1;    
% $$$ yind = 7;
% $$$ sax(end+1) = axes('Units','centimeters',...
% $$$                  'Position',[fig.page.xpos(xind),fig.page.ypos(yind)-fig.subplot.height/2,fig.subplot.width*1.9,fig.subplot.height],...
% $$$                  'FontSize', 8,...
% $$$                  'LineWidth',1);
% $$$ sesIds = [1,2,6,7,13:16];
% $$$ % ACCUMULATE slope within units' preferred behavior
% $$$ prefBhvId = [];
% $$$ prefBhvSlopes = [];
% $$$ prefBhvPhases = [];
% $$$ nprefBhvId = [];
% $$$ nprefBhvSlopes = [];
% $$$ nprefBhvPhases = [];
% $$$ for s = 1:5;    
% $$$     ind =   ismember(cluSessionMap,statesUSubs{s},'rows') ...
% $$$           & ismember(cluSessionMap(:,1),sesIds)           ...
% $$$           & rHRZall(:,1,1,[s+1])   > rthresh;
% $$$     %          & ppZscores(:,[s+1]) < zthresh                  ...    
% $$$     prefBhvId = cat(1,prefBhvId,s.*ones([sum(ind),1]));    
% $$$     prefBhvSlopes = cat(1,prefBhvSlopes,parmHRZall(ind,1,1,s+1));        
% $$$     prefBhvPhases = cat(1,prefBhvSlopes,parmHRZall(ind,2,1,s+1));
% $$$     nprefBhvId = cat(1,nprefBhvId,(statesIndCompPhz(s)-1).*ones([sum(ind),1]));
% $$$     nprefBhvSlopes = cat(1,nprefBhvSlopes,parmHRZall(ind,1,1,statesIndCompPhz(s)));
% $$$     nprefBhvPhases = cat(1,nprefBhvSlopes,parmHRZall(ind,2,1,statesIndCompPhz(s)));    
% $$$ end            
% $$$ boxplot([prefBhvSlopes;nprefBhvSlopes],      ...
% $$$         [prefBhvId*2-1;nprefBhvId*2],        ...
% $$$         'plotstyle',    'traditional',       ...
% $$$         'boxstyle',     'filled',            ...
% $$$         'colors',       'rrbbccggmm',        ...
% $$$         'symbol',       '.',                 ...
% $$$         'datalim',      [-10,7.5],           ...
% $$$         'labelorientation','inline',       ...
% $$$         'labels',       reshape([statesLabels(2:end);repmat({''},[1,numel(statesLabels)-1])],[],1)');        
% $$$ %        'labels',       reshape(repmat(statesLabels(2:end),[2,1]),[],1)');
% $$$ grid('on');
% $$$ sax(end).Position =[fig.page.xpos(xind),fig.page.ypos(yind)-fig.subplot.height/2,fig.subplot.width*1.9,fig.subplot.height];


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
% $$$ uind = ismember(cluSessionMap,cluSessionMap(unitSubset,:),'rows');
% $$$ for s = 1:numel(statesInds);
% $$$     xind = 2+s;
% $$$ % CREATE subplot axes
% $$$     sax(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[fig.page.xpos(xind),fig.page.ypos(yind+3)-fig.subplot.height/2,fig.subplot.width,fig.subplot.height],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$     plot(RmaxN_H(uind,s+1),PN_H(uind,s+1,1),'.')
% $$$     %scatter(RmaxN_H(uind,s),PN_H(uind,s,1),4,bsxfun(@minus,FSrC(uind(unitSubset),[2,3,1]),[0,0.2,0.2]),'Filled');
% $$$     xlim([0,1]);
% $$$     ylim([-15,15]);
% $$$     Lines([],0,'k');
% $$$     if s==1, ylabel({'entering','field'}); end    
% $$$ end
% $$$ 
% $$$ for s = 1:numel(statesInds);
% $$$     xind = 2+s;
% $$$ % CREATE subplot axes
% $$$     sax(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[fig.page.xpos(xind),fig.page.ypos(yind+4)-fig.subplot.height/2,fig.subplot.width,fig.subplot.height],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$     %scatter(RmaxP_H(uind,s),PP_H(uind,s,1),4,bsxfun(@minus,FSrC(uind(unitSubset),[2,3,1]),[0,0.2,0.2]),'Filled');    
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
% $$$     sax(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[fig.page.xpos(xind),fig.page.ypos(yind+5)-fig.subplot.height/2,fig.subplot.width,fig.subplot.height],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$     ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows');          
% $$$     hold('on');
% $$$ % $$$     plot(linspace(hrzBins(1),hrzBins(end),numel(hrzBins)-1),...
% $$$ % $$$          accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_mean)','+');
% $$$     errorbar(linspace(hrzBins(1),hrzBins(end),numel(hrzBins)-1)',...
% $$$              accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_mean)',...
% $$$              accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_var)'./2);    
% $$$     sax(end).XTickLabels = {};
% $$$     sax(end).YTickLabels = {};
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
% $$$     sax(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[fig.page.xpos(xind),fig.page.ypos(yind+5)-fig.subplot.height/2,fig.subplot.width,fig.subplot.height],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$     ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows');          
% $$$     hold('on');
% $$$ % $$$     plot(linspace(hrzBins(1),hrzBins(end),numel(hrzBins)-1),...
% $$$ % $$$          accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_mean)','+');
% $$$     errorbar(linspace(hrzBins(1),hrzBins(end),numel(hrzBins)-1)',...
% $$$              accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_mean)',...
% $$$              accumarray(discretize(spkhrz(ind),hrzBins), spkphz(ind),[numel(hrzBins)-1,1],@circ_var)'./2);    
% $$$     sax(end).XTickLabels = {};
% $$$     sax(end).YTickLabels = {};
% $$$     title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
% $$$     Lines([],pi,'k');    
% $$$ end





    
% $$$     sax(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[fig.page.xpos(xind),fig.page.ypos(yind+2)-fig.subplot.height/2,fig.subplot.width,fig.subplot.height],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$     ind = tind & any(spkstc(:,statesIndComp{s}),2) & ismember(spkmap,statesUSubs{s},'rows');
% $$$     hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
% $$$           linspace(-1,1,21),...
% $$$           linspace(-pi,3*pi,21));    
% $$$     sax(end).XTickLabels = {};
% $$$     sax(end).YTickLabels = {};
% $$$     %caxis([0,150]); 
% $$$     Lines([],pi,'k');    
% $$$     if s == 1, ylabel({'Complementary','behavior'}); end    
% $$$     if s == numel(statesInds),
% $$$         sax(end).YAxisLocation = 'right';
% $$$         sax(end).YTick = [0,pi,2*pi,3*pi];
% $$$         sax(end).YTickLabels = [0,180,360,540];
% $$$     else        
% $$$         sax(end).XTickLabels = {};
% $$$         sax(end).YTickLabels = {};
% $$$     end
% $$$     
% $$$     
% $$$     sax(end+1) = axes('Units','centimeters',...
% $$$                      'Position',[fig.page.xpos(xind),fig.page.ypos(yind+4),fig.subplot.width,fig.subplot.height],...
% $$$                      'FontSize', 8,...
% $$$                      'LineWidth',1);
% $$$     ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMapSubset_N,'rows');
% $$$     hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
% $$$           linspace(-1,1,21),...
% $$$           linspace(-pi,3*pi,21));    
% $$$     if s == numel(statesInds),
% $$$         sax(end).YAxisLocation = 'right';
% $$$         sax(end).YTick = [0,pi,2*pi,3*pi];
% $$$         sax(end).YTickLabels = [0,180,360,540];
% $$$     else        
% $$$         sax(end).XTickLabels = {};
% $$$         sax(end).YTickLabels = {};
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
    sax(end+1) = axes('Units','centimeters',...
                     'Position',[fig.page.xpos(xind),fig.page.ypos(yind+1)-fig.subplot.height/2,fig.subplot.width,fig.subplot.height],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,statesUSubs{s},'rows') & abs(spkego(:,2)<100);
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));
    sax(end).XTickLabels = {};
    sax(end).YTickLabels = {};
    title(['N = ',num2str(numel(unique(spkmap(ind,:),'rows'))/2)]);
    if s == 1, ylabel({'Behavior','Specific Units'}); end
    if s == numel(statesInds),
        sax(end).YAxisLocation = 'right';
        sax(end).YTick = [0,pi,2*pi,3*pi];
        sax(end).YTickLabels = [0,180,360,540];
        sax(end).XTickLabels = {};        
    else        
        sax(end).XTickLabels = {};
        sax(end).YTickLabels = {};
    end    
    Lines([],pi,'k');    

% NONPREFFERED state    
    sax(end+1) = axes('Units','centimeters',...
                     'Position',[fig.page.xpos(xind),fig.page.ypos(yind+2)-fig.subplot.height/2,fig.subplot.width,fig.subplot.height],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesIndComp{s}),2) & ismember(spkmap,statesUSubs{s},'rows') & abs(spkego(:,2)<100);
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));    
    sax(end).XTickLabels = {};
    sax(end).YTickLabels = {};
    Lines([],pi,'k');    
    if s == 1, ylabel({'Complementary','behavior'}); end    
    if s == numel(statesInds),
        sax(end).YAxisLocation = 'right';
        sax(end).YTick = [0,pi,2*pi,3*pi];
        sax(end).YTickLabels = [0,180,360,540];
    else        
        sax(end).XTickLabels = {};
        sax(end).YTickLabels = {};
    end
    
% THEREST state        
    sax(end+1) = axes('Units','centimeters',...
                     'Position',[fig.page.xpos(xind),fig.page.ypos(yind+4),fig.subplot.width,fig.subplot.height],...
                     'FontSize', 8,...
                     'LineWidth',1);
    ind = tind & any(spkstc(:,statesInds(s)),2) & ismember(spkmap,cluSessionMapSubset_N,'rows') & abs(spkego(:,2)<100);
    hist2([repmat(spkhrz(ind),[2,1]),reshape(bsxfun(@plus,repmat(spkphz(ind),[1,2]),[0,2*pi]),[],1)],...
          linspace(-1,1,21),...
          linspace(-pi,3*pi,21));    
    if s == numel(statesInds),
        sax(end).YAxisLocation = 'right';
        sax(end).YTick = [0,pi,2*pi,3*pi];
        sax(end).YTickLabels = [0,180,360,540];
    else        
        sax(end).XTickLabels = {};
        sax(end).YTickLabels = {};
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

% $$$ statesUSubs  = {cluSessionMap(rScore>-0.5 ,:), ...
% $$$                 cluSessionMap(hlScore<-0.1 & rScore<-0.2,:), ...
% $$$                 cluSessionMap(hlScore<-0.1 & rScore<-0.2,:), ...
% $$$                 cluSessionMap(hlScore>0.1  & rScore<-0.2,:), ...
% $$$                 cluSessionMap(hlScore>0.1  & rScore<-0.2,:)};

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

% $$$ uind = ismember(cluSessionMap,cluSessionMap(unitSubset,:),'rows') & ismember(cluSessionMap(:,1),[17:23]);
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
% $$$ uind = ismember(cluSessionMap,cluSessionMap(unitSubset,:),'rows');
% $$$ uind = ismember(cluSessionMap,cluSessionMap(unitSubset,:),'rows') & ismember(cluSessionMap(:,1),[17:23]);
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


% DEF Conditional Expectation: req20190527(): phase VS gdz (gaussian distance zone)

figure,
sp = tight_subplot(4,4,0,0);
sp = reshape(sp,4,4);
for u = 1:size(cluSessionMapRestricted,1),
    axes(sp(2,1));
    imagesc(sq(fsrcz(u,:,:)));
    title(num2str(cluSessionMapRestricted(u,:)));
    colorbar();
    axes(sp(1,1));
    plot(phaseBinCenters,LR(phzOrder,:));
    xlim([0,2*pi])
    title('erpPCA');
    axes(sp(1,2));
    plot(phaseBinCenters,sq(pftTPZap(phzOrder,u,:)));
    xlim([0,2*pi])
    legend({'rear','high','low'});
    t = cluSessionMapRestricted(u,1);
    unit = cluSessionMapRestricted(u,2);
    mrate = prctile(nonzeros(plot(pftHZTPD{s}{t},unit,1,'text',[],false)),99).*1.5;
    for s = 1:4,
        try,
            axes(sp(s,3));
            plot(pftHZTPD{s}{t},unit,1,'text',[0,mrate],false);
            title([num2str(u) '-' stateLabels{s} '-' num2str(mrate)]);
            axes(sp(s,4));        
            plot(pftHZTPD{s}{t},unit,1,'text',[0,mrate],false);
        end
    end
    waitforbuttonpress();
    cla(sp);
end


%%%<<< DIAGNOSTIC figures 
numComp = 3;
fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan(size(validDims{1}));
    fpc{i}(validDims{1}) = eigVecs{1}(:,i);
end

fpcMinMax = [min(cellfun(@min,fpc)),max(cellfun(@max,fpc))];
figure();
yind = yind + 1;
for i = 1:3,
    subplot(1,3,i)    
    imagescnan({pfd{1}.adata.bins{:},abs(reshape_eigen_vector(fpc{i},pfd(1,1)))},...
               fpcMinMax,'linear',false,nanColor,1,1);                % PRINT eigenvectors
    axis('xy');
    axis('tight');
    ylabel(['F',num2str(i)])
    xlim(pfd{1}.adata.bins{1}([1,end]));
    xlim([-2,nonzeros(xlim.*[0,1])-0.2])            
    ylim(pfd{1}.adata.bins{2}([1,end]));
end

%%%>>>