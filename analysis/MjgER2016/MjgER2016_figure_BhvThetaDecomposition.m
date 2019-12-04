 g% MjgER2016_figure5
% Behavior Theta Decomposition

% DECOMPOSITION of behavior space across theta phases
%
% SUBPLOT : examples of individual units
%     ELEMENTS : { pft(1), bhv(1), PhzDecomp(6), mTPZ(1), rTPZ(1), PCs(2), ScrPhz(1),BhvPhz(1) }
%
% SUBPLOT : distribution of mean phase for all single factor units
%
% SUBPLOT : distribution of mean phase for all double factor units
%     ELEMENTS : { mPhzPC1, mPhzPC2 }
%
% SUBPLOT : inter-factor angle vs factors' score ratio
%







% FIGURE : Behavioral theta phase decomposition Examples -------------------------------------------%

% bhvMulti  mod 52, 60, 61, 80,
% bhvSingle mod 81, 86, 103, 104 
% bhvProg   mod 63, 34, 83, 72, 145

MjgER2016_load_data();
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


% SET analysis args
MjgER2016_figure5_args('section 1');

overwrite = false;

% DEF place field rate maps
[pfts]          = cf(@(t,u) pfs_2d_theta(t,u), Trials,units);

% DEF behavior field rate maps
brm        = cf(@(t,u)   compute_bhv_ratemaps(t,u),                      Trials, units);
brmShuff   = cf(@(t,u)   compute_bhv_ratemaps_shuffled(t,u),             Trials, units);
[pfds]                = req20180123_ver5(Trials);
[pfss,metaData]       = req20181106(Trials,[],[],units,'overwrite',false);

% COMPUTE bhv ratemap nnmf
[eigVecs,eigScrs,eigVars,unitSubset,validDims] = compute_bhv_ratemaps_erpPCA(brm,units,'overwrite',overwrite);
% $$$ [eigVecs,eigScrs,dScale, unitSubset,validDims] = compute_bhv_ratemaps_nnmf(  brm,units,'overwrite',overwrite);
% $$$ [fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS] = ...
% $$$         compute_bhv_ratemaps_nnmf_scores(Trials,units,brm,brmShuff,eigVecs,validDims,unitSubset);
[fsrcz,FSrC,rmaps,FSCFr,FSrM,FSrS,fsrsMean,fsrsStd] = ...
    compute_bhv_ratemaps_erpPCA_scores(Trials,units,brm,brmShuff,eigVecs,validDims,unitSubset,overwrite);


% $$$ evec = nan([brm{1}.adata.binSizes']);
% $$$ figure
% $$$ for j = 1:size(eigVecs,2),
% $$$     subplot(1,size(eigVecs,2),j);
% $$$     evec(validDims) = eigVecs(:,j);
% $$$     imagesc(brm{1}.adata.bins{:},evec');
% $$$     axis('xy');
% $$$ end


    
% ANALYSIS theta phase rate maps -------------------------------------------------------------------


stateLabels = {'theta','rear','high','low'}
stateColors = 'brgk';

% DEF theta phase rate maps for each state
pftTPZ = cf(@(s) cf(@(T,u) MTAApfs(T,u,'tag',['tp-',s]), Trials, units), stateLabels);

figure();
phzOrder = [9:16,1:8];
t = 21
rmb = pftTPZ{1}{1}.adata.bins{1}(phzOrder)+2*pi*double(pftTPZ{1}{1}.adata.bins{1}(phzOrder)<0);
for u = units{t};
hold('on');    
for s = 1:4,
    %plot(pftTPZ{s}{t}.adata.bins{:},plot(pftTPZ{s}{t},u),stateColors(s),'LineWidth',1);
    %plot(pftTPZ{s}{t}.adata.bins{:},plot(pftTPZ{s}{t},u,'mean_circ_smooth'),stateColors(s),'LineWidth',1);
    rmp = plot(pftTPZ{s}{t},u,'mean_circ_smooth');
    %rmp = plot(pftTPZ{s}{t},u,'mean');    
    plot(rmb,rmp(phzOrder,:),stateColors(s),'LineWidth',1);    
    
end
xlim([0,2*pi])
title(num2str(u));
waitforbuttonpress();
clf();
end


% DEF shuffled theta phase rate maps for each state
pftTPZshuff = cf(@(s) cf(@(T,u) MTAApfs(T,u,'tag',['tp-shuff-',s]), Trials, units), stateLabels);

figure();
phzOrder = [9:16,1:8];
t = 21
rmb = pftTPZshuff{1}{1}.adata.bins{1}(phzOrder)+2*pi*double(pftTPZshuff{1}{1}.adata.bins{1}(phzOrder)<0);
for u = units{t};
hold('on');    
for s = 1:4,
    %plot(pftTPZ{s}{t}.adata.bins{:},plot(pftTPZ{s}{t},u),stateColors(s),'LineWidth',1);
    %plot(pftTPZ{s}{t}.adata.bins{:},plot(pftTPZ{s}{t},u,'mean_circ_smooth'),stateColors(s),'LineWidth',1);
    rmp = plot(pftTPZshuff{s}{t},u,'mean_circ_smooth');
    %rmp = plot(pftTPZ{s}{t},u,'mean');    
    plot(rmb,rmp(phzOrder,:),stateColors(s),'LineWidth',1);    
    
end
xlim([0,2*pi])
title(num2str(u));
waitforbuttonpress();
clf();
end
 
%---------------------------------------------------------------------------------------------------



% $$$ % HRZ vs PHZ - egocentric rate zones versus theta phaze given behavioral state and placefield proximity
% $$$ % req20190126();
% $$$ pftHZTP = cf(@(s) cf(@(T,u) MTAApfs(T,u,'tag',['hztp-',s]), Trials, units), stateLabels);
% $$$ figure,
% $$$ t = 19;
% $$$ sp = tight_subplot(2,4,0,0.05);
% $$$ sp = reshape(sp,4,2)';
% $$$ for u = units{t},
% $$$     mrate = prctile(nonzeros(plot(pftHZTP{s}{t},u,1,'text',[],false)),99);
% $$$     for s = 1:4,
% $$$         try,
% $$$         axes(sp(1,s));
% $$$         plot(pftHZTP{s}{t},u,1,'text',[0,mrate],false);
% $$$         title([num2str(u) '-' stateLabels{s} '-' num2str(mrate)]);
% $$$         axes(sp(2,s));
% $$$         plot(pftHZTP{s}{t},u,1,'text',[0,mrate],false);
% $$$         end
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end

% $$$ % ddr vs PHZ - egocentric rate zones versus theta phaze given behavioral state and placefield proximity
% $$$ % req20190127();
% $$$ pftHZTP = cf(@(s) cf(@(T,u) MTAApfs(T,u,'tag',['ddtp-',s]), Trials, units), stateLabels);
% $$$ figure,
% $$$ t = 20;
% $$$ sp = tight_subplot(2,4,0,0.05);
% $$$ sp = reshape(sp,4,2)';
% $$$ for u = units{t},
% $$$     mrate = prctile(nonzeros(plot(pftHZTP{s}{t},u,1,'text',[],false)),99).*1.5;
% $$$     for s = 1:4,
% $$$         try,
% $$$         axes(sp(1,s));
% $$$         plot(pftHZTP{s}{t},u,1,'text',[0,mrate],false);
% $$$         title([num2str(u) '-' stateLabels{s} '-' num2str(mrate)]);
% $$$         axes(sp(2,s));
% $$$         plot(pftHZTP{s}{t},u,1,'text',[0,mrate],false);
% $$$         end
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end




sesIds = [3:5,8:12,17:23];
sesIds = [8:12,17:23];
sesInds = ismember(cluSessionMap(:,1),sesIds);
uind = diag(rhoB)>0.5&diag(rhoS)>0.5; % SOURCE req20190

cluSessionMapSubset = cluSessionMap(unitSubset,:);
cluSessionMapRestricted = cluSessionMapSubset(uind&sesInds(unitSubset),:);

pftTPZap = pftTPZa(:,uind&sesInds(unitSubset),2:4);
npftTPZa = reshape(pftTPZap,16,[]);
%npftTPZa = reshape(permute(pftTPZap,[2,1,3]),[],16.*3);
npftTPZa(~nniz(npftTPZa(:))) = 0;
nind = sum(npftTPZa==0)==16;
% $$$ nind = sum(npftTPZa==0,2)==16*3;
%npftTPZa(nind,:) =[];
npftTPZa(:,nind) =[];

[LU,LR,FSr,VT] = erpPCA(npftTPZa',3);
fscr = zeros([numel(nind),3]);
fscr(~nind,:) = FSr;
fscr = reshape(fscr,[],3,3);


figure,plot(LR)
figure,plot(VT(1:10,4),'-+');

[W,H,D] = nnmf(npftTPZa,3);
Ha = zeros([numel(nind),3]);
Ha(~nind,:) = H';
Ha = reshape(Ha,[],3,3);

figure();plot(W)



figure,
for u = 1:size(cluSessionMapRestricted,1),
subplot(321);
plot(pftTPZ{1}{1}.adata.bins{1},LR);
title('erpPCA');
subplot(323);
imagesc(sq(fscr(u,:,:)));
title(num2str(cluSessionMapRestricted(u,:)));
colorbar();
subplot(322);
plot(pftTPZ{1}{1}.adata.bins{1},W);
title('nnmf');
subplot(324);
imagesc(sq(Ha(u,:,:)));
title(num2str(cluSessionMapRestricted(u,:)));
colorbar();
subplot(3,2,[5,6]);
plot(pftTPZ{1}{1}.adata.bins{1},sq(pftTPZap(:,u,:)));
legend({'rear','high','low'});
waitforbuttonpress();
clf();
end



figure,
plot(fscr(:,2,1),fscr(:,3,3),'.');
grid('on');
line([-3,3],[-3,3]);




% REDUCE placefield objects to 3D tensor
pftTPZa = cf(@(pfs) cf(@(p,u) mean(RectFilter(p.data.rateMap(:,ismember(p.data.clu,u),:),3,1,'circular'),3,'omitnan'),  pfs,units), pftTPZ);
clu     = cf(@(pfs) cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),                        pfs,units), pftTPZ);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
pftTPZa = cf(@(pfs) cat(2,pfs{:}), pftTPZa);
clu     = cf(@(c)   cat(2,c{:}),clu);
assert(all(cell2mat(cf(@(c) all(c), cf(@(c,o) c==o, clu,circshift(clu,1,2))))),'Clus do not match');
pftTPZa = cat(3,pftTPZa{:});
clu = clu{1};
tlu=cat(2,tlu{:});
clu=[tlu',clu'];
[clu,rind] = sortrows(clu);
clu     = clu(unitSubset,:);
pftTPZa = pftTPZa(:,rind,:);
pftTPZa = pftTPZa(:,unitSubset,:);
pftTPZa(isnan(pftTPZa)) = 0;


% COMPUTE mean theta phase for each state
meanStateThetaPhase = pftTPZa;
pftTPZcpx = sq(sum(reshape(bsxfun(@times,                                                     ...
                               reshape(pftTPZa,size(pftTPZa,1),[]),                        ...
                               exp(-i.*pftTPZ{1}{1}.adata.bins{1})),             ...
                        [size(pftTPZa,1),size(pftTPZa,2),size(pftTPZa,3)]),1)              ...
            ./sum(pftTPZa,1));

uind = diag(rhoB)>0.5&diag(rhoS)>0.5;
figure,
tang = angle(pftTPZcpx(:,1));
for s = 2:4,
    for j = 2:4,    
        if s~=j,
            subplot2(4,4,s,j);            
            hold('on');
            %rose(angle(pftTPZcpx(abs(pftTPZcpx(:,s))>0.1,s)),16);
            %histcirc(angle(pftTPZcpx(abs(pftTPZcpx(:,s))>0.1,s)),16,pi/2)        
            cpx = [pftTPZcpx(:,s),pftTPZcpx(:,j)];
            aind = all(abs(cpx)>0.1,2);
            dmax = [max(pftTPZa(:,:,s))',max(pftTPZa(:,:,j))'];
            
            ind = all(dmax>2,2);
            aind = aind&ind&uind&sesInds(unitSubset);

% $$$             scatter(circ_dist(tang(aind),angle(cpx(aind,1))),...
% $$$                     circ_dist(tang(aind),angle(cpx(aind,2))),20,diff(dmax(aind,:),1,2),'filled')
            
            
            scatter(angle(cpx(aind,1)),angle(cpx(aind,2)),20,-diff(dmax(aind,:),1,2),'filled')
            scatter(angle(cpx(aind,1))+2*pi,angle(cpx(aind,2))+2*pi,20,-diff(dmax(aind,:),1,2),'filled')
            scatter(angle(cpx(aind,1))+2*pi,angle(cpx(aind,2)),20,-diff(dmax(aind,:),1,2),'filled')
            scatter(angle(cpx(aind,1)),angle(cpx(aind,2))+2*pi,20,-diff(dmax(aind,:),1,2),'filled')            
            
% $$$             scatter(circ_dist(angle(cpx(aind,1)),angle(cpx(aind,2))),angle(cpx(aind,1)),20,diff(dmax(aind,:),1,2),'filled')
% $$$             scatter(circ_dist(angle(cpx(aind,1)),angle(cpx(aind,2))),angle(cpx(aind,1))+2*pi,20,diff(dmax(aind,:),1,2),'filled')

% $$$             scatter(diff(dmax(aind,:),1,2),angle(cpx(aind,1)),20,circ_dist(angle(cpx(aind,1)),angle(cpx(aind,2))),'filled')
% $$$             scatter(diff(dmax(aind,:),1,2),angle(cpx(aind,1))+2*pi,20,circ_dist(angle(cpx(aind,1)),angle(cpx(aind,2))),'filled')
            Lines([],pi,'k');
            Lines(pi,[],'k');
% $$$             Lines([],0,'k');
% $$$             Lines(0,[],'k');
            colormap('jet');            
            caxis([-10,10]);
            %colormap('hsv');            
% $$$             xlim([-pi,pi]);
% $$$             ylim([-pi,pi]);
            %caxis([-pi,pi]);
            ylim([0,2*pi]);
            xlim([0,2*pi]);
            line([0,2*pi],[0,2*pi]);

            title([stateLabels{s} ' ' stateLabels{j}]);

        end
    end
end



    
t = 20;
Trial = Trials{t};
unitSubset = units{t};

% LOAD placefields and theta restricted behavior fields
[pft]          = pfs_2d_theta(Trial,unitSubset);
[pfd]          = req20180123_ver5(Trial);
[pfs,metaData] = req20181106(Trial,[],[],unitSubset);
%brmNormal     = compute_bhv_ratemaps(Trial,unitSubset);

% LOAD group factor analysis of theta-space restricted behavior fields
[~,tags,eigVec,eigVar,eigScore,validDims,...
 unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(Trials,'loadPFDFlag',false);

% LOAD group factor analysis of theta-space theta decomposition of each unit
[eigVecsFact, eigScrsFact, eigVarsFact,eSpiFact, FSrBhvThpFact,metaDataFact]= req20181119('overwrite',overwrite);
[eigVecsNnmf, eigScrsNnmf,eigVarsNnmf, eSpi, FSrBhvThp, metaData] = req20181119_nnmf('overwrite',overwrite);

figure,
for u = unitSubsets{1},
    for v = 1:3,
        subplot(2,4,v),
        imagesc(sq(eigVecs(u,v,:,:))'),axis('xy');caxis(max(abs(reshape(eigVecs(u,v,:,:),[],1))).*[-1,1]);
        subplot(2,4,v+4),    
        imagesc(sq(eigVecsNnmf(u,v,:,:))'),axis('xy');
    end
    subplot(2,4,4),
    plot(eigVars(u,1:3),'+');
    subplot(2,4,8),
    plot(eigVarsNnmf(u,1:3).*100,'+');
    waitforbuttonpress();
end


% LOAD Factor analysis scores
MjgER2016_load_bhv_erpPCA_scores();
%     rmaps    - matrix[D x S](numeric); rate maps corresponding to the valid eigenvector dims
%     FSrC     - matrix[U x V](Numeric); fscores
%     fsrcz    - matrix[U x V](Numeric); normalized fscores
clear(varNonEss{:},varNonAux{:});


% GET behavioral information scores
bsi = cf(@(p,u) p.data.si(:,ismember(p.data.clu,u),:), pfds(:,1),units');
clu =  cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:), pfds(:,1),units');
tlu =  cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
bsi = cat(2, bsi{:});
clu = cat(2, clu{:});
tlu = cat(2, tlu{:});
clu = [tlu',clu'];
[~,rind] = sortrows(clu);
bsi = bsi(:,rind);


rmap = plot(pfs,unit,1,'',[],false);
rmapNind = sq(sum(reshape(~isnan(rmap(:)),size(rmap)),2))>=numPhzBins;
zmap = reshape(permute(rmap,[1,3,2]),[],numPhzBins);
zmap(repmat(~rmapNind(:),[1,numPhzBins])) = 0;            
zmap(zmap(:)<0) = 0;
rmapNind = rmapNind&reshape(sum(reshape(zmap,[],numPhzBins)','omitnan')~=0,...
                            size(rmap,1),size(rmap,3));
mrate = mean(reshape(zmap(:,p),[],1),'omitnan');
pmap = zmap(:,p);
pmap = pmap(rmapNind);

bsi = sum(rmaps./mean(rmaps,'omitnan')...
                               .*log2(rmaps./mean(rmaps,'omitnan')),'omitnan')./size(rmaps,1);


%bsi(isnan(bsi)) = 0;


 %req20190122();% -> bhvInfoTS
cluExamples = [17, 32]
 
cluExamples = [18, 20;... 
               18, 21;...
               18, 25;...               
               18, 63];
               
cluExamples = [19,  3;...
               19, 13;...
               19, 50;...
               19, 66;...               
               19, 67;...
               19,117;...
               19,141];

cluExamples = [21, 14;... 
               21, 22;... 
               21, 61;... 
               21, 63];
cluExamples = [22, 13;... 
               22, 18;...                
               22, 19;...                               
               22, 41;...                                              
               22, 50;...                                                             
               22, 61;...
               23, 67;...               
               23, 69];            
 

% $$$ unitsExampleSubset = [34, 83, 63, ... % rearing transitions
% $$$                       52, 61, 72, ... % nonselective transitions
% $$$                       81, 86,104];    % selective stationary
% $$$ unitsExampleSubset = [34, 83, 63,119, ... % rearing transitions
% $$$                       110, 79, 111,80, ... % nonselective transitions
% $$$                       35,139,104,103];    % selective stationary
% $$$ unitsExampleSubset = [34, 83, 63, ... % rearing transitions
% $$$                        79,80, ... % nonselective transitions
% $$$                       35,139,104,103];    % selective stationary

cluExampleSet  = [18, 24; ... all walk high low phase shift
                  18, 25; ... all walk high low phase dropout
                  18, 33; ... low walk
                  18, 42; ... all to low walk
                  18, 49; ... mid walk
                  18, 60; ... all sorts of crazy but only in walk 
                  ...
                  20, 34; ... rear 
                  20, 35; ... 
                  20, 63; ...
                  20, 79; ...                  
                  20, 80; ... % nonselective transitions
                  20, 83; ...                  
                  20,103; ...
                  20,104; ...              
                  20,110; ...                  
                  20,111; ...                  
                  20,119; ... % rearing transitions
                  20,139; ...
                  ...                  
                  21,  6; ... low walk dominant
                  21, 22; ... all to high walk
                  21, 32; ... mid walk
                  21, 37; ... low walk
                  ...
                  22,  6; ... rear
                  22, 11; ... rear
                  22, 13; ... walk slight spatial shift between high and low
                  22, 18; ... low walk
                  22, 22; ... low walk                  
                  22, 61; ... high walk
                  ...
                  23, 18; ... high walk                  
                  23, 35; ... double placefield high and low walk
                  23, 36; ... high walk                  
                  23, 40; ... walk high and low phase shift
                  23, 63; ... high walk
                  23, 70; ... all walk                                    
                  23, 71; ... low walk                  
                  23, 69];





cluExampleSet  = [18, 24; ...
                  18, 25; ...
                  19, 67; ...
                  19, 141; ...                            
                  20, 80; ...
                  20, 83; ...
                  21, 63; ...
                  22, 19; ...
                  22, 50; ...
                  23, 67; ...
                  23, 69];

                  
cluExampleSet = [20, 34; ...
                 20, 83; ...
                 20, 63; ...
                 20,119; ... % rearing transitions
                 20,110; ...
                 20, 79; ...
                 20,111; ...
                 20, 80; ... % nonselective transitions
                 20, 35; ...
                 20,139; ...
                 20,104; ...
                 20,103];    % selective stationary




nx = pfss{1}.adata.binSizes(2)/2+4;
%ny = numel(unitsExampleSubset);
ny = size(cluExampleSet,1);


cond_round = @(rate) max([round(rate,0),round(rate,1)].*[rate>=10,rate<10]);
nanColor = [0.25,0.25,0.25];


% SETUP Main figure
[hfig,fig,fax] = set_figure_layout(figure(666006),'A4','portrait');


sp = gobjects([1,0]);

% BEGIN Main figure
phzOrder = [9:16,1:8];
%phzOrder = [10,12,14,16,2,4,6,8];
phzLables = {'','90','','180','','270','',''};
clear('i');
for y = 1:ny,
    yind = y;

% GET space and theta phase jointly restricted bhv rate maps
    rmap = plot(pfss{cluExampleSet(y,1)},cluExampleSet(y,2),1,'colorbar',[],false);
    rmapNind = sq(sum(reshape(~isnan(rmap(:)),size(rmap)),2))>=pfs.adata.binSizes(2);
    zmap = reshape(permute(rmap,[1,3,2]),[],pfs.adata.binSizes(2));
    zmap(repmat(~rmapNind(:),[1,pfs.adata.binSizes(2)])) = 0;            
    zmap(zmap(:)<0) = 0;            
    rmapNind = rmapNind&reshape(sum(reshape(zmap,[],pfs.adata.binSizes(2))','omitnan')~=0,...
                                size(rmap,1),size(rmap,3));      
    mrate= prctile(rmap(:),99.9);    
    
% PLOT theta restricted spatial rate maps
    xind = 1;
    sp(end+1) = axes('Units','centimeters',...
                 'Position',[fig.page.xpos(xind),fig.page.ypos(yind),fig.subplot.width,fig.subplot.height],...
                 'FontSize', 8);
    plot(pfts{cluExampleSet(y,1)},cluExampleSet(y,2),'mean','text',mrate,'nanColor',nanColor,'colorMap',@jet);
    %ylabel(['Unit: ',num2str(unitsExampleSubset(y))]);
    ylabel([num2str(cluExampleSet(y,:))]);
    %sp(end).YLabel.Rotation = 0;
    sp(end).YTickLabels ={};
    sp(end).XTickLabels ={};    
    if y == 1,
        title({'PFS'});
    end
    
% PLOT spatially restricted bhv rate maps
    xind = 2;
    sp(end+1) = axes('Units','centimeters',...
                 'Position',[fig.page.xpos(xind),fig.page.ypos(yind),fig.subplot.width,fig.subplot.height],...
                 'FontSize', 8);
    plot(pfds{cluExampleSet(y,1)},cluExampleSet(y,2),...
         'mean','text',mrate,false,'nanColor',nanColor,'colorMap',@jet);
    set(sp(end),'YTickLabels',{});
    set(sp(end),'XTickLabels',{});    
    sp(end).Position(1) = sp(end).Position(1)+0.1;    
    if y == 1,
% $$$         ylabel('Body Pitch');
% $$$         xlabel('Head Pitch');
        title({'BHV'});
    end
    
% PLOT space and theta phase jointly restricted bhv rate maps
    for p = 1:(pfs.adata.binSizes(2)/2),
        xind = 2+p;        
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[fig.page.xpos(xind),fig.page.ypos(yind),fig.subplot.width,fig.subplot.height],...
                         'FontSize', 8);
        hold(sp(end),'on');
        imagescnan({pfd{1}.adata.bins{:},sq(sum(rmap(:,phzOrder([p*2-1,p*2]),:),2)./2)'},...
                   [0,mrate],[],false,nanColor,'colorMap',@jet);
        sp(end).Position(1) = sp(end).Position(1)+0.2;
        axis('xy');
        axis('tight');
        set(sp(end),'YTickLabels',{});
        set(sp(end),'XTickLabels',{});    
    end 
    text(pfs.adata.bins{1}(end)-0.45*diff(pfs.adata.bins{1}([1,end])),...
         pfs.adata.bins{3}(end)-0.10*diff(pfs.adata.bins{3}([1,end])),...
         sprintf('%2.0f',mrate),...
        'Color','w','FontWeight','bold','FontSize',8)
   
   
% COMPUTE circular statistics of rate changes within the theta cycle
% MEAN theta phase map within behavior space
    mask = double(sq(max(rmap,[],2))>2)';
    mask(mask==0) = nan;
    

% RESULTANT length of theta phase dependent rate changes
   xind = 4+p;        
   sp(end+1) = axes('Units','centimeters',                                                      ...
                    'Position',[fig.page.xpos(xind),                                                     ...
                                fig.page.ypos(yind),                                                     ...
                                fig.subplot.width,                                              ...
                                fig.subplot.height],                                            ...
                    'FontSize', 8);
   hold(sp(end),'on');
   imagescnan({pfd{1}.adata.bins{:},abs(rmapComplex)'.*mask},...
              [0,0.5],[],false,nanColor,'colorMap',@copper);
   sp(end).Position(1) = sp(end).Position(1)+0.3;
   axis('xy');
   axis('tight');
   set(sp(end),'YTickLabels',{});
   set(sp(end),'XTickLabels',{});    
   if y==1, title('r'); end

   

% GET unit row within population data
   uind = find(ismember(cluSessionMap,cluExampleSet(y,:),'rows'));   
   for e = 1:2
       if eigVars(uind,e)>15,
           xind = 4+p+e;        
% SETUP subplot
           sp(end+1) = axes('Units','centimeters',                                             ...
                            'Position',[fig.page.xpos(xind),                                            ...
                               fig.page.ypos(yind),                                                     ...
                               fig.subplot.width,                                              ...
                               fig.subplot.height],                                            ...
                            'FontSize', 8);
           hold(sp(end),'on');
% PLOT eigenvector
           eVecLims = repmat(max(abs([min(reshape(sq(eigVecsNnmf(uind,e,:,:)),[],1)),...
                                      max(reshape(sq(eigVecsNnmf(uind,e,:,:)),[],1))])),[1,2]).*[-1,1];
           imagescnan({pfd{1}.adata.bins{:},sq(eigVecsNnmf(uind,e,:,:))'},...
                      [eVecLims],[],false,nanColor,'colorMap',@parula);
% $$$            eVecLims = repmat(max(abs([min(reshape(sq(eigVecs(uind,e,:,:)),[],1)),...
% $$$                                       max(reshape(sq(eigVecs(uind,e,:,:)),[],1))])),[1,2]).*[-1,1];
% $$$            imagescnan({pfd{1}.adata.bins{:},sq(eigVecs(uind,e,:,:))'},...
% $$$                       [eVecLims],[],false,nanColor,'colorMap',@parula);
           sp(end).Position(1) = sp(end).Position(1)+0.4;
           axis('xy');
           axis('tight');
           set(sp(end),'YTickLabels',{});
           set(sp(end),'XTickLabels',{});    
           
           text(pfs.adata.bins{1}(end)-0.45*diff(pfs.adata.bins{1}([1,end])),...
                pfs.adata.bins{3}(end)-0.10*diff(pfs.adata.bins{3}([1,end])),...
                sprintf('%2.0f',eigVarsNnmf(uind,e)*100),...
                'Color','w','FontWeight','bold','FontSize',8)
% $$$            text(pfs.adata.bins{1}(end)-0.45*diff(pfs.adata.bins{1}([1,end])),...
% $$$                 pfs.adata.bins{3}(end)-0.10*diff(pfs.adata.bins{3}([1,end])),...
% $$$                 sprintf('%2.0f',eigVars(uind,e)),...
% $$$                 'Color','w','FontWeight','bold','FontSize',8)
       else
           sp(end+1) = gobjects([1,1]);
       end
   end%for e
end%for y


% CREATE reference bars
axes(fax)
xlim([0,fax.Position(3)]);
ylim([0,fax.Position(4)]);

% SUPERAXIS Theta phase 
supAxisThpOffset = -(ny+1)*fig.subplot.height-(ny).*fig.subplot.verticalPadding+0.7;
% $$$ line([sp(3).Position(1),sum(sp(nx).Position([1,3]))],...
% $$$      sum(sp(3).Position([2,4])).*[1,1]+0.5,'LineWidth',1)
text(mean([sp(3).Position(1),sum(sp(nx-2).Position([1,3]))]),...
     sum(sp(3).Position([2,4]))-0.1+supAxisThpOffset,'Theta Phase','HorizontalAlignment','center','FontSize',8)
% TICKS theta phase
text(sp(3).Position(1),sum(sp(3).Position([2,4]))+0.25+supAxisThpOffset,'0^{o}','HorizontalAlignment','center','FontSize',8)
text(mean([sp(3).Position(1),sum(sp(nx-2).Position([1,3]))]),...
     sum(sp(3).Position([2,4]))+0.25+supAxisThpOffset,'180^{o}','HorizontalAlignment','center','FontSize',8)
text(sum(sp(nx-2).Position([1,3])),...
     sum(sp(3).Position([2,4]))+0.25+supAxisThpOffset,'360^{o}','HorizontalAlignment','center','FontSize',8)

% ADD phase colorbar to the top of the examples
hcb = colorbar(sp(3),'NorthOutside');
hcb.Units = 'centimeters';
hcb.Position(2) = sp(1).Position(2) + fig.subplot.height + 0.5 + supAxisThpOffset;
hcb.Position(3) = fig.subplot.width*round(numel(phzOrder)/2);
colormap(hcb,'hsv');
caxis([0,360]);
hcb.YTick = [];

% REFERENCE bar behavior space
line([sp(end-nx).Position(1),sum(sp(end-nx).Position([1,3]).*[1,11/28])],...
     sp(end-nx).Position(2).*[1,1]-0.1,'LineWidth',1);
text(sp(end-nx).Position(1),sp(end-nx).Position(2)-0.3,'1 rad','FontSize',8);

% REFERENCE bar physical space
line([sp(end-nx-1).Position(1),sum(sp(end-nx-1).Position([1,3]).*[1,0.5])],...
     sp(end-nx-1).Position(2).*[1,1]-0.1,'LineWidth',1);
text(sp(end-nx-1).Position(1),sp(end-nx-1).Position(2)-0.3,'50 cm','FontSize',8);


% END FIGURE 5 ------------------------------------------------------------------------------------%







for t = 1:numel(Trials),
    for u = 1:numel(units{t}),
        rmap = plot(pfss{t},units{t}(u),1,'',[],false);
        rmap = plot(pfss{9},20,1,'',[],false);        
        rmapNind = sq(sum(reshape(~isnan(rmap(:)),size(rmap)),2)) >= pfs.adata.binSizes(2);
        zmap = reshape(permute(rmap,[1,3,2]),[],pfs.adata.binSizes(2));
        zmap(repmat(~rmapNind(:),[1,pfs.adata.binSizes(2)])) = 0;
        zmap(zmap(:)<0) = 0;            
        
    end
end

[W,H,D] = nnmf(zmap,5);

figure,
for w = 1:5,
subplot(1,6,w);
imagesc(reshape(W(:,w),[28,28])'),axis('xy');
colorbar
end
subplot(1,6,w+1);
plot(H(1:3,:)')
Lines(8.5,[],'k');

sum(bsxfun(@lt,dist,[300,200,100,50]),2)


rmapComplex = zeros([28,28,0]);
rmapComplex = sum(reshape(bsxfun(@times,                                                     ...
                                 reshape(permute(rmap,[1,3,2]),[],size(rmap,2)),             ...
                                 exp(-i.*pfs.adata.bins{2}(phzOrder))'),                     ...
                          [size(rmap,1),size(rmap,3),size(rmap,2)]),                         ...
                  3)                                                                         ...
              ./sq(sum(rmap,2));



mask = nan(pfd{1}.adata.binSizes');
mask(validDims{1}) = 1;

nEigVecs = eigVecsNnmf./repmat(max(max(abs(eigVecsNnmf),[],4),[],3),[1,1,pfs.adata.binSizes([1,3])']);
mnEigVecs = mask;
%eigVecsRng = 
pfig = figure(23032903);
clf();
for e = 1:numel(ind),    
    subplot2(2,numel(ind),1,e);
% $$$     plot(eigVars(ind{e},1),sum(reshape(nEigVecs(ind{e},1,:,:),sum(ind{e}),[]),2,'omitnan')...
% $$$          ./sqrt(sum(reshape(nEigVecs(ind{e},1,:,:),sum(ind{e}),[]).^2,2,'omitnan')),'.')
%    plot(eigVars(ind{e},1),...

    plot(eigVarsNnmf(ind{e},1),...
         max(reshape(eigVecs(ind{e},1,:,:),sum(ind{e}),[]),[],2),'.');
% $$$ 
% $$$     plot(eigVarsNnmf(ind{e},1),...
% $$$          log10(sum(reshape(nEigVecs(ind{e},1,:,:),sum(ind{e}),[])>0,2,'omitnan')...
% $$$          ./sum(reshape(nEigVecs(ind{e},1,:,:),sum(ind{e}),[])<0,2,'omitnan')),'.');
% $$$ 

% $$$     plot(sum(reshape(nEigVecs(ind{e},1,:,:),sum(ind{e}),[]),2,'omitnan')...
% $$$           ./sqrt(sum(reshape(nEigVecs(ind{e},1,:,:),sum(ind{e}),[]).^2,2,'omitnan')),...
% $$$            clip(prctile(reshape(nEigVecs(ind{e},1,:,:),sum(ind{e}),[]),98,2) ...
% $$$          ./prctile(reshape(nEigVecs(ind{e},1,:,:),sum(ind{e}),[]), 2,2),-30,30),'.');
% $$$ 
% $$$     plot(sum(reshape(nEigVecs(ind{e},1,:,:),sum(ind{e}),[]),2,'omitnan')...
% $$$           ./sqrt(sum(reshape(nEigVecs(ind{e},1,:,:),sum(ind{e}),[]).^2,2,'omitnan')),...
% $$$            clip(prctile(reshape(nEigVecs(ind{e},1,:,:),sum(ind{e}),[]),98,2) ...
% $$$          ./prctile(reshape(nEigVecs(ind{e},1,:,:),sum(ind{e}),[]), 2,2),-30,30),'.');
    
    %xlim([0,100]);
    %ylim([-31,31]);    
    %ylim([-15,15]);
% $$$     subplot2(2,numel(ind),2,e);    
% $$$     plot(eigVars(ind{e},1),sum(eigVars(ind{e},[1,2]),2),'.')
end

mfig = figure(2939203);
zfig = figure(2390294);
efig = figure(2039040);
chax = gobjects([1,1]);
while 1,
    figure(pfig);   
    waitforbuttonpress();    
    currentAxis = gca();
    cp = currentAxis.CurrentPoint(1,[1,2]);
    currentAxisIndex = find(flipud(ismember(get(pfig,'Children'),currentAxis)));
    currentCluMap = cluSessionMap(ind{currentAxisIndex},:);
    xydata = [currentAxis.Children(1).XData',currentAxis.Children(1).YData'];
    [~,cmin] = min(sqrt(sum(bsxfun(@minus,cp,xydata).^2,2)));
    currentClu = currentCluMap(cmin,:);
% $$$     hold(currentAxis,'on');    
% $$$     delete(chax);
% $$$     chax = circle(xydata(cmin,1),xydata(cmin,2),10);
% $$$     hold(currentAxis,'off');    
    
    figure(mfig);
    clf();
    plot(pfds{currentClu(1),1},currentClu(2),1,'text',[],false,[],[],[],[],mask);

    figure(zfig);
    clf();
    rmps = reshape(pfss{currentClu(1)}.data.rateMap(:,pfss{currentClu(1)}.data.clu==currentClu(2),1),...
                   pfss{1}.adata.binSizes');
    mrate = prctile(rmps(:),99);
    for s = 1:16,
        subplot(1,16,s);
        imagescnan({pfss{currentClu(1)}.adata.bins{1},...
                    pfss{currentClu(1)}.adata.bins{3},...
                    sq(rmps(:,s,:))'},[0,mrate],[],s==16);
        axis('xy');
    end
    
    figure(efig);
    clf();
    for s = 1:4,
        subplot(1,4,s);
        if s~=4,
        imagescnan({pfss{1}.adata.bins{1},...
                    pfss{1}.adata.bins{3},...
                    sq(eigVecsNnmf(ismember(cluSessionMap,currentClu,'rows'),s,:,:))'},[],[],true);
        axis('xy');
        else
            plot(eigVarsNnmf(ismember(cluSessionMap,currentClu,'rows'),1:3),'-+');
        end
        
    end
end

%fin










% COLLATE set of units base on behavioral scores
thresholdVariance =  20;
sesIds = [3:5,8:12,17:23];
sesIds = [8:12,17:23];
sesInds = ismember(cluSessionMap(:,1),sesIds);
reqInds = ismember(cluSessionMap,cluSessionMap(unitSubsets{1},:),'rows');
indR = ismember(cluSessionMap,cluSessionMapSubset_R,'rows');
indH = ismember(cluSessionMap,cluSessionMapSubset_H,'rows');
indL = ismember(cluSessionMap,cluSessionMapSubset_L,'rows');
indN = ismember(cluSessionMap,cluSessionMapSubset_N,'rows');
indC = ismember(cluSessionMap,cluSessionMapSubset_C,'rows');
eigenVectorVarianceRatios = [sum(  eigVars(indR&sesInds&reqInds,1)>thresholdVariance                    ...
                                 & eigVars(indR&sesInds&reqInds,2)>=thresholdVariance),                 ...
                             sum(  eigVars(indR&sesInds&reqInds,1)>thresholdVariance                    ...
                                 & eigVars(indR&sesInds&reqInds,2)<thresholdVariance);                  ...
                                                                                                ... 
                             sum(  eigVars(indH&sesInds&reqInds,1)>thresholdVariance                    ...
                                 & eigVars(indH&sesInds&reqInds,2)>=thresholdVariance),                 ...
                             sum(  eigVars(indH&sesInds&reqInds,1)>thresholdVariance                    ...
                                 & eigVars(indH&sesInds&reqInds,2)<thresholdVariance);                  ...
                                                                                                ...
                             sum(  eigVars(indL&sesInds&reqInds,1)>thresholdVariance                    ...
                                 & eigVars(indL&sesInds&reqInds,2)>thresholdVariance),                  ...
                             sum(  eigVars(indL&sesInds&reqInds,1)>thresholdVariance                    ...
                                 & eigVars(indL&sesInds&reqInds,2)<thresholdVariance);                  ...
                                                                                                ... 
                             sum(  eigVars(indN&sesInds&reqInds,1)>thresholdVariance                    ...
                                 & eigVars(indN&sesInds&reqInds,2)>thresholdVariance),                  ...
                             sum(  eigVars(indN&sesInds&reqInds,1)>thresholdVariance                    ...
                                 & eigVars(indN&sesInds&reqInds,2)<thresholdVariance)];



% INDEX something
validInds = ismember(cluSessionMap,cluSessionMapSubset,'rows');
reqInds = ismember(cluSessionMap,cluSessionMap(unitSubsets{1},:),'rows');
indS ={indR,indH,indL,indN};
ind = cf(@(i) i&reqInds&sesInds,indS);



indUnitFsrcz = cf(@(c)  ismember(cluSessionMapSubset(:,1),sesIds)  ...
                        &ismember(cluSessionMapSubset,c,'rows'),cluSessionMapSubsets);

rmapsMax = max(rmaps);
figure,
hold('on');
subplot(131);
scatter(bsi(sesInds(reqInds)),fsrcz(sesInds(reqInds),1),20,rmapsMax(sesInds(reqInds)),'filled');
colormap('jet');caxis([0,10]);title('LLoc');
subplot(132);
scatter(bsi(sesInds(reqInds)),fsrcz(sesInds(reqInds),2),20,rmapsMax(sesInds(reqInds)),'filled');
colormap('jet');caxis([0,10]);title('Rear');
subplot(133);
scatter(bsi(sesInds(reqInds)),fsrcz(sesInds(reqInds),3),20,rmapsMax(sesInds(reqInds)),'filled');
colormap('jet');caxis([0,10]);title('HLoc');

rmapsMax = max(rmaps);
figure,
hold('on');
subplot(131);
scatter(bsi(sesInds(reqInds)),FSrC(sesInds(reqInds),1),20,rmapsMax(sesInds(reqInds)),'filled');
colormap('jet');caxis([0,10]);title('LLoc');
subplot(132);
scatter(bsi(sesInds(reqInds)),FSrC(sesInds(reqInds),2),20,rmapsMax(sesInds(reqInds)),'filled');
colormap('jet');caxis([0,10]);title('Rear');
subplot(133);
scatter(bsi(sesInds(reqInds)),FSrC(sesInds(reqInds),3),20,rmapsMax(sesInds(reqInds)),'filled');
colormap('jet');caxis([0,10]);title('HLoc');

figure,
subplot(131);
scatter(FSrC(sesInds(reqInds),2),FSrC(sesInds(reqInds),3),20,rmapsMax(sesInds(reqInds)),'filled');
colormap('jet');caxis([0,10]);xlabel('Rear');ylabel('HLoc');title('max rate');
subplot(132);
scatter(FSrC(sesInds(reqInds),2),FSrC(sesInds(reqInds),3),20,fsrcz(sesInds(reqInds),1:3)+0.25 ,'filled');
colormap('jet');xlabel('Rear');ylabel('HLoc');title('sig');
subplot(133);
scatter(FSrC(sesInds(reqInds),2),FSrC(sesInds(reqInds),3),20,bsi(sesInds(reqInds)),'filled');
colormap('jet');caxis([0,1.2]),xlabel('Rear');ylabel('HLoc');title('bsi');


figure,
subplot(131);
scatter(FSrC(sesInds(reqInds),1),FSrC(sesInds(reqInds),3),20,rmapsMax(sesInds(reqInds)),'filled');
colormap('jet');caxis([0,10]);xlabel('LLoc');ylabel('HLoc');title('max rate');
subplot(132);
scatter(FSrC(sesInds(reqInds),1),FSrC(sesInds(reqInds),3),20,fsrcz(sesInds(reqInds),1:3)+0.25 ,'filled');
colormap('jet');xlabel('LLoc');ylabel('HLoc');title('sig');
subplot(133);
scatter(FSrC(sesInds(reqInds),1),FSrC(sesInds(reqInds),3),20,bsi(sesInds(reqInds)),'filled');
colormap('jet');caxis([0,1.2]),xlabel('LLoc');ylabel('HLoc');title('bsi');

figure();plot3(fsrcz(sesInds(reqInds)&fsrcz(:,2)<1,1),...
                  fsrcz(sesInds(reqInds)&fsrcz(:,2)<1,3),...
               bsi(sesInds(reqInds)&fsrcz(:,2)<1),'.b')

figure();plot3(FSrC(sesInds(reqInds)&FSrC(:,2)<0,1),...
               FSrC(sesInds(reqInds)&FSrC(:,2)<0,3),...
               bsi(sesInds(reqInds)&FSrC(:,2)<0),'.b')


%FSrC
%fsrcz
figure,
for f = 1:4
subplot2(4,2,f,1);imagesc(nunity(rmaps(:,indUnitFsrcz{f}&eigVars(reqInds,1)>40&eigVars(reqInds,2)<20))');
subplot2(4,2,f,2);imagesc(nunity(rmaps(:,indUnitFsrcz{f}&eigVars(reqInds,1)<50&eigVars(reqInds,2)>16))');
end


%fsrcz
figure,
for f = 1:4
subplot2(4,2,f,1);imagesc((rmaps(:,indUnitFsrcz{f}&eigVars(reqInds,1)>40&eigVars(reqInds,2)<20))');
subplot2(4,2,f,2);imagesc((rmaps(:,indUnitFsrcz{f}&eigVars(reqInds,1)<50&eigVars(reqInds,2)>16))');
end


figure();
hold('on');
subplot(141);scatter(eigVars(indR&sesInds&reqInds,1),...
                     eigVars(indR&sesInds&reqInds,2),...
                     25,...
                     FSrC(indUnitFsrcz{1},2),...
                     'filled');
colormap('jet');
caxis([-1,3]);
xlim([10,80]);ylim([0,45]);
subplot(142);scatter(eigVars(indH&sesInds&reqInds,1),...
                     eigVars(indH&sesInds&reqInds,2),...
                     25,...
                     FSrC(indUnitFsrcz{2},3),...
                     'filled');
colormap('jet');
caxis([-1,3]);
xlim([10,80]);ylim([0,45]);
subplot(143);scatter(eigVars(indL&sesInds&reqInds,1),...
                     eigVars(indL&sesInds&reqInds,2),...
                     25,...
                     FSrC(indUnitFsrcz{3},1),...
                     'filled');
colormap('jet');
caxis([-1,3]);
xlim([10,80]);ylim([0,45]);
subplot(144);
scatter(eigVars(indN&sesInds&reqInds,1),...
        eigVars(indN&sesInds&reqInds,2),...
        25,...
        FSrC(indUnitFsrcz{4},3),...
        'filled');
colormap('jet');
caxis([-1,3]);
xlim([10,80]);ylim([0,45]);

xlim([10,80]);ylim([0,45]);
subplot(142);scatter(eigVars(indH&sesInds&reqInds,1),eigVars(indH&sesInds&reqInds,2),5);
xlim([10,80]);ylim([0,45]);
subplot(143);scatter(eigVars(indL&sesInds&reqInds,1),eigVars(indL&sesInds&reqInds,2),5);
xlim([10,80]);ylim([0,45]);
subplot(144);scatter(eigVars(indN&sesInds&reqInds,1),eigVars(indN&sesInds&reqInds,2),5);
xlim([10,80]);ylim([0,45]);







figure();
hold('on');
subplot(141);
scatter(eigVars(indR&sesInds&reqInds,1)./(eigVars(indR&sesInds&reqInds,2)+eigVars(indR&sesInds&reqInds,1)),...
        bsi(indR(reqInds)&sesInds(reqInds))',...
        25,...
        fsrcz(indUnitFsrcz{1},2),...
        'filled');
colormap('jet');caxis([2,6])
subplot(142);
scatter(eigVars(indH&sesInds&reqInds,1)./(eigVars(indH&sesInds&reqInds,2)+eigVars(indH&sesInds&reqInds,1)),...
        bsi(indH(reqInds)&sesInds(reqInds))',...
        25,...
        fsrcz(indUnitFsrcz{2},3),...
        'filled');
colormap('jet');caxis([2,6])
subplot(143);
scatter(eigVars(indL&sesInds&reqInds,1)./(eigVars(indL&sesInds&reqInds,2)+eigVars(indL&sesInds&reqInds,1)),...
        bsi(indL(reqInds)&sesInds(reqInds))',...
        25,...
        fsrcz(indUnitFsrcz{3},1),...
        'filled');
colormap('jet');caxis([2,6])
subplot(144);
scatter(eigVars(indN&sesInds&reqInds,1)./(eigVars(indN&sesInds&reqInds,2)+eigVars(indN&sesInds&reqInds,1)),...
        bsi(indN(reqInds)&sesInds(reqInds))',...
        25,...
        fsrcz(indUnitFsrcz{4},[2,3,1])+0.25,...
        'filled');
%        FSrC(indUnitFsrcz{4},1:3)+1,...

figure
bar(bsxfun(@rdivide,eigenVectorVarianceRatios,sum(eigenVectorVarianceRatios,2)),'stacked')



%csmUnitSubset = cluSessionMap(unitSubsets{1},:);
% $$$ ind = indL;
% $$$ csmEigScrs = permute(eigScrsNnmf(sesInds&reqInds&ind,:,:),[1,3,2]);
% $$$ csmEigVecs = eigVecsNnmf(sesInds&reqInds&ind,:,:,:);
% $$$ csmEigVars = eigVars(sesInds&reqInds&ind,:);
% $$$ 
% $$$ csmEigScrs = eigScrs(sesInds&reqInds&ind,:,:);
% $$$ csmEigVecs = eigVecs(sesInds&reqInds&ind,:,:,:);
% $$$ csmEigVars = eigVars(sesInds&reqInds&ind,:);

csmEigScrs = eigScrsNnmf;
csmEigVecs = eigVecs;
csmEigVars = eigVars;

% $$$ csmEigScrs = permute(eigScrsNnmf,[1,3,2]);;
% $$$ csmEigVecs = eigVecsNnmf;%csmEigVars = eigVars;

clear('i');
[~,csmPhzPrefMax] = max(sq(csmEigScrs(:,:,:)),[],3);



% COMPUTE complex scores
csmEigScrComplex = sum(reshape(bsxfun(@times,reshape(bsxfun(@plus,csmEigScrs,abs(min(csmEigScrs,[],3))),[],size(csmEigScrs,3)),exp(-i.*pfs.adata.bins{2})'),size(csmEigScrs)),3)./sum(bsxfun(@plus,csmEigScrs,abs(min(csmEigScrs,[],3))),3);

csmEigScrComplex = sq(sum(csmEigScrs.*repmat(exp(-i.*pfs.adata.bins{2})',[size(csmEigScrs,1),1,size(csmEigScrs,3)]),2)./sum(csmEigScrs,2));


% COMPUTE each eigenvector's preferred theta phase
csmPhzPrefMean = angle(csmEigScrComplex);
csmPhzPrefMean(csmPhzPrefMean<0)  = csmPhzPrefMean(csmPhzPrefMean<0)+2*pi;
csmPhzPrefRes = abs(csmEigScrComplex);


% $$$ rose(csmPhzPrefMean( eigVars(indR&sesInds&reqInds,1)>thresholdVariance ...
% $$$                     &eigVars(indR&sesInds&reqInds,2)<thresholdVariance,1))

puInd =sesInds&reqInds&eigVars(:,1)>thresholdVariance&eigVars(:,2)<thresholdVariance; 
suInd =sesInds&reqInds&eigVars(:,1)>thresholdVariance&eigVars(:,2)>thresholdVariance; 
figure,plot(csmPhzPrefMean(suInd,1),csmPhzPrefMean(suInd,2),'.');

figure,
subplot(141);  rose(csmPhzPrefMean( indR&puInd,1));  title('Rear');
subplot(142);  rose(csmPhzPrefMean( indH&puInd,1));  title('High');
subplot(143);  rose(csmPhzPrefMean( indL&puInd,1));  title('Low');
subplot(144);  rose(csmPhzPrefMean( indN&puInd,1));  title('No Pref');

figure,
clf()
suptitle('Theta Phase Preferrence for Units with Single Factor');
sp = gobjects([0,1]);
sp(end+1)=subplot(141);  plot(csmPhzPrefMean( indR&puInd,1),csmPhzPrefRes( indR&puInd,1),'.');  title('Rear');
ylabel('Resultant Length');xlabel('Mean Phase');
sp(end+1)=subplot(142);  plot(csmPhzPrefMean( indH&puInd,1),csmPhzPrefRes( indH&puInd,1),'.');  title('High');
sp(end+1)=subplot(143);  plot(csmPhzPrefMean( indL&puInd,1),csmPhzPrefRes( indL&puInd,1),'.');  title('Low');
sp(end+1)=subplot(144);  plot(csmPhzPrefMean( indN&puInd,1),csmPhzPrefRes( indN&puInd,1),'.');  title('No Pref');
ForAllSubplots('xlim([0,2*pi]);ylim([0,0.75]);');
af(@(a) set(a,'XTick',[0,pi/2,pi,3*pi/2,2*pi]),sp);
%af(@(a) set(a,'XTickLabel',{'0','\pi/2','\pi','3*\pi/2','2*\pi'}),sp);
af(@(a) set(a,'XTickLabel',{'0','90','180','270','360'}),sp);
af(@(a) set(a,'Units','centimeters'),sp);
af(@(a) set(a,'Position',[a.Position(1:2),2.5,2.5]),sp);
ForAllSubplots('grid(''on'')');



figure,
clf()
suptitle('Theta Phase Preferrence for Units with Two Factors');
sp = gobjects([0,1]);
sp(end+1)=subplot(241);  plot(csmPhzPrefMean( indR&suInd,1),csmPhzPrefMean( indR&suInd,2),'.');  title('Rear');
ylabel({'Mean Phase','(2^{nd} factor)'});xlabel({'Mean Phase','(1^{st} factor)'});
sp(end+1)=subplot(242);  plot(csmPhzPrefMean( indH&suInd,1),csmPhzPrefMean( indH&suInd,2),'.');  title('High');
sp(end+1)=subplot(243);  plot(csmPhzPrefMean( indL&suInd,1),csmPhzPrefMean( indL&suInd,2),'.');  title('Low');
sp(end+1)=subplot(244);  plot(csmPhzPrefMean( indN&suInd,1),csmPhzPrefMean( indN&suInd,2),'.');  title('No Pref');
ForAllSubplots('xlim([0,2*pi]);ylim([0,2*pi]);');
af(@(a) set(a,'XTick',[0,pi/2,pi,3*pi/2,2*pi]),sp);
%af(@(a) set(a,'XTickLabel',{'0','\pi/2','\pi','3*\pi/2','2*\pi'}),sp);
af(@(a) set(a,'XTickLabel',{'0','90','180','270','360'}),sp);
af(@(a) set(a,'YTick',[0,pi/2,pi,3*pi/2,2*pi]),sp);
%af(@(a) set(a,'YTickLabel',{'0','\pi/2','\pi','3*\pi/2','2*\pi'}),sp);
af(@(a) set(a,'YTickLabel',{'0','90','180','270','360'}),sp);



sp(end+1)=subplot(245);  plot(csmPhzPrefRes( indR&suInd,1),csmPhzPrefRes( indR&suInd,2),'.');  title('Rear');
ylabel({'Resultant Length','(2^{nd} factor)'});xlabel({'Resultant Length', '(1^{st} factor)'});
sp(end+1)=subplot(246);  plot(csmPhzPrefRes( indH&suInd,1),csmPhzPrefRes( indH&suInd,2),'.');  title('High');
sp(end+1)=subplot(247);  plot(csmPhzPrefRes( indL&suInd,1),csmPhzPrefRes( indL&suInd,2),'.');  title('Low');
sp(end+1)=subplot(248);  plot(csmPhzPrefRes( indN&suInd,1),csmPhzPrefRes( indN&suInd,2),'.');  title('No Pref');
af(@(a) xlim(a,[0,0.75]), sp(5:end))
af(@(a) set(a,'XTick',[0,0.25,0.5,0.75]), sp(5:end))
af(@(a) ylim(a,[0,0.75]), sp(5:end))
af(@(a) set(a,'YTick',[0,0.25,0.5,0.75]), sp(5:end))
af(@(a) set(a,'Units','centimeters'),sp);
af(@(a) set(a,'Position',[a.Position(1:2),2.5,2.5]),sp);
ForAllSubplots('grid(''on'')');


ind = sesInds&reqInds&indH;
figure();
hist(csmPhzPrefMean(ind & eigVars(:,1)>thresholdVariance ...
                        & eigVars(:,2)<thresholdVariance,1),20)

figure,hist(csmPhzPrefMean(eigVars(indL&sesInds&reqInds,2)>thresholdVariance,1),20)

figure,hist(circ_dist(csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)),30)


figure,plot(pfs.adata.bins{2}(abs(csmPhzPrefMax(:,1)))+randn(size(csmPhzPrefMax(:,1)))/10,...
            csmPhzPrefMean(:,1),'.');
% $$$ figure,plot(sq(csmEigScrs(:,1,:))')

figure,plot(log10(csmEigVars(:,1)),log10(csmEigVars(:,2)),'.')

figure();
hold('on');
plot(circ_dist(csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)),csmEigVars(:,1),'.b')
plot(circ_dist(csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)),csmEigVars(:,2),'.r')


figure,
plot(circ_mean([csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)]')+pi,csmEigVars(:,1),'.b')

for j = 1:758,
    for k = 1:5,
        eigVarsNnmf(j,k) = sum(reshape(reshape(eigVecsNnmf(j,k,:,:),[],1)*eigScrsNnmf(j,:,k),[],1),'omitnan');
    end
end


figure,
for j = 600:758,
    clf();
    for k = 1:5,
        subplot(5,5,k);
        imagescnan(sq(eigVecs(j,k,:,:))');
        axis('xy');
        subplot(5,5,k+5);
        plot(sq(eigScrs(j,k,:)));        
        subplot(5,5,k+10);
        imagescnan(sq(eigVecsNnmf(j,k,:,:))');
        axis('xy');
        subplot(5,5,k+15);
        plot(eigScrsNnmf(j,:,k));        
    end
    subplot(5,5,21);    
    plot(pfds{cluSessionMap(j,1)},cluSessionMap(j,2),1,'text',[],false);
    subplot(5,5,22);        
    plot(eigVars(j,:),'-+');
    title(num2str(cluSessionMap(j,:)));
    subplot(5,5,23);
    plot(eigVarsNnmf(j,:),'-+');
    title(num2str(cluSessionMap(j,:)));
    waitforbuttonpress();
end

    
    

% $$$ 
% $$$ 
% $$$ p = pfs.adata.bins{2};
% $$$ for h = 1:size(rmap,1),
% $$$     for b = 1:size(rmap,3),
% $$$         rc = rmap(h,:,b)'.*exp(i*p);
% $$$         Rmat = rc*conj(rc)';
% $$$         R(h,b) = sum(sum(triu(Rmat,1)))./sum(sum(triu(rmap(h,:,b)'*rmap(h,:,b),1)));
% $$$     end
% $$$ end
% $$$ 
% $$$ figure();
% $$$ subplot(121);
% $$$ imagesc(abs(R)'),axis('xy');
% $$$ subplot(122);
% $$$ imagescnan({pfd{1}.adata.bins{:},abs(rmapComplex)'},...
% $$$            [],[],false,nanColor,'colorMap',@parula);
% $$$ axis('xy')
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure(), 
% $$$ pBins = linspace(0,2*pi,17);
% $$$ ks = [0.1:0.2:2];
% $$$ sp = tight_subplot(2,5,0.05,0.05);
% $$$ for k = 1:numel(ks),
% $$$     axes(sp(k));
% $$$     hold('on');
% $$$     for b = [10,100,1000],
% $$$         for c = 1:1000,
% $$$             th = VonMisesRnd(0,ks(k),b,1);
% $$$             pBinsInds = discretize(th,pBins);    
% $$$             thCount = accumarray(pBinsInds,ones(size(th)),[numel(pBins)-1,1],@sum);
% $$$             rc = thCount.*exp(i*p);
% $$$             Rmat = rc*conj(rc)';
% $$$             R = sum(sum(triu(Rmat,1)))./sum(sum(triu(thCount*thCount',1)));
% $$$             Rasm(c) = abs(R);
% $$$             Rppc(c) = PPC(th);
% $$$         end
% $$$         plot(Rppc,Rasm,'.');
% $$$     end
% $$$ end
% $$$ 
% $$$ af(@(s) xlim(s,[-0.2,1]),sp);
% $$$ af(@(s) ylim(s,[-0.2,1]),sp);
% $$$ af(@(s) grid(s,'on'),sp);




% req20190527();
sigma = 150;
pftHZTPD = cf(@(s) cf(@(T,u) MTAApfs(T,u,'tag',['ddtp-','s',num2str(sigma),'-',s]), Trials, units), stateLabels);
pftTPZDa = cf(@(pfs) ...
              cf(@(p,u) ...
                 p.data.rateMap(:,ismember(p.data.clu,u),:),  ...
                 pfs,units), ...
              pftHZTPD);
% $$$ pftTPZDa = cf(@(pfs) ...
% $$$               cf(@(p,u) ...
% $$$                  mean(reshape(permute(RectFilter(reshape(permute(p.data.rateMap(:,ismember(p.data.clu,u),:),[1,4,2,3]),[p.adata.binSizes',numel(u),p.parameters.numIter]),3,1,'circular'),[2,1,3,4]),size(p.data.rateMap,1),numel(u),p.parameters.numIter),3,'omitnan'),  ...
% $$$                  pfs,units), ...
% $$$               pftHZTPD);


u =3;figure;plot(sq(max(reshape(pftTPZDa{3}{20}(:,u,:),[40,16,101]))),'c');hold on, ...
   plot(sq(max(reshape(pftTPZDa{4}{20}(:,u,:),[40,16,101]))),'b');

u =3;figure;plot(mean(RectFilter(sq(max(reshape(pftTPZDa{3}{20}(:,u,:),[40,16,101]))),3,1,'circular'),2),'c');
hold('on'); plot(mean(RectFilter(sq(max(reshape(pftTPZDa{4}{20}(:,u,:),[40,16,101]))),3,1,'circular'),2),'b');

clu     = cf(@(pfs) cf(@(p,u) p.data.clu(:,ismember(p.data.clu,u),:),                        pfs,units), pftHZTPD);
tlu     = cf(@(i,u) repmat(i,size(u)), mat2cell(1:numel(units),1,ones([1,numel(units)])),units);
pftTPZDa = cf(@(pfs) cat(2,pfs{:}), pftTPZDa);
clu     = cf(@(c)   cat(2,c{:}),clu);
assert(all(cell2mat(cf(@(c) all(c), cf(@(c,o) c==o, clu,circshift(clu,1,2))))),'Clus do not match');
pftTPZDa = cat(4,pftTPZDa{:});
clu = clu{1};
tlu=cat(2,tlu{:});
clu=[tlu',clu'];
[clu,rind] = sortrows(clu);
clu     = clu(unitSubset,:);
pftTPZDa = pftTPZDa(:,rind,:,:);
pftTPZDa = pftTPZDa(:,unitSubset,:,:);
pftTPZDa = reshape(permute(pftTPZDa,[1,5,2,3,4]),[pftHZTPD{1}{1}.adata.binSizes',numel(unitSubset),pftHZTPD{1}{1}.parameters.numIter,numel(stateLabels)]);
%pftTPZDa = reshape(permute(pftTPZDa,[1,5,2,3,4]),[pftHZTPD{1}{1}.adata.binSizes',numel(unitSubset),1,numel(stateLabels)]);

% For each unit, for each state, for each phase, find the maximum rate and position along the gdz

[grate,gpos] = max(permute(sq(mean(RectFilter(permute(pftTPZDa(:,:,sesInds(unitSubset),:,:),[2,1,3,4,5]),3,1,'circular'),4,'omitnan')),[2,1,3,4,5]));
% $$$ [grate,gpos] = max(permute(sq(mean(permute(pftTPZDa(:,:,sesInds(unitSubset),:,:),[2,1,3,4,5]),4,'omitnan')),[2,1,3,4,5]));
grate = sq(grate(:,:,:,2:4));
gpos = sq(gpos);
[prate,pphz] = max(grate);
clusub = clu(sesInds(unitSubset),:);


npftTPZa = reshape(grate,16,[]);
npftTPZa(~nniz(npftTPZa(:))) = 0;
nind = sum(npftTPZa==0)==16;
npftTPZa(:,nind) =[];

% $$$ [LU,LR,FSr,VT] = erpPCA(npftTPZa',3);
% $$$ fscr = zeros([numel(nind),3]);
% $$$ fscr(~nind,:) = FSr;
% $$$ fscr = reshape(fscr,[],3,3);
% $$$ figure,plot(LR(phzOrder,:))
% $$$ figure,plot(VT(1:10,4),'-+');


[W,H,D] = nnmf(npftTPZa,3);
fscr = zeros([numel(nind),3]);
fscr(~nind,:) = H';
fscr = reshape(fscr,[],3,3);

figure,plot(W(phzOrder,:))


msr = sq(max(sq(grate(phzOrder,:,:))));
[msr,msi] = sort(msr,2,'descend');
dsr = sq(mean(sq(grate(phzOrder,:,:))));
[dsr,dsi] = sort(dsr,2,'descend');






% REVIEW examples 
figure(293023);
for u = 1:144;
clf();
% F scores
subplot(411);
imagesc(sq(fscr(u,:,:)));
% phase rates
subplot(412);
plot(sq(grate(phzOrder,u,1:3)));
legend(stateLabels(2:4));
% Eigenvectors
subplot(413);
plot(W(phzOrder,:));
% Delta F scores (dominate - subordinate)
subplot(414);
%sq(sum(fscr(u,:,:).*double(fscr(u,:,:)>0),3),'k-+')
%plot(bsxfun(@minus,sq(fscr(u,dsi(u),:)),sq(fscr(u,~ismember([1:3],dsi(u)),:))'));
% GET index of domiate state
waitforbuttonpress();
end


% DOMINATE state 
stateLabelSubset = stateLabels(2:4);
figure,
for s = 1:3,
ind = dsi==s;
dfscr = bsxfun(@minus,fscr(ind,~ismember(1:size(fscr,2),s),:),fscr(ind,s,:));
dgrate = bsxfun(@minus,grate(:,ind,[false,~ismember(2:4,s+1),false,false,false,false]),grate(:,ind,s+1));
rdsr = dsr(ind,s)./mean(dsr(ind,~ismember(1:size(fscr,2),s)),2);
subplot(3,2,2*s-1);
hold('on');
grid('on');
%plot3(dfscr(:,1,1),dfscr(:,1,2),dfscr(:,1,3),'.');
plot(sq(dfscr(:,1,:))');
%plot(sq(dgrate(phzOrder,:,1)));
%plot(mean(sq(dgrate(phzOrder,:,1)),2),'k--','LineWidth',3);
title(stateLabelSubset(~ismember(1:size(fscr,2),s)))
subplot(3,2,2*s);
hold('on');
grid('on');
%plot3(dfscr(:,2,1),dfscr(:,2,2),dfscr(:,2,3),'.');
plot(sq(dfscr(:,2,:))');
%plot(sq(dgrate(phzOrder,:,2)));
%plot(mean(sq(dgrate(phzOrder,:,2)),2),'k--','LineWidth',3);
end


% DOMINATE state rate difference
stateLabelSubset = stateLabels(2:4);
figure,
for s = 1:3,
ind = dsi==s;
%dfscr = bsxfun(@minus,fscr(ind,~ismember(1:size(fscr,2),s),:),fscr(ind,s,:));
dgrate = bsxfun(@minus,grate(:,ind,[false,~ismember(2:4,s+1),false,false,false,false]),grate(:,ind,s+1));

subplot(3,2,2*s-1);
hold('on');
grid('on');
%plot3(dfscr(:,1,1),dfscr(:,1,2),dfscr(:,1,3),'.');
%plot(sq(dfscr(:,1,:))');
plot(sq(dgrate(phzOrder,:,1)));
plot(mean(sq(dgrate(phzOrder,:,1)),2),'k--','LineWidth',3);
title(stateLabelSubset(~ismember(1:size(fscr,2),s)))

subplot(3,2,2*s);
hold('on');
grid('on');
%plot3(dfscr(:,2,1),dfscr(:,2,2),dfscr(:,2,3),'.');
%plot(sq(dfscr(:,2,:))');
plot(sq(dgrate(phzOrder,:,2)));
plot(mean(sq(dgrate(phzOrder,:,2)),2),'k--','LineWidth',3);
end


% Compare erpPCA fscores for all cells within dominate state
stateLabelSubset = stateLabels(2:4);
figure,
for s = 1:3,
ind = dsi==s;
dfscr = bsxfun(@minus,fscr(ind,~ismember(1:size(fscr,2),s),:),fscr(ind,s,:));
subplot(3,3,3*s-2);
plot3(msr(ind,3),msr(ind,2),msr(ind,1),'.');
xlim([0,50]);ylim([0,50]);zlim([0,50]);daspect([1,1,1])
subplot(3,3,3*s-1);
hold('on');
grid('on');
plot3(dfscr(:,1,1),dfscr(:,1,2),dfscr(:,1,3),'.');
title(stateLabelSubset(~ismember(1:size(fscr,2),s)))
subplot(3,3,3*s);
hold('on');
grid('on');
plot3(dfscr(:,2,1),dfscr(:,2,2),dfscr(:,2,3),'.');
end
ForAllSubplots('ylim([-0.2,0.2]);xlim([-0.2,0.2])')
%

% CONSOLIDATE the fscore differences 
%  color by dominate state
%  order fscore differences by mean rate
% 
figure();
stsColor = 'rgb';
sp = gobjects([0,1]);
sp(end+1) = subplot(121);
sp(end+1) = subplot(122);
af(@(h) hold(h,'on'), sp);
for s = 1:3
    ind = dsi(:,1)==s;
    dfscr = bsxfun(@minus,fscr(ind,~ismember(1:size(fscr,2),s),:),fscr(ind,s,:));
    dispOpts = {'.',...
                'MarkerFaceColor',stsColor(s),...
                'MarkerEdgeColor',stsColor(s),...
                'MarkerSize',10};
    for j = 1:size(dfscr,1),
% $$$         if dsi(j,2)>dsi(j,3)
% $$$             %if dsr(j,2)>4,
% $$$                 axes(sp(1));
% $$$                 plot(dfscr(j,2,1),dfscr(j,2,2),dispOpts{:});
% $$$                 %end
% $$$                 %if dsr(j,3)>4,
% $$$                 %axes(sp(2));
% $$$                 plot(dfscr(j,1,1),dfscr(j,1,2),dispOpts{:});
% $$$                 %end
% $$$         else
% $$$             %if dsr(j,2)>4,
% $$$                 axes(sp(1));
% $$$                 plot(dfscr(j,1,1),dfscr(j,1,2),dispOpts{:});
% $$$                 %end
% $$$                 %if dsr(j,3)>4,
% $$$                 %axes(sp(2));
% $$$                 plot(dfscr(j,2,1),dfscr(j,2,2),dispOpts{:});
% $$$                 %end
% $$$         end
        if dsi(j,2)>dsi(j,3)
            %if dsr(j,2)>4,
                axes(sp(1));
                plot3(dfscr(j,2,1),dfscr(j,2,2),dfscr(j,2,3),dispOpts{:});
                %end
                %if dsr(j,3)>4,
                %axes(sp(2));
                plot3(dfscr(j,1,1),dfscr(j,1,2),dfscr(j,1,3),dispOpts{:});
                %end
        else
            %if dsr(j,2)>4,
                axes(sp(1));
                plot3(dfscr(j,1,1),dfscr(j,1,2),dfscr(j,1,3),dispOpts{:});
                %end
                %if dsr(j,3)>4,
                %axes(sp(2));
                plot3(dfscr(j,2,1),dfscr(j,2,2),dfscr(j,2,3),dispOpts{:});
                %end
        end
        
    end
end
ForAllSubplots('ylim([-0.2,0.1]);xlim([-0.2,0.1])')


