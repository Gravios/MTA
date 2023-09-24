
%%%<<< SET default arguments -----------------------------------------------------------------------
%MjgER2016_figure_BhvPlacefields_args();
configure_default_args();
% MjgER2016_figure_BhvPlacefields_args:
%
% Default arument overrides:
%     fet_HB_pitchB
%     compute_bhv_ratemaps
%     compute_bhv_ratemaps_shuffled
%     compute_bhv_ratemap_erpPCA
%%%>>>

%%%<<< LOAD MjgER2016 data -------------------------------------------------------------------------
MjgER2016_load_data();
%  MjgER2016_load_data:
%
%  Variables:
%      Trials      units      cluSessionMap      pitchReferenceTrial      FigDir
%      sessionListName        sessionList        states                   numStates      
%      interpParPfsp          interpParDfs
%      
%  Functions:
%      reshape_eigen_vector
%%%>>>

%%%<<< LOAD data -----------------------------------------------------------------------------------
% EXAMPLE Trial : jg05-20120312
tind = 20;
Trial = Trials{tind};
stc = Trials{tind}.stc.copy();

% LOAD rate maps -----------------------------------------------------------------------------------
% LOAD place restricted behavior fields
bfs      = cf(@(t,u)  compute_bhv_ratemaps(t,u),          Trials, units);
bfsEx    = bfs{tind};
bfsShuff = cf(@(t,u)  compute_bhv_ratemaps_shuffled(t,u), Trials, units);
% LOAD all place fields
pfst = cf(@(t,u)  pfs_2d_theta(t,u),    Trials,units);
pfss = cf(@(t,u)  pfs_2d_states(t,u),   Trials,units);
pfsr = cf(@(t,u)  pfs_2d_states(t,u,'',{'rear&theta','hbhv&theta','lbhv&theta'}),   Trials,units);
pfsa = cf(@(s,t)  cat(2,{t},s),       pfss,pfst);
% --------------------------------------------------------------------------------------------------


% COMPUTE bfs erpPCA -------------------------------------------------------------------------------
[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfs, units, [], [], false);
numComp = size(eigVecs,2);
fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan(size(validDims));
    fpc{i}(validDims) = eigVecs(:,i);
end
fpcLims = [min(cellfun(@min,fpc)),max(cellfun(@max,fpc))];
% RESHAPE eigenvectors
bhvMask = false(size(validDims));
bhvMask(validDims) = true;
bhvMask = reshape_eigen_vector(bhvMask,bfs)';
% SET behavior space visual limits
bhvLims = [-1.7, 0.5; ...
           -0.5, 1.7];
% --------------------------------------------------------------------------------------------------


% LOAD example theta place fields
pft = pfs_2d_theta(Trial);
pft.parameters.states = 'theta-groom-sit';
% LOAD behavior place field
pfs = cat(2,{pft},pfs_2d_states(Trial));
% SORT place field states to match states
pfStates = cf(@(s) ['^',s,'$'],cf(@(p) p.parameters.states,pfs));
for s = 1:numel(pfStates),
    pfStates{s} = strrep(pfStates{s},'&','[&]');
    %pfStates{s} = strrep(pfStates{s},'-','[-]');
end    
for s = 1:numStates,
    psi(s) = find(~cellfun(@isempty,regexp(repmat(states(s),size(pfStates)),pfStates)));
end
pfs = pfs(psi);







[hfig,fig,fax,sax] = set_figure_layout(figure(666087),'A4','portrait',[],1.75,1.75,0.05,0.05);
% $$$ % Multisession examples
% $$$ % selected units session x unit id
% $$$ clumap = [17,147;...
% $$$            3,171;...
% $$$            3,197;...
% $$$           16, 51;...
% $$$            1, 15;...
% $$$           17,165;...
% $$$           18,129;...
% $$$          ];

% two patches [63, 89, 133


cluMap = [20,44;... 74;...
          20,31;... %20,73;...
          20,34;... double
          20,81;... %20,85;...
          20,83;...
          20,141;...
          20,79;...
          20,103];
%20,104];
% ...          20,109;...
%          20,59];

% {accg,BFS{HP,BP},BFS{HP,BS},theta,rear,hloc,lloc,hpause,lpause}


%labels = {{'Behavior','Rate Map'},{'Theta','Rate Map'},'Rear','H Loc','H Pause','L Loc','L Pause'};
labels = {{'Theta','    '},'Rear','H Loc','H Pause','L Loc','L Pause',{'Posture'}};

numStates = numel(states);
% PLOT 
xind = 1;
yind = 1;
yinit = yind;
for u = cluMap',
    uid = find(ismember(cluSessionMap,u','rows'));
    xind = 7;
% SET color scale max for each row -----------------------------------------------------------------
    maxPfsRate = max([cell2mat(cf(@(p,u,m)                                                        ...
                                 p.maxRate(u,true&numel(m)==1,1,[],[],[],m),                     ...
                                 [pfs,{bfsEx}],                                                  ...
                                 repmat({u(2)},[1,numel(pfs)+1]),                                ...
                                 [repmat({true},[1,numel(pfs)]),{bhvMask(:)}])                   ...
                              ),max(nonzeros(rmapP(validDims,uid,:)))]);
    pfsMaxRatesMean = cell2mat(cf(@(p,u) ...
                                  max(p.maxRate(u,true,'mean')),...
                                  pfs,repmat({u(2)},[1,numel(pfs)])));
    bfsMaxRatesMean = max(bfsEx.maxRate(u(2),false,1,[],[],[],bhvMask(:)));


% BHV FIELDS ---------------------------------------------------------------------------------------

    patchInds = find(~isnan(patchCntrF(uid,:,1)));
    for patchId = patchInds
        if xind > 8
            break;
        end
        sax(end+1) = axes('Units','centimeters',                                                 ...
                          'Position',[fig.page.xpos(xind)+0.2*(numStates-3),                     ...
                                      fig.page.ypos(yind),                                       ...
                                      fig.subplot.width,                                         ...
                                      fig.subplot.height],                                       ...
                          'FontSize', 8);
        %bfsEx.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet,bhvMask,nanColor);
        rmapEx = rmapP(:,ismember(cluSessionMap,u','rows'),patchId);
        rmapEx(~validDims) = nan;
        imagescnan({bfsEx.adata.bins{:},reshape(rmapEx,bfsEx.adata.binSizes')'},                 ...
                   [0,maxPfsRate],                                                               ...
                   'linear',                                                                     ...
                   false,                                                                        ...
                   nanColor,                                                                     ...
                   [],[],@jet);
        axis('xy');
        sax(end).YTickLabel = {};
        sax(end).XTickLabel = {};        
        % ANNOTATE plot with max rate
        text(0.45,1.475,num2str(round(max(rmapEx))),'FontSize',8,'Color',[1,1,1],'HorizontalAlignment','right');
        xlim(bhvLims(1,:));
        ylim(bhvLims(2,:));
        sax(end).Color = nanColor;
        
        patchColors = 'mw';        
        rectangle(fax,'Position',sax(end).Position+[0.02,0.02,-0.04,-0.04],'EdgeColor',patchColors(patchId));
        uistack(fax,'top');                  
        xind = xind + 1;        

        if yind == yinit && patchId == 1
            line([-1.6,-1.6],...
                 [1.1,1.6],...
                 'LineWidth',2,...
                 'Color','w');
            line([-1.6,-1.1],...
                 [1.6,1.6],...
                 'LineWidth',2,...
                 'Color','w');
        end
    end%for patchId 

% $$$     if yind == yinit, 
% $$$         title(labels{7});
% $$$     end 
    if yind==round(size(cluMap,1)),
        ylabel({'Body Pitch'});
        sax(end).YAxisLocation = 'right';
    end

% DRAW the scale bars for placefields (spatial ratemaps)
    if yind==size(cluMap,1),
        xlabel({'Head-Body','Pitch'});
% $$$         % y-axis scale bar
% $$$         line(fax,...
% $$$              [sax(end).Position([1])].*[1,1]-0.15,...     
% $$$              [sax(end).Position([2]),sum(sax(end).Position([2,4]).*[1,0.5])],...     
% $$$              'LineWidth',1,...
% $$$              'Color',[0,0,0]);
% $$$         % y-axis scale bar label        
% $$$         text(fax,...
% $$$              sax(end).Position(1)-0.35,...
% $$$              sum(sax(end).Position([2,4]).*[1,0.25]),...
% $$$              '1 rad',...
% $$$              'HorizontalAlignment','center',...
% $$$              'VerticalAlignment','middle',...
% $$$              'Rotation',90,...
% $$$              'FontSize',8);
    end
    
    xind = 1;
    for s = 1:numStates
% PLACEFIELDS MTAApfs ------------------------------------------------------------------------------
        sax(end+1) = axes('Units','centimeters',                                                 ...
                         'Position',[fig.page.xpos(xind)+0.2+0.2.*double(s>=2),                  ...
                                     fig.page.ypos(yind),                                        ...
                                     fig.subplot.width,                                          ...
                                     fig.subplot.height],                                        ...
                         'FontSize', 8);
        plot(pfs{s},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpParPfs,@jet,[],nanColor);
        sax(end).YTickLabel = {};
        sax(end).XTickLabel = {};
        if yind == yinit, 
            tax = title(labels{xind});
            if s==1,
                %tax.Position(2) = 800;
            end        
        end
                
        xind = xind+1;
        text(490,380,num2str(round(pfsMaxRatesMean(s))),'FontSize',8,'Color',[1,1,1],'HorizontalAlignment','right');
        %text(200,380,num2str(cond_round(pfsMaxRatesMean(s))),'FontSize',8,'Color',[1,1,1]);
        if yind == yinit && s == 1
            line([-480,-480],...
                 [280,480],...
                 'LineWidth',2,...
                 'Color','w');
            line([-480,-280],...
                 [480,480],...
                 'LineWidth',2,...
                 'Color','w');
        end
        
        
        if xind ==2
            hold(sax(end),'on');
            % Circles
            uid = find(ismember(cluSessionMap,u','rows'));
            gpatches = find(~isnan(patchCntrF(uid,:,1)));
            patchColors = 'mw';
            if numel(gpatches)>2,
                gpatches = gpatches(1:2);
            end
            for gp = gpatches
                circle(patchCntrF(uid,gp,1),patchCntrF(uid,gp,2),100,['-',patchColors(gp)]);
            end
% CONTURES 
% $$$             rmap = plot(pfs{1},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpParPfs);
% $$$             [~,mxp] = pft.maxRate(u(2));
% $$$             binGrid = cell([1,2]);
% $$$             [binGrid{:}] = ndgrid(interpParPfs.bins{:});
% $$$             binGrid = cat(3,binGrid{:});
% $$$             binDist = sqrt(sum((binGrid-repmat(permute(mxp,[1,3,2]),[cellfun(@numel,interpParPfs.bins),1])).^2,3));
% $$$             [cont,conHax] = contour(binGrid(:,:,1),binGrid(:,:,2),(rmap./max(rmap(:)))>0.2&binDist<250,[0.5,0.5],'m','LineWidth',1);
% $$$ 
% $$$             % REMOVE first coordinates for this example
% $$$             cont(:,1) = [];
% $$$             
% $$$             [contYMax,contYMaxInd] = max(reshape(binGrid(:,:,2).*double((rmap./max(rmap(:)))>0.2&binDist<250),[],1));
% $$$             [contYMin,contYMinInd] = min(reshape(binGrid(:,:,2).*double((rmap./max(rmap(:)))>0.2&binDist<250),[],1));
% $$$             contMaxInd = zeros([1,2]);
% $$$             [contMaxInd(1),contMaxInd(2)] = ind2sub(cellfun(@numel,interpParPfs.bins),contYMaxInd);
% $$$             contMinInd = zeros([1,2]);
% $$$             [contMinInd(1),contMinInd(2)] = ind2sub(cellfun(@numel,interpParPfs.bins),contYMinInd);
% $$$             contMaxPos = sq(binGrid(contMaxInd(1),contMaxInd(2),:));
% $$$             contMinPos = sq(binGrid(contMinInd(1),contMinInd(2),:));
% $$$             
% $$$             [contMaxYord, contMaxInd] = max(cont(2,:));
% $$$             contMaxXord = cont(1,contMaxInd);
% $$$             [contMinYord, contMinInd] = min(cont(2,:));
% $$$             contMinXord = cont(1,contMinInd);
% $$$             
% $$$             [~,contMinXind] = NearestNeighbour(interpParPfs.bins{1},contMinXord);
% $$$             [~,contMinYind] = NearestNeighbour(interpParPfs.bins{2},contMinYord);            
% $$$             [~,contMaxXind] = NearestNeighbour(interpParPfs.bins{1},contMaxXord);
% $$$             [~,contMaxYind] = NearestNeighbour(interpParPfs.bins{2},contMaxYord);            

% $$$             if yind == 1,
% $$$              % TOP line                
% $$$              line(fax,...
% $$$                   [sax(end-1).Position(1)+fig.subplot.width, sax(end).Position(1)+(contMaxXind./numel(interpParPfs.bins{1})).*fig.subplot.width],...
% $$$                   [sax(end-1).Position(2)+fig.subplot.height,sax(end).Position(2)+(contMaxYind./numel(interpParPfs.bins{2})).*fig.subplot.height],...
% $$$                   'Color','m','LineWidth',1);
% $$$             % BOTTOM line            
% $$$              line(fax,...
% $$$                   [sax(end-1).Position(1)+fig.subplot.width, sax(end).Position(1)+(contMinXind./numel(interpParPfs.bins{1})).*fig.subplot.width],...
% $$$                   [sax(end-1).Position(2),sax(end).Position(2)+(contMinYind./numel(interpParPfs.bins{2})).*fig.subplot.height],...
% $$$                   'Color','m','LineWidth',1);
% $$$              rectangle(fax,'Position',sax(end-1).Position+[-0.02,-0.02,0.04,0.04],'EdgeColor','m');
% $$$             end% if yind==1
% $$$             uistack(fax,'top');
        end% if xind==2
% $$$         if  yind == 1 && s == numStates,
% $$$             cax = colorbar(sax(end));
% $$$             colormap(cax,'jet');
% $$$             cax.Units = 'centimeters';
% $$$             cax.Position = [sum(sax(end).Position([1,3]))+0.15,sax(end).Position(2),0.15,sax(end).Position(4)];
% $$$             cax.XTick = [0,1];
% $$$             cax.XTickLabel = {'0','Max'};
% $$$             cax.Label.String = 'Hz';
% $$$             cax.Label.Units = 'centimeters';
% $$$             cax.Label.FontSize = 8;
% $$$             cax.Label.Position(1) = 0.3;
% $$$         end
    end% for s
    yind = yind+1;
end


% DRAW line and arrows for theta to behavioral state decomposition
s = 7*size(cluMap,1)-3;
%sax(end-s:end-s+5)
% HORIZONTAL line
line(fax,                                                                                         ...
     [sax(3).Position(1)+fig.subplot.width-0.2,sax(8).Position(1)+fig.subplot.width/2+0.07],                  ...
      repmat(sax(3).Position(2)+sax(3).Position(4)+0.575,[1,2]),                            ...
     'Color','k',                                                                                 ...
     'LineWidth',1);
for a = 1:5,
    patch(fax,                                                                                    ...
          repmat(sax(3+a).Position(1)+fig.subplot.width/2,[1,3])+[-0.07,0.07,0],                        ...
          repmat(sax(3+a).Position(2)+sax(3+a).Position(4)+0.475,[1,3])+[0.1,0.1,-0.09],        ...
          'k');
end


