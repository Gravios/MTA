% MjgER2016_figure_singleUnitEPP
%
% Definitions:
%    EPP - Egocentric Phase Precession  
%        FEPP - Forward
%        LEPP - Lateral
%    RTM - Realative To Maze
%    HFV - Head Forward Velocity
%    HLV - Head Lateral Velocity
%    HPM - Head Pitch RTM
%    HRM - Roll RTM
%    HBA - Head Body Angle
%
% Relations:
%    FEPP ( HFV, HPM, HRM, HRM )
%    LEPP ( HBA, HFL, HVA )
%
% TODO : compute FWE threshold from Null distributions
% 

p2z = @(p) icdf('normal',p,0,1);
goodTrialInds = find(~cellfun(@isempty,units));

%%%<<< egohba

pfs = cf(@(t,u,x,s,p,rt,hc,ch,pc)                              ... Egocentric ratemap given theta phase and head body angle.
         compute_egohba_ratemap(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),   ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed
             num2cell(phzCorrection)                           ... theta phase offset
);

pfsh = cf(@(t,u,x,s,p,rt,hc,ch,pc)                             ... Egocentric ratemap given theta phase and head body angle.
         compute_egohba_ratemap_shuffled(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),   ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed
             num2cell(phzCorrection)                           ... theta phase offset
);
% egohba

%%%>>>

%%%<<< egohbahvl

pfl = cf(@(t,u,x,s,p,rt,hc,ch,pc)                              ... Egocentric ratemap given (TP, HBA, HVL gt 0)
         compute_egohbahvl_ratemap(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed 
             num2cell(phzCorrection)                           ... theta phase offset
);

pflh = cf(@(t,u,x,s,p,rt,hc,ch,pc)                             ... Egocentric ratemap given (TP, HBA, HVL gt 0)
         compute_egohbahvl_ratemap_shuffled(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed 
             num2cell(phzCorrection)                           ... theta phase offset
);

%%%>>>

%%%<<< egohvf

t = 1:28;

pfv = cf(@(t,u,x,s,p,rt,hc,ch,pc)                              ... Egocentric ratemap given theta phase and head body angle.
         compute_egohvf_ratemap(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),   ...
             Trials(t),                                           ... MTATrial
             units(t),                                            ... Unit subset, placefields away from the maze walls
             xyz(t),                                              ... MTADxyz object, head position
             spk(t),                                              ... MTASpk object, spike time and id collection 
             pft(t),                                              ... MTAApfs object, theta state placefields 
             num2cell(rot(t)),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection(t)),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan(t)),                              ... lfp channel from which phase is computed
             num2cell(phzCorrection(t))                           ... theta phase offset
);

pfvh = cf(@(t,u,x,s,p,rt,hc,ch,pc)                             ... Egocentric ratemap given theta phase and head body angle.
         compute_egohvf_ratemap_shuffled(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),   ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed
             num2cell(phzCorrection)                           ... theta phase offset
);

%%%>>>

%%%<<< egohrl

pfr = cf(@(t,u,x,s,p,rt,hc,ch,pc)                              ... Egocentric ratemap given theta phase and head body angle.
         compute_egohrl_ratemap(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),   ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hrlCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed
             num2cell(phzCorrection)                           ... theta phase offset
);

pfrh = cf(@(t,u,x,s,p,rt,hc,ch,pc)                             ... Egocentric ratemap given theta phase and head body angle.
         compute_egohrl_ratemap_shuffled(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),   ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hrlCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed
             num2cell(phzCorrection)                           ... theta phase offset
);
% egohrl
%%%>>>

t = 20;

hang = transform_origin(Trials{t},filter(copy(xyz{t}),'ButFilter',4,20),'hcom','nose',{'hcom','head_right'});

hbaCorrection = [-0.25,-0.25,                               ... er01
                 0.2,0.2,0.2,                       ... ER06
                 0,0,                               ... Ed10
                 0,0,0,0,0,0,0,0,0,                 ... jg04
                 -0.25.*ones([1,9]),                ... jg05
                0.2,0,-0.25]; % new units - jg05, jg05, ER06, Ed10, er01
hrlCorrection =[0,0,                                ... er01
                -0.27.*ones([1,3]),                 ... ER06
                -0.42,-0.42,                        ... Ed10
                -0.05.*ones([1,9]),                 ... jg04
                -0.48.*ones([1,9]),                 ... jg05
                -0.18,-0.42,0]; %                   ... ER06, Ed10, er01



headBodyAng = [xyz{t}(:,'spine_upper',[1,2])-xyz{t}(:,'bcom',[1,2]),...
               xyz{t}(:,'nose',[1,2])-xyz{t}(:,'hcom',[1,2])];
headBodyAng = sq(bsxfun(@rdivide,headBodyAng,sqrt(sum(headBodyAng.^2,3))));
headBodyAng = cart2pol(headBodyAng(:,:,1),headBodyAng(:,:,2));
headBodyAng = circ_dist(headBodyAng(:,2),headBodyAng(:,1));
headBodyAng = MTADfet.encapsulate(Trials{t},-(headBodyAng+hbaCorrection(t)),sampleRate,'hba','hba','h');

% $$$ 
hrl = copy(xyz{t});
hrl.data = -(hang.roll+hrlCorrection(t));
;
% $$$ figure;
% $$$ bar(linspace(-1,1,50),histc(hrl(Trials{t}.stc{'w'},:)+offset,linspace(-1,1,50)),'histc');
% $$$ 
% $$$ figure;
wper = Trials{t}.stc{'w'};
wper.cast('TimeSeries');
wper.resample(xyz{t});
%aind = (headBodyAng.data-hbaCorrection(t))<0.6 & (headBodyAng.data-hbaCorrection(t))>-0.6 & logical(wper.data);
aind = logical(wper.data) ;

figure();
subplot(131);
bar(linspace(-1,1,50),histc(hrl(aind,:)+offset,linspace(-1,1,50)),'histc');
subplot(132);
hist2([hrl(aind,:)+offset,headBodyAng(aind)],linspace(-1,1,30),linspace(-1.8,1.8,30))
subplot(133);
bar(linspace(-1.2,1.2,50),histc(headBodyAng(aind,:),linspace(-1.2,1.2,50)),'histc');

hrollBinInd = discretize((dhroll-0.18-0.2*double(~ismember(dtind,[3,4,5]))),hrollBinEdges);
cdhroll = (dhroll-0.18-0.2*double(~ismember(dtind,[3,4,5])));
hrollBinEdges = linspace(-0.5,0.5,8);
% egohrl
%%%>>>

%%%<<< SCRATCH
figure();
for a = 1:size(Smoother,3),
    subplot(size(Smoother,3),1,a);
    imagesc(Smoother(:,:,a)');
end

sum(nonzeros(Smoother(:,:,3)))

    

vxy = xyz{20}.vel('hcom',[1,2]);

figure,hist(log10(vxy([Trials{20}.stc{'walk+pause&theta'}])),100);

figure,hist((hvf([Trials{20}.stc{'walk+pause&theta'}],1)),100);

% SCRATCH 
%%%>>>


 
%%%<<< GENERATE spatial mask 

lims = {[-250,250],[-250,250]};
maskBinInd = cf(@(p,l) l(1) < p & p < l(2),  pfs{goodTrialInds(1)}{1}.adata.bins([1:2]),lims);
mask = maskBinInd{1};
for dim = 2:numel(maskBinInd),
     mask = bsxfun(@and,mask,permute(maskBinInd{dim},[2:dim,1]));
end
maskDimVal = cf(@(p,i) p(i), pfs{goodTrialInds(1)}{1}.adata.bins(1:2),maskBinInd);
maskBinPos = cell(size(maskDimVal));
[maskBinPos{:}] = ndgrid(maskDimVal{:});
maskBinPos = reshape(cat(ndims(maskBinPos)+1,maskBinPos{:}),[],ndims(maskBinPos));
%%%>>>

%%%<<< COMPUTE rate-weighted positions

rateMapNormal = pfs;
rateMapShuffled = pfsh;

rateMapNormal = pfl;
rateMapShuffled = pflh;

rateMapNormal = pfv;
rateMapShuffled = pfvh;

rateMapNormal = pfr;
rateMapShuffled = pfrh;

% COMPUTE the rate weighted field postion
rweightedPos = {};
trialAnatomicalGroup = {};

goodTrialInds = find(~cellfun(@isempty,units));

nAng = 5;
for trialId = goodTrialInds,
    rweightedPos{trialId} =                                                               ...
        nan([numel(units{trialId}),                                                       ...
             2,                                                                           ...
             numel(rateMapNormal{trialId}),                                               ...
             rateMapNormal{trialId}{1}.adata.binSizes(end)]);
    
    for unitId = 1:numel(units{trialId}),
        for thetaPhase = 1:numel(rateMapNormal{trialId}),
            rmap = plot(rateMapNormal{trialId}{thetaPhase},units{trialId}(unitId));
            rmap(rmap<1)=0;
            % Ratemap weighted position of egocetric placefield
            %     maskBinPos - coordinates of each bin
            %     ratemap normalized to 1
            for hbAngle = 1:nAng
                rweightedPos{trialId}(unitId,:,thetaPhase,hbAngle) =                      ...
                    sum(bsxfun(@times,                                                    ...
                               maskBinPos,                                                ...
                               reshape(rmap(maskBinInd{:},hbAngle),[],1)                  ...
                               ./sum(reshape(rmap(maskBinInd{:},hbAngle),[],1),'omitnan'))...
                        ,'omitnan');
            end%for hbAngle
        end%for thetaPhase
    end%for unitId
    trialAnatomicalGroup{trialId} = repmat(ismember(trialId,[1,2,6,7,26,27,28]),[numel(units{trialId}),1]);
end%for trialId

% rwpa
rateWeightedPosition = cat(1,rweightedPos{:});
trialAnatomicalGroup = cat(1,trialAnatomicalGroup{:});

rweightedPosSH = {};
nAng  = 5;
nIter = size(rateMapShuffled{goodTrialInds(1)}{1}.data.rateMap,3);
for trialId = goodTrialInds,
    rweightedPosSH{trialId} =                                                             ...
        nan([numel(units{trialId}),                                                       ...
             2,                                                                           ...
             numel(rateMapShuffled{trialId}),                                                   ...
             rateMapShuffled{trialId}{1}.adata.binSizes(end),                                   ...
             nIter]);
    for unitId = 1:numel(units{trialId}),
        for thetaPhase = 1:numel(rateMapShuffled{trialId}),
            clusterIndex = units{trialId}(unitId)==pfsh{trialId}{thetaPhase}.data.clu;
            rmap = reshape(pfsh{trialId}{thetaPhase}                                      ... MTAApfs object
                             .data                                                        ... data
                               .rateMap(:,                                                ... rate values
                                        clusterIndex,                                     ... cluster Index
                                        :),                                               ... samples
                           [],nAng,nIter);
            rmap(rmap<1)=0;            
            % Ratemap weighted position of egocetric placefield
            %     maskBinPos - coordinates of each bin
            %     ratemap normalized to 1
            rweightedPosSH{trialId}(unitId,:,thetaPhase,:,:) =                            ...
                permute(sum(bsxfun(@times,                                                ...
                                repmat(maskBinPos,[1,1,nAng,nIter]),                      ...
                                permute(rmap(mask(:),:,:),[1,4,2,3])                      ...
                                   ./ sum(permute(rmap(mask(:),:,:),[1,4,2,3]),           ...
                                          'omitnan')),                                    ...
                            'omitnan'),                                                   ...
                        [1,2,5,3,4]);
        end%for thetaPhase
    end%for unitId
end%for trialId



unitSessionMap = cf(@(u,i) u*0+i, units,num2cell(1:numel(units)));
unitSessionMap = cat(2,cat(2,unitSessionMap{:})',cat(2,units{:})');

%rwpaSH = cat(1,rweightedPosSH{:});
rateWeightedPositionShuffled = cat(1,rweightedPosSH{:});

zscr = (bsxfun(@minus,                                                                    ...
               rateWeightedPosition,                                                      ...
               mean(rateWeightedPositionShuffled,5)))                                     ...
       ./std(rateWeightedPositionShuffled,[],5);


zscrFWE = bsxfun(@minus,                                                                    ...
               rateWeightedPositionShuffled,                                              ...
               repmat(mean(rateWeightedPositionShuffled,5),[1,1,1,1,size(rateWeightedPositionShuffled,5)]))                                     ...
       ./repmat(std(rateWeightedPositionShuffled,[],5),[1,1,1,1,size(rateWeightedPositionShuffled,5)]);


p  = 1;
x  = 1;
al = 2;
ar = 4;
tida   = logical(trialAnatomicalGroup);
rwpa   = rateWeightedPosition;
rwpaSH = rateWeightedPositionShuffled;


%%%<<< FIGURES epp group stats

figDir = create_directory('/storage/share/Projects/EgoProCode2D/eppGroupStats');


nPhz = 5;

x  = 1;
nfn = @(x) ~x;
figure();
for p  = 1:nPhz;
    for a = 1:nAng,
        subplot2(5,5,a,p);
        hold('on');
        %bar(linspace(-150,200,30),hist(rwpa(tida,x,p,a)-25,linspace(-150,200,30)),'histc');
        plot(rwpa(nfn(tida),2,p,a),...
             rwpa(nfn(tida),1,p,a),...
             '.',...
             'MarkerSize',8);
        plot(mean(rwpa(nfn(tida),2,p,a)),...
             mean(rwpa(nfn(tida),1,p,a)),...
             '*m',...
             'MarkerSize',8);
        xlim([-100,100]);
        ylim([-150,200]);
        grid('on');
        %Lines(mean(rwpa(tida,x,p,a)-25),[],'r');
        if p == 1,
            if a == 3,
                %ylabel({'Forward Head Speed', [num2str(hvfBinCtr(a)),' mm/s']});
                %ylabel({'Head Roll', [num2str(hrlBinCtr(a)),' mm/s']});
            else
                %ylabel([num2str(hvfBinCtr(a)),' mm/s']);
                %ylabel([num2str(hrlBinCtr(a)),' rad']);
            end
        end
        if a == 1,
            if p == 3,
                title({'theta phase',num2str(binPhzc(p))});
            else,
                title(num2str(binPhzc(p)));
            end
        end
    end
end


% MEAN RostCaud erm position 
hfig = figure();
cmap = cool(5);
anaLabels = {'CA3','CA1'};

for i = 0:1;
hfig = figure();
cmap = cool(5);
anaLabels = {'CA3','CA1'};
    
    clf(hfig);
    grid('on');
    if logical(i),
        nfn = @(x) ~x;
    else,
        nfn = @(x) x;
    end;%if logical(i),
    
    lind = nfn(0)+1;

    hold('on');
    for p= 1:5,
        plot(hvfBinCtr,sq(mean(rwpa(nfn(tida),x,p,:))),'-*','Color',cmap(p,:));
    end
    ylabel('mm')
    xlabel('cm/s')
    legend(cf(@(t) ['\theta=',num2str(round(t,2))], num2cell(binPhzc(1:5))));
    xlim([-20,100]);ylim([-40,80]);
    
    title({'Mean Rostrocaudal Position of',...
           [anaLabels{lind},' Egofields as a function of'],...
           ['Theta phase and Forward Head Speed. n = ',num2str(sum(nfn(tida)))]});
    
% $$$     print(hfig,                                                               ...% figure handle
% $$$           '-dpng',                                                            ...% image format
% $$$           fullfile(figDir,['mean_RC_pos_Fn_thetaPhase_headVelocityFront_Hpc_',anaLabels{lind},'.png']));     % image path
end



% MEAN RostCaud erm position 
hfig = figure();
cmap = cool(5);
anaLabels = {'CA3','CA1'};

for i = 0:1;
    clf(hfig);
    grid('on');
    if logical(i),
        nfn = @(x) ~x;
    else,
        nfn = @(x) x;
    end;%if logical(i),
    
    lind = nfn(0)+1;

    hold('on');
    for p= 1:5,
        plot(hvfBinCtr,sq(mean(rwpa(nfn(tida),2,p,:))),'-*','Color',cmap(p,:));
    end
    ylabel('mm')
    xlabel('cm/s')
    legend(cf(@(t) ['\theta=',num2str(round(t,2))], num2cell(binPhzc(1:5))));
    xlim([-20,100]);ylim([-40,80]);
    
    title({'Mean Rostrocaudal Position of',...
           [anaLabels{lind},' Egofields as a function of'],...
           ['Theta phase and Forward Head Speed. n = ',num2str(sum(nfn(tida)))]});
    
    print(hfig,                                                               ...% figure handle
          '-dpng',                                                            ...% image format
          fullfile(figDir,['mean_ML_pos_Fn_thetaPhase_HeadBodyAngle_Hpc_',anaLabels{lind},'.png']));     % image path
end

%%%>>>




x = 1;
ny = 1;
figure();
%for p = 1:5,
ar=2;
ar=4;
for p = 2,
    y = 1;
    subplot2(ny,2,y,1);
    hold('on');
    grid('on');    
    % CA1
    sigThresh = p2z(1-(1-0.95)^(1/sum(~tida)));
    plot(zscr(~tida,x,p,al),zscr(~tida,x,p,ar),'.b','MarkerSize',8);
    plot(zscrFWE(~tida,x,p,al,1),zscrFWE(~tida,x,p,ar,1),'.m','MarkerSize',8);    
    plot(mean(zscr(~tida,x,p,al)),mean(zscr(~tida,x,p,ar)),'ob','MarkerSize',8);
    [P_CA1,S_CA1] = polyfit(zscr(~tida,x,p,al),zscr(~tida,x,p,ar),1);
    line([-6,6],polyval(P_CA1,[-6,6]),'Color','b','MarkerSize',8);
    [R_CA1, Pval_CA1] = corrcoef(zscr(~tida,x,p,al),zscr(~tida,x,p,ar));
    % Subplot Opts    
    xlim([-6,6]);
    ylim([-6,6]);
    Lines([],-sigThresh,'k','--','LineWidth',2);
    Lines(sigThresh,[],'k','--','LineWidth',2);
    line([-6,6],[6,-6],'Color','c');
    xlabel('leftward HBA (z-score)')
    ylabel('rightward HBA (z-score)')
    title('CA1 ');
    
    subplot2(ny,2,y,2);
    hold('on');
    grid('on');        
    % CA3
    sigThresh = p2z(1-(1-0.95)^(1/sum(tida)));
    plot(zscr(tida,x,p,al),zscr(tida,x,p,ar),'.g','MarkerSize',8);
    plot(zscrFWE(tida,x,p,al,1),zscrFWE(tida,x,p,ar,1),'.m','MarkerSize',8);        
    plot(mean(zscr(tida,x,p,al)),mean(zscr(tida,x,p,ar)),'og','MarkerSize',8);
    [P_CA3, S_CA3   ] = polyfit(zscr(tida,x,p,al),zscr(tida,x,p,ar),1);
    line([-6,6],polyval(P_CA3,[-6,6]),'Color','g');
    [R_CA3, Pval_CA3] = corrcoef(zscr(tida,x,p,al),zscr(tida,x,p,ar));
    % Subplot Opts        
    xlim([-6,6]);
    ylim([-6,6]);
    Lines([],-sigThresh,'k','--','LineWidth',2);
    Lines(sigThresh,[],'k','--','LineWidth',2);
    line([-6,6],[6,-6],'Color','c');
    xlabel('leftward HBA (z-score)')
    ylabel('rightward HBA (z-score)')
    title('CA3 ');
end

