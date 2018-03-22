% MjgER2016 Figure3
%
% Subplots:
%    A. autoccgs
%    B. Placefield examples for each state
%    C. 
%    D. Feature Matrix
%    E. Behavioral labels - hand labels and neural network labels



sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';

FigDir = '/storage/gravio/figures/placefields';
mkdir(FigDir);

% SET behvaioral states
states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',    ...
          'pause&theta','lpause&theta','hpause&theta',           ...
          'theta-groom-sit'};
statesCcg = {'loc','lloc','hloc','rear',    ...
          'pause','lpause','hpause',           ...
          'theta-groom-sit'};

 
% LOAD Trials
% COMPUTE placefield statistics
Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
          cf(@(t)  t.load('nq'),           Trials);
pfstats = cf(@(t)  compute_pfstats_bs(t),       Trials);


% $$$ cf(@(t)  MjgER2016_drzfields(t,true), Trials);




patchComparison = [];

% FOR each Trial -------------------------------------------------------------
for t = 1:20,

% LOAD Trial
    Trial = Trials{t};
    units = pfstats{t}.cluMap;    
    disp(['Processing trial: ' Trial.filebase]);

% LOAD placefields and subsampled estimate
    for sts = 1:numel(states),        
        defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
        defargs.units = units;
        %defargs.overwrite = true;
        defargs.states = states{sts};
        defargs = struct2varargin(defargs);        
        pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
    end

% COMPARE patches
    for u = 1:numel(units),
        unit = units(u);
        for s = 1:8,
            for o = 1:8,
                rmap = plot(pfkbs{o},unit,'mean');
                patchInd = sq(pfstats{t}.pfmstats(s,unit==units).patchRateInd(1,1,1,:,:));
                oPatchRates = rmap(sub2ind(size(rmap),patchInd(1,nniz(patchInd(1,:)'))',...
                                           patchInd(2,nniz(patchInd(2,:)'))'));
                if s == 1 && o ==1, pcuInd = size(patchComparison,1)+1; end
                
                if sum(nniz(oPatchRates))>0.9*numel(oPatchRates), 
                    patchComparison(pcuInd,s,o) = mean(oPatchRates,'omitnan');
                    patchComparison(pcuInd,s,o) = mean(oPatchRates,'omitnan');                    
                else,
                    patchComparison(pcuInd,s,o) = nan;
                end
                    
            end%for o
        end%for s 
    end%for u
end%for t

        
figure();
subplot(141); plot(patchComparison(:,4,1),patchComparison(:,4,4),'.'); daspect([1,1,1]);xlim([0,35]);ylim([0,35]);grid('on');
xlabel('mean loc rate in rear');ylabel('mean rear rate in rear');
subplot(142); plot(patchComparison(:,1,1),patchComparison(:,1,4),'.'); daspect([1,1,1]);xlim([0,35]);ylim([0,35]);grid('on');
subplot(143); plot(patchComparison(:,2,2),patchComparison(:,2,3),'.'); daspect([1,1,1]);xlim([0,35]);ylim([0,35]);grid('on');
subplot(144); plot(patchComparison(:,3,2),patchComparison(:,3,3),'.'); daspect([1,1,1]);xlim([0,35]);ylim([0,35]);grid('on');

figure();
subplot(121); plot(patchComparison(:,8,1),patchComparison(:,8,4),'.'); daspect([1,1,1]);xlim([0,35]);ylim([0,35]);grid('on');
subplot(122); plot(patchComparison(:,8,2),patchComparison(:,8,3),'.'); daspect([1,1,1]);xlim([0,35]);ylim([0,35]);grid('on');
title;xlabel('low loc');ylabel('high loc');


figure();
subplot(
hist((patchComparison(:,8,2)-patchComparison(:,8,3))./(patchComparison(:,8,2)+patchComparison(:,8,3)),50);
        
% $$$ 
% $$$     
% $$$ 
% $$$     
% $$$     
% $$$ % LOAD DRZ fields
% $$$     dfs = cell([1,3]);
% $$$     [dfs{:}] = MjgER2016_drzfields(Trial,true);
% $$$     dfst = {'pitch','height','rhm'};
% $$$ 
% $$$ 
% $$$ 
% $$$     %[accg,tbins] = autoccg(Trial);
% $$$ 
% $$$     % PLACEFIELD figure now located in placefield_summary.m
% $$$ end
% $$$ 
% $$$ pr = reshape(cell2mat(cf(@(p) [p.pfmstats.peakFR],pfstats)),8,[]);
% $$$ 
% $$$ 
% $$$ figure();
% $$$ 
% $$$ clf();
% $$$ % loc x rear
% $$$ subplot(221);  hold('on');  plot(pr(1,:),pr(4,:),'.');  daspect([1,1,1]);
% $$$                line([0,40],[0,40]);  line([0,40],[0,20]); line([0,20],[0,40]);
% $$$ subplot(222);  hold('on');  plot(pr(2,:),pr(4,:),'.');   daspect([1,1,1]);
% $$$                line([0,40],[0,40]);  line([0,40],[0,20]); line([0,20],[0,40]);
% $$$                
% $$$ subplot(224);  hold('on');  plot(pr(3,:),pr(4,:),'.');   daspect([1,1,1]);
% $$$                line([0,40],[0,40]);  line([0,40],[0,20]); line([0,20],[0,40]);               





% Create parts for figure 3 place fields





% selected units session x unit id
clumap = [17,147;...
           3,171;...
           3,197;...
          16, 51;...
           1, 15;...
          17,165;...
          18,129;...
         ];

FigDir = create_directory('/storage/gravio/figures/analysis/parts/MjgER2016/');

% LOAD session list
sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);


% SET states to plot
states = {'rear&theta','hloc&theta','lloc&theta',     ...
          'hpause&theta','lpause&theta',            ...
          'theta-groom-sit'};
numStates = numel(states);


pfdVersion = '7';


interpPar = struct('bins',{{linspace(-500,500,100),linspace(-500,500,100)}},             ...
                   'nanMaskThreshold', 0,                                                ...
                   'methodNanMap',     'linear',                                         ...
                   'methodRateMap',    'linear');


% SET figure opts
hfig = figure(666001);
hfig.Units = 'centimeters';
hfig.Position = [1, 1, 40, 5];
hfig.PaperPositionMode = 'auto';

% SET number of subplots
ny = 12;


for u = clumap',
    
    clf();    
    
    sp = gobjects([1,0]);
    
    Trial = MTATrial.validate(sessionList(u(1)));


% LOAD place fields
    pft = pfs_2d_theta(Trial,u(2));
    pft.parameters.states = 'theta-groom-sit';
    pfs = cat(2,pfs_2d_states(Trial,u(2)),{pft});

    
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
    


% LOAD DRZ fields
    dfs = req20180123_ver5(Trial,[],pfdVersion);
    dfst = {'HPITCHxBPITCH','HPITCHxBSPEED','BPITCHxBSPEED','BPITCHxHSPEED','HPITCHxRHM'};

% LOAD accgs
    [accg,tbins] = autoccg(Trial);


% SET color scale max
        maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                                      [pfs,dfs],repmat({u(2)},[1,numel(pfs)+numel(dfs)]))));

% $$$     maxPfsRate = max([maxRate(pfs{1},u(2),false,'prctile99',0.5),...
% $$$                       maxRate(pfs{2},u(2),false,'prctile99',0.5),...
% $$$                       maxRate(pfs{3},u(2),false,'prctile99',0.5),...
% $$$                       maxRate(pfs{4},u(2),false,'prctile99',0.5),...
% $$$                       maxRate(pfs{5},u(2),false,'prctile99',0.5),...
% $$$                       maxRate(pfs{6},u(2),false,'prctile99',0.5),...
% $$$                       maxRate(dfs{1},u(2),false,'prctile99',0.5),...                          
% $$$                       maxRate(dfs{2},u(2),false,'prctile99',0.5),...
% $$$                       maxRate(dfs{3},u(2),false,'prctile99',0.5),...                          
% $$$                       maxRate(dfs{4},u(2),false,'prctile99',0.5),...
% $$$                       maxRate(dfs{5},u(2),false,'prctile99',0.5)]);


    pfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,true,'mean')),pfs,repmat({u(2)},[1,numel(pfs)])));
    dfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,false,'mean')),dfs,repmat({u(2)},[1,numel(dfs)])));


% ACCG 
    sp(end+1) = subplot(1,ny,1);
    bar(tbins,accg(:,u(2)));axis tight;
    set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});
    
    for s = 1:numStates,
% PLACEFIELDS MTAApfs
        sp(end+1) = subplot(1,ny,s+1);
        plot(pfs{s},u(2),'mean',false,[0,maxPfsRate],true,0.5,false,interpPar,@jet);
        set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});
        title(sprintf('Max Rate: %3.2f',pfsMaxRatesMean(s)));
    end

% DRZFIELDS 
    for s = 1:numel(dfs),
        sp(end+1) = subplot(1,ny,s+numStates+1);
        dfs{s}.plot(u(2),'mean',false,[0,maxPfsRate],false,0.5,false,[],@jet);
        title(sprintf('Max Rate: %3.2f',dfsMaxRatesMean(s)));    
    end

    af(@(h) set(h,'FontSize',4), sp);
    af(@(h) set(h,'Units','centimeters'), sp);
    af(@(h) set(h,'Position',[h.Position(1:2),1,1]), sp);

    FigName = ['MjgER2016_figure3_parts_pfs_dfs','_',Trial.filebase,'_unit_',num2str(u(2))];
    print(hfig,'-depsc2',fullfile(FigDir,[FigName,'.eps']));        
    print(hfig,'-dpng',  fullfile(FigDir,[FigName,'.png']));
    
end


fetSets{end+1}     = 'fet_HB_pitchB';
fetSets{end+1}     = 'fet_HB_HPS';
fetSets{end+1}     = 'fet_HB_HPS';
fetSets{end+1}     = 'fet_HB_HPS';
fetSets{end+1}     = 'fet_HB_HPR';