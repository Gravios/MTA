function req20160411(stateIndex)
%function req20160411(Trial,varargin)
% tSNE over all jg05 sessions excluding jg05-20120317
% jg05-20120317 features are then mapped onto tsne space

% Figure Settings ----------------------------------------------------------------------
OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'MjgEd2016/manuscript/Figures/Figure_2';
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8);
% $$$       'defaultuicontrolunits','centimeters',...
% $$$       'defaultfigureunits','centimeters',...
% $$$       'defaultaxesunits','centimeters');

% --------------------------------------------------------------------------------------

% PARAMETERS ---------------------------------------------------------------------------
% stateIndex = 1;
iteration = '';
iteration = '_1';
states = {'rear','sit','groom','pause','walk','turn'};
Trial = MTATrial.validate('jg05-20120317.cof.all');
% jg05-20120317
% rear  opt 5
% sit   ori 5
% groom opt 6
% pause opt 5
% walk  opt 3
featureSet = 'fet_all';
featureSet = 'fet_mis';

mfilename = 'req20160411';
SAVEFIG = true;         % boolean flag, Save figures along with png's 
HOST_PATH = fullfile(getenv('PROJECT'),'figures'); % Where reportfig should
overwrite = false;


dsd = load(fullfile(Trial.spath,'req20160310_1_preproc-afet.mat'));
bs = load(fullfile(Trial.spath,'req20160310_5_genfigs.mat'));

sortedFeatureMaxIndex   = [5,5,6,5,3];
useOriginalFeatureOrder = [0,1,0,0,0];


ori = load(fullfile(Trial.spath,['req20160310_4_accumStats',num2str(stateIndex),'.mat']));    
opt = load(fullfile(Trial.spath,['req20160310_7_accumOptStats',num2str(stateIndex),'.mat']));
states = dsd.stateOrd(stateIndex:end);


if useOriginalFeatureOrder(stateIndex),
    featureSubsetInds = ori.sbind(1:sortedFeatureMaxIndex(stateIndex));
else
    featureSubsetInds = opt.sbind(1:sortedFeatureMaxIndex(stateIndex));
end

%stateIndex = 0;
switch stateIndex
  case 0, 
    featureSubsetInds = 1:18;
  case 1,
    featureSubsetInds = [3,5,14,10,16];
  case 2,
    featureSubsetInds = [7,8,9,10,4];
  case 3,
    featureSubsetInds = [7,10,4,1,2,18];
  case 4,
    featureSubsetInds = [7,11,12,14,6];
  case 5,
    featureSubsetInds = [17,13,6];
end


clear parameters
parameters.featureSet = featureSet; % set of features 
parameters.featureSubsetInds = featureSubsetInds; %rear
parameters.initial_dims = 3;       % number of dimensions of the feature matrix used
                                   % in preliminary rounds of t-SNE
parameters.perplexity = 80;        % arbitrary scale use to control the malability
                                   % of the distributions during dimension
                                   % reduction
parameters.newSampleRate = 12;     % Reduce the sample rate of the feature matrix
parameters.states = states;
parameters.normalize = true;       % Normalize features 
parameters.map2reference = true;   % Adjust features to "fit" the data between subjects
parameters.sessionSet = 'hand_labeled'; % List of sessions to include 
parameters.nSamples = 4000;        % Number of samples per state to pass to t-SNE
parameters.refTrial = 'jg05-20120317.cof.all';

% --------------------------------------------------------------------------------------

tag = DataHash(parameters);

% SET file name
fileLoc = fullfile(MTASession([]).path.project,'analysis',[mfilename,'-',tag,iteration,'.mat']);

% POPULATE variable from parameter struct
[featureSet,featureSubsetInds,initial_dims,perplexity,newSampleRate,states,...
normalize,map2reference,sessionSet,nSamples,refTrial] = DefaultArgs({},parameters,'--struct');



%Reference Trial Stuff
RefTrial = MTATrial.validate(refTrial);


if ~exist(fileLoc,'file')||overwrite,

    sessionList = get_session_list(sessionSet);

    cfet = [];
    Stc = {};
    sts = [];
    for s = 1:numel(sessionList),
        Trial = MTATrial.validate(sessionList(s));
        Stc = Trial.stc.copy;
        % Load Feature matrix of the session    
        
        % Load features
        if strcmp(Trial.filebase,RefTrial.filebase)&&map2reference,
            [tfet] = feval(featureSet,Trial,newSampleRate); ...
            [~,rMean,rStd] = unity(tfet);                
        else
            if strcmp('fet_all',featureSet)
                rt = RefTrial;
            else
                rt = [];
            end
            [tfet] = feval(featureSet,Trial,newSampleRate,rt);
            tfet.map_to_reference_session(Trial,RefTrial);
        end

        if normalize&&~strcmp(Trial.filebase,RefTrial.filebase)
            tfet.unity([],rMean,rStd);
        end
        
        % Resample features
        [sStc,~,sfet] = resample_whole_state_bootstrap_noisy_trim(Stc,tfet,states,90,nSamples);

        % Concatenate features
        if s == 1,
            fet = sfet.copy;
        else
            fet.data = cat(1,fet.data,sfet.data);
        end

        
        tsts = [];
        for t = states,
            tper = sStc{t{1}};
            tper.cast('TimeSeries');
            tper.resample(fet.sampleRate);
            tsts = cat(2,tsts,tper.data);
        end     
        % Inefficient ... I don't care
        sts = cat(1,sts,tsts);
        sts(~nniz(fet),:) = [];        
        fet.data(~nniz(fet),:) = [];
        cfet(end+1) = length(fet.data);
    end





    % Add colors to states
    asmat = MTADfet('data',sts,'sampleRate',fet.sampleRate);
    [~,asmat.data] = max(asmat.data,[],2);
    c = jet(numel(states));
    c = [1,0,0;...
         1,1,0;...
         1,0,1;...
         0,1,1;...
         0,0,1;...
         0,1,0];
    c = c(stateIndex:end,:);
    csmat = asmat.copy; 
    csmat.data = c(csmat.data,:);


    start = 1;
    skip = 3;
    stop = size(fet,1);
    no_dims = 2;
    rind = randperm(stop);
    rind(rind<start)=[];
    ind = rind(1:skip:end);
    

    mappedX = tsne(fet(ind,featureSubsetInds), csmat(ind,:), no_dims, initial_dims, perplexity); 

    
    save(fileLoc);
else
    load(fileLoc);
end



%% Figures (t-SNE) ----------------------------------------------------------------------------------

% PLOT t-SNE for given states

FigName = ['mis_tsne_' num2str(stateIndex) '_' states{1} '-' tag];

% figure properties
hfig = figure(3923924);clf;
hfig.PaperPositionMode = 'auto';
hfig.Units = 'centimeters';
hfig.Position = [0,0,8,8];
hold('on');

% loop over states and plot each with desiganated color
osts = numel(states);
mc = csmat(ind,:);
for nc = 1:osts,
    nind = all(bsxfun(@eq,c(nc,:),mc),2);
    h = scatter(mappedX(nind,1),mappedX(nind,2),2,mc(nind,:));
    try,h.MarkerFaceColor = h.CData(1,:);end
end

% axis properties
xlim([-110,110]),ylim([-110,110])
daspect([1,1,1])

legend(states,'location','SouthOutside','Orientation','horizontal');

% REPORT figure to analysis folder
print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% REPORT figure to standard analysis folder
% $$$ reportfig(HOST_PATH,           ...
% $$$           'FigHandle',hfig,    ...
% $$$           'FileName','tsne-stateXsession',   ...
% $$$           'FigDir','req',      ...
% $$$           'Preview',false,     ...
% $$$           'Tag',[sessionSet,'-',featureSet],...
% $$$           'Comment',['tsne-all'],    ...
% $$$           'Resolution',100,    ...
% $$$           'SaveFig', SAVEFIG,  ...
% $$$           'format','png',      ...
% $$$           'width',12,           ...
% $$$           'height',8);


%% Figures (t-SNE State vs Session) ------------------------------------------------------
FigName = ['mis_tsne_' num2str(stateIndex) '_' states{1} '_SessionMapping' '-' tag];
nsts = numel(states);
nses = numel(sessionList);

hfig = figure(3923925);clf;
hfig.PaperPositionMode = 'auto';
hfig.Units = 'centimeters';
hfig.Position = [0,0,10*nsts,8];
hold('on');



cses = [0,0,1;...
        1,0,0;...
        0,1,0;...
        1,0,1;...
        0,1,1;...
        1,1,0;];
        
sesIdMat = MTADfet('data',reshape(bsxfun(@times,ones([nSamples*nsts,nses]),1:nses),[],1),...
                   'sampleRate',fet.sampleRate);

cSesMat = sesIdMat.copy; 
cSesMat.data = cses(cSesMat.data,:);
mSesC = cSesMat(ind,:);
mc = csmat(ind,:);
for stsId = 1:nsts;
    subplot(1,nsts,stsId);
    hold on;    
    for nc = 1:nses,
        nind = all(bsxfun(@eq,cses(nc,:),mSesC),2)&all(bsxfun(@eq,c(stsId,:),mc),2);    
        h = scatter(mappedX(nind,1),mappedX(nind,2),2,mSesC(nind,:));
        h.MarkerFaceColor = h.CData(1,:);
    end
    title(states(stsId))
    daspect([1,1,1]);
    xlim([-150,150])
    ylim([-125,125])
end

hl = legend({sessionList.sessionName},'location','SouthOutside','Orientation','horizontal');
hl.Position(2) = 0.01;

% REPORT figure to project folder
print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));



%% Figures (t-SNE feature map) ------------------------------------------------------.
FigName = ['mis_tsne_' num2str(stateIndex) '_' states{1} '_FeatureMapping' '-' tag];

mtfet =  MTADxyz('data',fet(ind,:),'sampleRate',fet.sampleRate);
mtpos =  MTADxyz('data',mappedX,'sampleRate',fet.sampleRate);


[~,fett,fetd] = fet_mis(Trial);

hfig = figure(38381);
hfig.Units = 'centimeters';
hfig.Position = [1,1,14,8];
for i = featureSubsetInds,
    mtfet = MTADxyz('data',fet(ind,i),'sampleRate',fet.sampleRate);
[meanMap,Bins,distdw,distIndw]= PlotKNNPF(Trial,mtfet,mtpos,[5,5],20,5,'xy',[],[],[-110,125;-125,110],@nanmean);
mmap = reshape(meanMap,numel(Bins{1}),numel(Bins{2}),[]);
[stdMap,Bins,distdw,distIndw]= PlotKNNPF(Trial,mtfet,mtpos,[5,5],20,5,'xy',[],[],[-110,125;-125,110],@nanstd);
smap = reshape(stdMap,numel(Bins{1}),numel(Bins{2}),[]);
    
    clf
    hax = axes;
    hax.Units = 'centimeters';
    hax.Position = [2,2,4,4];
    [ha,hc] = imagescnan({Bins{1},Bins{2},mmap(:,:)},[],false, ...
                        true,[0,0,0]);
    ylabel(hc,'mean z-score');
    axis xy
    title(fett{i})
    daspect([1,1,1])
    
    hax = axes;
    hax.Units = 'centimeters';
    hax.Position = [8,2,4,4];    
    [ha,hc] = imagescnan({Bins{1},Bins{2},smap(:,:)},[],false, ...
                        true,[0,0,0]);
    ylabel(hc,'std z-score');
    axis xy
    title(fett{i})
    daspect([1,1,1])

    
% REPORT figure to project folder
print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'_',num2str(i),'.eps']));
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'_',num2str(i),'.png']));
   
%reportfig([], hfig, 'tsne', 'req',false,fett{i},fetd{i},[],SAVEFIG,'png',8,8);
end



%% JPDF with contours


[~,featureTitles] = fet_mis(RefTrial);

csmat.data(round(linspace(nSamples+1,size(csmat,1),numel(states)))',:) = [];
targetState = sStc.states{1};
targetState.data = ThreshCross(all(bsxfun(@eq,csmat.data,c(1,:)),2),0.5,10);

remainingStates = sStc.states{end};                      
remainingStates.data =  ThreshCross(~all(bsxfun(@eq,csmat.data,c(1,:)),2),0.5,10);
remainingStates.data(end,2) = size(fet,1);
nbin = 50;                                        % number of bins per axis in JPDF


for f = 1:numel(featureSubsetInds),
    for k = f+1:numel(featureSubsetInds),
        x = featureSubsetInds(f);                         % fet_mis feature index
        y = featureSubsetInds(k);                         % fet_mis feature index
        stateThreshold = 10e-4;
        astateThreshold = 10e-4;
        patchPercentState = 0;
        patchPercentAState = 0;
        display = [true,true];
        FigName = ['mis_jpdf_' num2str(stateIndex) '_' states{1} '_FeatureJPDF_'...
                   num2str(x) '-' num2str(y) '-' tag];



        for i = 1:2
            display = ~display;
            hfig = figure(383812);clf;hold on;delete(gca);    % Create and clear the figure
            set(hfig,'Units','centimeters',...
                     'PaperPositionMode', 'auto')
            hfig.Position(3:4) = [6,6];
            hind = nniz(fet);                                % Index for not zero, inf or nan           

            % SET bins: try refined edgs See (featureAxLim) Definitions below or percentile
            edgs    = {linspace([prctile(fet(hind,x),[2,98])+[-1,1],nbin])};    
            edgs(2) = {linspace([prctile(fet(hind,y),[2,98])+[-1,1],nbin])};    
            edc = cell(size(edgs));
            [edc{:}] = get_histBinCenters(edgs);               


            % PLOT JPDF of features
            hax = axes('Units','centimeters','Position',[2,2,2,2]); hold on;
            [histMat,~,~,histBinsX,histBinsY] = ...
                histcounts2(fet(hind,x),fet(hind,y),...
                            edgs{1},edgs{2},...
                            'Normalization','probability');
            imagesc(edgs{1},edgs{2},histMat');                     % Plot JPDF of hand labeled data
            axis xy;
            axis tight;

            % COMPUTE JPDF for target state
            while patchPercentState < 0.95 || display(1),
                [stateOverlay,~,~,stateOverlayBinsX,stateOverlayBinsY] = ...
                    histcounts2(fet(targetState,x),fet(targetState,y),...
                                edgs{1},edgs{2},...
                                'Normalization','probability');
                stateOverlay = imgaussfilt(stateOverlay,2);
                stateOverlayThresholdMatrix = stateOverlay>stateThreshold;
                [stateOverlayBoundaries,statePatchMatrix] = bwboundaries(stateOverlayThresholdMatrix);
                for b = stateOverlayBoundaries'
                    plot(edc{1}(stateOverlayBoundaries{1}(:,1)),edc{2}(stateOverlayBoundaries{1}(:,2)),'-r');
                end    
                nzb = stateOverlayBinsX&stateOverlayBinsY;
                patchPercentState = sum(statePatchMatrix(...
                    sub2ind(size(statePatchMatrix),...
                            stateOverlayBinsX(nzb),...
                            stateOverlayBinsY(nzb)))...
                                        ==1)/sum(nzb);
                if patchPercentState < 0.95,
                    stateThreshold = stateThreshold - 0.5e-4;
                else
                    display(1) = false;
                end
            end


            % COMPUTE JPDF for all states excluding target state
            while patchPercentAState < 0.95 || display(2),
                [astateOverlay,~,~,astateOverlayBinsX,astateOverlayBinsY] = ...
                    histcounts2(fet(remainingStates,x),fet(remainingStates,y),...
                                edgs{1},edgs{2},...
                                'Normalization','probability');
                astateOverlay = imgaussfilt(astateOverlay,2);
                astateOverlayThresholdMatrix = astateOverlay>astateThreshold;
                [astateOverlayBoundaries,astatePatchMatrix] = bwboundaries(astateOverlayThresholdMatrix);
                for b = astateOverlayBoundaries'
                    plot(edc{1}(astateOverlayBoundaries{1}(:,1)),edc{2}(astateOverlayBoundaries{1}(:,2)),'-c');
                end    
                nzb = astateOverlayBinsX&astateOverlayBinsY;
                patchPercentAState = sum(astatePatchMatrix(...
                    sub2ind(size(astatePatchMatrix),...
                            astateOverlayBinsX(nzb),...
                            astateOverlayBinsY(nzb)))...
                                         ==1)/sum(nzb);
                if patchPercentAState < 0.95,
                    astateThreshold = astateThreshold - 0.5e-4;
                else
                    display(2) = false;
                end
            end

            caxis([0,0.003])
            xlabel(featureTitles{x});
            ylabel(featureTitles{y});
        end

        print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
        print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
    end
end







%% PART 2

% compute t-SNE mapping on full feature space
pind = all(bsxfun(@eq,csmat.data,[1,0,1]),2);
subMappedX = tsne(fet(pind,:), csmat(pind,:), no_dims, initial_dims, perplexity); 


% Load features of a target session
Trial = MTATrial.validate('jg05-20120317.cof.all');
xyz = Trial.load('xyz');
rt = []; %if strcmp('fet_all',featureSet), rt = RefTrial; else, rt = []; end

[pfet] = feval(featureSet,Trial,newSampleRate,rt);
if map2reference, pfet.map_to_reference_session(Trial,RefTrial); end
if normalize,     pfet.unity([],rMean,rStd); end


gind = Trial.stc{'m'};
%gInd.cast('TimeSeries',pfet);

hfig = figure(gen_figure_id);
plot(subMappedX(:,1),subMappedX(:,2),'.');
cIds = ClusterPP(hfig);


distMat = sqrt(sum(bsxfun(@minus,...
                          permute(pfet(gind,:),[1,3,2]),...
                          permute(fet(pind,:),[3,1,2])).^2,3));
[~,closestDistInd] = sort(distMat,2);

groomCluId = gind.copy;
groomCluId.cast('TimeSeries',pfet);
groomCluId.data = double(groomCluId.data);
groomCluId.data(groomCluId.data==1) = mode(cIds(closestDistInd(:,1:100)),2);
groomCluId.resample(xyz);



% $$$ reportfig(HOST_PATH,           ...
% $$$               'FigHandle',hfig,    ...
% $$$               'FileName','tsne-stateXsession',   ...
% $$$               'FigDir','req',      ...
% $$$               'Preview',false,     ...
% $$$               'Tag',[sessionSet,'-',featureSet],...
% $$$               'Comment',['tsne-',states{stsId}],    ...
% $$$               'Resolution',100,    ...
% $$$               'SaveFig', SAVEFIG,  ...
% $$$               'format','png',      ...
% $$$               'width',12,           ...
% $$$               'height',8);



% $$$ osts = numel(states);
% $$$ hfig = figure(3923924);clf
% $$$ hold on;
% $$$ mc = csmat(ind,:);
% $$$ for nc = 1:osts,
% $$$     nind = all(bsxfun(@eq,c(nc,:),mc),2);
% $$$     h = scatter(mappedX(nind,1),mappedX(nind,2),2,mc(nind,:));
% $$$     try,h.MarkerFaceColor = h.CData(1,:);end
% $$$ end
% $$$ legend(states,'location','south','Orientation','horizontal');
% $$$ reportfig(HOST_PATH,           ...
% $$$           'FigHandle',hfig,    ...
% $$$           'FileName','tsne',   ...
% $$$           'FigDir','req',      ...
% $$$           'Preview',false,     ...
% $$$           'Tag',[sessionSet,'-',featureSet,'-nomap'],...
% $$$           'Comment','tsne',    ...
% $$$           'Resolution',100,    ...
% $$$           'SaveFig', SAVEFIG,  ...
% $$$           'format','png',      ...
% $$$           'width',8,           ...
% $$$           'height',8);
% $$$ 
% $$$ 

%mkdir('/storage/gravio/figures/req/req20151203');