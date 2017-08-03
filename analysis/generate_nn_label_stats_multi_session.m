function [out] = generate_nn_label_stats_multi_session(mode,varargin)
%function [out] = generate_nn_label_stats_multi_session(varargin)
%
% DEFARGS ----------------------------------------------------------------------
%
%    'sessionList',                 '',                                            
%    'featureSet',                  'fet_bref_rev7',                               
%    'states',                      {{'walk','rear','turn','pause','groom','sit'}},
%    'model',                       [],                                            
%    'sampleRate',                  10,                                            
%    'nNeurons',                    25,                                            
%    'nIter',                       10,                                            
%    'randomizationMethod',         'WSBNT',                                       
%    'map2reference',               true,                                          
%    'normalize',                   true,                                          
%    'referenceTrial',              'jg05-20120317.cof.all',                       
%    'trainingSessionList',         'hand_labeled',                                
%    'normalizationSessionList',    'hand_labeled',                                
%    'dropIndex',                   [],                                            
%    'prctTrain',                   90,                                            
%    'stcTag',                      'msnn'                                         
%    'postProcessingTag'            '' 

global MTA_PROJECT_PATH

MODEL_TYPE = 'multiSesPatNet';

% DEFARGS ----------------------------------------------------------------------------------------
defargs = struct('sessionList',                 '',                                            ...
                 'featureSet',                  'fet_bref_rev7',                               ...
                 'states',                      {{'walk','rear','turn','pause','groom','sit'}},...
                 'keys',                        {{'w','r','n','p','m','s'}},                   ...
                 'model',                       [],                                            ...
                 'sampleRate',                  10,                                            ...
                 'nNeurons',                    25,                                            ...
                 'nIter',                       10,                                            ...
                 'randomizationMethod',         'WSBNT',                                       ...
                 'map2reference',               true,                                          ...
                 'normalize',                   true,                                          ...
                 'referenceTrial',              'jg05-20120317.cof.all',                       ...
                 'trainingSessionList',         'hand_labeled',                                ...
                 'normalizationSessionList',    'hand_labeled',                                ...
                 'dropIndex',                   [],                                            ...
                 'prctTrain',                   90,                                            ...
                 'postProcessingTag',           ''                                             ...
);
[sessionList,featureSet,states,keys,model,sampleRate,nNeurons,nIter,randomizationMethod,       ...
 map2reference,normalize,referenceTrial,trainingSessionList,normalizationSessionList,          ...
 dropIndex,prctTrain,postProcessingTag] = DefaultArgs(varargin,defargs,'--struct');
% ------------------------------------------------------------------------------------------------    



% MAIN -------------------------------------------------------------------------------------------

switch mode
% TRAIN ------------------------------------------------------------------------------------------
  case 'train' 
    trainingList = get_session_list(trainingSessionList); 
    for dropIndex = 1:numel(trainingList),
        bhv_nn_multi_session_patternnet(...
            sessionList,         featureSet, states, keys, model,                             ...
            sampleRate,          nNeurons,   nIter,  randomizationMethod,                     ...
            map2reference,       normalize,  referenceTrial,                                  ...
            trainingSessionList, normalizationSessionList,                                    ...
            dropIndex,           prctTrain                                                    ...
        );
    end


    
% COMPUTE ----------------------------------------------------------------------------------------
  case 'compute'


    if isempty(sessionList),
        sessionList = trainingSessionList; 
        sesList = get_session_list(sessionList);
        stc           = cell([1,numel(sesList)]); 
        networkOutput = cell([1,numel(sesList)]);
        labelingStats = cell([1,numel(sesList)]);
        for dropIndex = 1:numel(sesList),
            [stc(dropIndex),labelingStats{dropIndex},networkOutput{dropIndex}] =             ...
                bhv_nn_multi_session_patternnet(                                             ...
                    sesList(dropIndex),     featureSet,   states,   keys, model,             ...
                    sampleRate,             nNeurons,     nIter,    randomizationMethod,     ...
                    map2reference,          normalize,    referenceTrial,                    ...
                    trainingSessionList,    normalizationSessionList,                        ...
                    dropIndex,              prctTrain                                        ...
            );
        end
    else
        sesList = get_session_list(sessionList);
        stc           = cell([1,numel(sesList)]); 
        networkOutput = cell([1,numel(sesList)]);
        labelingStats = cell([1,numel(sesList)]);        
        for s = 1:numel(sesList),
            [stc(s),labelingStats{s},networkOutput(s)] = bhv_nn_multi_session_patternnet(    ...
                sessionList,            featureSet,   states,   keys,  model,                ...
                sampleRate,             nNeurons,     nIter,    randomizationMethod,         ...
                map2reference,          normalize,    referenceTrial,                        ...
                trainingSessionList,    normalizationSessionList,                            ...
                dropIndex,              prctTrain                                            ...
            );   
        end
    end

% CREATE statitics filename
    statsName = ['MTAC_STATS+TRN+' trainingSessionList                                       ...
                '+LBS+' sessionList '+'                                                      ...
                featureSet                                                                   ...
                '+SR'  num2str(sampleRate)                                                   ...
                'NN'  num2str(nNeurons)                                                      ...
                'NI'   num2str(nIter)                                                        ...
                'M'    num2str(map2reference)                                                ...
                'MREF+' referenceTrial                                                       ...
                '+N'    num2str(normalize)                                                   ...
                'NREF+' normalizationSessionList                                             ...
                '+RND' randomizationMethod                                                   ...
                '+PRCT' num2str(prctTrain)                                                   ...
                '+STS+' strjoin(keys,'')                                                     ...
                '-'    MODEL_TYPE];
 
% SAVE labeling statitics
    save(fullfile(MTA_PROJECT_PATH,'analysis',[statsName,'.mat']),                           ...
         '-v7.3','nNeurons','nIter','sampleRate','featureSet',                               ...
         'randomizationMethod','states','stc','labelingStats',                               ...
         'networkOutput','trainingSessionList','sessionList','model'                         ...
    );

    
    
% OPTIMIZE --------------------------------------------------------------------------------------
  case 'optimize'
    sessionListWasEmpty = false;
    if isempty(sessionList),
        sessionList = trainingSessionList; 
        sessionListWasEmpty = true;
    end
    
% LOAD neural network labeled state collections
    statsName = ['MTAC_STATS+TRN+' trainingSessionList                                       ...
                '+LBS+' sessionList '+'                                                      ...
                featureSet                                                                   ...
                '+SR'   num2str(sampleRate)                                                  ...
                'NN'    num2str(nNeurons)                                                    ...
                'NI'    num2str(nIter)                                                       ...
                'M'     num2str(map2reference)                                               ...
                'MREF+' referenceTrial                                                       ...
                '+N'    num2str(normalize)                                                   ...
                'NREF+' normalizationSessionList                                             ...
                '+RND'  randomizationMethod                                                  ...
                '+PRCT' num2str(prctTrain)                                                   ...
                '+STS+' strjoin(keys,'')                                                     ...
                '-'     MODEL_TYPE];
    disp(statsName);
    stsRaw = load(fullfile(MTA_PROJECT_PATH,'analysis',[statsName,postProcessingTag,'.mat']));

    
    if sessionListWasEmpty,
        sessionList = trainingSessionList; 
        sesList = get_session_list(sessionList);
        stc           = cell([1,numel(sesList)]);         
        for dropIndex = 1:numel(sesList),
            stc(dropIndex) = optimize_stc_transition(stsRaw.stc(dropIndex),                  ...
                    sesList(dropIndex),     featureSet,   states,   keys, model,             ...
                    sampleRate,             nNeurons,     nIter,    randomizationMethod,     ...
                    map2reference,          normalize,    referenceTrial,                    ...
                    trainingSessionList,    normalizationSessionList,                        ...
                    dropIndex,              prctTrain                                        ...
            );
        end
    else
        sesList = get_session_list(sessionList);
        stc           = cell([1,numel(sesList)]); 
        networkOutput = cell([1,numel(sesList)]);
        labelingStats = cell([1,numel(sesList)]);        
        for s = 1:numel(sesList),
            [Stc(s),labelingStats{s},networkOutput(s)] = bhv_nn_multi_session_patternnet(    ...
                sessionList,            featureSet,   states,   keys,  model,                ...
                sampleRate,             nNeurons,     nIter,    randomizationMethod,         ...
                map2reference,          normalize,    referenceTrial,                        ...
                trainingSessionList,    normalizationSessionList,                            ...
                dropIndex,              prctTrain                                            ...
            );   
        end
    end

% REMOVE states disjoint to target state set
    cf(@(stc,states) set(stc,'states',stc(states{:})), stc,repmat({states},[1,numel(stc)]));
% COMPUTE cross label statistics between optimized and hand labeled state collections
    labelingStats = compute_inter_stc_stats(sessionList,stc,states,sampleRate);

% SAVE labeling statitics    
    save(fullfile(MTA_PROJECT_PATH,'analysis',[statsName,'_pp.mat']),                        ...
         '-v7.3','stc','labelingStats','nNeurons','nIter','sampleRate','featureSet',         ...
         'randomizationMethod','states','trainingSessionList','sessionList'                  ...
     );


    
    
    
    
% DISPLAY --------------------------------------------------------------------------------------
  case 'display'

% DISPLAY - Load Data and Set Parameters --------------------------------------------------------    
% SET figure save paths
    OwnDir = '/storage/gravio/ownCloud/';
    FigDir = 'Shared/Behavior Paper/Figures/Figure_1/parts';
    try,mkdir(fullfile(OwnDir,FigDir));end    

% SET figure dimensions
    fig.w = 6;
    fig.h = 4;
    fig.units = 'centimeters';
    ax.w = 4;
    ax.h = 3;
    ax.units = 'centimeters';

% LOAD labeling statistics
    if isempty(sessionList), sessionList = trainingSessionList;  end
    statsName = ['MTAC_STATS+TRN+' trainingSessionList                                       ...
                '+LBS+' sessionList '+'                                                      ...
                featureSet                                                                   ...
                '+SR'   num2str(sampleRate)                                                  ...
                'NN'    num2str(nNeurons)                                                    ...
                'NI'    num2str(nIter)                                                       ...
                'M'     num2str(map2reference)                                               ...
                'MREF+' referenceTrial                                                       ...
                '+N'    num2str(normalize)                                                   ...
                'NREF+' normalizationSessionList                                             ...
                '+RND'  randomizationMethod                                                  ...
                '+PRCT' num2str(prctTrain)                                                   ...
                '+STS+' strjoin(keys,'')                                                     ...
                '-'     MODEL_TYPE];
    stsRaw = load(fullfile(MTA_PROJECT_PATH,'analysis',[statsName,postProcessingTag,'.mat']));
    stsOpt = load(fullfile(MTA_PROJECT_PATH,'analysis',[statsName,postProcessingTag,'_pp.mat']));

% DISPLAY - Plot Accuracy -----------------------------------------------------------------------

% CREATE figure        
    hfig = figure();
    clf();
    hold('on');
    
% APPLY figure properties
    set(hfig,'PaperPositionMode','auto');
    set(hfig,'units',fig.units);
    set(hfig,'Position',[0,0,fig.w,fig.h]);
    
% FORMAT statistical data for boxplot
    dat = [cell2mat(cf(@(x)x.accuracy,stsRaw.labelingStats));...
           cell2mat(af(@(x)x.accuracy,stsOpt.labelingStats))]'.*100;
    grp = bsxfun(@times,ones(size(dat)),[1:size(dat,2)]);
    
% PLOT boxplot
    boxplot(dat(:),grp(:));
    h = findobj(gca,'Tag','Box');
    
% OVERLAY each boxplot with patch
    for j=1:length(h)
            p = patch(get(h(j),'XData'),...
                      get(h(j),'YData'),...
                      [0.5,1,0.75].*double([rem(j,3)==1,rem(j,3)==2,rem(j,3)==0]),...
                      'FaceAlpha',.5);
% REORDER patches behind boxplots
            uistack(p,'bottom');
    end

% FORMAT figure
    title('accuracy')
    ylim([50,100]);
    grid('on');
    set(gca,'gridlinestyle',':');
    set(gca,'units',ax.units);
    set(gca,'Position',[hfig.Position(3)/2-ax.w/2,hfig.Position(4)/2-ax.h/2,ax.w,ax.h]);
    set(gca,'XTick',[1:3]);
    set(gca,'YTick',[50:5:100]);
    set(gca,'box','on');
    
% PRINT figure
    print(hfig,'-depsc2',fullfile(OwnDir,FigDir,['fig1_labeling_accuracy-',featureSet,postProcessingTag,'.eps']))
    print(hfig,'-dpng',fullfile(OwnDir,FigDir,['fig1_labeling_accuracy-',featureSet,postProcessingTag,'.png']))    

    
    
% DISPLAY - Plot Precision ----------------------------------------------------------------

% CREATE figure        
    hfig = figure();
    clf();
    hold('on');
    
% APPLY figure properties
    set(hfig,'PaperPositionMode','auto');
    set(hfig,'units',fig.units);
    set(hfig,'Position',[0,0,fig.w,fig.h]);
    
% FORMAT statistical data for boxplot
    dat = [cell2mat(cf(@(x)x.precision,stsRaw.labelingStats));...
           cell2mat(af(@(x)x.precision,stsOpt.labelingStats))]';
    grp = bsxfun(@plus,repmat(repmat([1:numel(states)]',[numel(stsRaw.labelingStats),1]),...
                              [1,size(dat,2)]),[0:1].*numel(states));
    
% PLOT boxplot
    boxplot(dat(:),grp(:));
    h = findobj(gca,'Tag','Box');
    
% OVERLAY each boxplot with patch
    pcolor = [repmat([1,0,0],[numel(states),1]);repmat([0,0,1],[numel(states),1])];
    for j=1:length(h)

        p = patch(get(h(j),'XData'),...
                  get(h(j),'YData'),...
                  pcolor(j,:),...
                  'FaceAlpha',.5);

% REORDER patches behind boxplots
        uistack(p,'bottom');
    end

% FORMAT figure
    title('precision')
    ylim([10,100]);
    grid('on');
    set(gca,'gridlinestyle',':');
    set(gca,'units',ax.units);
    set(gca,'Position',[hfig.Position(3)/2-ax.w/2,hfig.Position(4)/2-ax.h/2,ax.w,ax.h]);
    set(gca,'XTick',[1:numel(states)*2]);
    set(gca,'XTickLabel',[states(:),states(:)])
    set(gca,'YTick',[10:10:100]);
    set(gca,'box','on');
    
% PRINT figure
    print(hfig,'-depsc2',fullfile(OwnDir,FigDir,['fig1_labeling_precision-',featureSet,postProcessingTag,'.eps']))
    print(hfig,'-dpng',fullfile(OwnDir,FigDir,['fig1_labeling_precision-',featureSet,postProcessingTag,'.png']))    
    
% DISPLAY - Plot Precision ----------------------------------------------------------------

% CREATE figure        
    hfig = figure();
    clf();
    hold('on');
    
% APPLY figure properties
    set(hfig,'PaperPositionMode','auto');
    set(hfig,'units',fig.units);
    set(hfig,'Position',[0,0,fig.w,fig.h]);
    
% FORMAT statistical data for boxplot
    dat = [cell2mat(cf(@(x)x.sensitivity,stsRaw.labelingStats));...
           cell2mat(af(@(x)x.sensitivity,stsOpt.labelingStats))]';
    grp = bsxfun(@plus,repmat(repmat([1:numel(states)]',[numel(stsRaw.labelingStats),1]),...
                              [1,size(dat,2)]),[0:1].*numel(states));
    
% PLOT boxplot
    boxplot(dat(:),grp(:));
    h = findobj(gca,'Tag','Box');
    
% OVERLAY each boxplot with patch
    pcolor = [repmat([1,0,0],[numel(states),1]);repmat([0,0,1],[numel(states),1])];
    for j=1:length(h)

        p = patch(get(h(j),'XData'),...
                  get(h(j),'YData'),...
                  pcolor(j,:),...
                  'FaceAlpha',.5);

% REORDER patches behind boxplots
        uistack(p,'bottom');
    end

% FORMAT figure
    title('sensitivity')
    ylim([10,100]);
    grid('on');
    set(gca,'gridlinestyle',':');
    set(gca,'units',ax.units);
    set(gca,'Position',[hfig.Position(3)/2-ax.w/2,hfig.Position(4)/2-ax.h/2,ax.w,ax.h]);
    set(gca,'XTick',[1:numel(states)*2]);
    set(gca,'XTickLabel',[states(:),states(:)])
    set(gca,'YTick',[10:10:100]);
    set(gca,'box','on');
    
% PRINT figure
    print(hfig,'-depsc2',fullfile(OwnDir,FigDir,['fig1_labeling_sensitivity-',featureSet,postProcessingTag,'.eps']))
    print(hfig,'-dpng',fullfile(OwnDir,FigDir,['fig1_labeling_sensitivity-',featureSet,postProcessingTag,'.png']))
    
end




