% req20170726 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: ctx units
%  Bugs: NA
 



%% SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sessionList = get_session_list('ctx');
for i = 1:3%numel(sessionList);
s = MTASession.validate(sessionList(i));
s.spk.create(s);
s.save
%headMarkers = regexpi(s.xyz.model.ml,'^head_[A-Za-z]*','match');
%headMarkers = cellfun(@(x) cellstr(x),headMarkers(~cellfun('isempty',headMarkers)));
%ERCOR_fillgaps_RidgidBody(s,'MANUAL','BEST_SWAP_PERMUTATION',[],headMarkers);
end





QuickTrialSetup(s,'all','cof',[15,-15]);

Trial = MTATrial.validate(sl(i));
figure,
subplot2(2,2,1,1),
pXY(Trial),
subplot2(2,2,1,2),
pZ(Trial)
subplot2(2,2,2,[1,2]),
PlotSessionErrors(Trial,gcf)

Trials = af(@(s) MTATrial.validate(s),sessionList);



% SET neural network labeling args
varargin = {};
defargs = struct('sessionList',                 'ctx',                                         ...
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
                 'dropIndex',                   2,                                             ...
                 'prctTrain',                   90,                                            ...
                 'postProcessingTag',           ''                                             ...
);
[sessionList,featureSet,states,keys,model,sampleRate,nNeurons,nIter,randomizationMethod,       ...
 map2reference,normalize,referenceTrial,trainingSessionList,normalizationSessionList,          ...
 dropIndex,prctTrain,postProcessingTag] = DefaultArgs(varargin,defargs,'--struct');


% COMPOSE neural network model name
dropIndexTag = '';
if ~isempty(dropIndex), dropIndexTag = num2str(dropIndex);
    dropIndexTag(dropIndexTag==[' '])='_';
end
model = ['MTAC_BATCH+' trainingSessionList '+'            ...
         dropIndexTag '+'                                 ...
         featureSet                                       ...
         '+SR'   num2str(sampleRate)                      ...
         'NN'    num2str(nNeurons)                        ...
         'NI'    num2str(nIter)                           ...
         'M'     num2str(map2reference)                   ...
         'MREF+' referenceTrial                           ...                          
         '+N'    num2str(normalize)                       ...             
         'NREF+' normalizationSessionList                 ...
         '+RND'  randomizationMethod                      ...
         '+PRCT' num2str(prctTrain)                       ...
         '+STS+' strjoin(keys,'')                         ...
         '+'     'multiSesPatNet'];

% LABEL Behavior
[Stc,labelingStats,networkOutput] =                                              ...
    bhv_nn_multi_session_patternnet(                                             ...
        sessionList,            featureSet,   states,   keys,  model,            ...
        sampleRate,             nNeurons,     nIter,    randomizationMethod,     ...
        map2reference,          normalize,    referenceTrial,                    ...
        trainingSessionList,    normalizationSessionList,                        ...
        dropIndex,              prctTrain                                        ...
        );




%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sessionList = get_session_list('ctx');

% LOAD Trials
Trials = af(@(T)  MTATrial.validate(T),                 get_session_list(sessionList));
         cf(@(T)  T.load('nq'),                         Trials);
         cf(@(T)  T.load('stc','msnnN0+hand_labeled'),  Trials);
Stc    = cf(@(T)  T.stc.copy(),                         Trials);



states     = {'walk','turn','rear','pause','groom','sit','theta','gper-theta'};
durThresh  = {    1 ,  0.5 ,    1 ,     1 ,     1 ,   1 ,     1 ,          1 };


% COMPUTE unit x state transitions ccg
Bccg = cell([1,numel(states)]);
for s = 1:numel(states),
    Bccg{s} = cf(@(Trial,state,durThresh) gen_bhv_ccg(Trial,state,durThresh,[],1),...
                 Trials,...
                 repmat({states{s}},    [1,numel(Trials)]),...
                 repmat({durThresh{s}},[1,numel(Trials)]));
end

pfs = cf(@(T,s) pfs_2d_states(T,'all','msnnN0+hand_labeled',s), ...
         Trials, ...
         repmat({states},size(Trials)));
[accg,tbins] = cf(@(T) autoccg(MTASession.validate(T.filebase)),  Trials);

% figure save paths
OwnDir = '/storage/gravio/figures/';
FigDir = 'analysis/ctx_bhv_state_transition';
create_directory(fullfile(OwnDir,FigDir));

units = {};
units{1} = [3,4,5,11,22,23,28,38,42,43,46,54,57];
units{2} = [9,19,21,26,29,58];
units{3} = [9,12,17,21,34,43,51,52,60,61];
units{4} = [2,3,11,13,14,16,19,21,28,29,40,49,52,64,65,67,72,77,87];

pfs = cf(@(t,u) pfs_2d_states(t,u','msnnN0+hand_labeled',{'rear','turn','walk'}),Trials,units);

pfs = cf(@(ps,pt) cat(2, ps,{pt}), pfs,pft);

s = 4;
u = units{s}(6);

hfig = figure();
set(hfig,'PaperPositionMode','auto');
nx = 4;
ny = 4;


for tind = 1:numel(Trials),
    for u = Bccg{tind}{1}.cluMap(:)'
        clf();
        %maxPfsRate = max(cell2mat(cf(@(p,u) p.maxRate(u), pfs{s},repmat({u},1,numel(pfs{s})))));
        maxCcgRate = max(reshape([sq(Bccg{1}{tind}.ccg(:,u,:)),...
                                  sq(Bccg{2}{tind}.ccg(:,u,:)),...
                                  sq(Bccg{3}{tind}.ccg(:,u,:)),...
                                  sq(Bccg{4}{tind}.ccg(:,u,:))],[],1));
        if maxCcgRate==0, continue,end
        % ACCG
        subplot2(ny,nx,1,1);
        bar(tbins{tind},accg{tind}(:,u));axis tight;
        title({Trials{tind}.filebase,['autoccg unit: ',num2str(u)]});    

        % PFD-PP
        subplot2(ny,nx,1,2);
        plot(pfd{tind},u,'mean',false,maxCcgRate,false,0.9,true);
        
        % THETA PFS
        for s = 1:numel(states),
            subplot2(ny,nx,2,s);
            plot(pfs{tind}{s},u,'mean',false,maxCcgRate,false,0.9,false);
            title(pfs{tind}{s}.parameters.states);
        end


        for s = 1:numel(states),
            subplot2(ny,nx,3,s);
            Bccg{s}{tind}.plot(u,1);xlim([-2,2]);ylim([0,maxCcgRate]); % PLOT state onset ccg
            title(['unit: ',num2str(u),' ',states{s},' Onset']);
            subplot2(ny,nx,4,s);
            Bccg{s}{tind}.plot(u,2);xlim([-2,2]);ylim([0,maxCcgRate]); % PLOT state offset ccg
            title(['unit: ',num2str(u),' ',states{s},' Offset']);
        end
        
        
        drawnow();    
        pause(0.5);

        FigName = ['stateTransition_unit_CCG_',Trials{s}.filebase,'_unit_',num2str(u)];
        print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
        print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
    end
end




%pitchReferenceTrial = 'Ed05-20140529.ont.all';


xyz  = cf(@(T) preproc_xyz(T),             Trials);
pch  = cf(@(T) fet_HB_pitchB(T),           Trials);
tper = cf(@(T) [T.stc{'theta-groom-sit'}], Trials);
       cf(@(p) p{1}.resample(xyz), tper);

overwrite = true;

for tind = 1:numel(Trials),
pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units   = [];%units{tind};
pargs.numIter = 101;
pargs.halfsample = true;
pargs.tag            = 'DRZxHBPITCHxBPITCH_v2';
pargs.boundaryLimits = [-pi/2,pi/2;-2,pi/2];
pargs.binDims        = [0.1,0.1];
pargs.SmoothingWeights = [2,2];
if overwrite,  
    pargs.overwrite = true;    
    pargs.xyzp = MTADxyz('data',pch{tind}.data,'sampleRate',xyz{tind}.sampleRate);
    pargs.states = tper{tind};
    pfsArgs = struct2varargin(pargs);
    MTAApfs(Trials{tind},pfsArgs{:});

end
pargs.states    = tper{tind};
pargs.units     = [];%units{tind};
pargs.overwrite = false;
pfsArgs = struct2varargin(pargs);
pfd{tind} = MTAApfs(Trials{tind},pfsArgs{:});

end


figure,for u = pfd{3}.data.clu,
    clf();    
    subplot(241);
    % place fields

    plot(pfd{3},u,'mean',true,[],false,0.9,true);waitforbuttonpress();end
end