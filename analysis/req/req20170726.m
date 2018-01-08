% req20170726 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: ctx units
%  Bugs: NA
 

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




[Stc,labelingStats,networkOutput] =                                              ...
    bhv_nn_multi_session_patternnet(                                             ...
        sessionList,            featureSet,   states,   keys,  model,            ...
        sampleRate,             nNeurons,     nIter,    randomizationMethod,     ...
        map2reference,          normalize,    referenceTrial,                    ...
        trainingSessionList,    normalizationSessionList,                        ...
        dropIndex,              prctTrain                                        ...
        );


Trials = af(@(Trial) MTATrial.validate(Trial), get_session_list(sessionList));
cf(@(t) t.load('stc','msnnN0+hand_labeled'), Trials);
Stc = cf(@(t) t.stc.copy(), Trials);

state = 'walk';
durThresh  = 1;
BccgW = cf(@(Trial,state,durThresh) gen_bhv_ccg(Trial,state,durThresh),...
          Trials,...
          repmat({state},    [1,numel(Trials)]),...
          repmat({durThresh},[1,numel(Trials)]));

state = 'turn';
durThresh  = 0.5;
BccgN = cf(@(Trial,state,durThresh) gen_bhv_ccg(Trial,state,durThresh),...
          Trials,...
          repmat({state},    [1,numel(Trials)]),...
          repmat({durThresh},[1,numel(Trials)]));

state = 'rear';
durThresh  = 1;
BccgR = cf(@(Trial,state,durThresh) gen_bhv_ccg(Trial,state,durThresh),...
          Trials,...
          repmat({state},    [1,numel(Trials)]),...
          repmat({durThresh},[1,numel(Trials)]));

figure,
s = 1;
for u = 1:size(BccgW{s}.ccg,2)
    cla;
    BccgW{s}.plot(u,1);
    title(num2str(u));
    %pause(.4);
    waitforbuttonpress;
end



%nq = cf(@(t) NeuronQuality(t,[],[],[],false),  Trials);
cf(@(t) t.load('nq'),  Trials);
nq = cf(@(t) t.nq,     Trials);
pfs = cf(@(t,s) pfs_2d_states(t,s,{'rear','turn','walk'}), Trials, Stc);
[accg,tbins] = cf(@(t) autoccg(MTASession.validate(t.filebase)),  Trials);

% figure save paths
OwnDir = '/storage/gravio/nextcloud/';
FigDir = 'MjgER2016/figures/figure3/parts/state_transition_ctx';

units = {};
units{1} = [3,4,5,11,22,23,28,38,42,43,46,54,57];
units{2} = [9,19,21,26,29,58];
units{3} = [9,12,17,21,34,43,51,52,60,61];
units{4} = [2,3,11,13,14,16,19,21,28,29,40,49,52,64,65,67,72,77,87];

pfs = cf(@(t,u) pfs_2d_states(t,u','msnnN0+hand_labeled',{'rear','turn','walk'}),Trials,units);


s = 4;
u = units{s}(6);

figure();
set(gcf,'PaperPositionMode','auto');
nx = 4;
ny = 4;


for s = 1:4,
for u = units{s}

    clf();

    %maxPfsRate = max(cell2mat(cf(@(p,u) p.maxRate(u), pfs{s},repmat({u},1,numel(pfs{s})))));
    maxCcgRate = max(reshape([sq(BccgW{s}.ccg(:,u,:)),...
                              sq(BccgN{s}.ccg(:,u,:)),...
                              sq(BccgR{s}.ccg(:,u,:))],[],1));
    
    
    % accg
    subplot2(ny,nx,2,1);    
    %bar(tbins{s},accg{s}(:,u));axis tight;            
    title({Trials{s}.filebase,['autoccg unit: ',num2str(u)]});    
    
    % THETA 
    subplot2(ny,nx,1,1);
    %pfs{s}{1}.plot(u,[],[],maxPfsRate);
    title(['Theta']);
    
    % REAR
    subplot2(ny,nx,1,2);    
    %pfs{s}{2}.plot(u,[],[],maxPfsRate);
    title(['Rear']);
    subplot2(ny,nx,2,2);
    BccgR{s}.plot(u,1);xlim([-2,2]);ylim([0,maxCcgRate]);
    title(['unit: ',num2str(u),' Rear Onset']);
    subplot2(ny,nx,3,2);
    BccgR{s}.plot(u,2);xlim([-2,2]);ylim([0,maxCcgRate]);
    title(['unit: ',num2str(u),' Rear Offset']);
    %ForAllSubplots('xlim([-2,2])');

    % TURN
    subplot2(ny,nx,1,3);    
    %pfs{s}{3}.plot(u,[],[],maxPfsRate);
    title(['Turn']);    
    subplot2(ny,nx,2,3);
    BccgN{s}.plot(u,1);xlim([-2,2]);ylim([0,maxCcgRate]);
    title(['unit: ',num2str(u),' Turn Onset']);
    subplot2(ny,nx,3,3);
    BccgN{s}.plot(u,2);xlim([-2,2]);ylim([0,maxCcgRate]);
    title(['unit: ',num2str(u),' Turn Offset']);
    %ForAllSubplots('xlim([-2,2])');
    
    
    % Walk
    subplot2(ny,nx,1,4);
    %pfs{s}{4}.plot(u,[],[],maxPfsRate);
    title(['Walk']);        
    subplot2(ny,nx,2,4);
    BccgW{s}.plot(u,1);xlim([-2,2]);ylim([0,maxCcgRate]);
    title(['unit: ',num2str(u),' Walk Onset']);
    subplot2(ny,nx,3,4);
    BccgW{s}.plot(u,2);xlim([-2,2]);ylim([0,maxCcgRate]);
    title(['unit: ',num2str(u),' Walk Offset']);


    pause(0.5)    

    
    FigName = ['stateTransition_unit_CCG_',Trials{s}.filebase,'_unit_',num2str(u)];
    print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
    
end

end

