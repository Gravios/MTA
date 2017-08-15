% req20170726 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Labeling statistics 
%  Bugs: NA
 

sl = get_session_list('ctx');
for i = 1:3%numel(sessionList);
s = MTASession.validate(sl(i));
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

state = 'walk';
durThresh  = 1.5;
BccgW = cf(@(Trial,state,durThresh) gen_bhv_ccg(Trial,state,durThresh),...
          Trials,...
          repmat({state},    [1,numel(Trials)]),...
          repmat({durThresh},[1,numel(Trials)]));

state = 'rear';
durThresh  = 1.5;
BccgR = cf(@(Trial,state,durThresh) gen_bhv_ccg(Trial,state,durThresh),...
          Trials,...
          repmat({state},    [1,numel(Trials)]),...
          repmat({durThresh},[1,numel(Trials)]));

figure,
s = 1;
for u = 1:size(BccgW{s}.ccg,2)
    cla;
    Bccg{s}.plot(u,1);
    title(num2str(u));
    %pause(.4);
    waitforbuttonpress;
end





% figure save paths
OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'Shared/Behavior Paper/Figures/Figure_3/parts';

figure,
s = 1;
units = [3,4,5,11,22,23,28,38,42,43,46,54,57];

s = 2;
units = [9,19,21,26,29,58];

s = 3;
units = [9,12,17,21,34,43,51,52,60,61];

s = 4
units = [2,3,11,13,14,16,19,21,28,29,40,49,52,64,65,67,72,77,87];
for u = units

    clf;
    subplot2(2,2,1,1);
    BccgW{s}.plot(u,1);
    title(['unit: ',num2str(u),' Walk Onset']);
    subplot2(2,2,2,1);
    BccgW{s}.plot(u,2);
    title(['unit: ',num2str(u),' Walk Offset']);
    subplot2(2,2,1,2);
    BccgR{s}.plot(u,1);
    title(['unit: ',num2str(u),' Rear Onset']);
    subplot2(2,2,2,2);
    BccgR{s}.plot(u,2);
    title(['unit: ',num2str(u),' Rear Offset']);
    ForAllSubplots('xlim([-2,2])')

    fc = get(gcf,'Children');
    ForAllSubplots(['ylim([0,',num2str(max([fc(:).YLim])),'])'])

    FigName = ['stateTransition_unit_CCG_',Trials{s}.filebase,'_unit_',num2str(u)];
    print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
    
end



