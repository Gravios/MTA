function stats = req20160310_4_accumStats(Trial,s)

Trial = MTATrial.validate(Trial);
RefTrial = Trial;
%RefTrial = MTATrial.validate('jg05-20120317.cof.all');
%RefTrial.load('stc','hand_labeled_rev3_jg');


%'states','fetInds','stateOrd','afet','nNeurons','nIter','rndMethod'
ds = load(fullfile(Trial.spath,'req20160310_1_preproc-afet.mat'));


accum_acc = zeros([round(ds.afet.size(2)/2),1]);
accum_pre = zeros([round(ds.afet.size(2)/2),1]);
accum_sen = zeros([round(ds.afet.size(2)/2),1]);


oind = [repmat([1:59],1,2)',zeros([118,1])];
aind = oind(:,1);
for sh = 1:117,
    oind = [oind;[circshift(aind,-sh),aind]];
end
slind = oind(ds.fetInds{s},:);
ofet =reshape(slind,[],1);
best_inds = histc(ofet,1:59);
[~,sbind] = sort(best_inds,'descend');

gStates = ds.states(cellfun(@isempty,...
                         regexp(ds.states,...
                                ['(^',strjoin(ds.stateOrd(1:s-1),'$)|(^'),'$)'])...
                         )...
                 );


pobj = parpool(10);

%for f = 1:round(numel(sbind)/2),    
parfor f = 1:round(numel(sbind)/2),    
    opn = struct;
    sub_fet = ds.afet.copy;
    sub_fet.data = ds.afet(:,sbind(1:f));
    sub_fet.label = [ds.afet.label '-req20160310-' ds.stateOrd{s} '-' num2str(f)];
    sub_fet.key = 'x';
    sub_fet.updateFilename(Trial);
    
    model = ['MTAC_BATCH-' ds.stateOrd{s} sub_fet.label ...
             '_SR_'  num2str(sub_fet.sampleRate) ...
             '_NORM_' num2str(0) ...             
             '_REF_' RefTrial.filebase ...
             '_STC_' RefTrial.stc.mode ...
             '_NN_'  num2str(ds.nNeurons)...
             '_NI_'  num2str(ds.nIter)...
             '_'     'NN_multiPN'...
             '_'     'RAND_' ds.rndMethod];

    [opn.Stc,opn.d_state,opn.labelingStats, ...
     opn.labelingStatsMulti,opn.model,opn.p_state] = ...
        bhv_nn_multi_patternnet(Trial,gStates,Trial.stc,sub_fet,...
                                [], model,ds.nNeurons,ds.nIter,...
                                ds.rndMethod,'targetState',ds.stateOrd{s});

    accum_acc(f) = opn.labelingStats.accuracy;
    accum_pre(f) = opn.labelingStats.precision(1);
    accum_sen(f) = opn.labelingStats.sensitivity(2);
end

delete(pobj);

save(fullfile(Trial.spath,['req20160310_4_accumStats',num2str(s),'.mat']),'sbind','accum_acc','accum_pre','accum_sen','-v7.3');