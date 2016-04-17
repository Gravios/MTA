function stats = req20160310_4_accumStats(Trial,s)

Trial = MTATrial.validate(Trial);

load(fullfile(Trial.spath,'req20160310_1_preproc.mat'),...
     'states','fetInds','stateOrd','afet','nNeurons','nIter','rndMethod');


accum_acc = zeros([round(afet.size(2)/2),1]);
accum_pre = zeros([round(afet.size(2)/2),1]);
accum_sen = zeros([round(afet.size(2)/2),1]);


oind = [repmat([1:59],1,2)',zeros([118,1])];
aind = oind(:,1);
for sh = 1:117,
    oind = [oind;[circshift(aind,-sh),aind]];
end
slind = oind(fetInds{s},:);
ofet =reshape(slind,[],1);
best_inds = histc(ofet,1:59);
[~,sbind] = sort(best_inds,'descend');

pobj = parpool('local');

parfor f = 1:round(numel(sbind)/2),    
    opn = struct;
    sub_fet = afet.copy;
    sub_fet.data = afet(:,sbind(1:f));
    sub_fet.label = [afet.label '-req20160310-' stateOrd{1} '-' num2str(f)];
    sub_fet.key = 'x';
    sub_fet.updateFilename(Trial);
    
    model = ['MTAC_BATCH-' stateOrd{s} sub_fet.label ...
             '_SR_'  num2str(sub_fet.sampleRate) ...
             '_NORM_' num2str(0) ...             
             '_REF_' Trial.filebase ...
             '_STC_' Trial.stc.mode ...
             '_NN_'  num2str(nNeurons)...
             '_NI_'  num2str(nIter)...
             '_'     'NN_multiPN'...
             '_'     'RAND_' rndMethod];

    [opn.Stc,opn.d_state,opn.labelingStats, ...
     opn.labelingStatsMulti,opn.model,opn.p_state] = ...
        bhv_nn_multi_patternnet(Trial,states,Trial.stc,sub_fet,...
                                [], model,nNeurons,nIter,...
                                rndMethod,'targetState',stateOrd{s});

    accum_acc(f) = opn.labelingStats.accuracy;
    accum_pre(f) = opn.labelingStats.precision(1);
    accum_sen(f) = opn.labelingStats.sensitivity(2);
end

delete(pobj)

save(fullfile(Trial.spath,'req20160310_4_accumStats',num2str(s),'.mat'),'accum_acc','accum_pre','accum_sen','-v7.3');