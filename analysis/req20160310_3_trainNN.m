function req20160310_3_trainNN(Trial,s)

Trial = MTATrial.validate(Trial);

load(fullfile(Trial.spath,'req20160310_1_preproc.mat'),...
     'states','fetInds','stateOrd','afet','nNeurons','nIter','rndMethod');


oind = [repmat([1:59],1,2)',zeros([118,1])];
aind = oind(:,1);
for sh = 1:117,
    oind = [oind;[circshift(aind,-sh),aind]];
end
slind = oind(fetInds{s},:);
ofet =reshape(slind,[],1);
best_inds = histc(ofet,1:59);
[~,sbind] = sort(best_inds,'descend');


gStates = states(cellfun(@isempty,...
                         regexp(states,...
                                ['(^',strjoin(stateOrd(1:s-1),'$)|(^'),'$)'])...
                         )...
                 );

pobj = parpool('local');

parfor f = 1:round(numel(sbind)/2),%numel(sbind),    
    mpn = struct;
    sub_fet = afet.copy;
    sub_fet.data = afet(:,sbind(1:f));
    sub_fet.label = [afet.label '-req20160310-' stateOrd{s} '-' num2str(f)];
    sub_fet.key = 'x';
    sub_fet.updateFilename(Trial);
    
    [mpn.Stc,mpn.d_state,mpn.labelingStats, ...
     mpn.labelingStatsMulti,mpn.model,mpn.p_state] = ...
        bhv_nn_multi_patternnet(Trial,gStates,Trial.stc,sub_fet,...
                                [],[],nNeurons,nIter, ...
                                'WSBNT','targetState',stateOrd{s});
end

delete(pobj)

save(fullfile(Trial.spath,['req20160310_3_trainNN',num2str(s),'.mat']),'gStates','-v7.3');

