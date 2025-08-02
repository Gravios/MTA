function req20160310_3_trainNN(Trial,s)

Trial = MTATrial.validate(Trial);

%'states','fetInds','stateOrd','afet','nNeurons','nIter','rndMethod'
ds = load(fullfile(Trial.spath,'req20160310_1_preproc-afet.mat'));

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

%pobj = parpool(4);

for f = 1:round(numel(sbind)/2),%numel(sbind),    
    mpn = struct;
    sub_fet = ds.afet.copy;
    sub_fet.data = ds.afet(:,sbind(1:f));
    sub_fet.label = [ds.afet.label '-req20160310-' ds.stateOrd{s} '-' num2str(f)];
    sub_fet.key = 'x';
    sub_fet.updateFilename(Trial);
    
    [mpn.Stc,mpn.d_state,mpn.labelingStats, ...
     mpn.labelingStatsMulti,mpn.model,mpn.p_state] = ...
        bhv_nn_multi_patternnet(Trial,gStates,Trial.stc,sub_fet,...
                                [],[],ds.nNeurons,ds.nIter, ...
                                ds.rndMethod,'targetState',ds.stateOrd{s});
end

%delete(pobj)

save(fullfile(Trial.spath,['req20160310_3_trainNN',num2str(s),'.mat']),'gStates','-v7.3');

