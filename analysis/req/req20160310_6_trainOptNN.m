function req20160310_6_trainOptNN(Trial,s)

Trial = MTATrial.validate(Trial);

%'states','fetInds','stateOrd','afet','nNeurons','nIter','rndMethod'
ds = load(fullfile(Trial.spath,'req20160310_1_preproc-afet.mat'));
bs = load(fullfile(Trial.spath,'req20160310_5_genfigs.mat'));


gStates = ds.states(cellfun(@isempty,...
                         regexp(ds.states,...
                                ['(^',strjoin(ds.stateOrd(1:s-1),'$)|(^'),'$)'])...
                         )...
                 );

pobj = parpool(12);

parfor f = 1:numel(bs.bFetInds{s}),
    mpn = struct;
    sub_fet = ds.afet.copy;
    sub_fet.data = ds.afet(:,bs.bFetInds{s}(1:f));
    sub_fet.label = [ds.afet.label '-req20160310-optfet' ds.stateOrd{s} '-' num2str(f)];
    sub_fet.key = 'x';
    sub_fet.updateFilename(Trial);
    
    [mpn.Stc,mpn.d_state,mpn.labelingStats, ...
     mpn.labelingStatsMulti,mpn.model,mpn.p_state] = ...
        bhv_nn_multi_patternnet(Trial,gStates,Trial.stc,sub_fet,...
                                [],[],ds.nNeurons,ds.nIter, ...
                                ds.rndMethod,'targetState',ds.stateOrd{s});
end

delete(pobj)

save(fullfile(Trial.spath,[mfilename,num2str(s),'.mat']),'gStates','-v7.3');

