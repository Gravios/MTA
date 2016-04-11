function req20160310_2_tsne(Trial)

load(fullfile(Trial.spath,'req20160310_1_preproc.mat'),...
     'states','stateOrd','tfet','tstc','fet');

out = struct;
for sind = 1:numel(stateOrd),
    % State list which excludes previous states
    gStates = states(cellfun(@isempty,...
                             regexp(states,...
                                    ['(^',strjoin(stateOrd(1:sind-1),'$)|(^'),'$)'])...
                             )...
                     );
    % List of states {composite, target}
    tstates = {[strjoin({gStates{find(cellfun(@isempty,regexp(gStates,['(',strjoin(stateOrd(1:sind),')|('),')'])))}},'+'),'&gper'],stateOrd{sind}};
    
    sfet = fet.copy;
    sfet.data = tfet(:,fetInds{sind});
    out(sind) = mta_tsne(Trial,sfet,12,tstc,tstates,5,2,80,'ifReportFig',false,'overwrite',false);    
end

save(fullfile(Trial.spath,'req20160310_2_tsne'),'out','-v7.3');

