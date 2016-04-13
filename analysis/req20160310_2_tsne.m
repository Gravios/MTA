function req20160310_2_tsne(Trial,s)

load(fullfile(Trial.spath,'req20160310_1_preproc.mat'),...
     'states','fetInds','stateOrd','tfet','tstc','fet');

out = struct;

% State list which excludes previous states
gStates = states(cellfun(@isempty,...
                         regexp(states,...
                                ['(^',strjoin(stateOrd(1:s-1),'$)|(^'),'$)'])...
                         )...
                 );
% List of states {composite, target}
tstates = {[strjoin({gStates{find(cellfun(@isempty,regexp(gStates,['(',strjoin(stateOrd(1:s),')|('),')'])))}},'+'),'&gper'],stateOrd{s}};

sfet = fet.copy;
sfet.data = tfet(:,fetInds{s});
out = mta_tsne(Trial,sfet,12,tstc,tstates,5,2,80,'ifReportFig',false,'overwrite',false);    


save(fullfile(Trial.spath,'req20160310_2_tsne',num2str(s),'.mat'),'out','-v7.3');

