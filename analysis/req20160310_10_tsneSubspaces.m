function req20160310_10_tsneSubspaces(Trial,s)

%RefTrial = 'jg05-20120317.cof.all';

Trial = MTATrial.validate(Trial);
RefTrial = Trial;
%RefTrial = MTATrial.validate(RefTrial);

load(fullfile(Trial.spath,'req20160310_1_preproc-afet.mat'));
bs = load(fullfile(RefTrial.spath,'req20160310_5_genfigs.mat'));
ts = load(fullfile(RefTrial.spath,'req20160310_8_genOptfigs.mat'));


% State list which excludes previous states
gStates = states(cellfun(@isempty,...
                         regexp(states,...
                                ['(^',strjoin(stateOrd(1:s-1),'$)|(^'),'$)'])...
                         )...
                 );
% List of states {composite, target}
tstates = {[strjoin({gStates{find(cellfun(@isempty,regexp(gStates,['(',strjoin(stateOrd(1:s),')|('),')'])))}},'+'),'&gper'],stateOrd{s}};

sfet = afet.copy;
sfet.data = afet(:,ts.bfets{s});
out = mta_tsne(Trial,sfet,12,Trial.stc.copy,tstates,3,2,80,'ifReportFig',true,'overwrite',true);    


save(fullfile(Trial.spath,['req20160310_10_tsneSubspaces',num2str(s),'.mat']),'out','-v7.3');


