%% SupFig State Transition Statistics
%

% DEF figure variables -----------------------------------------------------------------
% figure save paths
OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'Shared/Behavior Paper/Figures/Suplementary/state_transition_stats';
% END figure variables -----------------------------------------------------------------


param = struct('sessionList',            'hand_labeled',...
               'referenceTrial',         'jg05-20120317.cof.all',...
);


sessionList     = get_session_list(param.sessionList);                    
numSessions     = numel(sessionList);                    
sampleRate = repmat({param.sampleRate},1,numSessions);
embeddingWindow = repmat({64},1,numSessions);

                    
% LOAD Trial objects
Trials = af(@(Trial) MTATrial.validate(Trial)  , sessionList);

% LOAD State Collections
Stc    = cf(@(Trial)  Trial.load('stc')         , Trials);
StcNN  = cf(@(Trial) Trial.load('stc','NN0317'), Trials);

states = {'walk','rear','turn','pause','groom','sit'};
nsts = numel(states);






% SUPFIG STSTRANSPROB ----------------------------------------------------------------------------
% project: MjgEd2016
% title: State Transition probabilities
% parent: figure 3
% subplots:
%    subplot 1: State Transition probabilities
%    subplot 2: probability of subsequent state given antecedent state
%    subplot 3: probability of antecedent state given subsequent state
% location: MjgEd2016_figure3_final.m


stpa = cell([nsts,nsts]);
for t = 1:nsts,
    for o = 1:nsts,
        if t~=o,
            % SELECT state to state transitions
            stpa{t,o} = cell2mat(cf(@(stc,Trial,state,x) ...
                                    stc.get_state_transitions(Trial,{states{t},states{o}},[],x),...
                                    Stc,Trials,repmat(states(t),1,numSessions),xyz)');
        end
    end
end
stf = cell2mat(cf(@length,stpa));

figure
subplot(221);
imagesc(round(stf./sum(stf(:)),4));
title('P(antecedent,subsequent)')
subplot(222);
imagesc(round(bsxfun(@rdivide,stf,sum(stf)),2));
caxis([0,1])
title('P(antecedent|subsequent)')
subplot(223);
imagesc(round(bsxfun(@rdivide,stf,sum(stf,2)),2))
caxis([0,1])
title('P(subsequent|antecedent)')

ForAllSubplots(['xlabel(''Subsequent State'');'...
                'ylabel(''Antecedent State'');'...
                'colorbar;'...
                'set(gca,''XTick'',[1:6]);'...
                'set(gca,''YTick'',[1:6]);'...                
                'set(gca,''XTickLabel'',{''' strjoin(states,''',''') '''});'...
                'set(gca,''YTickLabel'',{''' strjoin(states,''',''') '''});'...
]);
colormap('jet');

suptitle('State Transition Statistics')

FigName = ['stateTransitionMatrix_hand_labeled'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']),);
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));                


% END SUPFIG STSTRANSPROB ----------------------------------------------------------------------------
     


% SUPFIG sorted transition probs----------------------------------------------------------------------


figure,
plot(sort(log10(astf(:))));
xlabel('sorted transition probability');
ylabel('log10 probability');
FigName = ['stateTransition_sorted_log_prob'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END SUPFIG sorted transition probs------------------------------------------------------------------


% SUPFIG state occupancy -----------------------------------------------------------------------------


% END SUPFIG state occupancy -----------------------------------------------------------------------------
