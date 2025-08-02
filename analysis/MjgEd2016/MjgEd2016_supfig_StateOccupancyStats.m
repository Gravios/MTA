slist = get_session_list('hand_labeled');

states = {'walk','rear','turn','pause','groom','sit'};
socc = zeros([numel(slist),numel(states)]);

si = 1;
for s = slist
    Trial = MTATrial(s.sessionName,s.mazeName,s.trialName);
    stc = Trial.load('stc',s.stcMode);
    sti = 1;
    for state = [stc{states{:}}]
        state = state{1};
        socc(si,sti) = sum(diff(state.data,1,2)/state.sampleRate);
        sti = sti+1;
    end    
    si = si+1;    
end


set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)

hfig = figure(2016070501);
boxplot(socc/60)

title('State Occupancy of Manual Annotations');
ylabel('Time (minutes)')
set(gca,'xticklabels',states)

print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_2',...
                     'F_sup_state_occupancy.eps'))
