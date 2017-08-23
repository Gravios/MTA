


sessionList = 'MjgER2016';
Trials = af(@(t) MTATrial.validate(t), get_session_list(sessionList));
cf(@(t) t.load('stc','msnn_ppsvd'), Trials)

% CHECK that all sessions have theta periods labeled
tper = cf(@(t) t.stc{'t'}, Trials);



for Trial = Trials,
    Trial = Trial{1};
    xyz = Trial.load('xyz');
    hfig = figure;
    sp = gobjects(0);    
    sp(end+1) = subplot2(10,10,1:2,1:10);ylim([0,325]);
    pZ(Trial);
    sp(end+1) = subplot2(10,10,3:4,1:10);        
    PlotSessionErrors(Trial,hfig);    
    sp(end+1) = subplot2(10,10,5:6,1:10);
    pXYV(Trial,'w');
    linkaxes(sp,'x')
    subplot2(10,10,7:9,1:3);        
    pXY(Trial);
    subplot2(10,10,7:9,5:7);
    histogram(xyz(Trial.stc{w+n},'head_front',3),)
    
end



