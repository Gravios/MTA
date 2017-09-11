


sessionList = 'MjgER2016';
Trials = af(@(t) MTATrial.validate(t), get_session_list(sessionList));
cf(@(t) t.load('stc','msnn_ppsvd'), Trials)


slist = get_session_list(sessionList);

% CHECK that all sessions have theta periods labeled
tper = cf(@(t) t.stc{'t'}, Trials);






for Trial = Trials,
    Trial = Trial{1};
% $$$     for t = 1:numel(slist)
% $$$     Trial = MTATrial.validate(slist(t));
% $$$     Trial = labelTheta(Trial);
% $$$ end
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
e =    histogram(xyz([Trial.stc{'w+n'}],'head_front',3),linspace(10,150,100),'Orientation','horizontal')
    
% $$$ cd(Trial.spath)
% $$$ 
% $$$         CheckEegStates(Trial.name,'theta',[],[],[],[],'display',0);
% $$$         CheckEegStates(Trial.name,[],[],[],20,[],'display',1);
end






for t = 1:numel(slist)
    Trial = MTATrial.validate(slist(t));
    s = MTASession(Trial.name,Trial.maze.name);
    s.spk.create(s);
    s.save();
end    


Trials = af(@(t) MTATrial.validate(t), get_session_list(sessionList));
cf(@(t) t.load('nq'),Trials);

e = cf(@(t) size(t.nq.eDist,1)==size(t.spk.map,1), Trials);



for t = 1:numel(slist)
    Trial = MTATrial.validate(slist(t));
    try
        pft = pfs_2d_theta(Trial,[],[],true);    
    catch err
        disp(err)
    end
end


pft = cf(@(t) pfs_2d_theta(t), Trials);

t = 10;
figure(),
for u=1:size(Trials{t}.nq.eDist,1);
    clf(),plot(pft{t},u); colorbar();title(num2str(u));waitforbuttonpress();
end
close('all');