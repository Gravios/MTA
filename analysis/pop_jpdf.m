
Sessions = get_session_list('test_grp',...
                '/storage/gravio/data/processed/xyz/',...
                '/storage/gravio/data/processed/nlx/');




%% Speed
pvel = [];

for s = Sessions
    Trial = MTATrial(s.name,s.trialName,s.mazeName);    
    xyz = Trial.load('xyz');
    vl = vel(xyz,1:8,[1,2]);
    vl.data = ButFilter(vl.data,3,4/(vl.sampleRate/2),'low');
    ind = nniz(vl);
    pvel = cat(1,pvel,[median(vl(ind,1:2),2),median(vl(ind,5:8),2)]);
end


figure,
hist2(log10(pvel),linspace(-.5,2,75),linspace(-.5,2,75));
xlabel('log10 body speed (cm/s)');
ylabel('log10 head speed (cm/s)');
title('JPDF of log10 head and body speeds');




%% Speed

Sessions = get_session_list('jg05',...
                '/storage/gravio/data/processed/xyz/',...
                '/storage/gravio/data/processed/nlx/');


pvel = [];
for s = Sessions
    Trial = MTATrial(s.name,s.trialName,s.mazeName);    
    xyz = Trial.load('xyz');
    ang = create(MTADang,Trial,xyz);
    vl = vel(xyz,1,[1,2]);
    vl.data = ButFilter(vl.data,3,4/(vl.sampleRate/2),'low');
    ind = nniz(vl);
    %ind = vl.data>3;
    pvel = cat(1,pvel,[ang(ind,1,7,3),ang(ind,1,7,2)]);
end


figure,
hist2(pvel,linspace(100,280,75),linspace(-.25,1.5,75));
xlabel('Distance tail to head');
ylabel('Pitch tail to head');
title('JPDF of body length and pitch');



