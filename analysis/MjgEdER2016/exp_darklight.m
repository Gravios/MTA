

overWriteSessions = false;
overWriteTrials   = false;
overWriteStc      = false;


sessionList = 'Ed10_cof';
  trialList = 'Ed10_DarkLight';


  
Slist = SessionList(sessionList);
if overWriteSessions, 
    QuickSessionSetup(Slist,[],[],false); 
    QuickTrialSetup(Slist);
end



Tlist = SessionList(trialList,...
                '/storage/gravio/data/processed/xyz/Ed10/',...
                '/storage/eduardo/data/processed/nlx/Ed10/');
if overWriteTrials,   
    QuickTrialSetup(Tlist,'overwrite',true); 
end


s = MTASession('Ed10-20140817');

xyz = s.load('xyz');
pXY(s)
pZ(s)
PlotSessionErrors(s)



Trial = MTATrial.validate('Ed10-20140814.cof.all');

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
vxy = fxyz.vel([1,7],[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);


ang = create(MTADang,Trial,xyz);



ind = vxy(:,1)>0.5;
figure,hist(ang(ind,'head_back','head_front',2),100)

figure,hist(xyz(ind,'head_back',3),100)

