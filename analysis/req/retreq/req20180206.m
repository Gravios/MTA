

sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';


% LOAD Trials
% COMPUTE placefield statistics
Trials  = af(@(t)  MTATrial.validate(t),   sessionList([4,18]));

% $$$ 
% $$$ cf(@(t) transform_rigidBody(t,true,false), Trials);
% $$$ cf(@(t) preproc_xyz(t,'SPLINE_SPINE_HEAD_EQD',true), Trials);


numTrials = numel(Trials);

stc  = cf(@(t)   t.stc.copy(),         Trials);
xyz  = cf(@(t)   preproc_xyz(t,'trb'),  Trials);
xyzs = cf(@(t)   preproc_xyz(t,'SPLINE_SPINE_HEAD_EQD'),  Trials);
xyzi = cf(@(t)   preproc_xyz(t,'SPLINE_SPINE_HEAD_EQI'),  Trials);
ang  = cf(@(t,x) create(MTADang,t,x),  Trials,xyz);
angs = cf(@(t,x) create(MTADang,t,x),  Trials,xyzs);
angi = cf(@(t,x) create(MTADang,t,x),  Trials,xyzi);
fetm = cf(@(t)   fet_bref_rev12(t),    Trials);
fetm{1}.map_to_reference_session(Trials{1},Trials{2});

states = {'loc','rear','pause','groom','sit'};
numStates = numel(states);

% CREATE feature
fet  = cf(@(t,a) MTADfet.encapsulate(t,[a(:,'spine_middle','hcom',3)],...
                                    a.sampleRate,'tfet','tfet','t'), Trials, ang);
fets = cf(@(t,a) MTADfet.encapsulate(t,[a(:,'spine_middle','hcom',3)],...
                                    a.sampleRate,'tfet','tfet','t'), Trials, angs);
feti = cf(@(t,a) MTADfet.encapsulate(t,[a(:,'spine_middle','hcom',3)],...
                                    a.sampleRate,'tfet','tfet','t'), Trials, angi);

fet  = cf(@(t,a) MTADfet.encapsulate(t,[a(:,'spine_middle','hcom',2)],...
                                    a.sampleRate,'tfet','tfet','t'), Trials, ang);
fets = cf(@(t,a) MTADfet.encapsulate(t,[a(:,'spine_middle','hcom',2)],...
                                    a.sampleRate,'tfet','tfet','t'), Trials, angs);
feti = cf(@(t,a) MTADfet.encapsulate(t,[a(:,'spine_middle','hcom',2)],...
                                    a.sampleRate,'tfet','tfet','t'), Trials, angi);


fet  = cf(@(t,a) MTADfet.encapsulate(t,[a(:,'pelvis_root',3)],...
                                    a.sampleRate,'tfet','tfet','t'), Trials, xyz);
fets = cf(@(t,a) MTADfet.encapsulate(t,[a(:,'pelvis_root',3)],...
                                    a.sampleRate,'tfet','tfet','t'), Trials, xyzs);
feti = cf(@(t,a) MTADfet.encapsulate(t,[a(:,'pelvis_root',3)],...
                                    a.sampleRate,'tfet','tfet','t'), Trials, xyzi);

fet = cf(@(t,a,s) MTADfet.encapsulate(t,[a(:,'pelvis_root','spine_upper',2),...
                                         s(:,'pelvis_root','spine_upper',2)],...
                                    a.sampleRate,'tfet','tfet','t'), Trials, ang,angs);
edx = linspace(-pi/2,pi/2,50);
edy = linspace(-pi/2,pi/2,50);        


fet = cf(@(t,x) MTADfet.encapsulate(t,[x(:,'spine_middle' ,3),...
                                       x(:,'spine_upper',3)],...
                                    x.sampleRate,'tfet','tfet','t'), Trials, xyz);
edx = linspace(40,140,50);
edy = linspace(40,170,50);        

fet = cf(@(t,x) MTADfet.encapsulate(t,[x(:,10),...
                                       x(:,11)],...
                                    x.sampleRate,'tfet','tfet','t'), Trials, fetm);
edx = linspace(40,140,50);
edy = linspace(40,170,50);        



figure();
for t = 1:numTrials,
    for s = 1:numStates,
        subplot2(numStates,numTrials,s,t);
        ind = [stc{t}{states{s}}];
        hist2([fet{t}(ind,1),fet{t}(ind,2)],edx,edy);
        caxis([0,2400])
        grid('on');
    end
end


i = 1;
edy = linspace(10,180,150);
figure();
ind = [stc{1}{'x+p'}];
subplot(231);
hist(fet{1}(ind,i),edy);
ind = [stc{t}{'x+p'}];
subplot(234)
hist(fet{2}(ind,i),edy);

ind = [stc{1}{'x+p'}];
subplot(232);
hist(feti{1}(ind,i),edy);
ind = [stc{t}{'x+p'}];
subplot(235)
hist(feti{2}(ind,i),edy);

ind = [stc{1}{'x+p'}];
subplot(233);
hist(fets{1}(ind,i),edy);
ind = [stc{t}{'x+p'}];
subplot(236)
hist(fets{2}(ind,i),edy);

ForAllSubplots('xlim([30,120])')


i = 1;
edy = linspace(-pi/2,pi/2,150);
figure();
ind = [stc{1}{'x+p'}];
subplot(231);
hist(fet{1}(ind,i),edy);
ind = [stc{t}{'x+p'}];
subplot(234)
hist(fet{2}(ind,i),edy);

ind = [stc{1}{'x+p'}];
subplot(232);
hist(feti{1}(ind,i),edy);
ind = [stc{t}{'x+p'}];
subplot(235)
hist(feti{2}(ind,i),edy);

ind = [stc{1}{'x+p'}];
subplot(233);
hist(fets{1}(ind,i),edy);
ind = [stc{t}{'x+p'}];
subplot(236)
hist(fets{2}(ind,i),edy);

ForAllSubplots('xlim([80,170])')