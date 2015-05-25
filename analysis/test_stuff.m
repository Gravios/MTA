
SessionName = 'Ed05-20140529'; 
MazeName = 'ont';
TrialName = 't2'; 


s = MTASession(SessionName,MazeName);

%startStopShift = [20,-1]; 
%QuickTrialSetup(s,'all',startStopShift);

ignoredViconTrials = [1,3:10];
startStopShift = [20,-1]; 
QuickTrialSetup(s,TrialName,startStopShift,ignoredViconTrials,[],false);
Trial = MTATrial(SessionName,TrialName,MazeName);


%% Visually check the Filtering of xyz with Buterworth filter
Trial = MTATrial('jg05-20120317');

xyz = Trial.load('xyz');
ind = 1:3000;
figure,hold on
plot(xyz(ind,1,1),xyz(ind,1,2),'r')
xyz.filter('ButFilter',3,1,'low');
plot(xyz(ind,1,1),xyz(ind,1,2),'b')

xyz = Trial.load('xyz');
xyz.addMarker('f_spine_lower',[.7,1,.7],{{'spine_lower','pelvis_root',[0,0,1]}},ButFilter(xyz(:,'spine_lower',:),3,1./(xyz.sampleRate/2),'low'));
xyz.addMarker('f_spine_upper',[.7,1,.7],{{'spine_lower','pelvis_root',[0,0,1]}},ButFilter(xyz(:,'spine_upper',:),3,1./(xyz.sampleRate/2),'low'));
ang = create(MTADang,Trial,xyz);


figure,hold on
plot(ang(:,1,10,3))
Lines(Trial.stc{'w'}(:),[],'b');
Lines(Trial.stc{'n'}(:),[],'g');

figure,hold on
plot(diff(unwrap(circ_dist(ang(:,1,10,1),ang(:,10,11,1)))))



fet = Trial.ang.copy;
fet.data = ang(:,1,10,3);
[ys,fs,ts] = fet_spec(Trial,fet,'mtchglong',true);

figure
imagesc(ts,fs,log10(ys.data)'),axis xy
colormap jet
caxis([-8,-1])
Lines(Trial.stc{'w',1}(:),[],'b');

[U,S,V] = svd(cov(ys(Trial.stc{'a'},:)));

m_walk = mean(ys(Trial.stc{'w'},:));

figure,
plot(m_walk*V);


wf = V([1,4],:)*ys.data';
figure
plot(wf')
Lines(Trial.stc{'w',ys.sampleRate}(:),[],'b');

wf = MTADxyz('data',wf','sampleRate',ys.sampleRate);
ind = Trial.stc{'w'};
figure,
hist2(log10(abs(wf(ind,:))),linspace(-8,-1,100),linspace(-8,-1,100))

figure,hold on
ind = Trial.stc{'a-r-w'};
bar(linspace(-8,-1,100),histc(log10(abs(wf(ind,2))),linspace(-8,-1,100)),'histc')
ind = Trial.stc{'w'};
h = bar(linspace(-8,-1,100),histc(log10(abs(wf(ind,2))),linspace(-8,-1,100)),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

figure,hist(log10(abs(wf(:,nniz(wf'))))',linspace(-8,-1,100))


Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,40,'low');

figure,plot(sqrt(sum(sum(circshift(xyz(:,1:3,[1,2]),-30)-xyz(:,1:3,[1,2]),2).^2,3)))
Lines(Trial.stc{'w'}(:),[],'b');
Lines(Trial.stc{'n'}(:),[],'g');

figure,
ni = 1:2:15;
c = 1;
for i= ni,
    subplotfit(c,numel(ni)),hold on
    c = c+1;
wf = sqrt(sum(sum(circshift(xyz(:,1:4,[1,2]),-i)-xyz(:,1:4,[1,2]),2).^2,3));
wf = MTADxyz('data',wf,'sampleRate',xyz.sampleRate);
edgs = linspace(-2,4,100);
ind = Trial.stc{'a-r-w-n'};
bar(edgs,histc(log10(abs(wf(ind,1))),edgs),'histc')
ind = Trial.stc{'w'};
h = bar(edgs,histc(log10(abs(wf(ind,1))),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;
end


figure,hold on
edgs = linspace(-2,3,100);
ind = Trial.stc{'a-r-w-n'};
bar(edgs,histc(log10(abs(wf(ind,1))),edgs),'histc')
ind = Trial.stc{'w'};
h = bar(edgs,histc(log10(abs(wf(ind,1))),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;



hcom = xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'}));
wf = [0;sqrt(sum(diff(hcom(:,1,[1,2])).^2,3))];
wf = MTADxyz('data',wf,'sampleRate',xyz.sampleRate);
wf.filter('ButFilter',3,4,'low');
ind = Trial.stc{'a'};
ind = Trial.stc{'w'};
ind = JoinRanges(Trial.stc{'w'}.data,Trial.stc{'n'}.data);
edgs = linspace(-2,2,100);
figure,hist2(log10(abs([vxy(ind,2),circshift(wf(ind,1),-0)])),edgs,edgs)

ind = JoinRanges(Trial.stc{'w'}.data,Trial.stc{'n'}.data);
ind = Trial.stc{'w'};
edgs = linspace(-2,3,100);
adgs = linspace(-6,-1,100);
figure,hist2(log10(abs([wf(ind,1),abs(circ_dist(ang(ind,1,4,1),circshift(ang(ind,1,4,1),1)))])),edgs,adgs)

figure,plot(abs(circ_dist(ang(:,1,4,1),circshift(ang(:,1,4,1),1))))



