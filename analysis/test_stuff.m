
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



dsp = struct('nFFT',2^9,'Fs',fet.sampleRate,...
             'WinLength',2^7,'nOverlap',2^7*.875,...
             'FreqRange',[1,30]);

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,10,'low');
fet = Trial.ang.copy;
fet.data = [diff([0;diff(fxyz(:,1,3))]);0];
fet.data = [diff([0;diff(fxyz(:,4,3))]);0];
fet.data = [0;diff(ang(:,1,10,3))];
[ys,fs,ts] = fet_spec(Trial,fet,'mtchglong',true,'defspec',dsp);

figure
imagesc(ts,fs,log10(ays.data)'),axis xy
colormap jet
caxis([-8,-1])
Lines(Trial.stc{'w',1}(:),[],'b');

edgs = linspace(-8,-1,100);
figure,hold on
ind = Trial.stc{'a-w-n-r'};
bar(edgs,histc(log10(ays(ind,20)),edgs),'histc');
ind = Trial.stc{'w'};
%ind = JoinRanges(Trial.stc{'w',ays.sampleRate}.data,Trial.stc{'n',ays.sampleRate}.data);
%ind = [Trial.stc{'w&a',ys.sampleRate}]+[-.5,.5];
h = bar(edgs,histc(log10(ys(ind,20)),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

-[Trial.stc{'w&a',ys.sampleRate}]+[-.5,.5];

svxy = vxy.copy;svxy.resample(ys);
svxy.filter('ButFilter',3,2,'low');
swf = wf.copy;swf.resample(ys);
swf.filter('ButFilter',3,2,'low');

ind = Trial.stc{'a'};
ind = plus([Trial.stc{'w&a',ys.sampleRate}],[-.25,.25]);
edgs = linspace(-9,-2,70);
vdgs = linspace(-2,1,70);
figure
hist2(log10(abs([ys(ind,20),swf(ind,1)])),edgs,vdgs)

figure, hold on
ind = Trial.stc{'a-n-r-w'};
bar(linspace(-4,20,100),histc(log10(abs(ys(ind,20))).*log10(abs(swf(ind,1))),linspace(-4,20,100)),'histc');
ind = Trial.stc{'w'};
h = bar(linspace(-4,20,100),histc(log10(abs(ys(ind,20))).*log10(abs(swf(ind,1))),linspace(-4,20,100)),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;


dsx = Trial.load('xyz');
vxy = dsx.vel([1,7],[1,2]);
vxy.filter('ButFilter',3,6,'low');

edgs = linspace(-1,2,100);
figure,hold on
ind = Trial.stc{'a-w-n-r'};
bar(edgs,histc(log10(abs(svxy(ind,1))),edgs),'histc');
ind = Trial.stc{'w'};
%ind = JoinRanges(Trial.stc{'w',ays.sampleRate}.data,Trial.stc{'n',ays.sampleRate}.data);
%ind = [Trial.stc{'w&a',ys.sampleRate}]+[-.5,.5];
h = bar(edgs,histc(log10(abs(svxy(ind,1))),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

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


% Grooming stuff



dsp = struct('nFFT',2^9,'Fs',fet.sampleRate,...
             'WinLength',2^7,'nOverlap',2^7*.875,...
             'FreqRange',[1,50]);

fet = Trial.ang.copy;
fet.data = abs(circ_dist(ang(:,'spine_lower','pelvis_root',2),ang(:,'spine_lower','spine_middle',2)));
[gys,fs,ts] = fet_spec(Trial,fet,'mtchglong',true,'defspec',dsp);

figure
imagesc(ts,fs,log10(gys.data)'),axis xy
colormap jet
%caxis([-8,-1])
Lines(Trial.stc{'w',1}(:),[],'b');

edgs = linspace(-8,-1,100);
figure,hold on
ind = Trial.stc{'a-w-n-r'};
bar(edgs,histc(log10(ays(ind,20)),edgs),'histc');


[U,S,V] = svd(cov(gys(Trial.stc{'a-s'},:)));

m_walk = mean(gys(Trial.stc{'m'},:));

figure,
plot(m_walk*V);

figure,imagesc(log10(abs(gys(Trial.stc{'a'},:)*V))')


dp = (mean(log10(abs(gys(Trial.stc{'m'},:)*V)))-mean(log10(abs(gys(Trial.stc{'a'},:)*V))))./sqrt(.5*(var(log10(abs(gys(Trial.stc{'m'},:)*V)))+var(log10(abs(gys(Trial.stc{'a'},:)*V)))));


edgs = linspace(-10,-1,100);
figure,hold on
bar(edgs,histc(log10(abs(gys(Trial.stc{'a-m'},:)*V(:,9))),edgs),'histc');
h = bar(edgs,histc(log10(abs(gys(Trial.stc{'m'},:)*V(:,9))),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;


% More Change point detection stuff
dsx = Trial.load('xyz');
vxy = dsx.vel([1,7],[1,2]);
vxy.filter('ButFilter',3,2,'low');

vz = vxy.copy;
vz.data = dsx(:,7,3);
vz.filter('ButFilter',3,2,'low');

dsx.filter('ButFilter',3,2,'low');
ang = create(MTADang,Trial,dsx);
%circ_dist(ang(:,1,4,1),circshift(ang(:,1,4,1),-1);
% 
% figure,plot(log10(abs(diff(abs(vxy(:,1))))))
% hold on,plot(log10(abs(diff(diff(abs(vxy(:,1)))).*50)),'r')
% Lines(Trial.stc{'n'}(:),[],'g');Lines(Trial.stc{'w'}(:),[],'b');Lines(Trial.stc{'r'}(:),[],'r');

%fet = {abs(diff(abs(vxy(:,1))))};

ofet = {};
ofet = cat(1,ofet,{MTADxyz('data',vxy(:,1),'sampleRate',dsx.sampleRate)});
ofet = cat(1,ofet,{MTADxyz('data',vz(:,1),'sampleRate',dsx.sampleRate)});
ofet = cat(1,ofet,{MTADxyz('data',abs(circ_dist(ang(:,1,4,1),circshift(ang(:,1,4,1),-1))),'sampleRate',dsx.sampleRate)});


fet = {};
for f = 1:numel(ofet),
    fet = cat(1,fet,{abs(diff(abs(ofet{f}(:,1))))});
end

f = 3;

trv =10;
vpnts = LocalMinima(-log10(abs(fet{f})),round(.20*dsx.sampleRate),trv);
vprs = [vpnts,circshift(vpnts,-1)];
vprs(end,:) = [];

vstd = [];
vtoc = diff(vprs,1,2);
vprs(vtoc>round(10*dsx.sampleRate),:) = [];
vtoc(vtoc>round(10*dsx.sampleRate)) = [];
for i = vprs'    
    vstd(end+1,1) = std(abs(ofet{f}(i(1):i(2),1)));
end
% 
% i = [vprs(1,1),vprs(2,2)];
% vstd(end+1,1) = std(vxy(i(1):i(2),1));

Lsc = sum(vtoc.*log(vstd.^2));

Lsc = vxy.size(1)*sum(log(vtoc.*vstd.^2));

mpoint = zeros([size(vprs,1)-1,1]);
%for t = 1:size(vprs,1)-1,
lchng = [];
Lsc = Lsc(1);

t = 1;
while t<size(vprs,1),
    tvtoc = vtoc;
    tvstd = vstd;
    tvtoc(t+1) = vprs(t+1,2)-vprs(t,1);
    tvtoc(t) = [];
    tvstd(t+1) = std(abs(ofet{f}(vprs(t,1):vprs(t+1,2),1)));
    tvstd(t) = [];
    tLsc = vxy.size(1)*sum(log(tvtoc.*tvstd.^2));
    %tLsc = sum(tvtoc.*log(tvstd.^2));
    if Lsc(end)>tLsc;
        Lsc(end+1) = Lsc(end);
        t = t+1;
    else
        Lsc(end+1) = tLsc;
        vtoc = tvtoc;
        vstd = tvstd;
        vprs(t+1,:) = [vprs(t,1),vprs(t+1,2)];
        vprs(t,:) = [];
        
    end
end

figure,plot(abs(ofet{f}(:,1)))
Lines(vprs(1:100),[],'g');


f = 2;

    vpnts = LocalMinima(-log10(abs(fet{f})),round(.20*dsx.sampleRate),thr_fet(f));
    vprs = [vpnts,circshift(vpnts,-1)];
    vprs(end,:) = [];
    
%vxy = dsx.vel([1,7],[1,2]);
vsc = zeros([size(ofet,1),1]);
for i = vprs'    
    vsc(i(1):i(2),1) = median(ofet(i(1):i(2),f));
end

wvsc = [];
for i = Trial.stc{'w'}.data'    
    %wvsc(end+1) = median(unique(ofet{f}(i(1):i(2),1)));
    wvsc = [wvsc;unique(vsc(i(1):i(2),1))];
end

svsc = [];
for i = Trial.stc{'r'}.data'    
    %svsc(end+1) = median(unique(ofet{f}(i(1):i(2),1)));
    svsc = [svsc;unique(vsc(i(1):i(2),1))];
    %svsc(end+1) = median(unique(vsc(i(1):i(2),1)));
end

edgs = linspace(-5,2,100);
%edgs = linspace(1.5,2.6,100);
figure,hold on
bar(edgs,histc(log10(svsc),edgs),'histc');
h = bar(edgs,histc(log10(wvsc),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

figure,plot(vsc)
Lines(Trial.stc{'n'}(:),[],'g');Lines(Trial.stc{'w'}(:),[],'b');Lines(Trial.stc{'r'}(:),[],'r');

vsc = MTADxyz('data',vsc,'sampleRate',vxy.sampleRate);

edgs = linspace(-2,2,100);
figure,hold on
bar(edgs,histc(log10(vsc(Trial.stc{'a-w'},1)),edgs),'histc');
h = bar(edgs,histc(log10(vsc(Trial.stc{'w'},1)),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

edgs = linspace(-2,2,100);
figure,hold on
bar(edgs,histc(unique(log10(vsc(Trial.stc{'s'},1))),edgs),'histc');
h = bar(edgs,histc(unique(log10(vsc(Trial.stc{'w'},1))),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;



vs = vxy.segs(1:vxy.size(1),30,nan);
vs = sq(std(vs));
figure,plot(vs(:,1))
Lines(Trial.stc{'w'}(:),[],'b');

vs = MTADxyz('data',vs,'sampleRate',vxy.sampleRate);

edgs = linspace(-1,2,100);
figure,hold on
bar(edgs,histc(log10(vs(Trial.stc{'a-w'},1)),edgs),'histc');
h = bar(edgs,histc(log10(vs(Trial.stc{'s'},1)),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;


f = 3;
edgs = linspace(-2,2,100);
edgs = linspace(1.5,2.6,100);
edgs = linspace(-8,0,100);

cind = 1;
c = [];
for i = Trial.stc{'n'}.data',
    c(cind,:) = histc(log10(ofet{f}(i',1)),edgs);
    cind = cind + 1;
end
c = bsxfun(@rdivide,c,sum(c,2));
figure,bar(edgs,mean(c),'histc')


cind = 1;
s = [];
for i = Trial.stc{'m'}.data',
    s(cind,:) = histc(log10(ofet{f}(i',1)),edgs);
    cind = cind + 1;
end
s = bsxfun(@rdivide,s,sum(s,2));
figure,bar(edgs,mean(s),'histc')

figure,hold on
bar(edgs,mean(s),'histc');
h = bar(edgs,mean(c),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

ms = sum(mean(s).*edgs);
ss = cumsum(mean(s));
ss = abs(ms-edgs(find(.98>fliplr(ss),1,'first'))/2);

mc = sum(mean(c).*edgs);
sc = cumsum(mean(c));
sc = abs(mc-edgs(find(.98>fliplr(sc),1,'first'))/2);

(mc-ms)/sqrt(.5*(ss.^2+sc.^2))

edgs = linspace(-1,1,100);
m = histc(circ_dist(ang(Trial.stc{'s'},2,3,2),ang(Trial.stc{'s'},3,4,2)),edgs);
m = m/sum(m);
n = histc(circ_dist(ang(Trial.stc{'n'},2,3,2),ang(Trial.stc{'n'},3,4,2)),edgs);
%n = histc(log10(ofet{f}(Trial.stc{'n'},1)),edgs);
n = n/sum(n);

figure,hold on,
bar(edgs,m,'histc');
h = bar(edgs,n,'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;


% Had a Thought, How does the sitting periods distributions relate to their
% durations
Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
vxy = xyz.vel(1,[1,2]);
vxy.filter('ButFilter',3,4,'low');

ms = [];
vs = [];
ds = diff(Trial.stc{'w'}.data,1,2);
for s = Trial.stc{'w'}.data',
    ms(end+1) = mean(vxy(s'));
    vs(end+1) = var(vxy(s'));
end

figure,plot(log10(ms),log10(ds),'.')
figure,plot(log10(ms),log10(vs),'.')
figure,plot(log10(ds),log10(ms.*vs),'.')

figure,hist2(log10([ds,vs']),50,50)



%% Checking out Pfs stuff to see if it still works
pfs = {};
pfs{1} = MTAApfs(Trial,[],'gper');
pfs{2} = MTAApfs(Trial,[],'rear');
pfs{3} = MTAApfs(Trial,[],'groom');
pfs{4} = MTAApfs(Trial,[],'walk');


units = pfs{1}.data.clu;
unit = 1;
spo = [1,2,3,4];
hfig = figure(13939);
while unit~=-1,    
    for s = 1:numel(pfs),
        subplot(2,2,spo(s));
        pfs{s}.plot(unit);
        title([pfs{s}.parameters.states ' :' num2str(unit)]);
    end
    unit = figure_controls(hfig,unit,units);
end    

        
Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
vl = xyz.vel(1,[1,2]);
vl.filter('ButFilter',3,2,'low');
vl.data(nniz(vl)) = log10(vl(nniz(vl)));

ed = linspace(-.5,2,100);
figure,hold on
ind = Trial.stc{'a-w'};
bar(ed,histc(vl(ind),ed),'histc')
ind = Trial.stc{'w'};
h = bar(ed,histc(vl(ind),ed),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

wper = ThreshCross(vl.data,.5,25);
figure,hist(log10(diff(wper,1,2)),100)
figure,hist(log10(diff(Trial.stc{'w'}.data,1,2)),100)

wmax = [];
for w = wper',wmax(end+1,1)=max(vl(w'));end
wdur = abs(log10(diff(wper,1,2)));
figure,hist2([wdur,wmax],100,100)




Trial = MTATrial('jg05-20120317');

xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');
fang = create(MTADang,Trial,fxyz);

figure,plot(circ_dist(ang(:,1,4,1),fang(:,1,4,1)))
Lines(Trial.stc{'w'}(:),[],'b');



wf = MTADxyz('data',wf,'sampleRate',trajSampleRate);

figure,hold on,
ind = Trial.stc{'a-r-w'};
ha = bar(linspace(-6,6,500),histc(wf(ind),linspace(-6,6,500)),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .4;
ha.EdgeColor = 'c';
ha.EdgeAlpha = .4;
ind = Trial.stc{'w'};
hs = bar(linspace(-6,6,500),histc(wf(ind),linspace(-6,6,500)),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeColor = 'r';
hs.EdgeAlpha = .4;


load(fullfile(Trial.spath,...
    [Trial.filebase '-walk_fet_ppc.mat']))

man = Trial.xyz.copy;
man.data = mag;
man.filter('ButFilter',3,1,'low');

man.resample(trajSampleRate);

fvel = xyz.vel([],[1,2]);
fvel.filter('ButFilter',3,2.5,'low');
fvel.resample(trajSampleRate);
fvel.data(fvel.data<0)=.1;
fvel.data = log10(fvel.data);



figure
subplot(121)
ind = Trial.stc{'a-w-r-n'};
hist2([man(ind),wf(ind)],linspace(-.2,1,80),linspace(-6,6,80));
subplot(122)
ind = Trial.stc{'w+n'};
hist2([man(ind),wf(ind)],linspace(-.2,1,80),linspace(-6,6,80));

figure
subplot(121)
ind = Trial.stc{'a-w-r'};
hist2([fvel(ind,1),wf(ind)],linspace(-.7,2,80),linspace(-6,6,80));
subplot(122)
ind = Trial.stc{'w'};
hist2([fvel(ind,1),wf(ind)],linspace(-.7,2,80),linspace(-6,6,80));


ind = Trial.stc{'w'};
hist2([fvel(ind,1),wf(ind)],linspace(-.7,2,80),linspace(-6,6,80));





