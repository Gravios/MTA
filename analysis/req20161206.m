

%% Walk golani 1980
%QuickSessionSetup(get_session_list('jg03'));

s = MTASession('jg03-20110428');
s = MTASession('jg03-20110429');
s = MTASession('jg03-20110430');
s = MTASession('jg03-20110501');
s = MTASession('gs04-20110921');
s = MTASession('jg05-20120317');

Trial = MTATrial(s);
Trial.load('stc','NN0317');
Trial.load('stc','hand_labeled_rev3_jg');

% Select Walking periods
wper = Trial.stc{'w'};
% Compute the duration of each walk period
wdur = diff(wper.data,1,2);
% remove periods which are shorter than 1 second
wper.data(wdur<120,:)=[];


xyz = Trial.load('xyz');
rb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
hcom = xyz.com(rb);
xyz.addMarker('fbcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[1.5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('bcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
xyz.addMarker('fsl',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},xyz(:,'spine_lower',:));
rb = xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('hcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[1.5]./(xyz.sampleRate/2),'low'));
clear('hcom');

ang = create(MTADang,Trial,xyz);

% 30Hz low pass filtered xyz 
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,30,'low'); % Position
fvxy = fxyz.vel([],[1,2]);           % Speed xy
fvz = fxyz.vel([],[3]);              % Speed z
fang = create(MTADang,Trial,fxyz);   % Angles inter-marker

ffxyz = xyz.copy;
ffxyz.filter('ButFilter',3,2.4,'low'); % Position
ffvxy = ffxyz.vel([],[1,2]);           % Speed xy
ffvxy.data(ffvxy.data<1e-3)=1e-3;
ffvxy.data = log10(ffvxy.data);

% 1Hz low pass filtered xyz
flxyz = xyz.copy;
flxyz.filter('ButFilter',3,1,'low'); % Position
% 1Hz low pass filtered ang
flang = create(MTADang,Trial,flxyz);

% 1Hz High pass filtered speed
hfvxy = fvxy.copy;
hfvxy.filter('ButFilter',3,1,'high'); 

hfxyz = fxyz.copy;
hfxyz.filter('ButFilter',3,1,'high'); 

% traj NORM body vector projection
shft = 3;
tmar = {'spine_lower','pelvis_root','bcom','spine_middle','spine_upper', ...
        'head_back','head_front'};
for m = 1:numel(tmar),
    tvec(:,m,:) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,tmar{m},[1,2]),shft);
    cvec(:,m,:) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,tmar{m},[1,2]),shft);
end


% traj NORM body vector projection
bvec = xyz(:,'fbcom',[1,2])-xyz(:,'fsl',[1,2]);
ubvec = bsxfun(@rdivide,bvec,sqrt(sum(bvec.^2,3)));

% traj ORTHO body vector projection
mvec = fxyz(:,'fbcom',[1,2])-fxyz(:,'fsl',[1,2]);
umvec = bsxfun(@rdivide,bsxfun(@times,permute([1,-1],[1,3,2]),mvec(:,1,[2,1])),sqrt(sum(mvec.^2,3)));

mvec = fxyz(:,'bcom',[1,2])-fxyz(:,'spine_lower',[1,2]);
fumvec = bsxfun(@rdivide,bsxfun(@times,permute([1,-1],[1,3,2]),mvec(:,1,[2,1])),sqrt(sum(mvec.^2,3)));




nind = nniz(tvec);
for m = 1:numel(tmar),
    walkFetFilt(nind,m) = ButFilter(nunity(dot(tvec(nind,m,:),umvec(nind,:,:),3)),3,1/(xyz.sampleRate/2),'high');
    walkFet(nind,m) = nunity(dot(tvec(nind,m,:),umvec(nind,:,:),3));
    walkFetBody(nind,m) = nunity(dot(cvec(nind,m,:),ubvec(nind,:,:),3));
    walkFetComp(nind,m) = nunity(dot(cvec(nind,m,:),umvec(nind,:,:),3));
end


% traj NORM hip vector projection
% $$$ pmar = 'pelvis_root';
% $$$ pvec = circshift(xyz(:,pmar,[1,2]),-shft)...
% $$$        -circshift(xyz(:,pmar,[1,2]),shft);
% $$$ hvec = xyz(:,'hip_left',[1,2])-xyz(:,'hip_right',[1,2]);
% $$$ uhvec = bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3)));
% $$$ walkFetHips(nind) = nunity(dot(pvec(nind,:,:),uhvec(nind,:,:),3));


figure,hold on,whitebg('k')
plot(walkFet,'w-')
plot(walkFet(:,5),'g-')
plot(walkFet(:,4),'b-')
plot(walkFet(:,3),'c-')
plot(walkFet(:,2),'r-')

plot(var(walkFet,[],2),'r-')
plot(nunity(circ_dist(circshift(flang(:,'spine_lower','spine_upper',1),-3),circshift(flang(:,'spine_lower','spine_upper',1),3))),'g');
Lines(Trial.stc{'w'}(:)-.5,[],'w');
Lines(Trial.stc{'n'}(:),[],'g');


plot(walkFetBody,'m-')
plot(var(walkFet,[],2),'r-')

plot(walkFetBody,'c')
plot(walkFetComp,'m-.')
%plot(walkFetFilt,'m-.')
plot(walkFetHips,'r')
plot(nunity(fvxy(:,1)),'y')


figure,hold on,whitebg('k')
plot(walkFet,'w-')
plot(walkFetBody+5,'c')
plot(nunity(hfvxy(:,'spine_lower')),'m')
plot(nunity(hfxyz(:,'spine_lower',3)),'r');
Lines(Trial.stc{'w'}(:)-.5,[],'w');
Lines(Trial.stc{'n'}(:),[],'g');


Data = xyz.copy;
Data.data = walkFet;
Data.label = 'wfet';
Data.key = 'w';
Data.name = 'Walk Feature';
Data.ext = 'fet';
parspec = struct('nFFT',2^8,...
                 'Fs',  Data.sampleRate,...
                 'WinLength',2^7,...
                 'nOverlap',2^7*.875,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[1,20]);
[ys,fs,ts,phi,fst] = fet_spec(Trial,Data,'mtchglong',true,[],parspec,true);

% PHASE stuff
sp = [];
figure,
sp(1) = subplot(211);
imagesc(ts,fs,phi(:,:,2,5)')
Lines(Trial.stc{'w',1}(:),[],'m');
axis xy
caxis([-pi,pi])
colormap hsv;
sp(end+1) = subplot(212);
hold on;
plot(td,Data.data)
Lines(Trial.stc{'w',1}(:),[],'m');
Lines(Trial.stc{'n',1}(:),[],'g');
linkaxes(sp,'x');


% COHERENCE stuff
sp = [];
figure,
sp(1) = subplot(211);
imagesc(ts,fs,ys(:,:,2,5)')
Lines(Trial.stc{'w',1}(:),[],'m');
axis xy
caxis([0,1])
colormap jet;
sp(end+1) = subplot(212);
hold on;
plot(td,Data.data)
Lines(Trial.stc{'w',1}(:),[],'m');
Lines(Trial.stc{'n',1}(:),[],'g');
linkaxes(sp,'x');


yp = ys.copy;
yp.data = phi;


figure
states = {'a-r-m-n-s','r','m','p','n','w'};
for s = 1:numel(states)
subplot(2,3,s)
ind = Trial.stc{states{s}};
hist2([mean(log10(ys(ind,fs<4,2,2)),2),circ_mean(yp(ind,fs<4,2,5),[],2)],linspace(-6,-1,40),linspace(-pi,pi,40))
title(ind.label)
caxis([0,40])
end

td = [1:size(Data,1)]/Data.sampleRate;

sp = [];
figure,
sp(1) = subplot(211);
imagesc(ts,fs,log10(ys(:,:))')
Lines(Trial.stc{'w',1}(:),[],'m');
axis xy
caxis([-6,-1])
colormap jet;
sp(end+1) = subplot(212);
hold on;
plot(td,Data.data)
plot(td,envelope(Data.data,30,'rms'))
Lines(Trial.stc{'w',1}(:),[],'m');
linkaxes(sp,'x');

% COMPUTE envelope of walk feature
DataEnv = Data.copy;
DataEnv.data = envelope(Data.data,30,'rms');
DataEnv.name = 'Walk Feature Envelop';
DataEnv.label = 'wfete';
DataEnv.key = 'e';





wfb = xyz.copy;
wfb.data = log10(walkFetBody(:,2)+min(abs(walkFetBody(:,2)))+1e-3);

winds = LocalMinima(-abs(Data.data(:,2)),15,-.5);
windur = diff(winds)+1;
windur(end+1) = size(xyz,1)-winds(end)+1;

wdd = [];
wdist = [];
for i = 1:numel(windur)
    ind = winds(i:i+1);
    wdd(i) = sqrt(sum(diff(xyz(ind,1,[1,2])).^2,3));
    wdist(i) = wdd(i)./windur(i);
end


wdfet = [];
wddist = [];
for i = 1:numel(windur)
    inds = winds(i):winds(i)+windur(i);
    wdfet(inds) = wdist(i).*ones(size(inds));
    wddist(inds) = wdd(i).*ones(size(inds));
end

%figure,hold on,plot(td,[wdfet,1e-3,1e-3,1e-3]')

wdf = xyz.copy;
wdf.data = log10([wdfet,1e-3,1e-3,1e-3]'+randn([size(xyz,1),1])/16+1);
wd = xyz.copy;
wd.data = log10([wddist,1e-3,1e-3,1e-3]'+randn([size(xyz,1),1])/16+1);
%figure,plot(wdf.data)



edx = linspace(-0.2,1,150);
edy = linspace(-3,2,150);
figure
states = {'a-r-m-n-s','r','m','p','n','w'};
for i = 1:numel(states),
    subplot(2,3,i)
    ind =  Trial.stc{states{i}};
    %hist2([wfb(ind,1),log10(DataEnv(ind))],edx,edy)
    %hist2([wdf(ind,1),wfb(ind)],edx,edy)
    hist2([wdf(ind,1),ffvxy(ind)],edx,edy)    
    caxis([0,50])
    title(ind.label);
    grid on
end





figure,hold on
ind =  Trial.stc{states{end-2}};
h = bar(edx,histc(wdf(ind,1),edx),'histc');
h.FaceColor = 'c';h.EdgeColor = 'c';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;
ind =  Trial.stc{states{end}};
h = bar(edx,histc(wdf(ind,1),edx),'histc');
h.FaceColor = 'r';h.EdgeColor = 'r';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;

figure,hold on
ind =  Trial.stc{states{end-2}};
h = bar(edy,histc(ffvxy(ind,1),edy),'histc');
h.FaceColor = 'c';h.EdgeColor = 'c';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;
%ind =  Trial.stc{states{end}};
ind =  wper;
h = bar(edy,histc(ffvxy(ind,1),edy),'histc');
h.FaceColor = 'r';h.EdgeColor = 'r';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;




