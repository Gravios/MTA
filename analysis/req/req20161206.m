

%% Walk golani 1980
%QuickSessionSetup(get_session_list('jg03'));

s = MTASession('jg03-20110428');
s = MTASession('jg03-20110429');
s = MTASession('jg03-20110430');
s = MTASession('jg03-20110501');
s = MTASession('gs04-20110921');
s = MTASession('jg05-20120317');
%s = MTASession('Ed03-20140624');

Trial = MTATrial(s);
Trial.load('stc','NN0317');
Trial.load('stc','hand_labeled_rev3_jg');
%Trial.load('stc','hand_labeled_rev1_Ed');



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
% $$$ fxyz = xyz.copy;
% $$$ fxyz.filter('ButFilter',3,30,'low'); % Position
% $$$ fvxy = fxyz.vel([],[1,2]);           % Speed xy
% $$$ fvz = fxyz.vel([],[3]);              % Speed z
% $$$ fang = create(MTADang,Trial,fxyz);   % Angles inter-marker
% $$$ 
% $$$ ffxyz = xyz.copy;
% $$$ ffxyz.filter('ButFilter',3,2.4,'low'); % Position
% $$$ ffvxy = ffxyz.vel([],[1,2]);           % Speed xy
% $$$ ffvxy.data(ffvxy.data<1e-3)=1e-3;
% $$$ ffvxy.data = log10(ffvxy.data);
% $$$ 
% 1Hz low pass filtered xyz
flxyz = xyz.copy;
flxyz.filter('ButFilter',3,1,'low'); % Position
% 1Hz low pass filtered ang
flang = create(MTADang,Trial,flxyz);
% $$$ 
% $$$ % 1Hz High pass filtered speed
% $$$ hfvxy = fvxy.copy;
% $$$ hfvxy.filter('ButFilter',3,1,'high'); 
% $$$ 
% $$$ hfxyz = fxyz.copy;
% $$$ hfxyz.filter('ButFilter',3,1,'high'); 

% traj NORM body vector projection
shft = 3;
tmar = {'spine_lower','pelvis_root','bcom','spine_middle','spine_upper', ...
        'hcom'};
tvec = [];cvec = [];
for m = 1:numel(tmar),
    tvec(:,m,:) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,tmar{m},[1,2]),shft);
    cvec(:,m,:) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,tmar{m},[1,2]),shft);
end


% traj NORM body vector projection
bvec = xyz(:,'fbcom',[1,2])-xyz(:,'fsl',[1,2]);
ubvec = bsxfun(@rdivide,bvec,sqrt(sum(bvec.^2,3)));

% traj ORTHO body vector projection
mvec = xyz(:,'fbcom',[1,2])-xyz(:,'fsl',[1,2]);
umvec = bsxfun(@rdivide,bsxfun(@times,permute([1,-1],[1,3,2]),mvec(:,1,[2,1])),sqrt(sum(mvec.^2,3)));

mvec = xyz(:,'bcom',[1,2])-xyz(:,'spine_lower',[1,2]);
fumvec = bsxfun(@rdivide,bsxfun(@times,permute([1,-1],[1,3,2]),mvec(:,1,[2,1])),sqrt(sum(mvec.^2,3)));




nind = nniz(tvec);
for m = 1:numel(tmar),
    walkFetFilt(nind,m) = ButFilter(nunity(dot(tvec(nind,m,:),umvec(nind,:,:),3)),5,[1,6]/(xyz.sampleRate/2),'bandpass');
    walkFetFiltH(nind,m) = ButFilter(nunity(dot(tvec(nind,m,:),umvec(nind,:,:),3)),5,[.5]/(xyz.sampleRate/2),'high');
    walkFetFiltL(nind,m) = ButFilter(nunity(dot(tvec(nind,m,:),umvec(nind,:,:),3)),5,[1]/(xyz.sampleRate/2),'low');
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

plot(var(walkFet,[],2)+10,'r-')
plot(nunity(circ_dist(circshift(flang(:,'spine_lower','spine_upper',1),-3),circshift(flang(:,'spine_lower','spine_upper',1),3)))+10,'g');
Lines([],10,'w');
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
%Data.data = walkFetFilt;
Data.data = walkFet;
Data.label = 'wfet';
Data.key = 'w';
Data.name = 'Walk Feature';
Data.ext = 'fet';
parspec = struct('nFFT',2^9,...
                 'Fs',  Data.sampleRate,...
                 'WinLength',2^8,...
                 'nOverlap',2^8*.875,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[0.1,8]);
[ys,fs,ts] = fet_spec(Trial,Data,'mtcsdglong',true,[],parspec,true);

rs = ys.copy();
rs.clear();
for i = 1:6
for j = 1:6
    rs.data(:,:,i,j) = ys(:,:,i,j)./(sqrt(ys(:,:,i,i)).*sqrt(ys(:,:,j,j)));
end
end

hsp = rs.copy();
hsp.data = atan(imag(rs(:,10,1,3))./real(rs(:,10,1,3)));

edx = linspace(-pi/2,pi/2,100);
figure,hold on
ind =  Trial.stc{'a-w'};
h = bar(edx,histc(hsp(ind),edx),'histc');
h.FaceColor = 'c';h.EdgeColor = 'c';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;
ind =  Trial.stc{'w'};
h = bar(edx,histc(hsp(ind),edx),'histc');
h.FaceColor = 'r';h.EdgeColor = 'r';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;


figure
states = {'a','s','m','p','n','w','r'};
for s = 1:numel(states)
subplot(3,3,s)
ind = Trial.stc{states{s}};
hist2([log10(ys(ind,10,2,2).^2),hsp(ind)],linspace(-20,-4,40),linspace(-pi/2,pi/2,40))
title(ind.label)
caxis([0,40])
grid on
set(gca,'GridColor',[1,1,1])
end



figure,plot(hsp)
Lines(Trial.stc{'w',rs.sampleRate}(:),[],'m');
Lines(Trial.stc{'n',rs.sampleRate}(:),[],'g');


% POWER stuff

td = [1:size(Data,1)]/Data.sampleRate;

m = 2;
sp = [];
figure,whitebg(gcf,'k')
sp(1) = subplot(211);
title('Power Spectrum of body sway')
imagesc(ts,fs,log10(ys(:,:,m,m))')
Lines(Trial.stc{'w',1}(:),[],'m');
axis xy
caxis([-6,-1])
colormap jet;
sp(end+1) = subplot(212);
title('Normalized body sway during turn and walk');
hold on;
plot(td,Data.data)
%plot(td,envelope(Data.data,30,'rms'))
Lines(Trial.stc{'w',1}(:),[],'m');
Lines(Trial.stc{'n',1}(:),[],'g');
linkaxes(sp,'x');
legend(tmar)



% PHASE stuff
sp = [];
figure,
sp(1) = subplot(211);
imagesc(ts,fs,phi(:,:,1,2)')
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

winds = LocalMinima(-abs(Data.data(:,2)),15,-1);

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




% Can I use the peaks of 
WData = Data.copy;
WData.data = zeros(size(Data));
winds = {};
for i = 1:numel(tmar)
    winds{i} = LocalMinima(-abs(Data.data(:,i)),15,-.5);
    for j = winds{i}',
        try
            WData.data(j-5:j+5,i) = 1;
        end
    end
end

figure,imagesc(td,1:4,WData.data(:,1:4)')
Lines(Trial.stc{'w',1}(:),[],'m');
Lines(Trial.stc{'n',1}(:),[],'g');
Lines(Trial.stc{'r',1}(:),[],'c');
xlim([940,970]);


figure,hold on,
plot(nunity(ang(:,'spine_lower','fbcom',2))*5)
plot(walkFet(:,2))
Lines(Trial.stc{'w'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');



states = {'a-r-m-n-s','r','m','p','n','w'};
edx = linspace(-0.2,1,150);
figure,hold on
ind =  Trial.stc{states{end-2}};
h = bar(edx,histc(ang(ind,'spine_lower','fbcom',2),edx),'histc');
h.FaceColor = 'c';h.EdgeColor = 'c';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;
ind =  Trial.stc{states{end}};
h = bar(edx,histc(ang(ind,'spine_lower','fbcom',2),edx),'histc');
h.FaceColor = 'r';h.EdgeColor = 'r';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;




edx = linspace(-0.2,1.5,150);
edy = linspace(-3,2,150);
figure
states = {'a-r-m-n-s','r','m','p','n','w'};
for i = 1:numel(states),
    subplot(2,3,i)
    ind =  Trial.stc{states{i}};
    %hist2([wfb(ind,1),log10(DataEnv(ind))],edx,edy)
    %hist2([wdf(ind,1),wfb(ind)],edx,edy)
    %hist2([wdf(ind,1),ffvxy(ind)],edx,edy)    
    hist2([flang(ind,'spine_lower','bcom',2),ffvxy(ind)],edx,edy)    
    caxis([0,50])
    title(ind.label);
    grid on
end


segLength = 20;
dp = Data.phase([1,3]);

ddp = dp.copy;
ddp.data = bsxfun(@circ_dist,dp.data,permute(dp.data,[1,3,2]));

figure,hold on
plot(ddp.data(:,1,2))
plot(ddp.data(:,1,4))
plot(ddp.data(:,1,5))
plot(ddp.data(:,2,4))
plot(ddp.data(:,2,5))
plot(ddp.data(:,4,5))

plot(Data.data(:,[1,2]))
plot(Data.data(:,[2,3]))
Lines(Trial.stc{'w'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');

figure,hold on
plot(circ_dist(circshift(ddp.data(:,1,2),1),circshift(ddp.data(:,1,2),-1)))
plot(Data.data(:,[1,2]))
Lines(Trial.stc{'w'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');


figure,hold on
dpps = [ddp.data(:,1,2),ddp.data(:,1,4),ddp.data(:,1,5),ddp.data(:,2,4),ddp.data(:,2,5),ddp.data(:,4,5)];

dppcv = circ_var(dpps,[],[],2);


figure,hold on
plot(dppcv)
Lines(Trial.stc{'w'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');


ddps = ddp.segs([],segLength);

rinds = reshape(triu(reshape([1:size(dp,2)^2]',[size(dp,2),size(dp,2)]),1),[],1);
rinds(rinds==0)=[];

ddps = reshape(ddps,segLength,[],size(dp,2)^2);
ddps = ddps(:,:,rinds);

ddpsppc = zeros([size(ddps,2),size(ddps,3)]);
for j = 10:numel(rinds),
for i = 1:size(dp,1),
    ddpsppc(i,j) = PPC(ddps(:,i,j));
end
end



dp = Data.phase([1,3]);
ddp = dp.copy;
ddp.data = bsxfun(@circ_dist,dp.data,permute(dp.data,[1,3,2]));
ddp.data = ddp(:,1:4,1:4);

rinds = reshape(triu(reshape([1:16]',[4,4]),1),[],1);
rinds(rinds==0)=[];

ddps = reshape(ddps,segLength,[],16);
ddps = ddps(:,:,rinds);


ddprppc = zeros([size(ddps,2),1]);
for i = 1:size(dp,1),
    ddprppc(i) = PPC(reshape(ddps(:,i,:),[],1));
end


figure,plot(ddprppc)
Lines(Trial.stc{'w'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');

figure,hold on,
ind = Trial.stc{'a-w'};
plot(Data(ind,1),Data(ind,2))
ind = Trial.stc{'w'};
plot(Data(ind(1:10,:),1),Data(ind(1:10,:),2),'r')





s = MTASession.validate('JM11-20170112.sof.all');
Trial = s;


xyz = Trial.load('xyz');
rb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
hcom = xyz.com(rb);
xyz.addMarker('fbcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[1]./(xyz.sampleRate/2),'low'));
xyz.addMarker('bcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
xyz.addMarker('fsl',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_lower',:),4,[1.5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('fsm',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_middle',:),4,[1.5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('fsu',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_upper',:),4,[1.5]./(xyz.sampleRate/2),'low'));
rb = xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('hcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[1.5]./(xyz.sampleRate/2),'low'));
clear('hcom');

ang = create(MTADang,Trial,xyz);

fxyz = xyz.copy();
%fxyz.filter('ButFilter',3,2.4);
vxy = fxyz.vel([],[1,2]);
%vxy.data(vxy.data<=1e-3)=1e-3;

% traj NORM body vector projection
shft = 3;
tmar = {'spine_lower','pelvis_root','spine_middle','spine_upper','hcom','knee_left','knee_right'};
tvec = [];
for m = 1:numel(tmar),
    tvec(:,m,:) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,tmar{m},[1,2]),shft);
end
nind = nniz(tvec);

% body vector
mvec = xyz(:,'spine_upper',[1,2])-xyz(:,'fsl',[1,2]);
umvec = bsxfun(@rdivide,bsxfun(@times,permute([1,-1],[1,3,2]),mvec(:,1,[2,1])),sqrt(sum(mvec.^2,3)));

unvec = [];
rotationAngles = deg2rad([0,45,90,135]);
for theta = rotationAngles,
    rotMat = repmat(permute([cos(theta),-sin(theta);sin(theta),cos(theta)],[3,1,2]),[size(mvec,1),1,1]);
    unvec(:,end+1,:) = bsxfun(@rdivide,multiprod(mvec,rotMat,[2,3],[2,3]),sqrt(sum(mvec.^2,3)));
end

walkFetRot = [];
for t = rotationAngles;
    for m = 1:numel(tmar),
        walkFetRot(nind,t==rotationAngles,m) = dot(tvec(nind,m,:),unvec(nind,t==rotationAngles,:),3);
        %walkFetRot(nind,t==rotationAngles,m) = nunity(dot(tvec(nind,m,:),unvec(nind,t==rotationAngles,:),3));
    end
end


% DISPLAY spine sway along spine
OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'MjgEd2016/manuscript/Figures/Figure_3';
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)

ts = [1:size(walkFetRot,1)]./180;

hfig = figure(gen_figure_id);
hfig.Units = 'centimeters';
hfig.Position = [2,2,10,6];
hax = axes('Units','centimeters','Position',[1,1,8,4])
hold('on')
plot(ts,walkFetRot(:,3,1)/6*180/10),
plot(ts,walkFetRot(:,3,2)/6*180/10),
plot(ts,walkFetRot(:,3,3)/6*180/10),
plot(ts,walkFetRot(:,3,4)/6*180/10),
xlim([66360,67267]./180)
ylim([-35,35]);
FigName = 'rat_walk_bodyVec90_markerTraj_feature_demo_';
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));


hfig = figure(gen_figure_id);
hfig.Units = 'centimeters';
hfig.Position = [2,2,10,6];
hax = axes('Units','centimeters','Position',[1,1,8,4])
hold('on')
plot(ts,walkFetRot(:,1,1)/6*180/10,'LineWidth',1),
plot(ts,walkFetRot(:,1,2)/6*180/10,'LineWidth',1),
plot(ts,walkFetRot(:,1,3)/6*180/10,'LineWidth',1),
plot(ts,walkFetRot(:,1,4)/6*180/10,'LineWidth',1),
xlim([66360,67267]./180)
ylim([-10,75]);
Lines([],0,'k');
FigName = 'rat_walk_bodyVec0_markerTraj_feature_demo_';
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));


hfig = figure(gen_figure_id);
hfig.Units = 'centimeters';
hfig.Position = [2,2,10,6];
hax = axes('Units','centimeters','Position',[1,1,8,4])
hold('on')
plot(ts,walkFetRot(:,1,1)/6*180/10,'LineWidth',1),
plot(ts,walkFetRot(:,1,2)/6*180/10,'LineWidth',1),
plot(ts,walkFetRot(:,1,3)/6*180/10,'LineWidth',1),
plot(ts,walkFetRot(:,1,4)/6*180/10,'LineWidth',1),
xlim([66360,67267]./180)
ylim([-10,75]);
Lines([],0,'k');
FigName = 'rat_walk_bodyVec0_markerTraj_feature_demo_';
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));



fslf = ButFilter(ang(:,'spine_lower','fsl',3),3,[10]./(xyz.sampleRate/2),'low');
fsuf = ButFilter(ang(:,'spine_upper','fsu',3)*2,3,[10]./(xyz.sampleRate/2),'low');
plot(ButFilter(walkFetRot(:,3,1),3,[0.5,6]./(xyz.sampleRate/2),'bandpass').*5+35);
plot(fslf,'b')
plot(fsuf,'g')
Lines(Trial.stc{'w'}(:)-.5,[],'b');
Lines(Trial.stc{'n'}(:),[],'g');
Lines(Trial.stc{'m'}(:),[],'m');
Lines(Trial.stc{'r'}(:)-.5,[],'r');


% DISPLAY spine sway along spine
figure,hold on,
plot(walkFetRot(:,3,1).*5+40),
plot(walkFetRot(:,3,2).*5+50),
plot(walkFetRot(:,3,3).*5+60),
plot(walkFetRot(:,3,4).*5+70),
fslf = ButFilter(ang(:,'spine_lower','fsl',3),3,[10]./(xyz.sampleRate/2),'low');
fsuf = ButFilter(ang(:,'spine_upper','fsu',3)*2,3,[10]./(xyz.sampleRate/2),'low');
plot(ButFilter(walkFetRot(:,3,1),3,[0.5,6]./(xyz.sampleRate/2),'bandpass').*5+35);
plot(fslf,'b')
plot(fsuf,'g')
Lines(Trial.stc{'w'}(:)-.5,[],'b');
Lines(Trial.stc{'n'}(:),[],'g');
Lines(Trial.stc{'m'}(:),[],'m');
Lines(Trial.stc{'r'}(:)-.5,[],'r');

% PRINCOMP of spine sway during walk and turns
nz = nniz(xyz);
embeddingWindow = 64;
wfet = xyz.copy;
wfet.data= zeros([size(xyz,1),size(walkFetRot,2)*size(walkFetRot,3)]);
wfet.data(nz,:) = [reshape(walkFetRot(nz,:),[],size(walkFetRot,2)*size(walkFetRot,3))];
wfs = wfet.segs([],embeddingWindow);
wfs = circshift(wfs,embeddingWindow/2,2);
wfs = MTADxyz('data',reshape(permute(wfs,[2,1,3]),size(wfs,2),[]),'sampleRate',xyz.sampleRate);
wfs.data(isnan(wfs.data(:)))=0;
[Uw,Sw,Vw] = svd(wfs(Trial.stc{'w+n'},:),0);

% DISPLAY eigen vectors
figure,
for i = 1:25,
    subplot(5,5,i);imagesc(reshape(Vw(:,i),[],size(wfet,2))'),
    axis xy;
end

fetW = MTADxyz('data',wfs.data*Vw(:,1),'sampleRate',xyz.sampleRate);
for i = 1:10,fetW.data(:,i) = wfs.data*Vw(:,i);end




figure,hold on
plot(fetW(:,4)),
plot(nunity(ang(:,'fsl','knee_right',3),[],70,10)+10)
plot(nunity(ang(:,'fsl','knee_left',3),[],57,10)+10)
plot(nunity(ang(:,'fsm','shldr_right',3),[],83,5)+20)
plot(nunity(ang(:,'fsm','shldr_left',3),[],83,5)+20)
% $$$ 
% $$$ plot(ang(:,'spine_lower','knee_right',3))
% $$$ plot(ang(:,'spine_lower','knee_left',3))
plot(ButFilter(fetW(:,4),3,[1,4]./(xyz.sampleRate/2),'bandpass').*-fetW(:,1)./10,'m');
plot(bsxfun(@plus,sq(walkFetRot(:,3,[1,2,3,4])),[-1,1,2,3]).*10)



%plot(xyz(:,'spine_lower',3)-30)

plot(ang(:,'spine_middle','shldr_right',3)+20)
plot(ang(:,'spine_middle','shldr_left',3)+20)


fslf = ButFilter(ang(:,'spine_lower','fsl',3),3,[8]./(xyz.sampleRate/2),'low');
fsuf = ButFilter(ang(:,'spine_upper','fsu',3)*2,3,[8]./(xyz.sampleRate/2),'low');
plot(fslf)
plot(fsuf)

%xlm = [65400,67500]
xlm = [0,7000];
xlim([xlm])
%plot(ang(:,'spine_lower','knee_left',2)*10)
%usMin = LocalMinima(-fslf,10,-2);
usMin = LocalMinima(-fsuf,10,-2);
usMin = usMin(usMin>xlm(1)&usMin<xlm(2)) ;
Lines(usMin,[],'k');

videofile = ['/storage/share/exchange/gravio/optitrack/Session\ 2017-01-12/Take\ 2017-01-12\ 02.52.14\ PM-Camera\ 13029.avi'];
videofile = '/storage/gravio/data/project/general/JM11-20170112/Trial001_C13029.avi';
videofile_side = '/storage/gravio/data/project/general/JM11-20170112/Trial001_C12881.avi';
vdr = VideoReader(videofile);
vdr_side = VideoReader(videofile_side);



% Figure Settings ----------------------------------------------------------------------
OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'MjgEd2016/manuscript/Figures/Figure_3';
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)
% --------------------------------------------------------------------------------------


hfig = figure(20170131);clf
set(hfig,'Units','centimeters',...
         'PaperPositionMode', 'auto',...
         'Position',[0,0,21,30])
index = 66600;
timePoints = [0,15,45,75,85];
imageCloseUpBox = [195,600;170,340]
for i = 1:5;
    vdr_side.CurrentTime = (index+timePoints(i))/vdr_side.FrameRate;
    vs = readFrame(vdr_side);
    hax = axes('Units','centimeters',...
               'Position',[2+3*(i-1)+i/5,20,3,1.515]);
    imagesc(rot90(vs(imageCloseUpBox(2,1):imageCloseUpBox(2,2),...
                     imageCloseUpBox(1,1):imageCloseUpBox(1,2),1))');
    hax.XTickLabel = {};
    hax.YTickLabel = {};
end


vsTimeStamps = (1:size(ang,1))./vdr_side.FrameRate;


hax = axes('Units','centimeters',...
           'Position',[2,15.5,16,3.5]);

hold on 
%plot(fetW(:,4)),
plot(vsTimeStamps,nunity(ang(:,'fsl','knee_right',3),[],70,5)+15,'LineWidth',1)
plot(vsTimeStamps,nunity(ang(:,'fsl','knee_left',3),[],57,5)+15,'LineWidth',1)
plot(vsTimeStamps,nunity(ang(:,'fsm','shldr_right',3),[],83,2.5)+30,'LineWidth',1)
plot(vsTimeStamps,nunity(ang(:,'fsm','shldr_left',3),[],83,2.5)+30,'LineWidth',1)
% $$$ 
% $$$ plot(ang(:,'spine_lower','knee_right',3))
% $$$ plot(ang(:,'spine_lower','knee_left',3))
plot(vsTimeStamps,ButFilter(fetW(:,4),3,[1,4]./(xyz.sampleRate/2),'bandpass').*-fetW(:,1)./20,'b','LineWidth',1);
xlim([index-50,index+450]./vdr_side.FrameRate);
ylim([-10,40])


hold on,
plot(vsTimeStamps,ang(:,'spine_lower','knee_right',3))
plot(vsTimeStamps,ang(:,'spine_lower','knee_left',3))
xlim([index-50,index+450]./vdr_side.FrameRate);
ylim([45,90])
Lines((index+timePoints)./vdr_side.FrameRate,[],'k');
hax.XTickLabel = {};
ylabel('Distance');

hax = axes('Units','centimeters',...
           'Position',[2,18.9,16,2.5]);
hold on,
plot(vsTimeStamps,walkFet(:,[1,2,3]))
xlim([index-50,index+450]./vdr_side.FrameRate);
ylim([-1.5,1.5])
hax.XTickLabel = {};
ylabel('Spine Sway (AU)');

hax = axes('Units','centimeters',...
           'Position',[2,16.3,16,2.5]);
plot(vsTimeStamps,nunity(xyz(:,{'spine_lower','pelvis_root','spine_middle'},3)))
xlim([index-50,index+450]./vdr_side.FrameRate);
ylim([-0.5,1.75])
ylabel('Height (AU)');
hl = legend({'Body Lower','Body Pelvis','Body Middle'},'location','SouthOutside')
hl.Units = 'Centimeters';
hl.Position(2)=13;

xlabel('Time (s)');
FigName = 'walk_feature_demo_';

print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));


% $$$ figure,
% $$$ plot(circ_dist(ang(:,'bcom','hip_right',1),ang(:,'bcom','knee_right',1)))
% $$$ hold on,
% $$$ plot(circ_dist(ang(:,'bcom','hip_left',1),ang(:,'bcom','knee_left',1)))


%% Not sure what I'm doing here, written while writing.
[xyz,ss] = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD');    


figure,
hold on
for i = 1001:100:10000,
    mx = mean(ss(i,:,1));
    my = mean(ss(i,:,2));
    xy = [sq(ss(i,:,1))'-mx,sq(ss(i,:,2))'-my];
    theta = atan2(xy(1,2),xy(1,1));
% $$$     theta = sign(theta-pi/2)*theta;
    rotMat = [cos(theta),-sin(theta);sin(theta),cos(theta)];
     for j = 1:size(xy,1),
         xy(j,:) = xy(j,:)*rotMat;
     end    
    plot(xy(:,1),xy(:,2),'.')
end



figure,
hold on
for s = Trial.stc{'w'}.data'
for i = s(1):50:s(2),
    mx = mean(ss(i,55,1));
    my = mean(ss(i,55,2));
    xy = [sq(ss(i,:,1))'-mx,sq(ss(i,:,2))'-my];
    theta = atan2(xy(1,2),xy(1,1));
    rotMat = [cos(theta),-sin(theta);sin(theta),cos(theta)];
     for j = 1:size(xy,1),
         xy(j,:) = xy(j,:)*rotMat;
     end    
    plot3(xy(:,1),xy(:,2),sq(ss(i,:,3))','-')
end
end

daspect([1,1,1])

figure
sp = [];
sp(end+1) = subplot(311);plot(walkFet);
Lines(Trial.stc{'w'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');
sp(end+1) = subplot(312);plot(walkFetFiltH);
Lines(Trial.stc{'w'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');
sp(end+1) = subplot(313);plot(walkFetFiltL);
Lines(Trial.stc{'w'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');
linkaxes(sp,'xy');



signal = walkFet(:,1);

lev   = 16;
wname = 'db1';
nbcol = 64;
[c,l] = wavedec(signal,lev,wname);

len = length(signal);
cfd = zeros(lev,len);
for k = 1:lev
    d = detcoef(c,l,k);
    d = d(:)';
    d = d(ones(1,2^k),:);
    cfd(k,:) = wkeep1(d(:)',len);
end
cfd =  cfd(:);
I = find(abs(cfd)<sqrt(eps));
cfd(I) = zeros(size(I));
cfd    = reshape(cfd,lev,len);
cfd = wcodemat(cfd,nbcol,'row');

figure,hax = [];
h211 = subplot(2,1,1);
h211.XTick = [];
plot(signal,'r');
title('Analyzed signal.');
hax(1) = gca;
hax(1).XLim = [1 length(signal)];
subplot(2,1,2);
colormap%(cool(128));
image(cfd);
tics = 1:lev;
labs = int2str(tics');
hax(2) = gca;
hax(2).YTickLabelMode = 'manual';
hax(2).YDir = 'normal';
hax(2).Box = 'On';
hax(2).YTick = tics;
hax(2).YTickLabel = labs;
title('Discrete Transform, absolute coefficients.');
ylabel('Level');
linkaxes(hax,'x');



figure;
[cfs,f] = cwt(signal,ones(size(signal)));
hp = pcolor(1:length(signal),f,abs(cfs)); hp.EdgeColor = 'none';
set(gca,'YScale','log');
xlabel('Sample'); ylabel('log10(f)');c