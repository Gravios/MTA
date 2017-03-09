

%function transition(Stc,varargin)

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('transitionWindow',0.2,...
                 'ReferenceData',   Trial.load('xyz'));

[transitionWindow, ReferenceData] = DefaultArgs(varargin,defargs,'--struct')
Trial = MTATrial.validate('jg05-20120317.cof.all');
xyz = Trial.load('xyz');
smat = stc2mat(Stc,xyz,{'walk','rear'});
 = 0.2;%seconds 
indexShift = round(transitionWindow * xyz.sampleRate / 2);

stsTransitionMat = reshape(permute(cat(3,...
                                       circshift(smat,indexShift),...
                                       circshift(smat,-indexShift)),...
                                   [1,2,3]),...
                           size(smat,1),[]);

stsTargetMat     = repmat([1,0,0,1],size(smat,1),1);

stsTransitionPeriods = ThreshCross(all(stsTransitionMat==stsTargetMat,2),0.5,0);

figure,hold on;
plot(all(stsTransitionMat==stsTargetMat,2),'b')

Lines(Trial.stc{'r'}(:),[],'r');

cat(3,repmat(reshape[1,0,0,0,0,0;0,1,0,0,0,0],size(smat,1),1)
% $$$ s = MTASession('jg03-20110428');
% $$$ s = MTASession('jg03-20110429');
% $$$ s = MTASession('jg03-20110430');
% $$$ s = MTASession('jg03-20110501');
% $$$ s = MTASession('gs04-20110921');
%Trial = MTATrial(s);

Trial = MTATrial.validate('jg05-20120317.cof.all');
Trial.load('stc','hand_labeled_rev3_jg');
%Trial.load('stc','hand_labeled_rev1_Ed');


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
fxyz.filter('ButFilter',3,2.4);
vxy = fxyz.vel([],[1,2]);
vxy.data(vxy.data<=1e-3)=1e-3;

% traj NORM body vector projection
shft = 3;
tmar = {'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'};
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
        walkFetRot(nind,t==rotationAngles,m) = nunity(dot(tvec(nind,m,:),unvec(nind,t==rotationAngles,:),3));
    end
end

% DISPLAY spine sway along spine
figure,hold on,
%plot(sqrt(sum((xyz(:,'fsl',[1,2])-xyz(:,'hip_right',[1,2])).^2,3))-40)
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
% svd
bhvPeriods = Trial.stc{'w+n'};
embeddingWindow = 64;
wfet = xyz.copy;
wfet.data= zeros([size(xyz,1),size(walkFetRot,2)*size(walkFetRot,3)]);
wfet.data(nz,:) = [reshape(walkFetRot(nz,:),[],size(walkFetRot,2)*size(walkFetRot,3))];
wfs = wfet.segs([],embeddingWindow);
wfs = circshift(wfs,embeddingWindow/2,2);
wfs = MTADxyz('data',reshape(permute(wfs,[2,1,3]),size(wfs,2),[]),'sampleRate',xyz.sampleRate);
wfs.data(isnan(wfs.data(:)))=0;
[Uw,Sw,Vw] = svd(wfs(bhvPeriods,:),0);
[LU,LR,FSr,VT] = erpPCA(wfs(bhvPeriods,:),10);



% DISPLAY eigen vectors
figure,
for i = 1:10,
    subplot(2,5,i);imagesc(reshape(LR(:,i),[],size(wfet,2))'),
    axis xy;
    colorbar
end
% DISPLAY eigen vectors
figure,
for i = 1:4,
    subplot(6,6,i*6-);imagesc(reshape(Vw(:,i),[],size(wfet,2))'),
    axis xy;
    colorbar
end

fetWlr = MTADxyz('data',wfs.data*LR(:,1),'sampleRate',xyz.sampleRate);
for i = 1:10,fetWlr.data(:,i) = wfs.data*LR(:,i);end

fetW = MTADxyz('data',wfs.data*Vw(:,1),'sampleRate',xyz.sampleRate);
for i = 1:10,fetW.data(:,i) = wfs.data*Vw(:,i);end


sclr = [0,0,1;...
        0,1,0;...
        1,0,0;...
        1,0,1];
sclr = 'bgcr';

ts = [1:size(xyz,1)]./xyz.sampleRate;

figure,sp = [];
sp(end+1)=subplot2(10,4,[1:8],[2:4]);
hold on;
plot(ts,walkFetRot(:,3,1).*5+60), % spine sway 'spine_lower'
plot(ts,walkFetRot(:,3,2).*5+70), % spine sway 'pelvis_root'
plot(ts,walkFetRot(:,3,3).*5+80), % spine sway 'spine_middle'
plot(ts,walkFetRot(:,3,4).*5+90), % spine sway 'spine_upper'
% $$$ plot(ts,-fetWRot(:,1)*10,'b')            % Walk comp
% $$$ plot(ts,fetWRot(:,2)*10,'r')            % Walk comp
% $$$ plot(ts,fetWRot(:,3)*10,'g')            % Turn comp
% $$$ plot(ts,-fetWlr(:,1)/0.5e2,'b')            % Walk comp
% $$$ plot(ts,fetWlr(:,2)/0.5e2,'r')            % Walk comp
% $$$ plot(ts,fetWlr(:,4)/0.5e2,'g')            % Turn comp
plot(ts,fetW(:,1),'b')            % Walk comp
plot(ts,fetW(:,2),'r')            % Walk comp
plot(ts,fetW(:,3),'g')            % Turn comp
plot(ts,fetW(:,5),'m')            % Turn comp
Lines([],0,'k');
Lines([],5,'r');
Lines([],5,'r');
sp(end+1)=subplot2(10,4,[9,10],[2:4]);
plotSTC(Trial.stc,1,'text',{'walk','turn','pause','rear'},sclr,[],false);
linkaxes(sp,'x');
for i = 2:2:10,
sp(end+1)=subplot2(10,4,i-1:i,1);
imagesc(reshape(Vw(:,i/2),[],size(wfet,2))'),
    axis xy;
    colorbar
end

walkOn = fetW.segs(Trial.sts{'w'}(:,1),100);

B = rotatefactors(V, 'Method','varimax');

fetWRot = MTADxyz('data',wfs.data*B(:,1),'sampleRate',xyz.sampleRate);
for i = 1:10,fetWRot.data(:,i) = wfs.data*B(:,i);end


% DISPLAY eigen vectors
figure,
for i = 1:4,
    subplot(1,4,i);imagesc(reshape(B(:,end-i),[],size(wfet,2))'),
    axis xy;
    colorbar
end




for s = 1:numel(steps),
    line(ts(repmat(steps(s),1,2)),repmat(swayMagnitude(s),1,2)+[-2,2],'Color','m');
end

[steps,swayMagnitude] = LocalMinima(-abs(ButFilter(fetW(:,3),3,[1,6]./(xyz.sampleRate/2),'bandpass')),0,-8);
%[steps,swayMagnitude] = LocalMinima(-abs(wfs.data*Vw(:,3)),10,0);


stepsWalk = SelectPeriods(steps,[Trial.stc{'w'}.data],'d',1,0);
stepsNotWalk = SelectPeriods(steps,[Trial.stc{'w'}.data],'d',0,0);
figure,
edx = linspace(-50,50,100);
h = bar(edx,histc(fetW(   stepsWalk,3),edx),'histc');
h.FaceColor = 'r';h.EdgeColor = 'r';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;
hold on
h = bar(edx,histc(fetW(stepsNotWalk,3),edx),'histc');
h.FaceColor = 'c';h.EdgeColor = 'c';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;



figure,hold on
stepsW = swayMagnitude<-8;
plot(stepsW)
stepsW = stepsW...
         |ismember([stepsW,circshift(stepsW,-1)],[0,1],'rows')...
         |ismember([stepsW,circshift(stepsW,1)],[0,1],'rows');
plot(stepsW+.1,'r')

stepsW = [1;steps(find(stepsW));size(xyz,1)];



dstep = zeros([size(xyz,1),1]);
wstep = zeros([size(xyz,1),1]);
tstep = zeros([size(xyz,1),1]);
xstep = zeros([size(xyz,1),1]);
for i = 1:numel(stepsW)-2,
    dstep(stepsW(i)+1:stepsW(i+1)) = sqrt(sum((xyz(stepsW(i+1),1,[1,2])-xyz(stepsW(i),1,[1,2])).^2,3))./(stepsW(i+1)-stepsW(i));
    wstep(stepsW(i)+1:stepsW(i+1))  = sum(walkFetRot(stepsW(i:i+1),1,1))./(stepsW(i+1)-stepsW(i));
    tstep(stepsW(i)+1:stepsW(i+2))  = abs(circ_dist(ang(stepsW(i+1),'fsl','hcom',1),ang(stepsW(i),'fsl','hcom',1)))./(stepsW(i+1)-stepsW(i));
    xstep(stepsW(i)+1:stepsW(i+2))  = abs(circ_dist(ang(stepsW(i+1),'fsl','hcom',1),ang(stepsW(i),'fsl','hcom',1)));
end


dstep = MTADxyz('data',dstep,'sampleRate',xyz.sampleRate);
wstep = MTADxyz('data',wstep,'sampleRate',xyz.sampleRate);
tstep = MTADxyz('data',tstep,'sampleRate',xyz.sampleRate);
xstep = MTADxyz('data',xstep,'sampleRate',xyz.sampleRate);

figure,
hold on,
plot(dstep.data)
plot(wstep.data)
plot(tstep.data)
Lines(Trial.stc{'w'}(:)-.5,[],'b');
Lines(Trial.stc{'n'}(:),[],'g');
Lines(Trial.stc{'m'}(:),[],'m');
Lines(Trial.stc{'r'}(:)-.5,[],'r');


bfet = MTADxyz('data',wfs.data(:,32:64:size(wfs,2))*Vw(32:64:size(wfs,2),1),'sampleRate',xyz.sampleRate);
tfet = MTADxyz('data',var(walkFetRot(:,3,:),[],3),'sampleRate',xyz.sampleRate);
lfet = MTADxyz('data',abs(prod(walkFetRot(:,1,:),3)),'sampleRate',xyz.sampleRate);
lfet.data(lfet.data<1e-5) = 1e-6;
tprinc = reshape(Vw(:,i),[],size(wfet,2));
xfet = MTADxyz('data',sum(bsxfun(@times,...
                                 reshape(walkFetRot,size(wfet)),tprinc(25,:)),2),...
               'sampleRate',xyz.sampleRate);
xfet.data = abs(xfet.data-3);
xfet.data(lfet.data<1e-5) = 1e-6;

dstep.filter('ButFilter',3,2.5,'low');
wstep.filter('ButFilter',3,2.5,'low');
tstep.filter('ButFilter',3,2.5,'low');
xstep.filter('ButFilter',3,2.5,'low');

fet = wstep.copy;
fet = tstep.copy;
fet.data = log10(fet.data);
edx = linspace(-4,1,100);
edx = linspace(-6,-1,100);
%edx = linspace(-.1,.25,100);
figure,hold on
ind = Trial.stc{'a'};
h = bar(edx,histc(fet(ind),edx),'histc');
h.FaceColor = [0.1,0.1,0.1];h.EdgeColor = [0.1,0.1,0.1];h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;

ind = Trial.stc{'w'};
h = bar(edx,histc(fet(ind),edx),'histc');
h.FaceColor = 'r';h.EdgeColor = 'r';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;

ind = Trial.stc{'n'};
h = bar(edx,histc(fet(ind),edx),'histc');
h.FaceColor = 'g';h.EdgeColor = 'g';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;

ind = Trial.stc{'a-w-n-r'};
h = bar(edx,histc(fet(ind),edx),'histc');
h.FaceColor = 'c';h.EdgeColor = 'c';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;

ind = Trial.stc{'p'};
h = bar(edx,histc(fet(ind),edx),'histc');
h.FaceColor = 'm';h.EdgeColor = 'm';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;



lfet = MTADxyz('data',abs(prod(walkFetRot(:,1,:),3)./var(walkFetRot(:,1,:),[],3)),'sampleRate',xyz.sampleRate);
lfet.data(lfet.data<1e-10) = 1e-10;
vfet = MTADxyz('data',abs(prod(walkFetRot(:,1,:),3)),'sampleRate',xyz.sampleRate);
vfet.data(vfet.data<1e-10) = 1e-10;

figure
states = {'a','s','m','p','n','w','r'};
for s = 1:numel(states)
subplot(3,3,s)
ind = Trial.stc{states{s}};
%hist2([log10(abs(fetWlr(ind,2))),log10(abs(fetWlr(ind,4)+200))],linspace(0,4,100),linspace(0,4,100));
%hist2([log10(tstep(ind)),log10(abs(fetWlr(ind,4)+200))],linspace(-6,-1,100),linspace(0,4,100));
%hist2([log10(dstep(ind)),log10(tstep(ind))],linspace(-4,1,100),linspace(-6,-1,100))
%hist2([log10(dstep(ind)),log10(tstep(ind))],linspace(-4,1,50),linspace(-6,-1,50))
%hist2([log10(bfet(ind)),log10(wstep(ind))],linspace(-4,1,100),linspace(-4,1,100))
%hist2([log10(vxy(ind,1)),log10(lfet(ind))],linspace(-3,2,100),linspace(-9,4,100))
%hist2([log10(dstep(ind,1)),log10(lfet(ind))],linspace(-4,1,100),linspace(-9,4,100))
%hist2([log10(wstep(ind,1)),log10(vxy(ind,1))],linspace(-6,0,100),linspace(-2.9,2,100))
%hist2([log10(wstep(ind,1)),log10(tfet(ind,1))],linspace(-6,0,100),linspace(-4,2,100))
%hist2([log10(tstep(ind,1)),log10(lfet(ind))],linspace(-6,-1,100),linspace(-9,4,100))
%hist2([log10(tfet(ind)),log10(lfet(ind))],linspace(-4,2,100),linspace(-9,4,100))
%hist2([log10(vxy(ind,1)),log10(lfet(ind))],linspace(-3,2,100),linspace(-9,4,100))
%hist2([log10(dstep(ind,1)),log10(lfet(ind))],linspace(-4,1,100),linspace(-9,4,100))
%hist2([log10(dstep(ind)),log10(tfet(ind))],linspace(-4,1,100),linspace(-4,2,100))
hist2([log10(vxy(ind,'spine_lower')),log10(tfet(ind))],linspace(-2.9,2,100),linspace(-4,2,100))
%hist2([log10(vxy(ind,'spine_lower')),log10(tstep(ind))],linspace(-2.9,2,100),linspace(-6,-1,100))
%hist2([log10(vxy(ind,1)),log10(tstep(ind))],linspace(-2.9,2,100),linspace(-6,-1,100))
%hist2([log10(xfet(ind)),log10(wstep(ind))],linspace(.3,.8,100),linspace(-6,0,100))
%hist2([log10(vxy(ind)),log10(wstep(ind))],linspace(-2.9,2,50),linspace(-6,0,50))
title(ind.label)
caxis([0,150])
grid on
set(gca,'GridColor',[1,1,1])
end






%plot(ButFilter(walkFetRot(:,3,2),3,[0.5,6]./(xyz.sampleRate/2),'bandpass'));
Lines([],1,'k');
Lines([],-1,'k');
Lines(Trial.stc{'w'}(:)-.5,[],'b');
Lines(Trial.stc{'n'}(:),[],'g');
Lines(Trial.stc{'m'}(:),[],'m');


Lines(Trial.stc{'w'}(:)-.5,[],'b');
Lines(Trial.stc{'n'}(:),[],'g');
Lines(Trial.stc{'m'}(:),[],'m');

figure,hold on
ind = Trial.stc{'w'};
plot(fetW(ind,1),fetW(ind,2),'.b');
ind = Trial.stc{'n'};
plot(fetW(ind,1),fetW(ind,2),'.g');
ind = Trial.stc{'p'};
plot(fetW(ind,1),fetW(ind,2),'.c');




figure,hold on,whitebg('k')
plot(walkFet,'w-')
plot(walkFet(:,5),'g-')
plot(walkFet(:,4),'b-')
plot(walkFet(:,2),'r-')
plot(walkFetBody+10,'w-')
plot(walkFetBody(:,5)+10,'g-')
plot(walkFetBody(:,4)+10,'b-')
plot(walkFetBody(:,2)+10,'r-')
Lines([],10,'w');
Lines(Trial.stc{'w'}(:)-.5,[],'w');
Lines(Trial.stc{'n'}(:),[],'g');

figure,hold on
plot(xyz(:,'spine_lower',3)-xyz(:,'fsl',3))
Lines(Trial.stc{'w'}(:)-.5,[],'m');
Lines(Trial.stc{'n'}(:),[],'g');


parspec = struct('nFFT',2^9,...
                 'Fs',  Trial.xyz.sampleRate,...
                 'WinLength',2^8,...
                 'nOverlap',2^8*.5,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[.1,6]);

sessionList = get_session_list('hand_labeled');
for s = sessionList
    Trial = MTATrial.validate(s);
    [Feature,fs] = fet_hips(Trial,'xyz','mtchglong',false,parspec);
    behavior_feature_spectrum_binned_by_second_feature(Trial,Feature,fs,[],'walk',s.stcMode);
end









Data = xyz.copy;
%Data.data = walkFetFilt;
Data.data = [sq(walkFetRot(:,3,:))];
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
                 'FreqRange',[0.1,20]);
[ys,fs,ts] = fet_spec(Trial,Data,'mtcsdglong',true,[],parspec,true);

%[ys,fs,ts] = fet_spec(Trial,efet,'mtcsdglong',true,[],parspec,true);
rs = ys.copy();
rs.clear();
for i = 1:size(ys,3)
for j = 1:size(ys,4)
    rs.data(:,:,i,j) = ys(:,:,i,j)./(sqrt(ys(:,:,i,i)).*sqrt(ys(:,:,j,j)));
end
end


figure,sp = [];
sp(end+1)=subplot(211);imagesc(ts,fs,log10(ys(:,:,1,1))');axis xy;colormap jet;
sp(end+1)=subplot(212);imagesc(ts,fs,log10(ys(:,:,7,7))');axis xy;colormap jet;
linkaxes(sp,'x');

fxyz = xyz.copy();
fxyz.filter('ButFilter',3,2.4,'low');
vxy = fxyz.vel([],[1,2]);
vxy.resample(ys);
vxy.data(vxy.data<=1e-3) = 1e-3;
vxy.data(nniz(vxy),:)=log10(vxy(nniz(vxy),:));


m1 = 2;
m2 = 10;
f1 = 6;
f2 = 6;
fs(f1)
fs(f2)

hsp = rs.copy();
hsp.data = atan(imag(rs(:,f1,m1,m2))./real(rs(:,f1,m1,m2)));

mm1 = 1;
mm2 = 4;
hsp2 = rs.copy();
hsp2.data = atan(imag(rs(:,f2,mm1,mm2))./real(rs(:,f2,mm1,mm2)));


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
%hist2([log10(ys(ind,f2,m1,m1)),log10(ys(ind,f2,m2,m2))],linspace(-10,1,30),linspace(-10,1,30))
%hist2([log10(ys(ind,f2,m1,m1)),hsp(ind)],linspace(-10,-1,30),linspace(-pi/2,pi/2,30))
%hist2([log10(ys(ind,f2,m1,m1)),vxy(ind,1)],linspace(-10,-1,30),linspace(-2.9,2,30))
hist2([vxy(ind,1),hsp2(ind)],linspace(-2.9,2,35),linspace(-pi/2,pi/2,35))
%hist2([hsp(ind),hsp2(ind)],linspace(-pi/2,pi/2,35),linspace(-pi/2,pi/2,35))
title(ind.label)
caxis([0,50])
grid on
set(gca,'GridColor',[1,1,1])
end



figure,plot(hsp)
Lines(Trial.stc{'w',rs.sampleRate}(:),[],'m');
Lines(Trial.stc{'n',rs.sampleRate}(:),[],'g');


figure,hold on
plot(Data(:,[1,7]))
Lines(Trial.stc{'w'}(:)-.5,[],'m');
Lines(Trial.stc{'n'}(:),[],'g');




pfet = Data.phase([1,6]);

sclr = [0,0,1;...
        1,0,1;...
        1,0,0;...
        0,1,0];

figure,sp = [];
sp(end+1)=subplot2(8,1,[1:7],1);imagesc(pfet.data');
sp(end+1)=subplot2(8,1,8,1);plotSTC(Trial.stc,pfet.sampleRate,[],{'walk','pause','rear','turn'},sclr,[],false);
linkaxes(sp,'x');


wper = Trial.stc{'w'};
wdur = diff(wper.data,1,2);
phdmat = zeros([numel(wdur),max(wdur)]);
[~,wporder] = sort(wdur);

m1 = 1;
m2 = 2;
for w = 1:size(wdur,1)
    twp = wper(wporder(w),:);
    phdmat(w,1:wdur(wporder(w))+1) = circ_dist(pfet(twp,m1),pfet(twp,m2));
end
figure,imagesc(phdmat),colormap hsv

sxyz = xyz.copy();
sxyz.filter('ButFilter',3,40);

mvxy = sxyz.vel([],[1,2]);
fvxy = mvxy.copy();
fvxy.filter('ButFilter',3,1);

figure,hold on,
plot(nunity(mvxy(:,1)-fvxy(:,1)).*3,'w')
plot(Data(:,2),'c')
plot(nunity(ang(:,'spine_lower','fbcom',3)),'r');
Lines(Trial.stc{'w'}(:)-.5,[],'m');
Lines(Trial.stc{'n'}(:),[],'g');

figure,hold on,
plot(nunity(ang(:,'spine_lower','fbcom',3)));
plot(nunity(ang(:,'fsl','bcom',3)));
Lines(Trial.stc{'w'}(:)-.5,[],'m');
Lines(Trial.stc{'n'}(:),[],'g');




Data.data = Data(:,[1,2,4,5]);
dsegs = Data.segs([],20);

figure,hold on
plot(walkFet(:,[1,2,4,5]),'w-')
plot(var(Data.data,[],2)*2+5)
plot(var(reshape(permute(dsegs,[2,1,3]),size(Data,1),[]),[],2)*2+2.5)
plot(nunity(circ_dist(circshift(fang(:,'spine_lower','spine_upper',1),-3),circshift(fang(:,'spine_lower','spine_upper',1),3)))+10,'g');
Lines([],10,'w');
Lines(Trial.stc{'w'}(:)-.5,[],'b');
Lines(Trial.stc{'n'}(:),[],'g');


[U,S,V] = svd(pfet(Trial.stc{'w'},:),0);
[Ut,St,Vt] = svd(pfet(Trial.stc{'n'},:),0);


wper = Trial.stc{'w'};
wdur = diff(wper.data,1,2);


figure,subplot(221),plot(diag(S)),subplot(222);plot(U(:,1),U(:,2),'.');subplot(223);imagesc(V');

figure,hold on
dim1 = 4;
dim2 = 5;
plot(U(:,dim1),U(:,dim2),'.');
shiftInds = 1;
shiftInds = [];


shiftInds = find(wdur>120,numel(shiftInds)+1,'first')+1;
shiftInd = shiftInds(end);
hold on,plot(U(sum(wdur(1:shiftInd)+1)+1:sum(wdur(1:shiftInd+1)+1),dim1),...
             U(sum(wdur(1:shiftInd)+1)+1:sum(wdur(1:shiftInd+1)+1),dim2),'g')

figure,hold on,
shiftInd = 4;
ind = sum(wdur(1:shiftInd)):sum(wdur(1:shiftInd+1));
plot(U(ind,1).*100)
plot(U(ind,2).*100)
plot(walkFet(ind,:))


wvxy =vxy(wper,1);

c = jet(25);
eds = 0:2.5:50;
[~,cind] = histc(wvxy(:),eds);
cind(cind==0) = 1;

dim1 = 2;
dim2 = 3;
figure,
plot(U(:,dim1),U(:,dim2),'.')


wper = Trial.stc{'w'};
wdur = diff(wper.data,1,2)+1;


figure,hold on
shiftInds = 1;
shiftInds = [];
plot(U(:,dim1),U(:,dim2),'.')
wdurInds = find(wdur>180);
i = 1;
for shiftInd = wdurInds'
    shiftInd = wdurInds(i);
    ind = sum(wdur(1:shiftInd-1)):sum(wdur(1:shiftInd));
    plot(U(ind,dim1),U(ind,dim2),'g-')
    i = i+1;
 % $$$     hold on,scatter3(U(ind,dim1),...
% $$$                   U(ind,dim2),...
% $$$                   wvxy(ind),10,c(cind(ind),:),'filled')
end
%zlim([0,100])

figure,plot3(U(:,2),U(:,3),vxy(wper,1),'.')

figure,plot()
hold on,plot(efs.data*V(:,3))


figure,
for i = 1:5,
    subplot(1,5,i);imagesc(reshape(Vt(:,i),[],18)'), ...
end

fet1 = MTADxyz('data',efs.data*Vt(:,1),'sampleRate',xyz.sampleRate);
fet2 = MTADxyz('data',efs.data*Vt(:,2),'sampleRate',xyz.sampleRate);
figure,hold on,plot(fet1.data),plot(fet2.data),
Lines(Trial.stc{'w'}(:)-.5,[],'b');
Lines(Trial.stc{'n'}(:),[],'g');


wdurInds = find(wdur>180);
wperLong = wper(wdurInds,:);

figure,hold on
plot(fet1(wper),fet2(wper),'.')
for w = wperLong'
    ind = w(1):w(2);
    plot(fet1(ind),fet2(ind),'g-')
 % $$$     hold on,scatter3(U(ind,dim1),...
% $$$                   U(ind,dim2),...
% $$$                   wvxy(ind),10,c(cind(ind),:),'filled')
end





figure,plot(fetW(:,[9,10]))
Lines(Trial.stc{'w'}(:)-.5,[],'b');
Lines(Trial.stc{'n'}(:),[],'g');



fet1.filter('ButFilter',3,1,'high');
fet2.filter('ButFilter',3,1,'high');
figure,hold on,plot(fet1.data),plot(fet2.data),
Lines(Trial.stc{'w'}(:)-.5,[],'b');
Lines(Trial.stc{'n'}(:),[],'g');
Lines([],5,'k'),Lines([],-5,'k')
plot(fetb.data),
plot(fetC.data),


fxyz = xyz.copy();
fxyz.filter('ButFilter',3,2.4,'low');
vxy = fxyz.vel([],[1,2]);
vxy.data(vxy.data<=1e-3) = 1e-3;
vxy.data(nniz(vxy),:)=log10(vxy(nniz(vxy),:));


figure
states = {'a','s','m','p','n','w','r'};
for s = 1:numel(states)
subplot(3,3,s)
ind = Trial.stc{states{s}};
hist2([vxy(ind,1),log10(fetb(ind))],linspace(-2.9,2,35),linspace(-2.9,2,35))
title(ind.label)
caxis([0,500])
grid on
set(gca,'GridColor',[1,1,1])
end





nz = nniz(xyz);
% svd
bhvPeriods = Trial.stc{'r'};
embeddingWindow = 64;
rfet = xyz.copy;
rfet.data= zeros([size(xyz,1),numel(tmar)]);
rfet.data(nz,:) = xyz(nz,tmar,3);
rfs = rfet.segs([],embeddingWindow);
rfs = circshift(rfs,embeddingWindow/2,2);
rfs = reshape(permute(rfs,[2,1,3]),size(rfs,2),[]);
rfs = MTADxyz('data',rfs,'sampleRate',xyz.sampleRate);
[Ur,Sr,Vr] = svd(rfs(bhvPeriods,:),0);

figure,
subplot(221),plot(diag(Sr)),
subplot(222);plot(Ur(:,1),Ur(:,2),'.');
subplot(223);imagesc(Vr');

figure,
for i = 1:10,
    subplot(2,5,i);imagesc(reshape(Vr(:,i),[],size(rfet,2))');
end




fetR1 = MTADxyz('data',rfs.data*Vr(:,1),'sampleRate',xyz.sampleRate);
fetR2 = MTADxyz('data',rfs.data*Vr(:,2),'sampleRate',xyz.sampleRate);
fetR3 = MTADxyz('data',rfs.data*Vr(:,3),'sampleRate',xyz.sampleRate);
fetRC = MTADxyz('data',max(abs([fetR1.data,fetR2.data]),[],2),'sampleRate',xyz.sampleRate);



figure,hold on,plot(fetR1.data),plot(fetR2.data),
Lines(Trial.stc{'w'}(:)-.5,[],'b');
Lines(Trial.stc{'n'}(:),[],'g');
Lines([],5,'k'),Lines([],-5,'k')


figure,hold on,
plot(fetR1.data);
plot(fetR2.data);
plot(fetR3.data);
Lines(Trial.stc{'r'}(:),[],'r');
Lines(Trial.stc{'w'}(:),[],'b');
