%% FIG4: Examples of feature dynamics and head body independence ---|
% A  VIDEO snapshots of rat gait                                    |
% B  Examples of limb movements along side spine sway               |
% C   |
% __________________________________________________________________|




% Figure Settings ----------------------------------------------------------------------
OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'MjgEd2016/manuscript/Figures/Figure_3';
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)
% --------------------------------------------------------------------------------------




% LOAD Session with Side Video Tracking
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
        walkFetRot(nind,t==rotationAngles,m) = nunity(dot(tvec(nind,m,:),unvec(nind,t==rotationAngles,:),3));
    end
end




%% Fig.3.A Video snapshots 
% LOAD videos 
videofile = ['/storage/share/exchange/gravio/optitrack/Session\ 2017-01-12/Take\ 2017-01-12\ 02.52.14\ PM-Camera\ 13029.avi'];
videofile = '/storage/gravio/data/project/general/JM11-20170112/Trial001_C13029.avi';
videofile_side = '/storage/gravio/data/project/general/JM11-20170112/Trial001_C12881.avi';
vdr = VideoReader(videofile);
vdr_side = VideoReader(videofile_side);

% Plot Timepoints
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

%% Fig.3.B Spine Sway, Step Detection 
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


xlabel('Time (s)');
FigName = 'walk_feature_demo_';

print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));




%% Fig.3.C SVD of embedded marker trajectories relative to body
s = 0;
sessionList = get_session_list('hand_labeled');



s = s+1;
Trial = MTATrial.validate(sessionList(s));


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

%ang = create(MTADang,Trial,xyz);
%fxyz = xyz.copy();
%fxyz.filter('ButFilter',3,2.4);
%vxy = fxyz.vel([],[1,2]);
%vxy.data(vxy.data<=1e-3)=1e-3;

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


% SVD on XY translational vectors
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
wts = (1:embeddingWindow)./wfet.sampleRate;



% $$$ % DISPLAY eigen vectors
% $$$ hfig = figure;
% $$$ hfig.Units = 'centimeters';
% $$$ hfig.Position(3:4) = [30,4];
% $$$ hfig.PaperPositionMode = 'auto';
% $$$ for i = 1:10,
% $$$     subplot(1,10,i);imagesc(wts,1:size(wfet,2),reshape(Vw(:,i),[],size(wfet,2))'),
% $$$     caxis([-0.08,0.08]);
% $$$     axis xy
% $$$ end
% $$$ FigName = ['SVD_walkturn_',Trial.filebase];
% $$$ print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
% $$$ print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));




% COMPUTE timeseries score for first 10 eigenvectors
fetW = MTADxyz('data',wfs.data*Vw(:,1),'sampleRate',xyz.sampleRate);
for i = 1:10,fetW.data(:,i) = wfs.data*Vw(:,i);end

% $$$ 
% $$$ 
% $$$ sclr = 'bgcr';
% $$$ ts = [1:size(xyz,1)]./xyz.sampleRate;
% $$$ hfig = figure;
% $$$ hfig.Units = 'centimeters';
% $$$ hfig.Position(3:4) = [25,16];
% $$$ hfig.PaperPositionMode = 'auto';
% $$$ sp = [];
% $$$ sp(end+1)=subplot2(10,4,[1:8],[2:4]);
% $$$ hold on;
% $$$ plot(ts,walkFetRot(:,3,1).*5+60,'LineWidth',1); % spine sway 'spine_lower'
% $$$ plot(ts,walkFetRot(:,3,2).*5+70,'LineWidth',1); % spine sway 'pelvis_root'
% $$$ plot(ts,walkFetRot(:,3,3).*5+80,'LineWidth',1); % spine sway 'spine_middle'
% $$$ plot(ts,walkFetRot(:,3,4).*5+90,'LineWidth',1); % spine sway 'spine_upper'
% $$$ % $$$ plot(ts,-fetWRot(:,1)*10,'b')             % Walk comp
% $$$ % $$$ plot(ts,fetWRot(:,2)*10,'r')              % Walk comp
% $$$ % $$$ plot(ts,fetWRot(:,3)*10,'g')              % Turn comp
% $$$ % $$$ plot(ts,-fetWlr(:,1)/0.5e2,'b')           % Walk comp
% $$$ % $$$ plot(ts,fetWlr(:,2)/0.5e2,'r')            % Walk comp
% $$$ % $$$ plot(ts,fetWlr(:,4)/0.5e2,'g')            % Turn comp
% $$$ plot(ts,fetW(:,1),'b','LineWidth',1);           % Walk comp
% $$$ plot(ts,fetW(:,2),'r','LineWidth',1);           % Walk comp
% $$$ plot(ts,fetW(:,3),'y','LineWidth',1);           % Turn comp
% $$$ plot(ts,fetW(:,4),'g','LineWidth',1);           % Turn comp
% $$$ plot(ts,fetW(:,5),'m','LineWidth',1);           % Turn comp
% $$$ legend({'spine lower','pelvis root','spine middle','spine upper','PC1','PC2','PC3','PC4','PC5'})
% $$$ Lines([],0,'k');
% $$$ Lines([],5,'r');
% $$$ Lines([],-5,'r');
% $$$ sp(end+1)=subplot2(10,4,[9,10],[2:4]);
% $$$ plotSTC(Trial.stc,1,'text',{'walk','turn','pause','rear'},sclr,[],false);
% $$$ linkaxes(sp,'x');
% $$$ for i = 2:2:10,
% $$$ sp(end+1)=subplot2(10,4,i-1:i,1);
% $$$ imagesc(wts,1:size(wfet,2),reshape(Vw(:,i/2),[],size(wfet,2))');
% $$$     axis xy;
% $$$     caxis([-0.08,0.08]);
% $$$     ylabel(['PC',num2str(i/2)])
% $$$ end
% $$$ 
% $$$ 
% $$$ FigName = ['SVD_walkturn_timeseries_',Trial.filebase];
% $$$ print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
% $$$ print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));





% REDUCED pc 
rVw = Vw;
for i = 1:16,
    rVw(i:embeddingWindow:end,:) = 0;
end
for i = 48:64,
    rVw(i:embeddingWindow:end,:) = 0;
end

rfetW = MTADxyz('data',wfs.data*rVw(:,1),'sampleRate',xyz.sampleRate);
for i = 1:10,rfetW.data(:,i) = wfs.data*rVw(:,i);end


j = 3;
pcLabels = {'PC1','PC2','PC3','PC4','PC5'};
sclr = 'bgcr';
ts = [1:size(xyz,1)]./xyz.sampleRate;
hfig = figure;
hfig.Units = 'centimeters';
hfig.Position(3:4) = [25,16];
hfig.PaperPositionMode = 'auto';
sp = [];
sp(end+1)=subplot2(10,4,[1:8],[2:4]);
hold on;
plot(ts,walkFetRot(:,3,1).*5+60,'LineWidth',1); % spine sway 'spine_lower'
plot(ts,walkFetRot(:,3,2).*5+70,'LineWidth',1); % spine sway 'pelvis_root'
plot(ts,walkFetRot(:,3,3).*5+80,'LineWidth',1); % spine sway 'spine_middle'
plot(ts,walkFetRot(:,3,4).*5+90,'LineWidth',1); % spine sway 'spine_upper'
% $$$ plot(ts,fetW(:,1),'b','LineWidth',1);           % Walk comp
% $$$ plot(ts,fetW(:,2),'r','LineWidth',1);           % Walk comp
% $$$ plot(ts,fetW(:,3),'y','LineWidth',1);           % Turn comp
% $$$ plot(ts,fetW(:,j),'r','LineWidth',1);           % Turn comp
% $$$ plot(ts,rfetW(:,j),'m','LineWidth',1);           % Turn comp
plot(ts,fetW(:,j),'g','LineWidth',1);           % Turn comp
plot(ts,rfetW(:,j),'b','LineWidth',1);           % Turn comp
% $$$ plot(ts,fetW(:,4),'g','LineWidth',1);           % Turn comp
% $$$ plot(ts,fetW(:,5),'m','LineWidth',1);           % Turn comp
legend({'spine lower','pelvis root','spine middle','spine upper',pcLabels{j},[pcLabels{j},'R']})
Lines([],0,'k');
Lines([],5,'r');
Lines([],-5,'r');
sp(end+1)=subplot2(10,4,[9,10],[2:4]);
plotSTC(Trial.stc,1,'text',{'walk','turn','pause','rear'},sclr,[],false);
linkaxes(sp,'x');
for i = 2:2:10,
sp(end+1)=subplot2(10,4,i-1:i,1);
imagesc(wts,1:size(wfet,2),reshape(rVw(:,i/2),[],size(wfet,2))');
    axis xy;
    caxis([-0.08,0.08]);
    ylabel(['PC',num2str(i/2)])
end


FigName = ['SVD_walkturn_timeseries_',pcLabels{j},'R_MID_',Trial.filebase];
print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));





pcLabels = {'PC1','PC2','PC3','PC4','PC5'};
sclr = 'bgcr';
ts = [1:size(xyz,1)]./xyz.sampleRate;
hfig = figure;
hfig.Units = 'centimeters';
hfig.Position(3:4) = [25,16];
hfig.PaperPositionMode = 'auto';
sp = [];

for i = 2:2:10,
    subplot2(10,4,i-1:i,1);
    imagesc(wts,1:size(wfet,2),reshape(rVw(:,i/2),[],size(wfet,2))');
    axis xy;
    caxis([-0.08,0.08]);
    ylabel(['PC',num2str(i/2)])
    if i ~= 10,
        sp(end+1)=subplot2(10,4,i-1:i,[2:4]);    hold on
        plot(ts,fetW(:,i/2),'g','LineWidth',1);           % Turn comp
        plot(ts,rfetW(:,i/2),'b','LineWidth',1);           % Turn comp
        hold on;
        Lines([],0,'k');
        Lines([],5,'r');
        Lines([],-5,'r');

        legend({pcLabels{i/2},[pcLabels{i/2},'R']})
    end
end

sp(end+1)=subplot2(10,4,[9,10],[2:4]);
plotSTC(Trial.stc,1,'text',{'walk','turn','pause','rear'},sclr,[],false);

linkaxes(sp,'x');

FigName = ['SVD_walkturn_timeseries_Reduced_MID_',Trial.filebase];
print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

