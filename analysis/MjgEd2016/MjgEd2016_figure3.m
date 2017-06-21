%% FIG3: Examples of feature dynamics and head body independence ---|
% A  VIDEO snapshots of rat gait                                    |
% B  Examples of limb movements along side spine sway               |
% C  SVD decomposition of movements in horizontal plane relative    |
%        to the direction of the body.                              |
% __________________________________________________________________|


% FIG3FETMAT ---------------------------------------------------------------------------
% parent: figure 3
% subplots:
%    subplot 1: Selected timeperiod of feature matrix
%    subplot 2: Corresponds to subplot 1, contains state labels
% location: MjgEd2016_figure3_svd_walk_alt.m
%




% Figure Settings ----------------------------------------------------------------------
OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'MjgEd2016/figures/Figure_3';
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












