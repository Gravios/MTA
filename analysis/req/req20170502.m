

% Figure Settings ----------------------------------------------------------------------
OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'MjgEd2016/figures/Figure_3';
% --------------------------------------------------------------------------------------



sessionListTag  = 'hand_labeled';%'BHV_S4H5';
sessionList     = get_session_list(sessionListTag);
numSessions     = numel(sessionList);
%sampleRate      = 119.881035;
%embeddingWindow = 64;
%rotationAngles  = [0,45,90,135];




% START data processing ------------------------------------------------------------------
sampleRate = repmat({119.881035},1,numSessions);
embeddingWindow = repmat({64},1,numSessions);
rotationAngles = repmat({deg2rad([0,90])},1,numSessions);
%stcMode = 'NN0317';

% LOAD Sessions and basic data
Trials = af(@(Trial) MTATrial.validate(Trial), sessionList);
Stc    = cf(@(Trial) Trial.load('stc')       , Trials);
xyz    = cf(@(Trial) Trial.load('xyz')       , Trials);
% RESAMPLE to common sampleRate
cf(@(x,s) x.resample(s),xyz,sampleRate);

% ADD virtual markers for body
rb    = cf(@(x)   x.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'}), xyz);
bcom  = cf(@(x,r) x.com(r), xyz,rb);
cf(@(x,b) x.addMarker('fbcom','data',ButFilter(b,4,[0.5]./(x.sampleRate/2),'low')),xyz,bcom);
cf(@(x,b) x.addMarker('bcom','data',b), xyz,bcom);
cf(@(x,b) x.addMarker('fsl','data',ButFilter(x(:,'spine_lower',:),4,[0.5]./(x.sampleRate/2),'low')),xyz,bcom);
clear('bcom');
% ADD virtual markers for head
rb    = cf(@(x)   x.model.rb({'head_back','head_left','head_front','head_right'}),xyz);
hcom  = cf(@(x,r) x.com(r),xyz,rb);
cf(@(x,h) x.addMarker('fhcom','data',ButFilter(h,4,[0.8]./(x.sampleRate/2),'low')),xyz,hcom);
cf(@(x,h) x.addMarker('hcom','data',h), xyz,hcom);
clear('hcom');
% ADD virtual markers for all markers
rb    = cf(@(x)   x.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper',...
                              'head_back','head_front'}),xyz);
acom  = cf(@(x,r) x.com(r),xyz,rb);
cf(@(x,h) x.addMarker('facom','data',ButFilter(h,4,[0.8]./(x.sampleRate/2),'low')),xyz,acom);
cf(@(x,h) x.addMarker('acom','data',h), xyz,acom);
clear('acom');


% COMPUTE Features ------------------------------------------------------------------------------
% @walkFetRot
% NORM body vector projection
tmar = repmat({{'spine_lower','pelvis_root','spine_middle','spine_upper','hcom','acom','bcom'}},1,numSessions);
mxyz = cf(@(x) x.copy('empty'),xyz);
cf(@(m,x,t) set(m,'data',x(:,t(:),:)),mxyz,xyz,tmar);
tvec = cf(@(x) markerDiffMatrix(x), mxyz);
tvec = cf(@(x,a) repmat(reshape(x,[size(x,1),size(x,3).^2,size(x,4)]),[1,1,1,numel(a)]),tvec,rotationAngles);

% BODY vector
bodyUnitVector = cf(@(x,a) repmat(x(:,'fbcom',[1,2])-x(:,'fsl',[1,2]),1,1,1,numel(a)),xyz,rotationAngles);
rotMat = cf(@(t,m) reshape(repmat(permute([cos(t),-sin(t);...
                                           sin(t),cos(t)],...
                                          [3,1,2]),...
                                  [size(m,1),1,1]),...
                           size(m,1),2,2,numel(t)),rotationAngles,bodyUnitVector);
unvec  = cf(@(m,r,t) repmat(bsxfun(@rdivide,multiprod(m,r,[2,3],[2,3]),sqrt(sum(m.^2,3))),...
                            [1,numel(t)^2,1,1]),bodyUnitVector,rotMat,tmar);
% BODY oriented trajectory translation vector
walkFetRot = cf(@(t,u)  cat(3,sq(dot(t(:,:,[1,2],:),u,3)),t(:,:,3,1)),tvec,unvec);



figure;
hold('on');
s = 3;
for i = 2:5,
    plot(walkFetRot{s}(:,sub2ind([7,7],1,i),3));
end


figure;
hold('on');
s = 3;
for i = 2:5,
    plot(walkFetRot{s}(:,sub2ind([7,7],1,i),2));
end


figure;
hold('on');
s = 3;
m = 2;
for i = setdiff(1:5,m)
    plot(walkFetRot{s}(:,sub2ind([7,7],m,i),1));
end


wfrCat = cat(1,walkFetRot{:});
wfrMean = nanmean(wfrCat(nniz(wfrCat),:,:));
wfrStd = nanstd(wfrCat(nniz(wfrCat),:,:));
walkFetRot = cf(@(w,m,s) nunity(w,[],m,s),...
                walkFetRot,...
                repmat({wfrMean},1,numSessions),...
                repmat({wfrStd},1,numSessions));
clear('wfrCat','wfrMean','wfrStd')




% Set up Features
fets = cf(@(T) MTADfet(T.spath,...
              [],...
              [],...
              newSampleRate,...
              T.sync.copy,...
              T.sync.data(1),...
              [],'TimeSeries',[],'feature_req20170502','fet_rear_desc','m'),...
          Trials);

