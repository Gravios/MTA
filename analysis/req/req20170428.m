% req20170428 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: secondary exploration of marker movements in space
%               relative to body. Based on work in MjgEd2016_Figure3.m
%               Req 1.0: add z axis movements to svd analysis 
%               Req 1.1: stats between expert labeler classifier (ELC) labels 
%               Req 1.2: stats between expert labeler classifier (ELC) labels 
%                        and neural network classifier (NNC) labels 
%  Bugs: NA
 


% Figure Settings ----------------------------------------------------------------------
OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'MjgEd2016/figures/Figure_3';
% --------------------------------------------------------------------------------------


%% Fig.3.C SVD of embedded marker trajectories relative to body %%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis of horizontal motion relative of to the body
%
%

sessionListTag  = 'hand_labeled';%'BHV_S4H5';
sessionList     = get_session_list(sessionListTag);
numSessions     = numel(sessionList);
%sampleRate      = 119.881035;
%embeddingWindow = 64;
%rotationAngles  = [0,45,90,135];



% START data processing ------------------------------------------------------------------
sampleRate = repmat({119.881035},1,numSessions);
embeddingWindow = repmat({64},1,numSessions);
rotationAngles = repmat({deg2rad([0,45,90,135])},1,numSessions);
%stcMode = 'NN0317';

% LOAD Sessions and basic data
Trials = af(@(Trial) MTATrial.validate(Trial), sessionList);
Stc    = cf(@(Trial) Trial.load('stc')       , Trials);
xyz    = cf(@(Trial) Trial.load('xyz')       , Trials);
% RESAMPLE to common sampleRate
cf(@(x,s) x.resample(s),xyz,sampleRate);

% ADD virtual markers for body
rb    = cf(@(x)   x.model.rb({'spine_lower','pelvis_root','spine_middle'}), xyz);
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


% COMPUTE Features ----------------------------------------------------------------------------
% @walkFetRot
% NORM body vector projection
tmar = repmat({{'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'}},1,numSessions);
tvec = cf(@(x,m,t) repmat(circshift(x(:,m,[1,2]),-(round(x.sampleRate*0.05)/2))...
                          -circshift(x(:,m,[1,2]), (round(x.sampleRate*0.05)/2)),...
                          [1,1,1,numel(rotationAngles{1})]),...
          xyz,tmar,rotationAngles);
% BODY vector
mvec = cf(@(x)  repmat(x(:,'fbcom',[1,2])-x(:,'fsl',[1,2]),1,1,1,4),xyz);
%mvec = cf(@(x)  repmat(x(:,'spine_upper',[1,2])-x(:,'spine_lower',[1,2]),1,1,1,4),xyz);
rotMat = cf(@(t,m) reshape(repmat(permute([cos(t),-sin(t);...
                                           sin(t),cos(t)],...
                                          [3,1,2]),...
                                  [size(m,1),1,1]),...
                           size(m,1),2,2,numel(t)),rotationAngles,mvec);
unvec  = cf(@(m,r,t) repmat(bsxfun(@rdivide,multiprod(m,r,[2,3],[2,3]),sqrt(sum(m.^2,3))),...
                            [size(t),1,1,1]),mvec,rotMat,tmar);
% BODY oriented trajectory translation vector
walkFetRot = cf(@(t,u)  sq(dot(t,u,3)),tvec,unvec);

wfrCat = cat(1,walkFetRot{:});
wfrMean = nanmean(wfrCat(nniz(wfrCat),:,:));
wfrStd = nanstd(wfrCat(nniz(wfrCat),:,:));
walkFetRot = cf(@(w,m,s) nunity(w,[],m,s),...
                walkFetRot,...
                repmat({wfrMean},1,numSessions),...
                repmat({wfrStd},1,numSessions));
clear('wfrCat','wfrMean','wfrStd')



% DECOMPOSE features ------------------------------------------------------------------------------
% @wfs
% EMBED wfet 
wfet = cf(@(x) x.copy('empty'), xyz);
cf(@(w,wfr) set(w,'data',[reshape(wfr,[],size(wfr,2)*size(wfr,3))]),wfet,walkFetRot);
wfs  = cf(@(w,e) w.segs([],e),wfet,embeddingWindow);
wfs =  cf(@(w,e) circshift(w,e/2,2),wfs,embeddingWindow);
wfs =  cf(@(w,x) MTADxyz('data',reshape(permute(w,[2,1,3]),size(w,2),[]),...
              'sampleRate',x.sampleRate),wfs,xyz);
for i = 1:numel(wfs), wfs{i}.data(isnan(wfs{i}.data(:)))=0; end

% DECOMPOSE wfet with svd for walk and turn periods for each session
[~,Sw,Vw] = cf(@(w,s) svd(w([s{'w+n'}]+[0,-0.25],:),0), wfs, Stc);

% DECOMPOSE wfet with svd for walk and turn periods within all sessions
wfw = cf(@(w,s) w([s{'w+n'}]+[0,-0.25],:), wfs, Stc);
[~,Sww,Vww] = svd(cat(1,wfw{:}),0);

% COMPUTE eigenvector loadings for each session's eigen vectors
afetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:10)),...
                          'sampleRate',w.sampleRate),...
           wfs,repmat({Vww},1,numSessions));
cf(@(f,x) set(f,'sync',x.sync.copy), afetW, xyz); cf(@(f,x) set(f,'origin',x.origin), afetW, xyz);

% COMPUTE eigenvector loadings for each session's eigen vectors
fetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:10)),'sampleRate',w.sampleRate),wfs,Vw);
cf(@(f,x) set(f,'sync',x.sync.copy), fetW, xyz); cf(@(f,x) set(f,'origin',x.origin), fetW, xyz);

% CREATE mask for eigenvectors for each session's eigen vectors
maskEigVec = ones(size(Vw{1}));
for i = [1:16,48:64], maskEigVec(i:embeddingWindow{1}:end,:) = 0;end
maskEigVec = repmat({maskEigVec},1,numSessions);
rVw = cf(@(v,m) v.*m,Vw,maskEigVec);

% REDUCED eigenvector loadings for  each session's eigen vectors
rfetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:10)),'sampleRate',w.sampleRate),wfs,rVw);
cf(@(r,x) set(r,'sync',x.sync.copy), rfetW, xyz); cf(@(r,x) set(r,'origin',x.origin), rfetW, xyz);


% COMPUTE features with reference session ----------------------------------------------------------
referenceSessionIndex = 1;
% REFERENCED compute eigenvector loadings
fetWsvd = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:10),[1,2],[1,2]),...
                            'sampleRate',w.sampleRate),...
             wfs,...
             repmat(rVw(referenceSessionIndex),1,numSessions));
cf(@(f,x) set(f,'sync',x.sync.copy),fetWsvd,xyz); cf(@(f,x) set(f,'origin',x.origin) ,fetWsvd,xyz);

fetWsvdPhase = cf(@(f) f.phase([1,5]),fetWsvd);
cf(@(f,x) set(f,'sync',x.sync.copy),fetWsvdPhase,xyz); cf(@(f,x) set(f,'origin',x.origin) ,fetWsvdPhase,xyz);

% @rfetWsvd
% REFERENCED reduced eigenvector loadings
rfetWsvd = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:10),[1,2],[1,2]),...
                            'sampleRate',w.sampleRate),...
              wfs,...
              repmat(rVw(referenceSessionIndex),1,numSessions));
cf(@(f,x) set(f,'sync',x.sync.copy),rfetWsvd,xyz);
cf(@(f,x) set(f,'origin',x.origin) ,rfetWsvd,xyz);
rfetWsvdPhase = cf(@(f) f.phase([1,5]),rfetWsvd);
cf(@(f,x) set(f,'sync',x.sync.copy),rfetWsvdPhase,xyz);
cf(@(f,x) set(f,'origin',x.origin) ,rfetWsvdPhase,xyz);

% @ts
% TIME vectors
wts = cf(@(e,s)  [1:e]./s,                        embeddingWindow,sampleRate);
ts =  cf(@(x)    [1:size(x,1)]./x.sampleRate,     xyz);

% @ang
% COMPUTE  intermarker angles 
ang = cf(@(t,x) create(MTADang,t,x), Trials, xyz);
for s = 1:numSessions, ang{s}.data(~nniz(xyz{s}),:,:,:) = 0;end





% END data processing ------------------------------------------------------------------



%% FIG.3
% DISPLAY eigen vectors for all sessions
% D svd
% $$$ hfig = figure;
% $$$ hfig.Units = 'centimeters';
% $$$ hfig.Position(3:4) = [30,4];
% $$$ hfig.PaperPositionMode = 'auto';
% $$$ for i = 1:20,
% $$$     subplot(2,10,i);imagesc(wts{1},1:size(wfet{1},2),reshape(Vww(:,i),[],size(wfet{1},2))'),
% $$$     caxis([-0.08,0.08]);
% $$$     axis xy
% $$$ end
% $$$ 
% $$$ TrialName = [sessionList(s).sessionName,'.',sessionList(s).mazeName,'.',sessionList(s).trialName];
% $$$ FigName = ['SVD_walkturn_',TrialName];
% $$$ print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
% $$$ print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
