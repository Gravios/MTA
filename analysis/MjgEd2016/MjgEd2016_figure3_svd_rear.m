
%% Fig.3 of Rearing Characterization %%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis of horizontal motion relative of to the body
%
%

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
%rotationAngles  = [0,45,90,135];
sampleRate   = repmat({119.881035}            ,1,numSessions);
featureName  = repmat({'MjgEd2016-F3-SVDREAR'},1,numSessions);
featureLabel = repmat({'fet_bref'}            ,1,numSessions);
featureKey   = repmat({'r'}                   ,1,numSessions);
embeddingWindow = repmat({64},1,numSessions);
%stcMode = 'NN0317';

% LOAD Sessions and basic data
Trials = af(@(Trial) MTATrial.validate(Trial), sessionList);
Stc    = cf(@(Trial) Trial.load('stc')       , Trials);
xyz    = cf(@(Trial) Trial.load('xyz')       , Trials);
% RESAMPLE to common sampleRate
cf(@(x,s) x.resample(s),xyz,sampleRate);
fet    = cf(@(Trial) fet_bref(Trial), Trials);
cf(@(f,t) f.map_to_reference_session(t,'jg05-20120317.cof.all'),fet,Trials);

% $$$ % ADD virtual markers for body
% $$$ rb    = cf(@(x)   x.model.rb({'spine_lower','pelvis_root','spine_middle'}), xyz);
% $$$ bcom  = cf(@(x,r) x.com(r), xyz,rb);
% $$$ cf(@(x,b) x.addMarker('fbcom','data',ButFilter(b,4,[0.5]./(x.sampleRate/2),'low')),xyz,bcom);
% $$$ cf(@(x,b) x.addMarker('bcom','data',b), xyz,bcom);
% $$$ cf(@(x,b) x.addMarker('fsl','data',ButFilter(x(:,'spine_lower',:),4,[0.5]./(x.sampleRate/2),'low')),xyz,bcom);
% $$$ clear('bcom');
% $$$ % ADD virtual markers for head
% $$$ rb    = cf(@(x)   x.model.rb({'head_back','head_left','head_front','head_right'}),xyz);
% $$$ hcom  = cf(@(x,r) x.com(r),xyz,rb);
% $$$ cf(@(x,h) x.addMarker('fhcom','data',ButFilter(h,4,[0.8]./(x.sampleRate/2),'low')),xyz,hcom);
% $$$ cf(@(x,h) x.addMarker('hcom','data',h), xyz,hcom);
% $$$ clear('hcom');
% $$$  
% $$$ tmar = repmat({{'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'}},1,numSessions);
% $$$ 
% $$$ xyzSubset  = cf(@(x) x.copy,xyz);
% $$$ cf(@(x,m) set(x,'data',x(:,m,:)) ,xyzSubset,tmar);
% $$$ xyzSubDiff = cf(@(x) markerDiffMatrix(x),xyzSubset);
% $$$ 
% $$$ xyd = cf(@(x) reshape(sqrt(sum(x(:,:,:,[1,2]).^2,4)),size(x,1),[]),xyzSubDiff);
% $$$ zd  =  cf(@(x) reshape(x(:,:,:,[1,2]),size(x,1),[]),xyzSubDiff);
% $$$ 
% $$$ uindMat = repmat({logical(triu(ones([numel(tmar{1}),numel(tmar{1})]),1))},1,numSessions);
% $$$ xyvec = cf(@(x,u) reshape(x(reshape(repmat(permute(u,[3,1,2]),[size(x,1),1,1]),[],1)),size(x,1),[]),xyd,uindMat);
% $$$ zvec  = cf(@(x,u) reshape(x(reshape(repmat(permute(u,[3,1,2]),[size(x,1),1,1]),[],1)),size(x,1),[]),zd,uindMat);
% $$$ 
% $$$ zvec = cf(@(x,z) cat(2,x,z) ,xyvec,zvec);
% $$$ 
% $$$ 
% $$$ MTADfet.encapsulate(Trial,data,sampleRate,featureName,featureLabel,featureKey)

zfrCat = cf(@(f) get(f,'data'),fet);
zfrCat = cat(1,zfrCat{:});
zfrMean = nanmean(zfrCat(nniz(zfrCat),:,:));
zfrStd = nanstd(zfrCat(nniz(zfrCat),:,:));
fet    = cf(@(w,m,s) nunity(w,[],m,s),...
                fet,...
                repmat({zfrMean},1,numSessions),...
                repmat({zfrStd},1,numSessions));
clear('zfrCat','zfrMean','zfrStd')



% DECOMPOSE features ------------------------------------------------------------------------------
% @wfs
% EMBED wfet 
wfet = cf(@(x) x.copy('empty'), xyz);
cf(@(w,z) set(w,'data',z),wfet,fet);
wfs  = cf(@(w,e) w.segs([],e),wfet,embeddingWindow);
wfs =  cf(@(w,e) circshift(w,e/2,2),wfs,embeddingWindow);
wfs =  cf(@(w,x) MTADxyz('data',reshape(permute(w,[2,1,3]),size(w,2),[]),...
              'sampleRate',x.sampleRate),wfs,xyz);
for i = 1:numel(wfs), wfs{i}.data(isnan(wfs{i}.data(:)))=0; end

% DECOMPOSE wfet with svd for walk and turn periods for each session
[~,Sw,Vw] = cf(@(w,s) svd(w([s{'r'}]+[0,-0.25],:),0), wfs, Stc);

% DECOMPOSE wfet with svd for walk and turn periods within all sessions
wfw = cf(@(w,s) w([s{'r'}]+[0,-0.25],:), wfs, Stc);
[~,Sww,Vww] = svd(cat(1,wfw{:}),0);
