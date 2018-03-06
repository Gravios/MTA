function [projectionScore, model] = decompose_xy_motion_wrt_body(Trial,varargin)


% TESTING VARS
% $$$ sessionListName = 'hand_labeled';
% $$$ newSampleRate   = 119.881035;
% $$$ windowLength = 0.53386259;
% $$$ orientations = [0,90];
% $$$ maxNumComponents  = 5;
% $$$ mode = 'COMPUTE';
% $$$ tag = '';


global MTA_PROJECT_PATH;


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('sessionListName', 'hand_labeled',      ...   String
                 'newSampleRate',    119.881035,         ...   Hertz
                 'windowLength',     0.53386259,         ...   Seconds
                 'xyzPreprocOpts',   {{'LOAD_TRB_XYZ'}}, ...   Cellstr
                 'orientations',     [0,90],             ...   Degrees
                 'maxNumComponents', 5,                  ...   Numeric Int
                 'mode',             'RUN',              ...   String
                 'tag',              ''                  ...   String
);

[sessionListName, newSampleRate, windowLength, ...
 xyzPreprocOpts,  orientations,  maxNumComponents,...
 mode,            tag] = DefaultArgs(varargin,defargs,'--struct');

%--------------------------------------------------------------------------------------------------



% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create hash tag
metadata = struct('sessionListName',  sessionListName,     ...   String
                  'newSampleRate',    newSampleRate,      ...   Hertz
                  'windowLength',     windowLength,       ...   Seconds
                  'xyzPreprocOpts',   xyzPreprocOpts,     ...   Cellstr
                  'orientations',     orientations        ...   Degrees
);
if isempty(tag),
    tag = DataHash(metadata);
end
modelPath = fullfile(MTA_PROJECT_PATH,'analysis','models',['decompose_xy_motion_wrt_body-',tag,'.mat']);
%---------------------------------------------------------------------------------------------------


% ASSERTIONS ---------------------------------------------------------------------------------------
assert(~isempty(regexp(mode,'(^COMPUTE$)|(^RUN$)')),...
       'MTA:transformations:decompose_xy_motion_wrt:UnknownMode');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------
% LOAD Session(s) ----------------------------------------------------------------------------------

switch mode
  case 'COMPUTE',
    Trials = af(@(Trial)  MTATrial.validate(Trial),  get_session_list(sessionListName));
  case 'RUN',
    try,
        Trials = {MTATrial.validate(Trial)};
    catch 
        Trials = af(@(Trial)  MTATrial.validate(Trial),  get_session_list(Trial));
    end
    model = load(modelPath);
end    
numSessions= numel(Trials);

% ENCAPSULATE non cell vars into cell arrayss
sampleRate      = repmat({newSampleRate},1,numSessions);
embeddingWindow = repmat({round(windowLength.*newSampleRate)},1,numSessions);
rotationAngles  = repmat({deg2rad(orientations)},1,numSessions);%[0,45,90,135]

% LOAD basic data session data
xyz    = cf(@(Trial)  preproc_xyz(Trial,'trb'),  Trials);
cf(@(x,s) x.resample(s),xyz,sampleRate);

% ADD virtual markers for body
rb    = cf(@(x)   x.model.rb({'spine_lower','pelvis_root','spine_middle'}), xyz);
bcom  = cf(@(x,r) x.com(r), xyz,rb);
cf(@(x,b) x.addMarker('fbcom','data',ButFilter(b,4,[0.5]./(x.sampleRate/2),'low')),xyz,bcom);
cf(@(x,b) x.addMarker('bcom','data',b), xyz,bcom);
cf(@(x,b) x.addMarker('fsl','data',ButFilter(x(:,'spine_lower',:),4,[0.5]./(x.sampleRate/2),'low')),xyz,bcom);
clear('bcom');



% COMPUTE XY motion projections onto orthonormal body basis -------------------------------------------

% NORM body vector projection
tmar = repmat({{'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'}},1,numSessions);
tvec = cf(@(x,m,t) repmat(circshift(x(:,m,[1,2]),-(round(x.sampleRate*0.05)/2))...
                          -circshift(x(:,m,[1,2]), (round(x.sampleRate*0.05)/2)),...
                          [1,1,1,numel(rotationAngles{1})]),...
          xyz,tmar,rotationAngles);

% CREATE orthonormal body basis
bodyUnitVector = cf(@(x,r) repmat(x(:,'spine_upper',[1,2])-x(:,'spine_lower',[1,2]),1,1,1,numel(r)),xyz,rotationAngles);
rotMat = cf(@(t,m) reshape(repmat(permute([cos(t),-sin(t);...
                    sin(t),cos(t)],...
                                          [3,1,2]),...
                                  [size(m,1),1,1]),...
                           size(m,1),2,2,numel(t)),rotationAngles,bodyUnitVector);
unvec  = cf(@(m,r,t) repmat(bsxfun(@rdivide,multiprod(m,r,[2,3],[2,3]),sqrt(sum(m.^2,3))),...
                            [1,numel(t),1,1]),bodyUnitVector,rotMat,tmar);


% PROJECT marker trajectories onto orthonormal body basis
walkFetRot = cf(@(t,u)  sq(dot(t,u,3))  ,tvec,unvec);
switch mode
  case 'COMPUTE'
% COMPUTE mean and standard deviation for normalization
    wfrCat = cat(1,walkFetRot{:});
    model.mean = nanmean(wfrCat(nniz(wfrCat),:,:));
    model.std  = nanstd(wfrCat(nniz(wfrCat),:,:));
    clear('wfrCat')
end


% NORMALIZE projections
walkFetRot = cf(@(w,m,s) nunity(w,[],m,s),...
                walkFetRot,...
                repmat({model.mean},1,numSessions),...
                repmat({model.std },1,numSessions));



% DECOMPOSE features ------------------------------------------------------------------------------
% EMBED wfet 
wfet = cf(@(x) x.copy('empty'), xyz);
cf(@(w,wfr) set(w,'data',[reshape(wfr,[],size(wfr,2)*size(wfr,3))]),wfet,walkFetRot);
wfs  = cf(@(w,e) w.segs(1:size(w,1),e),wfet,embeddingWindow);
wfs =  cf(@(w,e) circshift(w,e/2,2),wfs,embeddingWindow);
wfs =  cf(@(w,x) MTADxyz('data',reshape(permute(w,[2,1,3]),size(w,2),[]),...
                         'sampleRate',x.sampleRate),wfs,xyz);
for i = 1:numel(wfs), wfs{i}.data(isnan(wfs{i}.data(:)))=0; end

switch mode
  case 'COMPUTE';
% LOAD default StateCollections
    Stc  = cf(@(Trial) Trial.load('stc'), Trials);
    model.eigenValues  = [];
    model.eigenVectors = [];
% DECOMPOSE wfet with svd for walk and turn periods within all sessions
    wfw = cf(@(w,s) w([s{'w+n'}]+[0,-0.25],:), wfs, Stc);
    [~,model.eigenValues,model.eigenVectors] = svd(cat(1,wfw{:}),0);
% SAVE model
    model.metadata = metadata;
    save(modelPath,'-struct','model');
    if nargout==0, return, end
end

% COMPUTE eigenvector loadings for each session's eigen vectors
projectionScore = cf(@(w,v,nc) MTADxyz('data',multiprod(w.data,v(:,1:nc)),...
                          'sampleRate',w.sampleRate),...
           wfs,...
           repmat({model.eigenVectors},1,numSessions),...
           repmat({maxNumComponents},1,numSessions));
cf(@(f,x) set(f,'sync',x.sync.copy), projectionScore, xyz); 
cf(@(f,x) set(f,'origin',x.origin),  projectionScore, xyz);
cf(@(f,m) f.update_hash(m),  projectionScore, repmat({metadata},size(projectionScore)));

if numel(projectionScore)==1,
    projectionScore = projectionScore{1};
end





% VISUALISATION 
% @ts
% TIME vectors
% $$$ wts = cf(@(e,s)  [1:e]./s,                        embeddingWindow,sampleRate);
% $$$ ts =  cf(@(x)    [1:size(x,1)]./x.sampleRate,     xyz);
% $$$ 
% $$$ hfig = figure;
% $$$ hfig.Units = 'centimeters';
% $$$ hfig.Position(3:4) = [30,4];
% $$$ hfig.PaperPositionMode = 'auto';
% $$$ for i = 1:20,
% $$$     subplot(2,10,i);imagesc(wts{1},1:size(wfet{1},2),reshape(Vww(:,i),[],size(wfet{1},2))'),
% $$$     caxis([-0.08,0.08]);
% $$$     axis xy
% $$$ end
% $$$ s = 1;
% $$$ stepFeature = ButFilter(afetW{s}(:,4),3,[1,6]./(afetW{s}.sampleRate/2),'bandpass');
% $$$ [steps,swayMagnitude] = LocalMinima(-abs(stepFeature),0,-4);
% $$$ states = {'walk','turn','pause','rear'};
% $$$ sclr = 'bgcr';
% $$$ sp = [];
% $$$ figure,
% $$$ sp(end+1) = subplot2(10,4,1:8,1:4); hold on
% $$$ plot(ts{s},stepFeature);
% $$$ plot(ts{s}(steps),stepFeature(steps),'or')
% $$$ plot(ts{s},-afetW{s}(:,1)./2+20);
% $$$ plot(ts{s},-reshape(walkFetRot{s},[],size(wfet{s},2))*mean(reshape(Vww(:,1),[],size(wfet{s},2)))'.*20+20);
% $$$ plot(ts{s},-abs(afetW{s}(:,2))./2-20);
% $$$ plot(ts{s},-abs(reshape(walkFetRot{s},[],size(wfet{s},2))*mean(reshape(Vww(:,2),[],size(wfet{s},2)))').*20-20);
% $$$ Lines([],+20,'k');
% $$$ legend({'Step Feature','Putative Steps','PC1 Forward displacement','PC1 raw','PC2 turning','PC2 raw'})
% $$$ sp(end+1)=subplot2(10,4,[9,10],[1:4]);
% $$$ %plotSTC(StcNN{s},1,'text',states,sclr,[],false);
% $$$ plotSTC(Stc{s},1,'text',states,sclr,[],false);
% $$$ linkaxes(sp,'x');

% END MAIN -----------------------------------------------------------------------------------------