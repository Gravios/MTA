function [Stc] = label_bhv_shake(Stc,Trial,varargin)
% function [Stc] = label_bhv_shake(Trial,Stc,varargin);
% feature set: fet_bref

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct(...
    'eigenVectorIndices',     2:4,...    
    'param', struct('sessionList',            'hand_labeled',...
                    'referenceTrial',         'jg05-20120317.cof.all',...
                    'featureSet',             'fet_bref',...
                    'svdState',               'shake',...
                    'sampleMode',             'centered',...
                    'sampleRate',             119.881035,...
                    'embeddingWindow',        64 ),...
    'display', false ...
);

[eigenVectorIndices,param,display] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------

OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'Shared/Behavior Paper/Figures/Figure_3/parts';

states = {'walk','rear','turn','pause','groom','sit','shake'};
sclr = 'brgcmyk';



if isempty(Trial),
    sessionList = get_session_list(param.sessionList);                    
elseif ischar(Trial),
    sessionList = get_session_list(Trial);                    
else
    sessionList = {Trial}; 
end

numSessions     = numel(sessionList);                    
sampleRate = repmat({param.sampleRate},1,numSessions);
embeddingWindow = repmat({param.embeddingWindow},1,numSessions);
trimWindow = repmat({[ param.embeddingWindow./4./param.sampleRate,...
                      -param.embeddingWindow./4./param.sampleRate]},1,numSessions);

% normalizationParameters = struct('sessionList',   param.sessionList   ,...
%                                  'referenceTrial',param.referenceTrial,...
%                                  'featureSet',    param.featureSet);

svdParameters = struct('sessionList',            param.sessionList,            ...
                       'referenceTrial',         param.referenceTrial,         ...
                       'featureSet',             param.featureSet,             ...
                       'svdState',               param.svdState,               ...
                       'sampleMode',             param.sampleMode,             ...
                       'sampleRate',             param.sampleRate,             ...
                       'embeddingWindow',        param.embeddingWindow         ...
);


% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create hash tag
%normalizationTag = DataHash(normalizationParameters);
svdTag           = DataHash(svdParameters);
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% LOAD Trial objects
if iscell(sessionList)
    Trials = cf(@(t)  MTATrial.validate(t),  sessionList);
else
    Trials = af(@(t)  MTATrial.validate(t),  sessionList);
end

% LOAD State Collections
if isempty(Stc),
    Stc = cf(@(t)  t.load('stc'),  Trials);
elseif ~iscell(Stc) && isa(Stc,'MTAStateCollection'),
    Stc = {Stc};
elseif ischar(Stc) && numSessions>1,
    Stc = cf(@(t,s)  t.load('stc',s), Trials,repmat({Stc},[1,numSessions]));
elseif ischar(Stc),    
    Stc = {Trial.load('stc',Stc)};
end

    
    
% LOAD Position data
xyz    = cf(@(t)   preproc_xyz(t)         , Trials);
         cf(@(x,s) x.resample(s)          , xyz,sampleRate);
fxyz   = cf(@(x)   x.copy()               , xyz);
         cf(@(f)   f.filter('ButFilter',5,2.4,'low'),  fxyz);

% LOAD and MAP Features to reference session
features = cf(@(t,p)   feval(p.featureSet,t),  Trials,repmat({param},[1,numSessions]));
           cf(@(f,t,p) f.map_to_reference_session(t,p.referenceTrial),...
              features,Trials,repmat({param},[1,numSessions]));
for s = 1:numSessions, features{s}.data(~nniz(xyz{s}),:,:) = 0;end
           cf(@(f,s)   f.resample(s),          features,sampleRate);

% NORMALIZE feature matrix
[refMean,refStd] = load_normalization_parameters_unity(param.featureSet,...
                                                       param.referenceTrial,...
                                                       param.sessionList);
cf(@(w,m,s) set(w,'data',nunity(w,[],m,s)),...
            features,...
            repmat({refMean},1,numSessions),...
            repmat({refStd},1,numSessions));
clear('featureCat','featureMean','featureStd')


% EMBBED feature
wfs  = cf(@(w,e) w.segs(1:size(w,1),e),features,embeddingWindow);
wfs =  cf(@(w,e) circshift(w,e/2,2),wfs,embeddingWindow);
wfs =  cf(@(w,x) MTADxyz('data',reshape(permute(w,[2,1,3]),size(w,2),[]),...
              'sampleRate',x.sampleRate),wfs,xyz);
for i = 1:numel(wfs), wfs{i}.data(isnan(wfs{i}.data(:)))=0; end


% EMBEDDED singular value decomposition
svdParameterFile = ...
    fullfile(MTASession([]).path.project,'analysis',...
             ['svd_parameters-',svdTag,'.mat']);
if ~exist(svdParameterFile),
% DECOMPOSE features with svd for a behavioral state within all Trials
    switch param.sampleMode
      case 'trimmed'
        wfw = cf(@(w,s,tw,sr) w([s{param.svdState,sr}]+tw,:), ...
                 wfs, Stc, trimWindow,sampleRate);    
      case 'centered'
        wfw = cf(@(w,s,sr) w(round(mean([s{param.svdState,sr}.data],2)),:), ...
                 wfs, Stc,sampleRate);
    end
    [~,Sww,Vww] = svd(cat(1,wfw{:}),0);
    save(svdParameterFile,'svdParameters','Sww','Vww');
    
    return
else
% LOAD SVD
    load(svdParameterFile);
end




if display
% FIGEIGVECT ------------------------------------------------------------------------------
    % project: MjgEd2016
    % parent: 
    % subplots: 1->N eigen vectors for all sessions during shake periods
    % location: label_bhv_shake.m
    wts = [1:param.embeddingWindow]./param.sampleRate;
    hfig = figure;
    hfig.Units = 'centimeters';
    hfig.Position(3:4) = [30,16];
    hfig.PaperPositionMode = 'auto';
    for i = 1:20,
        pc = reshape(Vww(:,i),[],size(features{1},2));
        pc = pc(:,[1:2:9,2:2:10,11:15,16:2:24,17:2:25,26:30]);
        subplot(2,10,i);imagesc(wts,1:size(features{1},2),pc'),
        caxis([-0.08,0.08]);
        %caxis([-.5,.5]);    
        axis xy
    end

    FigName = ['SVD_eigVect_' param.svdState];
    print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
% END FIGEIGVECT ---------------------------------------------------------------------------
end


% COMPUTE eigenvector loadings for each session's eigen vectors
% contains mask to select feature subset important to turning
sfet = cf(@(x) x.copy('empty'), xyz);
for i = eigenVectorIndices,
    cf(@(r,w,v) set(r,'data',[get(r,'data'),multiprod(w.data,v)]),...
       sfet,wfs,repmat({Vww(:,i)},1,numSessions));
end
cf(@(f,x) set(f,'sync',x.sync.copy), sfet, xyz); 
cf(@(f,x) set(f,'origin',x.origin), sfet, xyz);
cf(@(f,x) f.filter('ButFilter',5,[10,20],'bandpass'), sfet);

cf(@(f) set(f,'data',sqrt(sum((sum(circshift(f.segs(1:size(f,1),round(f.sampleRate/4),0),...
                                 round(f.sampleRate/8),2).^2,3))))'),...
   sfet);


if display,
    % TIME vectors
    ts =  cf(@(x)    [1:size(x,1)]./x.sampleRate,     xyz);
    
    figure,plot(sfet{1}.data),Lines(Stc{1}{'k'}(:),[],'r');
    skon = cf(@(f,s) f(s{'k',f.sampleRate}(:,1)),sfet,Stc);
    skoff = cf(@(f,s) f(s{'k',f.sampleRate}(:,2)),sfet,Stc);
    skdia = cf(@(f,s) f(round(mean(s{'k',f.sampleRate}.data,2))),sfet,Stc);
    figure,
    subplot(131)
    hist(vertcat(skon{:}),100)
    subplot(132)
    hist(vertcat(skdia{:}),100)
    subplot(133)
    hist(vertcat(skoff{:}),100)
end

% GENERATE shake periods



% REGRESSION TIME
% $$$ 
% $$$ shift = [-round(param.sampleRate),round(2*param.sampleRate)];
% $$$ % OPTIMIZATION of state transitions
% $$$ sts = cf(@(s,f,t) bsxfun(@plus,[s{'k',f.sampleRate}(:,1)],ones([1,2]).*round(0.2*f.sampleRate).*[-1,1]),...
% $$$          Stc, sfet);
% $$$ 
% $$$ 
% $$$ % TRIM timepoints at ends if necessary
% $$$ for s = 1:numSessions, 
% $$$     sts{s}(sum([sts{s}+repmat(shift,size(sts{s},1),1)<=0, ...
% $$$         sts{s}+repmat(shift,size(sts{s},1),1)>size(sfet{s},1)],1)>0,:)=[];
% $$$ end        
% $$$ 
% $$$ % COLLECT segments
% $$$ wfw = cf(@(w,s,t,i) w.segs(round(mean(s,2))-round(1.*w.sampleRate),round(3.*w.sampleRate)),...
% $$$          sfet, sts, Trials);        
% $$$ % COCATENATE trial data
% $$$ onf = cat(2,wfw{:});
% $$$ fetSegs = onf(:,:,:);
% $$$ 
% $$$ windowIndices = 101:181;
% $$$ 
% $$$ % COMPUTE mean trajectory
% $$$ mfs = repmat({sq(nanmean(fetSegs(windowIndices,:,:),2))'},1,numSessions); %rear
% $$$ 
% $$$ rof = cf(@(f) f.copy(),sfet);
% $$$ cf(@(f,m) set(f,'data',circshift(f.segs(1:size(f,1),size(m,2)),round(size(m,2)/2),2)),...
% $$$          rof,mfs);
% $$$ cf(@(f,m) set(f,'data',repmat(permute(f.data,[4,1,2,3]),3,1,1)-...
% $$$                      repmat(linspace(-.5,.5,3)',[1,size(m,2),size(f,2),size(f,3)])),...
% $$$   rof,mfs);
% $$$ 
% $$$ csw = cf(@(f) f.copy(),rof);
% $$$ cf(@(c,f,m) set(c,'data',sq(sum((f.data-repmat(permute(m,[3,2,4,1]),size(f,1),1,size(f,3))).^2,2))),...
% $$$    csw,rof,mfs);
% $$$ 
% $$$ csn = cf(@(c) c.copy(), csw);
% $$$ cf(@(c) set(c,'data',nunity(c(2,:,1)')),csn);

sfetThreshold = 20;
nkper = cf(@(f,t) ThreshCross(abs(f.data),t,round(0.05.*f.sampleRate)), ...
           sfet,repmat({sfetThreshold},1,numSessions));
nkper = cf(@(p,f) bsxfun(@plus,p,round([0.042,-0.042].*f.sampleRate)),...
           nkper,sfet);
for s = 1:numSessions,
    kind = diff(nkper{s},1,2)<round(0.17.*sfet{s}.sampleRate);
    nkper{s}(kind,:) = bsxfun(@plus,nkper{s}(kind,:),round([-0.01,0.01].*sfet{s}.sampleRate));
end



if display,
s = 1;
sp = [];
figure,
sp(end+1)=subplot2(6,1,[1:4],1);
hold('on');
plot(ts{s},sfet{s}(:))
%plot(ts{s},nansum(csn{s}(:),3)),
Lines(mean(sts{s},2)./sampleRate{s},[],'g');
Lines(nkper{s}(:)./sampleRate{s},[],'m');
sp(end+1)=subplot2(6,1,5:6,1);
plotSTC(Stc{s},1,'text',states,sclr,[],false);
linkaxes(sp,'x');
end


% $$$ [nsmins,nsvals] = cf(@(c) LocalMinima(sum(c(2,:,1),3),60,regressionThreshold),csw);
% $$$ %nsmins = cf(@(m,s) SelectPeriods(m,bsxfun(@plus,s,[-10,10]),'d',1),nsmins,sts);
% $$$ [nsmins,nsinds] = cf(@(m,s) SelectPeriods(m,bsxfun(@plus,[s{'shake'}(:,1)],[-60,60]),'d',1),...
% $$$                      nsmins,Stc);
% $$$ [ssmins,ssinds] = cf(@(m,s) SelectPeriods([s{'shake'}(:,1)],bsxfun(@plus,m,[-60,60]),'d',1),...
% $$$                      nsmins,Stc);
% $$$ csegs = cf(@(a,m) a.segs(m-120,360),sfet,nsmins);
% $$$ 
% $$$ [nkonind] = cf(@(m,s) SelectPeriods(m(:,2),bsxfun(@plus,[s{'shake'}(:,2)],[-60,60]),'d',1),...
% $$$                      nkper,Stc);
% $$$ [skonind] = cf(@(m,s) SelectPeriods([s{'shake'}(:,2)],bsxfun(@plus,m(:,2),[-60,60]),'d',1),...
% $$$                      nkper,Stc);
% $$$ csegs = cf(@(a,m) a.segs(m-120,360),sfet,nkonind);
% $$$ 
% $$$ 
% $$$ 
% $$$ [sccg,txx,pxx] = cf(@(s,n,sr) CCG([s;n],[ones(size(s));2*ones(size(n))],...
% $$$                                   2,40,sr,[1,2],'count'),...
% $$$                     skonind,nkonind,sampleRate);
% $$$ accg = sum(cat(4,sccg{:}),4);

%medianCorrectionOffset = median(cat(1,nsmins{:})-cat(1,ssmins{:}))./sampleRate{1};
% $$$ medianCorrectionOffset = median(cat(1,nkonind{:})-cat(1,skonind{:}))./sampleRate{1};
% $$$ 
% $$$ figure();
% $$$ hist(cat(1,nsvals{:}),100)
% $$$ 
% $$$ kdur = cf(@(s) diff(s{'k'}.data,1,2),Stc);
% $$$ kdur = cat(1,kdur{:});
% $$$ figure,hist(kdur,100)
% $$$ 
% $$$ std(kdur)
% $$$ kdur*2

if display,
    % SUPFIG Rear
    figure();
    bar(txx{1},accg(:,1,2));
    Lines(medianCorrectionOffset.*1000,[],'r');
    title(['CCG between shake onsets and local minima of onset regression'])
    xlabel('Time shift(ms) centered on local minima')
    ylabel('count')

    FigName = 'State_transition_regression_realignment_PC1_all_to_rearOn';
    print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
end


aper = cf(@(t,f) resample(t.sync.copy(),f.sampleRate),Trials,sfet);
cf(@(p) set(p,'data',p.data-p.data(1)+1),aper);
cf(@(p) p+[2,-2],aper);

cf(@(s,p) cf(@(sts,per) set(sts,'data',get(sts-bsxfun(@plus,per,[-1,1]),'data')),...
             s.states,repmat({p},[1,numel(s.states)])),...
   Stc,nkper);

cf(@(t,s,f,p) s.addState(t.spath,t.filebase,p,f.sampleRate,t.sync.copy,t.sync.data(1),'shake','k'),...
   Trials,Stc,sfet,nkper);
cf(@(s,a) set(s.states{s.gsi('k')},'data',get(s.states{s.gsi('k')}&a.data,'data')), ...
   Stc,aper);


for s = 1:numel(Stc),
    [stcMatrix] = stc2mat(Stc{s},xyz{s},states);
    stateVectors = bsxfun(@times,eye(numel(states)),1:numel(states));
    for sts = states(1:end-1),
        sts = sts{1};
        [stcMatrix] = reassign_low_duration_state_to_neighboring_states(stcMatrix,...
                                                  [Stc{s}{sts}],0.18*sfet{s}.sampleRate);
    end
    Stc{s} = mat2stc(stcMatrix,Stc{s},xyz{s},Trials{s});
end




cf(@(s) s.save(1), Stc);

if numel(Stc)==1,
    Stc = Stc{1};
end

    

% END MAIN -----------------------------------------------------------------------------------------

