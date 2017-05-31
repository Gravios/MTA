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
Trials = af(@(Trial) MTATrial.validate(Trial)  , sessionList);
Stc    = cf(@(Trial) Trial.load('stc')         , Trials);
StcNN  = cf(@(Trial) Trial.load('stc','NN0317'), Trials);
xyz    = cf(@(Trial) Trial.load('xyz')         , Trials);
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


% COMPUTE Features ------------------------------------------------------------------------------
% @walkFetRot
% NORM body vector projection
tmar = repmat({{'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'}},1,numSessions);
tvec = cf(@(x,m,t) repmat(circshift(x(:,m,[1,2]),-(round(x.sampleRate*0.05)/2))...
                          -circshift(x(:,m,[1,2]), (round(x.sampleRate*0.05)/2)),...
                          [1,1,1,numel(rotationAngles{1})]),...
          xyz,tmar,rotationAngles);
% BODY vector
bodyUnitVector = cf(@(x)  repmat(x(:,'fbcom',[1,2])-x(:,'fsl',[1,2]),1,1,1,4),xyz);
%bodyUnitVector = cf(@(x)  repmat(x(:,'spine_upper',[1,2])-x(:,'spine_lower',[1,2]),1,1,1,4),xyz);
rotMat = cf(@(t,m) reshape(repmat(permute([cos(t),-sin(t);...
                                           sin(t),cos(t)],...
                                          [3,1,2]),...
                                  [size(m,1),1,1]),...
                           size(m,1),2,2,numel(t)),rotationAngles,bodyUnitVector);
unvec  = cf(@(m,r,t) repmat(bsxfun(@rdivide,multiprod(m,r,[2,3],[2,3]),sqrt(sum(m.^2,3))),...
                            [1,numel(t),1,1]),bodyUnitVector,rotMat,tmar);
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
hfig = figure;
hfig.Units = 'centimeters';
hfig.Position(3:4) = [30,4];
hfig.PaperPositionMode = 'auto';
for i = 1:20,
    subplot(2,10,i);imagesc(wts{1},1:size(wfet{1},2),reshape(Vww(:,i),[],size(wfet{1},2))'),
    caxis([-0.08,0.08]);
    axis xy
end

% $$$ TrialName = [sessionList(s).sessionName,'.',sessionList(s).mazeName,'.',sessionList(s).trialName];
% $$$ FigName = ['SVD_walkturn_',TrialName];
% $$$ print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
% $$$ print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% $$$ 
% $$$ % DISPLAY eigen vectors for specific session
% $$$ % D svd
% $$$ hfig = figure;
% $$$ hfig.Units = 'centimeters';
% $$$ hfig.Position(3:4) = [30,4];
% $$$ hfig.PaperPositionMode = 'auto';
% $$$ for i = 1:10,
% $$$     subplot(1,10,i);imagesc(wts{s},1:size(wfet{s},2),reshape(Vw{s}(:,i),[],size(wfet{s},2))'),
% $$$     caxis([-0.08,0.08]);
% $$$     axis xy
% $$$ end
% $$$ 
% $$$ TrialName = [sessionList(s).sessionName,'.',sessionList(s).mazeName,'.',sessionList(s).trialName];
% $$$ FigName = ['SVD_walkturn_',TrialName];
% $$$ print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
% $$$ print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
% $$$ 
% $$$ % DISPLAY eigen vectors with loadings and features
% $$$ s = 1;
% $$$ states = {'walk','turn','pause','rear'};
% $$$ sclr = 'bgcr';
% $$$ pcLabels = {'PC1','PC2','PC3','PC4','PC5'};
% $$$ hfig = figure;
% $$$ hfig.Units = 'centimeters';
% $$$ hfig.Position(3:4) = [25,16];
% $$$ hfig.PaperPositionMode = 'auto';
% $$$ sp = [];
% $$$ sp(end+1)=subplot2(10,4,[1:8],[2:4]);
% $$$ hold on;
% $$$ plot(ts{s},walkFetRot{s}(:,3,1).*5+60,'LineWidth',1); % spine sway 'spine_lower'
% $$$ plot(ts{s},walkFetRot{s}(:,3,2).*5+70,'LineWidth',1); % spine sway 'pelvis_root'
% $$$ plot(ts{s},walkFetRot{s}(:,3,3).*5+80,'LineWidth',1); % spine sway 'spine_middle'
% $$$ plot(ts{s},walkFetRot{s}(:,3,4).*5+90,'LineWidth',1); % spine sway 'spine_upper'
% $$$ plot(ts{s},fetW{s}(:,1),'b','LineWidth',1);           % Walk comp
% $$$ plot(ts{s},fetW{s}(:,2),'r','LineWidth',1);           % Walk comp
% $$$ plot(ts{s},fetW{s}(:,3),'y','LineWidth',1);           % Turn comp
% $$$ plot(ts{s},fetW{s}(:,4),'g','LineWidth',1);           % Turn comp
% $$$ % $$$ plot(ts{s},fetW{s}(:,5),'m','LineWidth',1);           % Turn comp
% $$$ legend({'spine lower','pelvis root','spine middle','spine upper','PC1','PC2','PC3','PC4','PC5'})
% $$$ Lines([],0,'k');
% $$$ Lines([],5,'r');
% $$$ Lines([],-5,'r');
% $$$ sp(end+1)=subplot2(10,4,[9,10],[2:4]);
% $$$ plotSTC(Stc{s},1,'text',states,sclr,[],false);
% $$$ linkaxes(sp,'x');
% $$$ for i = 2:2:10,
% $$$     sp(end+1)=subplot2(10,4,i-1:i,1);
% $$$     imagesc(wts{s},1:size(wfet{s},2),reshape(Vw{s}(:,i/2),[],size(wfet{s},2))');
% $$$     axis xy;
% $$$     caxis([-0.08,0.08]);
% $$$     ylabel(['PC',num2str(i/2)])
% $$$ end
% $$$ 
% $$$ TrialName = [sessionList(s).sessionName,'.',sessionList(s).mazeName,'.',sessionList(s).trialName];
% $$$ FigName = ['SVD_walkturn_timeseries_',TrialName];
% $$$ print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
% $$$ print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ % DISPLAY eigen vectors with loadings for each PC and reduced PC
% $$$ s = 1;
% $$$ states = {'walk','turn','pause','rear'};
% $$$ sclr = 'bgcr';
% $$$ pcLabels = {'PC1','PC2','PC3','PC4','PC5'};
% $$$ hfig = figure;clf
% $$$ hfig.Units = 'centimeters';
% $$$ hfig.Position(3:4) = [25,16];pause(0.1);
% $$$ hfig.PaperPositionMode = 'auto';
% $$$ sp = [];
% $$$ 
% $$$ for i = 2:2:10,
% $$$     subplot2(10,4,i-1:i,1);
% $$$     imagesc(wts{s},1:size(wfet{s},2),reshape(rVw{s}(:,i/2),[],size(wfet{s},2))');
% $$$     axis xy;
% $$$     caxis([-0.08,0.08]);
% $$$     ylabel(['PC',num2str(i/2)])
% $$$     if i ~= 10,
% $$$         sp(end+1)=subplot2(10,4,i-1:i,[2:4]);    hold on
% $$$         plot(ts{s},fetW{s}(:,i/2),'g','LineWidth',1);           % Turn comp
% $$$         plot(ts{s},rfetW{s}(:,i/2),'b','LineWidth',1);           % Turn comp
% $$$         hold on;
% $$$         Lines([],0,'k');
% $$$         Lines([],5,'r');
% $$$         Lines([],-5,'r');
% $$$ 
% $$$         legend({pcLabels{i/2},[pcLabels{i/2},'R']})
% $$$     end
% $$$ end
% $$$ 
% $$$ sp(end+1)=subplot2(10,4,[9,10],[2:4]);
% $$$ plotSTC(Stc{s},1,'text',states,sclr,[],false);
% $$$ 
% $$$ linkaxes(sp,'x');
% $$$ 
% $$$ TrialName = [sessionList(s).sessionName,'.',sessionList(s).mazeName,'.',sessionList(s).trialName];
% $$$ FigName = ['SVD_walkturn_timeseries_Reduced_MID_',TrialName];
% $$$ print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
% $$$ print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ % DISPLAY eigen vectors with loadings for each PC and reduced PC
% $$$ s = 6;
% $$$ referenceSessionIndex = 1;
% $$$ states = {'walk','turn','pause','rear'};
% $$$ sclr = 'bgcr';
% $$$ pcLabels = {'PC1','PC2','PC3','PC4','PC5'};
% $$$ 
% $$$ hfig = figure;clf
% $$$ hfig.Units = 'centimeters';
% $$$ hfig.Position(3:4) = [25,16];pause(0.1);
% $$$ hfig.PaperPositionMode = 'auto';
% $$$ sp = [];
% $$$ 
% $$$ for i = 2:2:10,
% $$$     subplot2(10,5,i-1:i,1);
% $$$     imagesc(wts{s},1:size(wfet{s},2),reshape(Vw{referenceSessionIndex}(:,i/2),[],size(wfet{s},2))');
% $$$     axis xy;
% $$$     caxis([-0.08,0.08]);
% $$$     ylabel(['PC Ref',num2str(i/2)])
% $$$     subplot2(10,5,i-1:i,2);
% $$$     imagesc(wts{s},1:size(wfet{s},2),reshape(Vw{s}(:,i/2),[],size(wfet{s},2))');
% $$$     axis xy;
% $$$     caxis([-0.08,0.08]);
% $$$     ylabel(['PC',num2str(i/2)])
% $$$     if i ~= 10,
% $$$         sp(end+1)=subplot2(10,5,i-1:i,[3:5]);    hold on
% $$$         plot(ts{s},fetW{s}(:,i/2),'g','LineWidth',1);           % Turn comp
% $$$         plot(ts{s},fetWsvd{s}(:,i/2),'b','LineWidth',1);           % Turn comp
% $$$         hold on;
% $$$         Lines([],0,'k');
% $$$         Lines([],5,'r');
% $$$         Lines([],-5,'r');
% $$$ 
% $$$         legend({pcLabels{i/2},[pcLabels{i/2},'R']})
% $$$     end
% $$$ end
% $$$ sp(end+1)=subplot2(10,5,[9,10],[3:5]);
% $$$ plotSTC(Stc{s},1,'text',states,sclr,[],false);
% $$$ linkaxes(sp,'x');
% $$$ 
% $$$ RefTrialName = [sessionList(referenceSessionIndex).sessionName,'.',...
% $$$                 sessionList(referenceSessionIndex).mazeName,'.',...
% $$$                 sessionList(referenceSessionIndex).trialName];
% $$$ TrialName = [sessionList(s).sessionName,'.',sessionList(s).mazeName,'.',sessionList(s).trialName];
% $$$ FigName = ['SVD_walkturn_timeseries_REF_',RefTrialName,'_TARGET_',TrialName];
% $$$ print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
% $$$ print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ % DISPLAY Skeleton trajectory of spine sway
% $$$ %s = 6;i = 3;%2,11,13,22;tsubset = [60,183];[30,143];
% $$$ %s = 2;i = 1,17,20;
% $$$ %s = 3;i = 2,3,9,
% $$$ %s = 4;i = 4,5,13;
% $$$ %s = 5;i = 30,32,35
% $$$ 
% $$$ s = 2;
% $$$ i = 1;
% $$$ tsubset = [50,120];%,320];
% $$$ tsubset = [80,180];%,320];
% $$$ statePeriods = Stc{s}{'w'}.data;
% $$$ statePeriods(diff(statePeriods,1,2)<2*xyz{s}.sampleRate,:) = [];
% $$$ %index = 66600;statePeriods = [index-50,index+450];
% $$$ %for i = 1:size(statePeriods,1)
% $$$ hfig = figure(20170321);clf
% $$$ hfig.Units = 'centimeters';
% $$$ hfig.Position(3:4) = [25,16];pause(0.1);
% $$$ hfig.PaperPositionMode = 'auto';
% $$$ ts = [1:diff(statePeriods(i,:))+1]/fetWsvd{s}.sampleRate;
% $$$ subplot2(4,1,[1],1); hold on
% $$$ plot(ts{s},fetWsvd{s}(statePeriods(i,:),3),'LineWidth',1)
% $$$ plot(ts(tsubset(1):tsubset(2)),fetWsvd{s}(statePeriods(i,[1,1])+tsubset,3)+10,'LineWidth',2)
% $$$ % $$$ scatter(ts(tsubset(1):tsubset(2)),...
% $$$ % $$$         fetWsvd{s}(statePeriods(i,[1,1])+tsubset,3),...
% $$$ % $$$         10,...
% $$$ % $$$         mc(phaseBin(statePeriods(i,1)+[tsubset(1):tsubset(2)]),:),...
% $$$ % $$$         'filled')
% $$$ 
% $$$ colorbar
% $$$ subplot2(4,1,[2:4],1);
% $$$ hold on
% $$$ phaseBins = linspace(-pi,pi,100);
% $$$ [~,~,phaseBin] = histcounts(fetWsvdPhase{s}(:,3),phaseBins);
% $$$ mc = hsv(100);
% $$$ for m = {'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'}
% $$$     scatter(xyz{s}(statePeriods(i,[1,1])+tsubset,m{1},1),...
% $$$             xyz{s}(statePeriods(i,[1,1])+tsubset,m{1},2),...
% $$$             20,...
% $$$             mc(phaseBin(statePeriods(i,1)+[tsubset(1):tsubset(2)]),:),...
% $$$             'filled')
% $$$     scatter(xyz{s}(statePeriods(i,[1,1])+tsubset(1),m{1},1),...
% $$$             xyz{s}(statePeriods(i,[1,1])+tsubset(1),m{1},2),...
% $$$             30,...
% $$$             'k',...
% $$$             'filled')
% $$$ 
% $$$ end
% $$$ cax = colorbar;
% $$$ colormap(mc)
% $$$ daspect([1,1,1])
% $$$ cax.Ticks = linspace(0,1,5);
% $$$ cax.TickLabels = {'-180','-90','0','90','180'};
% $$$ 
% $$$ RefTrialName = [sessionList(referenceSessionIndex).sessionName,'.',...
% $$$                 sessionList(referenceSessionIndex).mazeName,'.',...
% $$$                 sessionList(referenceSessionIndex).trialName];
% $$$ TrialName = [sessionList(s).sessionName,'.',sessionList(s).mazeName,'.',sessionList(s).trialName];
% $$$ FigName = ['SVD_walk_sway_fet_skeleton_ALT_REF_',RefTrialName,'_TARGET_',TrialName];
% $$$ print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
% $$$ print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));







rhm = cf(@(t) fet_rhm(t),Trials);
cf(@(h,x) set(h,'sync',x.sync.copy),rhm,xyz); cf(@(h,x) set(h,'origin',x.origin) ,rhm,xyz);
prhm = cf(@(a) a.phase([6,12]),rhm);
cf(@(h,x) set(h,'sync',x.sync.copy),prhm,xyz); cf(@(h,x) set(h,'origin',x.origin) ,prhm,xyz);

edp = linspace(-pi,pi,50); edr = linspace(-pi/2,pi/2,50);
outPhase         = cf(@(p,f,s) hist2([p(s{'w&a'},1),f(s{'w&a'},3)],edp,edp), prhm, fetWsvdPhase, Stc);
outHang          = cf(@(p,f,s) hist2([p(s{'w&a'},1),f(s{'w&a'},3)],edr,edp), rhm,  fetWsvdPhase, Stc);
%outPhaseFiltered = cf(@(p,f,s) hist2([p(s{'w&a'},2),f(s{'w&a'},3)],edp,edp), phang, fetWsvdPhase, Stc);
%outHangFiltered  = cf(@(p,f,s) hist2([p(s{'w&a'},2),f(s{'w&a'},3)],edr,edp), hang,  fetWsvdPhase, Stc);





hang = cf(@(a) MTADang('data',[circ_dist(circshift(a(:,'head_back','head_front',1),-10),...
                                        circshift(a(:,'head_back', 'head_front',1),10)),...
                               circ_dist(a(:,'head_back','head_front',1),...
                                         a(:,'fsl','fbcom',1)),...
                               circ_dist(a(:,'head_back','head_front',1),...
                                         a(:,'spine_lower','spine_upper',1))],...
                       'sampleRate',a.sampleRate),ang);


cf(@(h,x) set(h,'sync',x.sync.copy),hang,xyz); cf(@(h,x) set(h,'origin',x.origin) ,hang,xyz);
phang = cf(@(a) a.phase([1,5]),hang);
cf(@(h,x) set(h,'sync',x.sync.copy),phang,xyz); cf(@(h,x) set(h,'origin',x.origin) ,phang,xyz);

edp = linspace(-pi,pi,50); edr = linspace(-pi/2,pi/2,50);
outPhase         = cf(@(p,f,s) hist2([p(s{'w&a'},3),f(s{'w&a'},3)],edp,edp), phang, fetWsvdPhase, Stc);
outHang          = cf(@(p,f,s) hist2([p(s{'w&a'},3),f(s{'w&a'},3)],edr,edp), hang,  fetWsvdPhase, Stc);
outPhaseFiltered = cf(@(p,f,s) hist2([p(s{'w&a'},2),f(s{'w&a'},3)],edp,edp), phang, fetWsvdPhase, Stc);
outHangFiltered  = cf(@(p,f,s) hist2([p(s{'w&a'},2),f(s{'w&a'},3)],edr,edp), hang,  fetWsvdPhase, Stc);


% Display for each animal
% $$$ figure,
% $$$ for s = 1:numSessions,
% $$$     subplot(1,numSessions,s);
% $$$     imagesc(linspace(-pi/2,pi/2,50),linspace(-pi,pi,50),sum(cat(3,outp{s}),3)')
% $$$ end

%figure, for s = 1:numSessions, subplot(1,numSessions,s); imagesc(out{s}'); end

% DISPLAY relationship between lateral spine sway (lss) and head-body yaw
sp = [];
hfig = figure(gen_figure_id);
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
hfig.Position = [1,1,16,16];
sp(end+1) = subplot(221); imagesc(edr,edp,sum(cat(3,outHang{:})         ,3)');
sp(end+1) = subplot(222); imagesc(edp,edp,sum(cat(3,outPhase{:})        ,3)');
sp(end+1) = subplot(223); imagesc(edr,edp,sum(cat(3,outHangFiltered{:}) ,3)');
sp(end+1) = subplot(224); imagesc(edp,edp,sum(cat(3,outPhaseFiltered{:}),3)');
% ADD colorbar
hcax = af(@(a) colorbar(a),sp); cf(@(a) set(a,'Position',a.Position+[0.1,0,0,0]),hcax);
% caxis units are count
% x&y axis units are radians
suptitle('lateral spine sway and head-body yaw; row1 - unfiltered body angle; row2 - filtered')
FigName = ['LSS_head_body_' sessionListTag];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));



% DISPLAY step detection 
s = 1
stepFeature = ButFilter(afetW{s}(:,4),3,[1,6]./(afetW{s}.sampleRate/2),'bandpass');
[steps,swayMagnitude] = LocalMinima(-abs(stepFeature),0,-4);
states = {'walk','turn','pause','rear'};
sclr = 'bgcr';
sp = [];
figure,
sp(end+1) = subplot2(10,4,1:8,1:4); hold on
plot(ts{s},stepFeature);
plot(ts{s}(steps),stepFeature(steps),'or')
plot(ts{s},-afetW{s}(:,1)./2+20);
plot(ts{s},-reshape(walkFetRot{s},[],size(wfet{s},2))*mean(reshape(Vww(:,1),[],size(wfet{s},2)))'.*20+20);
plot(ts{s},-abs(afetW{s}(:,2))./2-20);
plot(ts{s},-abs(reshape(walkFetRot{s},[],size(wfet{s},2))*mean(reshape(Vww(:,2),[],size(wfet{s},2)))').*20-20);
Lines([],+20,'k');
legend({'Step Feature','Putative Steps','PC1 Forward displacement','PC1 raw','PC2 turning','PC2 raw'})
sp(end+1)=subplot2(10,4,[9],[1:4]);
%plotSTC(StcNN{s},1,'text',states,sclr,[],false);
plotSTC(Stc{s},1,'text',states,sclr,[],false);
sp(end+1)=subplot2(10,4,[10],[1:4]);
plotSTC(StcNN{s},1,'text',states,sclr,[],false);
linkaxes(sp,'x');



stepsWalk = [];
stepsDistWalk = [];
stepsDispWalk = [];
stepsHeadPitch = [];
stepsPC1 = [];
% $$$ stepsWalkR = {};
% $$$ stepsWalkL = {};
markers = {'spine_lower','pelvis_root','spine_middle','spine_upper','fbcom','hcom','fhcom'};
for s = 1:numSessions,
    %stepFeature = ButFilter(fetWsvd{s}(:,3),3,[1,6]./(fetWsvd{s}.sampleRate/2),'bandpass');
    stepFeature = ButFilter(afetW{s}(:,3),3,[1,6]./(afetW{s}.sampleRate/2),'bandpass');
    [steps,swayMagnitude] = LocalMinima(-abs(stepFeature),0,-4);
% $$$     [stepsL,swayMagnitudeL] = LocalMinima(-stepFeature,0,-4);
% $$$     [stepsR,swayMagnitudeR] = LocalMinima( stepFeature,0,-4);

%   SELECT walk peirods with a duration greater that 0.5 seconds
    statePeriods = [Stc{s}{'w&a'}.data];
    statePeriods(diff(statePeriods,1,2)<1*xyz{s}.sampleRate,:) = [];
%   SELECT steps with in walk peirods
    stepsWalkTemp = SelectPeriods(steps,statePeriods,'d',1,0);
    stepsWalk = [stepsWalk;stepsWalkTemp];
% $$$     stepsWalkR{s} = SelectPeriods(stepsR,statePeriods,'d',1,0);
% $$$     stepsWalkL{s} = SelectPeriods(stepsL,statePeriods,'d',1,0);

%   COMPUTE Distance of each inter step interval 
    %stepsDistTurn = sqrt(sum((xyz{s}(stepsTurn-10,:,[1,2])-xyz{s}(stepsTurn+10,:,[1,2])).^2,3));
    stepsDispWalk = cat(1,stepsDispWalk,...
                        sqrt(sum((xyz{s}(circshift(stepsWalkTemp,-1),markers,[1,2])-...
                                  xyz{s}(stepsWalkTemp,markers,[1,2])).^2,3)));
    stepsDistWalk = cat(1,stepsDistWalk,sqrt(sum((xyz{s}(stepsWalkTemp-20,markers,[1,2])-xyz{s}(stepsWalkTemp+20,markers,[1,2])).^2,3))./40*xyz{s}.sampleRate/10);
    stepsHeadPitch = cat(1,stepsHeadPitch,ang{s}(stepsWalkTemp,'head_back','head_front',2));
    stepsPC1 = cat(1,stepsPC1,rfetWsvd{s}(stepsWalkTemp,1));
end




% COMPUTE inter step interval 
stepsWalkISI = circshift(stepsWalk,-1)-stepsWalk;

ind = stepsWalkISI<60&stepsWalkISI>0;

abins = linspace(-pi/2,0.8,100);
[~,hpc] = histc(stepsHeadPitch,abins);
ac = jet(numel(abins));

tss = (stepsWalkISI(ind)+randn([sum(ind==1),1])/4)/sampleRate{1};

figure,
subplot(331);
%plot(tss,lstepsDistWalk(ind,1),'.b'); 
plot(tss,log10(stepsDistWalk(ind,1)),'.b'); 
%plot(stepsDistWalk(ind,1),stepsPC1(ind,1),'.b'); 
%plot(tss,stepsPC1(ind,1),'.b'); 
ylabel('swing speed (cm/s)')
xlabel('time (s)')

subplot(332);
plot(tss,stepsDispWalk(ind,1),'.b');
ylabel('step distance (mm)')
xlabel('time (s)')

subplot(333);
plot(log10(stepsDispWalk(ind,1)),log10(stepsDistWalk(ind,1)),'.b');xlim([1,2.2]),ylim([0.5,2.2])
ylabel('step distance (mm)')
xlabel('swing speed (cm/s)')

subplot(334)
hist2([tss,log10(stepsDistWalk(ind,1))],linspace(0,0.5,30),linspace(0.5,1.8,30));
subplot(335)
hist2([tss,log10(stepsDispWalk(ind,1))],linspace(0,0.5,30),linspace(1,2.2,30));
subplot(336)
hist2([log10(stepsDispWalk(ind,1)),log10(stepsDistWalk(ind,1))],linspace(1,2.2,30),linspace(0.5,1.8,30));

sp = [];
sp(end+1) = subplot(337);
plot3(tss, log10(stepsDistWalk(ind,1)), stepsHeadPitch(ind),'.')
sp(end+1) = subplot(338);
plot3(tss, log10(stepsDispWalk(ind,1)), stepsHeadPitch(ind),'.')
sp(end+1) = subplot(339);
plot3(log10(stepsDistWalk(ind,1)),...
      log10(stepsDispWalk(ind,1)),...
      stepsHeadPitch(ind),'.')

linkprop(sp,'View');




figure,
subplot(121);
scatter(stepsWalkISI(ind)+randn([sum(ind==1),1])*2, ...
        log10(stepsDistWalk(ind,1)),10,...
        ac(hpc(ind),:),'filled')
xlim([10,64]),ylim([0,2.2],

subplot(122);
plot(stepsWalkISI(ind)+randn([sum(ind==1),1])*2,log10(stepsDispWalk(ind,1)),'.b'); xlim([10,64]),ylim([0,2.2])

figure,
plot(stepsWalkISI(ind),stepsDistWalk(ind,1),'.b'); xlim([10,64]),

figure,
plot(log10(stepsDispWalk(ind,1)),log10(stepsDistWalk(ind,1)),'.b');xlim([0,2.2]),ylim([0,2.2])





% PAD ends with 
stepsW = swayMagnitude<-4;
stepsW = stepsW...
         |ismember([stepsW,circshift(stepsW,-1)],[0,1],'rows')...
         |ismember([stepsW,circshift(stepsW,1)],[0,1],'rows');
stepsW = [1;steps(find(stepsW));size(xyz{s},1)];



figure,
fetIndex = 1;
transWindow = 0.5;
tsts = {'walk','turn','pause','rear','groom','sit'};
nsts = numel(tsts);
sp = [];
for t = 1:nsts,
    for o = 1:nsts, 
        sp(t,o) = subplot2(nsts,nsts,t,o);
        hold('on');
    end
end
for s = 1:numSessions,
    for t = 1:nsts,
        for o = 1:nsts,
            if t~=o,
                axes(sp(t,o));
                %   SELECT state to state transitions
                wp = Stc{s}.get_state_transitions(Trials{s},{tsts{t},tsts{o}},transWindow,xyz{s});
                wp(diff(wp,1,2)<round(transWindow*xyz{s}.sampleRate),:)=[];
                if ~isempty(wp)
                    %fws = abs(reshape(fetWsvd{s}(wp,fetIndex),round(transWindow*xyz{s}.sampleRate)+1,[]));
                    fws = reshape(fetWsvd{s}(wp,fetIndex),round(transWindow*xyz{s}.sampleRate)+1,[]);
                    plot(fws(:,:,1))
                    Lines([],0,'k');
                    ylim([-80,10]);
                    xlim([1,size(fws,1)]);
                    title([tsts{t},' -> ' tsts{o}]);
                end
            end
        end
    end
end
suptitle(['stateTransition_PC' num2str(fetIndex)]);
FigName = ['stateTransition_PC' num2str(fetIndex)];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));



figure,
axes;hold on;
transWindow = 0.8;
fetIndex = 2;
tsts = {'walk','turn','pause','rear','groom','sit'};
t = 3;
o = 2;
cfws = [];
for s = 1:numSessions,
    wp = Stc{s}.get_state_transitions(Trials{s},{tsts{t},tsts{o}},transWindow,xyz{s});
    %wp = StcNN{s}.get_state_transitions(Trials{s},{tsts{t},tsts{o}},transWindow,xyz{s});
    %wp = bsxfun(@plus,StcNN{s}{'n'}(:,1),[-60,60]);
%wp(diff(wp,1,2)<round(transWindow*xyz{s}.sampleRate),:)=[];
if ~isempty(wp)
    %fws = abs(reshape(fetWsvd{s}(wp,fetIndex),round(transWindow*xyz{s}.sampleRate)+1,[]));
    fws = abs(reshape(afetW{s}(wp,fetIndex),round(transWindow*xyz{s}.sampleRate)+1,[]));
    fws = abs(reshape(cang{s}(wp,1),round(transWindow*xyz{s}.sampleRate)+1,[]));
    %fws = diff(reshape(rfetWsvd{s}(wp,fetIndex),round(transWindow*xyz{s}.sampleRate)+1,[]));
    cfws = cat(2,cfws,fws);
    plot(fws(:,:,1))
    Lines([],0,'k');
    xlim([1,size(fws,1)]);
    title([tsts{t},' -> ' tsts{o}]);
end
end
axis tight


% COMPUTE Filtered xyz
fxyz = cf(@(x) x.copy(), xyz);
cf(@(x) x.filter('ButFilter',3,1,'low'),fxyz);
% COMPUTE Filtered ang
fang = cf(@(t,x) create(MTADang,t,x), Trials, fxyz);
for s = 1:numSessions, ang{s}.data(~nniz(xyz{s}),:,:,:) = 0;end


cang = cf(@(a) a.copy('empty'), ang);
cf(@(c,a) set(c,'data',circ_dist(circshift(a(:,1,4,1),-10),circshift(a(:,1,4,1),10))),cang,fang);
cf(@(c,a) set(c,'data',circ_dist(circshift(a(:,3,7,1),-10),circshift(a(:,3,7,1),10))),cang,fang);

StcNN    = cf(@(Trial) Trial.load('stc','NN0317'), Trials);

s = 1;
figure,
transWindow = repmat({0.4},1,numSessions);
tsts = {'pause','walk','rear','groom'};
for t = 1:4
    stpa = cf(@(stc,Trial,state,tw,x) stc.get_state_transitions(Trial,{state,'turn'},tw,x),...
              Stc,Trials,repmat(tsts(t),1,numSessions),transWindow,xyz);
    stpa = cf(@(s,tw,sr) s(diff(s,1,2)==round(tw.*sr),:),stpa,transWindow,sampleRate);

    stpaNN = cf(@(stc,Trial,state,tw,x) stc.get_state_transitions(Trial,{state,'turn'},tw,x),...
              StcNN,Trials,repmat(tsts(t),1,numSessions),transWindow,xyz);
    stpaNN = cf(@(s,tw,sr) s(diff(s,1,2)==round(tw.*sr),:),stpaNN,transWindow,sampleRate);

    subplot(2,4,t);
    cfws = cf(@(a,s,tw,sr) reshape(a(s,2),round(tw*sr)+1,[]),afetW,   stpa,transWindow,sampleRate);
    %cfws = cf(@(a,s,tw,sr) reshape(a(s,2),round(tw*sr)+1,[]),rfetWsvd,stpa,transWindow,sampleRate);
    %cfws = cf(@(a,s,tw,sr) reshape(a(s,1),round(tw*sr)+1,[]),cang,    stpa,transWindow,sampleRate);
    [~,dind] = cf(@max,cf(@diff,cf(@abs,cfws)));
    %[~,dind] = cf(@max,cf(@abs,cfws));    
    %[~,dind] = cf(@max,cf(@diff,cf(@diff,cf(@abs,cfws))));
    hist(cat(2,dind{s}),24), axis tight
    title([tsts{t},' -> turn']);

    subplot(2,4,t+4);
    cfws = cf(@(a,s,tw,sr) reshape(a(s,2),round(tw*sr)+1,[]),afetW,   stpaNN,transWindow,sampleRate);
    %cfws = cf(@(a,s,tw,sr) reshape(a(s,2),round(tw*sr)+1,[]),rfetWsvd,stpaNN,transWindow,sampleRate);
    %cfws = cf(@(a,s,tw,sr) reshape(a(s,1),round(tw*sr)+1,[]),cang,    stpaNN,transWindow,sampleRate);
    [~,dind] = cf(@max,cf(@diff,cf(@abs,cfws)));
    %[~,dind] = cf(@max,cf(@abs,cfws));    
    %[~,dind] = cf(@max,cf(@diff,cf(@diff,cf(@abs,cfws))));
    hist(cat(2,dind{s}),24),  axis tight

% $$$     plot(bsxfun(@rdivide,abs(fws(:,:,1)),max(abs(fws(:,:,1)))))
% $$$     Lines([],0,'k');
% $$$     ylim([0,1]);
% $$$     xlim([1,size(fws,1)]);

end

%
figure,
transWindow = repmat({0.2},1,numSessions);
tsts = {'pause','turn','rear','groom'};
for t = 1:4
    subplot(1,4,t);
    stpa = cf(@(stc,Trial,state,tw,x) stc.get_state_transitions(Trial,{state,'walk'},tw,x),Stc,Trials,repmat(tsts(t),1,numSessions),transWindow,xyz);
    %stpa = cf(@(stc,Trial,state,tw,x) stc.get_state_transitions(Trial,{state,'walk'},tw,x),StcNN,Trials,repmat(tsts(t),1,numSessions),transWindow,xyz);
    stpa = cf(@(s,tw,sr) s(diff(s,1,2)==round(tw.*sr),:),stpa,transWindow,sampleRate);

    cfws = cf(@(a,s,tw,sr) -reshape(a(s,1),round(tw*sr)+1,[]),afetW,stpa,transWindow,sampleRate);
    %cfws = cf(@(a,s,tw,sr) -reshape(a(s,1),round(tw*sr)+1,[]),rfetWsvd,stpa,transWindow,sampleRate);
    %cfws = cf(@(a,s,tw,sr) reshape(a(s,1),round(tw*sr)+1,[]),cang,stpa,transWindow,sampleRate);
    %[~,dind] = cf(@max,cf(@diff,cf(@abs,cfws)));
    [~,dind] = cf(@max,cf(@diff,cfws));
    %[~,dind] = cf(@max,cf(@abs,cfws));
    %hist(cat(2,dind{:}),48)
    %figure
    hist(cat(2,dind{:}),24);

% $$$     plot(bsxfun(@rdivide,abs(fws(:,:,1)),max(abs(fws(:,:,1)))))
% $$$     Lines([],0,'k');
% $$$     ylim([0,1]);
% $$$     xlim([1,size(fws,1)]);
    title([tsts{t},' -> walk']);
end


% PLOT state epoch duration distributions
figure,
sts = 'wnprms';
nsts = numel(sts);
edx = -1:0.1:3;
for s = 1:nsts,
    subplot(1,nsts,s);
    aStsPer = cf(@(x) x{sts(s)}.data, Stc);
    aStsPer = cat(1,aStsPer{:});
    bar(edx,histc(log10(diff(aStsPer,1,2)./sampleRate{1}),edx),'histc')
    title(Stc{1}{sts(s)}.label);
    xlim([edx(1),edx(end)])
end

cummulativeIndex = cell2mat(cf(@size,xyz));
cummulativeIndex = cumsum(cummulativeIndex(1:3:end));
cummulativeIndex(end)=[];
cummulativeIndex = mat2cell([0,cummulativeIndex],1,ones([1,numel(xyz)]));


stsTag ={' onset',' offset'};
figure,
for s = 1:nsts;
    shl = cf(@(x,i) x{sts(s)}.data+i, Stc  ,cummulativeIndex); shl = cat(1,shl{:});%./sampleRate{1};
    snn = cf(@(x,i) x{sts(s)}.data+i, StcNN,cummulativeIndex); snn = cat(1,snn{:});%./sampleRate{1};
    for e = 1:2,
        [cxg,t,pairs] = CCG([shl(:,e);snn(:,e)],...                 Trains
                [ones([size(shl,1),1]);ones([size(snn,1),1])*2],... Groups
                12,              ... BinSize
                30,             ... HalfBins
                sampleRate{1},  ... sampleRate
                [1,2],             ... GSubset
                'count',        ... Normalization
                []              ... Epochs
        );
        subplot2(2,nsts,e,s);
        bar(t/1000,cxg(:,1,2))

        title({[Stc{1}{sts(s)}.label stsTag{e}],['HL Vs NN']});
    end
end




figure,
stsTag ={' onset',' offset'};
figure,
for s = 1:nsts;  
    shl = cf(@(x,i) x{sts(s)}.data+i, Stc  ,cummulativeIndex); shl = cat(1,shl{:});%./sampleRate{1};
    snn = cf(@(x,i) x{sts(s)}.data+i, StcNN,cummulativeIndex); snn = cat(1,snn{:});%./sampleRate{1};
    dshl = log10(diff(shl,1,2));
    dsnn = log10(diff(snn,1,2));
    [~,dshlBins] = histc(dshl,edx);
    [~,dsnnBins] = histc(dsnn,edx);
    for e = 1:2,
        [cxg,t,pairs] = CCG([shl(:,e);snn(:,e)],...                 Trains
                [dshlBins,dsnnBins+numel(edx)-1],...
        ...%[ones([size(shl,1),1]);ones([size(snn,1),1])*2],... Groups
                12,              ... BinSize
                30,             ... HalfBins
                sampleRate{1},  ... sampleRate
                [1,2],             ... GSubset
                'count',        ... Normalization
                []              ... Epochs
        );
        subplot2(2,nsts,e,s);
        bar(t/1000,cxg(:,1,2))

        title({[Stc{1}{sts(s)}.label stsTag{e}],['HL Vs NN']});
    end
end







% COMPUTE State transition matrix
figure,
transWindow = repmat({0.2},1,numSessions);
tsts = {'walk','pause','turn','rear','groom'};
for t = 1:4
    subplot(1,4,t);
    stpa{} = cf(@(stc,Trial,state,tw,x) stc.get_state_transitions(Trial,{state,'walk'},tw,x),Stc,Trials,repmat(tsts(t),1,numSessions),transWindow,xyz);



    title([tsts{t},' -> walk']);
end







transWindow = repmat({0.2},1,numSessions);
tsts = {'walk','rear','turn','pause','groom','sit'};
nsts = numel(tsts);
stpa = {};
for t = 1:nsts,
    for o = 1:nsts,
        if t~=o,
            %   SELECT state to state transitions
            stpa{t,o} = cell2mat(cf(@(stc,Trial,state,tw,x) stc.get_state_transitions(Trial,{tsts{t},tsts{o}},tw,x),Stc,Trials,repmat(tsts(t),1,numSessions),transWindow,xyz)');
        end
    end
end

stf = cell2mat(cf(@length,stpa));
round(bsxfun(@rdivide,stf,sum(stf)),2)
round(bsxfun(@rdivide,stf,sum(stf,2)),2)
round(stf./sum(stf(:)),4)

%         0    0.0224    0.0391   0.2591  0.0022  0.0088
%    0.0099         0    0.0177   0.0188  0.0011       0
%    0.0893    0.0011         0   0.0757  0.0011  0.0004
%    0.2157     0.019    0.1012        0  0.0237  0.0291
%    0.0041    0.0015    0.0019   0.0201       0       0
%     0.011    0.0004    0.0041   0.0207  0.0006       0



for t = 1:nsts
    for o = 1:nsts
        
    end
end


%---------------------------------------------------------------

% $$$ figure
% $$$ %hist2([swayPhase(ind,4),circshift(hswPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$ %hist2([swayPhase(ind,4),circshift(rhmExpPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$ hist2([swayPhase(ind,4),circshift(rhmOriPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$ 
% $$$ 
% $$$ % NCP 
% $$$ figure,
% $$$ subplot(131);
% $$$ ind = Stc{s}{'w'};
% $$$ hist2([swayPhase(ind,4),circshift(rhmPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$ subplot(132);
% $$$ hist2([swayPhase(ind,4),circshift(ncpPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$ subplot(133);
% $$$ hist2([rhmPhase(ind,1),circshift(ncpPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$ 
% $$$ 
% $$$ % LFP 
% $$$ chan = 3;
% $$$ shift = 0;
% $$$ figure,
% $$$ for chan = 1:8,
% $$$     %subplot2(8,3,chan,1);
% $$$ ind = Stc{s}{'w'};
% $$$ hist2([swayPhase(ind,3),circshift(rhmPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$ subplot2(8,3,chan,2);
% $$$ %hist2([swayPhase(ind,3),circshift(lfpPhase(ind,chan),shift)],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$ %subplot2(8,3,chan,3);
% $$$ hist2([rhmPhase(ind,1),circshift(lfpPhase(ind,chan),shift)],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ % LFP 
% $$$ chan = 3;
% $$$ shift = -180:10:-20;
% $$$ shift = -150:30:150;
% $$$ 
% $$$ ind = Stc{s}{'w'};
% $$$ figure,
% $$$ for chan = 1:8,
% $$$     for shiftInd = 1:numel(shift)
% $$$         subplot2(8,numel(shift),chan,shiftInd);
% $$$ %hist2([swayPhase(ind,3),circshift(rhmPhase(ind),shift(shiftInd))],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$ %hist2([swayPhase(ind,3),circshift(lfpPhase(ind,chan),shift(shiftInd))],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$ hist2([rhmPhase(ind,1),circshift(lfpPhase(ind,chan),shift(shiftInd))],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$     end
% $$$ end
% $$$ 
% $$$ figure,chan = 5;shiftInd = 6;
% $$$ %hist2([swayPhase(ind,4),circshift(rhmPhase(ind),shift(shiftInd))],linspace(-pi,pi,30),linspace(-pi,pi,30))
% $$$ hist2([rhmPhase(ind,1),circshift(lfpPhase(ind,chan),shift(shiftInd))],linspace(-pi,pi,30),linspace(-pi,pi,30));
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ rhmlpf = rhm.copy;
% $$$ rhmlpf.filter('ButFilter',3,5,'low');
% $$$ 
% $$$ ind = Stc{s}{'w'};
% $$$ figure,plot(rhmlpf(ind),fetW{s}(ind,3),'.');
% $$$ 
% $$$ figure,plot(nunity(rhmlpf(:)))
% $$$ hold on,plot(fetW{s}(:,3))




ind = Stc{1}{'n'};
ind.data(diff(ind.data,1,2)<20,:)=[]; 
[x,xlag] = xcorr(circ_dist(circshift(ang{1}(ind,5,7,1),-3),...
                           circshift(ang{1}(ind,5,7,1), 3)),...
                 afetW{1}(ind,2),120);
figure,plot(xlag/120,x);
[~,xm] = max(x);




