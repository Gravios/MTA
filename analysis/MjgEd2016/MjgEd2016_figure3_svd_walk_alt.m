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
%stcMode = 'NN0317';

% LOAD Sessions and basic data
Trials = af(@(Trial) MTATrial.validate(Trial)  , sessionList);
Stc    = cf(@(Trial)  Trial.load('stc')         , Trials);
StcNN  = cf(@(Trial) Trial.load('stc','NN0317'), Trials);
xyz    = cf(@(Trial) preproc_xyz(Trial)         , Trials);
% RESAMPLE to common sampleRate
cf(@(x,s) x.resample(s),xyz,sampleRate);
fxyz = cf(@(x) x.copy(), xyz);
cf(@(f) f.filter('ButFilter',5,[2.4],'low'),fxyz);

fet    = cf(@(Trial) fet_bref(Trial), Trials);
%fet    = cf(@(Trial) fet_bref(Trial,[],'LOAD_TRB_XYZ'), Trials);
cf(@(f,t) f.map_to_reference_session(t,'jg05-20120317.cof.all'),fet,Trials);
for s = 1:numSessions, fet{s}.data(~nniz(xyz{s}),:,:) = 0;end


zfrCat = cf(@(f) get(f,'data'),fet);
zfrCat = cat(1,zfrCat{:});
zfrMean = nanmean(zfrCat(nniz(zfrCat),:,:));
zfrStd = nanstd(zfrCat(nniz(zfrCat),:,:));
cf(@(w,m,s) set(w,'data',nunity(w,[],m,s)),...
            fet,...
            repmat({zfrMean},1,numSessions),...
            repmat({zfrStd},1,numSessions));
clear('zfrCat','zfrMean','zfrStd')


ffet = cf(@(f) f.copy, fet);
cf(@(f) f.filter('ButFilter',5,[1.2,10],'bandpass'),ffet);

% @wfs
% EMBED fet 


wfs  = cf(@(w,e) w.segs(1:size(w,1),e),fet,embeddingWindow);
wfs =  cf(@(w,e) circshift(w,e/2,2),wfs,embeddingWindow);
wfs =  cf(@(w,x) MTADxyz('data',reshape(permute(w,[2,1,3]),size(w,2),[]),...
              'sampleRate',x.sampleRate),wfs,xyz);
for i = 1:numel(wfs), wfs{i}.data(isnan(wfs{i}.data(:)))=0; end

% @svd
% DECOMPOSE fet with svd for walk and turn periods within all sessions
wfw = cf(@(w,s) w([s{'w+n'}]+[0.125,-0.125],:), wfs, Stc);
%wfw = cf(@(w,s) w([s{'w'}]+[0.5,-0.5],:), wfs, Stc);
[~,Sww,Vww] = svd(cat(1,wfw{:}),0);

% @afetW
% COMPUTE eigenvector loadings for each session's eigen vectors
afetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:20)),...
                          'sampleRate',w.sampleRate),...
           wfs,repmat({Vww},1,numSessions));
cf(@(f,x) set(f,'sync',x.sync.copy), afetW, xyz); cf(@(f,x) set(f,'origin',x.origin), afetW, xyz);

% @tfetW
% COMPUTE eigenvector loadings for each session's eigen vectors
tVww = reshape(Vww(:,2),[],size(fet{1},2));
tVww(:,[1:16,18:2:24,26:30]) = zeros;
tVww = reshape(tVww,[],1);
tfetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1)),...
                          'sampleRate',w.sampleRate),...
           wfs,repmat({tVww},1,numSessions));
cf(@(f,x) set(f,'sync',x.sync.copy), afetW, xyz); cf(@(f,x) set(f,'origin',x.origin), afetW, xyz);

% @ts
% TIME vectors
wts = cf(@(e,s)  [1:e]./s,                        embeddingWindow,sampleRate);
ts =  cf(@(x)    [1:size(x,1)]./x.sampleRate,     xyz);

% @ang
% COMPUTE  intermarker angles 
ang = cf(@(t,x) create(MTADang,t,x), Trials, fxyz);
for s = 1:numSessions, ang{s}.data(~nniz(xyz{s}),:,:,:) = 0;end


% END data processing ------------------------------------------------------------------





%% FIG.3
% DEF figure variables -----------------------------------------------------------------

% figure save paths
OwnDir = '/storage/gravio/ownCloud/';
%FigDir = 'MjgEd2016/figures/Figure_3';
FigDir = 'Shared/Behavior Paper/Figures/Figure_3/parts';

states = {'walk','rear','turn','pause','groom','sit'};
sclr = 'brgcmy';

s = 1;                           % jg05-20120317.cof.all
exampleTimePeriod = [2170,2198]; % seconds
%exampleTimePeriod = [2183,2191]; % seconds
exampleTimePeriodStr = num2str(exampleTimePeriod);
exampleTimePeriodStr(exampleTimePeriodStr==' ') = '_';

% END figure variables -----------------------------------------------------------------




% FIG3FETMAT ---------------------------------------------------------------------------
% project: MjgEd2016
% parent: figure 3
% subplots:
%    subplot 1: Selected timeperiod of feature matrix
%    subplot 2: Corresponds to subplot 1, contains state labels
% location: MjgEd2016_figure3_svd_walk_alt.m
%

hfig = figure(gen_figure_id);
hfig.Units = 'centimeters';
hfig.Position(3:4) = [15,10];
% subplot 1 - Feature Matrix
sp = subplot2(4,1,1:3,1); 
imagesc(ts{s},1:30,fet{s}(:,[1:2:9,2:2:10,11:15,16:2:24,17:2:25,26:30])');
axis xy
caxis([-5,5]);
sp(1).YTick = [2.5,7.5,12.5,17.5,22.5,27.5];
sp(1).YTickLabels = {'Dist2COM_F','Dist2COM_L','Dist2MAZE_Z',...
                    'dF/dt','dL/dt','dZ/dt'};
hcb = colorbar;
hcb.Position(1) = hcb.Position(1)+0.1;
% subplot 2 - State Labels
sp(2) = subplot2(4,1,4,1);
plotSTC(Stc{s},1,'text',states,sclr,[],false);
linkaxes(sp,'x');

% preprint formating
ForAllSubplots(['xlim([',num2str(exampleTimePeriod),'])'])

% save figure
TrialName = [sessionList(s).sessionName,'.',sessionList(s).mazeName,'.',sessionList(s).trialName];
FigName = ['featureMatrix_',TrialName,'_',exampleTimePeriodStr];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIG3FETMAT ---------------------------------------------------------------------------





% FIG3EIGVECT ------------------------------------------------------------------------------
% project: MjgEd2016
% parent: figure 3 
% subplots: 1->N eigen vectors for all sessions during walking periods
% location: MjgEd2016_figure3_svd_walk_alt.m

hfig = figure;
hfig.Units = 'centimeters';
hfig.Position(3:4) = [30,16];
hfig.PaperPositionMode = 'auto';
for i = 1:20,
    pc = -reshape(Vww(:,i),[],size(fet{1},2));
    pc = reshape(LR(:,i),[],size(fet{1},2));
    pc = pc(:,[1:2:9,2:2:10,11:15,16:2:24,17:2:25,26:30]);
    subplot(2,10,i);imagesc(wts{1},1:size(fet{1},2),pc'),
    caxis([-0.08,0.08]);
    %caxis([-.5,.5]);    
    axis xy
end

FigName = ['SVD_eigVect_walkturn'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIG3EIGVECT ---------------------------------------------------------------------------








% FIG3EIGVALUES -----------------------------------------------------------------------------
% project: MjgEd2016
% parent: figure 3
hfig = figure;
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
subplot(211),plot(log10(diag(Sww)));
subplot(212);plot(log10(diag(Sww)));
xlim([0,50])


FigName = ['SVD_eigVal_walkturn'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIG3EIGVALUES -------------------------------------------------------------------------


% FIG3EIGTS ---------------------------------------------------------------------------------
% DISPLAY eigen vectors with loadings and features
s = 1;
pcLabels = {'PC1','PC2','PC3','PC4','PC5'};
hfig = figure;
hfig.Units = 'centimeters';
hfig.Position(3:4) = [25,16];
hfig.PaperPositionMode = 'auto';
sp = [];
sp(end+1)=subplot2(10,4,[1:8],[1:4]);
hold on;
plot(ts{s},ffet{s}(:,17)*5+450,'g','LineWidth',1);   Lines([], 450,'k');
plot(ts{s},fet{s}(:,17)*5+350,'b','LineWidth',1);   Lines([], 350,'k');
plot(ts{s},ffet{s}(:,2)*10+250,'g','LineWidth',1);    Lines([], 250,'k');
plot(ts{s},fet{s}(:,2)*10+150,'b','LineWidth',1);    Lines([], 150,'k');
plot(ts{s},afetW{s}(:,1),'LineWidth',1);         Lines([],   0,'k');
plot(ts{s},afetW{s}(:,2)-100,'LineWidth',1);     Lines([],-100,'k');
plot(ts{s},afetW{s}(:,3)-200,'LineWidth',1);    Lines([],-200,'k');
plot(ts{s},afetW{s}(:,5)-300,'LineWidth',1);     Lines([],-300,'k');
plot(ts{s},afetW{s}(:,6)-300,'LineWidth',1);     Lines([],-300,'k');
plot(ts{s},afetW{s}(:,7)-400,'LineWidth',1);     Lines([],-400,'k');
legend({'FF17','F17','FF2','F2','PC1','PC2','PC3','PC5','PC6','PC7'})

sp(end+1)=subplot2(10,4,[9,10],[1:4]);
plotSTC(Stc{s},1,'text',states,sclr,[],false);
linkaxes(sp,'x');

xlim(sp(1),exampleTimePeriod)
xlim(sp(2),exampleTimePeriod)

TrialName = [sessionList(s).sessionName,'.',sessionList(s).mazeName,'.',sessionList(s).trialName];
FigName = ['SVD_walkturn_timeseries_',TrialName];
print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIG3EIGTS -----------------------------------------------------------------------------



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
% $$$     imagesc(wts{s},1:size(fet{s},2),reshape(rVw{s}(:,i/2),[],size(fet{s},2))');
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
% $$$     imagesc(wts{s},1:size(fet{s},2),reshape(Vw{referenceSessionIndex}(:,i/2),[],size(fet{s},2))');
% $$$     axis xy;
% $$$     caxis([-0.08,0.08]);
% $$$     ylabel(['PC Ref',num2str(i/2)])
% $$$     subplot2(10,5,i-1:i,2);
% $$$     imagesc(wts{s},1:size(fet{s},2),reshape(Vw{s}(:,i/2),[],size(fet{s},2))');
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




% RHM stuff -----------------------------------------------------------------------------

rhm = cf(@(t) fet_rhm(t),Trials);
cf(@(h,x) set(h,'sync',x.sync.copy),rhm,xyz); cf(@(h,x) set(h,'origin',x.origin) ,rhm,xyz);
prhm = cf(@(a) a.phase([6,12]),rhm);
cf(@(h,x) set(h,'sync',x.sync.copy),prhm,xyz); cf(@(h,x) set(h,'origin',x.origin) ,prhm,xyz);


pfet = cf(@(a) a.phase([1,12]),fet);
cf(@(h,x) set(h,'sync',x.sync.copy),pfet,xyz); cf(@(h,x) set(h,'origin',x.origin) ,pfet,xyz);

hang = cf(@(a) MTADang('data',[circ_dist(circshift(a(:,'head_back','head_front',1),-10),...
                                        circshift(a(:,'head_back', 'head_front',1),10)),...
                               circ_dist(a(:,'head_back','head_front',1),...
                                         a(:,'spine_lower','spine_upper',1))],...
                       'sampleRate',a.sampleRate),ang);


cf(@(h,x) set(h,'sync',x.sync.copy),hang,xyz); cf(@(h,x) set(h,'origin',x.origin) ,hang,xyz);
phang = cf(@(a) a.phase([1,5]),hang);
cf(@(h,x) set(h,'sync',x.sync.copy),phang,xyz); cf(@(h,x) set(h,'origin',x.origin) ,phang,xyz);

edp = linspace(-pi,pi,50); edr = linspace(-pi/2,pi/2,50);
outPhase         = cf(@(p,f,s) hist2([p(s{'w&a'},2),f(s{'w&a'},17)],edp,edp), phang, pfet, Stc);
outHang          = cf(@(p,f,s) hist2([p(s{'w&a'},2),f(s{'w&a'},17)],edr,edp), hang,  pfet, Stc);
%outPhaseFiltered = cf(@(p,f,s) hist2([p(s{'w&a'},2),f(s{'w&a'},3)],edp,edp), phang, pfet, Stc);
%outHangFiltered  = cf(@(p,f,s) hist2([p(s{'w&a'},2),f(s{'w&a'},3)],edr,edp), hang,  pfet, Stc);


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


% END RHM stuff -----------------------------------------------------------------------------




% FIG3STEPS ---------------------------------------------------------------------------------
% project: MjgEd2016
% parent: figure 3 
% subplots: total(3)
%     subplot 1: first two principle componets and step feature
%     overlaid with circles indicating steps.
%     subplot 2: display hand labels
%     subplot 3: display neural network labels
% location: MjgEd2016_figure3_svd_walk_alt.m

% DISPLAY step detection 
stepFeature = ffet{s}(:,17);
[steps,swayMagnitude] = LocalMinima(-abs(stepFeature),8,-1);
% $$$ stepFeature = ButFilter(afetW{s}(:,5),3,[1,10]./(afetW{s}.sampleRate/2),'bandpass');
% $$$ [steps,swayMagnitude] = LocalMinima(-abs(stepFeature),0,-4);
states = {'walk','turn','pause','rear'};
sclr = 'bgcr';
sp = [];
figure,
sp(end+1) = subplot2(10,4,1:8,1:4); hold on
plot(ts{s},stepFeature*4);
plot(ts{s}(steps),stepFeature(steps)*4,'or')
plot(ts{s},-afetW{s}(:,1)./2+20);
plot(ts{s},-abs(afetW{s}(:,2))./2-20);
Lines([],+20,'k');
legend({'Step Feature','Putative Steps','PC1 Forward displacement','PC2 turning'})
sp(end+1)=subplot2(10,4,[9],[1:4]);
%plotSTC(StcNN{s},1,'text',states,sclr,[],false);
plotSTC(Stc{s},1,'text',states,sclr,[],false);
sp(end+1)=subplot2(10,4,[10],[1:4]);
plotSTC(StcNN{s},1,'text',states,sclr,[],false);
linkaxes(sp,'x');

xlim(sp(1),exampleTimePeriod)
xlim(sp(2),exampleTimePeriod)

FigName = ['Step_detection_walk'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIG3STEPS ---------------------------------------------------------------------------






% COLLECT steps across sessions
stepsWalk = [];
stepsDistWalk = [];
stepsDispWalk = [];
stepsHeadPitch = [];
stepsPC1 = [];
% $$$ stepsWalkR = {};
% $$$ stepsWalkL = {};
markers = {'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'};
for s = 1:numSessions,
    %stepFeature = ButFilter(fetWsvd{s}(:,3),3,[1,6]./(fetWsvd{s}.sampleRate/2),'bandpass');
    stepFeature = ffet{s}(:,17);
    [steps,swayMagnitude] = LocalMinima(-abs(stepFeature),15,-1.5);
    %stepFeature = ButFilter(afetW{s}(:,3),3,[1,6]./(afetW{s}.sampleRate/2),'bandpass');
    %[steps,swayMagnitude] = LocalMinima(-abs(stepFeature),0,-4);
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
    stepsDistWalk = cat(1,stepsDistWalk,sqrt(sum((xyz{s}(stepsWalkTemp-10,markers,[1,2])-xyz{s}(stepsWalkTemp+10,markers,[1,2])).^2,3))./20*xyz{s}.sampleRate/10);
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



% FIG3STEPSSTATS ----------------------------------------------------------------------------
% project: MjgEd2016
% parent: figure 3 
% subplots: total(3)
%     subplot 1: first two principle componets and step feature
%     overlaid with circles indicating steps.
%     subplot 2: display hand labels
%     subplot 2: display neural network labels
% location: MjgEd2016_figure3_svd_walk_alt.m

figure,
subplot(231);hold on,       
%plot(tss,lstepsDistWalk(ind,1),'.b'); 
plot(tss,log10(stepsDistWalk(ind,1)),'.b'); 
[b,bs] = robustfit(tss,log10(stepsDistWalk(ind,1)),'bisquare')
t = [0.1,0.5];
plot(t,b(1)+b(2).*t,'r-')
%plot(stepsDistWalk(ind,1),stepsPC1(ind,1),'.b'); 
%plot(tss,log10(stepsPC1(ind,1)),'.b'); 
title({['corr coeff: ',num2str(bs.coeffcorr(2))],...
       ['b: ' num2str(b')]})
xlabel('Inter-Step Interval (s)')
ylabel('log10 swing speed (cm/s)')
xlim(t)
ylim([1,2])
 
subplot(232); hold on
plot(tss,log10(stepsDispWalk(ind,1)),'.b');
[b,bs] = robustfit(tss,log10(stepsDispWalk(ind,1)),'bisquare')
t = [0.1,0.5];
plot(t,b(1)+b(2).*t,'r-')
xlabel('Inter-Step Interval (s)')
ylabel('log10 step distance (mm)')
title({['corr coeff: ',num2str(bs.coeffcorr(2))],...
       ['b: ' num2str(b')]})
xlim(t)
ylim([1.2,2.2])

subplot(233);hold on
plot(log10(stepsDispWalk(ind,1)),log10(stepsDistWalk(ind,1)),'.b');xlim([1,2.2]),ylim([0.5,2.2])
[b,bs] = robustfit(log10(stepsDispWalk(ind,1)),log10(stepsDistWalk(ind,1)),'bisquare')
t = [1.2,2.2];
plot(t,b(1)+b(2).*t,'r-')
xlabel('log10 step distance (mm)')
ylabel('log10 swing speed (cm/s)')
title({['corr coeff: ',num2str(bs.coeffcorr(2))],...
       ['b: ' num2str(b')]})
xlim(t)
ylim([1,2])

subplot(234)
hist2([tss,log10(stepsDistWalk(ind,1))],linspace(0.1,0.5,30),linspace(1,2,30));
subplot(235)
hist2([tss,log10(stepsDispWalk(ind,1))],linspace(0.1,0.5,30),linspace(1.2,2.2,30));
subplot(236)
hist2([log10(stepsDispWalk(ind,1)),log10(stepsDistWalk(ind,1))],linspace(1,2.2,30),linspace(1,2,30));

FigName = ['Step_stats_all'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));


% END FIG3STEPSSTATS ----------------------------------------------------------------------------







% SECTION -3- State transition optimization
% Motivation: state label boundaries obtained from neural network
% classiifiers may vary with respect to annotations by the
% experimenter. 
% Description: 

% brainstorming
% Detect peaks of pricipal component scores (PCS) and or derivatives.
% compute JPDF of xcorr and pc amplitude peaks

% FIG3STSTRIGAVE --------------------------------------------------------------------------------

afetWFilt = cf(@(f) f.copy, afetW);
cf(@(f) set(f,'data',-f.data), afetW);
%cf(@(f) f.filter('ButFilter',5,[0.2,5],'bandpass'),afetWFilt);
cf(@(f) f.filter('ButFilter',5,[5],'low'),afetWFilt);
zfrCat = cf(@(f) get(f,'data'),afetWFilt);
zfrCat = cat(1,zfrCat{:});
zfrMean = nanmean(zfrCat(nniz(zfrCat),:,:));
zfrStd = nanstd(zfrCat(nniz(zfrCat),:,:));
cf(@(w,m,s) set(w,'data',nunity(w,[],m,s)),...
            afetWFilt,...
            repmat({zfrMean},1,numSessions),...
            repmat({zfrStd},1,numSessions));
clear('zfrCat','zfrMean','zfrStd')


fetInd = 1;
transitionStatePre = {'pause','pause+turn','turn','gper-turn'};
transitionStatePost = {'walk'};
shift = [0,round(0.5*afetW{1}.sampleRate)];
%shift = [round(-0.5*afetW{1}.sampleRate),0];

for pre = 1:numel(transitionStatePre),
    for p = 1:numel(transitionStatePost),
        sts = cf(@(s,w,t) [s.get_state_transitions(t,{transitionStatePre{pre},...
                            transitionStatePost{p}},[],w)],...
                 StcNN, afetWFilt, Trials);
        for s = 1:numSessions, 
            sts{s}(sum([sts{s}+repmat(shift,size(sts{s},1),1)<=0, ...
                       sts{s}+repmat(shift,size(sts{s},1),1)>size(afetW{s},1)],2)>0,:)=[];
        end

        
        
        wfw = cf(@(w,s,t) w.segs(round(mean(s,2))-round(2.*w.sampleRate),round(4.*w.sampleRate)),...
                 afetWFilt, sts, Trials);

        stsSind = cf(@(s,a) sign(circ_dist(...
                                 a(round(mean(s,2)+shift(1)),'spine_lower','spine_upper',1),...
                                 a(round(mean(s,2)+shift(2)),'spine_lower','spine_upper',1)))==1,...
                  sts,ang);

        
        onf = cat(2,wfw{:});
        fetSegs = onf(:,:,fetInd);

        segTime = linspace(-2,2,round(sampleRate{1}*4));
        ind = { cat(1,stsSind{:})==1, cat(1,stsSind{:})==0 };

        hfig = figure();
        % Mean turning trace
        subplot2(1,3,1,[1,2]);
        hold('on');
        for i = ind,
            plot(segTime,nanmean(fetSegs(:,i{:}),2))
            plot(segTime,nanmean(fetSegs(:,i{:}),2)+...
                 nanstd(fetSegs(:,i{:}),[],2)*2.576/sqrt(size(fetSegs,2)),'r')
            plot(segTime,nanmean(fetSegs(:,i{:}),2)-...
                 nanstd(fetSegs(:,i{:}),[],2)*2.576/sqrt(size(fetSegs,2)),'r')
% $$$             plot(segTime,nanmean(fetSegs(:,i{:}),2)+nanstd(fetSegs(:,i{:}),[],2)*2,'c')
% $$$             plot(segTime,nanmean(fetSegs(:,i{:}),2)-nanstd(fetSegs(:,i{:}),[],2)*2,'c')
        end
        title([transitionStatePre{pre} ,' -> ',transitionStatePost{p}]);
        ylabel('z-score')
        xlabel('time (s)')
        ylim([-5,5]);

        subplot2(1,3,1,3);
        hold('on');
        hout = cf(@(f,s,i) histc(f(s{'a-n-s-m-r'},i),-5:.2:5),afetWFilt,StcNN,repmat({fetInd},[1,numSessions]));
        hout = sum(cat(2,hout{:}),2);
        hax = barh(-5:.2:5,hout./sum(hout),'histc');
        set(hax, 'FaceAlpha',0.4,'FaceColor',[0,0,1],...
                 'EdgeAlpha',0.4,'EdgeColor',[0,0,1]);
        hout = cf(@(f,s,g,i,shift) histc(f([round(mean(s(g==1,:),2)),...
                                            round(mean(s(g==1,:),2))]+repmat(shift,sum(g==1),1),i),...
                                 -5:.2:5),...
                  afetWFilt,sts,stsSind,repmat({fetInd},[1,numSessions]),repmat({shift},[1,numSessions]));
        hout = sum(cat(2,hout{:}),2);
        hax = barh(-5:.2:5,hout./sum(hout),'histc');
        set(hax, 'FaceAlpha',0.4,'FaceColor',[0,1,0],...
                 'EdgeAlpha',0.4,'EdgeColor',[0,1,0]);
% $$$         hout = cf(@(f,s,g,i,sr) histc(f([round(mean(s(g==0,:),2)),round(mean(s(g==0,:),2)+0.5*sr)],i),...
% $$$                                  -5:.2:5),...
% $$$                   afetWFilt,sts,stsSind,repmat({fetInd},[1,numSessions]),sampleRate);
        hout = cf(@(f,s,g,i,shift) histc(f([round(mean(s(g==0,:),2)),...
                                            round(mean(s(g==0,:),2))]+repmat(shift,sum(g==0),1),i),...
                                 -5:.2:5),...
                  afetWFilt,sts,stsSind,repmat({fetInd},[1,numSessions]),repmat({shift},[1,numSessions]));
        hout = sum(cat(2,hout{:}),2);
        hax = barh(-5:.2:5,hout./sum(hout),'histc');
        set(hax, 'FaceAlpha',0.4,'FaceColor',[0,1,1],...
                 'EdgeAlpha',0.4,'EdgeColor',[0,1,1]);
        ylim([-5,5]);
        ylabel('z-score')
        xlabel('prob')

        FigName = ['State_transition_',transitionStatePre{pre},...
                   '_to_',transitionStatePost{p},'_PC',num2str(fetInd)];
        print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
        print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
    end
end


% END FIG3STSTRIGAVE --------------------------------------------------------------------------------



figure,plot(fetSegs)

figure;hold on;
nHalf = round(size(fetSegs,2)/2);
fetSegsSE = nanstd(nanmean(reshape(fetSegs(:,randi([1,size(fetSegs,2)],nHalf)),[size(fetSegs,1),nHalf,nHalf]),2),[],3)./sqrt(nHalf);
plot(nanmean(fetSegs,2))
plot(nanmean(fetSegs,2)+fetSegsSE*2.576,'r')
plot(nanmean(fetSegs,2)-fetSegsSE*2.576,'r')




plot(nanmean(fetSegs(:,2))


%fitresult = fit(onf(:,:,1),pop,'poly2')
con = confint(onf(:,:,1),0.95);

onf = diff(onf);

t = 2;
figure;hold on;
plot(nanmean(onf(:,:,t),2))
plot(nanmean(onf(:,:,t),2)+nanstd(onf(:,:,t),[],2)*2,'r')
plot(nanmean(onf(:,:,t),2)-nanstd(onf(:,:,t),[],2)*2,'r')


figure,plot(onf(:,:,2))
figure,plot(diff(onf(:,:,1)))

% END FIG3STSXCORR ------------------------------------------------------------------------------



% FIG3PC2FDF ------------------------------------------------------------------------------------


fafetW = cf(@(x) x.copy(), afetW);
cf(@(w) w.filter('ButFilter',3,2,'low'),fafetW);

hang = cf(@(a) MTADang('data',[circ_dist(circshift(a(:,'spine_lower','spine_upper',1),-10),...
                                        circshift(a(:,'spine_lower', 'spine_upper',1),10)),...
                               circ_dist(circshift(a(:,'spine_lower','head_front',1),-10),...
                                        circshift(a(:,'spine_lower', 'head_front',1),10))],...
                       'sampleRate',a.sampleRate),ang);


cf(@(w) w.filter('ButFilter',5,6,'low'),tfetW);
dtfetW = cf(@(x) x.copy(), tfetW);
cf(@(w,z) set(w,'data',[0;diff(z.data)]),dtfetW,tfetW);
cf(@(w) w.filter('ButFilter',5,6,'low'),dtfetW);
atfetW = cf(@(x) x.copy(), tfetW);
cf(@(w,z) set(w,'data',[diff(z.data);0]),atfetW,dtfetW);
cf(@(w) w.filter('ButFilter',5,6,'low'),atfetW);

figure();
hold('on');
ind = Stc{1}{'a-n-m-s-r'};
ind.cast('TimeSeries');
ind = find(ind.data==1);
ind = ind(randi(size(ind,1),5000,1));
scatter(abs(tfetW{1}(ind+30)),abs(dtfetW{1}(ind)),20,[0,0,1],'Filled');
ind = Stc{1}{'n'};
scatter(abs(tfetW{1}(ind(:,1)+30)),abs(dtfetW{1}(ind(:,1))),20,[0,1,0],'Filled');



% @wfs
% EMBED tfetW
tfet = cf(@(x) x.copy('empty'), xyz);
cf(@(w,x,y,z) set(w,'data',[x.data,y.data,z.data]),tfet,tfetW,dtfetW,atfetW);

zfrCat = cf(@(f) get(f,'data'),tfet);
zfrCat = cat(1,zfrCat{:});
zfrMean = nanmean(zfrCat(nniz(zfrCat),:,:));
zfrStd = nanstd(zfrCat(nniz(zfrCat),:,:));
cf(@(w,m,s) set(w,'data',nunity(w,[],m,s)),...
            tfet,...
            repmat({zfrMean},1,numSessions),...
            repmat({zfrStd},1,numSessions));
clear('zfrCat','zfrMean','zfrStd')

embeddingWindow = repmat({32},1,numSessions);
wft  = cf(@(w,e) w.segs(1:size(w,1),e),tfet,embeddingWindow);
wft =  cf(@(w,e) circshift(w,e/2,2),wft,embeddingWindow);
wft =  cf(@(w,x) MTADxyz('data',reshape(permute(w,[2,1,3]),size(w,2),[]),...
              'sampleRate',x.sampleRate),wft,xyz);
for i = 1:numel(wft), wft{i}.data(isnan(wft{i}.data(:)))=0; end


%tft = cf(@(w,s) w([s{'n'}(:,1)],:), wft, Stc);
tft = cf(@(w,s) w(bsxfun(@plus,[s{'n&a'}(:,1)],round([-0.1,0.1].*w.sampleRate)),:), wft, Stc);

[~,Stt,Vtt] = svd(cat(1,tft{:}),0);

figure,
for i = 1:20,
subplot(2,10,i);
imagesc(reshape(Vtt(:,i),[],size(tfet{1},2))'),axis xy
caxis([-0.1,0.1])
end
% @afetW
% COMPUTE eigenvector loadings for each session's eigen vectors
ntfetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:20)),...
                          'sampleRate',w.sampleRate),...
            wft,repmat({Vtt},1,numSessions));



i = 10;
tshift=[-60,60];
figure();
hold('on');
turnPer = Stc{1}{'n'}(i,:);
ind = turnPer+tshift;
plot(tfet{1}  (ind,1));
plot(ntfetW{1}(ind,1));
Lines([60,diff(turnPer)+60],[],'r');


figure();
hold('on');
plot(tfet{1}(:,1));
plot(ntfetW{1}(:,1),'g');
plot(hang{1}(:,1)*120,'c');
plot(hang{1}(:,2)*120,'m');

dhang = diff(hang{1}(:,:).*120,1,2).^2;
dhang(dhang<10)=10;
figure,plot(mean(hang{1}(:,:)*120,2)./dhang,'c');


Lines(Stc{1}{'n'}(:),[],'r');

figure
(ntfewW{1}(:,1)),Lines(Stc{1}('n')


s = 1;
figure,hold on
ind = Stc{1}{'a-n-m-s-r'};
plot(sum(ntfetW{s}(ind,1:3),2),afetW{s}(ind,2),'.')
ind = Stc{s}{'n'};
plot(sum(ntfetW{s}(ind,1:3),2),afetW{s}(ind,2),'.g')


figure();
hold('on');
ind = Stc{1}{'a-n-m-s-r'};
ind.cast('TimeSeries');
ind = find(ind.data==1);
ind = ind(randi(size(ind,1),10000,1));
scatter(abs(tfetW{1}(ind,1)),abs(ntfetW{1}(ind,1)),20,[0,0,1],'Filled');
ind = Stc{1}{'n'};
ind = round(mean(ind(:,:),2));
scatter(abs(tfetW{1}(ind,1)),abs(ntfetW{1}(ind,1)),20,[0,1,0],'Filled');


figure();
hold('on');
ind = Stc{1}{'a-n-m-s-r'};
[tfinds,tfmins] = LocalMinima(-abs(tfetW{1}(ind,1)),0,0);
%hf = hang{1}(ind,1);
hf = ntfetW{1}(ind,1);
scatter(-tfmins,abs(hf(tfinds)),20,[0,0,1],'Filled');
ind = Stc{1}{'n'};
ind = round(mean(ind(:,:),2));
scatter(abs(tfetW{1}(ind,1)),abs(ntfetW{1}(ind,1)),20,[0,1,0],'Filled');

figure();
hold('on');
ind = Stc{1}{'a-n-m-s-r'};
[tfinds,tfmins] = LocalMinima(-abs(ntfetW{1}(ind,1)),0,0);
tfinds([1:3,end-3:end]) = [];
%hf = hang{1}(ind,1);
tf = ntfetW{1}(ind,1);
hf = ntfetW{1}(ind,2);
scatter(tf(tfinds+40),abs(hf(tfinds)),20,[0,0,1],'Filled');
ind = Stc{1}{'n'};
ind = round(mean(ind(:,:),2));
scatter(ntfetW{1}(ind+40,1),abs(ntfetW{1}(ind,2)),20,[0,1,0],'Filled');


figure,
ind = Stc{s}{'n'};
plot(afetW{s}(ind,1),fntfetW{s}(ind,1),'.')

fntfetW = cf(@(x) x.copy(), ntfetW);
cf(@(w) w.filter('ButFilter',3,2,'low'),fntfetW);


ind = Stc{s}{'n'};
turnSegs = ntfetW{s}.segs(ind(:,1)-40,80,nan);



% Tau from peak
turnPeak = [];
turnTau  = [];
turnFwdM  = [];
turnDang  = [];
turnDur  = [];
f = 2;
for s = 1:numSessions,
    ind = Stc{s}{'n&a'};
    for i = 1:size(ind,1),
        fetChunkAmp = -abs(tfet{s}(ind(i,:)+[-30,30],1));
        fetChunkDer = -abs(diff([nan;tfet{s}(ind(i,:)+[-30,30],1)]).*5);
        %fetChunkAmp = -abs(fntfetW{s}(ind(i,:)+[-30,30],1));
        %fetChunkDer = -abs(fntfetW{s}(ind(i,:)+[-30,30],2));
        %fetChunkDer = -abs(diff([nan;fntfetW{s}(ind(i,:)+[-30,30],1)]).*5);
        [turnTauAmpTemp,turnPeakAmpTemp] = LocalMinima(fetChunkAmp,0,0);
        [turnTauDerTemp,turnPeakDerTemp] = LocalMinima(fetChunkDer,0,0);
        
        if isempty(turnTauAmpTemp),
            [turnPeakAmpTemp,turnTauAmpTemp] = min(fetChunkAmp);
        end        
        [~,mind] = min(turnTauAmpTemp);
        if fetChunkAmp(1)<turnPeakAmpTemp(mind),
            turnPeak(end+1) = fetChunkAmp(1);
            turnTau(end+1) = 1;
        else
            turnPeak(end+1) = turnPeakAmpTemp(mind);
            turnTau(end+1) = turnTauAmpTemp(mind);
        end
        turnFwdM(end+1) = fafetW{s}(turnTau(end)+ind(i,1),1);
        turnDur(end+1) = diff(ind(i,:));
        turnDang(end+1) = circ_dist(ang{s}(ind(i,1),'spine_lower','head_front',1),...
                                    ang{s}(ind(i,2),'spine_lower','head_front',1));
        %[turnPeak(s,i),turnTau(s,i)] = max(abs(fntfetW{s}(ind(i,:),1)));
    end
end

% $$$ turnPeak(turnPeak==0) = [];
% $$$ turnTau(turnTau==0) = [];

figure,plot(turnPeak,turnTau-30,'.');
figure,plot(turnPeak,abs(turnDang),'.');
figure,plot(turnDur,abs(turnDang),'.');
figure,plot(-turnPeak,turnDur,'.');
figure,plot(turnDur,turnTau,'.');

figure,hold on
hout = cf(@(n,s) histc(n(s{'a-n-r-s-m'},1),-80:2:80),ntfetW,Stc);
hax = bar(-80:2:80,sum(cat(2,hout{:}),2),'histc')
set(hax, 'FaceAlpha',0.4,'FaceColor',[0,0,1],...
         'EdgeAlpha',0.4,'EdgeColor',[0,0,1]);
hout = cf(@(n,s) histc(n(s{'n'},1),-80:2:80),ntfetW,Stc);
hax = bar(-80:2:80,sum(cat(2,hout{:}),2),'histc')
set(hax, 'FaceAlpha',0.4,'FaceColor',[0,1,0],...
         'EdgeAlpha',0.4,'EdgeColor',[0,1,0]);



% FIG3TURNFETMAT ---------------------------------------------------------------------------
% project: MjgEd2016
% parent: figure 3
% subplots:
%    subplot 1: Selected timeperiod of feature matrix
%    subplot 2: Corresponds to subplot 1, contains state labels
% location: MjgEd2016_figure3_svd_walk_alt.m
%

s = 1;
hfig = figure(gen_figure_id);
hfig.Units = 'centimeters';
hfig.Position(3:4) = [15,10];
% subplot 1 - Feature Matrix
sp = subplot2(4,1,1:3,1); 
plot(ts{s},fntfetW{s}(:,1:3));
axis xy
caxis([-5,5]);
hcb = colorbar;
hcb.Position(1) = hcb.Position(1)+0.1;
% subplot 2 - State Labels
sp(2) = subplot2(4,1,4,1);
plotSTC(Stc{s},1,'text',states,sclr,[],false);
linkaxes(sp,'x');

% preprint formating
ForAllSubplots(['xlim([',num2str(exampleTimePeriod),'])'])

% save figure
TrialName = [sessionList(s).sessionName,'.',sessionList(s).mazeName,'.',sessionList(s).trialName];
FigName = ['featureMatrix_',TrialName,'_',exampleTimePeriodStr];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIGTURN3FETMAT ---------------------------------------------------------------------------





% FIG3XCORRHIST  --------------------------------------------------------------------------------

% COLLECT the center timepoint of each turn
nper = cf(@(s) round(mean(s{'w'}.data,2)),Stc);
nper = cf(@(s) round(mean(s{'n'}.data,2)),Stc);
% COLLECT the local minima of the 'turn' feature
%[nmin,nval] = cf(@(f) LocalMinima(-abs(f(:,2)),20,-5),afetW);
[nmin,nval] = cf(@(f,s) LocalMinima(-abs(ButFilter(f(:,2),5,0.1/[s/2],'high')),20,-5),...
                 afetW,sampleRate);
[nmin,nval] = cf(@(f,s) LocalMinima(-abs([0;diff(ButFilter(f(:,2),5,0.1/[s/2],'high'))]),20,-1),...
                 afetW,sampleRate);

[~,nind] = cf(@(n,s) SelectPeriods(n,s{'n'},'d',1,0),nmin,StcNN);
inind        = cf(@(n,s) ismember(1:size(n,1),s),nval,nind)

hout = cf(@(n,i) histc(abs(n(i)),10:2:70),nval,inind);
figure,hold on
hax = bar(10:2:70,sum(cat(2,hout{:}),2),'histc')
set(hax, 'FaceAlpha',0.4,'FaceColor',[0,1,0],...
         'EdgeAlpha',0.4,'EdgeColor',[0,1,0]);

[~,sind] = cf(@(n,s) SelectPeriods(n,s{'r+m+s+n'},'d',0,0),nmin,StcNN);
isind        = cf(@(n,s) ~ismember(1:size(n,1),s),nval,sind)
sind = cf(@(s) find(s),isind);
txx = cf(@(s) sum(s),isind)
nmin        = cf(@(n,s) n(s),nmin,isind);
nval        = cf(@(n,s) n(s),nval,isind);

hout = cf(@(n,i) histc(abs(n),10:2:70),nval);
hax = bar(10:2:70,sum(cat(2,hout{:}),2),'histc')
set(hax, 'FaceAlpha',0.4,'FaceColor',[0,0,1],...
         'EdgeAlpha',0.4,'EdgeColor',[0,0,1]);



% DIAGNOSTIC 
ttx =ButFilter(afetW{1}(:,2),5,0.1/[sampleRate{1}/2],'high');
figure,hold on
plot(ttx);
plot(nmin{1},nval{1}.*-sign(ttx(nmin{1})),'or');
Lines(nper{1},[],'g');
Lines(Stc{1}{'n'}(:,1),[],'k');
Lines(Stc{1}{'n'}(:,2),[],'r');
plot(nunity(circ_dist(circshift(ang{1}(:,'spine_lower','spine_upper',1),-3),...
               circshift(ang{1}(:,'spine_lower','spine_upper',1),3))).*10)


tta = nunity(circ_dist(circshift(ang{1}(:,'spine_lower','spine_upper',1),-3),...
               circshift(ang{1}(:,'spine_lower','spine_upper',1),3)));
figure,
hold on
ind = nmin{1}(nind{1});
plot(ttx(ind),tta(ind),'.g')
ind = nmin{1}(isind{1});
plot(ttx(ind),tta(ind),'.b')



















% END FIG3XCORRHIST  ------------------------------------------------------------------------------








% PAD ends with bi
stepsW = swayMagnitude<-1;
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
transWindow = .8;
fetIndex = 2;
tsts = {'walk','turn','pause','rear','groom','sit'};
t = 3;
o = 2;
cfws = [];
for s = 1:numSessions,
    wp = Stc{s}.get_state_transitions(Trials{s},{tsts{t},tsts{o}},transWindow,xyz{s});
    %wp = StcNN{s}.get_state_transitions(Trials{s},{tsts{t},tsts{o}},transWindow,xyz{s});
    %wp = bsxfun(@plus,StcNN{s}{'n'}(:,1),[-60,60]);
    wp(diff(wp,1,2)<round(transWindow*xyz{s}.sampleRate/2),:)=[];
    wp = [wp(:,1),wp(:,1)+round(transWindow*xyz{s}.sampleRate)];

    wp(end,:) = [];
    if ~isempty(wp)
        %fws = abs(reshape(fetWsvd{s}(wp,fetIndex),round(transWindow*xyz{s}.sampleRate)+1,[]));
        fws = reshape(afetW{s}(wp,fetIndex),round(transWindow*xyz{s}.sampleRate)+1,[]);
% $$$         fws = bsxfun(@times,...
% $$$                      sign(afetW{s}(wp(:,2),fetIndex))',...
% $$$                      reshape(afetW{s}(wp,fetIndex),round(transWindow*xyz{s}.sampleRate)+1,[]));
        %fws = abs(reshape(cang{s}(wp,1),round(transWindow*xyz{s}.sampleRate)+1,[]));
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
% COMPUTE angular speed
cang = cf(@(a) a.copy('empty'), ang);
cf(@(c,a) set(c,'data',circ_dist(circshift(a(:,1,4,1),-10),circshift(a(:,1,4,1),10))),cang,fang);
cf(@(c,a) set(c,'data',circ_dist(circshift(a(:,3,7,1),-10),circshift(a(:,3,7,1),10))),cang,fang);

% LOAD neural network labeled states
StcNN    = cf(@(Trial) Trial.load('stc','NN0317'), Trials);

% FIG3XCORSTSTRANS ------------------------------------------------------------------------------------------
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

% END FIG3XCORSTSTRANS --------------------------------------------------------------------------------------


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
        axis tight
        ylim([0,500])
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
    stpa{} = cf(@(stc,Trial,state,tw,x) ...
                stc.get_state_transitions(Trial,{state,'walk'},tw,x),...
                Stc,Trials,repmat(tsts(t),1,numSessions),transWindow,xyz);

    title([tsts{t},' -> walk']);
end







transWindow = repmat({0.2},1,numSessions);
tsts = {'walk','rear','turn','pause','groom','sit'};
nsts = numel(tsts);
stpa = {};
for t = 1:nsts,
    for o = 1:nsts,
        if t~=o,
            % SELECT state to state transitions
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

% FIG3XCORR

ind = Stc{1}{'w+n'};
ind.data(diff(ind.data,1,2)<120,:)=[]; 
% $$$ [x,xlag] = xcorr(circ_dist(circshift(ang{1}(ind,5,7,1),-3),...
% $$$                            circshift(ang{1}(ind,5,7,1), 3)),...
% $$$                  afetW{1}(ind,2),120);
[x,xlag] = xcorr(circ_dist(circshift(ang{1}(ind,5,7,1),-3),...
                           circshift(ang{1}(ind,5,7,1), 3)),...
                 circ_dist(circshift(ang{1}(ind,1,4,1),-3),...
                           circshift(ang{1}(ind,1,4,1), 3)),...
                 120);
figure,plot(xlag/120,x);
[~,xm] = max(x);




sp = []
s = 1;
figure,
sp(end+1)=subplot2(10,4,[1:8],[1:4]);hold on
plot(ts{s},circ_dist(circshift(ang{s}(:,1,4,1),10),circshift(ang{s}(:,1,4,1),-10)))
plot(ts{s},circ_dist(circshift(ang{s}(:,5,7,1),10),circshift(ang{s}(:,5,7,1),-10)))
sp(end+1)=subplot2(10,4,[9,10],[1:4]);
plotSTC(Stc{s},1,'text',states,sclr,[],false);
linkaxes(sp,'x');

f = 16;
sp = []
s = 1;
figure,
sp(end+1)=subplot2(10,4,[1:8],[1:4]);hold on
plot(fet{s}(:,f))
plot(ffet{s}(:,f))
sp(end+1)=subplot2(10,4,[9,10],[1:4]);
plotSTC(Stc{s},sampleRate{s},'text',states,sclr,[],false);
linkaxes(sp,'x');



sclr = 'brgcmy';
figure
for f = 1:15,%size(fet,2),
    subplot(3,5,f);hold on
    eds = linspace([repmat(nanmean(fet{1}(Stc{1}{'a'},f)),[1,2])+...
                    repmat(nanstd(fet{1}(Stc{1}{'a'},f)),[1,2]).*[-10,10],300]);
    for s = 1:numSessions,
        hax = bar(eds,histc(fet{s}(Stc{s}{'w'},f),eds),'histc');
        hax.FaceColor = sclr(s);
        hax.EdgeColor = sclr(s);
        hax.FaceAlpha = 0.4;
        hax.EdgeAlpha = 0.4;
    end
    axis tight
end


rfet = cf(@(f) f.copy, fet);
rfet = cf(@(r,f) set(r,'data',permute(rms(f.segs(1:size(f),round(f.sampleRate/2),0)),[2,3,1])),rfet,afetW);



sclr = 'brgcmy';
figure, hold on
s = 1;
mf = rfet{s};
ms = Stc{s};
for f = 2;
    %subplot(3,5,f);hold on
    eds = linspace([repmat(nanmean(mf(ms{'a'},f)),[1,2])+...
                    repmat(nanstd(mf(ms{'a'},f)),[1,2]).*[-200,200],300]);
    for sts = 1:numSessions,
        hax = bar(eds,histc(mf(ms{states{sts}},f),eds),'histc');
        hax.FaceColor = sclr(sts);
        hax.EdgeColor = sclr(sts);
        hax.FaceAlpha = 0.4;
        hax.EdgeAlpha = 0.4;
    end
    axis tight
end





mfet = cf(@(a) a.copy(), fet);
cf(@(f,t) f.map_to_reference_session(t,'jg05-20120317.cof.all'),mfet,Trials);


figure,
x = 16;
y = 5;
edx = linspace(-30,30,100);
edy = linspace(-300,300,100);
stsk = 'amrwnps'
for s = 1:6,
    p = 1;
    for sts = stsk
        sp(s,p,1) = subplot2(6,7,s,p);
        ind = Stc{s}{sts};
        hist2([fet{s}(ind,x),fet{s}(ind,y)],edx,edy);
        title(ind.label)
        p=p+1;
    end
end
ForAllSubplots('caxis([0,100]),grid on')

figure,
for s = 1:6,
    p = 1;
    for sts = stsk
        sp(s,p,2) = subplot2(6,7,s,p);
        ind = Stc{s}{sts};
        hist2([mfet{s}(ind,x),mfet{s}(ind,y)],edx,edy);
        title(ind.label)
        p=p+1;
    end
end
linkaxes(sp(:),'x');
ForAllSubplots('caxis([0,100]),grid on')





% DISPLAY eigen vectors with loadings and features
s = 1;
states = {'walk','rear','turn','pause','groom','sit'};
sclr = 'brgcmy';
pcLabels = {'PC1','PC2','PC3','PC4','PC5'};
hfig = figure;
hfig.Units = 'centimeters';
hfig.Position(3:4) = [25,16];
hfig.PaperPositionMode = 'auto';
sp = [];
sp(end+1)=subplot2(10,4,[1:8],[2:4]);
hold on;
plot(ts{s},fet{s}(:,17).*20+2600,'LineWidth',1); % spine sway 'spine_lower'
plot(ts{s},fet{s}(:,19).*20+2700,'LineWidth',1); % spine sway 'pelvis_root'
plot(ts{s},fet{s}(:,21).*20+2800,'LineWidth',1); % spine sway 'spine_middle'
plot(ts{s},fet{s}(:,23).*20+2900,'LineWidth',1); % spine sway 'spine_upper'
plot(ts{s},efetW{s}(:,1),'r','LineWidth',1);           % Walk comp
plot(ts{s},efetW{s}(:,2),'b','LineWidth',1);           % Walk comp
plot(ts{s},efetW{s}(:,3),'c','LineWidth',1);           % Walk comp
plot(ts{s},efetW{s}(:,4)-2000,'g','LineWidth',1);           % Walk comp
plot(ts{s},efetW{s}(:,5)-2000,'m','LineWidth',1);           % Walk comp
plot(ts{s},efetW{s}(:,6)-2000,'r','LineWidth',1);           % Walk comp
legend({'F1','F2','F3','F4','PC1','PC2','PC4','PC6'})
Lines([],0,'k');
Lines([],5,'r');
Lines([],-5,'r');
sp(end+1)=subplot2(10,4,[9,10],[2:4]);
plotSTC(Stc{s},1,'text',states,sclr,[],false);
linkaxes(sp,'x');
for i = 1:1:6,
    sp(end+1)=subplot2(10,4,i,1);
    imagesc(wts{s},1:size(fet{s},2),reshape(LR(:,i),[],size(fet{s},2))');
    axis xy;
    caxis([-1,1]);
    ylabel(['PC',num2str(i)])
end












%% SCRAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ fet = cf(@(x) x.copy('empty'), xyz);
% $$$ cf(@(w,wfr) set(w,'data',[reshape(wfr,[],size(wfr,2)*size(wfr,3))]),fet,walkFetRot);
% $$$ wfs  = cf(@(w,e) w.segs([],e),fet,embeddingWindow);
% $$$ wfs =  cf(@(w,e) circshift(w,e/2,2),wfs,embeddingWindow);
% $$$ wfs =  cf(@(w,x) MTADxyz('data',reshape(permute(w,[2,1,3]),size(w,2),[]),...
% $$$               'sampleRate',x.sampleRate),wfs,xyz);
% $$$ for i = 1:numel(wfs), wfs{i}.data(isnan(wfs{i}.data(:)))=0; end

% DECOMPOSE fet with svd for walk and turn periods for each session
%[~,Sw,Vw] = cf(@(w,s) svd(w([s{'w+n'}]+[0.1,-0.1],:),0), wfs, Stc);

% $$$ % DECOMPOSE fet with svd for walk and turn periods within all sessions
% $$$ wfw = cf(@(w,s) w([s{'w+n'}]+[0.125,-0.125],:), wfs, Stc);
% $$$ %wfw = cf(@(w,s) w([s{'w'}]+[0.5,-0.5],:), wfs, Stc);
% $$$ [~,Sww,Vww] = svd(cat(1,wfw{:}),0);


% VARIMAX pca
% $$$ [LU,LR,FSr,VT] = erpPCA(cat(1,wfw{:}),20);
% $$$ 
% COMPUTE eigenvector loadings for each session's eigen vectors
efetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:20)),...
                          'sampleRate',w.sampleRate),...
           wfs,repmat({LR},1,numSessions));
cf(@(f,x) set(f,'sync',x.sync.copy), efetW, xyz); cf(@(f,x) set(f,'origin',x.origin), efetW, xyz);


% COMPUTE eigenvector loadings for each session's eigen vectors
% $$$ afetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:20)),...
% $$$                           'sampleRate',w.sampleRate),...
% $$$            wfs,repmat({Vww},1,numSessions));
% $$$ cf(@(f,x) set(f,'sync',x.sync.copy), afetW, xyz); cf(@(f,x) set(f,'origin',x.origin), afetW, xyz);

% COMPUTE eigenvector loadings for each session's eigen vectors
%fetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:20)),'sampleRate',w.sampleRate),wfs,Vw);
%cf(@(f,x) set(f,'sync',x.sync.copy), fetW, xyz); cf(@(f,x) set(f,'origin',x.origin), fetW, xyz);

% CREATE mask for eigenvectors for each session's eigen vectors
% $$$ maskEigVec = ones(size(Vww));
% $$$ for i = [1:16,48:64], maskEigVec(i:embeddingWindow{1}:end,:) = 0;end
% $$$ %for i = [1:32], maskEigVec(i:embeddingWindow{1}:end,:) = 0;end
% $$$ rVww = Vww.*maskEigVec;

% REDUCED eigenvector loadings for  each session's eigen vectors
% $$$ rfetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:20)),'sampleRate',w.sampleRate),...
% $$$            wfs,repmat({rVww},1,numSessions));
% $$$ cf(@(r,x) set(r,'sync',x.sync.copy), rfetW, xyz); cf(@(r,x) set(r,'origin',x.origin), rfetW, xyz);


% COMPUTE features with reference session ----------------------------------------------------------
% REFERENCED compute eigenvector loadings
% $$$ fetWsvd = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:20),[1,2],[1,2]),...
% $$$                             'sampleRate',w.sampleRate),...
% $$$              wfs,...
% $$$              repmat(Vww,1,numSessions));
% $$$ cf(@(f,x) set(f,'sync',x.sync.copy),fetWsvd,xyz); cf(@(f,x) set(f,'origin',x.origin) ,fetWsvd,xyz);
% $$$ 
% $$$ fetWsvdPhase = cf(@(f) f.phase([1,5]),fetWsvd);
% $$$ cf(@(f,x) set(f,'sync',x.sync.copy),fetWsvdPhase,xyz); cf(@(f,x) set(f,'origin',x.origin) ,fetWsvdPhase,xyz);
% $$$ 
% $$$ % @rfetWsvd
% $$$ % REFERENCED reduced eigenvector loadings
% $$$ rfetWsvd = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:20),[1,2],[1,2]),...
% $$$                             'sampleRate',w.sampleRate),...
% $$$               wfs,...
% $$$               repmat({rVww},1,numSessions));
% $$$ cf(@(f,x) set(f,'sync',x.sync.copy),rfetWsvd,xyz);
% $$$ cf(@(f,x) set(f,'origin',x.origin) ,rfetWsvd,xyz);
% $$$ rfetWsvdPhase = cf(@(f) f.phase([1,5]),rfetWsvd);
% $$$ cf(@(f,x) set(f,'sync',x.sync.copy),rfetWsvdPhase,xyz);
% $$$ cf(@(f,x) set(f,'origin',x.origin) ,rfetWsvdPhase,xyz);
