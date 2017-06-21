
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
%sampleRate   = repmat({119.881035}            ,1,numSessions);
sampleRate   = repmat({60}                    ,1,numSessions);
featureName  = repmat({'MjgEd2016-F3-SVDREAR'},1,numSessions);
featureLabel = repmat({'fet_bref'}            ,1,numSessions);
featureKey   = repmat({'r'}                   ,1,numSessions);
embeddingWindow = repmat({60},1,numSessions);
%stcMode = 'NN0317';

% LOAD Sessions and basic data
Trials = af(@(Trial) MTATrial.validate(Trial), sessionList);
Stc    = cf(@(Trial) Trial.load('stc')       , Trials);
xyz    = cf(@(Trial) Trial.load('xyz')       , Trials);
% RESAMPLE to common sampleRate
cf(@(x,s) x.resample(s),xyz,sampleRate);
fet    = cf(@(Trial,x) fet_bref(Trial,x), Trials,xyz);
%fet    = cf(@(Trial) fet_bref(Trial,[],'LOAD_TRB_XYZ'), Trials);
cf(@(f,t) f.map_to_reference_session(t,'jg05-20120317.cof.all'),fet,Trials);


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


% DECOMPOSE features ------------------------------------------------------------------------------
% @wfs
% EMBED wfet 
wfet = cf(@(x) x.copy('empty'), xyz);
cf(@(w,z) set(w,'data',z.data),wfet,fet);
wfs  = cf(@(w,e) w.segs([],e),wfet,embeddingWindow);
wfs =  cf(@(w,e) circshift(w,e/2,2),wfs,embeddingWindow);
wfs =  cf(@(w,x) MTADxyz('data',reshape(permute(w,[2,1,3]),size(w,2),[]),...
              'sampleRate',x.sampleRate),wfs,xyz);
for i = 1:numel(wfs), wfs{i}.data(isnan(wfs{i}.data(:)))=0; end

% DECOMPOSE wfet with svd for walk and turn periods for each session
%[~,Sw,Vw] = cf(@(w,s) svd(w([s{'r'}]+[0.25,-0.25],:),0), wfs, Stc);

% DECOMPOSE wfet with svd for walk and turn periods within all sessions
wfw = cf(@(w,s) w([s{'r'}]+[0.25,-0.25],:), wfs, Stc);
[~,Sww,Vww] = svd(cat(1,wfw{:}),0);

% COMPUTE eigenvector loadings for each session's eigen vectors
afetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:10)),...
                          'sampleRate',w.sampleRate),...
           wfs,repmat({Vww},1,numSessions));
cf(@(f,x) set(f,'sync',x.sync.copy), afetW, xyz); cf(@(f,x) set(f,'origin',x.origin), afetW, xyz);


% DECOMPOSE wfet with svd for walk and turn periods within all sessions
% $$$ wfw = cf(@(w,s,r) w(bsxfun(@plus,[s{'r'}(:,1)],round([-0.2,0.2].*r)),:), wfs, Stc,sampleRate);
% $$$ [~,Son,Von] = svd(cat(1,wfw{:}),0);
% $$$ 
% $$$ % COMPUTE eigenvector loadings for each session's eigen vectors
% $$$ onfetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:10)),...
% $$$                            'sampleRate',w.sampleRate),...
% $$$            wfs,repmat({Von},1,numSessions));
% $$$ cf(@(f,x) set(f,'sync',x.sync.copy), onfetW, xyz); cf(@(f,x) set(f,'origin',x.origin), onfetW, xyz);
% $$$ 
% $$$ 
% $$$ 
% $$$ % COMPUTE eigenvector loadings for each session's eigen vectors
% $$$ fetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:10)),'sampleRate',w.sampleRate),wfs,Vw);
% $$$ cf(@(f,x) set(f,'sync',x.sync.copy), fetW, xyz); cf(@(f,x) set(f,'origin',x.origin), fetW, xyz);

% $$$ % CREATE mask for eigenvectors for each session's eigen vectors
% $$$ maskEigVec = ones(size(Vw{1}));
% $$$ for i = [1:16,48:64], maskEigVec(i:embeddingWindow{1}:end,:) = 0;end
% $$$ maskEigVec = repmat({maskEigVec},1,numSessions);
% $$$ rVw = cf(@(v,m) v.*m,Vw,maskEigVec);
% $$$ 
% $$$ % REDUCED eigenvector loadings for  each session's eigen vectors
% $$$ rfetW = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:10)),'sampleRate',w.sampleRate),wfs,rVw);
% $$$ cf(@(r,x) set(r,'sync',x.sync.copy), rfetW, xyz); cf(@(r,x) set(r,'origin',x.origin), rfetW, xyz);
% $$$ 
% $$$ 
% $$$ % COMPUTE features with reference session ----------------------------------------------------------
% $$$ referenceSessionIndex = 1;
% $$$ % REFERENCED compute eigenvector loadings
% $$$ fetWsvd = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:10),[1,2],[1,2]),...
% $$$                             'sampleRate',w.sampleRate),...
% $$$              wfs,...
% $$$              repmat(rVw(referenceSessionIndex),1,numSessions));
% $$$ cf(@(f,x) set(f,'sync',x.sync.copy),fetWsvd,xyz); cf(@(f,x) set(f,'origin',x.origin) ,fetWsvd,xyz);
% $$$ 
% $$$ fetWsvdPhase = cf(@(f) f.phase([1,5]),fetWsvd);
% $$$ cf(@(f,x) set(f,'sync',x.sync.copy),fetWsvdPhase,xyz); cf(@(f,x) set(f,'origin',x.origin) ,fetWsvdPhase,xyz);
% $$$ 
% $$$ % @rfetWsvd
% $$$ % REFERENCED reduced eigenvector loadings
% $$$ rfetWsvd = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:10),[1,2],[1,2]),...
% $$$                             'sampleRate',w.sampleRate),...
% $$$               wfs,...
% $$$               repmat(rVw(referenceSessionIndex),1,numSessions));
% $$$ cf(@(f,x) set(f,'sync',x.sync.copy),rfetWsvd,xyz);
% $$$ cf(@(f,x) set(f,'origin',x.origin) ,rfetWsvd,xyz);
% $$$ rfetWsvdPhase = cf(@(f) f.phase([1,5]),rfetWsvd);
% $$$ cf(@(f,x) set(f,'sync',x.sync.copy),rfetWsvdPhase,xyz);
% $$$ cf(@(f,x) set(f,'origin',x.origin) ,rfetWsvdPhase,xyz);

% @ts
% TIME vectors
wts = cf(@(e,s)  [1:e]./s,                        embeddingWindow,sampleRate);
ts =  cf(@(x)    [1:size(x,1)]./x.sampleRate,     xyz);

% @ang
% COMPUTE  intermarker angles 
ang = cf(@(t,x) create(MTADang,t,x), Trials, xyz);
for s = 1:numSessions, ang{s}.data(~nniz(xyz{s}),:,:,:) = 0;end



% DEF figure variables -----------------------------------------------------------------

% figure save paths
OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'Shared/Behavior Paper/Figures/Figure_3/parts';

states = {'walk','rear','turn','pause','groom','sit'};
sclr = 'brgcmy';

s = 1;                           % jg05-20120317.cof.all
exampleTimePeriod = [2170,2198]; % seconds
%exampleTimePeriod = [2183,2191]; % seconds
exampleTimePeriodStr = num2str(exampleTimePeriod);
exampleTimePeriodStr(exampleTimePeriodStr==' ') = '_';

% END figure variables -----------------------------------------------------------------




% FIG3EIGVECT ------------------------------------------------------------------------------
% project: MjgEd2016
% parent: figure 3 
% subplots: 1->N eigen vectors for all sessions during rearing periods
% location: MjgEd2016_figure3_svd_rear.m

hfig = figure;
hfig.Units = 'centimeters';
hfig.Position(3:4) = [30,16];
hfig.PaperPositionMode = 'auto';
for i = 1:40,
    pc = -reshape(Vww(:,i),[],size(wfet{1},2));
    pc = pc(:,[1:2:9,2:2:10,11:15,16:2:24,17:2:25,26:30]);
    subplot(4,10,i);imagesc(wts{1},1:size(wfet{1},2),pc'),
    caxis([-0.08,0.08]);
    axis xy
end

FigName = ['SVD_eigVect_rear'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIG3EIGVECT ---------------------------------------------------------------------------




% FIG3EIGVALUES -----------------------------------------------------------------------------
% project: MjgEd2016
% parent: figure 3
% subplots: 
%    subplot1: sorted eigen value
%    subplot2: first 50 sorted eigen value
% location: MjgEd2016_figure3_svd_rear.m

hfig = figure;
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
subplot(211),plot(log10(diag(Sww)));
subplot(212);plot(log10(diag(Sww)));
xlim([0,50])


FigName = ['SVD_eigVal_rear'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIG3EIGVALUES -------------------------------------------------------------------------





% FIG3EIGTS ---------------------------------------------------------------------------------
% project: MjgEd2016
% parent: figure 3
% description: eigen vectors with loadings and features
% subplots: 
%    subplot1: projections onto eigenvectors
%    subplot2: state labels
% location: MjgEd2016_figure3_svd_rear.m

s = 1;
states = {'walk','rear','turn','pause','groom','sit'};
sclr = 'brgcmy';
pcLabels = {'PC1','PC2','PC3','PC4','PC5'};
hfig = figure;
hfig.Units = 'centimeters';
hfig.Position(3:4) = [25,16];
hfig.PaperPositionMode = 'auto';
sp = [];
sp(end+1)=subplot2(10,4,[1:8],[1:4]);
hold on;
plot(ts{s},ffet{s}(:,30)*5+450,'g','LineWidth',1);   Lines([], 450,'k');
plot(ts{s},fet{s}(:,30)*5+350,'b','LineWidth',1);   Lines([], 350,'k');
plot(ts{s},ffet{s}(:,15)*10+250,'g','LineWidth',1);    Lines([], 250,'k');
plot(ts{s},fet{s}(:,15)*10+150,'b','LineWidth',1);    Lines([], 150,'k');
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
FigName = ['SVD_rear_timeseries_',TrialName];
print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIG3EIGTS -----------------------------------------------------------------------------




% FIG3STSXCORR ------------------------------------------------------------------------------
fetOnSet = cf(@(x) x.copy('empty'), xyz);
cf(@(nf,f,s) set(nf,'data',f(:,1:2)),fetOnSet,afetW,sampleRate);
%cf(@(nf,f,s) set(nf,'data',ButFilter(f(:,1:2),3,0.8/[s/2],'low')),fetOnSet,afetW,sampleRate);
wfw = cf(@(w,s,r) ...
         w.segs([s{'r',w.sampleRate}(:,2)]-round(2.*w.sampleRate),round(4.*w.sampleRate)),...
         fetOnSet, Stc);


onf = cat(2,wfw{:});
figure;hold on;
plot(nanmean(onf(:,:,1),2))
plot(nanmean(onf(:,:,1),2)+nanstd(onf(:,:,1),[],2)*2,'r')
plot(nanmean(onf(:,:,1),2)-nanstd(onf(:,:,1),[],2)*2,'r')

fitresult = fit(onf(:,:,1),pop,'poly2')
con = confint(onf(:,:,1),0.95);

onf = diff(onf);

t = 1;
figure;hold on;
plot(nanmean(onf(:,:,t),2))
plot(nanmean(onf(:,:,t),2)+nanstd(onf(:,:,t),[],2)*2,'r')
plot(nanmean(onf(:,:,t),2)-nanstd(onf(:,:,t),[],2)*2,'r')


figure,plot(onf(:,:,2))
figure,plot(diff(onf(:,:,1)))

% END FIG3STSXCORR ------------------------------------------------------------------------------




% FIG3FETPDF ---------------------------------------------------------------------------------
% project: MjgEd2016
% parent: figure 3
% description: 
% subplots: 
%    subplot1: projections onto eigenvectors
%    subplot2: state labels
% location: MjgEd2016_figure3_svd_rear.m

v = 1;
s = 1;
hfet = xyz{s}.copy;
%hfet.data = xyz{s}(:,'head_front',3)
hfet.data = ang{s}(:,'spine_middle','spine_upper',2);
figure,
edx = linspace(-pi/2,pi/2,100);
ind = Stc{s}{'r'}+[0.2,-0.2];
h = bar(edx,histc(hfet(ind,v),edx),'histc');
h.FaceColor = 'r';h.EdgeColor = 'r';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;
hold on
ind = Stc{s}{'a'}-(Stc{s}{'r'}+[0.2,-0.2]);
h = bar(edx,histc(hfet(ind,v),edx),'histc');
h.FaceColor = 'c';h.EdgeColor = 'c';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;

% END FIG3FETPDF------------------------------------------------------------------------------




% FIG3FETJPDF ---------------------------------------------------------------------------------
% project: MjgEd2016
% parent: figure 3
% description: 
% subplots: 
%    subplot1: 
%    subplot2: 
% location: MjgEd2016_figure3_svd_rear.m

s = 5
edx = linspace(-pi/2,pi/2,100);
edy = linspace(20,300,100);
hfet.data = [ang{s}(:,'spine_middle','spine_upper',2),xyz{s}(:,'head_back',3)];

%figure,
subplot(221);
ind = Stc{s}{'r'}+[0.2,-0.2];
hist2(hfet(ind,:),edx,edy);
subplot(223);
ind = Stc{s}{'a'}-(Stc{s}{'r'}+[0.2,-0.2]);
hist2(hfet(ind,:),edx,edy);
ForAllSubplots('caxis([0,10])')
subplot(222);
ind = Stc{s}{'p'}+[0.2,-0.2];
hist2(hfet(ind,:),edx,edy);
subplot(224);
ind = Stc{s}{'m'}+[0.2,-0.2];
hist2(hfet(ind,:),edx,edy);

% END FIG3FETJPDF------------------------------------------------------------------------------




s = 1;
figure,
hax = axes('Units','Centimeters',...
           'Position',[1,1,6,2.5],...
           'FontSize',8)
plot(ts{s},fet{s}(:,11:15),'LineWidth',1)


eds = linspace(-40,120,100)
sclr = 'brgcmy';
figure,hold on
for s = 1:numSessions,
    hax = bar(eds,histc(-afetW{s}(Stc{s}('r'),1),eds),'histc');
    hax.FaceColor = sclr(s);
    hax.EdgeColor = sclr(s);
    hax.FaceAlpha = 0.4;
    hax.EdgeAlpha = 0.4;
end



sclr = 'brgcmy';
figure
for f = 1:15,%size(fet,2),
    subplot(3,5,f);hold on
    eds = linspace([repmat(nanmean(fet{1}(Stc{1}{'a'},f)),[1,2])+...
                    repmat(nanstd(fet{1}(Stc{1}{'a'},f)),[1,2]).*[-10,10],300]);
    for s = 1:numSessions,
        hax = bar(eds,histc(fet{s}(Stc{s}{'m'},f),eds),'histc');
        hax.FaceColor = sclr(s);
        hax.EdgeColor = sclr(s);
        hax.FaceAlpha = 0.4;
        hax.EdgeAlpha = 0.4;
    end
    axis tight
end

