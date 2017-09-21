% MjgEd2016 Figure 4 Head Body kinematic organization
%
% A. Timeseries of head and body angular velocities
% B. Time lagged mutual information of head and body angular velocities
% C. 


sessionList = 'hand_labeled';
Trials = af(@(Trial) MTATrial.validate(Trial), get_session_list(sessionList));
stc = cf(@(t) t.load('stc'), Trials);
%stc = cf(@(t) t.load('stc','msnnN0+hand_labeled'), Trials);

s = 1
xyz = cf(@(Trial) Trial.load('xyz'),                         Trials);
vxy = cf(@(x)     x.vel({'spine_lower','head_front'},[1,2]), xyz   );
for s = 1:numel(Trials), vxy{s}.data(vxy{s}.data<1e-3,:) = 1e-3;end
for s = 1:numel(Trials), vxy{s}.data = log10(vxy{s}.data);end
cf(@(x) x.filter('ButFilter',3,5,'low'),xyz);;
ang = cf(@(t,x) create(MTADang,t,x),Trials,xyz);

ind = cf(@(stc) [stc{'a-s-m'}],stc);

%ind = cf(@(s,w,t) [s.get_state_transitions(t,{'pause+walk','rear'},1,w)],...
%         stc, xyz, Trials);

% $$$ shift = [0,30];
% $$$ for s = 1:numel(Trials), 
% $$$     ind{s}(sum([ind{s}+repmat(shift,size(ind{s},1),1)<=0, ...
% $$$                 ind{s}+repmat(shift,size(ind{s},1),1)>size(xyz{s},1)],2)>0,:)=[];
% $$$ end

s = 1;
tind = [277,283];
figure
sp = [];
pdh = cf(@(ang,ind) circ_dist(circshift(ang(ind,'head_back','head_front',2),-1),...
                              circshift(ang(ind,'head_back','head_front',2), 1)),...
         ang,ind);
adh = cf(@(ang,ind) circ_dist(circshift(ang(ind,'head_back','head_front',1),-1),...
                              circshift(ang(ind,'head_back','head_front',1), 1)),...
         ang,ind);

for marker = {'spine_lower','pelvis_root','spine_middle'},
    pdb = cf(@(ang,ind) circ_dist(circshift(ang(ind,marker,'spine_upper',2),-1),...
                                  circshift(ang(ind,marker,'spine_upper',2), 1)),...
             ang,ind);
    adb = cf(@(ang,ind) circ_dist(circshift(ang(ind,marker,'spine_upper',1),-1),...
                                  circshift(ang(ind,marker,'spine_upper',1), 1)),...
             ang,ind);

    sp(end+1) = subplot2(2,8,1,1:5);hold('on')
    plot([1:size(pdb{s},1)]./xyz{s}.sampleRate,pdb{s},'LineWidth',1);
    plot([1:size(pdh{s},1)]./xyz{s}.sampleRate,pdh{s},'LineWidth',1);
    sp(end+1) = subplot2(2,8,2,1:5);hold('on')
    plot([1:size(adb{s},1)]./xyz{s}.sampleRate,adb{s},'LineWidth',1);
    plot([1:size(adh{s},1)]./xyz{s}.sampleRate,adh{s},'LineWidth',1);
    
    p = xcorr(cat(1,pdb{:}),cat(1,pdh{:}),240,'coeff');
    a = xcorr(cat(1,adb{:}),cat(1,adh{:}),240,'coeff');

    subplot2(2,8,1,[7:8]); hold('on'); grid('on');
    plot([-240:240]./xyz{1}.sampleRate,a,'LineWidth',1)
    xlim([-2,2])

    subplot2(2,8,2,[7:8]); hold('on'); grid('on');
    plot([-240:240]./xyz{1}.sampleRate,p,'LineWidth',1)
    xlim([-2,2]);

end
linkaxes(sp,'x')
xlim(sp(1),tind);
xlim(sp(2),tind);
exampleTimePeriodStr = strjoin(num2str(tind),'-'
FigName = ['ccg_head_body_angular_velocity,'_',exampleTimePeriodStr];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));




% 
phaseFrequencyRangeLow = [1.2,6];
phaseFrequencyRangeHigh = [6,12];
% LOAD features
features = cf(@(t) fet_bref(t), Trials);
featuresPhaseLow  = cf(@(f,frq) f.phase(frq), features,repmat({phaseFrequencyRangeLow},[1,numel(Trials)]));
featuresPhaseHigh = cf(@(f,frq) f.phase(frq), features,repmat({phaseFrequencyRangeHigh},[1,numel(Trials)]));
% FILTER features
ffet = cf(@(f) f.copy(), features);
cf(@(f) f.filter('ButFilter',3,[1.2,6],'bandpass'),ffet);










%% Fig.4 SVD of embedded marker trajectories relative to body at head raises%%%%%%%%%%%%%%%%%%%
%
%


param = struct('sessionList',            'hand_labeled',...
               'referenceTrial',         'jg05-20120317.cof.all',...
               'svdState',               'walk+turn+pause',...
               'preState',               {{'walk,turn,pause'}},...
               'postState',              {{'walk,turn,pause'}},...
               'eigenVectorFeaturesMask',{{[2:15,17:2:25,26:30],[1:16,18:24,26:30]}},...
               'eigenVectorTemporalMask',[1:15,46:64],...
               'eigenVectorIndices',     [1,2],...
               'sampleRate',             119.881035,...
               'embeddingWindow',        64 ...
);




% DEF figure variables -----------------------------------------------------------------

% figure save paths
OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'Shared/Behavior Paper/Figures/Figure_4/parts';

states = {'walk','rear','turn','pause','groom','sit'};
nsts = numel(states);
sclr = 'brgcmy';

% $$$ s = 1;                           % jg05-20120317.cof.all
% $$$ exampleTimePeriod = [2170,2198]; % seconds
% $$$ %exampleTimePeriod = [2183,2191]; % seconds
% $$$ exampleTimePeriodStr = num2str(exampleTimePeriod);
% $$$ exampleTimePeriodStr(exampleTimePeriodStr==' ') = '_';

% END figure variables -----------------------------------------------------------------






% START data processing ------------------------------------------------------------------

sessionList     = get_session_list(param.sessionList);                    
numSessions     = numel(sessionList);                    
sampleRate = repmat({param.sampleRate},1,numSessions);
embeddingWindow = repmat({64},1,numSessions);

                    

% LOAD Trial objects
% $$$ Trials = af(@(Trial) MTATrial.validate(Trial)  , get_session_list('nn_labeled'));
% $$$ StcNN  = cf(@(Trial)  Trial.load('stc'), Trials);
% $$$ StcNNC = cf(@(Trial)  Trial.load('stc',[Trial.stc.mode,'_svdc']), Trials);
% LOAD State Collections


% LOAD Trial objects
Trials = af(@(Trial) MTATrial.validate(Trial)           ,sessionList);

% LOAD State Collections
Stc    = cf(@(Trial) Trial.load('stc')                  ,Trials);
%StcNN  = cf(@(Trial) Trial.load('stc','NN0317')         ,Trials);
StcHL  = cf(@(Trial)  Trial.load('stc')         , Trials);
StcHLC = cf(@(Trial)  Trial.load('stc',[Trial.stc.mode,'_SVDTRAJADJ']), Trials);


% LOAD Position data
xyz    = cf(@(Trial) preproc_xyz(Trial)                 ,Trials);
         cf(@(x,s)   x.resample(s)                      ,xyz,sampleRate);
fxyz   = cf(@(x)     x.copy()                           ,xyz);
         cf(@(f)     f.filter('ButFilter',5,[2.4],'low'),fxyz);

% LOAD and MAP Features to reference session
fet    = cf(@(Trial) fet_bref(Trial)                    ,Trials);
         cf(@(f,t,r) f.map_to_reference_session(t,r)    ,fet,Trials,...
            repmat({param.referenceTrial},1,numSessions));
for s = 1:numSessions, fet{s}.data(~nniz(xyz{s}),:,:) = 0;end

% NORMALIZE feature matrix along the columns 
zfrCat = cf(@(f) get(f,'data'),fet);
zfrCat = cat(1,zfrCat{:});
zfrMean = nanmean(zfrCat(nniz(zfrCat),:,:));
zfrStd = nanstd(zfrCat(nniz(zfrCat),:,:));
cf(@(w,m,s) set(w,'data',nunity(w,[],m,s)),...
            fet,...
            repmat({zfrMean},1,numSessions),...
            repmat({zfrStd},1,numSessions));
clear('zfrCat','zfrMean','zfrStd')

% FILTERED feature matrix
ffet   = cf(@(f)     f.copy                            ,fet);
         cf(@(f)     f.filter('ButFilter',5,[1.5],'low'),ffet);

headRaiseEvents = cf(@(f) ThreshCross(f(:,15),0,50),ffet);

ind = cf(@(s) [s{'w+n+p'}],StcHL);
headRaiseEvents = cf(@(h,i) SelectPeriods(h(:,1),i.data,'d',1),headRaiseEvents,ind);


% EMBED feature
wfs  = cf(@(w,e) w.segs(1:size(w,1),e),fet,embeddingWindow);
wfs =  cf(@(w,e) circshift(w,e/2,2),wfs,embeddingWindow);
wfs =  cf(@(w,x) MTADxyz('data',reshape(permute(w,[2,1,3]),size(w,2),[]),...
              'sampleRate',x.sampleRate),wfs,xyz);
for i = 1:numel(wfs), wfs{i}.data(isnan(wfs{i}.data(:)))=0; end

% @svd
% DECOMPOSE fet with svd for walk and turn periods within all sessions
wfw     = cf(@(w,h)    w(h,:), wfs, headRaiseEvents);
%wfw = cf(@(w,s) w([s{'w'}]+[0.5,-0.5],:), wfs, Stc);
[~,Sww,Vww] = svd(cat(1,wfw{:}),0);

% @afetW
% COMPUTE eigenvector loadings for each session's eigen vectors
afetW   = cf(@(w,v) MTADxyz('data',multiprod(w.data,v(:,1:20)),...
                          'sampleRate',w.sampleRate),...
             wfs,repmat({Vww},1,numSessions));
cf(@(f,x) set(f,'sync'  ,x.sync.copy)    ,afetW, xyz); 
cf(@(f,x) set(f,'origin',x.origin   )    ,afetW, xyz);


% COMPUTE eigenvector loadings for each session's eigen vectors
% contains mask to select feature subset important to turning
sfet = cf(@(x) x.copy('empty'), xyz);
for i= param.eigenVectorIndices
    eigenVector = reshape(Vww(:,i),[],size(fet{1},2));
    eigenVector(:,param.eigenVectorFeaturesMask{i}) = 0;
    eigenVector(param.eigenVectorTemporalMask,:) = 0;
    eigenVector = reshape(eigenVector,[],1);
    cf(@(r,w,v) set(r,'data',[get(r,'data'),multiprod(w.data,v)]),...
       sfet,wfs,repmat({eigenVector},1,numSessions));
end
cf(@(f,x) set(f,'sync'  ,x.sync.copy)   ,sfet, xyz); 
cf(@(f,x) set(f,'origin',x.origin   )   ,sfet, xyz);


% @ts
% TIME vectors
wts = cf(@(e,s)  [1:e]./s,                        embeddingWindow,sampleRate);
ts =  cf(@(x)    [1:size(x,1)]./x.sampleRate,     xyz);

% @ang
% COMPUTE  intermarker angles 
ang = cf(@(t,x) create(MTADang,t,x), Trials, fxyz);
for s = 1:numSessions, ang{s}.data(~nniz(xyz{s}),:,:,:) = 0;end

% FIG4EIGVECT ------------------------------------------------------------------------------
% project: MjgEd2016
% parent: figure 4
% subplots: 1->N eigen vectors for all sessions during walking periods
% location: MjgEd2016_figure3_svd_walk_alt.m

hfig = figure;
hfig.Units = 'centimeters';
hfig.Position(3:4) = [30,16];
hfig.PaperPositionMode = 'auto';
for i = 1:20,
    pc = -reshape(Vww(:,i),[],size(fet{1},2));
    %pc = reshape(LR(:,i),[],size(fet{1},2));
    pc = pc(:,[1:2:9,2:2:10,11:15,16:2:24,17:2:25,26:30]);
    subplot(2,10,i);imagesc(wts{1},1:size(fet{1},2),pc'),
    caxis([-0.08,0.08]);
    %caxis([-.5,.5]);    
    axis xy
end


FigName = ['SVD_eigVect_headRaises'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

% END FIG3EIGVECT ---------------------------------------------------------------------------


%figure,plot(fet{1}(:,[17,25,30]));



% FIG4PITCHANGVEL ---------------------------------------------------------------------------
OwnDir = '/storage/gravio/nextcloud/';
FigDir = 'Shared/Behavior Paper/Figures/Figure_4/parts';



sessionList = 'hand_labeled';
Trials = af(@(Trial)  MTATrial.validate(Trial),  get_session_list(sessionList));
stc = cf(@(t)  t.load('stc'),  Trials);
%features = cf(@(Trial)  fet_HB_pitchvel(Trial),  Trials);
features = cf(@(Trial)  fet_HB_angvel(Trial),  Trials);

s = 1;                           % jg05-20120317.cof.all
%exampleTimePeriod = [2170,2198]; % seconds
exampleTimePeriod = [2183,2191]; % seconds
exampleTimePeriodStr = num2str(exampleTimePeriod);
exampleTimePeriodStr(exampleTimePeriodStr==' ') = '_';



hfig = figure(gen_figure_id);
hfig.Units = 'centimeters';
hfig.Position(3:4) = [15,10];
% subplot 1 - Feature Matrix
subplot2(4,1,1:3,1); 
plot([1:size(features{s},1)]./features{s}.sampleRate,features{s}.data.*features{s}.sampleRate);
ylim([-15,15]);
legend({'body','head'});

% subplot 2 - State Labels
subplot2(4,1,4,1);
plotSTC(stc{s},1);
linkaxes(findobj(hfig,'Type','axes'),'x');

% preprint formating
ForAllSubplots(['xlim([',num2str(exampleTimePeriod),'])'])

% save figure
TrialName = [Trials{s}.name,'.',Trials{s}.maze.name,'.',Trials{s}.trialName];
FigName = ['F4_',features{s}.label,'_',TrialName,'_',exampleTimePeriodStr];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
% END FIG4PITCHANGVEL ------------------------------------------------------------------------