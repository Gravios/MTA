

param = struct('sessionList',            'hand_labeled',...
               'referenceTrial',         'jg05-20120317.cof.all',...
               'svdState',               'walk+turn',...
               'preState',               {{'walk,turn,pause'}},...
               'postState',              {{'walk,turn,pause'}},...
               'eigenVectorFeaturesMask',{{[6:10,16:25],[1:10,16:25]}},...
               'eigenVectorTemporalMask',[1:15,46:end],...
               'eigenVectorIndices',     [1,2],...
               'sampleRate',             119.881035,...
               'embeddingWindow',        64...                    
);



sessionList     = get_session_list(param.sessionList);                    
numSessions     = numel(sessionList);                    
sampleRate = repmat({param.sampleRate},1,numSessions);
embeddingWindow = repmat({64},1,numSessions);

                    
% LOAD Trial objects
Trials = af(@(Trial) MTATrial.validate(Trial)  , sessionList);

% LOAD State Collections
Stc    = cf(@(Trial)  Trial.load('stc')         , Trials);
StcNN  = cf(@(Trial) Trial.load('stc','NN0317'), Trials);

% LOAD Position data
xyz    = cf(@(Trial) preproc_xyz(Trial)         , Trials);
cf(@(x,s) x.resample(s),xyz,sampleRate);
fxyz = cf(@(x) x.copy(), xyz);
cf(@(f) f.filter('ButFilter',5,[2.4],'low'),fxyz);

% LOAD and MAP Features to reference session
fet    = cf(@(Trial) fet_bref(Trial), Trials);
cf(@(f,t) f.map_to_reference_session(t,param.referenceTrial),fet,Trials);
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
