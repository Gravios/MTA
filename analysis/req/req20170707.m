
param = struct('sessionList',            'hand_labeled',...
               'referenceTrial',         'jg05-20120317.cof.all',...
               'svdState',               'walk+turn',...
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
FigDir = 'Shared/Behavior Paper/Figures/Figure_2/parts';

states = {'walk','rear','turn','pause','groom','sit'};
nsts = numel(states);
sclr = 'brgcmy';


s = 1;                           % jg05-20120317.cof.all
exampleTimePeriod = [2170,2198]; % seconds
%exampleTimePeriod = [2183,2191]; % seconds
exampleTimePeriodStr = num2str(exampleTimePeriod);
exampleTimePeriodStr(exampleTimePeriodStr==' ') = '_';

% END figure variables -----------------------------------------------------------------



% END figure variables -----------------------------------------------------------------






% START data processing ------------------------------------------------------------------

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

dsfet = cf(@(f) f.copy(),fet);
cf(@(f) f.resample(12), dsfet);

fetSubsetGroom = cf(@(f,s) f(s{'m'},:),dsfet,Stc);
fetSubsetGroom = cat(1,fetSubsetGroom{:});

tsneMap = tsne(fetSubsetGroom,[],2,30,80);