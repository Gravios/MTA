
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
fetSubsetGroomId = cf(@(f,s) ones([size(f,1),1]).*s,fetSubsetGroom,mat2cell(1:6,1,ones([1,numSessions])));
fetSubsetGroomIdAll = cat(1,fetSubsetGroomId{:});
fetSubsetGroomAll = cat(1,fetSubsetGroom{:});

tsneMap_sub = tsne(fetSubsetGroomAll(1:4:end,:),fetSubsetGroomIdAll(1:4:end),2,5,80);




groomPerTs = cf(@(s) s{'m'}, Stc);
cf(@(p) p.cast('TimeSeries',12), groomPerTs);


figure,hold on
for i = 1:numSessions,
    ind = fetSubsetGroomIdAll==i;
    scatter(tsneMap(ind,1),tsneMap(ind,2),20,sclr(i))
end

hfig = figure();
plot(tsneMap(:,1),tsneMap(:,2),'.');
[cpnts]=ClusterPP(hfig);



fetSubsetGroomCumSize = mat2cell(cumsum([0,cellfun(@length, ...
                    fetSubsetGroom(1:end-1))]),1,ones([1,numSessions]));

groomPerTsIds = cf(@(p) p.copy(), groomPerTs);


cf(@(p,c,t) set(p,'data',t(c+1:c+sum(p.data==1))),...
   groomPerTsIds,fetSubsetGroomCumSize,repmat({cpnts},1,numSessions))
    
for i = 1:numSessions,
    groomPerTsIds{i}.data(groomPerTsIds{i}.data==1) =  ...
        cpnts(fetSubsetGroomCumSize{i}+1:fetSubsetGroomCumSize{i}+sum(groomPerTsIds{i}.data==1));
end


groomSubStateTags = cf(@(t) ['groom_ss_',t], ...
                       mat2cell(num2str(unique(cpnts),'%i'),1,ones([1,numel(unique(cpnts))])));




