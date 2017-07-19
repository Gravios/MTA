
param = struct('sessionList',            'hand_labeled',...
               'referenceTrial',         'jg05-20120317.cof.all',...
               'svdState',               'rear',...
               'preState',               {{'rear'}},...
               'postState',              {{'gper'}},...
               'eigenVectorFeaturesMask',{{[6:10,16:25],[1:10,16:25]}},...
               'eigenVectorTemporalMask',[1:15,46:64],...
               'eigenVectorIndices',     [1,2],...
               'sampleRate',             119.881035,...
               'embeddingWindow',        64 ...                    
);

sessionList     = get_session_list(param.sessionList);                    
numSessions     = numel(sessionList);                    
sampleRate = repmat({param.sampleRate},1,numSessions);
embeddingWindow = repmat({64},1,numSessions);
                    
% LOAD Trial objects
Trials = af(@(Trial) MTATrial.validate(Trial)  , sessionList);

% LOAD State Collections
StcHL  = cf(@(Trial)  Trial.load('stc')         , Trials);
StcHLC = cf(@(Trial)  Trial.load('stc',[Trial.stc.mode,'_SVDTRAJADJ']), Trials);

Trials = af(@(Trial) MTATrial.validate(Trial)  , get_session_list('nn_labeled'));
StcNN  = cf(@(Trial)  Trial.load('stc'), Trials);
StcNNC = cf(@(Trial)  Trial.load('stc',[Trial.stc.mode,'_svdc']), Trials);

% LOAD Position data
xyz    = cf(@(Trial) preproc_xyz(Trial)         , Trials);
cf(@(x,s) x.resample(s),xyz,sampleRate);
fxyz = cf(@(x) x.copy(), xyz);
cf(@(f) f.filter('ButFilter',5,[2.4],'low'),fxyz);

% LOAD and MAP Features to reference session
fet    = cf(@(Trial) fet_bref_rev5(Trial), Trials);
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




% DIAGNOSTIC of mapminmax for neural network normalization. 
% Need to add parameter struct to features which contains the
% domains of all features
[StcRnd,labelingEpochs,trainingFeatures] = cf(@(s,f,sts) ...
    resample_whole_state_bootstrap_noisy_trim(s,f,sts,[100],10000),...
    Stc,fet,repmat({states},1,numel(Trials)));


[x,ps] = cf(@(f) mapminmax(f.data'),trainingFeatures);

fdata = [];
for s = 1:numSessions,
    fdata = cat(1,fdata,trainingFeatures{s}.data);
end
[xa,psa] = mapminmax(fdata');

fet_bref_rev5.mapminmax_args.mat

psa = 

         name: 'mapminmax'
        xrows: 30
         xmax: [30x1 double]
         xmin: [30x1 double]
       xrange: [30x1 double]
        yrows: 30
         ymax: 1
         ymin: -1
       yrange: 2
    no_change: 0
         gain: [30x1 double]
      xoffset: [30x1 double]

mm = bsxfun(@plus,prctile(fdata,[0.1,99.9]),[-2;2]);

psa.name: 'mapminmax'      
psa.xrows = size(fet,2);
psa.xmax = mm(2,:)';
psa.xmin = mm(1,:)';
psa.xrange = diff(mm)';
psa.yrows = size(fet,2);
psa.ymax  = 1;
psa.ymin: =-1;
psa.yrange= 2;
psa.no_change = 0;
psa.gain = 2./psa.xrange;
psa.xoffset = psa.xmin;


figure,
eds = linspace(psa.xmin(15),psa.xmax(15),100);
bar(eds,histc(fet{6}(:,15),eds),'histc');


% CCG of state collections 


sto = cf(@(s) s{'r'},StcHL);
stn = cf(@(s) s{'r'},StcNN);
stoc = cf(@(s) s{'r'},StcHLC);
stnc = cf(@(s) s{'r'},StcNNC);


% PLOT CCG between bhv transition and residules minimas - used to determin medianCorrectionOffset
%medianCorrectionOffset = median(cat(1,nsmins{:})-cat(1,ssmins{:}))./sampleRate{1};
figure();
o = 1;
transType = {'onset','offset'};

subplot(221);
[sccg,txx,pxx] = cf(@(s,n,sr) CCG([s(:,o);n(:,o)],[ones([size(s,1),1]);2*ones([size(n,1),1])],...
                                  2,40,sr,[1,2],'count'),...
                    sto,stn,sampleRate);
accg = sum(cat(4,sccg{:}),4);

bar(txx{1},accg(:,1,2));
ylabel('count')


subplot(222);
[sccg,txx,pxx] = cf(@(s,n,sr) CCG([s(:,o);n(:,o)],[ones([size(s,1),1]);2*ones([size(n,1),1])],...
                                  2,40,sr,[1,2],'count'),...
                    sto,stnc,sampleRate);
accg = sum(cat(4,sccg{:}),4);
bar(txx{1},accg(:,1,2));
ylabel('count')


subplot(223);
[sccg,txx,pxx] = cf(@(s,n,sr) CCG([s(:,o);n(:,o)],[ones([size(s,1),1]);2*ones([size(n,1),1])],...
                                  2,40,sr,[1,2],'count'),...
                    stoc,stn,sampleRate);
accg = sum(cat(4,sccg{:}),4);
bar(txx{1},accg(:,1,2));
ylabel('count')


subplot(224);
[sccg,txx,pxx] = cf(@(s,n,sr) CCG([s(:,o);n(:,o)],[ones([size(s,1),1]);2*ones([size(n,1),1])],...
                                  2,40,sr,[1,2],'count'),...
                    stoc,stnc,sampleRate);
accg = sum(cat(4,sccg{:}),4);
bar(txx{1},accg(:,1,2));
ylabel('count')


