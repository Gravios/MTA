function Stc = adjust_state_boundaries_svd(Stc,Trial,param,varargin)
%function  adjust_state_boundaries_svd(Stc,Trial,param,varargin)
% param (struct):
%    sessionList             (string)                   'hand_labeled'
%    referenceTrial          (string)                   'jg05-20120317.cof.all'
%    featureSet              (string)                   'fet_bref'
%    sampleMode              (string)                   'centered','trimmed'
%    svdState                (string)                   'rear' 
%    antecedentState         (string)                   'rear'
%    subsequentState         (string)                   'gper-rear'
%    immutableStates         (cellArray{string})        {{}} 
%    sampleRate              (numeric)                  119.881035
%    eigenVectorFeaturesMask (cellArray{numericArray})  {{[6:10,16:25],[1:10,16:25]}}
%    eigenVectorTemporalMask (numericArray)             [1:15,46:64]
%    eigenVectorIndices      (numericArray)             [1,2]
%    sortTurnsIndex          (logicalArray)             [false,false]
%    embeddingWindow         (numeric)                  64
%    regressionWindow        (numericArray)             100:181,
%    regressionThreshold     (numeric)                  100
%    residualSearchWindow',  (numeric)                  0.25 
%    medianCorrectionOffset  (numeric)                  0.1333 



% TESTARGS ----------------------------------------------------------------------------------------
% $$$ param = struct('sessionList',            'hand_labeled',...
% $$$                'referenceTrial',         'jg05-20120317.cof.all',...
% $$$                'featureSet',             'fet_bref',...
% $$$                'sampleMode',             'centered',...
% $$$                'svdState',               'rear',...
% $$$                'antecedentState',        'gper-rear',...
% $$$                'subsequentState',        'rear',...
% $$$                'immutableStates',        {{}} 
% $$$                'eigenVectorFeaturesMask',{{[6:10,16:25],[1:10,16:25]}},...
% $$$                'eigenVectorTemporalMask',[1:15,46:64],...
% $$$                'eigenVectorIndices',     [1,2],...
% $$$                'sampleRate',             119.881035,...
% $$$                'embeddingWindow',        64, ...                    
% $$$                'regressionWindow',       100:181, ...
% $$$                'regressionThreshold',    5e4,...
% $$$                'residualSearchWindow',   0.25 
% $$$                'medianCorrectionOffset', 0.14);
% $$$ 
% $$$ varargin = {};
% $$$ Trial = [];
% $$$ Stc = [];
% $$$ states = {'walk','rear','turn','pause','groom','sit'};

%--------------------------------------------------------------------------------------------------



% DEFARGS -----------------------------------------------------------------------------------------
defargs = struct(                                                                               ...
    'states',                {{'walk','rear','turn','pause','groom','sit'}},                    ...
    'overwriteNORM',         false,                                                             ...
    'overwriteSVD',          false,                                                             ...
    'overwriteREG',          false                                                              ...
);
[states,overwriteNORM,overwriteSVD,overwriteREG] = DefaultArgs(varargin,defargs,'--struct');


normalizationParameters = struct('sessionList',            param.sessionList,                   ...
                                 'referenceTrial',         param.referenceTrial,                ...
                                 'featureSet',             param.featureSet                     ...
);

svdParameters =           struct('sessionList',            param.sessionList,                   ...
                                 'referenceTrial',         param.referenceTrial,                ...
                                 'featureSet',             param.featureSet,                    ...
                                 'svdState',               param.svdState,                      ...
                                 'sampleMode',             param.sampleMode,                    ...
                                 'sampleRate',             param.sampleRate,                    ...
                                 'embeddingWindow',        param.embeddingWindow                ...
);

regressionParameters =    struct('sessionList',            param.sessionList,                   ...
                                 'referenceTrial',         param.referenceTrial,                ...
                                 'featureSet',             param.featureSet,                    ...
                                 'svdState',               param.svdState,                      ...
                                 'sampleMode',             param.sampleMode,                    ...
                                 'sampleRate',             param.sampleRate,                    ...
                                 'embeddingWindow',        param.embeddingWindow,               ...
                                 'regressionWindow',       param.regressionWindow,              ...
                                 'regressionThreshold',    param.regressionThreshold            ...
);
%--------------------------------------------------------------------------------------------------



% TAG creation ------------------------------------------------------------------------------------
% ID Vars - create hash tag
normalizationTag = DataHash(normalizationParameters);
svdTag           = DataHash(svdParameters);
regressionTag    = DataHash(regressionParameters);
%--------------------------------------------------------------------------------------------------



% PREPROCESSING ----------------------------------------------------------------------------------


% WRAP basic vars into cells
if isempty(Trial)
    sessionList = get_session_list(param.sessionList);
elseif ischar(Trial)
    sessionList = get_session_list(Trial);
elseif iscell(Trial),
    sessionList = Trial; 
else
    sessionList = {Trial}; 
end
numSessions     = numel(sessionList);                    
sampleRate = repmat({param.sampleRate},1,numSessions);
embeddingWindow = repmat({param.embeddingWindow},1,numSessions);
trimWindow = repmat({[ param.embeddingWindow./4./param.sampleRate,...
                      -param.embeddingWindow./4./param.sampleRate]},1,numSessions);
residualSearchWindow = round(param.residualSearchWindow.*param.sampleRate);

% LOAD Trial objects
if isstruct(sessionList),
    Trials = af(@(Trial) MTATrial.validate(Trial)  , sessionList);
else
    Trials = cf(@(Trial) MTATrial.validate(Trial)  , sessionList);
end    

% LOAD State Collections
if isempty(Stc),
    Stc    = cf(@(Trial)  Trial.load('stc')         , Trials);
elseif ~iscell(Stc) && isa(Stc,'MTAStateCollection'),
    Stc = {Stc};
elseif ischar(Stc),
    Stc = {Trial.load('stc',Stc)};
end

% LOAD Position data
xyz    = cf(@(Trial) preproc_xyz(Trial)         , Trials);
cf(@(x,s) x.resample(s),xyz,sampleRate);
fxyz = cf(@(x) x.copy(), xyz);
cf(@(f) f.filter('ButFilter',5,[2.4],'low'),fxyz);

% LOAD and MAP Features to reference session
fet    = cf(@(Trial) fet_bref(Trial), Trials);
cf(@(f,t) f.map_to_reference_session(t,param.referenceTrial),fet,Trials);
for s = 1:numSessions, fet{s}.data(~nniz(xyz{s}),:,:) = 0;end



% NORMALIZE feature matrix
normalizationParameterFile = ...
    fullfile(MTASession([]).path.project,'analysis',...
             ['normalization_parameters-',normalizationTag,'.mat']);
featureCat = [];
if ~exist(normalizationParameterFile,'file')||overwriteNORM,
    featureCat = cf(@(f) get(f,'data'),fet);
    featureCat = cat(1,featureCat{:});
    featureMean = nanmean(featureCat(nniz(featureCat),:,:));
    featureStd = nanstd(featureCat(nniz(featureCat),:,:));
    save(normalizationParameterFile,'normalizationParameters','featureMean','featureStd');
else
    load(normalizationParameterFile);
end
cf(@(w,m,s) set(w,'data',nunity(w,[],m,s)),...
            fet,...
            repmat({featureMean},1,numSessions),...
            repmat({featureStd},1,numSessions));
clear('featureCat','featureMean','featureStd')



% EMBBED feature
wfs  = cf(@(w,e) w.segs(1:size(w,1),e),fet,embeddingWindow);
wfs =  cf(@(w,e) circshift(w,e/2,2),wfs,embeddingWindow);
wfs =  cf(@(w,x) MTADxyz('data',reshape(permute(w,[2,1,3]),size(w,2),[]),...
              'sampleRate',x.sampleRate),wfs,xyz);
for i = 1:numel(wfs), wfs{i}.data(isnan(wfs{i}.data(:)))=0; end



% EMBBEDDED singular value decomposition
svdParameterFile = ...
    fullfile(MTASession([]).path.project,'analysis',...
             ['svd_parameters-',svdTag,'.mat']);
if ~exist(svdParameterFile)||overwriteSVD,
% DECOMPOSE fet with svd for a behavioral state within all Trials
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
else
% LOAD SVD
    load(svdParameterFile);
end



% COMPUTE eigenvector loadings for each session's eigen vectors
% contains mask to select feature subset important to turning
sfet = cf(@(x) x.copy('empty'), xyz);
for i = 1:numel(param.eigenVectorIndices),
    eigenVector = reshape(Vww(:,i),[],size(fet{1},2));
    eigenVector(:,param.eigenVectorFeaturesMask{i}) = 0;
    eigenVector(param.eigenVectorTemporalMask,:) = 0;
    eigenVector = reshape(eigenVector,[],1);
    cf(@(r,w,v) set(r,'data',[get(r,'data'),multiprod(w.data,v)]),...
       sfet,wfs,repmat({eigenVector},1,numSessions));
end
cf(@(f,x) set(f,'sync',x.sync.copy), sfet, xyz); 
cf(@(f,x) set(f,'origin',x.origin), sfet, xyz);



% COMPUTE  intermarker angles 
ang = cf(@(t,x) create(MTADang,t,x), Trials, fxyz);
for s = 1:numSessions, ang{s}.data(~nniz(xyz{s}),:,:,:) = 0;end


% END PREPROCESSING -----------------------------------------------------------------------




% MAIN  -----------------------------------------------------------------------------------

shift = [-round(param.sampleRate),round(2*param.sampleRate)];


% OPTIMIZATION of state transitions
sts = cf(@(s,w,t) [s.get_state_transitions(t,...
                            {param.antecedentState,param.subsequentState},[],w)],...
         Stc, sfet, Trials);



% TRIM timepoints at ends if necessary
for s = 1:numSessions, 
    sts{s}(sum([sts{s}+repmat(shift,size(sts{s},1),1)<=0, ...
                sts{s}+repmat(shift,size(sts{s},1),1)>size(sfet{s},1)],2)>0,:)=[];
end        



% SEPARATE left and right trajectories for turns
if any(~cell2mat(cf(@isempty,regexp({param.antecedentState,param.subsequentState},'turn')))),
    if ~isempty(regexp(param.subsequentState,'turn'))
        shift = [0,round(0.5*param.sampleRate)];
    else
        shift = [round(0.5*param.sampleRate),-round(0.25*param.sampleRate)];
    end
    
    stsSind = cf(@(s,a) sign(circ_dist(...
        a(round(mean(s,2)+shift(1)),'spine_lower','spine_upper',1),...
        a(round(mean(s,2)+shift(2)),'spine_lower','spine_upper',1)))==1,...
                 sts,ang);
else
% OR DON'T
    stsSind = cf(@(s,a) true([size(s,1),1]),sts,ang);
end



% LOAD regression parameters
regressionParameterFile = ...
    fullfile(MTASession([]).path.project,'analysis',...
             ['regression_parameters-',regressionTag,'.mat']);
if ~exist(regressionParameterFile,'file')||overwriteREG,
% COLLECT segments
    for sid = 1:numel(param.sortTurnsIndex),
        tfet = cf(@(s) s.copy(), sfet);
        cf(@(t,sid) set(t,'data',t(:,sid)),tfet,repmat({sid},1,numSessions));
        if param.sortTurnsIndex(sid)
            wfw = cf(@(w,s,t,i)  w.segs(round(mean(s(i==1,:),2))-round(1.*w.sampleRate),round(3.*w.sampleRate)),...
                     tfet, sts, Trials, stsSind);
            wfn = cf(@(w,s,t,i) -w.segs(round(mean(s(i==0,:),2))-round(1.*w.sampleRate),round(3.*w.sampleRate)),...
                     tfet, sts, Trials, stsSind);
        else
            wfw = cf(@(w,s,t,i)  w.segs(round(mean(s(i==1,:),2))-round(1.*w.sampleRate),round(3.*w.sampleRate)),...
                     tfet, sts, Trials, stsSind);
            wfn = cf(@(w,s,t,i)  w.segs(round(mean(s(i==0,:),2))-round(1.*w.sampleRate),round(3.*w.sampleRate)),...
                     tfet, sts, Trials, stsSind);
        end
% COCATENATE segment data
        fetSegs(:,:,sid) = cat(2,wfw{:},wfn{:});
    end
% COMPUTE mean trajectory
    mfs = sq(nanmean(fetSegs(param.regressionWindow,:,:),2))';
    save(regressionParameterFile,'regressionParameters','mfs');
else
    load(regressionParameterFile);
end



rof = cf(@(a) a.copy(),sfet);
cf(@(f,m) set(f,'data',circshift(f.segs(1:size(f,1),size(m,2)),round(size(m,2)/2),2)),...
   rof,repmat({mfs},size(rof)));
cf(@(f,m) set(f,'data',repmat(permute(f.data,[4,1,2,3]),3,1,1)-...
                repmat(linspace(-.5,.5,3)',[1,size(m,2),size(f,2),size(f,3)])),...
   rof,repmat({mfs},size(rof)));

% COMPUTE regressed feature
csw = cf(@(f) f.copy(),rof);
cf(@(c,f,m) set(c,'data',sq(sum((f.data-repmat(permute(m,[3,2,4,1]),size(f,1),1,size(f,3))).^2,2))),...
   csw,rof,repmat({mfs},size(rof)));

% DETECT transitions based on regression score
nsmins = cf(@(c,t) LocalMinima(sum(c(2,:,:),3),0,t),...
            csw,repmat({param.regressionThreshold},1,numSessions));

% SELECT putative transitions around state label transitions 
%% THIS SECTION NEEDS REPAIRS
[nsmins,nsinds] = cf(@(m,s,w)                                                                     ...
                     SelectPeriods(m,JoinRanges(bsxfun(@plus,[round(nanmean(s,2))],[-w,w]),[1,2]),'d',1),...
                     nsmins,sts,repmat({residualSearchWindow},1,numSessions));
stsstc = cf(@(m,s,p,w) SelectPeriods([s{p}(:,1)],bsxfun(@plus,mean(m,2),[-5,5]),'d',1),...
                     sts,Stc,repmat({param.subsequentState},1,numSessions),...
                     repmat({residualSearchWindow},1,numSessions));

%nsmins = cf(@(n) unique(n),nsmins);
[ssmins,ssinds] = cf(@(m,s,w)                                                                ...
                     SelectPeriods(mean(s,2),JoinRanges(bsxfun(@plus,m,[-w,w]),[1,2]),'d',1),...
                     nsmins,stsstc,repmat({residualSearchWindow},1,numSessions));
[nsmins,nsinds] = cf(@(m,s,w)                                                                     ...
                     SelectPeriods(m,JoinRanges(bsxfun(@plus,[round(nanmean(s,2))],[-w,w]),[1,2]),'d',1),...
                     nsmins,ssmins,repmat({residualSearchWindow},1,numSessions));





% $$$ 
% $$$ 
% $$$ [ssmins,ssinds] = cf(@(m,s,w) SelectPeriods(s,bsxfun(@plus,m,[-w,w]),'d',1),...
% $$$                      nsmins,stsstc,repmat({residualSearchWindow},1,numSessions));
% $$$ ssmins = cf(@(s) unique(s),ssmins);
% $$$ 
% $$$ [nsmins,nsinds] = cf(@(m,s,p,w) SelectPeriods(m,bsxfun(@plus,[s{p}(:,1)],[-w,w]),'d',1),...
% $$$                      nsmins,Stc,repmat({param.subsequentState},1,numSessions),...
% $$$                      repmat({residualSearchWindow},1,numSessions));
% $$$ nsmins = cf(@(n) unique(n),nsmins);
% $$$ [ssmins,ssinds] = cf(@(m,s,p,w) SelectPeriods([s{p}(:,1)],bsxfun(@plus,m,[-w,w]),'d',1),...
% $$$                      nsmins,Stc,repmat({param.subsequentState},1,numSessions));
% $$$                      repmat({residualSearchWindow},1,numSessions));
%% THIS SECTION NEEDS REPAIRS


% REMOVE redundant onsets; keep nearest to Stc label
for s = 1:numSessions,
    i = 1;
    while i < numel(nsmins{s}),
        if abs(nsmins{s}(i+1) - ssmins{s}(i)) <= residualSearchWindow*2
            if abs(nsmins{s}(i)-ssmins{s}(i)) >= abs(nsmins{s}(i+1)-ssmins{s}(i))
                nsmins{s}(i) = [];
            else,
                nsmins{s}(i+1) = [];
            end
        else
            
            i = i+1;            
        end
    end
end

% LAZY bum
stsstc = cf(@(m,s,p,w) SelectPeriods([s{p}(:,1)],bsxfun(@plus,mean(m,2),[-5,5]),'d',1),...
                     sts,Stc,repmat({param.subsequentState},1,numSessions),...
                     repmat({residualSearchWindow},1,numSessions));
[ssmins,ssinds] = cf(@(m,s,w)                                                                ...
                     SelectPeriods(mean(s,2),JoinRanges(bsxfun(@plus,m,[-w,w]),[1,2]),'d',1),...
                     nsmins,stsstc,repmat({residualSearchWindow},1,numSessions));
[nsmins,nsinds] = cf(@(m,s,w)                                                                     ...
                     SelectPeriods(m,JoinRanges(bsxfun(@plus,[round(nanmean(s,2))],[-w,w]),[1,2]),'d',1),...
                     nsmins,ssmins,repmat({residualSearchWindow},1,numSessions));

% REMOVE redundant onsets; keep nearest to Stc label
for s = 1:numSessions,
    i = 1;
    while i < numel(nsmins{s}),
        if abs(nsmins{s}(i+1) - ssmins{s}(i)) <= residualSearchWindow*2
            if abs(nsmins{s}(i)-ssmins{s}(i)) >= abs(nsmins{s}(i+1)-ssmins{s}(i))
                nsmins{s}(i) = [];
            else,
                nsmins{s}(i+1) = [];
            end
        else
            
            i = i+1;            
        end
    end
end



% ADD temporal offset to detected transitions
if param.medianCorrectionOffset
    %medianCorrectionOffset = repmat({median(cat(1,nsmins{:})-cat(1,ssmins{:}))},1,numSessions);
    nsmins = cf(@(n,m) n-m, nsmins, repmat({round(param.medianCorrectionOffset.*param.sampleRate)},1,numSessions));
end


% APPLY adjusted state transitions to the state collection (Stc)
for s = 1:numSessions,
    [smat,labels,keys] = stc2mat(Stc{s},fet{s},states);
    %stateIndex = find(~cellfun(@isempty,regexp(labels,['^',param.subsequentState,'$'])));
    sper = Stc{s}.states{Stc{s}.gsi(param.subsequentState)};
    for i = 1:numel(nsmins{s}),
% FIND index of old period transition 
        oldSubPerIndex = find(ssmins{s}(i)==sper.data(:,1));
        shiftMagnitude = nsmins{s}(i)-ssmins{s}(i);
        stateIndex = smat(sper.data(oldSubPerIndex,1)+1,:);
% CLEAR region for transition adjustment
        if shiftMagnitude == 0,
            continue;
        elseif shiftMagnitude < 0,            
            ind = sper.data(oldSubPerIndex,1)+shiftMagnitude:sper.data(oldSubPerIndex,1)-1;
            backfill = false;
            smat(ind,~ismember(1:size(smat,2),stateIndex)) = 0;
            %[Stc{s}.gsi(param.immutableStates)]
            
        elseif shiftMagnitude > 0,            
            ind = sper.data(oldSubPerIndex,1):sper.data(oldSubPerIndex,1)+shiftMagnitude-1;
            backfill = true;
            smat([ind],logical(stateIndex)) = 0;                        
        end        

        if backfill,
% FIND antecedent states
            backSearchInd = sper.data(oldSubPerIndex,1)-round(0.25.*sampleRate{1}):sper.data(oldSubPerIndex,1);
            backSearchInd(backSearchInd<=0)=[];
            backStateInd = find(any(smat(backSearchInd,:)) & ~stateIndex);
            if isempty(backStateInd), continue; end
            smat(ind,backStateInd) = repmat(backStateInd,numel(ind),1);
        else
            smat(ind,:) = repmat(stateIndex,numel(ind),1);
        end    
    end    

    %StcN{s} = mat2stc(smat,Stc{s}.copy,fet{s},Trials{s},labels,keys);
    Stc{s} = mat2stc(smat,Stc{s},fet{s},Trials{s},labels,keys);
end


% END MAIN ---------------------------------------------------------------------------------