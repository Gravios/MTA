function StcCor = optimize_stc_transition(Stc,varargin)
% function optimize_stc_transition(varargin)
%
% Default Arguments 
%   varargin:
%     sessionList              (String,MTASession)           ''
%     stcMode                  (String,MTAStateCollection)   'msnnN0+hand_labeled'
%     featureSet               (String)                      'fet_bref_rev7'
%     states                   (CellStr)                     {'walk','rear','turn','pause','groom','sit'}
%     keys                     (CellStr)                     {'w',   'r',   'n',   'p',    'm',    's'  }
%     model                    (String)                      ''
%     sampleRate               (Numeric)                     10
%     nNeurons                 (Numeric)                     25
%     nIter                    (Numeric)                     10
%     randomizationMethod      (String)                      'WSBNT'
%     map2reference            (Logical)                     true
%     normalize                (Logical)                     true
%     referenceTrial           (String)                      'jg05-20120317.cof.all'
%     trainingSessionList      (String)                      'hand_labeled'
%     normalizationSessionList (String)                      'hand_labeled'
%     dropIndex                (NumericArray)                []
%     prctTrain                (Numeric)                     90
%     postprocessingTag        (String)                      '_ppsvd'
%

% global MTA_VERBOSITY_LEVEL

verbosePrefix = 'MTA:classifiers:optimize_stc_transition: ';

% DEFARGS ----------------------------------------------------------------------------------------
defArgs = struct('sessionList',                 '',                                            ...
                 'featureSet',                  'fet_mis',                                     ...
                 'states',                      {{'walk','rear','turn','pause','groom','sit'}},...
                 'keys',                        {{'w','r','n','p','m','s'}},                   ...
                 'model',                       [],                                            ...
                 'sampleRate',                  10,                                            ...
                 'nNeurons',                    25,                                            ...
                 'nIter',                       100                                            ...
                 'randomizationMethod',         'WSBNT',                                       ...
                 'map2reference',               true,                                          ...
                 'normalize',                   true,                                          ...
                 'referenceTrial',              'jg05-20120317.cof.all',                       ...
                 'trainingSessionList',         'hand_labeled',                                ...
                 'normalizationSessionList',    'hand_labeled',                                ...
                 'dropIndex',                   2,                                             ...
                 'prctTrain',                   90,                                            ...
                 'postprocessingTag',           '_ppsvd',                                      ...
                 'modelType',                   'multiSesPatNet'                               ...
                 );
[sessionList,featureSet,states,keys,model,sampleRate,nNeurons,nIter,...
 randomizationMethod,map2reference,normalize,referenceTrial,trainingSessionList,...
 normalizationSessionList,dropIndex,prctTrain,postprocessingTag,modelType,stcMode] = ...
    DefaultArgs(varargin,defArgs,'--struct');
%-----------------------------------------------------------------------------

% varargin = {};




% MAIN ------------------------------------------------------------------------------------------

% LOAD Trials which have the target state collections
Trials = af(@(Trial)  MTATrial.validate(Trial), get_session_list(sessionList));
numTrials = numel(Trials);


% LOAD Stc for opitimization
if isempty(Stc),
    if isempty(model),
        % CREATE model name            
        dropIndexTag = '';
        if ~isempty(dropIndex), 
            dropIndexTag = num2str(dropIndex);
            dropIndexTag(dropIndexTag==[' '])='_';
        end
        model = ['MTAC_BATCH+' trainingSessionList '+'            ...
                 dropIndexTag '+'                                 ...
                 featureSet                                       ...
                 '+SR'   num2str(sampleRate)                      ...
                 'NN'    num2str(nNeurons)                        ...
                 'NI'    num2str(nIter)                           ...
                 'M'     num2str(map2reference)                   ...
                 'MREF+' referenceTrial                           ...                          
                 '+N'    num2str(normalize)                       ...             
                 'NREF+' normalizationSessionList                 ...
                 '+RND'  randomizationMethod                      ...
                 '+PRCT' num2str(prctTrain)                       ...
                 '+STS+' strjoin(keys,'')                         ...
                 '+'     modelType];
    end
    % LOAD state collection based on the default bhv_nn_multi_session_patternnet model name
    StcCor = cf(@(t,m) t.load('stc',m),    Trials, repmat({model},[1,numTrials]));        
elseif ischar(Stc)
    StcCor = cf(@(t,m) t.load('stc',m),    Trials, repmat({Stc},[1,numTrials]));        
elseif iscell(Stc),    
    StcCor = cf(@(s)   s.copy(),         Stc);
elseif isa(Stc,'MTAStateCollection'),
    StcCor = { Stc.copy() };
else
    error('MTA:classifiers:optimize_stc_transition:UnknownStcSignature');
end


cf(@(stc,states) set(stc,'states',stc(states{:})),...
   StcCor,repmat({states},[1,numel(Trials)]));


%<--# REMOVE IN FUTURE VERSION
for s = 1:numTrials,
    if isempty(Trials{s}.fet),
        Trials{s}.fet = MTADfet(Trials{s}.spath,...
                                [],...
                                [],...
                                [],...
                                Trials{s}.sync.copy,...
                                Trials{s}.sync.data(1),...
                                []);                  
    end
end
%# REMOVE IN FUTURE VERSION -->

% LOAD features
features = cf(@(t)     fet_bref(t), Trials);
           cf(@(f,t,r) f.map_to_reference_session(t,r), features,Trials,repmat({referenceTrial},[1,numTrials]));

xyz = cf(@(t) t.load('xyz'),                      Trials);
      cf(@(x) x.filter('ButFilter',5,1.5,'low'),  xyz);

% $$$ vxy = cf(@(x) xyz.vel({'spine_lower','head_back'},[1,2]), xyz);
ang  = cf(@(t,x) create(MTADang,t,x), Trials,xyz);
dang = cf(@(x)   x.copy('empty'),xyz);
cf(@(d,ang) set(d,'data',circ_dist(circshift(ang(:,'spine_lower','spine_upper',1),-1),...
                      circshift(ang(:,'spine_lower','spine_upper',1),1))*ang.sampleRate),...
   dang,ang);

% $$$             StcNN = StcCor.copy();
% $$$             StcNN.updateMode(['nn_',refSession.filebase(1:4)]);
% $$$             StcNN.save(1);



% FILTER Rears 
% ADJUST Rears'  onsets from all periods
% ADJUST Rears' offsets from all periods

% FILTER Turns' duration lt 0.1 sec -> pause
% FILTER Turns  
% ADJUST Turns' onsets from Pause

% FILTER Walks 
% ADJUST Walks' onsets from Pause




% REAR remove periods which have low mean heights 



for s = 1:numTrials,
    disp([verbosePrefix,'Processing: ',Trials{s}.filebase]);



    [stcMatrix] = stc2mat(StcCor{s},xyz{s},states);
    stateVectors = bsxfun(@times,eye(numel(states)),1:numel(states));
    
    disp([verbosePrefix,'Assign all periods with low duration to neighboring states']);    
    [stcMatrix] = reassign_low_duration_state_to_neighboring_states(stcMatrix,StcCor{s}{'n'},0.25*xyz{s}.sampleRate);
    [stcMatrix] = reassign_low_duration_state_to_neighboring_states(stcMatrix,StcCor{s}{'w'},0.25*xyz{s}.sampleRate);
    
    disp([verbosePrefix,'Assign all periods with high mean heights to rear']);    
    for key = 'wnpms',
        keyVec  = stateVectors(find(key==cell2mat(keys)),:);
        rearVec = stateVectors(find('r'==cell2mat(keys)),:);
        rthresh = 135;
        try,
            for rp = StcCor{s}{key}.data',
                rhh = max(features{s}(rp',14));
                if rhh > rthresh,
% REASSIGN state to rear                                    
                    stcMatrix(rp(1):rp(2),:) = repmat(rearVec,[diff(rp)+1,1]);
                end
            end
        catch err
            disp(err);
        end
    end


    disp([verbosePrefix,'REAR remove periods which have too low mean heights']);
    % REAR remove periods which have low mean heights 
    key = 'r';
    keyVec  = stateVectors(find(key==cell2mat(keys)),:);
    pauseVec = stateVectors(find('p'==cell2mat(keys)),:);    
    rthresh = 125;
    try
        for rp = StcCor{s}{key}.data',
            rhh = max(features{s}(rp',14));
            if rhh < rthresh,
% REASSIGN rear to pause                
                stcMatrix(rp(1):rp(2),:) = repmat(pauseVec,[diff(rp)+1,1]);
            end
        end
    catch err
        disp(err);
    end


    % TURN angular displacement
    %try, StcCor{s} = reassign_state_by_duration(StcCor{s},'n','p',0.1,tds,tps,@lt); end
    disp([verbosePrefix 'Sort out turns with too low of ang displacement']);
    try
        key = 'n';
        keyVec  = stateVectors(find(key==cell2mat(keys)),:);
        pauseVec = stateVectors(find('p'==cell2mat(keys)),:);
        walkVec = stateVectors(find('w'==cell2mat(keys)),:);
        wthresh = 4;
        athresh = 1.4;
        for rp = StcCor{s}{key}.data',
            wd = mean(features{s}(rp',16));
            ad = mean(abs(dang{s}(rp')));
            if wd > wthresh && ad < athresh
% REASSIGN turn to walk
                stcMatrix(rp(1):rp(2),:) = repmat(walkVec,[diff(rp)+1,1]);
            end           
            if wd <= wthresh && ad <= athresh
% REASSIGN turn to pause
                stcMatrix(rp(1):rp(2),:) = repmat(pauseVec,[diff(rp)+1,1]);
            end            
        end
    catch err
        disp(err);
    end


    disp([verbosePrefix,'Assign pause to sit if too low and still']);
    % PAUSE speed
    try
        key = 'p';
        keyVec  = stateVectors(find(key==cell2mat(keys)),:);
        sitVec = stateVectors(find('s'==cell2mat(keys)),:);
        walkVec = stateVectors(find('w'==cell2mat(keys)),:);
        bthresh = 5;
        wthresh = 0.2;
        hthresh = 86;
        for rp = StcCor{s}{key}.data',
            wh = mean(features{s}(rp',13));
            wb = features{s}(round(sum(rp)/2),16);
            wba = mean(abs(features{s}(rp',16)));
            if wba < wthresh && wh < hthresh,
% REASSIGN pause to sit
                stcMatrix(rp(1):rp(2),:) = repmat(sitVec,[diff(rp)+1,1]);
            end
% $$$             if wb(end)>bthresh,
% $$$                 stcMatrix(rp(1):rp(2),:) = repmat(walkVec,[diff(rp)+1,1]);
% $$$             end
        end
    catch err
        disp(err);
    end


    disp([verbosePrefix,'Assign groom to sit if too low and still']);
    try
        key = 'm';
        keyVec  = stateVectors(find(key==cell2mat(keys)),:);
        sitVec = stateVectors(find('s'==cell2mat(keys)),:);
        pauseVec = stateVectors(find('p'==cell2mat(keys)),:);        
        hthresh = 86;
        wthresh = 0.2;
        gthresh = 140;
        for rp = StcCor{s}{key}.data',
            wh   = mean(features{s}(rp',13));                    
            wd   = mean(abs(features{s}(rp',16)));
            gd   = mean(sqrt(sum(features{s}(rp',[1,9]).^2,2)));
            if wd < wthresh && wh < hthresh,
% REASSIGN groom to sit                
                stcMatrix(rp(1):rp(2),:) = repmat(sitVec,[diff(rp)+1,1]);
            end
            if gd > gthresh,
                stcMatrix(rp(1):rp(2),:) = repmat(pauseVec,[diff(rp)+1,1]);
            end
        end
    catch err
        disp(err);
    end
    

    disp([verbosePrefix,'Assign sit to pause if too fast']);
    try
        key = 's';
        keyVec  = stateVectors(find(key==cell2mat(keys)),:);
        pauseVec = stateVectors(find('p'==cell2mat(keys)),:);
        wthresh = .3;
        %wdthresh = 2.*xyz{s}.sampleRate;
        for rp = StcCor{s}{key}.data',
            wd   = mean(abs(features{s}(rp',16)));
            if wd > wthresh, % & wdur < wdthresh, 
% REASSIGN sit to pause
                stcMatrix(rp(1):rp(2),:) = repmat(pauseVec,[diff(rp)+1,1]);
            end
        end
    catch err
        disp(err);
    end
    
    StcCor{s} = mat2stc(stcMatrix,StcCor{s},features{s},Trials{s},states,keys);
    [stcMatrix] = stc2mat(StcCor{s},xyz{s},states);    
    [stcMatrix] = reassign_low_duration_state_to_neighboring_states(stcMatrix,StcCor{s}{'s'},5*xyz{s}.sampleRate);    
    [stcMatrix] = reassign_low_duration_state_to_neighboring_states(stcMatrix,StcCor{s}{'m'},1.5*xyz{s}.sampleRate);        
    [stcMatrix] = reassign_low_duration_state_to_neighboring_states(stcMatrix,StcCor{s}{'p'},0.2*xyz{s}.sampleRate);    

    StcCor{s} = mat2stc(stcMatrix,StcCor{s},features{s},Trials{s},states,keys);
    [stcMatrix] = stc2mat(StcCor{s},xyz{s},states);    
    [stcMatrix] = reassign_low_duration_state_to_neighboring_states(stcMatrix,StcCor{s}{'s'},5*xyz{s}.sampleRate);    

    StcCor{s} = mat2stc(stcMatrix,StcCor{s},features{s},Trials{s},states,keys);
    
end



% $$$ 
% $$$ disp([verbosePrefix,'Adjust rear onset boundary']);
% $$$ % ALL to rear
% $$$ StcCor = adjust_state_boundaries_svd(StcCor, Trials,                             ...
% $$$                                      struct('sessionList',            'hand_labeled',                ...
% $$$                                             'referenceTrial',         'jg05-20120317.cof.all',       ...
% $$$                                             'featureSet',             'fet_bref',                    ...
% $$$                                             'sampleMode',             'trimmed',                     ...
% $$$                                             'svdState',               'rear',                        ...
% $$$                                             'antecedentState',        'gper-rear',                   ...
% $$$                                             'subsequentState',        'rear',                        ...
% $$$                                             'immutableStates',        {{}},                          ...
% $$$                                             'sampleRate',             119.881035,                    ...
% $$$                                             'eigenVectorFeaturesMask',{{[6:10,16:25],[1:10,16:25]}}, ...
% $$$                                             'eigenVectorTemporalMask',[1:15,46:64],                  ...
% $$$                                             'eigenVectorIndices',     [1,2],                         ...
% $$$                                             'sortTurnsIndex',         [false,false],                 ...
% $$$                                             'embeddingWindow',        64,                            ...
% $$$                                             'regressionWindow',       100:181,                       ...
% $$$                                             'regressionThreshold',    5e3,                           ...
% $$$                                             'residualSearchWindow',   0.5,                           ...
% $$$                                             'medianCorrectionOffset', 0.14),                         ...
% $$$                                      [],                                                             ...
% $$$                                      false,                                                          ...
% $$$                                      false,                                                          ...
% $$$                                      false                                                           ...
% $$$                                      );
% $$$ 
% $$$ disp([verbosePrefix,'Adjust rear offset boundary']);
% $$$ % REAR to all
% $$$ StcCor = adjust_state_boundaries_svd(StcCor,Trials,                              ...
% $$$                                      struct('sessionList',            'hand_labeled',                ...
% $$$                                             'referenceTrial',         'jg05-20120317.cof.all',       ...
% $$$                                             'featureSet',             'fet_bref',                    ...
% $$$                                             'sampleMode',             'trimmed',                     ...
% $$$                                             'svdState',               'rear',                        ...
% $$$                                             'antecedentState',        'rear',                        ...
% $$$                                             'subsequentState',        'gper-rear',                   ...
% $$$                                             'immutableStates',        {{}},                          ...
% $$$                                             'sampleRate',             119.881035,                    ...
% $$$                                             'eigenVectorFeaturesMask',{{[6:10,16:25],[1:10,16:25]}}, ...
% $$$                                             'eigenVectorTemporalMask',[1:15,46:64],                  ...
% $$$                                             'eigenVectorIndices',     [1,2],                         ...
% $$$                                             'sortTurnsIndex',         [false,false],                 ...
% $$$                                             'embeddingWindow',        64,                            ...
% $$$                                             'regressionWindow',       50:150,                        ...
% $$$                                             'regressionThreshold',    5e3,                           ...
% $$$                                             'residualSearchWindow',   0.5,                           ...
% $$$                                             'medianCorrectionOffset', -0.14),                        ...
% $$$                                      [],                                                             ...
% $$$                                      false,                                                          ...
% $$$                                      false,                                                          ...
% $$$                                      false                                                           ...
% $$$                                      );

% $$$ disp([verbosePrefix,'Adjust turn onset boundary from pause']);
% $$$ % PAUSE to TURN -------------------------------------------------------------------------------
% $$$ StcCor = adjust_state_boundaries_svd(StcCor,Trials,                              ...
% $$$                                      struct('sessionList',            'hand_labeled',                ...
% $$$                                             'referenceTrial',         'jg05-20120317.cof.all',       ...
% $$$                                             'featureSet',             'fet_bref',                    ...
% $$$                                             'sampleMode',             'trimmed',                     ...
% $$$                                             'svdState',               'walk+turn',                   ...
% $$$                                             'antecedentState',        'pause',                       ...
% $$$                                             'subsequentState',        'turn',                        ...
% $$$                                             'immutableStates',        {{}},                          ...
% $$$                                             'sampleRate',             119.881035,                    ...
% $$$                                             'eigenVectorFeaturesMask',{{[1:16,18:24,26:30]}},        ...
% $$$                                             'eigenVectorTemporalMask',[1:15,46:64],                  ...
% $$$                                             'eigenVectorIndices',     [2],                           ...
% $$$                                             'sortTurnsIndex',         [true],                        ...
% $$$                                             'embeddingWindow',        64,                            ...
% $$$                                             'regressionWindow',       100:150,                       ...
% $$$                                             'regressionThreshold',    5,                             ...
% $$$                                             'residualSearchWindow',   0.25,                          ...
% $$$                                             'medianCorrectionOffset', 0),                            ...
% $$$                                      [],                                                             ...
% $$$                                      false,                                                          ...
% $$$                                      false,                                                          ...
% $$$                                      false                                                           ...
% $$$                                      );
% $$$ %----------------------------------------------------------------------------------------------
% $$$ 
% $$$ 
% $$$ 
% $$$ disp([verbosePrefix,'Adjust walk onset boundary from pause']);
% $$$ % PAUSE to walk -------------------------------------------------------------------------------
% $$$ StcCor = adjust_state_boundaries_svd(StcCor,Trials,                              ...
% $$$                                      struct('sessionList',            'hand_labeled',                ...
% $$$                                             'referenceTrial',         'jg05-20120317.cof.all',       ...
% $$$                                             'featureSet',             'fet_bref',                    ...
% $$$                                             'sampleMode',             'trimmed',                     ...
% $$$                                             'svdState',               'walk+turn',                   ...
% $$$                                             'antecedentState',        'pause',                       ...
% $$$                                             'subsequentState',        'walk',                        ...
% $$$                                             'immutableStates',        {{}},                          ...
% $$$                                             'sampleRate',             119.881035,                    ...
% $$$                                             'eigenVectorFeaturesMask',{{[2:15,17:2:25,26:30]}},      ...
% $$$                                             'eigenVectorTemporalMask',[1:15,46:64],                  ...
% $$$                                             'eigenVectorIndices',     [1],                           ...
% $$$                                             'sortTurnsIndex',         [false],                       ...
% $$$                                             'embeddingWindow',        64,                            ...
% $$$                                             'regressionWindow',       110:150,                       ...
% $$$                                             'regressionThreshold',    1000,                          ...
% $$$                                             'residualSearchWindow',   0.25,                          ...
% $$$                                             'medianCorrectionOffset', 0.08),                         ...
% $$$                                      [],                                                             ...
% $$$                                      false,                                                          ...
% $$$                                      false,                                                          ...
% $$$                                      false                                                           ...
% $$$                                      );


% UPDATE stc mode
cf(@(s,tag) s.updateMode([s.mode,tag]), StcCor,repmat({postprocessingTag},[1,numTrials]));
%cf(@(s,t) s.load(t), StcCor,Trials);

% SAVE stc
cf(@(s) s.save(1), StcCor);


% END MAIN ----------------------------------------------------------------------------------------




% AUX METHODS --------------------------------------------------------------------------------------

function [stcMatrix] = reassign_low_duration_state_to_neighboring_states(stcMatrix,state,duration)
for rp = state.data'
    wdur = diff(rp);
    if wdur < duration,
        try
            startVec = stcMatrix(rp(1)-1,:);
            stopVec  = stcMatrix(rp(2)+1,:);
            midpoint = round(sum(rp)/2);
            stcMatrix(rp(1):midpoint,:) = repmat(startVec,[midpoint-rp(1)+1,1]);
            stcMatrix([midpoint+1]:rp(2),:) = repmat(stopVec, [rp(2)-midpoint,1]);
        catch err
            disp(err);
        end
    end
end

% END AUX METHODS ----------------------------------------------------------------------------------