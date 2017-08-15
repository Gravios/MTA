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


MODEL_TYPE = 'multiSesPatNet';
verbPrefix = 'MTA:classifiers:optimize_stc_transition: ';

% DEFARGS ----------------------------------------------------------------------------------------
defArgs = struct('sessionList',                 '',                                            ...
                 'featureSet',                  'fet_bref_rev7',                               ...
                 'states',                      {{'walk','rear','turn','pause','groom','sit'}},...
                 'keys',                        {{'w','r','n','p','m','s'}},                   ...
                 'model',                       [],                                            ...
                 'sampleRate',                  10,                                            ...
                 'nNeurons',                    25,                                            ...
                 'nIter',                       10,                                            ...
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
 normalizationSessionList,dropIndex,prctTrain,postprocessingTag,modelType] = ...
    DefaultArgs(varargin,defArgs,'--struct');
%-----------------------------------------------------------------------------

% varargin = {};




% MAIN ------------------------------------------------------------------------------------------

% LOAD Trials which have the target state collections
Trials = af(@(Trial)  MTATrial.validate(Trial), get_session_list(sessionList));
numTrials = numel(Trials);


% LOAD Stc for opitimization
if isempty(Stc),
% CREATE model name
    if isempty(model),
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
    StcCor = cf(@(t,m) t.load('stc',m),    Trials, repmat({model},[1,numTrials]));
elseif iscell(Stc),    
    StcCor = cf(@(s)   s.copy(),         Stc);
elseif isa(Stc,'MTAStateCollection'),
    StcCor = {Stc.copy()};
else
    error('MTA:classifiers:optimize_stc_transition:UnknownStcSignature');
end


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
    disp([verbPrefix,'Processing: ',Trials{s}.filebase]);

    disp([verbPrefix,'Assign all periods with high mean heights to rear']);

    [stcMatrix] = stc2mat(StcCor{s},xyz{s},states);
    
    for key = 'wnpms',
        keyInd = find(key==cell2mat(keys));
        keyVec = zeros([1,numel(keys)]);
        keyVec(keyInd) = keyInd;
        
        rearInd = find('r'==keys);
        rhh = [];
        rthresh = 140;
        try
            for rp = StcCor{s}{key}.data',
                rhh(end+1) = max(features{s}(rp',14));
                if rhh(end)>rthresh,
                    stcMatrix(rp(1):rp(2),:) = rearVec;

                    StcCor{s}.states{StcCor{s}.gsi('rear')}.data = [StcCor{s}.states{StcCor{s}.gsi('rear')}.data;rp'];
                end
            end
            StcCor{s}.states{StcCor{s}.gsi(key)}.data(rhh>rthresh,:) = [];                
            StcCor{s}.states{StcCor{s}.gsi(key)} = StcCor{s}{key}+[-0,0];
            % clear other states from the timeperiods assigned to the state
            for sts = StcCor{s}.states,
                if strcmp(sts{1}.key,key),continue,end
                StcCor{s}.states{StcCor{s}.gsi(sts{1}.key)} = sts{1}-StcCor{s}{key}.data; 
            end
        catch err
            disp(err)
        end
    end


    disp([verbPrefix,'REAR remove periods which have too low mean heights']);
    % REAR remove periods which have low mean heights 
    key = 'r';
    rhh = [];
    rthresh = 140;
    try
        for rp = StcCor{s}{key}.data',
            rhh(end+1) = max(features{s}(rp',14));
            if rhh(end)<rthresh,
                pind = rp(1):rp(2);
                % maybe make a cat function for MTADepoch 
                StcCor{s}.states{StcCor{s}.gsi('pause')}.data = ...
                    [StcCor{s}.states{StcCor{s}.gsi('pause')}.data;pind([1,end])];
            end
        end
        StcCor{s}.states{StcCor{s}.gsi(key)}.data(rhh<rthresh,:) = [];                
        StcCor{s}.states{StcCor{s}.gsi(key)} = StcCor{s}{key}+[-0,0];
        % clear other states from the timeperiods assigned
        % to the state
        for sts = StcCor{s}.states,
            if strcmp(sts{1}.key,key),continue,end
            StcCor{s}.states{StcCor{s}.gsi(sts{1}.key)} = sts{1}-StcCor{s}{key}.data; 
        end
    catch err
        disp(err)
    end


    % TURN angular displacement
    %try, StcCor{s} = reassign_state_by_duration(StcCor{s},'n','p',0.1,tds,tps,@lt); end
    disp([verbPrefix 'Sort out turns with too low of ang displacement']);
    try
        key = 'n';
        wd = [];
        ad = [];
        wthresh = 1;
        athresh = 1.2;
        tails = [0,0];
        for rp = StcCor{s}{key}.data',
            wd(end+1) = mean(features{s}(rp',16));
            ad(end+1) = mean(abs(dang{s}(rp')));
            if wd(end) > wthresh && ad(end) < athresh
% REASSIGN turn to walk
                StcCor{s}.states{StcCor{s}.gsi('walk')}.data = ...
                    [StcCor{s}.states{StcCor{s}.gsi('walk')}.data;rp'];
            end           
            if wd(end) <= wthresh && ad(end) <= athresh
% REASSIGN turn to pause
                StcCor{s}.states{StcCor{s}.gsi('pause')}.data = ...
                    [StcCor{s}.states{StcCor{s}.gsi('pause')}.data;rp'];
            end
        end
        StcCor{s}.states{StcCor{s}.gsi(key)}.data( (wd <= wthresh & ad <= athresh)...
                                                  |(wd >  wthresh & ad <  athresh),:) = [];
        StcCor{s}.states{StcCor{s}.gsi(key)} = StcCor{s}{key}+tails;
        for sts = StcCor{s}.states, sts{1}.clean; end
        for sts = StcCor{s}.states,
            if strcmp(sts{1}.key,key),continue,end
            StcCor{s}.states{StcCor{s}.gsi(sts{1}.key)} = sts{1}-StcCor{s}{key}.data; 
        end
    catch err
        disp(err)                
    end


    disp([verbPrefix,'Assign pause to sit if too low and still']);
    % PAUSE speed
    try
        key = 'p';
        wh = [];
        wb = [];
        wba = [];                
        bthresh = 4;
        wthresh = 0.4;
        hthresh = 86;
        tails = [-0.0,0.0];
        for rp = StcCor{s}{key}.data',
            wh(end+1) = mean(features{s}(rp',13));
            wb(end+1) = features{s}(round(sum(rp)/2),16);
            wba(end+1) = mean(abs(features{s}(rp',16)));
            if wba(end)<wthresh && wh(end)<hthresh,
                StcCor{s}.states{StcCor{s}.gsi('sit')}.data = ...
                    [StcCor{s}.states{StcCor{s}.gsi('sit')}.data;rp'];
            end
            if wb(end)>bthresh,
                StcCor{s}.states{StcCor{s}.gsi('walk')}.data = ...
                    [StcCor{s}.states{StcCor{s}.gsi('walk')}.data;rp'];
            end
        end
        StcCor{s}.states{StcCor{s}.gsi(key)}.data((wba<wthresh & wh<hthresh)|(wb>bthresh),:) = [];
        StcCor{s}.states{StcCor{s}.gsi(key)} = StcCor{s}{key}+tails;
        for sts = StcCor{s}.states, sts{1}.clean; end
        for sts = StcCor{s}.states,
            if strcmp(sts{1}.key,key),continue,end
            StcCor{s}.states{StcCor{s}.gsi(sts{1}.key)} = sts{1}-StcCor{s}{key}.data; 
        end
        
    end


    disp([verbPrefix,'Assign groom to sit if too low and still']);
    % GROOM to sit if low speed
    try
        key = 'm';
        wd = [];
        wh = [];
        wdur = [];        
        wdthresh = 1.*xyz{s}.sampleRate;        
        hthresh = 86;
        wthresh = 0.2;
        tails = [-0.0,0.0];
        for rp = StcCor{s}{key}.data',
            wdur(end+1) = diff(rp([1,end]));
            wh(end+1) = mean(features{s}(rp',13));                    
            wd(end+1) = mean(abs(features{s}(rp',16)));
            if wd(end) < wthresh && wh(end) < hthresh,
                pind = rp(1):rp(2);
                StcCor{s}.states{StcCor{s}.gsi('sit')}.data = ...
                    [StcCor{s}.states{StcCor{s}.gsi('sit')}.data;pind([1,end])];
            end
            if wdur(end) < wdthresh,
                StcCor{s}.states{StcCor{s}.gsi('pause')}.data = ...
                    [StcCor{s}.states{StcCor{s}.gsi('pause')}.data;pind([1,end])];

            end            
        end
        StcCor{s}.states{StcCor{s}.gsi(key)}.data((wd<wthresh & wh(end)<hthresh)|(wdur<wdthresh),:) = [];
        StcCor{s}.states{StcCor{s}.gsi(key)} = StcCor{s}{key}+tails;
        for sts = StcCor{s}.states, sts{1}.clean; end
        for sts = StcCor{s}.states,
            if strcmp(sts{1}.key,key),continue,end
            StcCor{s}.states{StcCor{s}.gsi(sts{1}.key)} = sts{1}-StcCor{s}{key}.data; 
        end
    catch err
        disp(err)
    end
    

    disp([verbPrefix,'Assign sit to pause if too fast']);
% SPEED sit
    try
        key = 's';
        wd = []; 
        wdur = [];
        wthresh = .4;
        wdthresh = 3.*xyz{s}.sampleRate;
        tails = [0,0];
        for rp = StcCor{s}{key}.data',
            wd(end+1) = mean(abs(features{s}(rp',16)));
            wdur(end+1) = diff(rp([1,end]));
            if (wd(end)>wthresh)&(wdur(end)<wdthresh),
                pind = rp(1):rp(2);
                StcCor{s}.states{StcCor{s}.gsi('pause')}.data = ...
                    [StcCor{s}.states{StcCor{s}.gsi('pause')}.data;pind([1,end])];
            end
        end
        StcCor{s}.states{StcCor{s}.gsi(key)}.data((wd<wthresh)&(wdur<wdthresh),:) = [];
        StcCor{s}.states{StcCor{s}.gsi(key)} = StcCor{s}{key}+tails;
        for sts = StcCor{s}.states, sts{1}.clean; end
        for sts = StcCor{s}.states,
            if strcmp(sts{1}.key,key),continue,end
            StcCor{s}.states{StcCor{s}.gsi(sts{1}.key)} = sts{1}-StcCor{s}{key}.data; 
        end
    catch err
        disp(err)
    end

    try
        key = 'w';
        wd = []; 
        wthresh = 0;
        tails = [0,0];
        for rp = StcCor{s}{key}.data',
            wd(end+1) = mean(abs(features{s}(rp',16)));
            if wd(end)<wthresh,
                pind = rp(1):rp(2);
                StcCor{s}.states{StcCor{s}.gsi('pause')}.data = ...
                    [StcCor{s}.states{StcCor{s}.gsi('pause')}.data;pind([1,end])];
            end
        end
        StcCor{s}.states{StcCor{s}.gsi(key)}.data((wd<wthresh),:) = [];
        StcCor{s}.states{StcCor{s}.gsi(key)} = StcCor{s}{key}+tails;
        for sts = StcCor{s}.states, sts{1}.clean; end
        for sts = StcCor{s}.states,
            if strcmp(sts{1}.key,key),continue,end
            StcCor{s}.states{StcCor{s}.gsi(sts{1}.key)} = sts{1}-StcCor{s}{key}.data; 
        end
    catch err
        disp(err)
    end
    
% $$$     try, StcCor{s} = reassign_state_by_duration(StcCor{s},'w','p',0.1,tds,tps,@lt); end
% $$$     try, StcCor{s} = reassign_state_by_duration(StcCor{s},'s','p',1.5,tds,tps,@lt); end
% $$$     try, StcCor{s} = reassign_state_by_duration(StcCor{s},'m','p',1.5,tds,tps,@lt); end

end


disp([verbPrefix,'Adjust rear onset boundary']);
% ALL to rear
StcCor = adjust_state_boundaries_svd(StcCor, Trials,                             ...
                                     struct('sessionList',            'hand_labeled',                ...
                                            'referenceTrial',         'jg05-20120317.cof.all',       ...
                                            'featureSet',             'fet_bref',                    ...
                                            'sampleMode',             'trimmed',                     ...
                                            'svdState',               'rear',                        ...
                                            'antecedentState',        'gper-rear',                   ...
                                            'subsequentState',        'rear',                        ...
                                            'immutableStates',        {{}},                          ...
                                            'sampleRate',             119.881035,                    ...
                                            'eigenVectorFeaturesMask',{{[6:10,16:25],[1:10,16:25]}}, ...
                                            'eigenVectorTemporalMask',[1:15,46:64],                  ...
                                            'eigenVectorIndices',     [1,2],                         ...
                                            'sortTurnsIndex',         [false,false],                 ...
                                            'embeddingWindow',        64,                            ...
                                            'regressionWindow',       100:181,                       ...
                                            'regressionThreshold',    5e3,                           ...
                                            'residualSearchWindow',   0.5,                           ...
                                            'medianCorrectionOffset', 0.14),                         ...
                                     [],                                                             ...
                                     false,                                                          ...
                                     false,                                                          ...
                                     false                                                           ...
                                     );

disp([verbPrefix,'Adjust rear offset boundary']);
% REAR to all
StcCor = adjust_state_boundaries_svd(StcCor,Trials,                              ...
                                     struct('sessionList',            'hand_labeled',                ...
                                            'referenceTrial',         'jg05-20120317.cof.all',       ...
                                            'featureSet',             'fet_bref',                    ...
                                            'sampleMode',             'trimmed',                     ...
                                            'svdState',               'rear',                        ...
                                            'antecedentState',        'rear',                        ...
                                            'subsequentState',        'gper-rear',                   ...
                                            'immutableStates',        {{}},                          ...
                                            'sampleRate',             119.881035,                    ...
                                            'eigenVectorFeaturesMask',{{[6:10,16:25],[1:10,16:25]}}, ...
                                            'eigenVectorTemporalMask',[1:15,46:64],                  ...
                                            'eigenVectorIndices',     [1,2],                         ...
                                            'sortTurnsIndex',         [false,false],                 ...
                                            'embeddingWindow',        64,                            ...
                                            'regressionWindow',       50:150,                        ...
                                            'regressionThreshold',    5e3,                           ...
                                            'residualSearchWindow',   0.5,                           ...
                                            'medianCorrectionOffset', -0.14),                        ...
                                     [],                                                             ...
                                     false,                                                          ...
                                     false,                                                          ...
                                     false                                                           ...
                                     );

% $$$ disp([verbPrefix,'Adjust turn onset boundary from pause']);
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
% $$$ disp([verbPrefix,'Adjust walk onset boundary from pause']);
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


% END MAIN -------------------------------------------------------------------------







