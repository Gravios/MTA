function optimize_stc_transition_svd(varargin)



% DEFARGS ----------------------------------------------------------------------
defargs = struct('rlist',  'training_hand_labeled',                          ...
                 'slist',  {{'hand_labeled_jg';'hand_labeled_Ed'}},          ...
                 'fetSet', 'fet_bref_exp',                                   ...
                 'tag_preprocessing', '+seh+',                               ...
                 'tag_postprocessing','_PPSVD',                              ...
                 'sampleRate', 15,                                           ...
                 'nNeurons',   200,                                          ...
                 'nIter',      100,                                          ...
                 'states',    {{'walk','rear','turn','pause','groom','sit'}},...
                 'rndMethod', 'WSBNT',                                       ...
                 'norm',       true,                                         ...
                 'mref',       true,                                         ...
                 'prctTrain', []                                             ...                 
);%-----------------------------------------------------------------------------

% varargin = {};

[rlist,slist,fetSet,tag_preprocessing,tag_postprocessing,sampleRate,...
 nNeurons,nIter,states,rndMethod,norm,mref,prctTrain,...
] = DefaultArgs(varargin,defargs,'--struct');

if ~isempty(prctTrain), 
    prctTrainTag = ['_PRT_',prctTrain]; 
else, 
    prctTrainTag = ''; 
end                    

rlist = SessionList(rlist);




% MAIN -------------------------------------------------------------------------

for rli = 1:numel(rlist),

    model = ['MTAC_BATCH-' tag_preprocessing fetSet                                           ...
             '_SR_' num2str(sampleRate)                                                       ...
             '_NORM_' num2str(norm)                                                           ...         
             '_REF_' rlist(rli).sessionName, '.' rlist(rli).mazeName '.' rlist(rli).trialName ...
             '_STC_' rlist(rli).stcMode                                                       ...
             '_NN_' num2str(nNeurons)                                                         ...
             '_NI_' num2str(nIter)                                                            ...
             prctTrainTag                                                                     ...
             '_NN_multiPN_RAND_' rndMethod];
    
    refSession = MTATrial.validate(rlist(rli));
    rfet = feval(fetSet,refSession,sampleRate,false);
    [~,refMean,refStd] = nunity(rfet(refSession.stc{'a'},:));

    for sli = 1:numel(slist),

        SesList = SessionList(slist{sli});
        mapped = '-map2ref';
        ds = load(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,mapped,'.mat']));

        %% Pre-process features
        for s = 1:numel(SesList);
            Trial = MTATrial.validate(SesList(s));


            if isempty(Trial.fet),
                Trial.fet = MTADfet(Trial.spath,...
                                    [],...
                                    [],...
                                    [],...
                                    Trial.sync.copy,...
                                    Trial.sync.data(1),...
                                    []);                  
            end
            

            features = fet_bref(Trial);
            features.map_to_reference_session(Trial,refSession);
            %features.filter('ButFilter',3,1,'low');
            %features.unity([],refMean,refStd);

% GET normalization parameters


            
            xyz = Trial.load('xyz');
            xyz.filter('ButFilter',5,1.5,'low');
            vxy = xyz.vel({'spine_lower','head_back'},[1,2]);

            ang = create(MTADang,Trial,xyz);
            dang = xyz.copy('empty');
            dang.data = circ_dist(circshift(ang(:,'spine_lower','spine_upper',1),-1),...
                                  circshift(ang(:,'spine_lower','spine_upper',1),1))...
                                  *features.sampleRate;
            
            StcHL = Trial.stc.copy;
            shl = MTADxyz('data',double(~~stc2mat(StcHL,features,states)),'sampleRate',features.sampleRate);
            ysm = MTADxyz('data',double(0<stc2mat(ds.stc{s},features,states)),'sampleRate',features.sampleRate); 

            StcCor = ds.stc{s}.copy();
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

            disp(['optimize_stc_transition_svd: ',Trial.filebase]);
            disp('Assign all periods with high mean heights to rear');
% REAR remove periods which have low mean heights 
            for key = 'wnpms',
                rhh = [];
                rthresh = 120;
                tds = ds.d_state{s};
                [tpv,tps] = sort(ds.p_state{s},2,'descend');
                try
                    for rp = ds.stc{s}{key}.data',
                        rhh(end+1) = max(features(rp',14));
                        if rhh(end)>rthresh,
                            pind = rp(1):rp(2);
                            tds(pind,2) = 0;
                            tps(pind,:) = StcCor.gsi('rear');
                            % maybe make a cat function for MTADepoch 
                            StcCor.states{StcCor.gsi('rear')}.data = ...
                                [StcCor.states{StcCor.gsi('rear')}.data;pind([1,end])];
                        end
                    end
                    StcCor.states{StcCor.gsi(key)}.data(rhh>rthresh,:) = [];                

                    StcCor.states{StcCor.gsi(key)} = StcCor{key}+[-0,0];

                    % clear other states from the timeperiods assigned
                    % to the state
                    for sts = StcCor.states,
                        if strcmp(sts{1}.key,key),continue,end
                        StcCor.states{StcCor.gsi(sts{1}.key)} = sts{1}-StcCor{key}.data; 
                    end
                catch err
                    disp(err)
                end
            end

            disp('REAR remove periods which have too low mean heights');
% REAR remove periods which have low mean heights 
            key = 'r';
            rhh = [];
            rthresh = 120;
            tds = ds.d_state{s};
            [tpv,tps] = sort(ds.p_state{s},2,'descend');
            try
                for rp = ds.stc{s}{key}.data',
                    rhh(end+1) = max(features(rp',14));
                    if rhh(end)<rthresh,
                        pind = rp(1):rp(2);
                        tds(pind,2) = 0;
                        tps(pind,:) = StcCor.gsi('pause');
                        % maybe make a cat function for MTADepoch 
                        StcCor.states{StcCor.gsi('pause')}.data = ...
                            [StcCor.states{StcCor.gsi('pause')}.data;pind([1,end])];
                    end
                end
                StcCor.states{StcCor.gsi(key)}.data(rhh<rthresh,:) = [];                

                StcCor.states{StcCor.gsi(key)} = StcCor{key}+[-0,0];

                % clear other states from the timeperiods assigned
                % to the state
                for sts = StcCor.states,
                    if strcmp(sts{1}.key,key),continue,end
                    StcCor.states{StcCor.gsi(sts{1}.key)} = sts{1}-StcCor{key}.data; 
                end
            catch err
                disp(err)
            end


            % TURN angular displacement
            %try, StcCor = reassign_state_by_duration(StcCor,'n','p',0.1,tds,tps,@lt); end
            disp('Sort out turns with too low of ang displacement');
            try
                key = 'n';
                wd = [];
                ad = [];
                wthresh = 3;
                athresh = 1.2;
                tails = [-0.0,0.1];
                for rp = StcCor{key}.data',
                    wd(end+1) = mean(features(rp',16));
                    ad(end+1) = mean(abs(dang(rp')));
                    if wd(end)>wthresh && ad(end) < athresh
                        pind = rp(1):rp(2);
                        tds(pind,2) = 0;
                        % maybe make a cat function for MTADepoch 
                        StcCor.states{StcCor.gsi('walk')}.data = ...
                            [StcCor.states{StcCor.gsi('walk')}.data;pind([1,end])];
                    end
                    if wd(end)<= wthresh && ad(end) <= athresh
                        pind = rp(1):rp(2);
                        tds(pind,2) = 0;
                        % maybe make a cat function for MTADepoch 
                        StcCor.states{StcCor.gsi('pause')}.data = ...
                            [StcCor.states{StcCor.gsi('pause')}.data;pind([1,end])];
                    end
                end
                StcCor.states{StcCor.gsi(key)}.data(wd <= wthresh & ad <= athresh,:) = [];
                StcCor.states{StcCor.gsi(key)}.data(wd >  wthresh & ad <  athresh,:) = [];
                StcCor.states{StcCor.gsi(key)} = StcCor{key}+tails;
                for sts = StcCor.states, sts{1}.clean; end
                for sts = StcCor.states,
                    if strcmp(sts{1}.key,key),continue,end
                    StcCor.states{StcCor.gsi(sts{1}.key)} = sts{1}-StcCor{key}.data; 
                end
                
            end


            disp('Assign pause to sit if too low and still');
% PAUSE speed
            try
                key = 'p';
                wh = [];
                wb = [];
                bthresh = 3;
                wthresh = 0.2;
                hthresh = 85;
                tails = [-0.0,0.0];
                for rp = StcCor{key}.data',
                    wh(end+1) = mean(features(rp',13));
                    wb(end+1) = features(round(sum(rp)/2),16);
                    wba(end+1) = mean(abs(features(rp',16)));
                    if wba(end)<wthresh && wh(end)<hthresh,
                        pind = rp(1):rp(2);
                        tds(pind,2) = 0;
                        % maybe make a cat function for MTADepoch 
                        StcCor.states{StcCor.gsi('sit')}.data = ...
                            [StcCor.states{StcCor.gsi('sit')}.data;pind([1,end])];
                    end
                    if wb(end)>bthresh,
                        pind = rp(1):rp(2);
                        tds(pind,2) = 0;
                        % maybe make a cat function for MTADepoch 
                        StcCor.states{StcCor.gsi('walk')}.data = ...
                            [StcCor.states{StcCor.gsi('walk')}.data;pind([1,end])];
                    end
                end
                StcCor.states{StcCor.gsi(key)}.data(wba<wthresh & wh<hthresh,:) = [];
                StcCor.states{StcCor.gsi(key)}.data(wb>bthresh,:) = [];
                StcCor.states{StcCor.gsi(key)} = StcCor{key}+tails;
                for sts = StcCor.states, sts{1}.clean; end
                for sts = StcCor.states,
                    if strcmp(sts{1}.key,key),continue,end
                    StcCor.states{StcCor.gsi(sts{1}.key)} = sts{1}-StcCor{key}.data; 
                end
                
            end
            

            disp('Assign groom to sit if too low and still');
% GROOM to sit if low speed
            try
                key = 'm';
                wd = [];
                wh = [];
                hthresh = 85;                
                wthresh = 0.2;
                tails = [-0.0,0.0];
                for rp = StcCor{key}.data',
                    wh(end+1) = mean(features(rp',13));                    
                    wd(end+1) = mean(abs(fet(rp',16)));
                    if wd(end) < wthresh && wh(end) < hthresh,
                        pind = rp(1):rp(2);
                        tds(pind,2) = 0;
                        StcCor.states{StcCor.gsi('sit')}.data = ...
                            [StcCor.states{StcCor.gsi('sit')}.data;pind([1,end])];
                    end
                end
                
                StcCor.states{StcCor.gsi(key)}.data(wd<wthresh && wh(end)<hthresh,:) = [];
                
                StcCor.states{StcCor.gsi(key)} = StcCor{key}+tails;

                for sts = StcCor.states, sts{1}.clean; end

                for sts = StcCor.states,
                    if strcmp(sts{1}.key,key),continue,end
                    StcCor.states{StcCor.gsi(sts{1}.key)} = sts{1}-StcCor{key}.data; 
                end
                
            end

            

            %% speed: walk
            try
                key = 'w';
                wd = [];
                wthresh = -0.2;
                tails = [-0.1,0.1];
                for rp = StcCor{key}.data',
                    wd(end+1) = mean(features(rp',16));
                    if wd(end)<wthresh,
                        pind = rp(1):rp(2);
                        tds(pind,2) = 0;
                        % maybe make a cat function for MTADepoch 
                        StcCor.states{StcCor.gsi('pause')}.data = ...
                            [StcCor.states{StcCor.gsi('pause')}.data;pind([1,end])];
                    end
                end
                StcCor.states{StcCor.gsi(key)}.data(wd<wthresh,:) = [];
                StcCor.states{StcCor.gsi(key)} = StcCor{key}+tails;
                for sts = StcCor.states, sts{1}.clean; end
                for sts = StcCor.states,
                    if strcmp(sts{1}.key,key),continue,end
                    StcCor.states{StcCor.gsi(sts{1}.key)} = sts{1}-StcCor{key}.data; 
                end
            end

            try, StcCor = reassign_state_by_duration(StcCor,'w','p',0.1,tds,tps,@lt); end
            try, StcCor = reassign_state_by_duration(StcCor,'s','p',2.5,tds,tps,@lt); end
            try, StcCor = reassign_state_by_duration(StcCor,'m','p',1,tds,tps,@lt); end

            disp('Adjust rear onset boundary');
% ALL to rear
            StcCor = adjust_state_boundaries_svd(StcCor, Trial,                             ...
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

            disp('Adjust rear offset boundary');            
% REAR to all
            StcCor = adjust_state_boundaries_svd(StcCor,Trial,                              ...
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

            disp('Adjust turn onset boundary from pause');            
% PAUSE to TURN -------------------------------------------------------------------------------
            StcCor = adjust_state_boundaries_svd(StcCor,Trial,                              ...
                            struct('sessionList',            'hand_labeled',                ...
                                   'referenceTrial',         'jg05-20120317.cof.all',       ...
                                   'featureSet',             'fet_bref',                    ...
                                   'sampleMode',             'trimmed',                     ...
                                   'svdState',               'walk+turn',                   ...
                                   'antecedentState',        'pause',                       ...
                                   'subsequentState',        'turn',                        ...
                                   'immutableStates',        {{}},                          ...
                                   'sampleRate',             119.881035,                    ...
                                   'eigenVectorFeaturesMask',{{[1:16,18:24,26:30]}},        ...
                                   'eigenVectorTemporalMask',[1:15,46:64],                  ...
                                   'eigenVectorIndices',     [2],                           ...
                                   'sortTurnsIndex',         [true],                        ...
                                   'embeddingWindow',        64,                            ...
                                   'regressionWindow',       100:150,                       ...
                                   'regressionThreshold',    5,                             ...
                                   'residualSearchWindow',   0.25,                          ...
                                   'medianCorrectionOffset', 0),                            ...
                            [],                                                             ...
                            false,                                                          ...
                            false,                                                          ...
                            false                                                           ...
            );
%----------------------------------------------------------------------------------------------



            disp('Adjust walk onset boundary from pause');            
% PAUSE to walk -------------------------------------------------------------------------------
            StcCor = adjust_state_boundaries_svd(StcCor,Trial,                              ...
                            struct('sessionList',            'hand_labeled',                ...
                                   'referenceTrial',         'jg05-20120317.cof.all',       ...
                                   'featureSet',             'fet_bref',                    ...
                                   'sampleMode',             'trimmed',                     ...
                                   'svdState',               'walk+turn',                   ...
                                   'antecedentState',        'pause',                       ...
                                   'subsequentState',        'walk',                        ...
                                   'immutableStates',        {{}},                          ...
                                   'sampleRate',             119.881035,                    ...
                                   'eigenVectorFeaturesMask',{{[2:15,17:2:25,26:30]}},      ...
                                   'eigenVectorTemporalMask',[1:15,46:64],                  ...
                                   'eigenVectorIndices',     [1],                           ...
                                   'sortTurnsIndex',         [false],                       ...
                                   'embeddingWindow',        64,                            ...
                                   'regressionWindow',       110:150,                       ...
                                   'regressionThreshold',    1000,                          ...
                                   'residualSearchWindow',   0.25,                          ...
                                   'medianCorrectionOffset', 0.08),                         ...
                            [],                                                             ...
                            false,                                                          ...
                            false,                                                          ...
                            false                                                           ...
            );

            StcCor = StcCor{1};
            
            csm = [];
            csm = MTADxyz('data',double(0<stc2mat(StcCor,features,states)),'sampleRate',features.sampleRate);
            labelingEpochs = Trial.stc{'a'}.cast('TimeSeries');

            errorEpochs = Trial.stc{'e'}.cast('TimeSeries');
            if isempty(errorEpochs)
                errorEpochs = labelingEpochs.copy;
                errorEpochs.data = ~errorEpochs.data;
            end
            
            %ind = cstc<0.05;
            ind = MTADepoch('data',any(csm.data,2)&any(shl.data,2),...
                            'sampleRate',features.sampleRate,...
                            'type','TimeSeries');
            ind.resample(labelingEpochs);
            errorEpochs.resample(labelingEpochs);
            ind = ind.data&labelingEpochs.data&~errorEpochs.data;

            tcm = confmat(shl(ind,:),csm(ind,:)); % #DEP: netlab
            ls{s}.confusionMatrix = round(tcm./features.sampleRate,2);
            ls{s}.precision = round(diag(tcm)./sum(tcm,2),4)'.*100;
            ls{s}.sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
            ls{s}.accuracy = sum(diag(tcm))/sum(tcm(:));
            
            
% $$$             stats(end+1,:) = [ls{s}.accuracy,ls{s}.precision,ls{s}.sensitivity];

            StcCor.updateMode([StcCor.mode,tag_postprocessing]);
            stc{s} = StcCor.copy;
% $$$             StcCor.updateMode(['nn_',refSession.filebase(1:4),'_svdc']);
% $$$             StcCor.save(1);

% $$$             sp(4) = subplot(414);
% $$$             plotSTC(StcCor,features,[],states,'brgcmy');
% $$$             linkaxes(sp,'x');


% $$$ figure();
% $$$ sp = [];
% $$$ sp(end+1)=subplot2(6,1,1:2,1);
% $$$ plotSTC(StcHL,1);
% $$$ sp(end+1)=subplot2(6,1,3:4,1);
% $$$ plotSTC(ds.stc{s},1);
% $$$ sp(end+1)=subplot2(6,1,5:6,1);
% $$$ plotSTC(StcCor,1);
% $$$ linkaxes(sp,'xy');

            
            
        end


        mapped = '-map2ref';

        save(fullfile(MTASession().path.data,'analysis',...
                      [slist{sli},'-',model,tag_postprocessing,mapped,'.mat']),...
             '-v7.3','slist','rlist','nNeurons','nIter','sampleRate','model',...
             'fetSet','rndMethod','states','stc','ls');


    end
end



% END MAIN -------------------------------------------------------------------------







