

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
                 'overwriteNormalizationParameters',false                    ...
);%-----------------------------------------------------------------------------


[rlist,slist,fetSet,tag_preprocessing,tag_postprocessing,sampleRate,...
 nNeurons,nIter,states,rndMethod,norm,mref,prctTrain,...
overwriteNormalizationParameters] = DefaultArgs(varargin,defargs,'--struct');

if ~isempty(prctTrain), 
    prctTrainTag = ['_PRT_',prctTrain]; 
else, 
    prctTrainTag = ''; 
end                    

rlist = SessionList(rlist);

if overwriteNormalizationParameters,
    fetNormParam = set_normalization_parameters();
else,
    fetNormParam = query_normalization_parameters();
end




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
            %features.unity([],refMean,refStd);

% GET normalization parameters


            
            xyz = Trial.load('xyz');
            xyz.filter('ButFilter',5,1,'low');
            vxy = xyz.vel({'spine_lower','head_back'},[1,2]);

            ang = create(MTADang,Trial,xyz);
            dang = xyz.copy('empty');
            dang.data = circ_dist(circshift(ang(:,'spine_lower','spine_upper',1),-1),...
                                  circshift(ang(:,'spine_lower','spine_upper',1),1))...
                                  *features.sampleRate;
            
            StcHL = Trial.stc.copy;
            shl = MTADxyz('data',double(~~stc2mat(StcHL,features,states)),'sampleRate',features.sampleRate);
            ysm = MTADxyz('data',double(0<stc2mat(ds.stc{s},features,states)),'sampleRate',features.sampleRate); 

            StcCor = ds.stc{s}.copy;

% REAR remove periods which have low mean heights 
            key = 'r';
            rhh = [];
            rthresh = 140;
            tds = ds.d_state{s};
            [tpv,tps] = sort(ds.p_state{s},2,'descend');
            try
                for rp = ds.stc{s}{key}.data',
                    rhh(end+1) = max(features(rp',15));
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

                StcCor.states{StcCor.gsi(key)} = StcCor{key}+[-0.15,0.0];

                % clear other states from the timeperiods assigned
                % to the state
                for sts = StcCor.states,
                    if strcmp(sts{1}.key,key),continue,end
                    StcCor.states{StcCor.gsi(sts{1}.key)} = sts{1}-StcCor{key}.data; 
                end
            catch err
                disp(err)
            end


% ADJUST rear offsets
            adjust_state_boundaries_svd(Trial,StcCor,...
                                        struct('sessionList',            'hand_labeled',...
                                               'referenceTrial',         'jg05-20120317.cof.all',...
                                               'featureSet',             'fet_bref',...
                                               'sampleMode',             'centered',...
                                               'svdState',               'rear',...
                                               'antecedentState',        'gper-rear',...
                                               'subsequentState',        'rear',...
                                               'eigenVectorFeaturesMask',{{[6:10,16:25],[1:10,16:25]}},...
                                               'eigenVectorTemporalMask',[1:15,46:64],...
                                               'eigenVectorIndices',     [1,2],...
                                               'sampleRate',             119.881035,...
                                               'embeddingWindow',        64, ...                    
                                               'regressionWindow',       100:181, ...
                                               'regressionThreshold',    100 )
            );
            
            
% ADJUST rear offsets
            adjust_state_boundaries_svd(Trial,StcCor,...
                                        struct('sessionList',            'hand_labeled',...
                                               'referenceTrial',         'jg05-20120317.cof.all',...
                                               'featureSet',             'fet_bref',...
                                               'sampleMode',             'centered',...
                                               'svdState',               'rear',...
                                               'antecedentState',        'rear',...
                                               'subsequentState',        'gper-rear',...
                                               'eigenVectorFeaturesMask',{{[6:10,16:25],[1:10,16:25]}},...
                                               'eigenVectorTemporalMask',[1:15,46:64],...
                                               'eigenVectorIndices',     [1,2],...
                                               'sampleRate',             119.881035,...
                                               'embeddingWindow',        64, ...                    
                                               'regressionWindow',       100:181, ...
                                               'regressionThreshold',    100 )
            );

            
            
            % 'pause' 'walk' evi = 1
            % 'pause' 'turn' evi = 2
            % 'turn'  'walk' evi = 1,2
            param = struct('sessionList',            'hand_labeled',...
                           'referenceTrial',         'jg05-20120317.cof.all',...
                           'sampleMode',             'centered',...
                           'svdState',               'walk+turn',...
                           'antecedentState',        'pause',...
                           'subsequentState',        'walk',...
                           'eigenVectorFeaturesMask',{{[2:15,17:2:25,26:30],[1:16,18:24,26:30]}},...
                           'eigenVectorTemporalMask',[1:15,46:64],...
                           'eigenVectorIndices',     [1],...
                           'sampleRate',             119.881035,...
                           'embeddingWindow',        64, ...
                           'regressionWindow',       91:151, ...
                           'regressionThreshold',    200 ...
            );
                        
            
            adjust_state_boundaries_svd(Trial,StcCor);


            
            

% TURN angular displacement
            try
                key = 'n';
                wd = [];
                ad = [];
                wthresh = 3;
                athresh = 1.2;
                tails = [-0.0,0.0];
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
            
            

            %% speed: pause
            try
                key = 'p';
                wd = [];
                wh = [];
                wb = [];
                bthresh = 4;
                wthresh = 0.3;
                hthresh = 90;
                tails = [-0.0,0.0];
                for rp = StcCor{key}.data',
                    wd(end+1) = mean(vxy(rp',1));
                    wh(end+1) = mean(features(rp',13));
                    wb(end+1) = mean(features(rp',16));                    
                    if wd(end)<wthresh && wh(end)<hthresh,
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
                
                StcCor.states{StcCor.gsi(key)}.data(wd<wthresh & wh<hthresh,:) = [];
                StcCor.states{StcCor.gsi(key)}.data(wb>bthresh,:) = [];
                
                StcCor.states{StcCor.gsi(key)} = StcCor{key}+tails;

                for sts = StcCor.states, sts{1}.clean; end

                for sts = StcCor.states,
                    if strcmp(sts{1}.key,key),continue,end
                    StcCor.states{StcCor.gsi(sts{1}.key)} = sts{1}-StcCor{key}.data; 
                end
                
            end
            
                        %% speed: pause
            try
                key = 'm';
                wd = [];
                wthresh = 0.5;
                tails = [-0.0,0.0];
                for rp = StcCor{key}.data',
                    wd(end+1) = mean(mean(vxy(rp(1):rp(2),2)));
                    if wd(end)<wthresh,
                        pind = rp(1):rp(2);
                        tds(pind,2) = 0;
                        % maybe make a cat function for MTADepoch 
                        StcCor.states{StcCor.gsi('sit')}.data = ...
                            [StcCor.states{StcCor.gsi('sit')}.data;pind([1,end])];
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

            

            %% speed: walk
            try
                key = 'w';
                wd = [];
                wthresh = 0;
                tails = [-0.0,0.1];
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
            
            

            try, StcCor = reassign_state_by_duration(StcCor,'n','w',0.2,tds,tps,@lt); end
            try, StcCor = reassign_state_by_duration(StcCor,'w','p',0.2,tds,tps,@lt); end
            try, StcCor = reassign_state_by_duration(StcCor,'s','p',1.5,tds,tps,@lt); end
            try, StcCor = reassign_state_by_duration(StcCor,'m','p',1,tds,tps,@lt); end
            %try, StcCor = reassign_state_by_duration(StcCor,'p','w',0.2,tds,tps,@lt); end
            
            
            
            
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

            StcCor.updateMode([StcCor.mode,tag_postprocessing]);
            stc{s} = StcCor.copy;

% $$$             sp(4) = subplot(414);
% $$$             plotSTC(StcCor,features,[],states,'brgcmy');
% $$$             linkaxes(sp,'x');
            
            
            
        end


        mapped = '-map2ref';

        save(fullfile(MTASession().path.data,'analysis',...
                      [slist{sli},'-',model,tag_postprocessing,mapped,'.mat']),...
             '-v7.3','slist','rlist','nNeurons','nIter','sampleRate','model',...
             'fetSet','rndMethod','states','stc','ls');


    end
end



% END MAIN -------------------------------------------------------------------------




% AUX METHODS ----------------------------------------------------------------------

function [StcCor,tempDState,tempPState] = reassign_state_by_duration(StcCor,key,defaultState,durationThreshold,tempDState,tempPState,logicFun)
pd = [];
states = StcCor.list_state_attrib;
pthresh = log10(durationThreshold.*StcCor{key}.sampleRate);
for rp = StcCor{key}.data',

    % collect period durations
    pd(end+1) = log10(abs(diff(rp)));

    if logicFun(pd(end),pthresh),
        pind = rp(1):rp(2);

        %???
        tempDState(pind,2) = 0;        
        % promote next best state 
        tempPState(pind,:) = circshift(tempPState(pind,:),-1,2);

        if isempty(defaultState), 
            % select next best state             
            msts = mode(tempPState(pind,1));
            if msts == 2; 
                tps(pind,:) = circshift(tempPState(pind,:),-1,2);
                msts = mode(tempPState(pind,1));
            end                    
        else % relabel state with a provided default state            
            msts = StcCor.gsi(defaultState);
        end
        
        
        % reassign state
        StcCor.states{StcCor.gsi(states{msts})}.data = ...
            [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
    end
end

% Resort and clean overlaps 
for sts = StcCor.states,
    sts{1}.clean; 
end

% Delete the reassigned periods within original state
StcCor.states{StcCor.gsi(key)}.data(logicFun(pd,pthresh),:) = [];






% END function set_svd_parameters ------------------------------------------------------------------




% END AUX METHODS ----------------------------------------------------------------------



