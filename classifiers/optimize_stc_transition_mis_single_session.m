function stc = optimize_stc_transition_mis_single_session(Trial,Stc,varargin)

% DEFARGS ----------------------------------------------------------------------
defargs = struct('RefTrial',  'jg05-20120317.cof.all',                       ...
                 'refStcMode','hand_labeled_rev3_jg',                        ...
                 'fetSet', 'fet_mis',                                        ...
                 'stcMode', '',                                              ...
                 'tag_preprocessing', '+seh+',                               ...
                 'tag_postprocessing','_PP',                                 ...
                 'sampleRate', 12,                                           ...
                 'nNeurons',   100,                                          ...
                 'nIter',      100,                                          ...
                 'states',    {{'walk','rear','turn','pause','groom','sit'}},...
                 'rndMethod', 'WSBNT',                                       ...
                 'norm',       true,                                         ...
                 'mref',    true,                                            ...
                 'prctTrain',  []                                            ...
);%-----------------------------------------------------------------------------

[RefTrial,refStcMode,fetSet,stcMode,tag_preprocessing,tag_postprocessing,sampleRate,...
 nNeurons,nIter,states,rndMethod,norm,mref,prctTrain] = DefaultArgs(varargin,defargs,'--struct');

RefTrial = MTATrial.validate(RefTrial);

if ~isempty(prctTrain), 
    prctTrainTag = ['_PRT_',prctTrain]; 
else, 
    prctTrainTag = ''; 
end                    


% MAIN -------------------------------------------------------------------------

model = ['MTAC_BATCH-' tag_preprocessing fetSet ...
         '_SR_' num2str(sampleRate)             ...
         '_NORM_' num2str(norm)                 ...         
         '_REF_' RefTrial.filebase              ...
         '_STC_' refStcMode                     ...
         '_NN_' num2str(nNeurons)               ...
         '_NI_' num2str(nIter)                  ...   
         prctTrainTag                           ...
         '_NN_multiPN_RAND_' rndMethod];

RefTrial.load('stc',refStcMode);

rfet = feval(fetSet,RefTrial,sampleRate,false);
[~,refMean,refStd] = nunity(rfet(RefTrial.stc{'a'},:));


Trial = MTATrial.validate(Trial);


if mref,
    mapped = '-map2ref';
else
    mapped = '';
end

ds = load(fullfile(Trial.spath,[Trial.filebase,'.','labelBhv_NN','-',model,mapped,'.mat']));

%% Pre-process features
% $$$ features = feval(fetSet,Trial,sampleRate,false);
% $$$ features.map_to_reference_session(Trial,RefTrial);
% $$$ features.unity([],refMean,refStd);

if isempty(Trial.fet),
    Trial.fet = MTADfet(Trial.spath,       ...
                        [],                ...
                        [],                ...
                        [],                ...
                        Trial.sync.copy,   ...
                        Trial.sync.data(1),...
                        []);                  
end


xyz = Trial.load('xyz','seh');
ss = Trial.load('fet','3dssh');
ss.resample(xyz);

xyz.addMarker('bcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));

fvelxy = xyz.copy;
fvelxy.filter('ButFilter',3,2.4,'low');
fvelxy = fvelxy.vel([],[1,2]);
fvelxy.data(fvelxy.data<=10e-5) = 10e-5;
fvelxy.data = log10(fvelxy.data);

ang = create(MTADang,Trial,xyz);

if ~isempty(stcMode),
    StcHL = Trial.load('stc',stcMode);
    shl = MTADxyz('data',double(~~stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
end

%% COMPOSITE: rear
StcCor = ds.stc.copy;
rhh = [];
rthresh = 20;
tds = ds.d_state;
[tpv,tps] = sort(ds.p_state,2,'descend');
try
    for rp = StcCor{'r'}.data',
        rhh(end+1) = max(xyz(rp','hcom',3).*ang(rp','spine_middle','spine_upper',2));
        if rhh(end)<rthresh,
            pind = rp(1):rp(2);
            tds(pind,2) = 0;
            tps(pind,:) = circshift(tps(pind,:),-1,2);
            % maybe make a cat function for MTADepoch 
            StcCor.states{StcCor.gsi(states{mode(tps(pind,1))})}.data = ...
                [StcCor.states{StcCor.gsi(states{mode(tps(pind,1))})}.data;pind([1,end])];
        end
    end
    for sts = StcCor.states, sts{1}.clean; end
    StcCor.states{StcCor.gsi('r')}.data(rhh<rthresh,:) = [];
end

%% DISTANCE: Walk
try
    key = 'w';
    wd = [];
    wthresh = 1.5;
    for rp = StcCor{key}.data',
        wd(end+1) = log10(sqrt(sum([xyz(rp(2),'spine_lower',[1,2])-xyz(rp(1),'spine_lower',[1,2])].^2,3)));
        if wd(end)<wthresh,
            pind = rp(1):rp(2);
            tds(pind,2) = 0;
            tps(pind,:) = circshift(tps(pind,:),-1,2);
            % maybe make a cat function for MTADepoch 
            msts = mode(tps(pind,1));
            if msts == 2; msts = 4; end
            StcCor.states{StcCor.gsi(states{msts})}.data = ...
                [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
        end
    end
    for sts = StcCor.states, sts{1}.clean; end
    StcCor.states{StcCor.gsi(key)}.data(wd<wthresh,:) = [];
end


%% DURATION: Sit
%try, StcCor = reassign_state_by_duration(StcCor,'s','p',5,tds,tps); end
%% DURATION: Groom
%try, StcCor = reassign_state_by_duration(StcCor,'m','p',3,tds,tps); end
%% DURATION: Pause
try, StcCor = reassign_state_by_duration(StcCor,'p', [],0.2,tds,tps); end
%% DURATION: Groom
try, StcCor = reassign_state_by_duration(StcCor,'m','p',3,tds,tps); end
%% DURATION: Sit
try, StcCor = reassign_state_by_duration(StcCor,'s','p',5,tds,tps); end



%% DURATION: Turn
try
    pd = [];
    ad = [];    
    key = 'n';
    pthresh = log10(0.1.*StcCor{key}.sampleRate);
    athresh = 0.2;
    for rp = StcCor{key}.data',
        pd(end+1) = log10(abs(diff(rp)));
        ad(end+1) = abs(circ_dist(ang(rp(2),'spine_lower','spine_upper',1),...
                                  ang(rp(1),'spine_lower','spine_upper',1)));
        if pd(end)<pthresh&ad(end)<athresh,
            pind = rp(1):rp(2);
            tds(pind,2) = 0;
            tps(pind,:) = circshift(tps(pind,:),-1,2);
            % maybe make a cat function for MTADepoch 
            msts = mode(tps(pind,1));
            if msts == 2; 
                tps(pind,:) = circshift(tps(pind,:),-1,2);
                msts = mode(tps(pind,1));
            end                    
            
            StcCor.states{StcCor.gsi(states{msts})}.data = ...
                [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
        end
    end
    
    tind = pd<pthresh&ad<athresh;
    if any(tind),StcCor.states{StcCor.gsi(key)}.data(tind,:) = [];end
    
    StcCor.states{StcCor.gsi(key)} = StcCor{key}+[-0.1,0.1];

    for sts = StcCor.states, sts{1}.clean; end

    for sts = StcCor.states,
        if strcmp(sts{1}.key,key),continue,end
        StcCor.states{StcCor.gsi(sts{1}.key)} = sts{1}-StcCor{key}.data; 
    end    

end




%% Speed: Walk 
try
    key = 'w';
    pd = [];
    ad = [];
    pthresh = log10(0.2.*StcCor{key}.sampleRate);
    athresh = 0.2;

    for rp = StcCor{key}.data',
        %for rp = StcHL{key}.data',
        pd(end+1) = log10(abs(diff(rp)));
        ad(end+1) = abs(circ_dist(ang(rp(2),1,4,1),ang(rp(1),1,4,1)));
        if pd(end)<pthresh&ad(end)>athresh,
            pind = rp(1):rp(2);
            tds(pind,2) = 0;
            tps(pind,:) = circshift(tps(pind,:),-1,2);
            % maybe make a cat function for MTADepoch 
            msts = mode(tps(pind,1));
            if msts == 2; 
                tps(pind,:) = circshift(tps(pind,:),-1,2);
                msts = mode(tps(pind,1));
            end                    
            
            StcCor.states{StcCor.gsi(states{msts})}.data = ...
                [StcCor.states{StcCor.gsi(states{msts})}.data;pind([1,end])];
        end
    end
    for sts = StcCor.states, sts{1}.clean; end
    StcCor.states{StcCor.gsi(key)}.data(pd<pthresh&ad>athresh,:) = [];
end


%% DURATION: Sit
try, StcCor = reassign_state_by_duration(StcCor,'s','p',5,tds,tps); end
%% DURATION: Groom 
try, StcCor = reassign_state_by_duration(StcCor,'m','p',5,tds,tps); end


ls = [];
if ~isempty(stcMode),
    csm = [];
    csm = MTADxyz('data',double(0<stc2mat(StcCor,xyz,states)),'sampleRate',xyz.sampleRate); 
    labelingEpochs = Trial.stc{'a'}.cast('TimeSeries');

    errorEpochs = Trial.stc{'e'}.cast('TimeSeries');
    if isempty(errorEpochs)
        errorEpochs = labelingEpochs.copy;
        errorEpochs.data = ~errorEpochs.data;
    end

    %ind = cstc<0.05;
    ind = MTADepoch('data',any(csm.data,2)&any(shl.data,2),...
                    'sampleRate',xyz.sampleRate,...
                    'type','TimeSeries');
    ind.resample(labelingEpochs);
    errorEpochs.resample(labelingEpochs);
    ind = ind.data&labelingEpochs.data&~errorEpochs.data;

    tcm = confmat(shl(ind,:),csm(ind,:)); % #DEP: netlab
    ls.confusionMatrix = round(tcm./xyz.sampleRate,2);
    ls.precision = round(diag(tcm)./sum(tcm,2),4)'.*100;
    ls.sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
    ls.accuracy = sum(diag(tcm))/sum(tcm(:));

end

StcCor.updateMode([StcCor.mode,tag_postprocessing]);
stc = StcCor.copy;
stc.save(1);


mapped = '-map2ref';

save(fullfile(Trial.spath,...
              [Trial.filebase,'.',mfilename,'-',model,tag_postprocessing,mapped,'.mat']),...
     '-v7.3','nNeurons','nIter','sampleRate','model',...
     'fetSet','rndMethod','states','stc','ls');





function [StcCor,tempDState,tempPState] = reassign_state_by_duration(StcCor,key,defaultState,durationThreshold,tempDState,tempPState)
pd = [];
states = StcCor.list_state_attrib;
pthresh = log10(durationThreshold.*StcCor{key}.sampleRate);
for rp = StcCor{key}.data',

    % collect period durations
    pd(end+1) = log10(abs(diff(rp)));

    if pd(end)<pthresh,
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
StcCor.states{StcCor.gsi(key)}.data(pd<pthresh,:) = [];
